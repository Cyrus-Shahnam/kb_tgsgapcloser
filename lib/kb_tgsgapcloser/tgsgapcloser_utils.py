# -*- coding: utf-8 -*-
"""Tool logic for kb_tgsgapcloser (TGS-GapCloser v1.2.1, bioconda build).

Pure functions (validation, command building, running, parsing, report files)
are kept apart from KBase I/O so they can be unit-tested and dry-run outside
the container (tools/dry_run_local.py). The SDK Impl only calls Runner.run().

Behaviour verified from the v1.2.1 driver script and binaries, not the README:
  * print_fatal() exits with status 0, so the log is checked for FATAL/Error
    lines and 'ALL DONE' instead of trusting the exit code.
  * done_step*_tag files in the cwd make a re-run silently skip steps, so
    every run gets a fresh directory and stale tags are refused.
  * outputs are <prefix>.scaff_seqs and <prefix>.gap_fill_detail (the README
    says .scaff_seq / .gap_fill_details); both spellings are accepted.
  * reads must be FASTA; minimap2 is called by bare name (PATH injection).
  * with --ne the candidate and correction steps are skipped, so min_nread,
    max_candidate, chunk and r_round have no effect there.
"""
import collections
import gzip
import html
import json
import logging
import os
import re
import shlex
import shutil
import subprocess
import time
import uuid
import zipfile

TOOL_ENV_BIN = os.environ.get('TGS_ENV_BIN', '/opt/conda/envs/tgsgapcloser/bin')
EXPECTED_TOOL_VERSION = '1.2.1'
OUT_PREFIX = 'tgs'

REQUIRED = ('workspace_name', 'assembly_ref', 'reads_ref', 'output_assembly_name')
DEFAULTS = collections.OrderedDict([
    ('read_type', 'ont'), ('error_correction', 'racon'), ('racon_rounds', 1),
    ('minimap2_preset', 'auto'), ('min_idy', None), ('min_match', None),
    ('min_nread', 1), ('max_candidate', 10), ('chunk', 3),
    ('gapsize_check', 0), ('threads', None),
])
TGSTYPE = {'ont': 'ont', 'pb_clr': 'pb', 'pb_hifi': 'pb'}
PRESETS = ('auto', 'asm5', 'asm10', 'asm20')
CORRECTIONS = ('racon', 'none')
INT_RANGES = {'racon_rounds': (1, 4), 'min_nread': (1, 1000), 'max_candidate': (1, 100),
              'chunk': (1, 32), 'threads': (1, 128), 'min_match': (1, 100000)}
OBJ_NAME_RE = re.compile(r'^[A-Za-z0-9|._-]+$')
GAP_RE = re.compile(r'[Nn]+')
READS_TYPES = ('KBaseFile.SingleEndLibrary', 'KBaseAssembly.SingleEndLibrary')
ASSEMBLY_TYPE = 'KBaseGenomeAnnotations.Assembly'
LOW_COVERAGE = 5.0
LOW_CLOSURE = 0.25


# ---------------------------------------------------------------- parameters
def _blank(v):
    return v is None or (isinstance(v, str) and v.strip() == '')


def _to_int(key, v):
    try:
        f = float(str(v).strip())
    except ValueError:
        raise ValueError('Parameter "%s" must be a whole number, got %r' % (key, v))
    if f != int(f):
        raise ValueError('Parameter "%s" must be a whole number, got %r' % (key, v))
    return int(f)


def _to_bool(v):
    s = str(v).strip().lower()
    if s in ('true', 'yes', 'on'):
        return True
    try:
        return float(s) != 0
    except ValueError:
        return False


def validate_params(params):
    """Validate and coerce Narrative params ("" for untouched fields, numbers
    as strings such as "4.0"). Raises ValueError before any I/O."""
    if not isinstance(params, dict):
        raise ValueError('params must be a mapping')
    for key in REQUIRED:
        if _blank(params.get(key)):
            raise ValueError('Required parameter "%s" is missing or blank' % key)
    p = {k: str(params[k]).strip() for k in REQUIRED}
    name = p['output_assembly_name']
    if not OBJ_NAME_RE.match(name) or name.isdigit():
        raise ValueError('output_assembly_name %r is not a valid object name: use letters, '
                         'digits and . _ - | only, and not a bare number' % name)
    for key, dflt in DEFAULTS.items():
        v = params.get(key)
        p[key] = dflt if _blank(v) else v
    for key, allowed in (('read_type', tuple(TGSTYPE)), ('error_correction', CORRECTIONS),
                         ('minimap2_preset', PRESETS)):
        p[key] = str(p[key]).strip().lower()
        if p[key] not in allowed:
            raise ValueError('Parameter "%s" must be one of %s, got %r'
                             % (key, ', '.join(allowed), p[key]))
    for key, (lo, hi) in INT_RANGES.items():
        if p[key] is None:
            continue
        p[key] = _to_int(key, p[key])
        if not lo <= p[key] <= hi:
            raise ValueError('Parameter "%s" must be between %d and %d, got %d'
                             % (key, lo, hi, p[key]))
    if p['min_idy'] is not None:
        try:
            p['min_idy'] = float(str(p['min_idy']).strip())
        except ValueError:
            raise ValueError('Parameter "min_idy" must be a number, got %r' % p['min_idy'])
        if not 0 < p['min_idy'] <= 1:
            raise ValueError('Parameter "min_idy" must be > 0 and <= 1, got %s' % p['min_idy'])
    p['gapsize_check'] = _to_bool(p['gapsize_check'])
    return p


# ------------------------------------------------------------------- command
def build_command(scaffold_fa, reads_fa, p, threads, reads_are_contigs=False):
    """Return (cmd, warnings, notes). Only flags that differ from the tool's
    own defaults are passed, except --thread (the tool's default is 16)."""
    warnings, notes = [], []
    tool = os.path.join(TOOL_ENV_BIN, 'tgsgapcloser')
    racon = os.path.join(TOOL_ENV_BIN, 'racon')
    scaff_arg = os.path.abspath(scaffold_fa)      # the tool runs with cwd=job dir
    reads_arg = os.path.abspath(reads_fa)
    cmd = [tool, '--scaff', scaff_arg, '--reads', reads_arg, '--output', OUT_PREFIX,
           '--thread', str(threads)]
    tgstype = TGSTYPE[p['read_type']]
    if tgstype != 'ont':
        cmd += ['--tgstype', tgstype]
    preset = p['minimap2_preset']
    if preset == 'auto' and (p['read_type'] == 'pb_hifi' or reads_are_contigs):
        preset = 'asm20'
        notes.append('minimap2 preset set to asm20 automatically for %s, as the '
                     'TGS-GapCloser README suggests for accurate long sequences.'
                     % ('pre-assembled contigs' if reads_are_contigs else 'PacBio HiFi reads'))
    if preset != 'auto':
        cmd += ['--minmap_arg', ' -x %s' % preset]
    else:
        notes.append('minimap2 preset: tool default (-x %s).'
                     % ('ava-pb' if tgstype == 'pb' else 'ava-ont'))
    candidate_flags = (('racon_rounds', '--r_round'), ('chunk', '--chunk'),
                       ('min_nread', '--min_nread'), ('max_candidate', '--max_candidate'))
    if p['error_correction'] == 'none':
        cmd.append('--ne')
        for key, flag in candidate_flags:
            if p[key] != DEFAULTS[key]:
                warnings.append('%s=%s was not passed: with error correction off (--ne) the '
                                'tool skips candidate selection and correction, so %s has '
                                'no effect.' % (key, p[key], flag))
    else:
        cmd += ['--racon', racon]
        for key, flag in candidate_flags:
            if p[key] != DEFAULTS[key]:
                cmd += [flag, str(p[key])]
        if reads_are_contigs:
            notes.append('Racon correction of candidates is usually unnecessary when filling '
                         'from pre-assembled contigs; error correction = none is faster.')
    if p['min_idy'] is not None:
        cmd += ['--min_idy', '%g' % p['min_idy']]
    if p['min_match'] is not None:
        cmd += ['--min_match', str(p['min_match'])]
    if p['gapsize_check']:
        cmd.append('--g_check')
    return cmd, warnings, notes


def available_threads(requested=None):
    try:
        n = len(os.sched_getaffinity(0))      # respects container cpusets
    except AttributeError:
        n = os.cpu_count() or 4
    return max(1, min(n, requested)) if requested else max(1, n)


# ------------------------------------------------------------------- running
def run_cmd(cmd, log_path, cwd, extra_env=None):
    """Run with the tool env's bin first on PATH (subprocess only), stream to
    a log, and raise with the last lines on a non-zero exit."""
    env = dict(os.environ, **(extra_env or {}))
    env['PATH'] = TOOL_ENV_BIN + os.pathsep + env.get('PATH', '')
    tail = collections.deque(maxlen=60)
    with open(log_path, 'w') as lf:
        proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, universal_newlines=True, bufsize=1)
        for line in proc.stdout:
            lf.write(line)
            tail.append(line.rstrip('\n'))
        rc = proc.wait()
    if rc != 0:
        raise RuntimeError('%s failed (exit %d):\n%s%s' % (os.path.basename(cmd[0]), rc,
                           '\n'.join(tail), _newest_log_tail(cwd, exclude=log_path)))
    return rc


def _newest_log_tail(job_dir, exclude=None, n=15):
    logs = [os.path.join(job_dir, f) for f in os.listdir(job_dir)
            if (f.endswith('.log') or f.endswith('.err'))
            and os.path.join(job_dir, f) != exclude]
    if not logs:
        return ''
    newest = max(logs, key=os.path.getmtime)
    with open(newest, errors='replace') as fh:
        lines = fh.read().splitlines()[-n:]
    return '\n--- tail of %s ---\n%s' % (os.path.basename(newest), '\n'.join(lines))


def _first_existing(job_dir, names):
    for n in names:
        path = os.path.join(job_dir, n)
        if os.path.exists(path):
            return path
    return None


def tool_version(name, args=('--version',)):
    try:
        out = subprocess.run([os.path.join(TOOL_ENV_BIN, name)] + list(args),
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             universal_newlines=True, timeout=60).stdout
        return out.strip().splitlines()[0] if out.strip() else 'unknown'
    except Exception:  # noqa: BLE001 - informational only
        return 'unknown'


def run_tgsgapcloser(cmd, job_dir):
    os.makedirs(job_dir, exist_ok=True)
    stale = sorted(f for f in os.listdir(job_dir) if f.startswith('done_step'))
    if stale:
        raise RuntimeError('%s already holds %s: TGS-GapCloser would silently skip those '
                           'steps. Use a fresh directory.' % (job_dir, ', '.join(stale)))
    log_path = os.path.join(job_dir, 'tgsgapcloser.log')
    t0 = time.time()
    run_cmd(cmd, log_path, cwd=job_dir)
    with open(log_path, errors='replace') as fh:
        text = fh.read()
    bad = [ln for ln in text.splitlines() if re.match(r'^\s*(FATAL|Error)\s*:', ln)]
    if bad or 'ALL DONE' not in text:
        raise RuntimeError('TGS-GapCloser did not complete. (v1.2.1 exits with status 0 on '
                           'fatal errors, so its log is checked instead.)\n%s%s'
                           % ('\n'.join(bad or text.splitlines()[-20:]),
                              _newest_log_tail(job_dir, exclude=log_path)))
    out_fa = _first_existing(job_dir, [OUT_PREFIX + '.scaff_seqs', OUT_PREFIX + '.scaff_seq'])
    details = _first_existing(job_dir, [OUT_PREFIX + '.gap_fill_detail',
                                        OUT_PREFIX + '.gap_fill_details'])
    if not out_fa or os.path.getsize(out_fa) == 0 or not details:
        raise RuntimeError('TGS-GapCloser reported success but %s.scaff_seqs or '
                           '.gap_fill_detail is missing or empty in %s' % (OUT_PREFIX, job_dir))
    m = re.search(r'Version\s*:\s*([\w.]+)', text)
    return {'scaff_fasta': out_fa, 'details': details, 'log': log_path,
            'version': m.group(1) if m else 'unknown', 'seconds': round(time.time() - t0, 1)}


# ------------------------------------------------------------------- parsing
def _open_text(path):
    with open(path, 'rb') as fh:
        magic = fh.read(2)
    if magic == b'\x1f\x8b':
        return gzip.open(path, 'rt', errors='replace')
    return open(path, 'r', errors='replace')


def iter_fasta(path):
    name, chunks = None, []
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                fields = line[1:].split()
                name, chunks = (fields[0] if fields else ''), []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, ''.join(chunks)


def n50(lengths):
    total, acc = sum(lengths), 0
    for x in sorted(lengths, reverse=True):
        acc += x
        if acc * 2 >= total:
            return x
    return 0


def fasta_stats(path):
    """Scaffold and scaftig statistics. Internal N-runs are the gaps the tool
    can fill; N-runs touching a sequence end are counted separately."""
    ids, lengths, scaftigs, gaps = [], [], [], []
    terminal = 0
    for name, seq in iter_fasta(path):
        ids.append(name)
        lengths.append(len(seq))
        for m in GAP_RE.finditer(seq):
            if m.start() == 0 or m.end() == len(seq):
                terminal += 1
            else:
                gaps.append((name, m.start() + 1, m.end()))
        scaftigs.extend(len(x) for x in GAP_RE.split(seq) if x)
    gap_bp = sum(e - s + 1 for _, s, e in gaps)
    return {'ids': ids, 'n_seqs': len(ids), 'total_length': sum(lengths),
            'internal_gaps': len(gaps), 'gap_bp': gap_bp, 'gaps': gaps,
            'short_gaps': sum(1 for _, s, e in gaps if e - s + 1 < 10),
            'terminal_n_runs': terminal, 'scaffold_n50': n50(lengths),
            'n_scaftigs': len(scaftigs), 'scaftig_n50': n50(scaftigs),
            'largest_scaftig': max(scaftigs) if scaftigs else 0}


def reads_to_fasta(in_path, out_path):
    """Stream FASTQ or FASTA (plain or gzip) to plain FASTA: TGS-GapCloser
    v1.2.1 only accepts FASTA reads."""
    n, bp, lens = 0, 0, []
    with _open_text(in_path) as fh:
        first = fh.readline()
    if first.startswith('>'):
        fmt = 'fasta'
        with open(out_path, 'w') as out:
            for name, seq in iter_fasta(in_path):
                out.write('>%s\n%s\n' % (name or 'read%d' % (n + 1), seq))
                n, bp = n + 1, bp + len(seq)
                lens.append(len(seq))
    elif first.startswith('@'):
        fmt = 'fastq'
        with _open_text(in_path) as fh, open(out_path, 'w') as out:
            while True:
                header = fh.readline()
                if not header:
                    break
                if not header.strip():
                    continue
                seq, plus = fh.readline().strip(), fh.readline()
                fh.readline()
                if not header.startswith('@') or not plus.startswith('+'):
                    raise ValueError('Malformed FASTQ near record %d of the reads' % (n + 1))
                fields = header[1:].split()
                out.write('>%s\n%s\n' % (fields[0] if fields else 'read%d' % (n + 1), seq))
                n, bp = n + 1, bp + len(seq)
                lens.append(len(seq))
    else:
        raise ValueError('Reads file is neither FASTQ nor FASTA (starts with %r)' % first[:40])
    if n == 0:
        raise ValueError('The reads file contains no reads')
    return {'format': fmt, 'n_reads': n, 'total_bp': bp,
            'mean_length': round(bp / float(n), 1), 'read_n50': n50(lens)}


def parse_gap_details(path):
    """One event per junction between consecutive input segments (S rows) in
    <prefix>.gap_fill_detail; S rows are: out_start out_end S in_start in_end."""
    events, scaffold, segs = [], None, []

    def flush():
        s_idx = [i for i, sg in enumerate(segs) if sg['type'] == 'S' and 'in_start' in sg]
        for j, (a, b) in enumerate(zip(s_idx, s_idx[1:]), 1):
            prev, nxt, between = segs[a], segs[b], segs[a + 1:b]
            types = set(sg['type'] for sg in between)
            size = {t: sum(sg['out_end'] - sg['out_start'] + 1
                           for sg in between if sg['type'] == t) for t in ('F', 'N')}
            status = ('merged' if not between else 'filled' if types == {'F'}
                      else 'unfilled' if types == {'N'} else 'partial')
            events.append({
                'scaffold': scaffold, 'junction': j, 'status': status,
                'input_gap_start': prev['in_end'] + 1, 'input_gap_end': nxt['in_start'] - 1,
                'input_gap_len': nxt['in_start'] - prev['in_end'] - 1,
                'fill_bp': size['F'], 'remaining_n_bp': size['N'],
                'output_start': between[0]['out_start'] if between else prev['out_end'] + 1,
                'output_end': between[-1]['out_end'] if between else prev['out_end']})

    with open(path, errors='replace') as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if scaffold is not None:
                    flush()
                scaffold, segs = line[1:].split()[0], []
                continue
            f = line.split()
            if len(f) < 3:
                continue
            seg = {'out_start': int(f[0]), 'out_end': int(f[1]), 'type': f[2]}
            if f[2] == 'S' and len(f) >= 5:
                seg['in_start'], seg['in_end'] = int(f[3]), int(f[4])
            segs.append(seg)
    if scaffold is not None:
        flush()
    return events


def summarize(before, after, events, reads_info, tool, reads_are_contigs):
    """Headline numbers plus the warnings that the numbers justify."""
    warnings = []
    counts = collections.Counter(e['status'] for e in events)
    closed = counts['filled'] + counts['merged']
    n_gaps = before['internal_gaps']
    rate = closed / float(n_gaps) if n_gaps else 0.0
    nongap_bp = max(1, before['total_length'] - before['gap_bp'])
    coverage = reads_info['total_bp'] / float(nongap_bp)
    summary = collections.OrderedDict([
        ('gaps_before', n_gaps), ('gaps_after', after['internal_gaps']),
        ('gaps_closed', closed), ('filled', counts['filled']), ('merged', counts['merged']),
        ('partial', counts['partial']), ('unfilled', counts['unfilled']),
        ('closure_rate', round(rate, 4)), ('gap_bp_before', before['gap_bp']),
        ('gap_bp_after', after['gap_bp']), ('filled_bp', sum(e['fill_bp'] for e in events)),
        ('scaftig_n50_before', before['scaftig_n50']),
        ('scaftig_n50_after', after['scaftig_n50']),
        ('read_coverage', round(coverage, 1)), ('tool_version', tool['version']),
        ('runtime_s', tool['seconds'])])
    if len(events) != n_gaps:
        warnings.append('gap_fill_detail lists %d junctions but the input has %d internal '
                        'N-gaps; per-gap results may be incomplete.' % (len(events), n_gaps))
    if n_gaps - after['internal_gaps'] != closed:
        warnings.append('Closed gaps by table (%d) differ from the drop in N-runs (%d -> %d).'
                        % (closed, n_gaps, after['internal_gaps']))
    missing = sorted(set(before['ids']) - set(after['ids']))
    extra = sorted(set(after['ids']) - set(before['ids']))
    if missing or extra:
        warnings.append('Sequence IDs changed: %d missing from output (%s), %d new (%s).'
                        % (len(missing), ', '.join(missing[:5]), len(extra),
                           ', '.join(extra[:5])))
    if closed == 0:
        warnings.append('No gaps were closed. Common causes: low read coverage, reads from '
                        'a different strain, wrong read type, or min_idy/min_match too strict.')
    elif rate < LOW_CLOSURE:
        warnings.append('Only %.0f%% of gaps were closed.' % (100 * rate))
    if not reads_are_contigs and coverage < LOW_COVERAGE:
        warnings.append('Estimated read coverage is %.1fx; the TGS-GapCloser paper used '
                        'about 10x.' % coverage)
    if not reads_are_contigs and reads_info['mean_length'] < 1000:
        warnings.append('Mean read length is %.0f bp: reads this short rarely bridge gaps. '
                        'Are these long reads?' % reads_info['mean_length'])
    if tool['version'] != EXPECTED_TOOL_VERSION:
        warnings.append('TGS-GapCloser version %s differs from the %s this wrapper was '
                        'validated with.' % (tool['version'], EXPECTED_TOOL_VERSION))
    return summary, warnings


def tech_warnings(tech, read_type):
    t = (tech or '').lower()
    if 'illumina' in t or 'ion' in t:
        return ['The reads object says "%s": short reads cannot bridge gaps.' % tech]
    if 'pacbio' in t and read_type == 'ont':
        return ['The reads object says "%s" but read type is Nanopore.' % tech]
    if ('nanopore' in t or 'minion' in t) and read_type != 'ont':
        return ['The reads object says "%s" but read type is PacBio.' % tech]
    return []


# ------------------------------------------------------------------- reports
def _e(x):
    return html.escape(str(x))


def _table(header, rows):
    head = ''.join('<th>%s</th>' % _e(h) for h in header)
    body = ''.join('<tr>%s</tr>' % ''.join('<td>%s</td>' % _e(c) for c in r) for r in rows)
    return '<table><tr>%s</tr>%s</table>' % (head, body)


def _fmt(x):
    return '{:,}'.format(x) if isinstance(x, int) else x


def write_report_files(out_dir, ctx):
    """Write index.html, the summary and per-gap TSVs and a zip of the tool
    logs. ctx keys: before, after, events, summary, warnings, notes, params,
    commands, versions, inputs, job_dir. Returns a dict of paths."""
    html_dir = os.path.join(out_dir, 'html')
    os.makedirs(html_dir, exist_ok=True)
    b, a, s = ctx['before'], ctx['after'], ctx['summary']
    metrics = [('Sequences', 'n_seqs'), ('Total length (bp)', 'total_length'),
               ('Internal N-gaps', 'internal_gaps'), ('Gap bases (N)', 'gap_bp'),
               ('Scaffold N50', 'scaffold_n50'), ('Scaftigs (contigs)', 'n_scaftigs'),
               ('Scaftig N50', 'scaftig_n50'), ('Largest scaftig', 'largest_scaftig')]
    summary_tsv = os.path.join(out_dir, 'gapfill_summary.tsv')
    with open(summary_tsv, 'w') as fh:
        fh.write('metric\tbefore\tafter\n')
        for label, key in metrics:
            fh.write('%s\t%s\t%s\n' % (label, b[key], a[key]))
        for k, v in s.items():
            fh.write('%s\t\t%s\n' % (k, v))
    cols = ['scaffold', 'junction', 'status', 'input_gap_start', 'input_gap_end',
            'input_gap_len', 'fill_bp', 'remaining_n_bp', 'output_start', 'output_end']
    per_gap_tsv = os.path.join(out_dir, 'gapfill_per_gap.tsv')
    with open(per_gap_tsv, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for ev in ctx['events']:
            fh.write('\t'.join(str(ev[c]) for c in cols) + '\n')
    zip_path = os.path.join(out_dir, 'tgsgapcloser_output.zip')
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zf:
        for f in sorted(os.listdir(ctx['job_dir'])):
            if f.endswith(('.log', '.err', '.gap_fill_detail', '.name_map')):
                zf.write(os.path.join(ctx['job_dir'], f), 'tgsgapcloser/' + f)
        zf.write(summary_tsv, 'gapfill_summary.tsv')
        zf.write(per_gap_tsv, 'gapfill_per_gap.tsv')
    headline = 'Closed %d of %d gaps (%.1f%%): %d filled, %d merged, %d partial, %d unfilled' % (
        s['gaps_closed'], s['gaps_before'], 100 * s['closure_rate'], s['filled'],
        s['merged'], s['partial'], s['unfilled'])
    parts = ['<h2>TGS-GapCloser gap filling</h2><p class="head">%s</p>' % _e(headline)]
    if ctx['warnings']:
        parts.append('<div class="warn"><b>Warnings</b><ul>%s</ul></div>'
                     % ''.join('<li>%s</li>' % _e(w) for w in ctx['warnings']))
    if ctx['notes']:
        parts.append('<div class="note"><b>Notes</b><ul>%s</ul></div>'
                     % ''.join('<li>%s</li>' % _e(n) for n in ctx['notes']))
    parts.append('<h3>Assembly before and after</h3>' + _table(
        ['Metric', 'Before', 'After'],
        [(label, _fmt(b[key]), _fmt(a[key])) for label, key in metrics]))
    shown = ctx['events'][:200]
    parts.append('<h3>Per-gap results</h3><p>%s</p>' % _e(
        'Showing %d of %d junctions; the full table is in gapfill_per_gap.tsv.'
        % (len(shown), len(ctx['events']))) + _table(cols, [[ev[c] for c in cols]
                                                            for ev in shown]))
    parts.append('<h3>Inputs and parameters</h3>' + _table(
        ['Name', 'Value'], list(ctx['inputs'].items()) +
        [(k, 'tool default' if v is None else v) for k, v in ctx['params'].items()
         if k not in REQUIRED]))
    parts.append('<h3>Commands executed</h3><pre>%s</pre>'
                 % _e('\n\n'.join(ctx['commands'])))
    parts.append('<h3>Versions and citation</h3>' + _table(
        ['Tool', 'Version'], list(ctx['versions'].items())) +
        '<p>Xu M, Guo L, Gu S, et al. (2020) TGS-GapCloser: a fast and accurate gap closer '
        'for large genomes with low coverage of error-prone long reads. GigaScience 9(9): '
        'giaa094. doi:10.1093/gigascience/giaa094 (GPL-3.0). Also cite minimap2 (Li 2018, '
        'Bioinformatics) and, when used, Racon (Vaser et al. 2017, Genome Research).</p>')
    style = ('body{font-family:sans-serif;margin:16px;color:#222}table{border-collapse:'
             'collapse;margin:8px 0 16px}td,th{border:1px solid #ccc;padding:3px 8px;'
             'font-size:13px;text-align:left}th{background:#f2f2f2}.head{font-size:16px;'
             'font-weight:bold}.warn{background:#fff4e5;border:1px solid #f0b36b;padding:'
             '6px 12px;margin:8px 0}.note{background:#eef5fb;border:1px solid #9cc3e4;'
             'padding:6px 12px;margin:8px 0}pre{background:#f7f7f7;padding:8px;'
             'white-space:pre-wrap;font-size:12px}')
    with open(os.path.join(html_dir, 'index.html'), 'w') as fh:
        fh.write('<!DOCTYPE html><html><head><meta charset="utf-8"><title>TGS-GapCloser'
                 '</title><style>%s</style></head><body>%s</body></html>'
                 % (style, ''.join(parts)))
    return {'html_dir': html_dir, 'summary_tsv': summary_tsv,
            'per_gap_tsv': per_gap_tsv, 'zip': zip_path}


def shell_join(cmd):
    return ' '.join(shlex.quote(c) for c in cmd)


# ------------------------------------------------------------ KBase runner
class TGSGapCloserRunner(object):
    """Orchestrates KBase I/O around the pure functions above."""

    def __init__(self, config, callback_url, token):
        self.scratch = config['scratch']
        self.ws_url = config['workspace-url']
        self.callback_url = callback_url
        self.token = token
        self.results = None

    def _ws(self):
        from installed_clients.WorkspaceClient import Workspace
        return Workspace(self.ws_url, token=self.token)

    def _obj_type(self, ref):
        info = self._ws().get_object_info3({'objects': [{'ref': ref}]})['infos'][0]
        return info[2].split('-')[0], info[1]

    def _assembly_fasta(self, ref, dest):
        from installed_clients.AssemblyUtilClient import AssemblyUtil
        path = AssemblyUtil(self.callback_url).get_assembly_as_fasta({'ref': ref})['path']
        shutil.copy(path, dest)
        return dest

    def _reads_fasta(self, ref, obj_type, obj_name, job_dir, notes):
        out = os.path.join(job_dir, 'reads.fasta')
        if obj_type == ASSEMBLY_TYPE:
            src = self._assembly_fasta(ref, os.path.join(job_dir, 'source_contigs.fa'))
            info = reads_to_fasta(src, out)
            info['sequencing_tech'] = 'pre-assembled contigs'
            notes.append('Gap-fill source is the Assembly "%s" (%d sequences).'
                         % (obj_name, info['n_reads']))
            return out, info, True
        if obj_type not in READS_TYPES:
            raise ValueError('Object "%s" is a %s. Use single-end long reads (Nanopore or '
                             'PacBio) or an Assembly of contigs as the gap-fill source.'
                             % (obj_name, obj_type))
        from installed_clients.ReadsUtilsClient import ReadsUtils
        res = ReadsUtils(self.callback_url).download_reads({'read_libraries': [ref]})
        rec = list(res['files'].values())[0]
        files = rec['files']
        if files.get('type') not in (None, 'single') or files.get('rev'):
            raise ValueError('Reads "%s" are paired-end; TGS-GapCloser needs single-end long '
                             'reads.' % obj_name)
        info = reads_to_fasta(files['fwd'], out)
        info['sequencing_tech'] = rec.get('sequencing_tech') or 'unknown'
        notes.append('Reads "%s" (%s) converted %s -> FASTA: %s reads, %s bp, mean %s bp, '
                     'N50 %s bp.' % (obj_name, info['sequencing_tech'], info['format'].upper(),
                                     _fmt(info['n_reads']), _fmt(info['total_bp']),
                                     info['mean_length'], _fmt(info['read_n50'])))
        return out, info, False

    def run(self, params):
        p = validate_params(params)
        job_dir = os.path.join(self.scratch, 'kb_tgsgapcloser_' + uuid.uuid4().hex[:12])
        os.makedirs(job_dir)
        warnings, notes = [], []
        asm_type, asm_name = self._obj_type(p['assembly_ref'])
        if asm_type != ASSEMBLY_TYPE:
            raise ValueError('"%s" is a %s; choose an Assembly.' % (asm_name, asm_type))
        draft = self._assembly_fasta(p['assembly_ref'], os.path.join(job_dir, 'draft.fa'))
        before = fasta_stats(draft)
        if before['internal_gaps'] == 0:
            raise ValueError('Assembly "%s" has no internal N-gaps, so there is nothing to '
                             'fill. TGS-GapCloser fills gaps inside scaffolds: run a '
                             'scaffolder (e.g. RagTag scaffold) first.' % asm_name)
        if before['short_gaps']:
            notes.append('%d N-runs are shorter than 10 bp; TGS-GapCloser treats every N-run '
                         'as a gap.' % before['short_gaps'])
        if before['terminal_n_runs']:
            notes.append('%d N-runs touch a sequence end and are not counted as gaps.'
                         % before['terminal_n_runs'])
        rd_type, rd_name = self._obj_type(p['reads_ref'])
        reads_fa, reads_info, is_contigs = self._reads_fasta(
            p['reads_ref'], rd_type, rd_name, job_dir, notes)
        if not is_contigs:
            warnings += tech_warnings(reads_info['sequencing_tech'], p['read_type'])
        threads = available_threads(p['threads'])
        cmd, w, n = build_command(draft, reads_fa, p, threads, is_contigs)
        warnings += w
        notes += n + ['Threads used: %d.' % threads,
                      'TGS-GapCloser keeps scaftig order and orientation but not the '
                      'estimated sizes of gaps it fills.']
        tool = run_tgsgapcloser(cmd, os.path.join(job_dir, 'run'))
        after = fasta_stats(tool['scaff_fasta'])
        events = parse_gap_details(tool['details'])
        summary, w = summarize(before, after, events, reads_info, tool, is_contigs)
        warnings += w
        out_fa = os.path.join(job_dir, p['output_assembly_name'] + '.fa')
        shutil.copy(tool['scaff_fasta'], out_fa)
        from installed_clients.AssemblyUtilClient import AssemblyUtil
        new_ref = AssemblyUtil(self.callback_url).save_assembly_from_fasta({
            'file': {'path': out_fa}, 'workspace_name': p['workspace_name'],
            'assembly_name': p['output_assembly_name']})
        versions = collections.OrderedDict([
            ('tgsgapcloser', tool['version']), ('minimap2', tool_version('minimap2')),
            ('racon', tool_version('racon') if p['error_correction'] == 'racon'
             else 'not used')])
        files = write_report_files(os.path.join(job_dir, 'report'), {
            'before': before, 'after': after, 'events': events, 'summary': summary,
            'warnings': warnings, 'notes': notes, 'params': p, 'versions': versions,
            'commands': ['# reads -> FASTA conversion done in Python (FASTA-only tool)',
                         shell_join(cmd)], 'job_dir': os.path.dirname(tool['log']),
            'inputs': collections.OrderedDict([('Assembly', asm_name), ('Reads', rd_name),
                                               ('Output', p['output_assembly_name'])])})
        from installed_clients.KBaseReportClient import KBaseReport
        report = KBaseReport(self.callback_url).create_extended_report({
            'workspace_name': p['workspace_name'],
            'report_object_name': 'kb_tgsgapcloser_report_' + uuid.uuid4().hex[:8],
            'objects_created': [{'ref': new_ref, 'description': 'Gap-filled assembly'}],
            'message': '', 'direct_html_link_index': 0,
            'html_links': [{'path': files['html_dir'], 'name': 'index.html',
                            'description': 'TGS-GapCloser report'}],
            'file_links': [
                {'path': files['summary_tsv'], 'name': 'gapfill_summary.tsv',
                 'description': 'Before/after metrics'},
                {'path': files['per_gap_tsv'], 'name': 'gapfill_per_gap.tsv',
                 'description': 'Status of every gap'},
                {'path': files['zip'], 'name': 'tgsgapcloser_output.zip',
                 'description': 'Tool logs, gap_fill_detail and TSVs'}]})
        self.results = {'params': p, 'before': before, 'after': after, 'events': events,
                        'summary': summary, 'warnings': warnings, 'notes': notes,
                        'command': cmd, 'reads_info': reads_info, 'job_dir': job_dir}
        logging.info('kb_tgsgapcloser summary: %s', json.dumps(summary))
        return {'report_name': report['name'], 'report_ref': report['ref'],
                'assembly_ref': new_ref, 'gapfill_summary': dict(summary)}
