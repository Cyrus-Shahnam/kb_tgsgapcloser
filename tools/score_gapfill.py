#!/usr/bin/env python3
"""Score a gap-filling run against a plant_gaps.py manifest.

Status of every planted gap comes from the report's gapfill_per_gap.tsv
(matched by draft coordinates). With --reference and --filled, each closed gap
is also checked for accuracy and placement: the fill plus FLANK bp either side
is cut from the filled assembly and aligned to the reference with minimap2.
The flanks are reference sequence, so the alignment's NM is the fill's edits.

  python3 tools/score_gapfill.py CFT073_gapped.manifest.json gapfill_per_gap.tsv \
      --reference GCF_000007445.1_ASM744v1_genomic.fna --filled filled.fa \
      --out score.tsv
minimap2 is taken from --minimap2, else $TGS_ENV_BIN/minimap2, else PATH.
"""
import argparse
import csv
import json
import os
import statistics
import subprocess
import sys
import tempfile


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.strip())
    return {n: ''.join(c) for n, c in seqs.items()}


def align(minimap2, reference, frags):
    """Best (longest) alignment per fragment: (target, tstart, tend, NM, block)."""
    with tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False) as fh:
        for fid, seq in frags.items():
            fh.write('>%s\n%s\n' % (fid, seq))
        path = fh.name
    try:
        paf = subprocess.run([minimap2, '-c', '-x', 'asm20', '--secondary=no', reference,
                              path], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                             universal_newlines=True, check=True).stdout
    finally:
        os.unlink(path)
    best = {}
    for line in paf.splitlines():
        f = line.split('\t')
        nm = sum(int(t[5:]) for t in f[12:] if t.startswith('NM:i:'))
        hit = (f[5], int(f[7]), int(f[8]), nm, int(f[10]))
        if f[0] not in best or hit[4] > best[f[0]][4]:
            best[f[0]] = hit
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('manifest')
    ap.add_argument('per_gap_tsv')
    ap.add_argument('--reference')
    ap.add_argument('--filled')
    ap.add_argument('--flank', type=int, default=1000)
    ap.add_argument('--minimap2')
    ap.add_argument('--out')
    a = ap.parse_args()
    manifest = json.load(open(a.manifest))
    with open(a.per_gap_tsv) as fh:
        events = list(csv.DictReader(fh, delimiter='\t'))
    rows = []
    for g in manifest['gaps']:
        hit = [e for e in events if e['scaffold'] == g['scaffold']
               and int(e['input_gap_start']) <= g['draft_end']
               and int(e['input_gap_end']) >= g['draft_start']]
        row = dict(g, status=hit[0]['status'] if len(hit) == 1 else 'unmatched(%d)' % len(hit),
                   fill_bp=int(hit[0]['fill_bp']) if len(hit) == 1 else '',
                   edits='', identity='', placed='')
        row['_event'] = hit[0] if len(hit) == 1 else None
        rows.append(row)
    closed = [r for r in rows if r['status'] in ('filled', 'merged')]
    if a.reference and a.filled:
        mm2 = a.minimap2 or os.path.join(os.environ.get('TGS_ENV_BIN', ''), 'minimap2')
        if not os.path.exists(mm2):
            mm2 = 'minimap2'
        filled = read_fasta(a.filled)
        frags = {}
        for r in closed:
            e, seq = r['_event'], filled.get(r['scaffold'], '')
            lo = max(0, int(e['output_start']) - 1 - a.flank)
            frags[r['id']] = seq[lo:int(e['output_end']) + a.flank]
        best = align(mm2, a.reference, frags)
        for r in closed:
            h = best.get(r['id'])
            if not h:
                r['placed'] = 'no alignment'
                continue
            tol = 100 + abs(r['fill_bp'] - r['length'])
            r['placed'] = 'yes' if (h[0] == r['scaffold']
                                    and abs(h[1] - (r['ref_start'] - 1 - a.flank)) <= tol
                                    and abs(h[2] - (r['ref_end'] + a.flank)) <= tol) else 'NO'
            r['edits'] = h[3]
            r['identity'] = round(1 - h[3] / float(max(r['length'], r['fill_bp'], 1)), 4)
    cols = ['id', 'scaffold', 'ref_start', 'ref_end', 'length', 'draft_start', 'draft_end',
            'status', 'fill_bp', 'edits', 'identity', 'placed']
    out = open(a.out, 'w') if a.out else sys.stdout
    out.write('\t'.join(cols) + '\n')
    for r in rows:
        out.write('\t'.join(str(r[c]) for c in cols) + '\n')
    if a.out:
        out.close()
    n = len(rows)
    print('\nplanted gaps           %d (%d bp)' % (n, sum(r['length'] for r in rows)),
          file=sys.stderr)
    print('closed                 %d (%.1f%%)' % (len(closed), 100.0 * len(closed) / n),
          file=sys.stderr)
    print('status counts          %s' % {s: sum(1 for r in rows if r['status'] == s)
                                        for s in sorted(set(r['status'] for r in rows))},
          file=sys.stderr)
    scored = [r for r in closed if r['edits'] != '']
    if scored:
        edits = sum(r['edits'] for r in scored)
        bp = sum(r['length'] for r in scored)
        print('correctly placed       %d / %d closed' % (sum(r['placed'] == 'yes'
                                                          for r in closed), len(closed)),
              file=sys.stderr)
        print('fill edits             %d over %d bp = %.1f per 100 kb (identity %.4f)'
              % (edits, bp, 1e5 * edits / bp, 1 - edits / float(bp)), file=sys.stderr)
        print('median gap identity    %.4f' % statistics.median(r['identity'] for r in scored),
              file=sys.stderr)
        print('mean |fill - true len| %.1f bp' % statistics.mean(
            abs(r['fill_bp'] - r['length']) for r in scored), file=sys.stderr)


if __name__ == '__main__':
    main()
