#!/usr/bin/env python3
"""Fake TGS-GapCloser v1.2.1 for dry runs (see tools/dry_run_local.py).

Mimics the real driver's CLI, cwd tag files, output names and formats, and
its quirk of exiting 0 after FATAL. Fills every internal N-run of >= 300 bp
with placeholder sequence; shorter runs stay as N. Installed under the names
tgsgapcloser, racon and minimap2 in a temp bin dir; only tgsgapcloser acts.
"""
import os
import re
import sys


def fatal(msg):
    print('\nFATAL :   %s ! exit ...' % msg)
    sys.exit(0)                                  # the real v1.2.1 exits 0 here


def main(argv):
    name = os.path.basename(sys.argv[0])
    if name != 'tgsgapcloser':
        print('%s fake 0.0' % name)
        return 0
    print('INFO  :   Run tgsgapcloser. from fake ;\n          Version : 1.2.1 ;')
    opts, i = {}, 0
    while i < len(argv):
        if argv[i] in ('--ne', '--g_check'):
            opts[argv[i]] = True
            i += 1
        else:
            opts[argv[i]] = argv[i + 1] if i + 1 < len(argv) else ''
            i += 2
    if os.environ.get('FAKE_TGS_FATAL'):
        fatal('File %s does not exist !!!' % opts.get('--scaff'))
    for key in ('--scaff', '--reads', '--output'):
        if not opts.get(key):
            fatal('Arg %s is not assigned !!!' % key)
    for key in ('--scaff', '--reads'):
        if not os.path.exists(opts[key]):
            fatal('File %s does not exist !!!' % opts[key])
    with open(opts['--reads']) as fh:
        if not fh.read(1) == '>':
            print('Segmentation fault (reads are not FASTA)')
            return 139
    prefix = opts['--output']
    seqs, name_ = [], None
    with open(opts['--scaff']) as fh:
        for line in fh:
            if line.startswith('>'):
                seqs.append([line[1:].split()[0], []])
            elif seqs:
                seqs[-1][1].append(line.strip())
    out_fa = open(prefix + '.scaff_seqs', 'w')
    det = open(prefix + '.gap_fill_detail', 'w')
    nmap = open(prefix + '.name_map', 'w')
    for k, (name_, chunks) in enumerate(seqs, 1):
        seq = ''.join(chunks)
        nmap.write('%d\t%s\n' % (k, name_))
        det.write('>%s\n' % name_)
        out, pos, opos = [], 0, 0
        runs = [m for m in re.finditer(r'[Nn]+', seq) if 0 < m.start() and m.end() < len(seq)]
        for m in runs + [None]:
            end = m.start() if m else len(seq)
            piece = seq[pos:end]
            det.write('%d\t%d\tS\t%d\t%d\t\n' % (opos + 1, opos + len(piece), pos + 1, end))
            out.append(piece)
            opos += len(piece)
            if m is None:
                break
            glen = m.end() - m.start()
            kind, fill = ('F', ('ACGT' * glen)[:glen - 5]) if glen >= 300 else ('N', 'N' * glen)
            det.write('%d\t%d\t%s\n' % (opos + 1, opos + len(fill), kind))
            out.append(fill)
            opos += len(fill)
            pos = m.end()
        out_fa.write('>%s\n%s\n' % (name_, ''.join(out)))
    for fh in (out_fa, det, nmap):
        fh.close()
    for tag in ('done_step1_tag', 'done_step_4.2_tag'):
        open(tag, 'w').close()
    print('\nINFO  :   ALL DONE !!! ')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
