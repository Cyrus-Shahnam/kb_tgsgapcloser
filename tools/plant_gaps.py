#!/usr/bin/env python3
"""Plant seeded N-gaps into a real reference genome (real-data known truth).

Two modes:
  default        mask each chosen region with N of the same length, so the
                 draft keeps the reference coordinates;
  --fixed-n 100  replace each region with exactly 100 N, like a scaffolder's
                 gap of unknown size, so the tool must fill more sequence than
                 the gap's N-count suggests. Draft coordinates shift; the
                 manifest records both draft and reference coordinates.

  python3 tools/plant_gaps.py GCF_000007445.1_ASM744v1_genomic.fna CFT073_gapped \
      --n-gaps 50 --seed 7 [--fixed-n 100]
writes CFT073_gapped.fa and CFT073_gapped.manifest.json.
"""
import argparse
import json
import random


def read_fasta(path):
    seqs, name = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                name = line[1:].split()[0]
                seqs.append([name, []])
            elif name:
                seqs[-1][1].append(line.strip().upper())
    return [(n, ''.join(c)) for n, c in seqs]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('reference')
    ap.add_argument('out_prefix')
    ap.add_argument('--n-gaps', type=int, default=50)
    ap.add_argument('--lengths', default='200,500,1000,2000,5000',
                    help='comma-separated gap lengths, cycled in seeded order')
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--min-spacing', type=int, default=20000,
                    help='minimum bp between gaps and from sequence ends')
    ap.add_argument('--fixed-n', type=int, default=0,
                    help='replace each region with exactly this many N (0 = same length)')
    a = ap.parse_args()
    rng = random.Random(a.seed)
    lengths = [int(x) for x in a.lengths.split(',')]
    seqs = read_fasta(a.reference)
    usable = [(n, s) for n, s in seqs if len(s) > 3 * a.min_spacing]
    if not usable:
        raise SystemExit('no sequence is longer than 3 x --min-spacing')
    chosen, tries = [], 0
    while len(chosen) < a.n_gaps and tries < 100000:
        tries += 1
        name, seq = rng.choices(usable, weights=[len(s) for _, s in usable])[0]
        glen = lengths[len(chosen) % len(lengths)]
        start = rng.randint(a.min_spacing, len(seq) - a.min_spacing - glen)
        if 'N' in seq[start:start + glen]:
            continue
        if all(c['scaffold'] != name or abs(c['ref_start'] - start) > a.min_spacing + glen
               for c in chosen):
            chosen.append({'scaffold': name, 'ref_start': start + 1, 'ref_end': start + glen,
                           'length': glen})
    if len(chosen) < a.n_gaps:
        raise SystemExit('placed only %d gaps; lower --n-gaps or --min-spacing' % len(chosen))
    out_seqs = []
    for name, seq in seqs:
        mine = sorted((c for c in chosen if c['scaffold'] == name), key=lambda c: c['ref_start'])
        parts, pos, shift = [], 0, 0
        for c in mine:
            parts.append(seq[pos:c['ref_start'] - 1])
            n_count = a.fixed_n or c['length']
            c['draft_start'] = c['ref_start'] + shift
            c['draft_end'] = c['draft_start'] + n_count - 1
            parts.append('N' * n_count)
            shift += n_count - c['length']
            pos = c['ref_end']
        parts.append(seq[pos:])
        out_seqs.append((name, ''.join(parts)))
    chosen.sort(key=lambda c: (c['scaffold'], c['ref_start']))
    for i, c in enumerate(chosen, 1):
        c['id'] = 'gap%03d' % i
    with open(a.out_prefix + '.fa', 'w') as fh:
        for name, seq in out_seqs:
            fh.write('>%s\n' % name)
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + '\n')
    manifest = {'reference': a.reference, 'seed': a.seed, 'n_gaps': len(chosen),
                'lengths': lengths, 'fixed_n': a.fixed_n, 'min_spacing': a.min_spacing,
                'planted_bp': sum(c['length'] for c in chosen), 'gaps': chosen}
    with open(a.out_prefix + '.manifest.json', 'w') as fh:
        json.dump(manifest, fh, indent=1)
    print('planted %d gaps (%d bp) in %d sequences -> %s.fa'
          % (len(chosen), manifest['planted_bp'], len(set(c['scaffold'] for c in chosen)),
             a.out_prefix))


if __name__ == '__main__':
    main()
