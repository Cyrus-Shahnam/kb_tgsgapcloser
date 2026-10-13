#!/usr/bin/env python3
"""Seeded known-truth fixture for kb_tgsgapcloser.

Writes into OUTDIR:
  truth.fa       gap-free truth: scaffolds s1..sN plus an ungapped control contig
  draft.fa       truth with planted gaps N-masked (coordinates unchanged)
  reads.fastq    ONT-like reads simulated from the TRUTH (substitutions, indels)
  manifest.json  seed, settings, and every planted gap with its bridging-read count

Every gap and read has its own seeded RNG, so a gap's position and a read's
sequence depend only on (seed, id). Stdlib only; deterministic.
Usage: python3 make_gap_fixture.py OUTDIR [--seed 20260929] [--coverage 25]
"""
import argparse
import json
import os
import random

BASES = 'ACGT'


def rand_seq(rng, n):
    return ''.join(rng.choice(BASES) for _ in range(n))


def mutate(rng, seq, sub, ins, dele):
    out = []
    for b in seq:
        r = rng.random()
        if r < dele:
            continue
        if r < dele + sub:
            out.append(rng.choice([x for x in BASES if x != b]))
        else:
            out.append(b)
        if rng.random() < ins:
            out.append(rng.choice(BASES))
    return ''.join(out)


def revcomp(s):
    return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def wrap(s, w=80):
    return '\n'.join(s[i:i + w] for i in range(0, len(s), w))


def make_fixture(outdir, seed=20260929, n_scaffolds=3, scaffold_len=60000,
                 control_len=20000, gap_lengths=(200, 500, 1000, 2000),
                 coverage=25.0, mean_read=8000, sub=0.03, ins=0.02, dele=0.02,
                 bridge_flank=1000):
    os.makedirs(outdir, exist_ok=True)
    truth = {}
    for i in range(1, n_scaffolds + 1):
        truth['s%d' % i] = rand_seq(random.Random('%s:truth:s%d' % (seed, i)), scaffold_len)
    truth['control'] = rand_seq(random.Random('%s:truth:control' % seed), control_len)

    gaps, draft = [], dict(truth)
    for i in range(1, n_scaffolds + 1):
        sid = 's%d' % i
        lengths = list(gap_lengths)
        random.Random('%s:gaplens:%s' % (seed, sid)).shuffle(lengths)
        slot = scaffold_len // (len(lengths) + 1)
        seq = list(draft[sid])
        for k, glen in enumerate(lengths, 1):
            rng = random.Random('%s:gap:%s:%d' % (seed, sid, k))
            centre = slot * k + rng.randint(-slot // 6, slot // 6)
            start = centre - glen // 2 + 1                     # 1-based inclusive
            end = start + glen - 1
            seq[start - 1:end] = ['N'] * glen
            gaps.append({'id': '%s_g%d' % (sid, k), 'scaffold': sid, 'start': start,
                         'end': end, 'length': glen, 'bridging_reads': 0})
        draft[sid] = ''.join(seq)

    total = sum(len(s) for s in truth.values())
    n_reads = int(coverage * total / mean_read)
    names = sorted(truth)
    weights = [len(truth[n]) for n in names]
    with open(os.path.join(outdir, 'reads.fastq'), 'w') as fq:
        for r in range(n_reads):
            rng = random.Random('%s:read:%d' % (seed, r))
            sid = rng.choices(names, weights=weights)[0]
            src = truth[sid]
            length = min(len(src), max(1500, int(rng.lognormvariate(0, 0.45) * mean_read)))
            start = rng.randint(0, len(src) - length)
            frag = src[start:start + length]
            for g in gaps:
                if (g['scaffold'] == sid and start < g['start'] - bridge_flank
                        and start + length > g['end'] + bridge_flank):
                    g['bridging_reads'] += 1
            read = mutate(rng, frag if rng.random() < 0.5 else revcomp(frag), sub, ins, dele)
            fq.write('@r%d %s:%d-%d\n%s\n+\n%s\n' % (r, sid, start + 1, start + length,
                                                  read, '5' * len(read)))
    for name, seqs in (('truth.fa', truth), ('draft.fa', draft)):
        with open(os.path.join(outdir, name), 'w') as fh:
            for sid in names:
                fh.write('>%s\n%s\n' % (sid, wrap(seqs[sid])))
    manifest = {'seed': seed, 'coverage': coverage, 'mean_read': mean_read,
                'error_rates': {'sub': sub, 'ins': ins, 'del': dele},
                'bridge_flank': bridge_flank, 'n_reads': n_reads, 'total_truth_bp': total,
                'control_ids': ['control'], 'gaps': gaps}
    with open(os.path.join(outdir, 'manifest.json'), 'w') as fh:
        json.dump(manifest, fh, indent=1)
    return manifest


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('outdir')
    ap.add_argument('--seed', type=int, default=20260929)
    ap.add_argument('--coverage', type=float, default=25.0)
    a = ap.parse_args()
    m = make_fixture(a.outdir, seed=a.seed, coverage=a.coverage)
    print('planted %d gaps (%d bp); %d reads; bridgeable gaps: %d' % (
        len(m['gaps']), sum(g['length'] for g in m['gaps']), m['n_reads'],
        sum(1 for g in m['gaps'] if g['bridging_reads'] > 0)))
