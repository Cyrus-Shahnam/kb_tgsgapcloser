#!/usr/bin/env python3
"""Dry-run the wrapper's orchestration without KBase services.

Runs: fixture -> reads FASTQ->FASTA -> command -> tool -> parsing -> report
files, then scores planted gaps against the manifest. Nothing touches KBase.

  # before the first Docker build, with the fake tool:
  python3 tools/dry_run_local.py --fake
  # inside the built image, with the real tool:
  docker run --rm --entrypoint bash -v $PWD:/kb/module <image> \
      -c "cd /kb/module && python3 tools/dry_run_local.py"
Add --fatal (with --fake) to check that exit-0 fatal errors are caught.
"""
import argparse
import json
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, 'lib'))
sys.path.insert(0, os.path.join(ROOT, 'test', 'data'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fake', action='store_true', help='use tools/fake_tgsgapcloser.py')
    ap.add_argument('--fatal', action='store_true', help='fake tool prints FATAL, exits 0')
    ap.add_argument('--ne', action='store_true', help='error_correction = none')
    ap.add_argument('--out', default=None, help='keep outputs here (default: temp dir)')
    a = ap.parse_args()
    work = a.out or tempfile.mkdtemp(prefix='tgs_dry_')
    os.makedirs(work, exist_ok=True)
    if a.fake:
        fakebin = os.path.join(work, 'fakebin')
        os.makedirs(fakebin, exist_ok=True)
        for name in ('tgsgapcloser', 'racon', 'minimap2'):
            dst = os.path.join(fakebin, name)
            if not os.path.exists(dst):
                shutil.copy(os.path.join(HERE, 'fake_tgsgapcloser.py'), dst)
                os.chmod(dst, 0o755)
        os.environ['TGS_ENV_BIN'] = fakebin
        if a.fatal:
            os.environ['FAKE_TGS_FATAL'] = '1'
    from kb_tgsgapcloser import tgsgapcloser_utils as U   # after TGS_ENV_BIN is set
    from make_gap_fixture import make_fixture

    fx = os.path.join(work, 'fixture')
    manifest = make_fixture(fx)
    p = U.validate_params({'workspace_name': 'dry', 'assembly_ref': '1/1/1',
                           'reads_ref': '1/2/1', 'output_assembly_name': 'dry_out',
                           'error_correction': 'none' if a.ne else '',
                           'racon_rounds': '', 'min_idy': '', 'threads': '4.0'})
    before = U.fasta_stats(os.path.join(fx, 'draft.fa'))
    reads_fa = os.path.join(work, 'reads.fasta')
    info = U.reads_to_fasta(os.path.join(fx, 'reads.fastq'), reads_fa)
    cmd, warnings, notes = U.build_command(os.path.join(fx, 'draft.fa'), reads_fa, p,
                                           U.available_threads(p['threads']))
    print('command:', U.shell_join(cmd))
    run_dir = os.path.join(work, 'run')
    if os.path.exists(run_dir):
        shutil.rmtree(run_dir)
    try:
        tool = U.run_tgsgapcloser(cmd, run_dir)
    except RuntimeError as err:
        print('RUN FAILED (expected with --fatal):\n%s' % err)
        return 0 if a.fatal else 1
    if a.fatal:
        print('ERROR: --fatal run was NOT detected as a failure')
        return 1
    after = U.fasta_stats(tool['scaff_fasta'])
    events = U.parse_gap_details(tool['details'])
    summary, w = U.summarize(before, after, events, info, tool, False)
    files = U.write_report_files(os.path.join(work, 'report'), {
        'before': before, 'after': after, 'events': events, 'summary': summary,
        'warnings': warnings + w, 'notes': notes, 'params': p, 'versions': {'tgs': 'dry'},
        'commands': [U.shell_join(cmd)], 'job_dir': run_dir,
        'inputs': {'Assembly': 'draft.fa', 'Reads': 'reads.fastq'}})
    print(json.dumps(summary, indent=1))
    ok = 0
    for g in manifest['gaps']:
        hit = [e for e in events if e['scaffold'] == g['scaffold']
               and e['input_gap_start'] <= g['end'] and e['input_gap_end'] >= g['start']]
        status = hit[0]['status'] if len(hit) == 1 else 'NOT FOUND (%d matches)' % len(hit)
        ok += status in ('filled', 'merged')
        print('%-7s %5d bp  bridging=%-3d -> %s' % (g['id'], g['length'],
                                                    g['bridging_reads'], status))
    ctrl = dict(U.iter_fasta(tool['scaff_fasta'])).get('control')
    truth = dict(U.iter_fasta(os.path.join(fx, 'truth.fa')))['control']
    print('closed %d/%d planted gaps; control contig unchanged: %s'
          % (ok, len(manifest['gaps']), ctrl == truth))
    print('warnings:', warnings + w)
    print('report:', os.path.join(files['html_dir'], 'index.html'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
