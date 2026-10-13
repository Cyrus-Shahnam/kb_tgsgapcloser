# -*- coding: utf-8 -*-
"""kb_tgsgapcloser tests.

Fast tier (test_f*): validation, command building, parsers. No KBase I/O.
Known-truth tier (test_kt*): a seeded fixture (test/data/make_gap_fixture.py)
plants 12 N-gaps in 3 random scaffolds and simulates ~7%-error ONT-like reads
from the gap-free truth. The tool runs once per mode, lazily, and several tests
assert on the cached result. Local reference run (1 CPU, v1.2.1): racon mode
closed 12/12, 13 edits over 11,100 planted bp; --ne mode 780 edits; truth
contigs as the source rebuilt the truth exactly.
"""
import gzip
import os
import subprocess
import sys
import tempfile
import time
import unittest
from configparser import ConfigParser

from kb_tgsgapcloser import tgsgapcloser_utils as U
from kb_tgsgapcloser.kb_tgsgapcloserImpl import kb_tgsgapcloser
from kb_tgsgapcloser.kb_tgsgapcloserServer import MethodContext
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.ReadsUtilsClient import ReadsUtils
from installed_clients.WorkspaceClient import Workspace

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
sys.path.insert(0, DATA_DIR)
from make_gap_fixture import make_fixture  # noqa: E402

# Loose first-run thresholds; tighten after the first real run (see README).
MIN_CLOSURE_OF_BRIDGEABLE = 0.75     # observed 1.00
MAX_EDITS_PER_PLANTED_BP = 0.02      # observed 0.0012 (13 / 11,100)
MIN_ALIGNED_FRACTION = 0.98          # observed 1.00
MIN_BRIDGEABLE_GAPS = 6              # fewer -> skip, the fixture can't test anything
BASE = {'workspace_name': 'w', 'assembly_ref': '1/2/3', 'reads_ref': '1/3/1',
        'output_assembly_name': 'out'}

README_DETAIL = """>scaffold_1
1\t1000\tS\t1000\t2000
1001\t1010\tN
1011\t1100\tS\t2201\t2290
1101\t1110\tF
1111\t1200\tS\t2301\t2390
1201\t1300\tS\t2350\t2449
>scaffold_2
1\t500\tS\t1\t500
"""


def _write(path, text, gz=False):
    with (gzip.open(path, 'wt') if gz else open(path, 'w')) as fh:
        fh.write(text)
    return path


class kb_tgsgapcloserTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.token = os.environ.get('KB_AUTH_TOKEN', None)
        config = ConfigParser()
        config.read(os.environ.get('KB_DEPLOYMENT_CONFIG'))
        cls.cfg = dict(config.items('kb_tgsgapcloser'))
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': cls.token, 'user_id': 'test', 'authenticated': 1,
                        'provenance': [{'service': 'kb_tgsgapcloser',
                                        'method': 'please_never_use_it_in_production',
                                        'method_params': []}]})
        cls.ws = Workspace(cls.cfg['workspace-url'], token=cls.token)
        cls.impl = kb_tgsgapcloser(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        cls.wsName = 'test_kb_tgsgapcloser_%d' % int(time.time() * 1000)
        cls.ws.create_workspace({'workspace': cls.wsName})
        cls._staged = None
        cls._runs = {}

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.ws.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    # ------------------------------------------------------------ helpers
    @classmethod
    def stage(cls):
        """Generate the fixture in scratch (the utility containers mount
        scratch, not test/) and upload it once."""
        if cls._staged:
            return cls._staged
        fx = os.path.join(cls.scratch, 'fixture_%d' % int(time.time() * 1000))
        manifest = make_fixture(fx)
        au = AssemblyUtil(cls.callback_url)
        refs = {}
        for name, fname in (('draft', 'draft.fa'), ('truth', 'truth.fa')):
            refs[name] = au.save_assembly_from_fasta({
                'file': {'path': os.path.join(fx, fname)},
                'workspace_name': cls.wsName, 'assembly_name': 'kt_%s' % name})
        refs['reads'] = ReadsUtils(cls.callback_url).upload_reads({
            'fwd_file': os.path.join(fx, 'reads.fastq'), 'sequencing_tech': 'Nanopore',
            'single_genome': 1, 'wsname': cls.wsName, 'name': 'kt_sim_ont_reads'})['obj_ref']
        cls._staged = {'dir': fx, 'manifest': manifest, 'refs': refs}
        return cls._staged

    @classmethod
    def run_mode(cls, mode):
        """Run once per mode and cache. 'racon' goes through the Runner to keep
        its detailed results; 'ne' and 'contigs' go through the Impl."""
        if mode in cls._runs:
            return cls._runs[mode]
        st = cls.stage()
        params = dict(workspace_name=cls.wsName, assembly_ref=st['refs']['draft'],
                      reads_ref=st['refs']['reads'], output_assembly_name='kt_' + mode)
        if mode == 'racon':
            runner = U.TGSGapCloserRunner(cls.cfg, cls.callback_url, cls.token)
            out = runner.run(params)
            res = runner.results
        else:
            params['error_correction'] = 'none'
            if mode == 'contigs':
                params['reads_ref'] = st['refs']['truth']
            out = cls.impl.run_tgsgapcloser(cls.ctx, params)[0]
            res = None
        fa = AssemblyUtil(cls.callback_url).get_assembly_as_fasta(
            {'ref': out['assembly_ref']})['path']
        cls._runs[mode] = {'out': out, 'results': res, 'fasta': fa}
        return cls._runs[mode]

    def edit_distance(self, query_fa):
        """(total NM, aligned fraction) of query vs truth via minimap2 NM tags."""
        truth = os.path.join(self.stage()['dir'], 'truth.fa')
        paf = subprocess.run([os.path.join(U.TOOL_ENV_BIN, 'minimap2'), '-c', '-x', 'asm20',
                              '--secondary=no', truth, query_fa], stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, universal_newlines=True,
                             check=True).stdout
        nm, aligned, lengths = 0, 0, {}
        for line in paf.splitlines():
            f = line.split('\t')
            lengths[f[0]] = int(f[1])
            aligned += int(f[3]) - int(f[2])
            nm += sum(int(t[5:]) for t in f[12:] if t.startswith('NM:i:'))
        return nm, aligned / float(max(1, sum(lengths.values())))

    def planted_bp(self):
        return sum(g['length'] for g in self.stage()['manifest']['gaps'])

    # ---------------------------------------------------------- fast tier
    def test_f01_missing_required(self):
        for key in U.REQUIRED:
            for bad in (None, '', '   '):
                p = dict(BASE)
                p[key] = bad
                with self.assertRaisesRegex(ValueError, key):
                    U.validate_params(p)

    def test_f02_blank_and_stringified_params(self):
        p = U.validate_params(dict(BASE, read_type='', racon_rounds='2.0', min_idy='',
                                   min_match='500', threads='4', gapsize_check='1',
                                   max_candidate=' '))
        self.assertEqual(p['read_type'], 'ont')
        self.assertEqual(p['racon_rounds'], 2)
        self.assertIsNone(p['min_idy'])
        self.assertEqual(p['min_match'], 500)
        self.assertEqual(p['max_candidate'], 10)
        self.assertTrue(p['gapsize_check'])
        self.assertFalse(U.validate_params(dict(BASE, gapsize_check=0))['gapsize_check'])

    def test_f03_bad_values(self):
        for key, val in (('read_type', 'illumina'), ('racon_rounds', '9'),
                         ('min_nread', '0'), ('min_idy', '1.5'), ('threads', '2.5'),
                         ('error_correction', 'pilon'), ('minimap2_preset', 'sr'),
                         ('output_assembly_name', '123'),
                         ('output_assembly_name', 'has space')):
            with self.assertRaises(ValueError, msg='%s=%s' % (key, val)):
                U.validate_params(dict(BASE, **{key: val}))

    def test_f04_command_racon_defaults(self):
        p = U.validate_params(dict(BASE))
        cmd, warnings, _ = U.build_command('d.fa', 'r.fa', p, 8)
        self.assertEqual(cmd[cmd.index('--racon') + 1], os.path.join(U.TOOL_ENV_BIN, 'racon'))
        self.assertTrue(os.path.isabs(cmd[cmd.index('--scaff') + 1]))
        self.assertTrue(os.path.isabs(cmd[cmd.index('--reads') + 1]))
        self.assertEqual(cmd[cmd.index('--thread') + 1], '8')
        for flag in ('--ne', '--tgstype', '--r_round', '--minmap_arg', '--min_nread'):
            self.assertNotIn(flag, cmd)
        self.assertEqual(warnings, [])

    def test_f05_ne_drops_candidate_flags(self):
        p = U.validate_params(dict(BASE, error_correction='none', min_nread=2,
                                   max_candidate=5, racon_rounds=2))
        cmd, warnings, _ = U.build_command('d.fa', 'r.fa', p, 2)
        self.assertIn('--ne', cmd)
        for flag in ('--racon', '--min_nread', '--max_candidate', '--r_round'):
            self.assertNotIn(flag, cmd)
        self.assertEqual(len(warnings), 3)

    def test_f06_hifi_and_contigs_presets(self):
        p = U.validate_params(dict(BASE, read_type='pb_hifi'))
        cmd, _, notes = U.build_command('d.fa', 'r.fa', p, 2)
        self.assertEqual(cmd[cmd.index('--tgstype') + 1], 'pb')
        self.assertEqual(cmd[cmd.index('--minmap_arg') + 1], ' -x asm20')
        cmd, _, notes = U.build_command('d.fa', 'r.fa', U.validate_params(dict(BASE)), 2,
                                        reads_are_contigs=True)
        self.assertEqual(cmd[cmd.index('--minmap_arg') + 1], ' -x asm20')
        self.assertTrue(any('contigs' in n for n in notes))

    def test_f07_parse_gap_details(self):
        d = tempfile.mkdtemp(dir=self.scratch)
        ev = U.parse_gap_details(_write(os.path.join(d, 'x.gap_fill_detail'), README_DETAIL))
        self.assertEqual([e['status'] for e in ev], ['unfilled', 'filled', 'merged'])
        self.assertEqual((ev[0]['input_gap_start'], ev[0]['input_gap_end']), (2001, 2200))
        self.assertEqual(ev[1]['fill_bp'], 10)
        self.assertEqual(ev[0]['remaining_n_bp'], 10)

    def test_f08_stats_and_reads_conversion(self):
        d = tempfile.mkdtemp(dir=self.scratch)
        st = U.fasta_stats(_write(os.path.join(d, 'a.fa'),
                                  '>a desc\nNNACGTNNNNACGT\nACNGT\n>b\nACGT\n'))
        self.assertEqual((st['ids'], st['internal_gaps'], st['terminal_n_runs']),
                         (['a', 'b'], 2, 1))
        self.assertEqual(st['gap_bp'], 5)
        fq = _write(os.path.join(d, 'r.fq.gz'), '@r1 x\nACGTAC\n+\nIIIIII\n\n'
                    '@r2\nAC\n+r2\nII\n', gz=True)
        info = U.reads_to_fasta(fq, os.path.join(d, 'r.fa'))
        self.assertEqual((info['format'], info['n_reads'], info['total_bp']), ('fastq', 2, 8))
        with open(os.path.join(d, 'r.fa')) as fh:
            self.assertEqual(fh.read(), '>r1\nACGTAC\n>r2\nAC\n')
        with self.assertRaises(ValueError):
            U.reads_to_fasta(_write(os.path.join(d, 'bad.txt'), 'hello\n'),
                             os.path.join(d, 'x.fa'))

    def test_f09_summary_warnings(self):
        before = {'internal_gaps': 4, 'gap_bp': 400, 'total_length': 10400,
                  'ids': ['a', 'b'], 'scaftig_n50': 1}
        after = dict(before, ids=['a'])
        tool = {'version': '1.2.1', 'seconds': 1}
        reads = {'total_bp': 20000, 'mean_length': 500}
        _, w = U.summarize(before, after, [], reads, tool, False)
        text = ' '.join(w)
        for needle in ('No gaps were closed', 'coverage', 'Mean read length',
                       'IDs changed', 'junctions'):
            self.assertIn(needle, text)

    def test_f10_stale_tags_refused(self):
        d = tempfile.mkdtemp(dir=self.scratch)
        open(os.path.join(d, 'done_step1_tag'), 'w').close()
        with self.assertRaisesRegex(RuntimeError, 'silently skip'):
            U.run_tgsgapcloser(['/bin/false'], d)

    # --------------------------------------------------- known-truth tier
    def test_kt0_racon_run_outputs(self):
        run = self.run_mode('racon')
        for key in ('report_name', 'report_ref', 'assembly_ref', 'gapfill_summary'):
            self.assertTrue(run['out'].get(key), key)
        self.assertEqual(run['results']['summary']['gaps_before'], 12)

    def test_kt1_planted_gaps_by_id(self):
        gaps = self.stage()['manifest']['gaps']
        bridgeable = [g for g in gaps if g['bridging_reads'] >= 3]
        if len(bridgeable) < MIN_BRIDGEABLE_GAPS:
            self.skipTest('only %d bridgeable gaps in fixture' % len(bridgeable))
        events = self.run_mode('racon')['results']['events']
        closed = 0
        for g in gaps:
            hit = [e for e in events if e['scaffold'] == g['scaffold']
                   and e['input_gap_start'] <= g['end'] and e['input_gap_end'] >= g['start']]
            self.assertEqual(len(hit), 1, 'planted gap %s maps to %d events' % (g['id'], len(hit)))
            closed += g in bridgeable and hit[0]['status'] in ('filled', 'merged')
        print('closed %d / %d bridgeable planted gaps' % (closed, len(bridgeable)))
        self.assertGreaterEqual(closed / float(len(bridgeable)), MIN_CLOSURE_OF_BRIDGEABLE)

    def test_kt2_control_contig_unchanged(self):
        out = dict(U.iter_fasta(self.run_mode('racon')['fasta']))
        truth = dict(U.iter_fasta(os.path.join(self.stage()['dir'], 'truth.fa')))
        self.assertEqual(sorted(out), sorted(truth))
        self.assertEqual(out['control'], truth['control'])
        for sid in truth:
            self.assertLess(abs(len(out[sid]) - len(truth[sid])), 0.01 * len(truth[sid]))

    def test_kt3_fill_accuracy_racon(self):
        nm, aligned = self.edit_distance(self.run_mode('racon')['fasta'])
        print('racon: %d edits over %d planted bp, aligned %.3f'
              % (nm, self.planted_bp(), aligned))
        self.assertGreaterEqual(aligned, MIN_ALIGNED_FRACTION)
        self.assertLessEqual(nm / float(self.planted_bp()), MAX_EDITS_PER_PLANTED_BP)

    def test_kt4_correction_beats_no_correction(self):
        ne = self.run_mode('ne')
        self.assertIn('gaps_closed', ne['out']['gapfill_summary'])
        nm_ne, _ = self.edit_distance(ne['fasta'])
        nm_racon, _ = self.edit_distance(self.run_mode('racon')['fasta'])
        print('edits vs truth: racon %d, no correction %d' % (nm_racon, nm_ne))
        self.assertLess(nm_racon, nm_ne)

    def test_kt5_contigs_source_rebuilds_truth(self):
        out = dict(U.iter_fasta(self.run_mode('contigs')['fasta']))
        truth = dict(U.iter_fasta(os.path.join(self.stage()['dir'], 'truth.fa')))
        nm, aligned = self.edit_distance(self.run_mode('contigs')['fasta'])
        print('contigs source: %d edits, aligned %.3f' % (nm, aligned))
        self.assertEqual(sorted(out), sorted(truth))
        self.assertLessEqual(nm, 10)

    def test_kt6_no_gaps_is_a_readable_error(self):
        st = self.stage()
        with self.assertRaisesRegex(ValueError, 'no internal N-gaps'):
            self.impl.run_tgsgapcloser(self.ctx, dict(
                workspace_name=self.wsName, assembly_ref=st['refs']['truth'],
                reads_ref=st['refs']['reads'], output_assembly_name='kt_nogaps'))

    def test_kt7_wrong_reads_type_is_readable(self):
        st = self.stage()
        report_ref = self.run_mode('racon')['out']['report_ref']
        with self.assertRaisesRegex(ValueError, 'single-end long reads'):
            self.impl.run_tgsgapcloser(self.ctx, dict(
                workspace_name=self.wsName, assembly_ref=st['refs']['draft'],
                reads_ref=report_ref, output_assembly_name='kt_badtype'))


if __name__ == '__main__':
    unittest.main()
