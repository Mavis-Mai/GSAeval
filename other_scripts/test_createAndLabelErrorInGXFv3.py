import importlib.util
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location('modifier', Path(__file__).with_name('createAndLabelErrorInGXFv3.py'))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def feature(kind, start, end, strand='+', phase='.'):
    return m.Feature(f'chr1\ttest\t{kind}\t{start}\t{end}\t.\t{strand}\t{phase}\tgene_id "g"; transcript_id "t";\n')


class ModifierTests(unittest.TestCase):
    def test_fast_alignment_matches_reference_including_ties(self):
        import itertools
        import random
        words = [''.join(chars) for size in range(1, 4)
                 for chars in itertools.product('MA', repeat=size)]
        pairs = [(a, b) for a in words for b in words]
        rng = random.Random(42)
        pairs += [(''.join(rng.choices('MACD', k=rng.randint(65, 100))),
                   ''.join(rng.choices('MACD', k=rng.randint(65, 100)))) for _ in range(100)]
        for a, b in pairs:
            expected = m._protein_alignment_python(a, b)
            self.assertEqual(m.protein_alignment(a, b, 'python'), expected)
            self.assertEqual(m.protein_alignment(a, b, 'auto'), expected)
            if m.np is not None:
                self.assertEqual(m._protein_alignment_numpy(a, b), expected)

    def test_alignment_without_numpy(self):
        from unittest.mock import patch
        a, b = 'MACD'*30, 'MCCD'*30
        expected = m._protein_alignment_python(a, b)
        m.protein_alignment.cache_clear()
        with patch.object(m, 'np', None):
            self.assertEqual(m.protein_alignment(a, b, 'auto'), expected)
            with self.assertRaisesRegex(ValueError, 'requires NumPy'):
                m.protein_alignment(a, b, 'numpy')
        m.protein_alignment.cache_clear()

    def test_added_exon_keeps_two_base_gaps(self):
        for strand in ('+', '-'):
            for right_start in (18, 19, 50):
                exons = [feature('exon', 1, 10, strand),
                         feature('exon', right_start, right_start+10, strand)]
                for seed in range(100):
                    result = m.FastGTFModifier(seed=seed)._mutate(exons, 'add_exon')
                    self.assertIsNotNone(result)
                    self.assertEqual(len(result), 3)
                    self.assertGreaterEqual(result[1].start-result[0].end-1, 2)
                    self.assertGreaterEqual(result[2].start-result[1].end-1, 2)
                    self.assertGreaterEqual(result[1].end-result[1].start+1, 3)
                    self.assertEqual((result[0].start, result[0].end), (1, 10))
                    self.assertEqual((result[2].start, result[2].end),
                                     (right_start, right_start+10))
            # Six intronic bases cannot fit 2 + 3 + 2 bases.
            self.assertIsNone(m.FastGTFModifier(seed=0)._mutate(
                [feature('exon', 1, 10, strand), feature('exon', 17, 30, strand)],
                'add_exon'))

    def make_gene(self, gid, tid, intervals, strand='+'):
        rows = [feature('transcript', min(a for a, b in intervals), max(b for a, b in intervals), strand)]
        rows += [feature('exon', a, b, strand) for a, b in intervals]
        for row in rows:
            row.attrs = {'gene_id': gid, 'transcript_id': tid, 'protein_id': 'stale', 'Name': 'stale'}
        return rows

    def test_prediction_rebuilds_phase_and_codons_on_both_strands(self):
        # Transcript ATGA + AATAA: CDS 6 bp, with a codon spanning the intron.
        for strand, dna in (('+', 'ATGACCCCCAATAA'), ('-', 'TTATTCCCCCTCAT')):
            modifier = m.FastGTFModifier(seed=1)
            modifier.genome = {'chr1': dna}
            modifier._used_ids = set()
            intervals = [(1, 4), (10, 14)] if strand == '+' else [(1, 5), (11, 14)]
            exons = [feature('exon', a, b, strand) for a, b in intervals]
            rows = modifier._predict(exons, 'newg', 'newt', [('g', 't')], 'event1')
            cds = m.ordered([f for f in rows if f.kind == 'CDS'])
            self.assertEqual(m.coding_dna(cds, modifier.genome), 'ATGAAA')
            self.assertEqual([f.fields[7] for f in cds], ['0', '2'])
            self.assertEqual(sum(f.end-f.start+1 for f in rows if f.kind == 'stop_codon'), 3)
            self.assertEqual(sum(f.end-f.start+1 for f in rows if f.kind == 'start_codon'), 3)
            self.assertEqual(rows[0].attrs['orf_status'], 'complete')
            self.assertTrue(all(f.attrs['gene_id'] == 'newg' and f.attrs['transcript_id'] == 'newt' for f in rows))
            self.assertTrue(all('Name' not in f.attrs for f in rows))

    def test_orf_absent_partial_and_minimum(self):
        modifier = m.FastGTFModifier()
        modifier._used_ids = set()
        for dna, expected in [('CCCCCC', 'no_orf'), ('ATGAAA', 'partial_3prime'), ('ATGAAATAA', 'complete')]:
            modifier.genome = {'chr1': dna}
            rows = modifier._predict([feature('exon', 1, len(dna))], 'g', 't', [('g', 't')], 'e')
            self.assertEqual(rows[0].attrs['orf_status'], expected)
        modifier.min_orf_aa = 3
        self.assertEqual(modifier._predict([feature('exon', 1, 9)], 'g', 't', [('g', 't')], 'e')[0].attrs['orf_status'], 'no_orf')

    def test_gene_fusion_all_isoforms_and_report(self):
        import csv
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, target = root/'source.gtf', root/'target.gtf'
            rows = []
            for gid, tid, start in [('g1', 't1', 1), ('g1', 't1b', 1), ('g2', 't2', 21), ('g2', 't2b', 21)]:
                rows += self.make_gene(gid, tid, [(start, start+8)])
                cds = feature('CDS', start, start+5, phase='0')
                cds.attrs = {'gene_id': gid, 'transcript_id': tid}
                rows.append(cds)
            source.write_text(''.join(f.line() for f in rows))
            modifier = m.FastGTFModifier(1, seed=1, error_types=['gene_fusion'])
            modifier.genome = {'chr1': 'ATGAAATAA'+'C'*11+'ATGCCCTAA'}
            modifier.modify_gtf(source, target)
            output = [m.Feature(line) for line in target.read_text().splitlines()]
            genes = [f for f in output if f.kind == 'gene']
            transcripts = [f for f in output if f.kind == 'transcript']
            self.assertEqual(len(genes), 1)
            self.assertEqual(len(transcripts), 4)
            self.assertNotIn(genes[0].attrs['gene_id'], ('g1', 'g2'))
            self.assertEqual((genes[0].start, genes[0].end), (1, 29))
            self.assertEqual(len({t.attrs['transcript_id'] for t in transcripts}), 4)
            with open(str(target)+'.error_report.tsv') as handle:
                report = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(len(report), 8)
            self.assertEqual({r['source_transcript_id'] for r in report}, {'t1', 't1b', 't2', 't2b'})
            self.assertTrue(all(r['error_type'] == 'gene_fusion' for r in report))
            self.assertTrue(all(t.attrs['comparison_status'] == 'multiple_sources_see_report' for t in transcripts))
            self.assertTrue(all(r['frameshift'] == 'NA' for r in report))
            self.assertTrue(all(r['coding_impact'] == 'protein_changed' for r in report))
            self.assertTrue(all(float(r['protein_identity_pct']) < 100 for r in report))

    def test_gene_split_and_negative_order(self):
        for strand in ('+', '-'):
            modifier = m.FastGTFModifier(1, seed=1, error_types=['gene_split'], require_protein_change=False)
            modifier.genome = {'chr1': 'ATGAAATAA'+'C'*11+'ATGCCCTAA'}
            modifier._used_ids = {'g', 't', 'g.split1'}
            groups = {'t': self.make_gene('g', 't', [(1, 9), (21, 29)], strand)}
            plans = modifier._plan_events(groups)
            self.assertEqual(len(plans), 2)
            self.assertEqual(len({p[0][0].attrs['gene_id'] for p in plans}), 2)
            self.assertTrue(all(p[2] == 'gene_split' and p[1] == ['t'] for p in plans))
            first, second = plans[0][0][0], plans[1][0][0]
            self.assertTrue(first.end < second.start if strand == '+' else first.start > second.end)
            self.assertNotEqual(first.attrs['gene_id'], 'g.split1')

    def test_split_roundtrip_removes_source_gene_and_rebuilds_hierarchy(self):
        import csv
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output = root/'in.gtf', root/'out.gtf'
            rows = self.make_gene('g', 't', [(1, 9), (21, 29)])
            gene = feature('gene', 1, 29)
            gene.attrs = {'gene_id': 'g', 'gene_name': 'old_name'}
            rows.append(gene)
            for start in (1, 21):
                rows.append(feature('CDS', start, start+8, phase='0'))
            source.write_text(''.join(f.line() for f in rows))
            modifier = m.FastGTFModifier(1, seed=1, error_types=['gene_split'], require_protein_change=False)
            modifier.genome = {'chr1': 'ATGAAATAA'+'C'*11+'ATGCCCTAA'}
            modifier.modify_gtf(source, output)
            result = [m.Feature(line) for line in output.read_text().splitlines()]
            genes = [f for f in result if f.kind == 'gene']
            self.assertEqual(len(genes), 2)
            self.assertTrue(all(f.attrs['gene_id'] != 'g' for f in result))
            for gene in genes:
                self.assertNotIn('transcript_id', gene.attrs)
                children = [f for f in result if f.kind != 'gene' and f.attrs['gene_id'] == gene.attrs['gene_id']]
                self.assertEqual(len([f for f in children if f.kind == 'transcript']), 1)
                self.assertEqual((gene.start, gene.end), (min(f.start for f in children), max(f.end for f in children)))
                self.assertEqual(gene.attrs['error_type'], 'gene_split')
                cds = [f for f in children if f.kind == 'CDS']
                self.assertEqual(len(m.coding_dna(cds, modifier.genome)), 6)
                self.assertEqual(cds[0].fields[7], '0')
                self.assertTrue(all('Name' not in f.attrs for f in children))
            with open(str(output)+'.error_report.tsv') as handle:
                report = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(len(report), 2)
            self.assertTrue(all(r['source_transcript_id'] == 't' for r in report))

    def test_fusion_rejects_incompatible_neighbors(self):
        for strand, start in [('-', 21), ('+', 8), ('+', 30)]:
            modifier = m.FastGTFModifier(1, seed=1, error_types=['gene_fusion'], max_fusion_gap=15, require_protein_change=False)
            modifier.genome = {'chr1': 'C'*50}
            modifier._used_ids = set()
            groups = {'t1': self.make_gene('g1', 't1', [(1, 9)]),
                      't2': self.make_gene('g2', 't2', [(start, start+8)], strand)}
            plans = modifier._plan_events(groups)
            self.assertTrue(all(p[3] == 'no_applicable_error' for p in plans))

    def test_split_cannot_cut_through_alternative_exon(self):
        modifier = m.FastGTFModifier(1, seed=1, error_types=['gene_split'], require_protein_change=False)
        modifier.genome = {'chr1': 'C'*40}
        modifier._used_ids = set()
        groups = {'t1': self.make_gene('g', 't1', [(1, 9), (21, 29)]),
                  't2': self.make_gene('g', 't2', [(1, 29)])}
        self.assertTrue(all(p[3] == 'no_applicable_error' for p in modifier._plan_events(groups)))

    def test_default_retries_utr_only_change_and_rolls_back_ids(self):
        from unittest.mock import patch
        modifier = m.FastGTFModifier(1, seed=1, error_types=['exon_truncation'], progress_every=0, max_attempts=3)
        modifier.genome = {'chr1': 'CCCATGAAACCCTAACCC'}
        original = self.make_gene('g', 't', [(1, 17)])+[feature('CDS', 4, 15, phase='0')]
        modifier._used_ids = {'g', 't'}
        with patch.object(modifier, '_mutate', side_effect=[[feature('exon', 2, 17)], [feature('exon', 1, 10)]]) as mutate:
            plans = modifier._plan_events({'t': original})
        self.assertEqual(mutate.call_count, 2)
        self.assertEqual(plans[0][3], 'applied')
        self.assertEqual(plans[0][0][0].attrs['transcript_id'], 'g.exon_truncation.t1')
        self.assertTrue(modifier._protein_changed(plans[0][0], ['t'], {'t': original}))

    def test_unchanged_candidate_is_not_marked_applied(self):
        from unittest.mock import patch
        modifier = m.FastGTFModifier(1, seed=1, error_types=['exon_truncation'], progress_every=0, max_attempts=3)
        modifier.genome = {'chr1': 'CCCATGAAACCCTAACCC'}
        original = self.make_gene('g', 't', [(1, 17)])+[feature('CDS', 4, 15, phase='0')]
        modifier._used_ids = {'g', 't'}
        with patch.object(modifier, '_mutate', return_value=[feature('exon', 2, 17)]) as mutate:
            plans = modifier._plan_events({'t': original})
        self.assertEqual(mutate.call_count, 3)
        self.assertEqual(plans[0][2:4], ('none', 'no_protein_changing_error'))
        self.assertEqual(modifier._used_ids, {'g', 't'})
        self.assertEqual([f.line() for f in plans[0][0]], [f.line() for f in original])

    def test_strict_mode_skips_unevaluable_sources(self):
        for dna, cds in [('CCCCCC', []), ('ATGNNN', [feature('CDS', 1, 6, phase='0')])]:
            modifier = m.FastGTFModifier(1, seed=1, progress_every=0)
            modifier.genome = {'chr1': dna}
            modifier._used_ids = {'g', 't'}
            groups = {'t': self.make_gene('g', 't', [(1, 6)])+cds}
            plans = modifier._plan_events(groups)
            self.assertEqual(plans[0][3], 'no_evaluable_source_protein')

    def test_strict_split_rejects_one_unchanged_product_atomically(self):
        modifier = m.FastGTFModifier(1, seed=1, error_types=['gene_split'], progress_every=0, max_attempts=2)
        modifier.genome = {'chr1': 'ATGAAATAA'+'C'*11+'ATGCCCTAA'}
        modifier._used_ids = {'g', 't'}
        original = self.make_gene('g', 't', [(1, 9), (21, 29)])
        original += [feature('CDS', 1, 9, phase='0'), feature('CDS', 21, 29, phase='0')]
        plans = modifier._plan_events({'t': original})
        self.assertEqual(len(plans), 1)
        self.assertEqual(plans[0][3], 'no_protein_changing_error')
        self.assertEqual(modifier._used_ids, {'g', 't'})

    def test_each_error_accepts_only_changed_protein(self):
        for error in m.FastGTFModifier.EXON_ERRORS+('gene_split',):
            accepted = 0
            for seed in range(5):
                modifier = m.FastGTFModifier(1, seed=seed, error_types=[error], progress_every=0)
                modifier.genome = {'chr1': 'C'*100+'ATGAAACCCGGG'+'C'*20+'AAATTTCCCTAA'+'C'*100}
                modifier._used_ids = {'g', 't'}
                original = self.make_gene('g', 't', [(101, 112), (133, 144)])
                original += [feature('CDS', 101, 112, phase='0'), feature('CDS', 133, 144, phase='0')]
                groups = {'t': original}
                for output, sources, actual_error, status, event_id in modifier._plan_events(groups):
                    if status == 'applied':
                        accepted += 1
                        self.assertEqual(actual_error, error)
                        labels = modifier._label(original, output)
                        self.assertEqual(labels['coding_impact'], 'protein_changed')
                        self.assertLess(float(labels['protein_identity_pct']), 100)
            self.assertGreater(accepted, 0, error)

    def test_identity(self):
        self.assertEqual(m.protein_identity('MAAA', 'MAA'), 75)
        self.assertEqual(m.protein_identity('MA', 'MT'), 50)
        self.assertEqual(m.protein_identity('MA', ''), 0)

    def test_phase_and_negative_strand(self):
        # ATG + AAA across an intron; phase on the second CDS must not be trimmed.
        plus = [feature('CDS', 1, 4, phase='0'), feature('CDS', 8, 9, phase='2')]
        self.assertEqual(m.coding_dna(plus, {'chr1': 'ATGACCCAA'}), 'ATGAAA')
        minus = [feature('CDS', 1, 2, '-', '2'), feature('CDS', 6, 9, '-', '0')]
        self.assertEqual(m.coding_dna(minus, {'chr1': 'TTCCCTCAT'}), 'ATGAAA')

    def test_coding_and_utr(self):
        mod = m.FastGTFModifier(seed=3)
        mod.genome = {'chr1': 'CCCATGAAACCCTAACCC'}
        original = [feature('exon', 1, 17), feature('CDS', 4, 15, phase='0')]
        mod._used_ids = set()
        utr = mod._predict([feature('exon', 2, 17)], 'g', 'utr', [('g', 't')], 'e')
        self.assertEqual(mod._label(original, utr)['coding_impact'], 'protein_unchanged')
        coding = mod._predict([feature('exon', 1, 7), feature('exon', 9, 17)], 'g', 'coding', [('g', 't')], 'e')
        self.assertEqual(mod._label(original, coding)['coding_impact'], 'protein_changed')
        self.assertEqual(mod._label(original, coding)['frameshift'], 'NA')
        mod.genome = None
        self.assertEqual(mod._label(original, coding)['coding_impact'], 'unknown')
        self.assertEqual(mod._label([], [])['coding_impact'], 'noncoding')

    def test_each_error_and_reproducibility(self):
        exons = [feature('exon', 101, 160), feature('exon', 201, 260)]
        for error in m.FastGTFModifier.EXON_ERRORS:
            mod = m.FastGTFModifier(1, seed=4, error_types=[error])
            mutated = mod._mutate(exons, error)
            self.assertIsNotNone(mutated, error)
            self.assertNotEqual([(f.start, f.end) for f in exons], [(f.start, f.end) for f in mutated])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root/'input.gtf'
            source.write_text(''.join(f.line() for f in [feature('transcript', 101, 260)] + exons +
                                     [feature('CDS', 110, 160, phase='0'), feature('CDS', 201, 250, phase='0')]))
            for name in ('a', 'b'):
                modifier = m.FastGTFModifier(1, seed=7)
                modifier.genome = {'chr1': 'ATG' * 100}
                modifier.modify_gtf(source, root/name)
            self.assertEqual((root/'a').read_text(), (root/'b').read_text())
            self.assertIn('\tCDS\t', (root/'a').read_text())
            self.assertIn('protein_', (root/'a.error_report.tsv').read_text())
            self.assertIn('>t', (root/'a.original_proteins.fa').read_text())
            self.assertTrue((root/'a.modified_proteins.fa').exists())

    def test_required_fasta_and_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output, fasta = root/'in.gtf', root/'out.gtf', root/'genome.fa'
            source.write_text(feature('CDS', 1, 9, phase='0').line())
            with self.assertRaisesRegex(ValueError, 'FASTA is required'):
                m.FastGTFModifier(0).modify_gtf(source, output)
            fasta.write_text('>chr1\nATGAAATAG\n')
            modifier = m.FastGTFModifier(0, genome=fasta)
            modifier.modify_gtf(source, output)
            self.assertIn('protein_unchanged', Path(str(output)+'.error_report.tsv').read_text())
            self.assertEqual(Path(str(output)+'.original_proteins.fa').read_text(), '>t\nMK\n')
            source.write_text(feature('CDS', 1, 10, phase='0').line())
            with self.assertRaisesRegex(ValueError, 'bounds'):
                modifier.modify_gtf(source, root/'bad.gtf')
            self.assertFalse((root/'bad.gtf').exists())
            fasta.write_text('>other\nATGAAATAG\n')
            with self.assertRaisesRegex(ValueError, 'missing from FASTA'):
                m.FastGTFModifier(0, genome=fasta).modify_gtf(source, root/'bad.gtf')

    def test_no_error_preserves_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory)/'in.gtf', Path(directory)/'out.gtf'
            cds = feature('CDS', 1, 9, phase='1')
            cds.fields[5] = '5.7'
            source.write_text(cds.line())
            modifier = m.FastGTFModifier(0)
            modifier.genome = {'chr1': 'ATGAAATAG'}
            modifier.modify_gtf(source, output)
            result = next(m.Feature(line) for line in output.read_text().splitlines() if '\tCDS\t' in line)
            self.assertEqual(result.fields[5:8], ['5.7', '+', '1'])


if __name__ == '__main__':
    unittest.main()
