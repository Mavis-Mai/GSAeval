import argparse
import importlib.util
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

spec=importlib.util.spec_from_file_location('csqcompare',Path(__file__).with_name('compare_CSQ_components.py'))
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def scores(ids, values):
    return pd.DataFrame({'BUSCO':ids, 'GeneA':['t'+str(i) for i in ids], 'GeneB':['r'+str(i) for i in ids],
                         **{metric:values for metric in m.METRICS}})


class ComparisonTests(unittest.TestCase):
    def test_signed_absolute_change_and_missing_pairs(self):
        before=scores(['1','2','3'],[.9,.5,.7])
        after=scores(['1','2','4'],[.7,.6,.1])
        pairs=m.pair_tables(before,after,['BUSCO'])
        summary=m.summarize(pairs,large_change=.15)
        row=summary.iloc[0]
        self.assertEqual(row.paired_valid_n,2)
        self.assertAlmostEqual(row.mean_delta,-.05)
        self.assertAlmostEqual(row.mean_abs_delta,.15)
        self.assertEqual(row.decreased_pct,50)
        self.assertEqual(row.large_decrease_pct,50)
        self.assertEqual(pairs.pair_status.value_counts().to_dict(),{'paired':2,'before_only':1,'after_only':1})

    def test_missing_is_not_zero_and_reference_switch(self):
        before=scores(['1','2'],[1,np.nan])
        after=scores(['1','2'],[np.nan,.5])
        after.loc[0,'GeneB']='different'
        pairs=m.pair_tables(before,after,['BUSCO'])
        row=m.summarize(pairs).iloc[0]
        self.assertEqual(row.paired_valid_n,0)
        self.assertEqual(row.paired_missing_n,2)
        self.assertTrue(np.isnan(row.mean_abs_delta))
        self.assertEqual(int(pairs.reference_changed.sum()),1)

    def test_source_mapping_preserves_fusion_edges(self):
        before=scores(['1','2'],[1,.8])
        after=scores(['1','2'],[.5,.4])
        after['GeneA']='fusion'
        mapping=pd.DataFrame({'transcript_id':['fusion','fusion'],'source_transcript_id':['t1','t2'],
                              'error_type':['gene_fusion']*2,'error_status':['applied']*2})
        pairs=m.pair_tables(before,after,['BUSCO'],mapping)
        self.assertEqual(len(pairs),4)
        matched=pairs.loc[pairs.pair_status.eq('paired')]
        self.assertEqual(len(matched),2)
        self.assertEqual(set(matched.GeneA_before),{'t1','t2'})
        self.assertTrue(np.allclose(sorted(matched[m.METRICS[0]+'_delta']),[-.5,-.4]))

    def test_headerless_preserves_first_row(self):
        with tempfile.TemporaryDirectory() as root:
            path=Path(root)/'a.tsv'
            path.write_text('BUSCO1\tt1\tr1\t2\t3\t0.8\t0.7\t0.6\t0.5\t0.9\t0.7\tGood\n')
            table=m.read_scores(path,['BUSCO'])
            self.assertEqual(len(table),1)
            self.assertEqual(table.iloc[0].BUSCO,'BUSCO1')
            self.assertEqual(table.iloc[0].CDS_number_score,.8)

    def test_reject_duplicate_keys_and_out_of_range(self):
        with tempfile.TemporaryDirectory() as root:
            path=Path(root)/'a.tsv'
            scores(['1','1'],[.5,.6]).to_csv(path,sep='\t',index=False)
            with self.assertRaisesRegex(ValueError,'duplicate'):
                m.read_scores(path,['BUSCO'])
            scores(['1'],[1.1]).to_csv(path,sep='\t',index=False)
            with self.assertRaisesRegex(ValueError,'range'):
                m.read_scores(path,['BUSCO'])

    def test_runner_single_and_paired(self):
        with tempfile.TemporaryDirectory() as root:
            root=Path(root)
            before,after=root/'before.tsv',root/'after.tsv'
            scores(['1','2'],[1,.8]).to_csv(before,sep='\t',index=False)
            scores(['1','2'],[.8,.7]).to_csv(after,sep='\t',index=False)
            args=argparse.Namespace(before=str(before),after=str(after),key=['BUSCO'],error_report=None,
                mapping_level='transcript',tolerance=1e-9,large_change=.1,output_dir=str(root/'comparison'),
                no_plots=True,before_label='Original',after_label='With errors')
            _,summary,audit=m.run(args)
            self.assertEqual(audit['pair_status_counts'],{'paired':2})
            self.assertTrue((root/'comparison/paired_scores.tsv').exists())
            self.assertTrue(np.allclose(summary.mean_abs_delta,.15))
            args.before=None
            args.output_dir=str(root/'single')
            _,summary,audit=m.run(args)
            self.assertEqual(len(summary),4)
            self.assertEqual(audit['after_rows'],2)


if __name__=='__main__':
    unittest.main()
