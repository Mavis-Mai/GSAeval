import argparse
import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

spec=importlib.util.spec_from_file_location('contribution',Path(__file__).with_name('detect_metric_contribution.py'))
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def table(values):
    count=len(values)
    return pd.DataFrame({'BUSCO':[f'b{i}' for i in range(count)],
                         'GeneA':[f'g{i}' for i in range(count)], 'GeneB':[f'r{i}' for i in range(count)],
                         **{name:np.asarray(values)[:,i] for i,name in enumerate(m.METRICS)},
                         'CSQ':['ignored_text']*count, 'label':['ignored_label']*count})


class ContributionTests(unittest.TestCase):
    def test_actual_change_summary_direction_magnitude_and_threshold(self):
        before=table(np.tile(np.array([.2,.8,.5,.6])[:,None],(1,5))).drop(columns=['CSQ','label'])
        after=table(np.tile(np.array([.6,.4,.5,.8])[:,None],(1,5))).drop(columns=['CSQ','label'])
        pairs=m.pair_tables(before,after,['BUSCO'])
        summary=m.change_distribution_summary(pairs,'all_pairs',large_change=.2)
        row=summary.iloc[0]
        self.assertAlmostEqual(row.gain_mean,.05)
        self.assertAlmostEqual(row.mean_abs_gain,.25)
        self.assertEqual(row.improved_pct,50)
        self.assertEqual(row.worsened_pct,25)
        self.assertEqual(row.unchanged_pct,25)
        self.assertEqual(row.large_improvement_n,2)
        self.assertEqual(row.large_absolute_change_n,3)
        self.assertAlmostEqual(row.gain_median,.1)

    def test_identical_marginal_distributions_can_hide_paired_changes(self):
        before=table(np.tile(np.array([.2,.8])[:,None],(1,5))).drop(columns=['CSQ','label'])
        after=table(np.tile(np.array([.8,.2])[:,None],(1,5))).drop(columns=['CSQ','label'])
        summary=m.change_distribution_summary(m.pair_tables(before,after,['BUSCO']),'all_pairs')
        self.assertTrue(np.allclose(summary.error_mean,summary.correct_mean))
        self.assertTrue(np.allclose(summary.gain_mean,0))
        self.assertTrue(np.allclose(summary.mean_abs_gain,.6))

    def test_ecdf_preserves_zero_mass_and_endpoints(self):
        x,y=m.ecdf_coordinates([-1,0,0,1],-1,1)
        self.assertEqual(y[np.flatnonzero(x<=0)[-1]],.75)
        self.assertEqual(y[np.flatnonzero(x<=-.1)[-1]],.25)
        self.assertEqual(y[-1],1)
        with self.assertRaises(ValueError):
            m.ecdf_coordinates([],0,1)
        with self.assertRaises(ValueError):
            m.ecdf_coordinates([1.1],0,1)

    def test_auc_direction_and_ties(self):
        self.assertEqual(m.auc([0,0],[1,1]),1)
        self.assertEqual(m.auc([1,1],[0,0]),0)
        self.assertEqual(m.auc([.5,.5],[.5,.5]),.5)
        self.assertEqual(m.auc([0,1],[1,1]),.75)
        x,y=m.roc(np.array([0,1]),np.array([1,1]))
        self.assertAlmostEqual(float(np.sum(np.diff(x)*(y[1:]+y[:-1])/2)),.75)
        self.assertEqual(m.auc([.1+.2],[.3]),.5)

    def test_unique_and_redundant_signal_ablation(self):
        error=np.full((4,5),.5)
        correct=error.copy()
        error[:,0]=.1; correct[:,0]=.9
        result=m.point_statistics(error,correct,1e-9)
        self.assertEqual(result['full_auc'],1)
        self.assertEqual(result['auc_loss_when_removed'][0],.5)
        self.assertTrue(np.all(result['auc_loss_when_removed'][1:]==0))
        result=m.point_statistics(np.zeros((4,5)),np.ones((4,5)),1e-9)
        self.assertTrue(np.all(result['single_auc']==1))
        self.assertTrue(np.all(result['auc_loss_when_removed']==0))
        expected=(result['mean_gain']-result['full_mean_gap'])/4
        self.assertTrue(np.allclose(result['mean_gap_loss_when_removed'],expected))

    def test_conflicting_metric_can_have_negative_contribution(self):
        error=np.tile([.2,.2,.2,.2,1],(4,1))
        correct=np.tile([.4,.4,.4,.4,0],(4,1))
        result=m.point_statistics(error,correct,1e-9)
        self.assertEqual(result['full_auc'],0)
        self.assertEqual(result['single_auc'][-1],0)
        self.assertEqual(result['auc_loss_when_removed'][-1],-1)

    def test_nonconstant_bootstrap_is_paired_and_reproducible(self):
        rng=np.random.default_rng(7)
        error=table(rng.uniform(.1,.5,size=(20,5))).drop(columns=['CSQ','label'])
        correct=error.copy()
        correct[m.METRICS]=np.clip(error[m.METRICS].to_numpy()+rng.normal(.1,.08,size=(20,5)),0,1)
        pairs=m.pair_tables(error,correct,['BUSCO'])
        first=m.analyze_scope(pairs,'all_pairs',['BUSCO'],None,40,42,1e-9,False)
        second=m.analyze_scope(pairs,'all_pairs',['BUSCO'],None,40,42,1e-9,False)
        pd.testing.assert_frame_equal(first[0],second[0])
        pd.testing.assert_frame_equal(first[1],second[1])
        self.assertTrue(first[0]['mean_gain_ci_high'].ge(first[0]['mean_gain_ci_low']).all())
        expected=(first[0]['mean_gain'].to_numpy()-first[2]['full_mean_gap'])/4
        self.assertTrue(np.allclose(first[1]['mean_gap_loss_when_removed'],expected))

    def test_headerless_and_ignored_csq_label(self):
        with tempfile.TemporaryDirectory() as root:
            path=Path(root)/'scores.tsv'
            path.write_text('b1\tg1\tr1\t2\t3\t0.1\t0.2\t0.3\t0.4\t0.5\tIGNORED\tUNKNOWN\n')
            frame=m.read_table(path,['BUSCO'])
            self.assertEqual(len(frame),1)
            self.assertEqual(frame.iloc[0]['Total_length_score'],.5)
            self.assertNotIn('CSQ',frame)
            self.assertNotIn('label',frame)

    def test_reject_duplicates_invalid_and_out_of_range(self):
        with tempfile.TemporaryDirectory() as root:
            path=Path(root)/'scores.tsv'
            frame=table(np.full((2,5),.5)); frame['BUSCO']='duplicate'
            frame.to_csv(path,sep='\t',index=False)
            with self.assertRaisesRegex(ValueError,'duplicate'):
                m.read_table(path,['BUSCO'])
            for value in [1.1,'bad','inf']:
                frame=table(np.full((2,5),.5))
                frame[m.METRICS[0]]=str(value)
                frame.to_csv(path,sep='\t',index=False)
                with self.assertRaises(ValueError):
                    m.read_table(path,['BUSCO'])

    def test_bootstrap_reproducibility_clusters_and_constant_correlation(self):
        error=table(np.zeros((4,5))).drop(columns=['CSQ','label'])
        correct=table(np.ones((4,5))).drop(columns=['CSQ','label'])
        error['event_id']=correct['event_id']=['e1','e1','e2','e2']
        paired=m.pair_tables(error,correct,['BUSCO'])
        first=m.analyze_scope(paired,'all_pairs',['BUSCO'],'event_id',20,42,1e-9,False)
        second=m.analyze_scope(paired,'all_pairs',['BUSCO'],'event_id',20,42,1e-9,False)
        pd.testing.assert_frame_equal(first[0],second[0])
        self.assertEqual(first[2]['bootstrap_units_n'],2)
        self.assertTrue(first[3].isna().all().all())
        self.assertTrue(first[1]['auc_loss_when_removed_ci_low'].eq(0).all())
        paired['event_id_before']=paired['event_id_after']='one'
        result=m.analyze_scope(paired,'all_pairs',['BUSCO'],'event_id',20,42,1e-9,False)
        self.assertTrue(result[0]['single_auc_ci_low'].isna().all())
        paired.loc[0,'event_id_after']='other'
        with self.assertRaisesRegex(ValueError,'agree'):
            m.cluster_labels(paired,'event_id',['BUSCO'])

    def test_complete_cases_and_stable_scope_end_to_end(self):
        with tempfile.TemporaryDirectory() as root:
            root=Path(root)
            before,after=root/'error.tsv',root/'correct.tsv'
            error=table(np.full((4,5),.2))
            correct=table(np.full((4,5),.8))
            correct.loc[3,m.METRICS[0]]=np.nan
            correct.loc[2,'GeneB']='different'
            error.to_csv(before,sep='\t',index=False)
            correct.to_csv(after,sep='\t',index=False)
            args=argparse.Namespace(before=str(before),after=str(after),output_dir=str(root/'out'),
                key=['BUSCO'],scope='both',bootstrap=20,seed=42,cluster_column=None,tolerance=1e-9,no_plots=True)
            response,ablation,audit=m.run(args)
            self.assertEqual(audit['complete_pairs'],3)
            self.assertEqual(audit['excluded_incomplete_pairs'],1)
            self.assertEqual(audit['stable_reference_complete_pairs'],2)
            self.assertEqual(set(response.loc[response.scope.eq('all_pairs'),'pairs_n']),{3})
            self.assertEqual(set(ablation.loc[ablation.scope.eq('stable_reference'),'pairs_n']),{2})
            self.assertTrue((root/'out/metric_response.tsv').exists())
            changes=pd.read_csv(root/'out/metric_change_distribution.tsv',sep='\t')
            self.assertEqual(len(changes),10)
            self.assertEqual(set(changes.loc[changes.scope.eq('all_pairs'),'pairs_n']),{3})
            self.assertTrue(np.allclose(changes.gain_mean,.6))
            self.assertNotIn('CSQ_before',pd.read_csv(root/'out/paired_metric_scores.tsv',sep='\t'))


if __name__=='__main__':
    unittest.main()
