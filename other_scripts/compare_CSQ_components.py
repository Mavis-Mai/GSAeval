#!/usr/bin/env python3
"""Compare four CSQ component scores in two tab-separated result tables."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

METRICS = ['CDS_number_score', 'Phase_match_score', 'CDS_length_score', 'Intron_site_score']
SHORT = ['CDS number', 'Phase match', 'CDS length', 'Intron site']
STANDARD_COLUMNS = ['BUSCO', 'GeneA', 'GeneB', 'cdsA_number', 'cdsB_number'] + METRICS + ['Total_length_score', 'CSQ', 'label']


def read_scores(path, keys):
    table = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    table.columns = table.columns.str.strip()
    required = list(dict.fromkeys(keys+METRICS))
    # Recognize the user-confirmed 12-column headerless CSQ layout only.
    if not set(METRICS).issubset(table.columns) and len(table.columns) == 12:
        table = pd.read_csv(path, sep='\t', names=STANDARD_COLUMNS, header=None,
                            dtype=str, keep_default_na=False)
    missing = set(required)-set(table.columns)
    if missing:
        raise ValueError(f'{path}: missing columns {sorted(missing)}; expected tab-separated input')
    for column in keys:
        table[column] = table[column].str.strip()
        if table[column].eq('').any():
            raise ValueError(f'{path}: blank key {column}')
    if table.duplicated(keys).any():
        raise ValueError(f'{path}: duplicate keys {keys}; choose a unique --key (no automatic aggregation)')
    for metric in METRICS:
        raw = table[metric].str.strip()
        empty = raw.str.lower().isin(['', 'na', 'nan', 'n/a', 'null', '.'])
        numeric = pd.to_numeric(raw.mask(empty), errors='coerce')
        if ((~empty & numeric.isna()) | (~np.isfinite(numeric) & numeric.notna())).any():
            raise ValueError(f'{path}: invalid numeric value in {metric}')
        if ((numeric < 0) | (numeric > 1)).any():
            raise ValueError(f'{path}: {metric} outside expected [0,1] range')
        table[metric] = numeric
    return table


def pair_tables(before, after, keys, mapping=None, mapping_level='transcript'):
    """Outer pairing keeps disappeared and new observations explicit.

    Optional mapping compares source/product edges, constrained by BUSCO when
    present in both tables. Direct pairing compares max-CSQ winners per key.
    """
    if mapping is None:
        left = before.rename(columns={c: c+'_before' for c in before if c not in keys})
        right = after.rename(columns={c: c+'_after' for c in after if c not in keys})
        pairs = left.merge(right, on=keys, how='outer', validate='one_to_one', indicator=True)
    else:
        if 'GeneA' not in before or 'GeneA' not in after:
            raise ValueError('Mapping requires GeneA in both score tables')
        current = 'transcript_id' if mapping_level == 'transcript' else 'gene_id'
        source = 'source_transcript_id' if mapping_level == 'transcript' else 'source_gene_id'
        required = {current, source, 'error_type', 'error_status'}
        if not required.issubset(mapping):
            raise ValueError(f'Mapping report missing {sorted(required-set(mapping))}')
        available = [current, source]+[c for c in ['error_type', 'error_status', 'event_id',
                    'protein_identity_loss_pct', 'coding_impact'] if c in mapping]
        links = mapping[available].copy().drop_duplicates()
        links = links.rename(columns={current: 'GeneA_after', source: 'GeneA_before'})
        if links.duplicated(['GeneA_before', 'GeneA_after']).any():
            raise ValueError('Ambiguous mapping edges; use --mapping-level transcript or resolve duplicate annotations')
        left = before.add_suffix('_before')
        right = after.add_suffix('_after')
        candidate = right.merge(links, on='GeneA_after', how='left', indicator='_mapping')
        on = ['GeneA_before']
        if 'BUSCO' in before and 'BUSCO' in after:
            candidate['BUSCO_before'] = candidate['BUSCO_after']
            on.append('BUSCO_before')
        if left.duplicated(on).any():
            raise ValueError('Source table is not unique on mapping key: '+str(on))
        pairs = candidate.merge(left, on=on, how='outer', indicator=True, validate='many_to_one')
    pairs['pair_status'] = pairs['_merge'].astype(str).map(
        {'both': 'paired', 'left_only': 'before_only', 'right_only': 'after_only'} if mapping is None else
        {'both': 'paired', 'left_only': 'after_only', 'right_only': 'before_only'})
    pairs = pairs.drop(columns='_merge')
    if 'GeneB_before' in pairs and 'GeneB_after' in pairs:
        known = pairs['GeneB_before'].notna() & pairs['GeneB_after'].notna()
        pairs['reference_changed'] = pd.Series(pd.NA, index=pairs.index, dtype='boolean')
        pairs.loc[known, 'reference_changed'] = pairs.loc[known, 'GeneB_before'].ne(pairs.loc[known, 'GeneB_after'])
    if 'GeneA_before' in pairs and 'GeneA_after' in pairs:
        known = pairs['GeneA_before'].notna() & pairs['GeneA_after'].notna()
        pairs['target_id_changed'] = pd.Series(pd.NA, index=pairs.index, dtype='boolean')
        pairs.loc[known, 'target_id_changed'] = pairs.loc[known, 'GeneA_before'].ne(pairs.loc[known, 'GeneA_after'])
    for metric in METRICS:
        pairs[metric+'_delta'] = pairs[metric+'_after']-pairs[metric+'_before']
    return pairs


def summarize(pairs, tolerance=1e-9, large_change=0.1):
    rows = []
    for metric in METRICS:
        delta = pairs.loc[pairs['pair_status'].eq('paired'), metric+'_delta'].dropna()
        count = len(delta)
        row = dict(metric=metric, paired_valid_n=count,
                   paired_missing_n=int(pairs['pair_status'].eq('paired').sum()-count))
        row.update(mean_delta=delta.mean(), median_delta=delta.median(),
                   mean_abs_delta=delta.abs().mean(), median_abs_delta=delta.abs().median(),
                   p90_abs_delta=delta.abs().quantile(.9) if count else np.nan,
                   rmse=float(np.sqrt(np.mean(np.square(delta)))) if count else np.nan)
        for label, condition in [('changed', delta.abs()>tolerance), ('decreased', delta < -tolerance),
                                 ('increased', delta > tolerance), ('unchanged', delta.abs()<=tolerance),
                                 ('large_change', delta.abs()>=large_change-tolerance),
                                 ('large_decrease', delta<=-large_change+tolerance)]:
            row[label+'_n'] = int(condition.sum())
            row[label+'_pct'] = 100*condition.mean() if count else np.nan
        rows.append(row)
    result = pd.DataFrame(rows)
    result['magnitude_rank'] = result['mean_abs_delta'].rank(ascending=False, method='min')
    return result


def distribution_summary(before, after, labels):
    rows = []
    for label, frame in zip(labels, (before, after)):
        for metric in METRICS:
            values = frame[metric].dropna()
            rows.append(dict(dataset=label, metric=metric, total_n=len(frame), valid_n=len(values),
                             missing_n=int(frame[metric].isna().sum()), mean=values.mean(),
                             median=values.median(), std=values.std(), q25=values.quantile(.25),
                             q75=values.quantile(.75), zero_pct=100*values.eq(0).mean(),
                             one_pct=100*values.eq(1).mean()))
    return pd.DataFrame(rows)


def plot_results(before, after, pairs, summary, output, labels, single=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.fonttype': 'none', 'pdf.fonttype': 42})
    colors = ['#3979A8', '#D07644']

    def save(figure, name):
        figure.savefig(output/(name+'.png'), dpi=220, bbox_inches='tight')
        figure.savefig(output/(name+'.pdf'), bbox_inches='tight')
        plt.close(figure)

    from matplotlib.ticker import FuncFormatter
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    score_bins = np.linspace(0, 1, 21)
    centers = (score_bins[:-1]+score_bins[1:])/2
    for ax, metric, title in zip(axes.flat, METRICS, SHORT):
        if single:
            datasets = [(labels[1], after[metric].dropna(), colors[1], 1)]
        else:
            valid = pairs.loc[pairs['pair_status'].eq('paired'), [metric+'_before', metric+'_after']].dropna()
            datasets = [(labels[0], valid[metric+'_before'], colors[0], 1),
                        (labels[1], valid[metric+'_after'], colors[1], -1)]
        for label, values, color, sign in datasets:
            if len(values):
                percentages = np.histogram(values, bins=score_bins)[0]*100/len(values)
                ax.plot(centers, sign*percentages, color=color, lw=2,
                        marker='o', markersize=2.6, label=f'{label} (n={len(values)})')
                ax.fill_between(centers, 0, sign*percentages, color=color, alpha=.08)
        ax.axhline(0, color='#666666', lw=.8)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, pos: f'{abs(value):g}'))
        ax.set(title=title, xlabel='Score', ylabel='BUSCOs per score bin (%)', xlim=(0, 1),
               ylim=(0, 100) if single else (-100, 100))
        if not single and len(valid):
            decrease = valid[metric+'_before']-valid[metric+'_after']
            ax.text(.03, .94, f'Mean decrease: {decrease.mean():.3f}',
                    transform=ax.transAxes, va='top', fontsize=10)
        ax.legend(loc='lower left', fontsize=9)
        ax.grid(axis='y', alpha=.13)
    fig.suptitle('Paired score distributions: original above, exon addition below' if not single else
                 'CSQ component score distributions')
    fig.supxlabel('Equal 0.05-wide bins; vertical sign separates datasets, not positive/negative scores', fontsize=10)
    save(fig, 'score_mirrored_lines')
    if single:
        return

    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    curve_colors = ['#3979A8', '#4B8D83', '#C66D3F', '#8E77A9']
    for metric, title, color in zip(METRICS, SHORT, curve_colors):
        decrease = -pairs.loc[pairs['pair_status'].eq('paired'), metric+'_delta'].dropna().to_numpy()
        if len(decrease):
            values = np.sort(decrease)
            x = np.r_[-1, values, 1]
            y = np.r_[0, np.arange(1, len(values)+1)/len(values), 1]
            ax.step(x, y, where='post', color=color, lw=2,
                    label=f'{title}: mean decrease {decrease.mean():.3f}')
    ax.axvline(0, color='#555555', ls='--', lw=1)
    ax.set(xlabel='Paired decrease = original score minus score after exon addition',
           ylabel='Cumulative proportion', xlim=(-.15, 1), ylim=(0, 1.02),
           title='Larger decreases shift a curve to the right')
    # Include increases beyond the usual small negative range when present.
    decreases = -pairs[[metric+'_delta' for metric in METRICS]].to_numpy(dtype=float)
    finite = decreases[np.isfinite(decreases)]
    if len(finite) and finite.min() < -.15:
        ax.set_xlim(max(-1, finite.min()-.03), 1)
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(alpha=.15)
    save(fig, 'paired_decrease_ecdf')


def run(args):
    if not 0 <= args.tolerance < args.large_change <= 1:
        raise ValueError('Require 0 <= tolerance < large-change <= 1')
    if not args.before:
        if args.error_report:
            raise ValueError('--error-report requires both --before and --after')
        after = read_scores(args.after, args.key)
        output = Path(args.output_dir)
        if Path(args.after).resolve() == (output/'distribution_summary.tsv').resolve():
            raise ValueError('Output would overwrite input')
        if not args.no_plots:
            try:
                import matplotlib
            except ImportError as exc:
                raise RuntimeError('Plotting requires matplotlib; install it or use --no-plots') from exc
        output.mkdir(parents=True, exist_ok=True)
        summary = distribution_summary(after.iloc[:0], after, (args.before_label, args.after_label))
        summary = summary.loc[summary['dataset'].eq(args.after_label)]
        summary.to_csv(output/'distribution_summary.tsv', sep='\t', index=False, na_rep='NA')
        (output/'analysis.md').write_text(
            f'# Single-table CSQ analysis\n\nRecords: {len(after)}. No before table supplied.\n'
            'Scores and distance from 1 describe the observed distribution; they are not paired changes.\n', encoding='utf-8')
        if not args.no_plots:
            empty = after.iloc[:0]
            pairs = pair_tables(empty, empty, args.key)
            plot_results(empty, after, pairs, summarize(pairs), output, (args.before_label, args.after_label), single=True)
        print(summary.to_string(index=False))
        return after, summary, {'mode': 'single_table', 'after_rows': len(after)}
    before, after = (read_scores(path, args.key) for path in (args.before, args.after))
    mapping = pd.read_csv(args.error_report, sep='\t', dtype=str, keep_default_na=False) if args.error_report else None
    pairs = pair_tables(before, after, args.key, mapping, args.mapping_level)
    output = Path(args.output_dir)
    # Source files must never be replaced by an output of this script.
    names = ['paired_scores.tsv', 'metric_change_summary.tsv', 'paired_decrease_summary.tsv', 'distribution_summary.tsv',
             'stable_reference_summary.tsv', 'by_error_type.tsv', 'pairing_audit.json', 'analysis.md']
    for source in [args.before, args.after, args.error_report]:
        if source and Path(source).resolve() in [(output/name).resolve() for name in names]:
            raise ValueError('An output path would overwrite an input file')
    if not args.no_plots:
        try:
            import matplotlib
        except ImportError as exc:
            raise RuntimeError('Plotting requires matplotlib: pip install matplotlib; or use --no-plots') from exc
    output.mkdir(parents=True, exist_ok=True)
    summary = summarize(pairs, args.tolerance, args.large_change)
    pairs.to_csv(output/'paired_scores.tsv', sep='\t', index=False, na_rep='NA')
    summary.to_csv(output/'metric_change_summary.tsv', sep='\t', index=False, na_rep='NA')
    decreases = summary.copy()
    decreases['mean_decrease'], decreases['median_decrease'] = -decreases['mean_delta'], -decreases['median_delta']
    for index, row in decreases.iterrows():
        metric = row['metric']
        valid = pairs.loc[pairs['pair_status'].eq('paired'), [metric+'_before', metric+'_after']].dropna()
        decreases.loc[index, 'mean_before'] = valid.iloc[:, 0].mean()
        decreases.loc[index, 'mean_after'] = valid.iloc[:, 1].mean()
    decreases.sort_values('mean_decrease', ascending=False)[
        ['metric', 'paired_valid_n', 'mean_before', 'mean_after', 'mean_decrease', 'median_decrease',
         'decreased_pct', 'large_decrease_pct']].to_csv(
             output/'paired_decrease_summary.tsv', sep='\t', index=False, na_rep='NA')
    distribution_summary(before, after, (args.before_label, args.after_label)).to_csv(
        output/'distribution_summary.tsv', sep='\t', index=False, na_rep='NA')
    audit = dict(before_rows=len(before), after_rows=len(after), pairing_mode='source_mapping' if mapping is not None else 'direct_keys',
                 keys=args.key, mapping_level=args.mapping_level if mapping is not None else None,
                 pair_status_counts=pairs['pair_status'].value_counts().to_dict(),
                 large_change=args.large_change, tolerance=args.tolerance)
    if 'reference_changed' in pairs:
        stable = pairs.loc[pairs['reference_changed'].eq(False).fillna(False)]
        summarize(stable, args.tolerance, args.large_change).to_csv(
            output/'stable_reference_summary.tsv', sep='\t', index=False, na_rep='NA')
        audit['reference_changed_pairs'] = int(pairs['reference_changed'].eq(True).sum())
        audit['stable_reference_pairs'] = int(pairs['reference_changed'].eq(False).sum())
    if 'target_id_changed' in pairs:
        audit['target_id_changed_pairs'] = int(pairs['target_id_changed'].eq(True).sum())
    if mapping is not None:
        audit['after_rows_without_mapping'] = int(pairs.get('_mapping', pd.Series(dtype=str)).eq('left_only').sum())
        groups = []
        eligible = pairs.loc[pairs['pair_status'].eq('paired') & pairs['error_status'].eq('applied')]
        for error, group in eligible.groupby('error_type'):
            frame = summarize(group, args.tolerance, args.large_change)
            frame.insert(0, 'error_type', error)
            groups.append(frame)
        by_error = pd.concat(groups, ignore_index=True) if groups else pd.DataFrame(columns=['error_type']+summary.columns.tolist())
        by_error.to_csv(output/'by_error_type.tsv', sep='\t', index=False, na_rep='NA')
    (output/'pairing_audit.json').write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding='utf-8')
    lines = ['# CSQ component comparison', '', f'Before: {args.before_label}; after: {args.after_label}.',
             f'Pairing mode: {audit["pairing_mode"]}. Raw rows: {len(before)} / {len(after)}.',
             '', 'Delta = after minus before. Scores have the same [0,1] scale. Rank by mean absolute delta;',
             'also inspect downward response, unchanged fraction, and the 90th percentile of absolute changes.',
             f'Large change threshold: {args.large_change:g} score units. Missing values are never filled with zero.', '',
             '| Metric | Valid pairs | Mean delta | Mean absolute delta | Decreased (%) | Large decrease (%) |',
             '|---|---:|---:|---:|---:|---:|']
    for row in summary.sort_values('magnitude_rank').itertuples():
        lines.append(f'| {row.metric} | {row.paired_valid_n} | {row.mean_delta:.4f} | {row.mean_abs_delta:.4f} | {row.decreased_pct:.1f} | {row.large_decrease_pct:.1f} |')
    lines += ['', '## Interpretation limits', '',
              '- Distribution plots use all valid rows in each table; paired summaries use only matched rows with both scores available.',
              '- max_CSQ tables select the highest TOTAL CSQ per BUSCO, not the maximum of each component. Winners can change.',
              '- Check target/reference ID switches and stable_reference_summary.tsv. Equal GeneB does not guarantee an identical CDS matching path.',
              '- Missing BUSCOs, lost ORFs and absent score rows are reported as unmatched, not assigned zero. Survivor-only summaries may underestimate damage.',
              '- Source mapping gives source/product comparisons; fusion/split edges are dependent observations and may overweight multi-product events.',
              '- IDs alone cannot identify unchanged biological objects after renaming. Use the source mapping and unsummarized score tables for controlled comparisons.',
              '- A large response measures sensitivity in this dataset, not annotation accuracy or specificity. Replicates/unchanged controls are needed for robustness.']
    (output/'analysis.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    if not args.no_plots:
        plot_results(before, after, pairs, summary, output, (args.before_label, args.after_label))
    print(summary[['metric', 'paired_valid_n', 'mean_delta', 'mean_abs_delta', 'decreased_pct', 'large_decrease_pct']].to_string(index=False))
    print('Output:', output.resolve())
    return pairs, summary, audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', help='Original scores; omit for single-table distribution summary')
    parser.add_argument('--after', required=True)
    parser.add_argument('-o', '--output-dir', required=True)
    parser.add_argument('--key', nargs='+', default=['BUSCO'])
    parser.add_argument('--before-label', default='Original')
    parser.add_argument('--after-label', default='With errors')
    parser.add_argument('--error-report', help='Optional error simulator TSV for source/product pairing')
    parser.add_argument('--mapping-level', choices=['transcript', 'gene'], default='transcript')
    parser.add_argument('--tolerance', type=float, default=1e-9)
    parser.add_argument('--large-change', type=float, default=.1)
    parser.add_argument('--no-plots', action='store_true')
    args = parser.parse_args()
    try:
        run(args)
    except (ValueError, RuntimeError) as exc:
        parser.error(str(exc))


if __name__ == '__main__':
    main()
