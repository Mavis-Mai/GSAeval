#!/usr/bin/env python3
"""Evaluate five annotation scores using paired erroneous/correct structures.

before = erroneous structures; after = corrected structures. CSQ and label are
never used. An equal-weight mean is only a fixed ablation baseline, not CSQ.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

METRICS = ['CDS_number_score', 'Phase_match_score', 'CDS_length_score',
           'Intron_site_score', 'Total_length_score']
DISPLAY = ['CDS number', 'Phase match', 'CDS length', 'Intron site', 'Total length']
STANDARD_COLUMNS = ['BUSCO', 'GeneA', 'GeneB', 'cdsA_number', 'cdsB_number']+METRICS+['CSQ', 'label']
MISSING = {'', '.', 'na', 'nan', 'n/a', 'null'}


def read_table(path, keys):
    """Headered TSV or exactly the user-confirmed 12-column headerless layout."""
    table = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, encoding='utf-8-sig')
    table.columns = table.columns.str.strip()
    if not set(METRICS).issubset(table.columns) and len(table.columns) == 12:
        table = pd.read_csv(path, sep='\t', header=None, names=STANDARD_COLUMNS,
                            dtype=str, keep_default_na=False, encoding='utf-8-sig')
    required = set(keys+METRICS)
    if required-set(table.columns):
        raise ValueError(f'{path}: missing columns {sorted(required-set(table.columns))}')
    for key in keys:
        table[key] = table[key].str.strip()
        if table[key].str.lower().isin(MISSING).any():
            raise ValueError(f'{path}: missing pairing key {key}')
    if table.duplicated(keys).any():
        raise ValueError(f'{path}: duplicate keys {keys}; use a unique --key, do not aggregate silently')
    for metric in METRICS:
        raw = table[metric].str.strip()
        empty = raw.str.lower().isin(MISSING)
        values = pd.to_numeric(raw.mask(empty), errors='coerce')
        if ((~empty & values.isna()) | (values.notna() & ~np.isfinite(values))).any():
            raise ValueError(f'{path}: invalid number in {metric}')
        if ((values < 0) | (values > 1)).any():
            raise ValueError(f'{path}: {metric} must be in [0, 1]')
        table[metric] = values
    # Exclude CSQ/label completely from pairing, scoring, plots and reporting.
    return table.drop(columns=[c for c in ('CSQ', 'label') if c in table])


def pair_tables(before, after, keys):
    left = before.rename(columns={c: c+'_before' for c in before if c not in keys})
    right = after.rename(columns={c: c+'_after' for c in after if c not in keys})
    pairs = left.merge(right, on=keys, how='outer', validate='one_to_one', indicator=True)
    pairs['pair_status'] = pairs.pop('_merge').astype(str).map(
        {'both': 'paired', 'left_only': 'error_only', 'right_only': 'correct_only'})
    columns = [metric+suffix for metric in METRICS for suffix in ('_before', '_after')]
    pairs['complete_five_metrics'] = pairs[columns].notna().all(axis=1) & pairs['pair_status'].eq('paired')
    pairs['reference_stable'] = pd.Series(pd.NA, index=pairs.index, dtype='boolean')
    for identity, flag in [('GeneB', 'reference_stable'), ('GeneA', 'target_id_equal')]:
        if identity in keys:
            pairs[flag] = pairs['pair_status'].eq('paired')
        elif identity+'_before' in pairs and identity+'_after' in pairs:
            a, b = pairs[identity+'_before'], pairs[identity+'_after']
            known = a.notna() & b.notna() & ~a.fillna('').str.lower().isin(MISSING) & ~b.fillna('').str.lower().isin(MISSING)
            pairs[flag] = pd.Series(pd.NA, index=pairs.index, dtype='boolean')
            pairs.loc[known, flag] = a.loc[known].eq(b.loc[known])
    for metric in METRICS:
        pairs[metric+'_gain'] = pairs[metric+'_after']-pairs[metric+'_before']
    for suffix in ('_before', '_after'):
        pairs['equal_mean'+suffix] = pairs[[m+suffix for m in METRICS]].mean(
            axis=1, skipna=False).where(pairs.complete_five_metrics)
        for metric in METRICS:
            pairs['without_'+metric+'_mean'+suffix] = pairs[
                [m+suffix for m in METRICS if m != metric]].mean(
                    axis=1, skipna=False).where(pairs.complete_five_metrics)
    return pairs


def average_ranks(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind='mergesort')
    sorted_values = values[order]
    starts = np.r_[0, np.flatnonzero(np.diff(sorted_values))+1]
    counts = np.diff(np.r_[starts, len(values)])
    ranks = np.empty(len(values), dtype=float)
    ranks[order] = np.repeat(starts+(counts+1)/2, counts)
    return ranks


def auc(error_scores, correct_scores):
    """P(correct > error) + 0.5 P(tie), across all cross-group comparisons.

    Higher scores are declared better in advance; never flip a below-0.5 AUC.
    The bootstrap, unlike the point estimator, resamples whole matched pairs.
    """
    # Equal means of decimal inputs can differ by floating-point roundoff.
    # Keep mathematical ties tied, rather than creating artificial AUC gains.
    error_scores = np.round(np.asarray(error_scores, dtype=float), 12)
    correct_scores = np.round(np.asarray(correct_scores, dtype=float), 12)
    if not len(error_scores) or not len(correct_scores):
        return float('nan')
    ranks = average_ranks(np.r_[error_scores, correct_scores])
    positive_n = len(correct_scores)
    return float((ranks[len(error_scores):].sum()-positive_n*(positive_n+1)/2)/
                 (len(error_scores)*positive_n))


def roc(error_scores, correct_scores):
    values = np.round(np.r_[error_scores, correct_scores].astype(float), 12)
    positive = np.r_[np.zeros(len(error_scores)), np.ones(len(correct_scores))]
    order = np.argsort(-values, kind='mergesort')
    ends = np.r_[np.flatnonzero(np.diff(values[order])), len(values)-1]
    true_positive = np.cumsum(positive[order])[ends]
    false_positive = (ends+1)-true_positive
    return np.r_[0, false_positive/len(error_scores)], np.r_[0, true_positive/len(correct_scores)]


def point_statistics(error, correct, tolerance):
    gains = correct-error
    full_error, full_correct = error.mean(axis=1), correct.mean(axis=1)
    full_auc = auc(full_error, full_correct)
    full_gap = (full_correct-full_error).mean()
    result = dict(mean_gain=gains.mean(axis=0), median_gain=np.median(gains, axis=0),
                  improved_pct=100*(gains>tolerance).mean(axis=0),
                  unchanged_pct=100*(np.abs(gains)<=tolerance).mean(axis=0),
                  worsened_pct=100*(gains < -tolerance).mean(axis=0),
                  paired_concordance=((gains>tolerance)+.5*(np.abs(gains)<=tolerance)).mean(axis=0),
                  single_auc=np.array([auc(error[:, i], correct[:, i]) for i in range(5)]),
                  full_auc=full_auc, full_mean_gap=full_gap)
    without_auc, without_gap = [], []
    for i in range(5):
        keep = [j for j in range(5) if j != i]
        removed_error, removed_correct = error[:, keep].mean(axis=1), correct[:, keep].mean(axis=1)
        without_auc.append(auc(removed_error, removed_correct))
        without_gap.append((removed_correct-removed_error).mean())
    result['without_auc'] = np.array(without_auc)
    result['auc_loss_when_removed'] = full_auc-result['without_auc']
    result['without_mean_gap'] = np.array(without_gap)
    result['mean_gap_loss_when_removed'] = full_gap-result['without_mean_gap']
    return result


def cluster_labels(frame, key, keys):
    if not key:
        return np.arange(len(frame)).astype(str)
    if key in keys:
        values = frame[key]
    else:
        a, b = key+'_before', key+'_after'
        if a not in frame and b not in frame:
            raise ValueError(f'Cluster column not found: {key}')
        if a in frame and b in frame:
            if not frame[a].eq(frame[b]).all():
                raise ValueError(f'Cluster IDs must agree across paired versions: {key}')
        values = frame[a] if a in frame else frame[b]
    if values.isna().any() or values.astype(str).str.strip().str.lower().isin(MISSING).any():
        raise ValueError(f'Missing cluster ID in {key}')
    return values.astype(str).to_numpy()


def change_distribution_summary(frame, scope, tolerance=1e-9, large_change=.1):
    """Describe actual score/gain distributions on the shared analysis cohort."""
    rows = []
    for metric in METRICS:
        before = frame[metric+'_before'].to_numpy(float)
        after = frame[metric+'_after'].to_numpy(float)
        gain = after-before
        row = dict(scope=scope, metric=metric, pairs_n=len(frame), large_change_threshold=large_change)
        for name, values in [('error', before), ('correct', after), ('gain', gain)]:
            row.update({name+'_mean': float(values.mean()), name+'_median': float(np.median(values)),
                        name+'_q25': float(np.quantile(values, .25)), name+'_q75': float(np.quantile(values, .75)),
                        name+'_std': float(np.std(values, ddof=1))})
        row.update(mean_abs_gain=float(np.abs(gain).mean()),
                   median_abs_gain=float(np.median(np.abs(gain))),
                   p90_abs_gain=float(np.quantile(np.abs(gain), .9)),
                   gain_min=float(gain.min()), gain_max=float(gain.max()))
        conditions = dict(improved=gain>tolerance, unchanged=np.abs(gain)<=tolerance,
                          worsened=gain < -tolerance, large_improvement=gain>=large_change-tolerance,
                          large_deterioration=gain<=-large_change+tolerance,
                          large_absolute_change=np.abs(gain)>=large_change-tolerance)
        for name, mask in conditions.items():
            row[name+'_n'], row[name+'_pct'] = int(mask.sum()), float(100*mask.mean())
        rows.append(row)
    return pd.DataFrame(rows)


def ecdf_coordinates(values, lower, upper):
    """Right-continuous ECDF, preserving ties and endpoint mass exactly."""
    values = np.asarray(values, dtype=float)
    if not len(values) or not np.isfinite(values).all() or values.min()<lower or values.max()>upper:
        raise ValueError('ECDF requires nonempty finite values within the axis bounds')
    unique, counts = np.unique(values, return_counts=True)
    return np.r_[lower, unique, upper], np.r_[0, np.cumsum(counts)/len(values), 1]


def plot_distributions(frames, output, tolerance, large_change):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    palette = ['#3979A8', '#4B8D83', '#C66D3F', '#8E77A9', '#A29238']

    def save(fig, filename):
        for extension in ('png', 'pdf'):
            fig.savefig(output/(filename+'.'+extension), dpi=200, bbox_inches='tight')
        plt.close(fig)

    for scope, frame in frames.items():
        fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), sharex=True, sharey=True, constrained_layout=True)
        for ax, metric, name in zip(axes.flat, METRICS, DISPLAY):
            for suffix, label, color in [('_before', 'Before: erroneous', '#3979A8'),
                                          ('_after', 'After: corrected', '#D07644')]:
                x, y = ecdf_coordinates(frame[metric+suffix], 0, 1)
                ax.step(x, y, where='post', color=color, lw=1.9, label=label)
            ax.set(title=name, xlabel='Score', ylabel='Cumulative proportion', xlim=(0, 1), ylim=(0, 1.02))
            ax.grid(alpha=.15)
        axes.flat[-1].axis('off')
        handles, labels = axes.flat[0].get_legend_handles_labels()
        axes.flat[-1].legend(handles, labels, loc='center', frameon=False, fontsize=12)
        fig.suptitle(f'Actual score distributions: {scope}, {len(frame)} complete matched pairs')
        fig.supxlabel('Both curves use the same paired objects; a rightward shift indicates higher scores after correction', fontsize=10)
        save(fig, scope+'.before_after_score_ecdf')

        fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
        for metric, name, color in zip(METRICS, DISPLAY, palette):
            x, y = ecdf_coordinates(frame[metric+'_gain'], -1, 1)
            ax.step(x, y, where='post', color=color, lw=1.9, label=name)
        ax.axvline(0, ls='--', lw=1, color='#666666')
        ax.axvline(large_change, ls=':', lw=1, color='#999999')
        ax.set(xlabel='Paired gain = corrected score minus erroneous score', ylabel='Cumulative proportion',
               title=f'Actual within-pair changes: {scope}, n={len(frame)}', xlim=(-1, 1), ylim=(0, 1.02))
        ax.legend(loc='lower right', fontsize=9)
        ax.grid(alpha=.15)
        fig.supxlabel('Rightward = larger gains; negative gains = deterioration. Exact zeros form a jump at zero.', fontsize=9)
        save(fig, scope+'.paired_gain_ecdf')

        fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
        for metric, name, color in zip(METRICS, DISPLAY, palette):
            x, y = ecdf_coordinates(np.abs(frame[metric+'_gain']), 0, 1)
            ax.step(x, y, where='post', color=color, lw=1.9, label=name)
        ax.set(xlabel='Absolute paired change |corrected minus erroneous|', ylabel='Cumulative proportion',
               title=f'Change magnitude regardless of direction: {scope}', xlim=(0, 1), ylim=(0, 1.02))
        ax.legend(loc='lower right', fontsize=9)
        ax.grid(alpha=.15)
        fig.supxlabel('Rightward = larger changes, including deterioration; this does not by itself mean better performance.', fontsize=9)
        save(fig, scope+'.absolute_gain_ecdf')


def analyze_scope(frame, scope, keys, cluster_column, bootstrap, seed, tolerance, progress=True):
    error = frame[[m+'_before' for m in METRICS]].to_numpy(float)
    correct = frame[[m+'_after' for m in METRICS]].to_numpy(float)
    if len(error) < 2:
        raise ValueError(f'{scope}: at least two complete pairs required')
    point = point_statistics(error, correct, tolerance)
    labels = cluster_labels(frame, cluster_column, keys)
    units = pd.unique(labels)
    clusters = [np.flatnonzero(labels == label) for label in units]
    rng = np.random.default_rng(seed)
    interval_keys = ['mean_gain', 'median_gain', 'improved_pct', 'paired_concordance', 'single_auc',
                     'without_auc', 'auc_loss_when_removed', 'mean_gap_loss_when_removed', 'full_auc', 'full_mean_gap']
    samples = {key: [] for key in interval_keys}
    if len(units) >= 2:
        for repeat in range(bootstrap):
            chosen = rng.integers(0, len(units), size=len(units))
            if cluster_column:
                indices = np.concatenate([clusters[index] for index in chosen])
            else:
                indices = chosen
            statistics = point_statistics(error[indices], correct[indices], tolerance)
            for key in interval_keys:
                samples[key].append(statistics[key])
            if progress and (repeat+1) % 250 == 0:
                print(f'[{scope}] paired bootstrap {repeat+1}/{bootstrap}', file=sys.stderr, flush=True)
    intervals = {}
    for key in interval_keys:
        intervals[key] = np.quantile(samples[key], [.025, .975], axis=0) if samples[key] else None
    response, ablation = [], []
    for i, metric in enumerate(METRICS):
        row = dict(scope=scope, metric=metric, pairs_n=len(frame), bootstrap_units_n=len(units),
                   mean_error_score=error[:, i].mean(), mean_correct_score=correct[:, i].mean())
        for key in ['mean_gain', 'median_gain', 'improved_pct', 'unchanged_pct', 'worsened_pct',
                    'paired_concordance', 'single_auc']:
            row[key] = float(point[key][i])
            if key in intervals:
                for suffix, q in [('_ci_low', 0), ('_ci_high', 1)]:
                    row[key+suffix] = float(intervals[key][q, i]) if intervals[key] is not None else np.nan
        response.append(row)
        row = dict(scope=scope, removed_metric=metric, pairs_n=len(frame), bootstrap_units_n=len(units),
                   full_auc=point['full_auc'], full_mean_gap=point['full_mean_gap'])
        for key in ['without_auc', 'auc_loss_when_removed', 'without_mean_gap', 'mean_gap_loss_when_removed']:
            row[key] = float(point[key][i])
            if key in intervals:
                for suffix, q in [('_ci_low', 0), ('_ci_high', 1)]:
                    row[key+suffix] = float(intervals[key][q, i]) if intervals[key] is not None else np.nan
        low, high = row['auc_loss_when_removed_ci_low'], row['auc_loss_when_removed_ci_high']
        row['interval_interpretation'] = ('unavailable' if pd.isna(low) else 'positive' if low>0 else
                                           'negative' if high<0 else 'includes_zero')
        ablation.append(row)
    model = dict(scope=scope, pairs_n=len(frame), bootstrap_units_n=len(units), full_auc=point['full_auc'],
                 full_mean_gap=point['full_mean_gap'], combination='equal_weight_mean_of_five')
    for key in ('full_auc', 'full_mean_gap'):
        model[key+'_ci_low'] = float(intervals[key][0]) if intervals[key] is not None else np.nan
        model[key+'_ci_high'] = float(intervals[key][1]) if intervals[key] is not None else np.nan
    ranked = np.column_stack([average_ranks(correct[:, i]-error[:, i]) for i in range(5)])
    correlation = np.full((5, 5), np.nan)
    for i in range(5):
        for j in range(5):
            if np.std(ranked[:, i]) > 0 and np.std(ranked[:, j]) > 0:
                correlation[i, j] = np.corrcoef(ranked[:, i], ranked[:, j])[0, 1]
    return pd.DataFrame(response), pd.DataFrame(ablation), model, pd.DataFrame(correlation, index=METRICS, columns=METRICS)


def plot_results(response, ablation, frames, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator, FormatStrFormatter
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    colors = {'all_pairs': '#3979A8', 'stable_reference': '#D07644'}
    positions = np.arange(5)
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 5.5), sharey=True, constrained_layout=True)
    definitions = [('mean_gain', response, 'metric', 'Response to correction', 'Mean gain (correct minus erroneous)'),
                   ('single_auc', response, 'metric', 'Single-metric discrimination', 'AUC'),
                   ('auc_loss_when_removed', ablation, 'removed_metric', 'Incremental contribution', 'AUC loss when metric removed')]
    for ax, (key, table, id_column, title, xlabel) in zip(axes, definitions):
        scopes = list(table.groupby('scope', sort=False))
        offsets = [0] if len(scopes) == 1 else [-.08, .08]
        for offset, (scope, group) in zip(offsets, scopes):
            group = group.set_index(id_column).loc[METRICS]
            values = group[key].to_numpy(float)
            color = colors[scope]
            low, high = group[key+'_ci_low'].to_numpy(float), group[key+'_ci_high'].to_numpy(float)
            valid = np.isfinite(low) & np.isfinite(high)
            ax.hlines((positions+offset)[valid], low[valid], high[valid], color=color, lw=2)
            ax.plot(low[valid], (positions+offset)[valid], '|', color=color)
            ax.plot(high[valid], (positions+offset)[valid], '|', color=color)
            ax.scatter(values, positions+offset, color=color, s=48, zorder=3,
                       label=f'{scope}, n={int(group.pairs_n.iloc[0])}')
            for i, value in enumerate(values):
                annotation = f'{value:.4f}' if key == 'auc_loss_when_removed' else f'{value:.3f}'
                ax.annotate(annotation, (value, positions[i]+offset), xytext=(0, 9),
                            textcoords='offset points', ha='center', fontsize=9, color=color)
        ax.axvline(.5 if key == 'single_auc' else 0, color='#999999', lw=1, ls='--')
        ax.set(title=title, xlabel=xlabel)
        if key == 'auc_loss_when_removed':
            ax.xaxis.set_major_locator(MaxNLocator(nbins=3, min_n_ticks=3))
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        ax.grid(axis='x', alpha=.13)
        if len(scopes)>1:
            ax.legend(fontsize=8)
    axes[0].set_yticks(positions, DISPLAY)
    axes[0].set_ylim(4.6, -.6)
    if len(frames) == 1:
        count = len(next(iter(frames.values())))
        restriction = 'no GeneB restriction' if 'all_pairs' in frames else 'same GeneB'
        fig.suptitle(f'Five-metric contribution analysis: {count} matched pairs, {restriction}')
    else:
        fig.suptitle('Five-metric contribution analysis')
    fig.supxlabel('Points are estimates; horizontal intervals are paired bootstrap 95% intervals', fontsize=10)
    for extension in ('png', 'pdf'):
        fig.savefig(output/('metric_contribution.'+extension), dpi=200, bbox_inches='tight')
    plt.close(fig)
    fig, axes = plt.subplots(1, len(frames), figsize=(6*len(frames), 5.5), squeeze=False, constrained_layout=True)
    palette = ['#3979A8', '#4B8D83', '#C66D3F', '#8E77A9', '#A29238']
    for ax, (scope, frame) in zip(axes.flat, frames.items()):
        error = frame[[m+'_before' for m in METRICS]].to_numpy(float)
        correct = frame[[m+'_after' for m in METRICS]].to_numpy(float)
        for i, (name, color) in enumerate(zip(DISPLAY, palette)):
            fpr, tpr = roc(error[:, i], correct[:, i])
            ax.plot(fpr, tpr, color=color, label=f'{name}: AUC={auc(error[:, i], correct[:, i]):.3f}')
        fpr, tpr = roc(error.mean(axis=1), correct.mean(axis=1))
        ax.plot(fpr, tpr, 'k--', lw=2, label=f'Equal mean: AUC={auc(error.mean(axis=1), correct.mean(axis=1)):.3f}')
        ax.plot([0, 1], [0, 1], ':', color='#999999')
        ax.set(xlabel='False positive rate (erroneous called correct)', ylabel='True positive rate (correct called correct)',
               title=f'{scope}, n={len(frame)} pairs', xlim=(0, 1), ylim=(0, 1.02), aspect='equal')
        ax.legend(loc='lower right', fontsize=8)
    for extension in ('png', 'pdf'):
        fig.savefig(output/('error_correct_roc.'+extension), dpi=200, bbox_inches='tight')
    plt.close(fig)


def run(args):
    if args.bootstrap < 0 or args.tolerance < 0 or args.seed < 0:
        raise ValueError('bootstrap, tolerance and seed must be nonnegative')
    large_change = getattr(args, 'large_change', .1)
    if not np.isfinite(args.tolerance) or not np.isfinite(large_change) or not 0<=args.tolerance<large_change<=1:
        raise ValueError('Require 0 <= tolerance < large-change <= 1, with finite values')
    output = Path(args.output_dir).resolve()
    for source in (args.before, args.after):
        if Path(source).resolve().parent == output:
            raise ValueError('Use a dedicated output directory separate from input tables')
    before, after = (read_table(path, args.key) for path in (args.before, args.after))
    pairs = pair_tables(before, after, args.key)
    eligible = pairs.loc[pairs.complete_five_metrics].copy()
    frames = {}
    if args.scope in ('both', 'all') and len(eligible) >= 2:
        frames['all_pairs'] = eligible
    stable = eligible.loc[eligible.reference_stable.eq(True).fillna(False)].copy()
    if args.scope in ('both', 'stable') and len(stable) >= 2:
        frames['stable_reference'] = stable
    if not frames:
        raise ValueError('No requested scope has at least two complete pairs (stable scope requires matching GeneB)')
    if not args.no_plots:
        try:
            import matplotlib
            if not hasattr(matplotlib, 'use'):
                raise ImportError('Incomplete matplotlib installation')
        except ImportError as exc:
            raise ValueError('Install matplotlib for plots, or use --no-plots') from exc
    audit = dict(before_role='erroneous', after_role='corrected', before_rows=len(before), after_rows=len(after),
                 pairing_keys=args.key, pair_status_counts=pairs.pair_status.value_counts().to_dict(),
                 complete_pairs=len(eligible), excluded_incomplete_pairs=int(
                     (pairs.pair_status.eq('paired') & ~pairs.complete_five_metrics).sum()),
                 stable_reference_complete_pairs=len(stable), scopes=list(frames),
                 metrics=METRICS, ignored_columns=['CSQ', 'label'], bootstrap=args.bootstrap, seed=args.seed,
                 resampling_unit=args.cluster_column or 'paired_row', scope_requested=args.scope,
                 large_change_threshold=large_change,
                 source_files=[str(Path(args.before).resolve()), str(Path(args.after).resolve())])
    responses, ablations, models, correlations = [], [], [], {}
    for scope, frame in frames.items():
        print(f'[{scope}] {len(frame)} complete pairs; before=erroneous, after=corrected', file=sys.stderr)
        response, ablation, model, correlation = analyze_scope(
            frame, scope, args.key, args.cluster_column, args.bootstrap, args.seed, args.tolerance)
        responses.append(response); ablations.append(ablation); models.append(model); correlations[scope] = correlation
    response, ablation = pd.concat(responses, ignore_index=True), pd.concat(ablations, ignore_index=True)
    changes = pd.concat([change_distribution_summary(frame, scope, args.tolerance, large_change)
                         for scope, frame in frames.items()], ignore_index=True)
    output.mkdir(parents=True, exist_ok=True)
    for name, table in [('paired_metric_scores.tsv', pairs), ('metric_response.tsv', response),
                        ('leave_one_out.tsv', ablation), ('equal_mean_baseline.tsv', pd.DataFrame(models))]:
        table.to_csv(output/name, sep='\t', index=False, na_rep='NA')
    changes.to_csv(output/'metric_change_distribution.tsv', sep='\t', index=False, na_rep='NA')
    for scope, correlation in correlations.items():
        correlation.to_csv(output/f'{scope}.gain_spearman.tsv', sep='\t', index_label='metric', na_rep='NA')
    (output/'audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')
    notes = ['# Five-metric error/correct contribution analysis', '',
             'before = erroneous; after = corrected. Gain = after minus before. CSQ and label are excluded.',
             'All five metrics are expected to be higher for correct structures; AUC is not automatically flipped.',
             'The fixed equal-weight mean is an ablation reference, not CSQ and not a trained classifier.', '',
             '## Reading the results', '',
             '- metric_response.tsv: response magnitude, coverage, within-pair concordance and single-metric AUC.',
             '- metric_change_distribution.tsv: actual before/after quartiles and mean/median gain, absolute change, directional fractions and large-change fractions on the same complete pairs.',
             f'- Large changes use threshold {large_change:g} score units (inclusive within tolerance). These distribution summaries are descriptive; bootstrap intervals for mean/median gain are in metric_response.tsv.',
             '- before_after_score_ecdf plots two raw-score ECDFs for each metric; paired_gain_ecdf compares the five actual within-pair gain distributions; absolute_gain_ecdf compares magnitude irrespective of direction.',
             '- ECDFs use actual unrounded gains; directional percentages use the stated tolerance. Distribution separation is not a paired change statistic.',
             '- leave_one_out.tsv: full AUC minus AUC without each metric. Positive supports incremental value in this fixed combination; negative means removal improved AUC here.',
             '- mean_gap_loss_when_removed: full correct/error mean gap minus the four-metric gap. Changing the denominator affects it; it is not independent importance.',
             '- gain_spearman.tsv: correlations of within-pair gains. Constant gains give NA, not zero correlation.',
             '- Intervals are paired/cluster percentile bootstrap 95% intervals. They are marginal intervals, not multiplicity-adjusted significance tests.',
             '- An interval including zero does not demonstrate no contribution. Correlated features can substitute for each other; these are not causal contributions or percentages that sum to 100%.', '',
             '## Scope and limitations', '',
             f'- Raw rows: {len(before)} erroneous, {len(after)} corrected; complete matched pairs: {len(eligible)}; stable reference: {len(stable)}.',
             '- Every metric and every ablation uses the same complete-case pairs within its scope. Missing/unmatched scores are never imputed as zero.',
             '- all_pairs uses all complete pairs; stable_reference uses pairs with identical nonempty GeneB. Reference grouping is optional; the default all_pairs makes no GeneB restriction. Neither proves a shared source isoform.',
             '- Stable BUSCO alone does not prove identical source genes. New IDs and fusion/split require a prepared source mapping/pair key; do not silently strip ID suffixes.',
             '- Input max-CSQ winner selection cannot be undone by ignoring CSQ at analysis time. Prefer unselected score tables with fixed comparison objects.',
             '- AUC compares all cross-group score pairs; within-pair concordance instead compares each object to its own corrected version.',
             '- Bootstrap resamples both versions together. Default unit is a paired row; use --cluster-column for repeated isoforms/events with a shared, explicit cluster identifier.',
             '- One bootstrap unit cannot give an informative resampling interval; intervals are NA. Zero bootstrap iterations also gives NA intervals.',
             '- No optimized weights or thresholds are fitted. This is descriptive evaluation on the provided data, not external validation; no accuracy/specificity claims at an optimized threshold are made.',
             '- Results are specific to the error type/data provided. Evaluate additional error types separately and use biological replicates for broader claims.']
    (output/'README.md').write_text('\n'.join(notes)+'\n', encoding='utf-8')
    if not args.no_plots:
        plot_results(response, ablation, frames, output)
        plot_distributions(frames, output, args.tolerance, large_change)
    print(ablation[['scope', 'removed_metric', 'full_auc', 'without_auc', 'auc_loss_when_removed', 'interval_interpretation']].to_string(index=False))
    return response, ablation, audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', required=True, help='Erroneous structure scores')
    parser.add_argument('--after', required=True, help='Corrected structure scores')
    parser.add_argument('-o', '--output-dir', required=True)
    parser.add_argument('--key', nargs='+', default=['BUSCO'], help='Unique shared pairing columns')
    parser.add_argument('--scope', choices=['both', 'all', 'stable'], default='all')
    parser.add_argument('--bootstrap', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--cluster-column', help='Shared gene/event cluster column for correlated pairs')
    parser.add_argument('--tolerance', type=float, default=1e-9)
    parser.add_argument('--large-change', type=float, default=.1, help='Threshold for a large gain/loss in score units')
    parser.add_argument('--no-plots', action='store_true')
    args = parser.parse_args()
    try:
        run(args)
    except ValueError as exc:
        parser.error(str(exc))


if __name__ == '__main__':
    main()
