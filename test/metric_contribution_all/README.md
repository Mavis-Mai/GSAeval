# Five-metric error/correct contribution analysis

before = erroneous; after = corrected. Gain = after minus before. CSQ and label are excluded.
All five metrics are expected to be higher for correct structures; AUC is not automatically flipped.
The fixed equal-weight mean is an ablation reference, not CSQ and not a trained classifier.

## Reading the results

- metric_response.tsv: response magnitude, coverage, within-pair concordance and single-metric AUC.
- leave_one_out.tsv: full AUC minus AUC without each metric. Positive supports incremental value in this fixed combination; negative means removal improved AUC here.
- mean_gap_loss_when_removed: full correct/error mean gap minus the four-metric gap. Changing the denominator affects it; it is not independent importance.
- gain_spearman.tsv: correlations of within-pair gains. Constant gains give NA, not zero correlation.
- Intervals are paired/cluster percentile bootstrap 95% intervals. They are marginal intervals, not multiplicity-adjusted significance tests.
- An interval including zero does not demonstrate no contribution. Correlated features can substitute for each other; these are not causal contributions or percentages that sum to 100%.

## Scope and limitations

- Raw rows: 208 erroneous, 211 corrected; complete matched pairs: 206; stable reference: 148.
- Every metric and every ablation uses the same complete-case pairs within its scope. Missing/unmatched scores are never imputed as zero.
- all_pairs uses all complete pairs; stable_reference uses pairs with identical nonempty GeneB. Reference grouping is optional; the default all_pairs makes no GeneB restriction. Neither proves a shared source isoform.
- Stable BUSCO alone does not prove identical source genes. New IDs and fusion/split require a prepared source mapping/pair key; do not silently strip ID suffixes.
- Input max-CSQ winner selection cannot be undone by ignoring CSQ at analysis time. Prefer unselected score tables with fixed comparison objects.
- AUC compares all cross-group score pairs; within-pair concordance instead compares each object to its own corrected version.
- Bootstrap resamples both versions together. Default unit is a paired row; use --cluster-column for repeated isoforms/events with a shared, explicit cluster identifier.
- One bootstrap unit cannot give an informative resampling interval; intervals are NA. Zero bootstrap iterations also gives NA intervals.
- No optimized weights or thresholds are fitted. This is descriptive evaluation on the provided data, not external validation; no accuracy/specificity claims at an optimized threshold are made.
- Results are specific to the error type/data provided. Evaluate additional error types separately and use biological replicates for broader claims.
