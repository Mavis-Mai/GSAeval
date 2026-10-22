# CSQ 四个小指标的分布与配对变化

脚本：`compare_CSQ_components.py`。依赖 pandas、numpy；绘图另需 matplotlib。

## 两张表比较

```sh
python other_scripts/compare_CSQ_components.py --before original_max_CSQ.txt --after CSQ2db-witherro_max_CSQ.txt -o CSQ_comparison
```

默认按 BUSCO 唯一配对；重复 BUSCO 会报错，不静默取第一行。未汇总表可用 `--key BUSCO GeneA` 等指定唯一键，修改后的 ID 需结合来源映射。

接受带表头 TSV，或用户确认的 12 列无表头格式（第一条数据不会丢失）：

```
BUSCO GeneA GeneB cdsA_number cdsB_number CDS_number_score Phase_match_score CDS_length_score Intron_site_score Total_length_score CSQ label
```

当前只有一张表时：

```sh
python other_scripts/compare_CSQ_components.py --after test/CSQ2db-witherro_addexon.CSQ.tsv -o test/CSQ_analysis
```

`--no-plots` 只输出统计表。脚本读取制表符，不接受以上示意中的空格作为实际分隔符。

## 图表

按用户偏好，默认只输出折线图（PNG 和 PDF），不画柱形分布图：

- score_mirrored_lines：四个指标的上下镜像频率折线。原始表在上方，错误后在下方；横轴得分 [0,1]，统一按 0.05 分箱，以箱中心连接频率点，不进行 KDE 平滑。纵轴显示频率绝对值，上下仅区分两组，并不是负得分。双表模式仅使用双方分数均有效的配对记录，两组样本数相同。四幅图采用相同坐标尺度。
- paired_decrease_ecdf：四个指标的配对下降量累积分布，下降量 = 原始 − 错误后。正值表示下降，负值表示上升。相同百分位下越靠右表示下降越大，曲线交叉时不代表某个指标在所有对象上都更敏感。

单表模式仅画 score_mirrored_lines，只有一条分布线，不生成虚假的配对差值图。

## 怎样衡量“变化剧烈”

所有分数都在 [0,1]，可直接比较绝对分值变化，不需要用容易在基线为零时失效的相对百分比。

- mean_abs_delta：平均 |Δ|，作为变化幅度排名的主要依据。
- p90_abs_delta：|Δ| 的第 90 百分位，检查一部分对象的大变化。
- mean_delta / median_delta：变化方向，正负可抵消，不能单独作为剧烈程度。
- changed_pct / decreased_pct / increased_pct / unchanged_pct：受影响范围和方向。
- large_decrease_pct：下降至少 0.1 分的比例；阈值可用 `--large-change 0.1` 修改。
- paired_valid_n / paired_missing_n：每个指标使用的配对样本数；缺失不填零。

## 输出与配对限制

`distribution_summary.tsv` 保存两组各自的均值、中位数、标准差、四分位数、0 分和 1 分比例。
`paired_scores.tsv` 保留完整配对和未配对行；`metric_change_summary.tsv` 保存变化统计。
`pairing_audit.json` 保存配对数、对象切换等；`analysis.md` 简述结果与限制。

max_CSQ 在项目中是按 BUSCO 选择综合 CSQ 最高的一条记录，而非四个分项各自最大。两个表可能选择不同 GeneA / GeneB。因此默认分析是 BUSCO 层面的最佳结果响应，不一定是同一转录本对相同参考的变化。脚本标记 GeneA / GeneB 是否切换，并输出 GeneB 未切换子集的 `stable_reference_summary.tsv`。GeneB 相同仍不能保证 CDS 动态匹配路径完全相同。

对于 ID 重写、融合和拆分，建议使用错误报告做来源映射：

```sh
python other_scripts/compare_CSQ_components.py --before original_max_CSQ.txt --after CSQ2db-witherro_max_CSQ.txt --error-report modified.gtf.error_report.tsv --mapping-level transcript -o CSQ_source_comparison
```

默认 GeneA 按 transcript_id/source_transcript_id 映射；若表中 GeneA 实际是 gene_id，可选择 `--mapping-level gene`，但歧义映射会显式报错。来源映射仍要求 BUSCO 相同（两表均有此列时）。若 `max_CSQ` 已经筛掉某个来源/产物，不能从聚合表恢复；更严格实验应保留未汇总结果，并固定同一 GeneB。

来源映射模式额外输出 `by_error_type.tsv`，只统计 error_status=applied 的已配对记录；一对多事件按边统计，不视为独立实验重复。融合/拆分多产物可能得到更高权重，应按事件另做重复采样评估，当前不提供错误的独立样本显著性检验。

原表有、错误后缺失的 BUSCO/ORF 不计零分，单独保留为 before_only。仅比较幸存对象可能低估损伤。指标对错误的响应幅度不等同于准确率、特异性或蛋白相似性。

测试：

```sh
python -m unittest discover -s other_scripts -p test_compare_CSQ_components.py
```
