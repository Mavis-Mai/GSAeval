# 五项指标的纠错响应与去除贡献分析

脚本：`detect_metric_contribution.py`。独立运行，不依赖此前的 compare_CSQ_components.py。

**方向固定：before = 错误结构，after = 修正后的正确结构。** 不根据文件名自动判断。CSQ 总分及 label 在读取后直接移除，不参与分析、标签定义、权重、筛选或阈值选择。

依赖：`numpy pandas`；绘图另需 `matplotlib`。

```sh
python -m pip install numpy pandas matplotlib

python other_scripts/detect_metric_contribution.py \
  --before test/CSQ2db-witherro_addexon.CSQ.tsv \
  --after test/CSQ2db-origin_addexon.CSQ.tsv \
  -o test/metric_contribution \
  --bootstrap 1000 --seed 42 --scope all
```

按用户最新要求，默认不限制 GeneB，使用全部有效配对；正式结果已生成于 test/metric_contribution_all，使用 2000 次重采样。

## 输入与配对

分析五项：CDS_number_score、Phase_match_score、CDS_length_score、Intron_site_score、Total_length_score。均要求数值在 [0,1]，预先规定越高越好，不根据结果翻转方向。

支持带表头 TSV，或以下已确认的 12 列无表头 TSV：

```
BUSCO GeneA GeneB cdsA_number cdsB_number CDS_number_score Phase_match_score CDS_length_score Intron_site_score Total_length_score CSQ label
```

默认按 BUSCO 唯一配对。重复键会报错，不静默取最高分或第一行。可用 `--key BUSCO GeneA` 等指定多列唯一键，但新旧 ID 不同的情况必须先提供共享的来源配对键；脚本不擅自去掉 ID 后缀，不猜测融合/拆分来源。

默认 `--scope all` 仅分析 all_pairs，不考虑 GeneB 是否相同。可选 `--scope both` 同时分析：

- all_pairs：五项都完整的全部配对。
- stable_reference：同一 GeneB、且五项完整的配对。GeneB 空白或缺失不算相同。

可选择 `--scope all` 或 `--scope stable`。每个分析范围至少需要 2 对。某个范围不足 2 对时不生成其结果；若所有请求范围都不足则报错。

所有指标、五项组合和五次去除，在每个范围内使用**完全相同的完整配对集合**。缺失不填零，也不会在去掉某项后悄悄纳入原本因该项缺失而被排除的对象。

BUSCO 相同或 GeneB 相同，并不能证明是同一个原始异构体。若输入已按最大 CSQ 选择了胜出记录，这个筛选影响无法通过忽略 CSQ 列消除；严格实验应从未筛选结果中固定相同来源及参考对象。

## 三类问题

### 1. 指标自己能响应纠错吗？

对每个配对计算 gain = 正确分数 − 错误分数。

- mean_gain / median_gain：平均和中位提升。
- improved_pct / unchanged_pct / worsened_pct：提升、不变、降低的对象比例，按 `--tolerance` 判断。
- paired_concordance：同一对象正确版本得分较高计 1，平分计 0.5，较低计 0，再求平均。它与下述跨组 AUC 不同。
- single_auc：仅凭这一项，正确结构是否能排在错误结构前面。正类为 after 的正确结构；AUC 比较所有跨组正确/错误组合，平分计半。低于 0.5 不会自动翻转。

### 2. 去掉这项，会损失多少区分能力？

先建立固定的五项等权平均 `S5`，然后每次去掉一项，对剩余四项求平均 `S_without`。它们只是消融实验基线，**不是 CSQ，也不是训练得到的新评分模型**。

- auc_loss_when_removed = AUC(S5) − AUC(S_without)。正值说明该指标在当前固定组合中有补充价值；接近 0 可能是冗余、饱和或证据不足；负值说明去掉它后当前数据上的 AUC 更高。
- mean_gap_loss_when_removed = 五项平均的“正确−错误”均值差 − 四项平均的对应均值差。

注意，均值分差的去除量受分母从 5 改为 4 影响。若第 i 项平均提升为 d_i，五项平均提升为 D，则该去除量严格等于 `(d_i − D)/4`。它主要表示相对五项平均响应的强弱，不能独立代表信息贡献。

不输出“贡献百分比”，各项去除量不要求相加为 100%。若多个指标都能很好地区分、又高度重复，可能出现各项单独 AUC 很高、去除损失几乎全为零；这不是程序错误。

### 3. 指标的信息是否重复？

对五项**配对提升量**计算 Spearman 相关，而不是对两组原始分数混在一起算相关。高度同步的提升可能提示冗余；相关不能证明因果。常量提升的相关性未定义，输出 NA。

## 不确定范围

默认 1000 次 95% 百分位 bootstrap：每次同时抽取一条配对的错误和正确版本，并对所有指标/所有去除组合使用相同抽样索引。

一个基因/事件有多个异构体时，可在两张输入表中准备相同的聚类列，例如 `event_id` 或 `source_gene_id`：

```sh
python other_scripts/detect_metric_contribution.py --before error.tsv --after correct.tsv \
  -o contribution --cluster-column source_gene_id --bootstrap 2000
```

此时按聚类整体抽样，保留组内全部配对。聚类列可只存在于一张表；若两张表都有，配对值必须一致。它必须真实表示共享单位，不能用已经重命名的新 gene_id 猜测来源。

仅有一个聚类单位不能得到有效区间，区间输出 NA。`--bootstrap 0` 可只计算点估计。少于数百次的重采样适合快速检查，不建议解读其区间稳定性。

`interval_interpretation` 中 positive / negative / includes_zero 仅描述单项区间是否跨零；不是经过多重比较校正的显著性检验。区间跨零不等于证明“没有贡献”。

AUC 和 ROC 内部将评分保留 12 位小数，以避免等权求和的浮点误差把数学上相等的评分当成不同值。配对提升的零判定另由 tolerance 控制。

## 输出

| 文件 | 用途 |
| --- | --- |
| metric_response.tsv | 单项响应、单指标 AUC、配对一致率及其区间 |
| leave_one_out.tsv | 逐项去除后的 AUC/均值分差及损失区间 |
| equal_mean_baseline.tsv | 五项等权组合的基准 AUC、均值分差及区间 |
| paired_metric_scores.tsv | 完整配对审计、五项提升、组合与去除后逐对象均值；不含 CSQ/label |
| all_pairs.gain_spearman.tsv | 全配对提升量相关性矩阵 |
| stable_reference.gain_spearman.tsv | 参考未变子集的提升量相关性矩阵（该范围可分析时） |
| metric_contribution.png/.pdf | 三列点区间图：平均提升、单项 AUC、去除后 AUC 损失；带 95% 区间，不连接指标点 |
| error_correct_roc.png/.pdf | 五项单指标与等权组合的 ROC 曲线 |
| audit.json / README.md | 样本范围、排除数、随机种子、方法与解释限制 |

不画柱形分布图。`--no-plots` 可只生成表格。

这些结果只解释提供的数据集和错误类型。没有拟合权重、分类阈值或预测器，也没有做外部验证；不能直接把结果当作新数据上的错误检测准确率。

## 验证

```sh
python -m unittest discover -s other_scripts -p test_detect_metric_contribution.py
```

测试包含已知 AUC、平分、方向不翻转、唯一信号、完全冗余、负贡献、完整配对范围、聚类、随机可复现性、缺失值和无表头输入。开发时另以独立人工数据完成绘图检查。

## 新增：得分与真实配对变化的分布

在原有贡献分析之外，运行同一条命令会自动新增下列文件（默认前缀为 all_pairs；可选范围分别单独输出）：

| 文件 | 看什么 |
| --- | --- |
| all_pairs.before_after_score_ecdf.png/.pdf | 五幅小图，每项两条累积分布曲线：before 错误、after 正确。横轴得分 0–1，正确曲线右移表示得分整体提高 |
| all_pairs.paired_gain_ecdf.png/.pdf | 五项放在一张图，展示每个对象的真实 Δ=after−before；负值是反向下降，0 是完全不变，正值是提升 |
| all_pairs.absolute_gain_ecdf.png/.pdf | 五项绝对变化 |Δ| 的分布，比较变化幅度，不将变化大直接解释为表现好 |
| metric_change_distribution.tsv | 两组实际分数的均值/中位数/四分位数/标准差，以及配对变化的统计 |

分布图使用与贡献分析一致的完整配对集合，不混用两组全部未配对样本。默认 all_pairs 不要求 GeneB 一致。分数缺失不填零，缺失配对对象不被隐式当成“没有变化”。

新增统计字段：

- error_mean/median/q25/q75/std、correct_mean/median/q25/q75/std：错误、正确结构的实际得分分布。
- gain_mean/median/q25/q75/std、gain_min/max：实际配对 Δ 的分布。
- mean_abs_gain、median_abs_gain、p90_abs_gain：变化幅度，用于识别正负抵消和一部分对象的大变化。
- improved/unchanged/worsened_n 及 _pct：提升、不变、下降的对象数与比例。
- large_improvement/large_deterioration/large_absolute_change_n 及 _pct：大幅改善、恶化、绝对变化的对象数与比例。

大幅变化阈值默认 0.1 分，可调整：

```sh
python other_scripts/detect_metric_contribution.py   --before error.tsv --after correct.tsv -o contribution   --scope all --large-change 0.1 --bootstrap 1000 --seed 42
```

阈值按 >= 判断，并考虑 tolerance 的浮点边界；要求 0 <= tolerance < large-change <= 1。统计表记录阈值。`--no-plots` 仍输出新分布统计表。

ECDF 的纵轴表示“不超过该横轴数值”的比例。两条原始得分分布分离，并不等于对象逐一发生相同幅度的变化，因此必须结合配对 Δ 曲线。相同百分位下 Δ 曲线更靠右通常表示提升更大；曲线交叉时不能断言一个指标对每个对象都更敏感。

ECDF 保留实际 Δ，不使用 KDE 平滑，不按 tolerance 将小变化强制归零。因此零点跳跃表示完全不变的比例，而表中的 unchanged_pct 按 |Δ|<=tolerance 统计。新的变化分布表为描述性统计；平均/中位提升的 bootstrap 区间仍在 metric_response.tsv。

目前仍需共享且唯一的配对键。融合/拆分产生重复 BUSCO 时，应先准备真实来源映射后的唯一配对记录；不能用删除重复行的方法制造配对。
