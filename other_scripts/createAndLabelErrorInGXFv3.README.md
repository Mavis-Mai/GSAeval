# 人工注释错误与蛋白影响标签

脚本：`createAndLabelErrorInGXFv3.py`。仅接受 GTF，GFF3 需先转换。
必须输入参考基因组 DNA FASTA（不是蛋白 FASTA），序列名须与 GTF 第一列一致。基础运行无需第三方 Python 依赖；建议安装 NumPy 以加速蛋白比较。

```sh
# 八种错误一起模拟
python other_scripts/createAndLabelErrorInGXFv3.py -i reference.gtf -g genome.fa -o perturbed.gtf -e 0.2 --seed 42

# 只模拟基因融合与基因拆分
python other_scripts/createAndLabelErrorInGXFv3.py -i reference.gtf -g genome.fa -o gene_errors.gtf -e 0.2 --seed 42 --error-types gene_fusion gene_split --max-fusion-gap 10000 --min-orf-aa 30
```

## 错误选择

**本版本改为按基因选择事件**：每个尚未参与事件的基因，以 `-e` 概率尝试引入错误。随机打乱允许的类型，按此顺序尝试，每轮每种允许类型最多构建一次候选，最多重试 `--max-attempts` 轮（默认 20）；采用第一个满足蛋白变化条件的事件。一个基因最多参与一个事件；融合会同时消费相邻基因，因此实际被改变的基因比例不保证等于 `-e`。`--seed` 保证相同输入和参数可复现。

| 类型 | 操作 |
| --- | --- |
| exon_deletion | 删除一个 exon；至少保留一个 exon |
| exon_truncation | 从一侧截短 exon |
| exon_shift | 整体移动 exon，拒绝重叠和越界结果 |
| add_exon | 在内含子中插入 3–90 bp exon，与两侧原 exon 均至少间隔 2 bp |
| merge_exon | 合并两个相邻 exon，保留中间的内含子序列 |
| split_exon | 移除一个碱基形成新内含子，将 exon 拆为两个 |
| gene_fusion | 合并同染色体、同方向的相邻基因；间隔需为 2 至 `--max-fusion-gap` bp，默认上限 10000 |
| gene_split | 在所有异构体共同的内含子间隙中拆开一个基因；间隙至少 2 bp，不允许跨越任何异构体的 exon |

融合不跳过中间基因，也不接受基因跨度重叠的情况。两个来源基因整体替换为一个新基因，转录本按两侧所有异构体的笛卡尔积生成。默认在两侧来源 CDS 内随机选择连接位置，保留转录方向上前一个基因的前段和后一个基因的后段，再拼接为新的 exon 链并预测 ORF。基因间序列作为内含子，不直接拼接两条原蛋白；异构体很多时可能产生大量融合转录本。

拆分会产生两个新基因。所有来源异构体按同一位置切开；只落在一侧的异构体仅在该侧输出，不生成空转录本。split1、split2 按转录方向编号，负链 split1 位于较高基因组坐标。单 exon 基因或没有共同间隙的基因当前不能实施 gene_split。

exon 层面的事件保留 gene_id，但为成功修改的每个转录本生成新的 transcript_id。不同异构体的随机 exon 操作独立；无法操作的异构体原样保留。基因融合、拆分生成新 gene_id 和 transcript_id。生成的 ID 会检查与已有 gene/transcript/exon/protein ID 的冲突。

## 默认要求实际蛋白变化

成功引入错误（error_status=applied）的条件是：最终预测蛋白与每个来源蛋白的氨基酸序列都不同。至少一个氨基酸的替换、插入、缺失，或整个 ORF 丢失，均算作变化；不要求达到人为设定的 identity 损失阈值。不能只凭修改了 CDS 坐标或 CDS 长度来判定。

- 仅改变 UTR、同义变化、或重新预测后仍得到原蛋白的候选不会标为 applied，会继续尝试。
- 每个来源必须原本有 CDS，翻译得到非空且不含 X 的蛋白；新蛋白也不能含 X。无法评价的来源保留原注释，通常标 no_evaluable_source_protein。
- 融合产物必须同时不同于两个来源蛋白。拆分的所有产物都必须不同于各自来源。任一个产物不满足时整个结构事件回退，不丢弃异构体来凑成功结果。
- exon 错误按转录本检查；同一基因中通过检查的转录本被修改，其余保留。
- 无可实施的结构操作时标 no_applicable_error；构建过候选但尝试上限内未找到满足蛋白变化条件的结果时标 no_protein_changing_error。均保留原始注释，error_type=none。
- 因此 `-e 1.0` 表示每个基因尝试修改，不保证每个基因一定能成功。有限次随机尝试失败也不证明数学上不存在可成功的候选。严格模式会产生蛋白改变样本的筛选偏向，不再包含不影响蛋白的负例。
- 如需旧的、不强制蛋白变化的行为，可使用 `--allow-protein-unchanged`；此时只尝试一轮，融合恢复为连接完整来源 exon 链。

重试阶段只比较翻译序列，不计算全局比对；仅最终输出才计算 identity，以控制运行时间。但增加尝试次数仍会增加结构构建与 ORF 预测耗时。

## CDS 重新预测与 phase

**所有成功修改的转录本，包括 exon 错误，均重新预测 ORF**；未修改转录本保留原注释。

1. 按转录方向拼接修改后的 exon 序列，负链取反向互补。
2. 在三个阅读框中扫描 ATG 起始、至第一个同框终止密码子的 ORF；没有 stop 的候选延伸至最后一个完整密码子。
3. 选择氨基酸长度最长的候选，同长度优先转录本中更靠前的起点。完整与不完整候选一起按长度排序；因此最长候选不一定完整。
4. 最短长度由 `--min-orf-aa` 控制，默认 1 aa；可提高至 30、50 等以减少短随机 ORF。
5. 将转录本坐标映射回基因组，重建 CDS、start_codon、stop_codon、five_prime_utr、three_prime_utr。无 stop 时不标注 three_prime_utr；没有合格 ORF 时保留 transcript/exon，不输出 CDS 或密码子特征。
6. CDS 不包含终止密码子，stop_codon 单独输出。密码子可以跨 exon。首段 CDS phase=0；后续按转录方向累计 CDS 长度 L，phase = (3 − L % 3) % 3。phase 表示距离下一个完整密码子的偏移，翻译拼接时不能逐段裁剪。

这是基础最长 ORF 规则，**不是经转录组/同源蛋白证据支持的基因预测器**。仅使用标准核遗传密码，不模拟 NMD、RNA 编辑、非 ATG 起始或翻译例外。重新选择 ORF 本身可能与原始注释所选 ORF 不同，因此 protein_changed 表示最终产物与原注释蛋白的差异，不证明差异完全由结构错误造成。

`cds_prediction=longest_ATG_ORF` 记录预测方法；`orf_status` 为 complete、partial_3prime 或 no_orf。原始未修改注释标为 original_annotation。

## ID 和层级

- 每个输出基因有一条 gene 行；其范围覆盖输出子特征。
- 修改后的 transcript/exon/CDS/密码子/UTR 使用统一的新 gene_id 和 transcript_id。
- exon 按转录方向重新编号，exon_id 和 protein_id 重新生成。旧 Name、Parent、旧 protein_id 等不复制到新产物。
- source_gene_ids、source_transcript_ids 和 event_id 保存来源关系；源基因的所有异构体参与结构事件，原基因行不会残留。

## 输出与比较

输出文件：

- `<output>`：GTF。
- `<output>.error_report.tsv`：每个输出转录本与每个来源转录本一行。融合通常两行，拆分产物分别与自己的来源比较。
- `<output>.original_proteins.fa`：原蛋白，以原 transcript_id 为标题。
- `<output>.modified_proteins.fa`：输出蛋白，以新 transcript_id 为标题；通过 TSV 配对，不再假定两份 FASTA 的 ID 相同。编码区丢失时可能存在空序列条目。

| 标签 | 含义 |
| --- | --- |
| error_type / error_status | 操作类型；applied、not_selected、no_applicable_error、no_protein_changing_error 或 no_evaluable_source_protein |
| gene_id / transcript_id | 输出 ID |
| source_gene_id / source_transcript_id / event_id | 来源和事件编号，支持多对多追踪 |
| coding_impact | protein_changed、protein_unchanged、protein_gained、noncoding 或 unknown |
| protein_identity_pct | 全局比对一致氨基酸数 / 总比对列数 × 100，gap 列计入分母 |
| protein_identity_loss_pct | 100 − identity |
| original_aligned_pct / modified_aligned_pct | 全局比对中两边均为氨基酸的列数 / 对应蛋白长度 × 100，不是局部同源覆盖率 |
| original_protein_length / modified_protein_length | 翻译至第一个 stop 前的氨基酸长度 |
| cds_length_delta | 新 CDS 总长 − 来源 CDS 总长，bp；新 CDS 排除 stop，原 CDS 则按输入定义统计 |
| frameshift | ORF 重新预测后为 NA：新 ORF 不具有可直接继承的旧阅读框锚点，不能简单用长度差判断移码 |
| comparison_status | ok、no_original_CDS、序列不确定等原因 |

融合产物的 GTF 中数值型比较标签为 NA，comparison_status=multiple_sources_see_report；每个来源的 identity、长度及比对比例以 TSV 为准。融合转录本可能只预测出其中一侧的 ORF，并不保证形成包含两个来源蛋白的嵌合蛋白。gene_fusion 是结构错误标签，coding_impact 是序列比较结果。

默认严格模式不修改原注释无 CDS 的来源。关闭严格模式后，原注释没有 CDS 时，新产物有 CDS 则标 protein_gained，否则 noncoding，不计算 identity。含不确定氨基酸 X 或原蛋白为空时标 unknown。无合格新 ORF、而原蛋白有效时，identity 为 0，蛋白长度为 0，标 protein_changed。

全局比对评分为匹配 +1、不匹配 −1、每个 gap 位 −1；并列优先 diagonal、deletion、insertion。与 BLAST 局部 identity 不等价。耗时 O(原蛋白长度 × 新蛋白长度)，内存线性；相同蛋白直接返回 100%。FASTA 整体载入内存，大基因组和大量长蛋白需注意资源消耗。

GTF phase 定义参考：[Ensembl GFF/GTF 格式说明](https://static.ensembl.org/info/website/upload/gff.html)。

验证：

```sh
python -m unittest discover -s other_scripts -p test_createAndLabelErrorInGXFv3.py
```

## 加速与进度日志

默认 `--alignment-backend auto`：环境有 NumPy 时，长蛋白使用 NumPy 向量化全局比对；没有则退回纯 Python。也可指定 `--alignment-backend numpy`，缺少 NumPy 时立即报错，避免无意使用慢速路径。安装：`python -m pip install numpy`。

加速保留原始评分、并列选择、identity 和比对比例，不采用近似比对。同序列直接返回；一条蛋白是另一条的子序列时使用可达到全局最优评分上界的线性算法；其他情况使用精确动态规划。重复序列的翻译与蛋白比较使用有容量上限的缓存。时间复杂度最坏仍为 O(n*m)。

本地人工基准：1500 aa 蛋白对、约每 11 位替换一次，纯 Python 约 2.73 s，NumPy 约 0.10 s；该样本约快 27 倍，不代表全流程或所有数据的提速。1000 组随机蛋白对与原算法的三项输出完全一致；自动测试也覆盖并列最优情况与无 NumPy 回退。

默认每 100 个基因/输出转录本记录一次进度；比较阶段也会在完成一条后检查是否已超过 30 秒。日志写到 stderr，并立即刷新。`--progress-every 500` 调整频率，`--progress-every 0` 关闭日志。每个阶段有开始/完成提示；大型输入读入、单个超长蛋白比较或写出期间可能暂时没有新行。

Linux 后台运行示例：

```sh
python -u createAndLabelErrorInGXFv3.py -i BUSCO.gtf -g Arabidopsis_thaliana.genome.fna -o modified.gtf -e 1.0 --seed 42 --alignment-backend numpy > modified.log 2>&1 &
tail -f modified.log
```

`&` 只负责后台执行，不会启用多核加速。`-e 1.0` 会让所有基因尝试引入错误，基因融合还可能产生多种异构体组合；大量组合的 ORF 预测、注释构建和输出也会增加耗时。
