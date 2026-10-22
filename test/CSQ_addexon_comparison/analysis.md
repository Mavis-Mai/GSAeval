# add_exon 前后四个 CSQ 小指标的配对比较

比较文件：`CSQ2db-origin_addexon.CSQ.tsv`（原始）与 `CSQ2db-witherro_addexon.CSQ.tsv`（错误后）。均按用户确认的 12 列读取，首行数据已保留。

原始表 211 个 BUSCO，错误后 208 个 BUSCO，成功配对 206 个；原始独有 5 个、错误后独有 2 个，没有将未配对对象计为 0 分。以下四个指标均使用这 206 对，且无分数缺失。

下降量定义为：原始分数 − 错误后分数，正值为下降。按平均下降量排名：

| 指标 | 原始均值 | 错误后均值 | 平均下降 | 下降对象比例 | 下降至少 0.1 的比例 |
| --- | ---: | ---: | ---: | ---: | ---: |
| CDS_length_score | 0.964 | 0.669 | 0.295 | 97.09% | 82.52% |
| Intron_site_score | 0.979 | 0.736 | 0.244 | 96.12% | 74.27% |
| CDS_number_score | 0.981 | 0.774 | 0.207 | 77.67% | 65.05% |
| Phase_match_score | 0.979 | 0.783 | 0.196 | 83.01% | 69.42% |

## 结论

这次 add_exon 实验中，CDS_length_score 的平均下降最大，且下降覆盖对象最多；其次为 Intron_site_score。CDS_number_score 的平均下降略大于 Phase_match_score，但 Phase_match_score 下降的对象比例更高。因此“变化幅度”与“受影响范围”应分别解释。

新增 exon 可能改变 CDS 长度、匹配关系和阅读框，但本实验还对成功修改的转录本重新预测了最长 ORF。以上是整个模拟流程的最终响应，不能只归因于新增 exon 这一个局部结构。

206 对中有 58 对 GeneB 参考对象切换。只保留 GeneB 相同的 148 对后，四项平均下降依次为：CDS_length 0.311、Intron_site 0.258、CDS_number 0.218、Phase_match 0.204，排名不变。这支持当前排序，但相同 GeneB 不保证原转录本异构体与动态 CDS 匹配路径完全一致。

## 图怎么读

- `score_mirrored_lines.png/.pdf`：上下镜像的频率折线。蓝色原始分布在上，橙色错误后分布在下。横轴是得分，统一以 0.05 分箱，连接箱中心频率点。下半图表示第二组数据，不是负分。每幅图采用相同坐标尺度；两组均为同一批 206 个配对 BUSCO。
- `paired_decrease_ecdf.png/.pdf`：四个指标真实配对下降量的累积分布。横轴为下降量，越靠右表示下降越大；纵轴为下降不超过该数值的对象比例。CDS_length 曲线在大部分中间百分位更靠右。少数曲线会交叉，不能说它对每个 BUSCO 都最敏感。

按偏好不绘制柱形图。两类折线图已完成视觉检查。

`paired_decrease_summary.tsv` 保存正向定义的下降统计；`paired_scores.tsv` 保存逐对象配对数据；`stable_reference_summary.tsv` 保存相同 GeneB 子集统计。
