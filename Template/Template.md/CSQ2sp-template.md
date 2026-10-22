## 2. CSQ2sp分析
* CSQ分析的简要介绍
<div style="padding: 15px; border: 1px solid #7a7a7a; border-radius: 5px; background-color: #f9f9f9;">
通过与同源物种（拟南芥 TAIR12）基因结构的剪切相似性，进而评估目标物种的剪切准确性。<br>
首先，通过双向最佳比对获取物种间的同源基因对。<br>
随后，计算同源基因对的CSQ得分，该得分越高，表明剪切模式相似性就越高。<br>
最后，统计所有同源基因对的CSQ分布，并提取该分布中的第一四分数（Q1）作为参考基准。<br>
（分布的Q1越接近1，说明两个物种的剪切模式越相似）
</div>

* 同源基因对数量为：10300
* 同源基因的CSQ分布

|                      图                      |               <div style="width:200px">  解析   </div>                |
| :------------------------------------------: | :-------------------------------------------------------------------: |
| ![Gene_Gma_venn](_v_images/CSQ2sp.png =500x) | <div style="width:200px">  同源基因对的分布情况。红色虚线标注Q1。 </div> |

* 列取CSQ得分最高的10个同源基因对情况（即同源基因具有相似的剪切模式）：

| GeneA | GeneB | cdsA_number | cdsB_number | CSQ |
|:---------:|:---------:|:--------:|:--------:|:--------:|
| rna-TAIR12_TAIR12_AT2G25810 | SAM6A429900.1 | 3 | 3 | 1.0 |
| rna-TAIR12_TAIR12_AT1G01230 | SAM1A272700.1 | 3 | 3 | 1.0 |
| rna-TAIR12_TAIR12_AT5G67590 | SAM2A182300.1 | 5 | 5 | 1.0 |
| rna-TAIR12_TAIR12_AT5G67500 | SAM2A180500.1 | 6 | 6 | 1.0 |
| rna-TAIR12_TAIR12_AT5G59890 | SAM1A315800.1 | 3 | 3 | 1.0 |
| rna-TAIR12_TAIR12_AT5G60220 | SAM5A419200.1 | 2 | 2 | 1.0 |
| rna-TAIR12_TAIR12_AT2G23940 | SAM4A129600.1 | 5 | 5 | 1.0 |
| rna-TAIR12_TAIR12_AT5G59840 | SAM3A142800.1 | 8 | 8 | 1.0 |
| rna-TAIR12_TAIR12_AT1G32200 | SAM7A240200.1 | 12 | 12 | 1.0 |
| rna-TAIR12_TAIR12_AT5G48760 | SAM7A375400.1 | 4 | 4 | 1.0 |


* 列取CSQ得分最低的10个同源基因对情况（即同源基因具有较大差异的剪切模式）：

| GeneA | GeneB | cdsA_number | cdsB_number | CSQ |
|:---------:|:---------:|:--------:|:--------:|:--------:|
| rna-TAIR12_TAIR12_AT5G24710 | SAM6A245800.1 | 21 | 2 | 0.193 |
| rna-TAIR12_TAIR12_AT1G24460 | SAM1A144000.1 | 8 | 2 | 0.229 |
| rna-TAIR12_TAIR12_AT3G28875 | SAM2A534800.1 | 2 | 12 | 0.23 |
| rna-TAIR12_TAIR12_AT1G64770 | SAM2A306900.1 | 3 | 21 | 0.249 |
| rna-TAIR12_TAIR12_AT1G13160 | SAM1A133800.1 | 14 | 3 | 0.255 |
| rna-TAIR12_TAIR12_AT1G49475 | SAM2A434600.1 | 3 | 11 | 0.259 |
| rna-TAIR12_TAIR12_AT1G66590-2 | SAM2A366500.1 | 4 | 18 | 0.263 |
| rna-TAIR12_TAIR12_AT3G29010-2 | SAM7A039000.1 | 4 | 17 | 0.264 |
| rna-TAIR12_TAIR12_AT3G30525 | SAM2A484200.1 | 2 | 9 | 0.267 |
| rna-TAIR12_TAIR12_AT2G28230 | SAM2A084000.1 | 3 | 16 | 0.269 |


* 结论
<div style="background: #f0f8ff; border-left: 4px solid #4682b4; padding: 10px;">
这里填写结论
</div>

---
