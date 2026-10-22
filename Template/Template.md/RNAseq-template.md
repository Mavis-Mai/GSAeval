## 3. RNA-seq分析
* RNA-seq分析的简要介绍
<div style="padding: 15px; border: 1px solid #7a7a7a; border-radius: 5px; background-color: #f9f9f9;">
通过高质量的转录组数据，对目标物种的结构注释准确性进行评估。<br>
* 评估基于转录组重构的转录本与已有基因结构的重叠情况<br>
* 评估转录组的剪切位点与已有基因的内含子剪切位点重合性<br>
* 基因结构出错的风险得分（Risk score）。
</div>

* 基本信息统计

|                     Labels                     |   Values (%)   |                                                         explaination                                                         |
| :--------------------------------------------: | :------------: | :--------------------------------------------------------------------------------------------------------------------------: |
|         genes without RNA-seq support          | 44689 (34.86%) |              <div style="width:300px"> 若缺乏转录组数据支持的基因结构比例过高，则提示注释结果中可能存在假基因。 </div>              |
|  RNA-seq support regions without gene models   | 11899 (12.5%)  |    <div style="width:300px"> 若转录组数据重构的基因结构大量缺失于原有注释中（即占比较高），则暗示当前注释结果可能不够完整。 </div>    |
| genes with conflicting splicing junctions (SJ) | 26067 (20.33%) |           <div style="width:300px"> 若与转录组的剪切位点具有冲突的基因结构比例过多，则暗示当前注释剪切准确性较差。 </div>           |
|              CDS distance (Q0.95)              |      0.24      | <div style="width:300px">重构转录本与当前注释基因在CDS的差异度。该数值越大，表明当前注释与转录组重构结果之间的整体差异越显著。 </div> |
|             Overlap score (Q0.05)              |      1.00      |  <div style="width:300px">重构转录本与当前注释基因在基因组上重叠情况。该数值越大，表明当前注释与转录组重构结果之间重叠性好。 </div>   |
|                SJ score (Q0.05)                |      1.00      |     <div style="width:300px">若与转录组的剪切位点与当前注释基因的匹配情况。该数值越大，表明当前注释与转录组的SJ匹配性好。 </div>     |
|               Risk score (Q0.95)               |      0.27      |               <div style="width:300px">风险值是综合以上指标的复合指标。该数值越大，表明当前注释错误风险越高。 </div>               |
*Ps: Q0.05是前0.05的阈值，而Q0.95是前0.95的阈值*

* 图表信息

|                              图                              |                                                                <div style="width:200px">  解析   </div>                                                                 |
| :----------------------------------------------------------: | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------: |
|       ![Gene_Gma_venn](_v_images/RNAseqvenn.png =500x)       |         <div style="width:200px">  基因结构与转录组数据相互重叠情况 <br>(左侧：转录组预测的结构与基因结构注释的重叠比例。右侧：基因结构注释具有转录组支持的比例）  </div>         |
|        ![SJ_Gma_pieplot](_v_images/SJlevel.png =500x)        |              <div style="width:200px">   基因结构的剪切位点是否具有转录组数据的支持 <br>（A: 完全一致的；B 存在与转录组冲突的剪切位点；NA 缺少转录组支持） </div>               |
| ![Gma_riskscore_distribution](_v_images/Riskscore.png =500x) | <div style="width:200px">基因的风险得分 <br>（Risk score）的分布（Risk score 是综合评估基因结构与转录组数据的准确性的指标。Risk score 越大说明基因的结构出错的可能越大）  </div> |

* 结论
<div style="background: #f0f8ff; border-left: 4px solid #4682b4; padding: 10px;">
这里填写结论
</div>

---

## 4. 典型错误案例展示
|                       图                       | <div style="width:200px">  错误解释   </div> |
| :--------------------------------------------: | :-----------------------------------------: |
| ![Gene_Gma_venn](_v_images/example1.png =500x) |     <div style="width:200px">    </div>     |
| ![Gene_Gma_venn](_v_images/example2.png =500x) |     <div style="width:200px">    </div>     |
| ![Gene_Gma_venn](_v_images/example3.png =500x) |     <div style="width:200px">    </div>     |