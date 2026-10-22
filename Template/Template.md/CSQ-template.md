# 基因结构准确性评价
这一部分通过同源参考基因集合（**Eudicots**）和转录组数据进行评估。

## 1. CSQ分析
* CSQ分析的简要介绍
<div style="padding: 15px; border: 1px solid #7a7a7a; border-radius: 5px; background-color: #f9f9f9;">
通过比对目标基因结构注释与参考集中同源基因的剪切模式相似度，可有效评估目标注释的准确性。<br>
在此基础上，通过计算同源基因对的CSQ得分，可对基因的剪切保守性进行量化，并将其划分为Good、Moderate和Poor三个等级：<br>
* Good：表示基因与同源基因的剪切模式非常保守；<br>
* Moderate：表示基因与同源基因的剪切模式具有差别；<br>
* Poor：表示基因与同源基因的剪切模式差异较大。
</div>

* 基本信息统计

|  lable   | count | frequency(%) |
| :------: | :---: | :----------: |
|   Good   |  1618 |     96.7     |
| Moderate |  43   |     2.6      |
|   Poor   |  13   |     0.8      |
| Missing  |  24   |     nan      |
*同源基因集合（Eudicots）*

|                    图                     |                   <div style="width:200px">  解析   </div>                    |
| :---------------------------------------: | :---------------------------------------------------------------------------: |
| ![Gene_Gma_venn](_v_images/CSQ.png =500x) | <div style="width:200px">  各类别（Good/Moderate/Poor）基因比例分布饼图  </div> |

* 结论
<div style="background: #f0f8ff; border-left: 4px solid #4682b4; padding: 10px;">
此处填写结论
</div>

---
