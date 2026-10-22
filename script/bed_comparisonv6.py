#!/usr/bin/env python3
import pybedtools
import matplotlib.pyplot as plt
import pandas as pd
import warnings
from time import time
import numpy as np
import sys
import seaborn as sns
from matplotlib_venn import venn2
import matplotlib

##要修改一下，按照表达高低进行---分类统计

warnings.filterwarnings("ignore", category=UserWarning)

def validate_inputs(bed1_path, bed2_path):
    """Validate input files exist and are non-empty"""
    import os
    errors = []
    
    if not bed1_path or not bed2_path:
        errors.append("Error: Both --bed1 and --bed2 must be specified")
    
    for path in [bed1_path, bed2_path]:
        if path and not os.path.exists(path):
            errors.append(f"Error: Input file not found - {path}")
        elif path and os.path.getsize(path) == 0:
            errors.append(f"Error: Empty input file - {path}")
    
    if errors:
        print("\n".join(errors) + "\n")
        sys.exit(1)

def timeit(stage):
    """Timing decorator"""
    def decorator(func):
        def wrapper(*args, **kwargs):
            start = time()
            result = func(*args, **kwargs)
            print(f"{stage} time: {time()-start:.2f}s")
            return result
        return wrapper
    return decorator


class BedAnalyzer:
    def __init__(self, bed1_path, bed2_path, output_prefix, min_overlap):
        self.bed1_path = bed1_path
        self.bed2_path = bed2_path
        self.output_prefix = output_prefix
        self.min_overlap = min_overlap
        self.a = None
        self.b = None
        self.df_A= None
        self.qualitydict={}


    @timeit("Initialize the BED file")
    def load_beds(self):
        #df_key=pybedtools.BedTool.from_dataframe(df_)
        self.df_A=pd.read_csv(self.bed1_path,sep="\t")
        #print(self.df_A)
        self.df_A=self.df_A.sort_values("Coverage", ascending=False)
        #self.df_B=pd.read_csv(self.bed2_path)
        #self.a = pybedtools.BedTool(self.bed1_path)
        df_A_=pd.read_csv(self.bed1_path,sep="\t",skiprows=1,usecols=range(6))
        #bambed
        self.a = pybedtools.BedTool.from_dataframe(df_A_).sort().merge()
        self.a_total= pybedtools.BedTool.from_dataframe(df_A_)
        #df_B_=pd.read_csv(self.bed1_path,sep="\t",usecols=[0,1,2,3,4,5])
        #gene
        self.b_total = pybedtools.BedTool(self.bed2_path)
        #self.b_total= pybedtools.BedTool.from_dataframe(df_B_)
        #存在基因区间的冗余，不利于计算
        self.b = pybedtools.BedTool(self.bed2_path).sort().merge()

    @timeit("Grade the gene regions of RNA-seq")
    def level_class(self):
        high_threshold = self.df_A["Coverage"].quantile(0.75)  # 前30%为High
        medium_threshold = self.df_A["Coverage"].quantile(0.25)  # 前60%为Moderate

        self.df_A["score_class"] = "Low"  # 默认Low
        self.df_A.loc[self.df_A["Coverage"] > medium_threshold, "score_class"] = "Moderate"
        self.df_A.loc[self.df_A["Coverage"] > high_threshold, "score_class"] = "High"
        self.df_A["depth"] = np.clip(self.df_A["Coverage"], a_min=None, a_max=high_threshold)
        print("Draw a histogram") 
        plt.figure(figsize=(10, 6))        
        # 绘制直方图
        sns.histplot(self.df_A["depth"], bins=100, kde=True, color="#ffae4c")
        # 标注分位数
        plt.axvline(high_threshold, color='r', linestyle='--', label=f"High (Top 25%): {high_threshold:.2f}")
        plt.axvline(medium_threshold, color='g', linestyle='--', label=f"Moderate (Top 75%): {medium_threshold:.2f}")
        
        plt.title("Distribution of RNA-seq scores (with High/Moderate/Low thresholds)")
        plt.xlabel("depth")
        plt.ylabel("Count")
        plt.legend()
        plt.ylim(0,5000)

        plt.savefig(f"Gene_{self.output_prefix}_cov_distribution.svg")
        plt.savefig(f"Gene_{self.output_prefix}_cov_distribution.png", dpi=300)


    @timeit("Calculate overlapping statistics")
    def calculate_overlaps(self):
        # 计算gene bed在不同表达水平RNA-seq中的分布
        self.qualitydict = {}
        total_b = len(self.b)
        
        # 计算gene bed与各表达水平RNA-seq的交集情况
        self.a_and_b_inB =0
        for key in ["High", "Moderate","Low"]:
            df_ = self.df_A.loc[self.df_A["score_class"] == key]
            df_key = pybedtools.BedTool.from_dataframe(df_)
            
            # Gene bed中与该表达水平RNA-seq有严格重叠的基因数量
            b_and_a_strict = len(df_key.intersect(self.b, f=self.min_overlap, u=True, r=True))
            # Gene bed中与该表达水平RNA-seq有重叠但不满足严格条件的基因数量
            b_and_a_loose =len(df_key.intersect(self.b, u=True))- b_and_a_strict
            # Gene bed中不与该表达水平RNA-seq重叠的基因数量
            a_only_level = len(df_) - b_and_a_strict - b_and_a_loose
            self.a_and_b_inB=b_and_a_strict+b_and_a_loose+self.a_and_b_inB
            self.qualitydict[key] = [
                b_and_a_strict,
                b_and_a_loose,
                a_only_level,
                len(df_)  # 该表达水平的RNA-seq区域总数(保持不变)
            ]
        print(self.qualitydict)
        # 保留原有的全集统计(供venn图使用)
        self.a_and_b_strict_inB = self.a.intersect(self.b, f=self.min_overlap, u=True, r=True)
        #self.b_only = self.b.intersect(self.a, v=True)
        self.b_only = len(self.b)-self.a_and_b_inB
        print(f"self.b_only:{self.b_only}")
        self.a_only = self.a.intersect(self.b, v=True)

    @timeit("Save the statistical results")
    def save_stats(self):
        rnaseq_support = len(self.a_only)
        overlap = self.a_and_b_inB
        overlap_pct_rnaseq = f'{overlap / (rnaseq_support + overlap) * 100:.1f}%' if (rnaseq_support + overlap) > 0 else "0%"
        total_b = len(self.b)
        with open(f"Gene_{self.output_prefix}_stats.tsv", "w+") as f:
            f.write(f"ab initio transcript regions overlap with gene:\n")
            f.write(f"total transcript region \t{rnaseq_support+overlap}\n")
            for k, v in self.qualitydict.items():
                f.write(f"{k} expression ab initio transcript regions overlap with gene region:\n")
                f.write(f"  strict overlap gene count\t{v[0]} ({(v[0]/total_b)*100:.2f} %)\n")
                f.write(f"  loose overlap gene count\t{v[1]} ({(v[1]/total_b)*100:.2f} %)\n")
                f.write(f"  no overlap gene count\t{v[2]} ({(v[2]/total_b)*100:.2f} %)\n")
                f.write(f"  total {k} ab initio regions\t{v[3]}\n")
                f.write(f"\n")
            
            f.write(f"Total overlap region:\t{overlap}({overlap_pct_rnaseq})\n")
            f.write(f"\n")
            f.write(f"Genes overlap with ab initio transcript regions:\n")
            f.write(f"total gene \t{len(self.b)}\n")
            f.write(f"overlap\t{self.a_and_b_inB}({(self.a_and_b_inB/len(self.b))*100:.2f})\n")
            f.write(f"without RNA-seq support\t{len(self.b)-self.a_and_b_inB}({(len(self.b)-self.a_and_b_inB)/len(self.b)*100:.2f})\n")


    @timeit("Create a quick search dictionary")
    def create_lookup_dicts(self):
        # 读取RNA-seq BED并创建坐标到行的映射
        df_a = pd.read_csv(self.bed1_path, sep='\t', skiprows=1, usecols=[0,1,2,3,4])
        df_a.columns = ["Chrom", "Start", "End", "Strand", "Coverage"]
        self.a_dict = {(r.Chrom, r.Start, r.End, r.Strand): r.Coverage
                      for r in df_a.itertuples()}

        # 读取Gene BED并创建坐标到GeneID的映射
        df_b = pd.read_csv(self.bed2_path, sep='\t', header=None, usecols=[0,1,2,3,4,5])
        print(df_b.head())
        df_b.columns = ["Chrom", "Start", "End", "GeneID","Coverage","Strand",]
        self.b_dict = {(r.Chrom, r.Start, r.End,r.Strand): r.GeneID
                      for r in df_b.itertuples()}

    @timeit("Handle all regions")
    def process_all_regions(self):
        records = []
        processed_a = set()
        processed_b = set()

        # 1. 先找出所有交集并计算重叠比例
        raw = self.a_total.intersect(self.b_total, wao=True).to_dataframe(
            names=['chrom_RNAseq', 'start_RNAseq', 'end_RNAseq','strand_RNAseq',"cov_RNAseq","lena",
               "chrom_gff","start_gff","end_gff","geneid",".","strand_gff","genename","overlap_lenb"])
        test = self.a_total.intersect(self.b_total, wao=True).to_dataframe()
        #print(raw.iloc[0])
        # 创建一个字典来保存每个A区域的最佳重叠记录
        best_overlap_for_gene = {}  
    
        for i in range(len(raw)):
            row = raw.iloc[i]
            gene_id = row["geneid"]  # 新增：提前获取gene_id
            key_a = (row["chrom_RNAseq"], row["start_RNAseq"], row["end_RNAseq"], row["strand_RNAseq"])
            
            # 计算重叠比例
            overlap_len = row["overlap_lenb"]
            genelen = row["end_gff"] - row["start_gff"]
            overlap_pct = 0 if genelen <= 0 else overlap_len / genelen
            
            if gene_id not in best_overlap_for_gene or overlap_pct > best_overlap_for_gene[gene_id]['overlap_pct']:
                best_overlap_for_gene[gene_id] = {
                    'row': row,
                    'overlap_pct': overlap_pct,
                    'key_a': key_a,
                    'key_b': (row["chrom_gff"], row["start_gff"], row["end_gff"], row["strand_gff"])
                }

        # 2. 处理最佳重叠记录
        for gene_id, data in best_overlap_for_gene.items():
            row = data['row']
            key_a = data['key_a']
            key_b = data['key_b']
            overlap_pct = data['overlap_pct']
            
            is_significant = "Yes" if overlap_pct >= self.min_overlap else "No"
            coverage = self.a_dict.get(key_a, "NA")
            strand = key_a[-1]
            
            records.append([
                gene_id,  # 直接使用外层循环的gene_id
                "Yes",    # In_RNAseq
                "Yes",    # In_Gene
                *key_a,
                coverage,
                f"{overlap_pct*100:.1f}%" if overlap_pct < 1 else "100%",
                is_significant
            ])
            
            processed_a.add(key_a)
            processed_b.add(key_b)
        # 3. 处理仅在A中的区域 (保持不变)
        for r in self.a_only:
            key_a = (r.chrom, r.start, r.end, r.strand)
            if key_a not in processed_a:
                coverage = self.a_dict.get(key_a, "NA")
                strand = key_a[-1]
                records.append([
                    "NA",   # GeneID
                    "Yes",  # In_RNAseq
                    "No",   # In_Gene
                    *key_a,
                    coverage,
                    "0%",   # Overlap_Percentage
                    "No"    # Significant_Overlap
                ])

        # 4. 处理仅在B中的区域 (保持不变)
        b_only = self.b.intersect(self.a, v=True)
        for r in b_only:
            key_b = (r.chrom, r.start, r.end, r.strand)
            if key_b not in processed_b:
                gene_id = self.b_dict.get(key_b, "NA")
                strand = key_b[-1]
                records.append([
                    gene_id,
                    "No",   # In_RNAseq
                    "Yes",  # In_Gene
                    *key_b,
                    "NA",
                    "0%",   # Overlap_Percentage
                    "No"    # Significant_Overlap
                ])

        return records


    @timeit("Save the detailed results")
    def save_results(self, records):
        df = pd.DataFrame(
            records,
            columns=[
                "GeneID","In_RNAseq", "In_Gene",
                "Chrom", "Start", "End","Strand",
                "Coverage_in_RNAseq",             
                "Gene_Overlap_Percentage",
                "Significant_Overlap"
            ]
        )

        df.to_csv(f"Gene_{self.output_prefix}_detailed_overlap.tsv", sep="\t", index=False)

    @timeit("Draw a stacked bar chart")
    def plot_stacked_bar(self):
        plt.figure(figsize=(7, 8))
        
        # Total stats
        total_a = len(self.a)
        total_b = len(self.b)
        
        # Gene bed overlaps (严格+宽松)
        b_overlap_high = self.qualitydict["High"][0]+ self.qualitydict["High"][1] # Genebed strictly overlaps with High RNA-seq
        b_overlap_moderate = self.qualitydict["Moderate"][0]+self.qualitydict["Moderate"][1]  # Genebed strictly overlaps with Moderate RNA-seq
        b_overlap_low = self.qualitydict["Low"][0]+self.qualitydict["Low"][1]
        #b_no_overlap = len(self.b_only) # Genebed with no RNA-seq overlap
        b_no_overlap = total_b-b_overlap_high-b_overlap_moderate-b_overlap_low# Genebed with no RNA-seq overlap
        
        # RNA-seq only (not overlapping any gene)
        a_only = len(self.a_only)
        
        # Categories & Heights
        categories = [
            "Genes",
            "RNA-seq\nevidence"
        ]
        
        heights = [
            b_overlap_high + b_overlap_moderate+b_overlap_low+b_no_overlap,  # Genebed overlaps RNA-seq (any level)
            len(self.a)
        ]
        
        # Bar colors (split Genebed-overlapping into High/Moderate)
        overlap_colors = ["#da5b10", "#f7a649","#f7d2a8"]  # Dark blue (High), Light blue (Moderate)
        no_overlap_color = "#67add6"  # Gray
        rnaseq_only_color = "#67add6"  # Orange
        
        # Plot stacked bars
        plt.bar(categories[0], b_no_overlap, color=no_overlap_color,label="no RNA-seq evi.", edgecolor="black")
        plt.bar(categories[0], b_overlap_low, color=overlap_colors[2],bottom=b_no_overlap,label="Low", edgecolor="black")
        plt.bar(categories[0], b_overlap_moderate, color=overlap_colors[1]
                ,bottom=b_no_overlap+b_overlap_low,label="Moderate", edgecolor="black")
        plt.bar(categories[0], b_overlap_high, color=overlap_colors[0],
                bottom=b_no_overlap+b_overlap_low+b_overlap_moderate,label="High", edgecolor="black")

        plt.bar(categories[1], a_only,color=rnaseq_only_color, 
                label="Noval genes", edgecolor="black", hatch=".")
        
        plt.bar(categories[1], len(self.a)-a_only,bottom=a_only, color=overlap_colors[2], 
                label="Annotated genes", edgecolor="black",hatch=".")      


        # Add total counts on top of each bar
        for i, h in enumerate(heights):
            plt.text(i, h + 0.01 * max(heights), str(h), ha='center', va='bottom', fontsize=15)

        # Legend & Labels
        plt.legend(loc="upper right")
        plt.ylabel("Number of Regions")
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        
        plt.tight_layout()
        plt.savefig(f"Gene_{self.output_prefix}_stacked_bar.svg")
        plt.savefig(f"Gene_{self.output_prefix}_stacked_bar.png", dpi=300)
        plt.close()

    #要不考虑堆一起？或者两个以gff出发，和以RNA-seq出发
    def plot_stacked_bar_ratio(self):
        plt.figure(figsize=(7, 8))
        
        # Total stats
        total_a = len(self.a)
        total_b = len(self.b)

        # Gene bed overlaps (严格+宽松)
        b_overlap_high = self.qualitydict["High"][0]+ self.qualitydict["High"][1] # Genebed strictly overlaps with High RNA-seq
        b_overlap_moderate = self.qualitydict["Moderate"][0]+self.qualitydict["Moderate"][1]  # Genebed strictly overlaps with Moderate RNA-seq
        b_overlap_low = self.qualitydict["Low"][0]+self.qualitydict["Low"][1]
        b_no_overlap = total_b-b_overlap_high-b_overlap_moderate-b_overlap_low# Genebed with no RNA-seq overlap
       
        # RNA-seq only (not overlapping any gene)
        a_only = len(self.a_only)
        
        # Categories & Heights
        categories = [
            "Genes",
            "RNA-seq\nevidence"
        ]
        # Calculate percentages
        #total = b_overlap_high + b_overlap_moderate + b_no_overlap + a_only
        
        # Split the overlapping portion into High/Moderate percentages
        overlap_high_pct = 100 * b_overlap_high / total_b
        overlap_moderate_pct = 100 * b_overlap_moderate / total_b
        overlap_low_pct = 100 * b_overlap_low / total_b
        b_no_overlap_pct =  100 * b_no_overlap/ total_b
        # Bar colors (split Genebed-overlapping into High/Moderate)
        overlap_colors = ["#da5b10", "#f7a649","#f7d2a8"]  # Dark blue (High), Light blue (Moderate)
        no_overlap_color = "#67add6"  # Gray
        rnaseq_only_color = "#67add6"  # Orange
        
        # Plot stacked bars as percentages
        plt.bar(categories[0], b_no_overlap_pct,color=no_overlap_color,label=f"No expression ({b_no_overlap})", edgecolor="black")
        plt.bar(categories[0],overlap_low_pct, bottom=b_no_overlap_pct,label=f"Low ({b_overlap_low})",edgecolor="black",color=overlap_colors[2])
        plt.bar(categories[0], overlap_moderate_pct, bottom=b_no_overlap_pct+overlap_low_pct, label=f"Moderate ({b_overlap_moderate})",
                    edgecolor="black",color=overlap_colors[1])
        plt.bar(categories[0], overlap_high_pct, color=overlap_colors[0], bottom=b_no_overlap_pct+overlap_low_pct+overlap_moderate_pct,
                label=f"High ({b_overlap_high})", edgecolor="black")

        a_only_pct=100*a_only/len(self.a)
        plt.bar(categories[1], a_only_pct, color=rnaseq_only_color, hatch=".",
                label=f"Noval genes ({a_only})", edgecolor="black")
        plt.bar(categories[1], 100-a_only_pct, bottom= a_only_pct, color=overlap_colors[2], hatch=".",
                label=f"Annotated genes ({len(self.a)-a_only})", edgecolor="black")

        text_bg_props = dict(
            boxstyle="round,pad=0.3",      # 圆角矩形，内边距 0.3
            facecolor="white",            # 背景色
            edgecolor="black",            # 边框色
            alpha=0.8,                    # 透明度
            linewidth=1                   # 边框粗细
        )

        plt.text(categories[0], b_no_overlap_pct/2,f"{b_no_overlap_pct:.1f}%", ha='center', va='center', fontsize=12)
        plt.text(categories[0], b_no_overlap_pct+overlap_low_pct/2, f"{overlap_low_pct:.1f}%", ha='center', va='center', fontsize=12)
        plt.text(categories[0], b_no_overlap_pct+overlap_low_pct+overlap_moderate_pct/2,
                f"{overlap_moderate_pct:.1f}%", ha='center', va='center', fontsize=12)
        plt.text(categories[0], b_no_overlap_pct+overlap_low_pct+overlap_moderate_pct+overlap_high_pct/2,
                f"{overlap_high_pct:.1f}%", ha='center', va='center', fontsize=12)
        plt.text(categories[1], a_only_pct/2,f"{a_only_pct:.1f}%", ha='center', va='center', fontsize=12, bbox=text_bg_props)
        plt.text(categories[1], a_only_pct+(100-a_only_pct)/2,f"{100-a_only_pct:.1f}%", ha='center', va='center', fontsize=12,bbox=text_bg_props)


        # Legend & Labels
        plt.legend(loc="upper right")
        plt.ylabel("Percentage of Regions (%)")
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        
        # Set y-axis limit to 100%
        plt.ylim(0, 100)
        
        plt.tight_layout()
        plt.savefig(f"Gene_{self.output_prefix}_stacked_bar_pct.svg")
        plt.savefig(f"Gene_{self.output_prefix}_stacked_bar_pct.png", dpi=300)
        plt.close()

    

    def plot_venn(self):
        """绘制RNA-seq支持和基因注释支持的韦恩图"""
        
        # 初始化图形
        plt.figure(figsize=(8, 8))
        
        try:
            # 计算各区域值
            rnaseq_support = len(self.a_only)
            gff_support = self.b_only
            overlap = self.a_and_b_inB
            
            # 计算重叠比例（考虑分母为零的情况）
            if rnaseq_support == 0 or gff_support == 0:
                overlap_ratio = 0
            else:
                overlap_ratio = overlap / min(rnaseq_support, gff_support)
                
            # 计算重叠百分比
            overlap_pct_rnaseq = f'{overlap / (rnaseq_support + overlap) * 100:.1f}%' if (rnaseq_support + overlap) > 0 else "0%"
            overlap_pct_gff = f'{overlap / (gff_support + overlap) * 100:.1f}%' if (gff_support + overlap) > 0 else "0%"
            
            # 绘制韦恩图
            set_colors = ('#4682B4', '#FF8C00')  # 钢蓝和深橙
            v = venn2(
                subsets=(rnaseq_support, gff_support, overlap), 
                set_labels=('RNA-seq support', 'Gene annotation support'),
                set_colors=set_colors,
                alpha=0.7
            )
            # 设置图形样式
            self._style_venn_diagram(v, rnaseq_support, gff_support)
            
            # 手动添加重叠区域百分比标签
            if overlap > 0 and hasattr(v, 'get_label_by_id'):
                try:
                    label = v.get_label_by_id('11')
                    if label:
                        original_text=label.get_text()
                        label.set_text(f'{original_text}\n{overlap_pct_rnaseq}//{overlap_pct_gff}')
                        label.set_fontsize(15)
                except:
                    pass
            
            # 设置标题
            plt.title("Gene loci supported by RNA-seq", fontsize=14, pad=15)
            
            # 保存图形
            plt.tight_layout()
            plt.savefig(f"Gene_{self.output_prefix}_venn.svg",  bbox_inches='tight')
            plt.savefig(f"Gene_{self.output_prefix}_venn.png", dpi=300, bbox_inches='tight')
            
        except Exception as e:
            print(f"Error plotting Venn diagram: {str(e)}")
        finally:
            plt.close()


    def _style_venn_diagram(self, venn_obj, rnaseq_support, gff_support):
        """设置韦恩图的样式"""
        
        # 设置所有区域的边框
        for patch in venn_obj.patches:
            if patch:  # 某些区域可能为0（比如无重叠）
                patch.set_edgecolor('black')
                patch.set_linewidth(1.5)
        
        # 设置圆环的边框
        if hasattr(venn_obj, 'circles'):
            for circle in venn_obj.circles:
                circle.set_edgecolor('black')
                circle.set_linewidth()

    def run(self):
        self.load_beds()
        self.level_class()
        self.calculate_overlaps()
        self.save_stats()
        self.create_lookup_dicts()
        results = self.process_all_regions()
        #print(results)
        self.plot_stacked_bar()
        self.plot_stacked_bar_ratio()
        self.plot_venn()
        
        self.save_results(results)

        print(f"Analysis completed! The result has been saved to Gene_{self.output_prefix}_*files")

def main(args_list=None):
    import argparse
    parser = argparse.ArgumentParser(description="Compare RNA-seq splicing sites with gene regions", 
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-r","--rnaseqbed", help="RNA-seq BED file (Chr start end strand coverage)")
    parser.add_argument("-g","--genebed", help="Gene BED file (Chr start end geneID)")
    parser.add_argument("-o", "--output", default="bed_comparison",
                       help="output prefiex")
    parser.add_argument("-m", "--min_overlap", type=float, default=0.5,
                       help="Minimum overlap ratio (0-1)")
    #parser.add_argument("-h", "--help", action="store_true",
    #                   help="Display help information")
    args = parser.parse_args(args_list) if args_list else parser.parse_args()

    
    validate_inputs(args.rnaseqbed, args.genebed)
    analyzer = BedAnalyzer(args.rnaseqbed, args.genebed, args.output, args.min_overlap)
    analyzer.run()


if __name__ == "__main__":
    main()