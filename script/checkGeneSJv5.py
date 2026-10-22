import argparse
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import sys
from intervaltree import IntervalTree
from matplotlib_venn import venn2
import numpy as np
import seaborn as sns
import matplotlib

#需要进行按照频次对SJ进行分级

def load_gene_regions(genebed_file):
    """加载基因区域信息 {GeneID: (chr, start, end, strand)}"""
    gene_regions = {}
    with open(genebed_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 6: continue
            try:
                chr, start, end,  geneID, _, strand,gene = fields[:7]
                gene_regions[geneID] = (chr, int(start), int(end), strand,gene)
            except ValueError:
                print(f"\n❌ erron in load gene region: {str(e)}", file=sys.stderr)
                pass
    return gene_regions

def load_gene_sjs(geneSJ_file):
    """加载基因SJ信息 {geneID: set((chr, start, end, strand))}"""
    gene_sjs = defaultdict(set)
    gene_map = {}
    with open(geneSJ_file) as f:
        for line in f:
            if line.startswith("GeneID") is False:
                fields = line.strip().split('\t')
                #print(fields)
                if len(fields) < 7: continue
                try:
                    geneID, gene, chr, start, end, _, strand = fields[:7]
                    gene_sjs[geneID].add((gene, chr, int(start), int(end), strand))
                except ValueError:
                    print(f"\n❌ error in load gene sjs: {str(e)}", file=sys.stderr)
                    pass
    return gene_sjs

def load_rnaseq_sjs(rnaseqSJ_file):
    """加载RNA-seq SJ信息 set((chr, start, end, strand,cov))"""
    rnaseq_sjs = set()
    with open(rnaseqSJ_file) as f:
        for line in f:
            fields =line.strip().split('\t')
            if len(fields) < 4: continue
            try:
                chr, start, end, strand = fields[:4]
                cov= int(fields[4])
                rnaseq_sjs.add((chr, int(start), int(end), strand,int(cov)))
            except ValueError:
                print(f"\n❌ erron in load rna_seq sjs: {str(e)}", file=sys.stderr)
                pass
    return rnaseq_sjs

def convert_sj_set_to_df(rnaseq_sjs):
    """将 rnaseq_sjs（集合）转为 DataFrame"""

    return pd.DataFrame(
        list(rnaseq_sjs),  # 集合→列表→DataFrame
        columns=["chr", "start", "end", "strand", "cov"]
    )

def level_class(df_A):
    df_A=df_A.sort_values("cov", ascending=False)

    high_threshold = df_A["cov"].quantile(0.75)  # 前25%为High
    medium_threshold = max(df_A["cov"].quantile(0.25),20)  # 前75%为Moderate

    df_A["depth_class"] = "Low"  # 默认Low
    df_A.loc[df_A["cov"] > medium_threshold, "depth_class"] = "Moderate"
    df_A.loc[df_A["cov"] > high_threshold, "depth_class"] = "High"
    df_A["depth"] = np.clip(df_A["cov"], a_min=None, a_max=high_threshold)
    
    print(df_A)
    plt.figure(figsize=(10, 6))        
    # 绘制直方图
    sns.histplot(df_A["depth"], bins=100, kde=True, color="#ffae4c")
    print("The histogram done.")
    # 标注分位数
    plt.axvline(high_threshold, color='r', linestyle='--', label=f"High (Top 25%): {high_threshold:.2f}")
    plt.axvline(medium_threshold, color='g', linestyle='--', label=f"Moderate (Top 75%): {medium_threshold:.2f}")
        
    plt.title("Distribution of RNA-seq scores (with High/Moderate/Low thresholds)")
    plt.xlabel("depth")
    plt.ylabel("Count")
    plt.legend()
    plt.ylim(0,10000)
    plt.savefig(f"SJ_cov_distribution.svg", dpi=300)
    plt.savefig(f"SJ_cov_distribution.png", dpi=300)
    return df_A

##评估每个基因SJ与转录组SJ的情况
def compare_sjs(gene_regions, gene_sjs, rnaseq_sjs):
    """优化后的基因SJ比较函数（改进 total_sj_matched 格式）"""
    results = []

    # 预处理：建立 chr->strand->{sjs} 的索引
    rnaseq_index = defaultdict(lambda: defaultdict(set))
    for sj in rnaseq_sjs:
        chr, start, end, strand, cov = sj
        #if cov >=5 : # 是否不要过滤比较好？
        rnaseq_index[chr][strand].add((start, end))  # 只存 (start, end) 便于查找

    # rnaseq_sj_chr_strand_cov 保留完整信息，用于恢复 total_sj_matched
    rnaseq_sj_chr_strand_cov = {
        (chr, start, end, strand): (cov) for chr, start, end, strand, cov in rnaseq_sjs
    }

    # 收集所有基因SJ用于找未匹配
    all_gene_sjs = set().union(*gene_sjs.values())
    unmatched_sjs = rnaseq_sjs - all_gene_sjs

    # 改用列表存储 matched_sjs（保留 chr, start, end, strand, cov）
    total_sj_matched = []  
    #print(gene_regions)
    for geneID, ( chr, g_start, g_end, strand,gene) in gene_regions.items():
        g_sjs = gene_sjs.get(geneID, set())
        total = len(g_sjs)
    
        # 获取该基因所在区域和链的所有RNA-seq SJ（仅 (start, end) ）
        chr_strand_sjs = rnaseq_index.get(chr, {}).get(strand, set())

        matched = 0
        novel = 0

        # 遍历基因 SJ，检查是否匹配 RNA-seq
        #print(g_sjs)
        for sj in g_sjs:
            _, sj_chr, sj_start, sj_end, sj_strand = sj  # 基因SJ格式（无cov）
            if (sj_start, sj_end) in chr_strand_sjs:
                matched += 1
                # 从 rnaseq_sjs 提取完整信息（chr, start, end, strand, cov）
                sj_key = (chr, sj_start, sj_end, strand)
                if sj_key in rnaseq_sj_chr_strand_cov:
                    cov = rnaseq_sj_chr_strand_cov[sj_key]
                    sj_matched = (chr, sj_start, sj_end, strand, cov)
                    total_sj_matched.append(sj_matched)

        # 检查基因区域内的 RNA-seq 特有 SJ（novel_sjs）
        for (sj_start, sj_end) in chr_strand_sjs:
            if (g_start <= sj_start and sj_end <= g_end) and \
               (sj_start, sj_end) not in {(s[1], s[2]) for s in g_sjs}:
                novel += 1

        results.append({
            'GeneID':geneID,
            'genename': gene,
            'Total_gene_sjs': total,
            'Matched_sjs': matched,
            'Novel_sjs': novel,
            'Status': "LEVEL A" if (matched == total and total > 0 and novel == 0) else
                      "LEVEL A" if (total == 0 and novel == 0) else  # 单外显子基因
                      "NA" if (total > 0 and matched+novel == 0) else      # 全部未匹配
                      "LEVEL A" if (total == matched and novel > 0) else  # 匹配但有新 SJ
                      "LEVEL B" if (total != matched) else
                      "Other situations"
        })

    # 转换 total_sj_matched 为 DataFrame
    total_sj_matched_df = pd.DataFrame(
        total_sj_matched,
        columns=["chr", "start", "end", "strand", "coverage"]
    ).drop_duplicates()  # 避免重复记录

    return pd.DataFrame(results), unmatched_sjs, total_sj_matched_df


def plot_sj_overlap(
    df_A_level: pd.DataFrame,
    total_sj_matched_df: pd.DataFrame,
    save_path: str,
    coverage_classes: list = ["High", "Moderate"],
    color_dict: dict = {"Mismatched": "#b8b8b8","Matched": "#ffaa4f", },
    figsize = (7, 8),
    ) -> plt.Axes:

    # Step 1. 数据预处理：标记匹配情况
    df_A_level = df_A_level.copy()
    total_sj_matched_df = total_sj_matched_df.copy()
    
    # 生成唯一键（chr, start, end, strand）
    df_A_level["sj_key"] = df_A_level.apply(
        lambda x: (x["chr"], x["start"], x["end"], x["strand"]), axis=1
    )
    total_sj_matched_df["sj_key"] = total_sj_matched_df.apply(
        lambda x: (x["chr"], x["start"], x["end"], x["strand"]), axis=1
    )
    
    # 标记是否匹配
    df_A_level["Matched"] = df_A_level["sj_key"].isin(total_sj_matched_df["sj_key"])
    df_filtered = df_A_level[df_A_level["depth_class"].isin(coverage_classes)]
    
    # Step 2. 统计分组计数
    count_data = (
        df_filtered.groupby(["depth_class", "Matched"])
        .size()
        .unstack(fill_value=0)
        .rename(columns={True: "Matched", False: "Mismatched"})
    )
    
    count_data_total = (
        df_A_level.groupby(["depth_class", "Matched"])
        .size()
        .unstack(fill_value=0)
        .rename(columns={True: "Matched", False: "Mismatched"})
    )

    print(count_data_total)
    count_data = count_data[["Mismatched","Matched"]]  
    count_data_total = count_data_total[["Mismatched","Matched"]]  

    unmatch_ratio=count_data["Mismatched"].sum()/(count_data["Mismatched"].sum()+count_data["Matched"].sum())*100
    match_ratio=count_data["Matched"].sum()/(count_data["Mismatched"].sum()+count_data["Matched"].sum())*100
    count_data_total.loc["Total ratio (%)"]=[unmatch_ratio,match_ratio]

    # Step 3. 绘制堆叠条形图
    fig, ax = plt.subplots(figsize=figsize)

    count_data.plot(
        kind="bar",
        stacked=True,
        color=[color_dict["Mismatched"],color_dict["Matched"] ],
        ax=ax,
        width=0.6
    )

    # 添加数值标签
    for i, (matched, unmatched) in enumerate(zip(count_data["Matched"], count_data["Mismatched"])):
        total = matched + unmatched
        ax.text(
            i, matched / 2 + unmatched,  
            f"{matched}\n({matched/total:.1%})", 
            ha="center", va="center", 
            color="white", fontweight="bold"
        )
        ax.text(
            i, unmatched / 2 , 
            f"{unmatched}\n({unmatched/total:.1%})", 
            ha="center", va="center", 
            color="black"
        )
    
    # 图表美化
    ax.set_title(
        "Overlap of Specified SJs with RNA-seq Data", 
        fontsize=15, pad=20
    )
    ax.set_xlabel("Coverage Class", fontsize=15)
    ax.set_ylabel("Count of SJs", fontsize=15)
    ax.legend(
        title="Overlap Status",
        bbox_to_anchor=(1.05, 1), 
        loc="upper left"
    )
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.xticks(fontsize=12)  # 调整 x 轴刻度字体大小
    plt.yticks(fontsize=12)  # 调整 y 轴刻度字体大小

    # 保存图表
    plt.savefig("SJ_"+save_path+"_stacked_bar.svg", dpi=300, bbox_inches="tight")
    plt.savefig("SJ_"+save_path+"_stacked_bar.png", dpi=300, bbox_inches="tight")
    
    percent_data = count_data.div(count_data.sum(axis=1), axis=0) * 100

    fig, ax = plt.subplots(figsize=figsize)
    percent_data.plot(
        kind="bar",
        stacked=True,
        color=[ color_dict["Mismatched"],color_dict["Matched"]],
        ax=ax,
        width=0.6
    )

     # 添加数值标签
    for i, (matched, unmatched) in enumerate(zip(percent_data["Matched"], percent_data["Mismatched"])):
        total = matched + unmatched
        ax.text(
            i, matched / 2 + unmatched, 
            f"{matched/total:.1%}", 
            ha="center", va="center", 
            color="white", fontweight="bold"
        )
        ax.text(
            i, unmatched / 2, 
            f"{unmatched/total:.1%}", 
            ha="center", va="center", 
            color="black"
        )
        # 图表美化
    ax.set_title(
        "Overlap of Specified SJs with RNA-seq Data", 
        fontsize=15, pad=20
    )
    ax.set_xlabel("Coverage Class", fontsize=15)
    ax.set_ylabel("precent(%)", fontsize=15)
    ax.legend(
        title="Overlap Status",
        bbox_to_anchor=(1.05, 1), 
        loc="upper left"
    )
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.xticks(fontsize=12)  # 调整 x 轴刻度字体大小
    plt.yticks(fontsize=12)  # 调整 y 轴刻度字体大小    
    # 保存图表
    plt.savefig("SJ_"+save_path+"_percentage_stacked_bar.svg", dpi=300, bbox_inches="tight")
    plt.savefig("SJ_"+save_path+"_percentage_stacked_bar.png", dpi=300, bbox_inches="tight")
    return count_data, count_data_total

def plot_venn(rnaseq_support,gff_support,overlap,output):
    """绘制RNA-seq支持和基因注释支持的韦恩图"""
    
    # 初始化图形
    plt.figure(figsize=(8, 8))
    
    try:
        
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
        _style_venn_diagram(v, rnaseq_support, gff_support)
        
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
        plt.title("SJs supported by RNA-seq", fontsize=15, pad=15)
        
        # 保存图形
        plt.tight_layout()
        plt.savefig(f"SJ_{output}_venn.svg", dpi=300, bbox_inches='tight')
        plt.savefig(f"SJ_{output}_venn.png", dpi=300, bbox_inches='tight')
        
    except Exception as e:
        print(f"Error plotting Venn diagram: {str(e)}")
    finally:
        plt.close()

def _style_venn_diagram(venn_obj, rnaseq_support, gff_support):
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


def drawpieplot(df, output):
    import matplotlib.pyplot as plt
    import numpy as np
    df = df.loc[df.groupby('genename')['Matched_sjs'].idxmax()]
    print(df.head())
    # 统计并排序确保顺序一致
    status_order = ['LEVEL A', 'LEVEL B', 'NA']
    status_counts = df['Status'].value_counts().reindex(status_order, fill_value=0)
    status_counts.to_csv(f"SJ_{output}_junction_match_gene_stat.tsv",sep="\t")
    # 定制颜色方案（红、橙、深灰、浅灰）
    custom_colors = [
        '#FF6B6B',  # 鲜艳红色
        #'#FFA500',  # 橙色
        '#7FB3D5',  # 深灰
        '#D3D3D3'   # 浅灰
    ]
    explode = {
            'LEVEL A': 0.05,
            'LEVEL B': 0.03,
        #    'LEVEL C': 0,
            'NA': 0
        }

    explode_values = [explode[cat] for cat in status_order]

    plt.figure(figsize=(8, 8))
    label_texts = [f'{label} ({count})' for label, count in zip(status_counts.index, status_counts)]
    wedges, texts, autotexts = plt.pie(
        status_counts,
        labels=label_texts,
        colors=custom_colors,
        autopct='%1.1f%%',
        explode=explode_values,
        startangle=0,
        textprops={'fontsize': 12, 'fontweight': 'bold'},
        wedgeprops={'linewidth': 1, 'edgecolor': 'white'}
    )
    
    # 设置百分比标签样式
    plt.setp(autotexts, size=12, color='black')
    
    
    plt.title('Gene SJ Status Distribution', fontsize=15, pad=20)
    plt.tight_layout()
    plt.savefig(f"SJ_{output}_pieplot.svg", dpi=300) 
    plt.savefig(f"SJ_{output}_pieplot.png", dpi=300) 
    plt.close()

def filter_intergenic_sjs(unmatched_sjs, gene_regions):
    """
    过滤掉位于已知基因区间内的剪接位点，只保留真正的基因间区剪接位点
    
    参数:
        unmatched_sjs: 剪接位点列表，格式为[(chr, start, end, strand), ...]
        gene_regions: 基因区间字典，格式为{geneID: (chr, start, end, strand), ...}
    
    返回:
        真正位于基因间区的剪接位点列表
    """
    # 1. 构建IntervalTree结构（按染色体和链方向索引）
    trees = {}
    #print(gene_regions)
    # 首先收集所有染色体的链方向信息
    for gene_info in gene_regions.values():
        chr, start, end, strand,gene = gene_info
        if chr not in trees:
            trees[chr] = {'+': IntervalTree(), '-': IntervalTree()}
    
    # 填充IntervalTree
    for geneID, (chr, start, end, strand,gene) in gene_regions.items():
        trees[chr][strand][start:end] = geneID
    
    # 2. 过滤剪接位点
    truly_intergenic = []
    
    for sj in unmatched_sjs:
        sj_chr, sj_start, sj_end, sj_strand,_ = sj
        
        # 检查该染色体是否存在
        if sj_chr not in trees:
            truly_intergenic.append(sj)
            continue
        
        # 检查该链方向是否存在
        if sj_strand not in trees[sj_chr]:
            truly_intergenic.append(sj)
            continue
        
        # 检查是否与任何基因重叠
        overlapping_genes = trees[sj_chr][sj_strand].overlaps(sj_start, sj_end)
        
        # 如果没有重叠，则保留
        if not overlapping_genes:
            truly_intergenic.append(sj)
    
    return truly_intergenic

def merge_sj_intervals(unmatched_sjs, window_size=100):
    # 按染色体和链分组
    chr_strand_intervals = defaultdict(list)
    for sj in unmatched_sjs:
        chr, start, end, strand ,_ = sj
        chr_strand_intervals[(chr, strand)].append((start, end))

    merged_intervals = []

    for (chr, strand), intervals in chr_strand_intervals.items():
        if not intervals:
            continue

        # 按起始位置排序
        intervals.sort()

        # 初始化合并
        merged = []
        current_start, current_end = intervals[0]

        for start, end in intervals[1:]:
            # 检查是否可以合并(考虑滑窗)
            if start <= current_end + window_size:
                # 扩展当前区间
                current_end = max(current_end, end)
            else:
                # 保存当前合并区间
                merged.append((current_start, current_end))
                current_start, current_end = start, end

        # 添加最后一个区间
        merged.append((current_start, current_end))

        # 转换为BED6格式
        for i, (start, end) in enumerate(merged):
            merged_intervals.append(
                {"Region":f"Merged_{i}",
                 "Chrom":chr,
                 "Start":start,
                 "End":end,
                 "Strand":strand})

    merged_intervals=pd.DataFrame(merged_intervals)
    return merged_intervals

def main(args_list=None):
    parser = argparse.ArgumentParser(
        description='''
        Assess gene annotation quality using splice junctions
        --------------------------
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
    Usage examples:
    # Basic usage (with default arguments):
    %(prog)s --genebed gene.bed --geneSJ geneSJ.bed --RNAseqSJ RNAseqSJ.bed  -o output

    # Output format:
    Generates two files: {output_prefix}.geneCheck.tsv and {output_prefix}.unmatched_RNAseqSJ.tsv
            '''
    )

    # 必需参数
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--genebed',
                        required=True,
                        metavar='FILE',
                        help='Gene annotation file (BED format)')
    required.add_argument('--geneSJ',
                        required=True,
                        metavar='FILE',
                        help='Gene splice junction information (BED format)')
    required.add_argument("-s",'--RNAseqSJ',
                        required=True,
                        metavar='FILE',
                        help='Transcriptome splice junction information (BED format)')

    optional = parser.add_argument_group('Optional arguments')
    required.add_argument('-o','--output',
                        default='output',
                        metavar='PREFIX',
                        help='Output file prefix (default: output)')
    required.add_argument('-w','--window_size',
                        default=100,
                        help='Sliding window size (default: 100)')
                        
    

    # 解析参数
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    #args = parser.parse_args(args_list) if args_list else parser.parse_args()
    args, _ = parser.parse_known_args(args_list) if args_list else parser.parse_args()
    try:
        print("Loading gene regions...")
        gene_regions = load_gene_regions(args.genebed)
        print("Gene regions loaded, count:", len(gene_regions))

        print("Loading gene SJs...")
        gene_sjs = load_gene_sjs(args.geneSJ)
        
        print("Gene SJs loaded, count:", sum(len(v) for v in gene_sjs.values()))

        print("Loading RNA-seq SJs...")
        rnaseq_sjs = load_rnaseq_sjs(args.RNAseqSJ)
        print("RNA-seq SJs loaded, count:", len(rnaseq_sjs))

        print("SJ comparing...")
        results,unmatched_rnaseq_sjs,totol_match_sjs = compare_sjs(gene_regions, gene_sjs, rnaseq_sjs)
        print(results)
        drawpieplot(results,args.output)
        results.to_csv("SJ_"+args.output+"_detail_comparing.tsv",sep="\t",index=False)

        df_A=convert_sj_set_to_df(rnaseq_sjs)
        df_A_level=level_class(df_A)
        count_data, count_data_total=plot_sj_overlap(df_A_level,totol_match_sjs,args.output)
        print("plot_sj_overlap")
        count_data_total.to_csv("SJ_"+args.output+"_junction_match_stat.tsv",sep="\t")

        overlap = count_data["Matched"].sum()
        rnaseq_support= count_data["Mismatched"].sum()
        gff_support=sum(len(v) for v in gene_sjs.values())-overlap

        plot_venn(rnaseq_support,gff_support,overlap,args.output)
        print("plot_sj_venn")
        

        print("unmatching SJ merging...")
        filterData=filter_intergenic_sjs(unmatched_rnaseq_sjs, gene_regions)
        merged_intervals=merge_sj_intervals(filterData,args.window_size)        
        merged_intervals.to_csv("SJ_"+args.output+".unmatch_RNAseqSJ.tsv",sep="\t",index=False)
    except Exception as e:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        tb_lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
        error_details = "".join(tb_lines)
        
        print(f"\n❌ Error details:\n{error_details}", file=sys.stderr)
        sys.exit(1)    
        #print(f"\n❌ 错误: {str(e)}", file=sys.stderr)
        #sys.exit(1)


if __name__ == "__main__":
    main()


