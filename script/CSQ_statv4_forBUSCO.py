import pandas as pd
import re
from collections import defaultdict
import matplotlib.pyplot as plt
import argparse
import sys

def extract_gene_prefix(gene_name):
    """安全提取基因名称前缀"""
    if pd.isna(gene_name) or not str(gene_name).strip():
        return 'Other'
    cleaned = re.sub(r'[^a-zA-Z0-9]', '', str(gene_name).strip())
    match = re.match(r'^([A-Za-z]+)[0-9]*', cleaned)
    return match.group(1).upper() if match else 'Other'

def calculate_prefix_stats(df,target_col):
    """按基因前缀统计"""
    prefix_stats = defaultdict(list)
    prefix_intron_stats = defaultdict(list)

    found_gene_cols = False
    for col in df.columns:
        if col == "GeneA" or col == "GeneB":
            found_gene_cols = True
            for _, row in df.iterrows():
                prefix = extract_gene_prefix(row[col])
                try:
                    value = float(row[target_col])
                    prefix_stats[prefix].append(value)

                    introna = float(row["cdsA_number"])
                    intronb = float(row["cdsB_number"])
                    prefix_intron_stats[prefix].append(intronb)
                    prefix_intron_stats[prefix].append(introna)

                except (ValueError, TypeError):
                    continue

    if not found_gene_cols:
        print("Warning: No gene columns detected (column names should contain 'gene')")
    elif not prefix_stats:
        print("Warning: No valid gene prefix data extracted")

    stats_data = []
    for prefix, values in prefix_stats.items():
        if values:
            s = pd.Series(values)
            i = pd.Series(prefix_intron_stats[prefix])
            stats_data.append({
                'Prefix': prefix,
                'Count': len(values),
                'CSQ_Mean': round(s.mean(), 4) if len(values) > 0 else 0,
                'intron_Mean': round(i.mean(), 4) if len(values) > 0 else 0,
            })

    return pd.DataFrame(stats_data).sort_values('Count', ascending=False)


def calculate_gene_stats(df, target_col, use_prefix=True, top_n=18):
    """计算基因统计
    
    Args:
        df: 数据框
        target_col: 目标列名
        use_prefix: 是否使用前缀分组统计
        top_n: 提取TOP N个基因
    """
    gene_stats = defaultdict(list)
    gene_prefix = {}
    gene_intron_stats = defaultdict(list)
    genelist = set(list(df["GeneA"])+list(df['GeneB']))
    #print(genelist)
    for key in genelist:
        tmp_df = df.loc[(df['GeneA'] == key) | (df['GeneB'] == key)]
        prefix = extract_gene_prefix(key)
        gene_prefix[key] = prefix
        
        for _, row in tmp_df.iterrows():
            value = float(row[target_col])
            gene_stats[key].append(value)
            introna = float(row["cdsA_number"])
            intronb = float(row["cdsB_number"])
            gene_intron_stats[key].append(intronb)
            gene_intron_stats[key].append(introna)

    stats_data = []
    for gene, values in gene_stats.items():
        if values:
            s = pd.Series(values)
            i = pd.Series(gene_intron_stats[gene])
            stats_data.append({
                'Prefix': gene_prefix[gene],
                'Gene': gene,
                'Count': len(values),
                'CSQ_Mean': round(s.mean(), 4) if len(values) > 0 else 0,
                'Intron_Mean': round(i.mean(), 4) if len(values) > 0 else 0,
            })
    #print(stats_data)
    stats_df = pd.DataFrame(stats_data).sort_values('CSQ_Mean', ascending=False)
    
    # 初始化Top_Flag列
    stats_df['Top_Flag'] = 'NA'
    
    if use_prefix:
        # 使用前缀分组的方式选择TOP基因
        top_prefixes = (
            stats_df.groupby('Prefix')['CSQ_Mean']
            .mean()
            .nlargest(top_n)
            .index
        )
        
        for prefix in top_prefixes:
            prefix_genes = stats_df[stats_df['Prefix'] == prefix]
            top_gene_row = prefix_genes.nlargest(1, 'CSQ_Mean')
            
            if not top_gene_row.empty:
                top_gene = top_gene_row.iloc[0]['Gene']
                stats_df.loc[
                    (stats_df['Prefix'] == prefix) & (stats_df['Gene'] == top_gene),
                    'Top_Flag'
                ] = 'Top'
    else:
        # 不使用前缀，直接提取TOP N个基因
        top_genes = stats_df.nlargest(top_n, 'CSQ_Mean')['Gene'].tolist()
        stats_df.loc[stats_df['Gene'].isin(top_genes), 'Top_Flag'] = 'Top'
    
    # 重新按Count排序
    stats_df = stats_df.sort_values('Count', ascending=False)
    
    return stats_df


def summarize_overall_data(df, output):
    """生成整体数据汇总（一行）
    
    Args:
        df: 输入数据框
        summary_cols: 需要汇总的列名列表，如果为None则汇总所有数值列（排除基因列）
    
    Returns:
        包含一行汇总数据的数据框
    """
    # 确定需要汇总的列
    numeric_cols = df.select_dtypes(include=['float64', 'int64']).columns.tolist()
        # 排除基因相关列
    exclude_patterns = ['Gene', 'gene']
    summary_cols = [col for col in numeric_cols 
                    if not any(pattern in col for pattern in exclude_patterns)]
    
    if not summary_cols:
        print("Warning: No numeric columns found to summarize")
        return pd.DataFrame()
    
    #print(f"\n汇总以下列: {', '.join(summary_cols)}")
    
    summary_info = {
        'BUSCO_id': output,
        'Total_Genes': len(set(df['GeneA'].dropna().tolist() + df['GeneB'].dropna().tolist())) if 'GeneA' in df.columns else 0
    }
    
    # 计算每列的均值
    for col in summary_cols:
        try:
            values = pd.to_numeric(df[col], errors='coerce').dropna()
            if len(values) > 0:
                summary_info[f'{col}'] = round(values.mean(), 4)
            else:
                summary_info[f'{col}'] = None
        except Exception as e:
            print(f"Warning: Error processing column {col}: {str(e)}")
            summary_info[f'{col}'] = None
    
    summary_df = pd.DataFrame([summary_info])
    
    return summary_df


def process_gene_data(input_data,  output, target_col="CSQ", use_prefix=True, top_n=18, 
                     visualize=False, do_summary=False):
    """主处理流程"""
    df = pd.read_csv(input_data, sep="\t", header=0)

    if df.empty:
        return pd.DataFrame(), None, None, None
    
    # 运行统计
    stats_df = calculate_prefix_stats(df, target_col)
    stats_df_gene = calculate_gene_stats(df, target_col, use_prefix=use_prefix, top_n=top_n)
    
    # 生成汇总文件（一行）
    summary_df = None
    if do_summary:
        summary_df = summarize_overall_data(df, output)
    
    fig = None

    if visualize and not stats_df.empty:
        try:
            fig, ax = plt.subplots(figsize=(10, 5))
            stats_df.plot.bar(x='Prefix', y='Count', ax=ax, color='skyblue')
            ax.set_title('Gene Prefix Distribution')
            ax.set_ylabel('Occurrence Count')
            plt.xticks(rotation=45)
            plt.tight_layout()
        except Exception as e:
            print(f"Failed to generate visualizations: {str(e)}")

    return stats_df, stats_df_gene, summary_df, fig


def main():
    parser = argparse.ArgumentParser(
        description='Gene prefix statistical analysis tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--input', required=True, help='Input file path')
    parser.add_argument('-o', '--output', help='Output file prefix', default='output')
    parser.add_argument('-t', '--target-col', type=str, default="CSQ",
                       help='Target data column name')
    parser.add_argument('--use-prefix', action='store_true', 
                       help='Group statistics by gene prefix (default: group directly by gene)')
    parser.add_argument('--top-n', type=int, default=10,
                       help='Extract the top N genes (default: 10)')
    parser.add_argument('--visualize', action='store_true', help='Generate visualizations')
    parser.add_argument('--summary', action='store_true', 
                       help='Generate an overall summary file (one row of means for all numeric columns)')
    args = parser.parse_args()

    try:
        stats_df, stats_df_gene, summary_df, fig = process_gene_data(
            args.input, args.output,
            target_col=args.target_col,
            use_prefix=args.use_prefix,
            top_n=args.top_n,
            visualize=args.visualize,
            do_summary=args.summary,
        )

        if stats_df.empty:
            print("No valid statistics generated. Possible reasons:")
            print("1. The input file is empty")
            print("2. No valid gene columns found")
            print("3. Invalid data in the target column")
            sys.exit(1)
        
        # 保存结果
        stats_df.to_csv(args.output + "_species_stat.tsv", index=False, sep="\t")
        print(f"Statistics saved to {args.output}_species_stat.tsv")

        stats_df_gene.to_csv(args.output + "_gene_stat.tsv", index=False, sep="\t")
        print(f"Statistics saved to {args.output}_gene_stat.tsv")
        
        # 保存汇总文件
        if summary_df is not None and not summary_df.empty:
            summary_df.to_csv(args.output + "_overall_summary.tsv", index=False, sep="\t")

     
        # 保存图表
        if fig:
            chart_path = args.output + '_prefix_distribution.png'
            fig.savefig(chart_path, dpi=300, bbox_inches='tight')
            print(f"\nVisualizations saved to {chart_path}")

    
        
    except Exception as e:
        print(f"Program execution error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
