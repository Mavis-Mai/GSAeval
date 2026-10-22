import pandas as pd
import re
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def merge_tables_by_geneid(
    file1_path: str,
    file2_path: str,
    AEDfile3: str,
    CDS_ratio_file: str,
    gfffile: str,
    output_path: str = None,
    merge_how: str = 'inner',
    drop_na_geneid: bool = True
    ) -> pd.DataFrame:
    """
    基于GeneID列合并两个表格文件
    
    参数:
        file1_path: str - 第一个表格文件路径
        file2_path: str - 第二个表格文件路径
        output_path: str - 可选，合并结果输出路径
        merge_how: str - 合并模式 ('inner', 'outer', 'left', 'right')
        drop_na_geneid: bool - 是否丢弃GeneID为NA的行
    
    返回:
        pd.DataFrame - 合并后的DataFrame
    """
    # 自动读取gene-transcript的关系
    def get_tx2gene_from_gtf(gtf_file):
        tx2gene_dict = {}
        gene_pattern = re.compile(r'gene_id\s+"([^"]+)"')
        transcript_pattern = re.compile(r'transcript_id\s+"([^"]+)"')
        f_open = open(gtf_file, 'r')
        with f_open as f:
            for line in f:
                if line.startswith('#'):
                    continue
                columns = line.strip().split('\t')
                if len(columns) < 9:
                    continue
                attributes = columns[8]
                if columns[2] != "transcript":
                    continue
                tx_match = transcript_pattern.search(attributes)
                if tx_match:
                    tx_id = tx_match.group(1)
                    if tx_id not in tx2gene_dict:
                        gene_match = gene_pattern.search(attributes)
                        if gene_match:
                            gene_id = gene_match.group(1)
                            tx2gene_dict[tx_id] = gene_id
        return tx2gene_dict


    # 自动检测分隔符读取文件
    def read_table(path):
        if path.endswith('.tsv') or path.endswith('.txt'):
            return pd.read_csv(path, sep='\t', dtype=str)
        elif path.endswith('.csv'):
            return pd.read_csv(path, sep=',', dtype=str)
        else:
            raise ValueError("Unsupported file format. Please use .tsv, .txt or .csv")
    
    # ------------------ 新增：用于记录分位数统计信息的字典 ------------------
    stats_data = {
        "Metric": [],
        "Quantile_0.90": [],
        "Quantile_0.95": []
    }
    
    def record_quantiles(metric_name, data_array):
        """计算并保存0.90和0.95的分位数"""
        series = pd.Series(data_array).dropna()
        q90 = series.quantile(0.90)
        q95 = series.quantile(0.95)
        stats_data["Metric"].append(metric_name)
        stats_data["Quantile_0.90"].append(q90)
        stats_data["Quantile_0.95"].append(q95)
        return q90, q95
    # --------------------------------------------------------------------------

    # 读取数据
    df1 = read_table(file1_path)
    df2 = read_table(file2_path)
    df3 = read_table(AEDfile3)
    df4 = read_table(CDS_ratio_file)
    
    if 'GeneID' not in df1.columns or 'GeneID' not in df2.columns:
        raise ValueError("Both tables must contain 'GeneID' column")
    
    if drop_na_geneid:
        df1 = df1[df1['GeneID'] != 'NA'].copy()
        df1 = df1[df1['GeneID'] != '.'].copy()
        df2 = df2[df2['GeneID'] != 'NA'].copy()
        df3 = df3[["Transcript1","Transcript2","AED"]]

    df1 = df1[["GeneID","Coverage_in_RNAseq","Gene_Overlap_Percentage","Significant_Overlap"]]
    merged = pd.merge(df1, df2, on='GeneID', how=merge_how, suffixes=('_1', '_2'))
    
    df3['AED'] = pd.to_numeric(df3['AED'], errors='coerce')
    merged_df = pd.merge(merged, df3, left_on='GeneID', right_on='Transcript1', how=merge_how, suffixes=('_1', '_2'))
    merged_df = pd.merge(merged_df, df4, on='GeneID', how=merge_how, suffixes=('_1', '_2'))
    
    # 提取基因ID信息：
    tx2gene_dict = get_tx2gene_from_gtf(gfffile)
    merged_df["gene"] = merged_df["GeneID"].map(tx2gene_dict)
    
    df3['AED'] = pd.to_numeric(df3['AED'], errors='coerce')
    merged_df['AED'].fillna(1.0, inplace=True) 
    print(merged_df.head())
    
    filtered_df = merged_df.loc[merged_df.groupby('gene')['AED'].idxmin()]

    filtered_df = filtered_df[["genename","Transcript1","Coverage_in_RNAseq","Gene_Overlap_Percentage","Significant_Overlap","Total_gene_sjs","Matched_sjs","Status","AED","CDS_ratio"]]
    filtered_df.columns = ["Gene","Transcript","Coverage_in_RNAseq","Gene_Overlap_Percentage","Significant_Overlap","Total_gene_sjs","Matched_sjs","SJ_level","CDS_distance","CDS_ratio"]
    
    # ---------------- 1. CDS_distance 画图及指标记录 ----------------
    plt.figure(figsize=(10, 6))
    sel_aed = filtered_df.groupby('Gene')['CDS_distance'].min().reset_index()
    sns.histplot(sel_aed['CDS_distance'], bins=50, kde=True, color='#ffae4c', stat='density')
    
    # 获取并记录分位数
    q90_cds_dist, q95_cds_dist = record_quantiles('CDS_distance', sel_aed['CDS_distance'])
    threshold_95_AED = q95_cds_dist # 原代码用于画图的其实是 0.95
    
    plt.axvline(threshold_95_AED, color='red', linestyle='--', label='AED ≤ 0.1 (High Quality)')
    plt.xlabel('CDS_distance')
    plt.ylabel('Density')
    plt.ylim(0,10)
    plt.text(
        x=threshold_95_AED + 0.01,
        y=plt.ylim()[1] * 0.9,
        s=f'Top 95%: {threshold_95_AED:.2f}',
        color='red',
        fontsize=12
        )
    if output_path:
        plt.savefig(f"{output_path}_CDS_distance_distribution.svg")
        plt.savefig(f"{output_path}_CDS_distance_distribution.png", dpi=300)

    # ---------------- 2. CDS_ratio 画图及指标记录 ----------------
    filtered_df['CDS_ratio'] = pd.to_numeric(filtered_df['CDS_ratio'], errors='coerce')
    
    plt.figure(figsize=(10, 6))
    sns.histplot(filtered_df['CDS_ratio']/100, bins=50, kde=True, color='#ffae4c', stat='density')
    
    # 获取并记录分位数
    q90_cds_ratio, q95_cds_ratio = record_quantiles('CDS_ratio', filtered_df['CDS_ratio']/100)
    threshold_95_cds_ratio = q95_cds_ratio
    
    plt.xlabel('CDS_ratio')
    plt.ylabel('Density')
    plt.axvline(threshold_95_cds_ratio, color='red', linestyle='--')
    plt.text(
        x=threshold_95_cds_ratio + 0.01,
        y=plt.ylim()[1] * 0.9,
        s=f'Top 95%: {threshold_95_cds_ratio:.2f}',
        color='red',
        fontsize=12
        )
    if output_path:
        plt.savefig(f"{output_path}_CDS_distribution.svg")
        plt.savefig(f"{output_path}_CDS_distribution.png", dpi=300)

    # 计算分数
    filtered_df.fillna(0)
    overlap_score = np.array([float(o.replace("%", "")) / 100 for o in filtered_df["Gene_Overlap_Percentage"]])  
    sj_score = np.where(
        filtered_df["Total_gene_sjs"].astype(float) != 0,
        (filtered_df["Total_gene_sjs"].astype(float) - filtered_df["Matched_sjs"].astype(float)) / filtered_df["Total_gene_sjs"].astype(float),
        np.nan
    )

    cds_score = filtered_df["CDS_ratio"]
    cdsdistance_score = filtered_df["CDS_distance"]
    cds_score_adjust = np.where(cds_score == 100, 0.8, cds_score/100)
                        
    filtered_df["Risk_score"] = np.where(
        ~pd.isna(sj_score),
        ((1 - overlap_score) + sj_score*2 + (1-cds_score_adjust) + cdsdistance_score*2) / 6,
        np.nan
    )

    filtered_df.loc[filtered_df["Coverage_in_RNAseq"].astype(float) < 1] = np.nan
                
    # ---------------- 3. Risk_score 画图及指标记录 ----------------
    filtered_df['Risk_score'] = pd.to_numeric(filtered_df['Risk_score'], errors='coerce')
    sel_rs = filtered_df.groupby('Gene')['Risk_score'].min().reset_index()
    
    plt.figure(figsize=(10, 6))
    sns.histplot(sel_rs['Risk_score'], bins=50, kde=True, color='#ffae4c', stat='density')
    
    # 获取并记录分位数
    q90_rs, q95_rs = record_quantiles('Risk_score', sel_rs['Risk_score'])
    threshold_95_rs = q95_rs
    
    plt.xlabel('Risk score')
    plt.ylabel('Density')
    plt.axvline(threshold_95_rs, color='red', linestyle='--')
    plt.text(
        x=threshold_95_rs + 0.01,
        y=plt.ylim()[1] * 0.9,
        s=f'Top 95%: {threshold_95_rs:.2f}',
        color='red',
        fontsize=12
    )
    if output_path:
        plt.savefig(f"{output_path}_riskscore_distribution.svg")
        plt.savefig(f"{output_path}_riskscore_distribution.png", dpi=300)

    # ---------------- 4. SJ_score 画图及指标记录 ----------------
    plt.figure(figsize=(10, 6))
    sj_score_array = 1 - sj_score  # 画图使用的是 1 - sj_score
    sns.histplot(sj_score_array, bins=50, kde=True, color='#ffae4c', stat='density')
    
    # 获取并记录分位数
    q90_sj, q95_sj = record_quantiles('SJ_score (1 - sj_score)', sj_score_array)
    threshold_95_sj = q95_sj
    
    plt.xlabel('sj score')
    plt.ylabel('Density')
    plt.axvline(threshold_95_sj, color='red', linestyle='--')
    plt.text(
        x=threshold_95_sj + 0.01,
        y=plt.ylim()[1] * 0.9,
        s=f'Top 95%: {threshold_95_sj:.2f}',
        color='red',
        fontsize=12
    )
    if output_path:
        plt.savefig(f"{output_path}_sjscore_distribution.svg")
        plt.savefig(f"{output_path}_sjscore_distribution.png", dpi=300)

    # ---------------- 5. Overlap_score 画图及指标记录 ----------------
    plt.figure(figsize=(10, 6))
    sns.histplot(overlap_score, bins=50, kde=True, color='#ffae4c', stat='density')
    
    # 获取并记录分位数
    q90_ov, q95_ov = record_quantiles('Overlap_score', overlap_score)
    threshold_95_overlap_score = q95_ov
    
    plt.xlabel('overlap ratio')
    plt.ylabel('Density')
    plt.axvline(threshold_95_overlap_score, color='red', linestyle='--')
    plt.text(
        x=threshold_95_overlap_score + 0.01,
        y=plt.ylim()[1] * 0.9,
        s=f'Top 95%: {threshold_95_overlap_score:.2f}',
        color='red',
        fontsize=12
    )
    if output_path:
        plt.savefig(f"{output_path}_overlapscore_distribution.svg")
        plt.savefig(f"{output_path}_overlapscore_distribution.png", dpi=300)

    # ================= 输出统计信息文件 =================
    if output_path:
        stats_df = pd.DataFrame(stats_data)
        stats_file_path = f"{output_path}_quantiles_statistics.tsv"
        stats_df.to_csv(stats_file_path, sep='\t', index=False)
        print(f"Quantile statistics for all metrics recorded and saved to: {stats_file_path}")
        print(stats_df)

    return filtered_df
