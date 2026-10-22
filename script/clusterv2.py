#!/usr/bin/env python3
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import argparse
from collections import defaultdict
import sys

def load_mash_matrix(dist_file):
    """加载Mash距离矩阵并验证完整性"""
    with open(dist_file) as f:
        lines = [line.strip() for line in f if line.strip()]

    samples = []
    dist_values = []

    # 首行应为样本数
    expected_samples = int(lines[0])
    for line in lines[1:]:
        parts = line.split()
        samples.append(parts[0])
        dist_values.extend([float(x) for x in parts[1:]])

    # 验证数据完整性
    if len(samples) != expected_samples:
        raise ValueError(f"The number of samples does not match. The expected number was {expected_samples}, but the actual number was {len(samples)}")
    
    dist_matrix = squareform(dist_values)
    np.fill_diagonal(dist_matrix, 0)  # 确保对角线为0
    
    return samples, dist_matrix

def find_representative_samples(cluster_df, dist_matrix):
    """为每个簇寻找代表性样品（中心点）"""
    cluster_dict = defaultdict(list)
    for idx, row in cluster_df.iterrows():
        cluster_dict[row['Cluster']].append(idx)
    
    representatives = {}
    intra_cluster_dists = {}
    
    for cluster, members in cluster_dict.items():
        submatrix = dist_matrix[np.ix_(members, members)]
        # 计算每个样品与簇内其他样品的平均距离
        avg_dists = submatrix.mean(axis=1)
        # 选择平均距离最小的作为代表
        rep_idx = members[np.argmin(avg_dists)]
        representatives[cluster] = {
            'sample': cluster_df.iloc[rep_idx]['Sample'],
            'index': rep_idx,
            'avg_distance': np.min(avg_dists),
            'cluster_size': len(members)
        }
        intra_cluster_dists[cluster] = avg_dists.mean()
    
    return representatives, intra_cluster_dists

def hierarchical_clustering(samples, dist_matrix, output_prefix, threshold=None):
    """
    执行层次聚类并生成完整报告
    """
    condensed_dist = squareform(dist_matrix)
    Z = linkage(condensed_dist, method='average')

    # 自动计算阈值（使用中位数的70%）
    if threshold is None:
        threshold = np.percentile(condensed_dist, 50)*0.9 
        print(f"\nAutomatically calculate the clustering threshold: {threshold:.4f} (distance median)*0.9")

    clusters = fcluster(Z, t=threshold, criterion='distance')
    cluster_df = pd.DataFrame({'Sample': samples, 'Cluster': clusters})
    
    # 寻找代表性样品
    representatives, intra_dists = find_representative_samples(cluster_df, dist_matrix)
    
    # 生成详细报告
    report = []
    for cluster in sorted(set(clusters)):
        rep_info = representatives[cluster]
        report.append({
            'Cluster': cluster,
            'Representative': rep_info['sample'],
            'RepresentativeAvgDist': f"{rep_info['avg_distance']:.4f}",
            'ClusterSize': rep_info['cluster_size'],
            'IntraClusterAvgDist': f"{intra_dists[cluster]:.4f}"
        })
    
    report_df = pd.DataFrame(report)
    report_df.to_csv(f"{output_prefix}_report.csv", index=False)
    
    # 在聚类结果中标记代表性样品
    cluster_df['IsRepresentative'] = cluster_df.apply(
        lambda r: r['Sample'] == representatives[r['Cluster']]['sample'], axis=1)
    cluster_df.to_csv(f"{output_prefix}_clusters.csv", index=False)

    # 获取代表性样品的索引
    rep_indices = [rep['index'] for rep in representatives.values()]
    
    # 绘制树状图
    plt.figure(figsize=(12, max(6, len(samples)*0.3)))
    dendro = dendrogram(Z, labels=samples, orientation='right', leaf_font_size=8)
    
    # 创建颜色映射字典
    color_dict = {idx: 'red' if idx in rep_indices else 'black' for idx in range(len(samples))}
    
    # 手动设置叶节点颜色 - 修正后的方法
    ax = plt.gca()
    xticks = ax.get_xticks()
    yticks = ax.get_yticks()
    
    for i, (label, ytick) in enumerate(zip(ax.get_yticklabels(), yticks)):
        if i < len(dendro['leaves']):
            idx = dendro['leaves'][i]
            label.set_color(color_dict[idx])
            if idx in rep_indices:
                label.set_fontweight('bold')
    
    plt.axvline(x=threshold, color='r', linestyle='--', 
               label=f'Cutoff: {threshold:.4f}')
    plt.legend()
    plt.title(f"Hierarchical Clustering (n={len(samples)})")
    plt.xlabel("Distance")
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_dendrogram.png", bbox_inches='tight', dpi=300)
    plt.close()

    return cluster_df, report_df


def print_cluster_stats(cluster_df, report_df):
    """打印美观的聚类统计信息"""
    from tabulate import tabulate
    
    print("\n📊 Cluster statistics:")
    print(tabulate(report_df.rename(columns={
        'Cluster': 'Cluster ID',
        'Representative': 'Representative samples',
        'RepresentativeAvgDist': 'Representative average distance',
        'ClusterSize': 'Cluster size',
        'IntraClusterAvgDist': 'Average distance within the cluster'
    }), headers='keys', tablefmt='grid', showindex=False))
    
    size_dist = cluster_df['Cluster'].value_counts().value_counts().sort_index()
    print("\n📈 Cluster size distribution:")
    print(tabulate(
        [(f"{size}Samples", count) for size, count in size_dist.items()],
        headers=['Cluster size', 'Frequency'],
        tablefmt='grid'
    ))

def main(args_list=None):
    """主程序调用入口"""
    parser = argparse.ArgumentParser(
        description='Hierarchical clustering analysis tool based on Mash distance',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-i", "--input", required=True,
                       help="input the distance matrix file (Mash triangle format)")
    parser.add_argument("-o", "--output", default="mash_clustering",
                       help="output prefix")
    parser.add_argument("-d", "--threshold", type=float,
                       help="Clustering distance threshold (0-1), if not set, it will be calculated automatically")
    parser.add_argument("--tablefmt", default="grid",
                       choices=["plain", "simple", "grid", "fancy_grid"],
                       help="Output format of statistical tables")
    
    args = parser.parse_args(args_list) if args_list else parser.parse_args()

    try:
        print("\n🔍 The distance matrix is being loaded...")
        samples, dist_matrix = load_mash_matrix(args.input)
        print(f"The distance matrix of {len(samples)} samples was successfully loaded")
        clusters, report = hierarchical_clustering(
            samples=samples,
            dist_matrix=dist_matrix,
            output_prefix=args.output,
            threshold=args.threshold
        )

        print_cluster_stats(clusters, report)
        
        print(f"\n💾 The result file has been generated:")
        print(f"  - Clustering details: {args.output}_clusters.csv")
        print(f"  - Analysis Report: {args.output}_report.csv")
        print(f"  - Tree diagram: {args.output}_dendrogram.png")
        
        return True
        
    except Exception as e:
        print(f"\n❌ Error: {str(e)}", file=sys.stderr)
        raise

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
