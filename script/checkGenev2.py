#!/usr/bin/env python3
"""
合并版 BAM/基因结构 可视化脚本
功能：输入特定 gene_id 或 transcript_id，自动提取区域(含可调边界)、提取BAM、解析GFF/GTF并按层绘制可视化图。
"""
import argparse
import os
import sys
import pysam
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from collections import defaultdict

def parse_attributes(attributes):
    """兼容解析 GTF 的 key "value" 和 GFF3 的 key=value 属性。"""
    attr_dict = {}
    for attr in attributes.split(';'):
        attr = attr.strip()
        if not attr:
            continue
        if '=' in attr:
            key, val = attr.split('=', 1)
        else:
            kv = attr.split(None, 1)
            if len(kv) != 2:
                continue
            key, val = kv
        attr_dict[key] = val.strip().strip('"')
    return attr_dict

def parse_gxf(ref_gxf):
    """解析GTF/GFF文件，建立 gene_id 和 transcript_id 坐标索引。"""
    gene_dict = defaultdict(list)
    transcript_dict = {}

    with open(ref_gxf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, _, strand, _, attributes = fields

            attr_dict = parse_attributes(attributes)

            gene_id = attr_dict.get("gene_id") or attr_dict.get("Parent") or attr_dict.get("geneID")
            if not gene_id:
                gene_id = "NA"

            transcript_id = attr_dict.get("transcript_id")
            if feature.lower() in ('mrna', 'transcript'):
                transcript_id = transcript_id or attr_dict.get("ID")
            elif not transcript_id:
                transcript_id = attr_dict.get("Parent")
            if not transcript_id:
                continue

            record = {
                "gene_id": gene_id,
                "transcript_id": transcript_id,
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "strand": strand
            }

            # gene 模式继续沿用转录本/mRNA条目，避免改变原有行为。
            if feature.lower() in ('mrna', 'transcript'):
                gene_dict[gene_id].append(record)

            # transcript 模式聚合所有带 transcript_id 的条目；即使 GTF 没有
            # transcript 行，也能由 exon/CDS 得到完整区间。
            for tid in transcript_id.split(','):
                tid = tid.strip()
                if not tid:
                    continue
                if tid not in transcript_dict:
                    transcript_dict[tid] = dict(record, transcript_id=tid)
                else:
                    old = transcript_dict[tid]
                    if old["chrom"] == chrom:
                        old["start"] = min(old["start"], int(start))
                        old["end"] = max(old["end"], int(end))
    return gene_dict, transcript_dict

def getGenomeRegionbyid(gene_id, gene_dict, extend=1000):
    """根据 gene_id 获取该基因的整体作图区域（合并所有转录本并向外扩展边界）"""
    if gene_id not in gene_dict or not gene_dict[gene_id]:
        return []
    
    # 取最远的起始点和终止点包含整个基因的全部转录本
    chrom = gene_dict[gene_id][0]["chrom"]
    strand = gene_dict[gene_id][0]["strand"]
    min_start = min(t["start"] for t in gene_dict[gene_id])
    max_end = max(t["end"] for t in gene_dict[gene_id])
    
    # 返回包含扩展边界的坐标（起始位置保证不小于1）
    return [(chrom, max(1, min_start - extend), max_end + extend, strand)]

def getGenomeRegionbyTranscript(transcript_id, transcript_dict, extend=1000):
    """根据 transcript_id 获取该转录本区间并向外扩展边界。"""
    transcript = transcript_dict.get(transcript_id)
    if not transcript:
        return []
    return [(
        transcript["chrom"],
        max(1, transcript["start"] - extend),
        transcript["end"] + extend,
        transcript["strand"]
    )]

def gettargetGXFbygff(region_list, ref_gxf, output_gff="cds_ref.gff3.tmp", transcript_id=None):
    """切取目标区域内相关的GFF注释以便于精简作图"""
    with open(ref_gxf, 'r') as fin, open(output_gff, 'w') as fout:
        fout.write("##gff-version 3\n")
        for line in fin:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom, _, feature, start, end, _, strand, _, attributes = fields
            start, end = int(start), int(end)

            if transcript_id:
                attr_dict = parse_attributes(attributes)
                line_transcript_id = attr_dict.get("transcript_id")
                if feature.lower() in ('mrna', 'transcript'):
                    line_transcript_id = line_transcript_id or attr_dict.get("ID")
                parent_ids = (attr_dict.get("Parent") or "").split(',')
                if line_transcript_id != transcript_id and transcript_id not in parent_ids:
                    continue

            for reg_chrom, reg_start, reg_end, reg_strand in region_list:
                if chrom == reg_chrom and start <= reg_end and end >= reg_start:
                    fout.write(line)
                    break
    return output_gff

def getbam(bamfile, region_list, output_bam="target_region.sorted.bam"):
    """使用pysam提取指定区域的bam文件及建立索引"""
    with pysam.AlignmentFile(bamfile, "rb") as in_bam:
        with pysam.AlignmentFile(output_bam, "wb", header=in_bam.header) as out_bam:
            for region in region_list:
                chrom, start, end, strand = region 
                try:
                    for read in in_bam.fetch(chrom, max(0, start - 1), end):
                        out_bam.write(read)
                except ValueError as e:
                    print(f"Unable to extract region {region}: {str(e)}")
  
    pysam.index(output_bam)
    print(f"Successfully extracted reads from the specified region to {output_bam}")
    return output_bam

def parse_gff3(gff_file, chrom, start, end):
    """解析临时GFF3文件，按转录本分组提取外显子和CDS特征"""
    transcripts = {}
  
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            curr_chrom, feat_type = parts[0], parts[2]
            if curr_chrom != chrom: continue
          
            attributes = parse_attributes(parts[8])

            transcript_id = attributes.get("transcript_id") or attributes.get("ID")
            if feat_type.lower() in ('exon', 'cds'):
                # GTF 子特征用 transcript_id 归组；GFF3 子特征用 Parent 归组。
                parent_id = attributes.get("transcript_id") or attributes.get("Parent")
            else:
                parent_id = attributes.get("Parent") or attributes.get("geneID") or attributes.get("gene_id")
          
            if feat_type.lower() in ('mrna', 'transcript'):
                transcripts[transcript_id] = {
                    'exon': [],
                    'cds': [],
                    'strand': parts[6],
                    'geneID': attributes.get('gene_id', parent_id)
                }
            elif feat_type.lower() in ('exon', 'cds'):
                feat_start, feat_end = int(parts[3]), int(parts[4])
                if feat_end < start or feat_start > end: continue

                # 部分 GTF 只有 exon/CDS 行，没有单独的 transcript 行。
                if parent_id and parent_id not in transcripts:
                    transcripts[parent_id] = {
                        'exon': [],
                        'cds': [],
                        'strand': parts[6],
                        'geneID': attributes.get('gene_id')
                    }
                if parent_id:
                    transcripts[parent_id][feat_type.lower()].append(
                        (max(feat_start, start), min(feat_end, end))
                    )
  
    for transcript in transcripts.values():
        for key in ['exon', 'cds']:
            transcript[key].sort(key=lambda x: x[0])
    return transcripts

def assign_reads_to_tracks(reads, min_gap=5):
    """分配reads到轨道 (从顶部开始)"""
    tracks = [[]]
    for read in reads:
        placed = False
        for i in range(len(tracks)):
            if not tracks[i] or (read.reference_start > tracks[i][-1].reference_end + min_gap):
                tracks[i].append(read)
                placed = True
                break
        if not placed:
            tracks.append([read])
    return tracks 

def plot_reads(ax, read, y_pos, read_height=0.4):
    """绘制单个read (颜色区分链方向)"""
    color = '#ebafaf' if read.is_reverse else '#afafeb'
    blocks = read.get_blocks()
  
    for start, end in blocks:
        ax.plot([start, end], [y_pos, y_pos], lw=5, color=color, alpha=0.7, solid_capstyle='butt')
  
    if len(blocks) > 1:
        for i in range(len(blocks)-1):
            ax.plot([blocks[i][1], blocks[i+1][0]], [y_pos, y_pos], lw=.5, linestyle=':', color='gray')

def plot_gene_structure(ax, features, y_center, color, label):
    """绘制基因结构 (带链方向和内含子连接线)"""
    strand = features.get('strand', '+')
    exons = features['exon']
    exon_height = 0.1

    for s, e in exons:
        ax.add_patch(plt.Rectangle((s, y_center-exon_height/2), e-s, exon_height, color=color, linewidth=0))
  
    for s, e in features.get('cds', []):
        ax.add_patch(plt.Rectangle((s, y_center-exon_height), e-s, exon_height*2, color=color, linewidth=0))
  
    if len(exons) > 1:
        for i in range(len(exons)-1):
            intron_start, intron_end = exons[i][1], exons[i+1][0]
            ax.plot([intron_start, intron_end], [y_center, y_center], color='gray', linestyle=':', linewidth=1, alpha=0.7)
            
            arrow_positions = np.linspace(intron_start, intron_end, 3)[1:-1]
            for pos in arrow_positions:
                arrow_char = '>' if strand == '+' else '<'
                ax.text(pos, y_center, arrow_char, ha='center', va='center', color='black', fontsize=10)

    ax.text(0.01, 0.9, label, transform=ax.transAxes, color='black', fontsize=8, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

def plot_bam(bam_file, region_list, gff_file1, output):
    """原汁原味地保留原有画图格式：基因结构 -> 覆盖度 -> reads 堆叠"""
    for chrom, start, end, _ in region_list:
        # 1. 初始化3层画布
        fig = plt.figure(figsize=(15, 5))
        gs = gridspec.GridSpec(3, 1, height_ratios=[1, 2, 3])
        ax_gene1 = fig.add_subplot(gs[0])  
        ax_cov = fig.add_subplot(gs[1])    
        ax_read = fig.add_subplot(gs[2])   

        # 2. 绘制基因结构
        y_pos = 0.5
        y_step = 0.3 
        gene1 = parse_gff3(gff_file1, chrom, start, end)
        for i, (transcript_id, features) in enumerate(gene1.items()):
            plot_gene_structure(ax_gene1, features, y_center=y_pos, color='#0319a9', label=transcript_id)
            y_pos += y_step 
        ax_gene1.set_xlim(start, end)
        ax_gene1.axis('off')

        bam = pysam.AlignmentFile(bam_file, "rb")
        
        # 3. 计算并绘制覆盖度
        coverage = bam.count_coverage(chrom, start, end, read_callback='all', quality_threshold=0)
        y_cov = np.sum(coverage, axis=0)
        x_cov = np.arange(start, end)

        ax_cov.fill_between(x_cov, y_cov, color='gray', alpha=0.6)
        ax_cov.set_ylabel("Coverage")
        ax_cov.set_xlim(start, end)
        ax_cov.set_xticks([])
      
        # 4. 绘制reads (从顶部向下)
        reads = [r for r in bam.fetch(chrom, max(0, start-1), end) if not r.is_unmapped]
        read_height = 0.5
        if reads:
            tracks = assign_reads_to_tracks(sorted(reads, key=lambda x: x.reference_start))
            total_tracks = len(tracks)

            for track_idx, track_reads in enumerate(tracks):
                y_pos = (total_tracks - track_idx - 1) * read_height  
                for read in track_reads:
                    plot_reads(ax_read, read, y_pos, read_height)
            ax_read.set_ylim(-read_height, total_tracks * read_height)
        else:
            ax_read.set_ylim(0, 1)

        ax_read.set_xlabel(f"Position on {chrom}")
        ax_read.grid(axis='x', linestyle='--', alpha=0.3)
        ax_read.set_xlim(start, end)
        ax_read.set_yticks([])
        plt.subplots_adjust(left=0.05, right=0.95)

        plt.tight_layout()
        output_file = f'{output}.png'
        plt.savefig(output_file, dpi=200, bbox_inches='tight')
        print(f"Figure saved to: {output_file}")


def main(args_list=None):
    parser = argparse.ArgumentParser(
        description="Gene BAM Visualization Script - Extracts region by gene_id and plots structural features, coverage, and reads.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    target_group = parser.add_mutually_exclusive_group(required=True)
    target_group.add_argument("-g", "--gene_id", help="Target gene ID (e.g. AT1G01010)")
    target_group.add_argument("-t", "--transcript_id", help="Target transcript ID in GTF/GFF")
    parser.add_argument("-b", "--bam_file", required=True, help="Input sorted and indexed BAM file")
    parser.add_argument("-r", "--ref_gff", required=True, help="Reference GFF/GTF file")
    parser.add_argument("-e", "--extend", type=int, default=1000, help="Custom boundary extension (bp) around the gene")
    parser.add_argument("-o", "--output", default="output", help="Output figure prefix")
  
    #args = parser.parse_args()
    if args_list:
        args, _ = parser.parse_known_args(args_list)
    else:
        args = parser.parse_args()

    # 1. 从 GTF/GFF 获取目标的坐标
    target_type = "transcript_id" if args.transcript_id else "gene_id"
    target_id = args.transcript_id or args.gene_id
    print(f"Searching in {args.ref_gff} for coordinates of {target_type}: {target_id}...")
    gene_dict, transcript_dict = parse_gxf(args.ref_gff)
    if args.transcript_id:
        region_list = getGenomeRegionbyTranscript(args.transcript_id, transcript_dict, extend=args.extend)
    else:
        region_list = getGenomeRegionbyid(args.gene_id, gene_dict, extend=args.extend)
    
    if not region_list:
        print(f"Error: In {args.ref_gff}, could not find {target_type}: '{target_id}'. Please check the name or file.")
        sys.exit(1)
        
    print(f"Plot region found: {region_list[0]}")

    # 2. 截取临时区域 GFF 文件
    tmp_gff = "cds_ref_tmp.gff3"
    gettargetGXFbygff(
        region_list,
        args.ref_gff,
        output_gff=tmp_gff,
        transcript_id=args.transcript_id
    )

    # 3. 截取临时 BAM 文件
    tmp_bam = "target_region.sorted.bam"
    getbam(args.bam_file, region_list, output_bam=tmp_bam)

    # 4. 绘图 (严格保留原有结构)
    print("Generating visualizations...")
    plot_bam(
        bam_file=tmp_bam,
        region_list=region_list,
        gff_file1=tmp_gff,
        output=args.output
    )

    # 5. 清理临时文件
    for f in [tmp_gff, tmp_bam, tmp_bam + ".bai"]:
        if os.path.exists(f):
            os.remove(f)

if __name__ == "__main__":
    main()
