import os
import argparse
import gffutils
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq

def run_analysis(fasta_file, gff_file, output_prefix="eval_result"):
    print(f"1. Parsing annotation file: {gff_file} (May take a few minutes for the first run to build DB)...")
    db_name = f"{gff_file}.db"
    
    # 强制清理旧缓存，避免被残留错误数据干扰
    if os.path.exists(db_name):
        print("   Removing old database cache...")
        os.remove(db_name)
        
    db = gffutils.create_db(gff_file, db_name, merge_strategy="create_unique", force=True, 
                            keep_order=True)
        
    print(f"   Loading genome sequence: {fasta_file} ...")
    try:
        genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    except Exception as e:
        print(f"Error loading FASTA: {e}")
        return
        
    # [新增的安全检查] 检查 GFF 中的染色体名是否和 FASTA 匹配
    fasta_seqids = set(genome_dict.keys())
    db_seqids = {f.seqid for f in db.all_features() if f.seqid}
    matched = db_seqids.intersection(fasta_seqids)
    if not matched:
        print("\n🚨 Warning: FASTA and GFF3 chromosome names do not match at all! Subsequent codon/intron checks may fail.")
        print(f"   Example names in FASTA: {list(fasta_seqids)[:3]}")
        print(f"   Example names in GFF3: {list(db_seqids)[:3]}\n")

    results = []
    total_classic_introns = 0
    total_non_classic_introns = 0
    
    print("2. Extracting features and calculating metrics...")
    
    # 【核心修改区】：不再遍历 'gene'，直接收集所有的 'mRNA' 或 'transcript'
    transcripts = list(db.features_of_type(['mRNA', 'transcript']))
    if not transcripts:
        print("🚨 Error: No 'mRNA' or 'transcript' features found in GFF3! Check the feature names in column 3.")
        print(f"Features present in the current database: {list(db.featuretypes())}")
        return
        
    for tx in transcripts:
        tx_id = tx.id
        
        # 【核心修改区】：从转录本的属性中反推所属的 Gene ID
        parents = tx.attributes.get('Parent', [])
        gene_id = parents[0] if parents else f"UnknownGene_{tx_id}"
        
        # ---------------- 1. Start/Stop Codon Intactness ----------------
        cds_feats = list(db.children(tx, featuretype='CDS', order_by='start'))
        start_stop_tag = 'C' # Default: Neither
        if cds_feats:
            cds_seq = ""
            for cds in cds_feats:
                if cds.seqid in genome_dict:
                    seq_frag = genome_dict[cds.seqid].seq[cds.start-1 : cds.end]
                    cds_seq += str(seq_frag)
            
            if tx.strand == '-':
                cds_seq = str(Seq(cds_seq).reverse_complement())
            
            if len(cds_seq) >= 3:
                has_start = cds_seq.upper().startswith('ATG')
                has_stop = cds_seq[-3:].upper() in ['TAA', 'TAG', 'TGA']
                if has_start and has_stop:
                    start_stop_tag = 'A'
                elif has_start or has_stop:
                    start_stop_tag = 'B'
    
        # ---------------- 2. Exons, Introns & Splice Sites ----------------
        exons = list(db.children(tx, featuretype='exon', order_by='start'))
        exon_count = len(exons)
        tx_length = sum([e.end - e.start + 1 for e in exons])
        
        avg_intron_length = 0
        has_non_classic_splice = "NA"
        
        if exon_count > 1:
            intron_lengths = []
            non_classic_flag = False
            
            for i in range(exon_count - 1):
                intron_start = exons[i].end + 1
                intron_end = exons[i+1].start - 1
                i_len = intron_end - intron_start + 1
                
                if i_len > 0:
                    intron_lengths.append(i_len)
                    
                    if exons[i].seqid in genome_dict and i_len >= 4:
                        intron_seq = str(genome_dict[exons[i].seqid].seq[intron_start-1 : intron_end])
                        if tx.strand == '-':
                            intron_seq = str(Seq(intron_seq).reverse_complement())
                        
                        if intron_seq[:2].upper() == 'GT' and intron_seq[-2:].upper() == 'AG':
                            total_classic_introns += 1
                        else:
                            total_non_classic_introns += 1
                            non_classic_flag = True
                    else:
                        total_non_classic_introns += 1
                        non_classic_flag = True

            if intron_lengths:
                avg_intron_length = sum(intron_lengths) / len(intron_lengths)
            
            has_non_classic_splice = "Yes" if non_classic_flag else "No"

        # ---------------- 3. UTR Ratio ----------------
        utr_feats = list(db.children(tx, featuretype=['five_prime_UTR', 'three_prime_UTR']))
        if not utr_feats:
            cds_length = sum([c.end - c.start + 1 for c in cds_feats])
            utr_length = tx_length - cds_length if tx_length > cds_length else 0
        else:
            utr_length = sum([u.end - u.start + 1 for u in utr_feats])
            
        utr_ratio = utr_length / tx_length if tx_length > 0 else 0

        # ---------------- 4. Overlap with Adjacent Genes ----------------
        is_overlap = "No"
        # 【核心修改区】：改为通过提取附近的 transcript 来判断是否属于不同的基因重叠
        overlapping_txs = list(db.region(region=(tx.seqid, tx.start, tx.end), featuretype=['mRNA', 'transcript']))
        for ot in overlapping_txs:
            ot_parents = ot.attributes.get('Parent', [])
            ot_gene = ot_parents[0] if ot_parents else "UnknownGene"
            if ot_gene == gene_id: 
                continue # 属于同一个基因的转录本不计入“重叠基因”
            
            overlap_len = min(tx.end, ot.end) - max(tx.start, ot.start) + 1
            if overlap_len > 0:
                ot_len = ot.end - ot.start + 1
                if overlap_len > 1000 or (overlap_len / ot_len) >= 0.5:
                    is_overlap = "Yes"
                    break

        # ---------------- Compile Data ----------------
        results.append({
            "Transcript_ID": tx_id,
            "Gene_ID": gene_id,
            "Start_Stop_Codon_Tag": start_stop_tag,
            "UTR_Ratio": round(utr_ratio, 4),
            "Exon_Count": exon_count,
            "Avg_Intron_Length": round(avg_intron_length, 2),
            "Has_Non_Classic_Splice_Site": has_non_classic_splice,
            "Is_Overlapping_Gene": is_overlap
        })
        
    print(f"   Successfully extracted {len(results)} transcript records.")
    
    # [新增的安全检查] 彻底拦截空数据触发 KeyError
    if not results:
        print("🚨 Fatal error: No results collected. Aborting plotting.")
        return

    df = pd.DataFrame(results)

    print("3. Generating plots...")
    
    # Plot 1: Codon Completeness
    plt.figure(figsize=(8, 6))
    codon_counts = df['Start_Stop_Codon_Tag'].value_counts()
    plt.pie(codon_counts, labels=[f"Type {k} ({v})" for k, v in codon_counts.items()], 
            autopct='%1.1f%%', colors=['#66b3ff','#99ff99','#ffcc99'])
    plt.title("Start/Stop Codon Completeness\n(A: Both, B: One, C: Neither)")
    plt.savefig(f"{output_prefix}_1_codon_completeness.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 2: Splice Sites (Classic vs Non-Classic)
    plt.figure(figsize=(8, 6))
    splice_counts = [total_classic_introns, total_non_classic_introns]
    splice_labels = [f"Classic GT-AG ({splice_counts[0]})", f"Non-Classic ({splice_counts[1]})"]
    if sum(splice_counts) > 0:
        plt.pie(splice_counts, labels=splice_labels, autopct='%1.1f%%', colors=['#77dd77','#ffb347'])
        plt.title("Intron Splice Sites (GT-AG vs Others)")
        plt.savefig(f"{output_prefix}_2_splice_sites.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 3a: UTR Presence
    plt.figure(figsize=(8, 6))
    utr_presence = ["With UTR" if r > 0 else "Without UTR" for r in df['UTR_Ratio']]
    utr_counts = pd.Series(utr_presence).value_counts()
    plt.pie(utr_counts, labels=[f"{k} ({v})" for k, v in utr_counts.items()], 
            autopct='%1.1f%%', colors=['#c2c2f0','#ffb3e6'])
    plt.title("Proportion of Transcripts with UTR")
    plt.savefig(f"{output_prefix}_3a_utr_presence.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 3b: UTR Ratio Density
    plt.figure(figsize=(8, 6))
    utr_data = df[df['UTR_Ratio'] > 0]['UTR_Ratio']
    if len(utr_data) > 1 and utr_data.var() > 0:
        sns.kdeplot(utr_data, fill=True, color="blue")
        plt.title("Density Distribution of UTR Ratio")
        plt.xlabel("UTR Ratio (UTR length / Transcript length)")
        plt.ylabel("Density")
        plt.savefig(f"{output_prefix}_3b_utr_ratio_density.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 4: Overlapping Genes
    plt.figure(figsize=(8, 6))
    overlap_counts = df['Is_Overlapping_Gene'].map({"Yes": "Overlapped", "No": "Not Overlapped"}).value_counts()
    plt.pie(overlap_counts, labels=[f"{k} ({v})" for k, v in overlap_counts.items()], 
            autopct='%1.1f%%', colors=['#ff9999','#99ffcc'])
    plt.title("Gene Overlap Proportion\n(Threshold: >1kb or >50% of adjacent gene)")
    plt.savefig(f"{output_prefix}_4_overlap_genes.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 5: Exon Count Distribution 
    plt.figure(figsize=(10, 6))
    max_exon = df['Exon_Count'].max()
    bins = [1, 2, 7, 12, 17, 22] + list(range(27, max_exon + 6, 5))
    labels = ["1"] + [f"{bins[i]}-{bins[i+1]-1}" for i in range(1, len(bins)-1)]
    binned_exons = pd.cut(df['Exon_Count'], bins=bins, labels=labels, right=False)
    binned_counts = binned_exons.value_counts(sort=False)
    
    binned_counts = binned_counts[binned_counts > 0]
    sns.barplot(x=binned_counts.index, y=binned_counts.values, color="skyblue")
    plt.title("Exon Count Distribution")
    plt.xlabel("Exon Count Groups")
    plt.ylabel("Frequency")
    
    single_exon_ratio = (df['Exon_Count'] == 1).mean() * 100
    plt.text(0.5, max(binned_counts.values)*0.9, f"Single-exon gene ratio: {single_exon_ratio:.2f}%", 
             fontsize=12, color='red', weight='bold')
    plt.xticks(rotation=45)
    plt.savefig(f"{output_prefix}_5_exon_distribution.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 6: Intron Length Density
    plt.figure(figsize=(8, 6))
    intron_data = df[df['Exon_Count'] > 1]['Avg_Intron_Length']
    if len(intron_data) > 1 and intron_data.var() > 0:
        q99 = intron_data.quantile(0.99)
        sns.kdeplot(intron_data[intron_data < q99], fill=True, color="purple")
        plt.title("Density Distribution of Average Intron Length (Top 1% excluded)")
        plt.xlabel("Average Intron Length (bp)")
        plt.ylabel("Density")
        plt.savefig(f"{output_prefix}_6_intron_density.png", dpi=300, bbox_inches='tight')
    plt.close()

    # ---------------- Export Consolidated Table ----------------
    print("4. Saving consolidated table...")
    output_tsv = f"{output_prefix}_summary.tsv"
    
    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Analysis Complete!\nTable saved as: {output_tsv}\nPlots are saved with prefix: {output_prefix}_")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Genome Annotation Quality Evaluation Tool")
    parser.add_argument("-f", "--fasta", required=True, help="Input Genome FASTA file")
    parser.add_argument("-g", "--gff", required=True, help="Input Annotation GFF3/GTF file")
    parser.add_argument("-o", "--output", default="eval_result", help="Output prefix (Default: eval_result)")

    args = parser.parse_args()
    run_analysis(args.fasta, args.gff, output_prefix=args.output)
