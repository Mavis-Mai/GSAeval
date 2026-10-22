import sys
import gzip
import argparse
from collections import defaultdict
import numpy as np

def print_debug(msg, args):
    """调试信息打印"""
    if args.debug:
        print(f"[DEBUG] {msg}", file=sys.stderr)

def parse_attributes(attr_str):
    """解析GTF文件的属性字段"""
    attr = {}
    for item in attr_str.split(';'):
        item = item.strip()
        if not item:
            continue
        if ' ' in item:  # 处理 key "value" 格式
            key, value = item.split(' ', 1)
            value = value.strip('"')
        elif '=' in item:  # 处理 key=value 格式
            key, value = item.split('=', 1)
            value = value.strip('"')
        else:
            continue
        attr[key] = [value] if key not in attr else attr[key] + [value]
    return attr

def gff3_to_gtf(gff3_file, gtf_file):
    import subprocess
    subprocess.run(["gffread", gff3_file, "-T", "-o", gtf_file])

def parse_gtf(gtf_file, args={}):
    """解析GTF文件，返回transcript、外显子和CDS信息"""
    transcripts = {}  # transcript_id -> (chr, start, end, strand)
    exons = defaultdict(list)  # transcript_id -> [(start, end)]
    cdses = defaultdict(list)  # transcript_id -> [(start, end, phase)]
    
    try:
        opener = gzip.open if gtf_file.endswith('.gz') else open
        with opener(gtf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                feature_type = fields[2]
                attributes = parse_attributes(fields[8])
                
                if 'transcript_id' not in attributes:
                    if 'gene_id' in attributes:
                        print_debug(f"Line has gene_id but no transcript_id: {line[:60]}...", args)
                    continue
                
                transcript_id = attributes['transcript_id'][0]
                
                try:
                    if feature_type == 'transcript':
                        transcripts[transcript_id] = (
                            fields[0],  # chr
                            int(fields[3]),  # start
                            int(fields[4]),  # end
                            fields[6]   # strand
                        )
                    elif feature_type == 'exon':
                        exons[transcript_id].append((
                            int(fields[3]),
                            int(fields[4])
                        ))
                    elif feature_type == 'CDS':
                        phase = int(fields[7]) if fields[7] != '.' else 0
                        cdses[transcript_id].append((
                            int(fields[3]),
                            int(fields[4]),
                            phase
                        ))
                except ValueError as e:
                    print_debug(f"Error parsing line: {str(e)}", args)
                    continue
        
        print_debug(f"Parsed {len(transcripts)} transcripts", args)
        return transcripts, exons, cdses
        
    except Exception as e:
        print(f"Error reading GTF file: {str(e)}", file=sys.stderr)
        #sys.exit(1)

def merge_exons_cds(exons, cdses):
    """合并外显子和CDS信息"""
    merged = defaultdict(list)
    if exons !={}:
        for tr_id in exons:
            exon_list = sorted(exons[tr_id])
            cds_list = sorted(cdses.get(tr_id, []))

            # 为每个外显子查找对应的CDS信息
            for exon in exon_list:
                exon_phase = "."
                for cds in cds_list:
                    if cds[0] == exon[0] and cds[1] == exon[1]:  # 完全重叠
                        exon_phase = cds[2]
                        break
                    elif exon[0] < cds[0] and cds[1] == exon[1]: #3utr
                        exon_phase = cds[2]
                        break
                    elif exon[0] == cds[0] and   exon[1] > cds[1]: #5utr
                        exon_phase = cds[2]
                        break
                    else: 
                        continue
                merged[tr_id].append((exon[0], exon[1], exon_phase))
    else:
        for tr_id in cdses:
            cds_list = sorted(cdses.get(tr_id, []))
            for cds in cds_list:
                merged[tr_id].append((cds[0], cds[1], cds[2]))
    
    return merged

def calculate_entry(exon, chrom_seq, isFirst, isEnd, strand):
    """计算单个外显子条目"""
    exon_start, exon_end, phase = exon
    exon_len = exon_end - exon_start + 1
    
    donor, acceptor = None, None
    if isFirst:
        if strand == '+':
            left =  None
            right = chrom_seq[exon_end:exon_end+2].upper() if exon_end+2 <= len(chrom_seq) else None
        else:
            left= None
            right = reverse_complement(chrom_seq[exon_start-3:exon_start-1].upper()) if exon_start-2 >=0 else None
    elif isEnd:
        if strand == '+':
            left = chrom_seq[exon_start-3:exon_start-1].upper() if exon_start-2 <= len(chrom_seq) else None
            right = None
        else:
            left = reverse_complement(chrom_seq[exon_end:exon_end+2].upper()) if exon_end+2 <= len(chrom_seq) else None
            right = None
    else:
        if strand == '+':
            left = chrom_seq[exon_start-3:exon_start-1].upper() 
            right = chrom_seq[exon_end:exon_end+2].upper() 
        else:
            left = reverse_complement(chrom_seq[exon_end:exon_end+2].upper()) 
            right = reverse_complement(chrom_seq[exon_start-3:exon_start-1].upper()) 
    return f"({phase if phase is not None else 'None'},{exon_len},{left if left else '.'},{right if right else '.'})"

def reverse_complement(seq):
    """获取反向互补序列"""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join([comp.get(base, base) for base in reversed(seq)])

def load_genome(genome_file):
    """加载基因组序列"""
    genome = {}
    current_chr = None
    seq = []
    
    opener = gzip.open if genome_file.endswith('.gz') else open
    
    with opener(genome_file, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_chr and seq:
                    genome[current_chr] = ''.join(seq)
                current_chr = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if current_chr and seq:
            genome[current_chr] = ''.join(seq)
    
    standardized_genome = {}
    for chr_name in genome:
        standardized_genome[chr_name] = genome[chr_name]
    
    return standardized_genome

def process_transcript_exons(exons, strand):
    """处理外显子顺序"""
    sorted_exons = sorted(exons, key=lambda x: (x[0], x[1]))
    if strand == '-':
        sorted_exons = sorted_exons[::-1]
    return sorted_exons

def process_transcript(transcripts,merged_exons,genome):
    # 处理每个transcript
    print("Processing transcripts...", file=sys.stderr)
    results = {}
    missing_chrom = set()
    for tr_id in transcripts:
        chrom, start, end, strand = transcripts[tr_id]
        
        if chrom not in genome:
            missing_chrom.add(chrom)
            continue
        
        tr_exons = merged_exons.get(tr_id, [])
        processed_exons = process_transcript_exons(tr_exons, strand)
        exon_count = len(processed_exons)
        
        entries = []
        cds_list=[]
        cds_length=[]
        if exon_count>1:
            for i in range(exon_count):
                exon = processed_exons[i]
                #print(exon)
                cds_list.append(genome[chrom][exon[0]-1:exon[1]])
                cds_length.append(exon[1]-exon[0]+1)
                try:
                    if i == 0:
                        entry = calculate_entry(exon, genome[chrom],True,False, strand)
                    elif i ==exon_count-1:
                        entry = calculate_entry(exon, genome[chrom],False,True, strand)
                    else:
                        entry = calculate_entry(exon, genome[chrom], False,False, strand)
                    entries.append(entry)
                except Exception as e:
                    print(f"Skipping exon due to error: {str(e)}", file=sys.stderr)
                    continue
            if strand=="-":
                cds_list=[reverse_complement(seq) for seq in cds_list]

            tr_cds=np.array(cds_length).sum()
            results[tr_id] = f"{tr_id}\t{exon_count}\t{tr_cds}\t{'|'.join(entries)}\t{'|'.join(cds_list)}"
    
    # 输出缺失染色体警告
    if missing_chrom:
        print(f"Warning: {len(missing_chrom)} chromosomes not found in genome: {', '.join(sorted(missing_chrom)[:5])}{'...' if len(missing_chrom)>5 else ''}", file=sys.stderr)
    return results

def write_tsv(output_file, results):
    """将结果写入TSV文件"""
    with open(output_file, 'w') as f:
      #  f.write("TranscriptID\tExonCount\tExonDetails\n")
        for tr_id in sorted(results.keys()):
            f.write(results[tr_id] + "\n")
    



def main(args_list=None):
    parser = argparse.ArgumentParser(description='Process GTF and genome files to extract splice site information.')
    parser.add_argument('-f', '--fasta', required=True, help='Input genome FASTA file (required)')
    parser.add_argument('-g', '--gtf', required=True, help='Input GTF annotation file (required)')
    parser.add_argument('-o', '--output', default='output_lib.tsv', help='Output TSV file (default: output_lib.tsv)')
    parser.add_argument('-d', '--type', default="CDS", help='mode:CDS or exon.')
    parser.add_argument('--debug', action='store_true', help='Enable debug mode')

    args = parser.parse_args(args_list) if args_list else parser.parse_args()
    print(f"Processing with parameters:", file=sys.stderr)
    print(f"  Genome FASTA: {args.fasta}", file=sys.stderr)
    print(f"  GTF file: {args.gtf}", file=sys.stderr)
    print(f"  Output file: {args.output}", file=sys.stderr)
    
    # 解析GTF文件
    print("Parsing GTF file...", file=sys.stderr)
    outputgtf=f"{args.gtf}.gtf"
    gff3_to_gtf(args.gtf,outputgtf)
    transcripts, exons, cdses = parse_gtf(outputgtf, args)
    #print(cdses)    
    
    if not transcripts:
        print("Error: No transcripts found in GTF file!", file=sys.stderr)
        sys.exit(1)
    
    # 合并外显子和CDS
    if args.type=="exon":
        merged_exons = merge_exons_cds(exons, cdses)
    else:
        merged_exons = merge_exons_cds({}, cdses)
    
    # 加载基因组
    print("Loading genome sequence...", file=sys.stderr)
    genome = load_genome(args.fasta)
    print(f"Loaded {len(genome)} chromosomes", file=sys.stderr)
    
    # 写输出文件
    results=process_transcript(transcripts,merged_exons,genome)
    print(f"Writing results to {args.output}...", file=sys.stderr)
    write_tsv(args.output, results)
    
    print(f"Done! Processed {len(results)} transcripts.", file=sys.stderr)
    print(f"Sample output (first 3):", file=sys.stderr)
    for tr_id in sorted(results.keys())[:3]:
        print(results[tr_id], file=sys.stderr)

if __name__ == '__main__':
    main()
