import sys
import re
import pandas as pd
from collections import defaultdict
from extarctGeneTranscirpt import extract_mappings
def extract_introns(gxf_file,parent_pattern,output,type="exon"):
    """从GFF3文件中提取内含子信息，并保存到intron.txt"""
    exons = defaultdict(list)
    strands = {}
    chrdict={}
    dict_map=extract_mappings(gxf_file,"gene_transcript.map")

    # 第一步：解析GFF3文件，按GeneID存储外显子和链信息
    with open(gxf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue  # 跳过注释行
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != type:
                continue
            attrs = fields[8]
            try:
                gene_id = parent_pattern.search(attrs).group(1)
                
                start, end = int(fields[3]), int(fields[4])
                exons[gene_id].append((start, end))
                strands[gene_id] = fields[6]  # 链信息（+/-）
                chrdict[gene_id]=fields[0]
            except AttributeError:
                print(line)
                pass

    # 第二步：计算内含子（按GeneID排序的外显子之间的区域）
    intron_data = []
    for gene_id in sorted(exons.keys()):
        exon_list = sorted(exons[gene_id], key=lambda x: x[0])  # 按start排序
        strand = strands[gene_id]
        chr=chrdict[gene_id]
        for i in range(len(exon_list) - 1):
            current_end = exon_list[i][1]
            next_start = exon_list[i + 1][0]
            
            # 内含子区间：[current_end + 1, next_start - 1]
            intron_start = current_end 
            intron_end = next_start - 1
            intron_length = intron_end - intron_start + 1
            
            intron_data.append({
                "GeneID": gene_id,
                "Gene_id":dict_map[gene_id]["GeneID"],
                "Chr": chr,
                "start": intron_start,
                "end": intron_end,
                "length": intron_length,
                "strand": strand
            })

    # 第三步：输出结果到文件
    df = pd.DataFrame(intron_data)
    outputfile=f'{output}_sjfromgff.bed'
    if len(df) > 0:
        df.to_csv(outputfile, sep="\t", index=False)
        avg_length = df["length"].mean()
        #print(f"Intron信息已保存到 {outputfile}，平均长度: {avg_length:.2f}")
    else:
        print("Warning: No introns detected!")

def main(args_list=None):
    import argparse
    parser = argparse.ArgumentParser(description='extarct splice junction from gxf file.')
    parser.add_argument('-g', '--gxf', required=True, help='gff3/gtf file')
    parser.add_argument('-o', '--output', default='output', 
                       help='output prefix')
    parser.add_argument('-r','--recordType', default='exon',
                       help='extract junction from feature type')
    parser.add_argument('--item', default="",
                       help='extract exon name from item. with --gxfform')                
    #args = parser.parse_args()
    #args = parser.parse_args(args_list) if args_list else parser.parse_args()
    if args_list:
        args, _ = parser.parse_known_args(args_list)  # 返回 (Namespace, remaining_args)
    else:
        args = parser.parse_args()  # 返回 Namespace

    if args.item!="":
        print(f"extract the feature: {args.item}")
        parent_pattern = re.compile(rf'{args.item}[=\s]*["\']?([A-Za-z0-9_.\-:|]+)["\']?[;\s]?') # 兼容更复杂的ID格式 
    else:
        print(f"extract the feature: exon")
        parent_pattern = re.compile(r'transcript_id[=\s]*["\']?([A-Za-z0-9_.\-:|]+)["\']?[;\s]?')

    
    extract_introns(args.gxf,parent_pattern,args.output,args.recordType)


if __name__ == "__main__":
    main()

