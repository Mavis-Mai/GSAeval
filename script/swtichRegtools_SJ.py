import sys

def convert_bed(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#') or line.strip() == '':
                continue  # 跳过注释行或空行
            
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue  # 跳过不合规的行
            
            chrom = fields[0]
            donor_ = int(fields[1])  # 供体位点（内含子开始）
            acceptor_ = int(fields[2])  # 受体位点（内含子结束）
            name = fields[3]
            strand = fields[5]
            cov=fields[4]
            donor_bias=int(fields[10].split(",")[0])
            acceptor_bias=int(fields[10].split(",")[1])
            
            # 计算内含子区间（保留原始数据）
            intron_start = donor_+donor_bias
            intron_end = acceptor_-acceptor_bias
            
            # 输出新的BED行（6列）
            new_line = "\t".join([
                chrom,
                str(intron_start),
                str(intron_end),
                strand,
                str(cov),
                name + ":intron",
            ]) + "\n"
            outfile.write(new_line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_introns.py input.bed output.bed")
        sys.exit(1)
    
    input_bed = sys.argv[1]
    output_bed = sys.argv[2]
    convert_bed(input_bed, output_bed)
    print(f"Introns successfully extracted and saved to {output_bed}")
