def convert_bed(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.strip() == "":
                continue  # 跳过空行
            fields = line.strip().split('\t')
            
            # (1) 转换 strand (第4列)
            strand_map = {"1": "+", "2": "-", "0": "."}
            strand = strand_map.get(fields[3], ".")  # 如果不是1/2/0，默认返回"."
            
            # (2) 计算 score = 第7列 + 第8列
            try:
                score = int(fields[6]) + int(fields[7])
            except (ValueError, IndexError):
                score = 0  # 如果转换失败，设为0
            
            # (3) 拼接新的BED行
            new_bed_line = "\t".join([
                fields[0],  # chrom
                fields[1],  # start
                fields[2],  # end
                strand,     # +/-
                str(score), # score
                "."         # 空列（通常BED6的第6列可忽略）
            ]) + "\n"
            
            outfile.write(new_bed_line)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python convert_bed.py input.bed output.bed")
        sys.exit(1)
    
    convert_bed(sys.argv[1], sys.argv[2])
    print(f"Conversion complete! Saved to {sys.argv[2]}")
