import gzip
import subprocess
import os
from collections import defaultdict
import sys
import pybedtools


from collections import defaultdict

def find_overlaps_with_cds(tx1_data, tx2_data, overlap_threshold=0.5):
    """
    基于转录本CDS总范围的快速重叠检测
    
    参数:
        tx1_data: 第一个GTF的解析结果 {tx_id: {'chrom': str, 'strand': str, 'exons': [...]}}
        tx2_data: 第二个GTF的解析结果（格式同上）
        overlap_threshold: 重叠比例阈值(default: 0.5)
        
    返回:
        [(tx1_id, tx2_id), ...] 重叠的转录本对列表
    """
    def extract_transcript_info(data_dict):
        """提取转录本的关键信息"""
        result = {}
        for tx_id, data in data_dict.items():
            # 关键修改: 优先使用 exon 计算重叠范围，如果没有 exon 则用 cds
            intervals = data['exon'] if data['exon'] else data['cds']
            if not intervals:
                continue
            min_start = min(s for s, e in intervals)
            max_end = max(e for s, e in intervals)
            total_length = sum(e - s for s, e in intervals)
            result[tx_id] = (
                data['chrom'],
                min_start,
                max_end,
                data['strand'],
                total_length
            )
        return result


    
    # 预处理转录本信息
    tx1_info = extract_transcript_info(tx1_data)
    tx2_info = extract_transcript_info(tx2_data)
    
    # 建立染色体索引
    tx1_chrom_index = defaultdict(list)
    tx2_chrom_index = defaultdict(list)
    
    for tx_id in tx1_info:
        tx1_chrom_index[tx1_info[tx_id][0]].append(tx_id)
    for tx_id in tx2_info:
        tx2_chrom_index[tx2_info[tx_id][0]].append(tx_id)
    
    # 结果存储
    results = []
    
    # 按染色体处理
    for chrom in tx1_chrom_index:
        if chrom not in tx2_chrom_index:
            continue
            
        # 处理当前染色体的转录本对
        for tx1_id in tx1_chrom_index[chrom]:
            chrom1, s1, e1, strand1, len1 = tx1_info[tx1_id]
            for tx2_id in tx2_chrom_index[chrom]:
                chrom2, s2, e2, strand2, len2 = tx2_info[tx2_id]
                
                # 只在同一条链上比较
                if strand1 != strand2:
                    continue
                
                # 计算重叠区间
                overlap_start = max(s1, s2)
                overlap_end = min(e1, e2)
                
                if overlap_end <= overlap_start:
                    continue
                
                # 计算重叠比例
                overlap_len = overlap_end - overlap_start
                overlap_ratio = max(
                    overlap_len / len1,
                    overlap_len / len2
                )
                
                if overlap_ratio >= overlap_threshold:
                    results.append((tx1_id, tx2_id))
    
    return results

def parse_exon_coords(gtf_file):
    """
    解析GTF，同时记录 exon 和 cds 的坐标
    """
    # 关键修改 1: 初始化时同时包含 exon 和 cds 的空列表
    transcripts = defaultdict(lambda: {'exon': [], 'cds': [], 'chrom': None, 'strand': None})
    
    try:
        opener = gzip.open if gtf_file.endswith('.gz') else open
        with opener(gtf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                feature_type = fields[2].lower()
                # 关键修改 2: 同时放行 exon 和 cds
                if feature_type not in ['exon', 'cds']:
                    continue
                
                chrom, strand = fields[0], fields[6]
                start, end = int(fields[3]), int(fields[4])
                
                attrs = {}
                for attr in fields[8].replace('";', '"; ').split(';'):
                    attr = attr.strip()          
                    if '=' in attr:
                        key, val = attr.split('=', 1)
                        attrs[key.strip()] = val.strip('"')
                    elif ' ' in attr:
                        key, val = attr.split(' ', 1)
                        attrs[key.strip()] = val.strip('"')
                
                tx_id = attrs.get('transcript_id', attrs.get('Parent'))
                if not tx_id:
                    continue
                
                if not transcripts[tx_id]['chrom']:
                    transcripts[tx_id]['chrom'] = chrom
                    transcripts[tx_id]['strand'] = strand
                
                # 关键修改 3: 根据特征类型分别存入对应的列表
                transcripts[tx_id][feature_type].append((start, end))
        
        return {
            tx: {
                'chrom': data['chrom'],
                'strand': data['strand'],
                'exon': sorted(data['exon']),
                'cds': sorted(data['cds'])
            }
            for tx, data in transcripts.items()
        }
    except Exception as e:
        print(f"Error parsing {gtf_file}: {str(e)}")
        raise



def calculate_aed(reference_exons, predicted_exons):
    # 完全相同检查（提前返回0）
    if reference_exons == predicted_exons:
        return 0.0

    # 1. 计算外显结构相似度
    def calculate_exon_similarity(ref, pred):
        # 对齐算法改进：考虑部分重叠
        match_count = 0
        #total_ref_len = sum(e[1]-e[0] for e in ref)
        #total_pre_len = sum(e[1]-e[0] for e in pred)
        total_ref_len=len(ref)
        total_pre_len=len(pred)

        i = j = 0
        while i < len(ref) and j < len(pred):
            ref_start, ref_end = ref[i]
            pred_start, pred_end = pred[j]
            
            #overlap = min(ref_end, pred_end) - max(ref_start, pred_start)
            #if overlap > 0:
            #    # 使用重叠比例和外显长度比例的几何平均
            #    overlap_ratio = overlap / (ref_end - ref_start)
            #    length_ratio = min(ref_end-ref_start, pred_end-pred_start) / max(ref_end-ref_start, pred_end-pred_start)
            #    total_similarity += (overlap_ratio * length_ratio) ** 0.5 * (ref_end - ref_start)
            
            if ref_start == pred_start:
                match_count += 1
            
            if ref_end == pred_end:
                match_count += 1
            # 移动指针
            if ref_end < pred_end:
                i += 1
            else:
                j += 1
                
        try:
            sn=match_count / (total_ref_len*2)
            sp=match_count / (total_pre_len*2)
            splice_similarity= (sn+sp)/2
            aed = 1.0 - splice_similarity
        except:
            return 1.0
        return  aed

    # 2. 剪接位点相似度（全部位点，包括5'/3'）
    def get_all_splice_sites(exons):
        sites = set()
        if len(exons) > 1:
            sites.update({exon[1] for exon in exons})  # 所有供体位点
            sites.update({exon[0] for exon in exons})   # 所有受体位点
        return sites
    
    
    ref_sites = get_all_splice_sites(reference_exons)
    pred_sites = get_all_splice_sites(predicted_exons)

    try:
        sn=len(ref_sites & pred_sites) / len(ref_sites)
        sp=len(ref_sites & pred_sites) / len(pred_sites)
        splice_similarity=(sn+sp)/2
        # 3. 综合计算（外显70% + 剪接30%）
        aed = 1.0 - splice_similarity
    except ZeroDivisionError:
        #如果其中一个参考是单外显子的，就比较长度。
        aed=calculate_exon_similarity(reference_exons, predicted_exons)

    #return aed
    return max(0.0, min(1.0, aed))

def compute_aed_for_overlapping_pairs(gtf1, gtf2, output_file="AED_results.txt"):
    """主函数：计算两个GTF间的AED""" 
    try:
        
        # 1. 解析外显子
        print("Parsing exon coordinates...")
        tx1_data = parse_exon_coords(gtf1)
        tx2_data = parse_exon_coords(gtf2)

        print("Finding overlapping feature...")
        overlaps = find_overlaps_with_cds(tx1_data, tx2_data)

        # 3. 计算AED
        print(f"Computing AED for {len(overlaps)} pairs...")
        results = defaultdict(list)
        for tx1, tx2 in overlaps:
            if tx1 in tx1_data and tx2 in tx2_data:
                aed_exon  = calculate_aed(tx1_data[tx1]["exon"], tx2_data[tx2]["exon"])
                aed_cds   =  calculate_aed(tx1_data[tx1]["cds"], tx2_data[tx2]["cds"])
                if aed_exon == None or aed_cds == None:
                    print(tx1)
                    print(tx2)
                best_aed = min(aed_exon, aed_cds)
                results[tx1].append((tx2, round(best_aed, 4)))

        # 5. 处理结果（修正部分）
        # 先转换results为(tx1,tx2)->aed的映射
        min_aed_results = {
            (tx1, min(pairs, key=lambda x: x[1])[0]): min(pairs, key=lambda x: x[1])[1]
            for tx1, pairs in results.items()
        }

        # 6. 写结果
        with open(output_file, 'w') as f:
            f.write("Transcript1\tTranscript2\tAED\tGeneID\n")
            for (tx1, tx2), aed in sorted(min_aed_results.items()):
                f.write(f"{tx1}\t{tx2}\t{aed}\n")

        print(f"Results saved to {output_file}")
        return min_aed_results  # 返回最终分组结果
        
    finally:
        # 清理临时文件
        pass
        #for f in [bed1, bed2]:
        #    if os.path.exists(f):
        #        os.remove(f)


if __name__ == "__main__":
    if len(sys.argv) <= 3:
        print("Usage: python script.py <gtf1> <gtf2>")
        sys.exit(1)

    pred_gtf = sys.argv[1]
    ref_gtf = sys.argv[2]
    output = sys.argv[3]

    if not os.path.exists(pred_gtf):
        print(f"File not found: {pred_gtf}")
        sys.exit(1)

    if not os.path.exists(ref_gtf):
        print(f"File not found: {ref_gtf}")
        sys.exit(1)

    compute_aed_for_overlapping_pairs(pred_gtf, ref_gtf,output)
