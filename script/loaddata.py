import sys
from collections import defaultdict

def load_junction_data(file_path):
    """
    简易版内含子数据加载函数
    输入: 文件路径
    输出: {gene_id: {'cds_count': int, 'cds': [cdsdict,...],"seq":[seq,...]}}
    """
    gene_data = {}
    with open(file_path) as f:
        for lineno, line in enumerate(f, 1):
            parts = line.strip().split('\t')
            if len(parts) != 5:
                if parts and parts[0]:
                    print(f"[loaddata] {file_path}:{lineno}: expected 5 tab-separated "
                          f"fields, skipped: {line[:60]!r}", file=sys.stderr)
                continue
            gene_id, cds_count, cds_length, cds_str, cds_seq = parts
            junctions = []

            for cds in cds_str.split('|'):
                phase, length ,donor, acceptor = cds.strip('()').split(',')
                junctions.append([int(phase) if phase != '.' else 0,int(length), donor, acceptor])

            if gene_id in gene_data:
                print(f"[loaddata] {file_path}:{lineno}: duplicate gene id "
                      f"{gene_id}, previous record overwritten", file=sys.stderr)
            gene_data[gene_id] = {
                'cds_count': int(cds_count),
                'cds_len':int(cds_length),
                'cds': junctions,
                'seq': cds_seq.strip().split('|')
            }
    return gene_data

# 读取内含子的数据库
def load_busco_junction_data(file_path):
    """
    简易版内含子数据加载函数
    输入: 文件路径
    输出: {busco, gene_id,  {'cds_count': int, 'cds': [cdsdict,...]}}
    """
    gene_data =defaultdict(dict)
    with open(file_path) as f:
        for lineno, line in enumerate(f, 1):
            parts = line.strip().split('\t')
            if len(parts) != 6:
                if parts and parts[0]:
                    print(f"[loaddata] {file_path}:{lineno}: expected 6 tab-separated "
                          f"fields, skipped: {line[:60]!r}", file=sys.stderr)
                continue
            busco_id, gene_id, cds_count, cds_length, cds_str,cds_seq = parts
            junctions = []

            for cds in cds_str.split('|'):
                phase, length ,donor, acceptor = cds.strip('()').split(',')
                junctions.append([int(phase) if phase != '.' else 0,int(length),  donor, acceptor])
            if gene_id in gene_data[busco_id]:
                print(f"[loaddata] {file_path}:{lineno}: duplicate gene id "
                      f"{gene_id} under {busco_id}, previous record overwritten", file=sys.stderr)
            gene_data[busco_id][gene_id] = {
                'cds_count': int(cds_count),
                'cds_len':int(cds_length),
                'cds': junctions ,
                'seq': cds_seq.strip().split('|')
                }
    return gene_data
