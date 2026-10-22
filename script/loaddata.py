from collections import defaultdict

def load_junction_data(file_path):
    """
    简易版内含子数据加载函数
    输入: 文件路径
    输出: {gene_id: {'cds_count': int, 'cds': [cdsdict,...],"seq":[seq,...]}}
    """
    gene_data = {}
    with open(file_path) as f:
        for line in f:
            gene_id, cds_count, cds_length, cds_str, cds_seq = line.strip().split('\t')
            junctions = []
            
            for cds in cds_str.split('|'):
                phase, length ,donor, acceptor = cds.strip('()').split(',')
                junctions.append([int(phase) if phase != '.' else 0,int(length), donor, acceptor])

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
        for line in f:
            busco_id, gene_id, cds_count, cds_length, cds_str,cds_seq = line.strip().split('\t')
            junctions = []
            
            for cds in cds_str.split('|'):
                phase, length ,donor, acceptor = cds.strip('()').split(',')
                junctions.append([int(phase) if phase != '.' else 0,int(length),  donor, acceptor])
            gene_data[busco_id][gene_id] = {
                'cds_count': int(cds_count),
                'cds_len':int(cds_length),
                'cds': junctions ,
                'seq': cds_seq.strip().split('|')
                }
    return gene_data