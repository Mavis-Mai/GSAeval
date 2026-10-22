from collections import defaultdict
import pandas as pd 

def parse_cds_coords(gtf_file):
    """解析GTF/GFF3，返回{transcript_id: [(start, end), ...]}"""
    cds_dict = defaultdict(list)
    exon_dict = defaultdict(list)
    
    print(gtf_file)
    try:
        opener = gzip.open if gtf_file.endswith('.gz') else open
        with opener(gtf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                feature=fields[2].lower()
                
                if len(fields) < 9 or feature not in ["cds",'exon']:
                    #print(feature)
                    continue
                
                # 解析属性
                attrs = {}
                for attr in fields[8].replace('";', '"; ').split(';'):
                    attr = attr.strip()
                    if '=' in attr:
                        key, val = attr.split('=', 1)
                        attrs[key.strip()] = val.strip('"')
                    elif ' ' in attr:
                        key, val = attr.split(' ', 1)
                        attrs[key.strip()] = val.strip('"')
                
                tx_id = attrs.get('transcript_id') or attrs.get('Parent')
                if not tx_id:
                    continue

                start, end = int(fields[3]), int(fields[4])
                if feature=="cds":          
                    cds_dict[tx_id].append(end-start)
                else:
                    exon_dict[tx_id].append(end-start)
        
        # 计算每个基因的cds比例
        #print(transcript_dict)
        #print(cds_dict)
        cds_ratio_dict={}
        for gene in exon_dict.keys():
            cds_ratio=sum(cds_dict[gene])/sum(exon_dict[gene])*100
            cds_ratio_dict[gene]=max(0.0,min(cds_ratio,100.0))
        
        
        df = pd.DataFrame.from_dict(
            cds_ratio_dict, 
            orient='index', 
            columns=["CDS_ratio"]
        ).reset_index().rename(columns={"index": "GeneID"})
        print(df)
        #df = pd.DataFrame.from_dict(cds_ratio_dict, orient='index', columns=["GeneID",'CDS_ratio'])
        return df
    except Exception as e:
        print(f"Error parsing {gtf_file}: {str(e)}")
        raise


if __name__ == "__main__":
    import sys
    gtf_file=sys.argv[1]
    output=sys.argv[2]
    result=parse_cds_coords(gtf_file)
    result.to_csv(f"{output}.cds.tbl",sep='\t')