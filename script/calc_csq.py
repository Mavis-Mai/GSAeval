import pandas as pd
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
# 检查是否导入成功，如果没有则显示帮助信息
try:
    from loaddata import load_busco_junction_data, load_junction_data
    from csq_calculation import calc_csq_,calc_csq
except ImportError as e:
    print(f"Error: The necessary modules were not found. \n\
    Please ensure that it is in the parent directory of the script's directory.")
    print(f"Error text: {e}")
    sys.exit(1)

def busco_table(busco_txt:"full_table.tsv"):
    busco_info=pd.read_csv(busco_txt,sep="\t",comment="#",names=["busco_id","status","gene_id","score","length"])
    genelist=busco_info[["gene_id","busco_id"]]
    genelist=genelist.fillna("0")
    return genelist

def main(args_list=None):
    import argparse
    parser = argparse.ArgumentParser(description='Quick Splice Conservation Score (qSCS) calculation',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-a', '--annot_a', required=True, 
                help='For the library of the species.')
    parser.add_argument('-A', '--annot_A', 
                help='For the library of the second species, please correspond to the gene pair file.')
    parser.add_argument( '-b','--busco_lib', 
                help='The library of high conservation BUSCO genes.')
    parser.add_argument('-p', '--pairs',  
                help='The gene pair file has two gene ids on each line, separated by tabs')
    parser.add_argument('-l','--busco_table', default="full_table.tsv", 
                help='Detailed information of the busco corresponding table for proteins')
    parser.add_argument('-w', '--weight_type', default='plant',
                help='Based on the species input, select plant or animal weight coefficients to calculate CSQ.[plant or animal]]')
    parser.add_argument('-o', '--output', default='output', 
                       help='output file')
    #args = parser.parse_args()
    args = parser.parse_args(args_list) if args_list else parser.parse_args()
    # 加载注释数据
    print("load library...")
    a_data = load_junction_data(args.annot_a)
    if args.annot_A:
        b_data = load_junction_data(args.annot_b)
    elif args.busco_lib:
        busco_data=load_busco_junction_data(args.busco_lib)
    else:
        b_data=a_data

    if args.pairs:
        # 处理基因对
        print("qSCS calculating ...")
        results = []
        dirfile=os.path.dirname( os.path.abspath(args.pairs))
        ouput_withoutIntron=open(f"{dirfile}/genepairs_withoutIntron.txt",'w+')
        with open(args.pairs) as f:
            for line in f:
                gene1, gene2 = line.strip().split('\t')[:2]
                if gene1 not in a_data :
                    ouput_withoutIntron.write(f'{gene1}\n')
                    continue
                    #print(f"警告: 基因 {gene1} 或 {gene2} 无注释信息，跳过")
                elif gene2 not in b_data :
                    ouput_withoutIntron.write(f'{gene2}\n')
                    continue
                try:
                    cdss_a,cdss_b=a_data[gene1],b_data[gene2]
                    score= calc_csq(cdss_a,cdss_b,args.weight_type)
                    #print(f'{gene1}_{gene2} finish qSCS')

                    results.append({
                        'GeneA': gene1,
                        'GeneB': gene2,
                        **score
                    })
                except TypeError:
                    continue


    elif args.busco_table and busco_data:
        #处理基因列表 最好为geneid\tbusco
        results = []
        print(os.path.abspath(args.busco_table))
        dirfile=os.path.dirname(os.path.abspath(args.busco_table))
        #print(dirfile)
        ouput_withoutIntron=open(f"{dirfile}/gene_withoutIntron.txt",'w+')
        genelist=busco_table(args.busco_table)
       
        print(genelist)
        for i in range(len(genelist)):
            gene1,busco=genelist.iloc[i]
            if gene1!="0":
                if gene1 not in a_data:
                    ouput_withoutIntron.write(f'{gene1}\n')
                    continue
                try:
                    busco_data_tmp=busco_data[busco]
                    intorn_a=a_data[gene1]
                    #print(intorn_a)
                    for gene_b in busco_data_tmp:            
                        score= calc_qscs(intorn_a,busco_data_tmp[gene_b])
                        #print(score)
                        results.append({
                                'BUSCO': busco,
                                'GeneA': gene1,
                                'GeneB': gene_b,
                                **score
                        })
                    
                except TypeError:
                    continue

        #这块可以考虑不放进去，因为不属于精确性评估的范畴
            #else:
                #if busco in busco_data:
                    #print(f'{busco} missing in the dataset')
                    #results.append({'BUSCO': busco,'GeneA':"Missing",'GeneB': "NA",
                                            #'cdsA number':0,
                                            #'cdsB number':0,
                                            #'CDS score': 0,
                                            #'Match score':0,
                                            #'Phase score':0,
                                            #'Length similarity': 0,
                                            #'Intron site score':0,
                                            #'Total length score':0,
                                            #'AED':1,
                                            #'qSCS': 0})
    else:
        print("Either the gene list or the gene pair needs to be Provided.")
    
    # 保存结果
    ouputfile=f'{args.output}.csq.tsv'
    pd.DataFrame(results).to_csv(ouputfile, sep='\t', index=False)
    print(f"Results save as: {ouputfile}")

if __name__ == '__main__':
    main()

