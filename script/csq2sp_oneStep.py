import sys
from pathlib import Path
import importlib.util
import argparse
import os 
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import shutil
import glob
import tempfile

sys.path.insert(0, str(Path(__file__).parent))
try:
    from build_junction_lib_v5 import *
    from loaddata import *
    from csq_calculation import *
    from csq_oneStep import *
    from generatereport import *
except ImportError as e:
    print(f"Error: Required module overlap_resolver was not found. Ensure it is in the parent directory of the script directory.")
    print(f"Error details: {e}")
    sys.exit(1)

def single_diamond_blast(seq_q, seq_r, threads, k=2):
    """Return each query's best forward hit, preserving complete gene IDs."""
    diamond_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "diamond")
    print("Commencing one-way comparison (Query -> Reference)...")
    with tempfile.TemporaryDirectory(prefix="diamond_single_") as tmp_dir:
        database = os.path.join(tmp_dir, "ref.dmnd")
        output = os.path.join(tmp_dir, "q2r.tsv")
        subprocess.run(
            [diamond_path, "makedb", "--in", os.path.abspath(seq_r),
             "--db", database, "--quiet"], check=True)
        subprocess.run(
            [diamond_path, "blastp", "--db", database,
             "--query", os.path.abspath(seq_q), "--threads", str(threads),
             "--outfmt", "6", "qseqid", "sseqid", "bitscore",
             "-o", output, "-k", str(k), "--evalue", "1e-5", "--quiet"],
            check=True)
        best_hits = {}
        with open(output) as handle:
            for line in handle:
                if not line.strip():
                    continue
                query, reference, score = line.rstrip().split("\t")
                score = float(score)
                if query not in best_hits or score > best_hits[query][1]:
                    best_hits[query] = (reference, score)
    return pd.DataFrame(
        [(query, reference, score) for query, (reference, score) in best_hits.items()],
        columns=["Gene_Q", "Gene_R", "Bitscore_Q2R"])


def reciprocal_diamond_blast(seq_q, seq_r, threads, k=1,  program="blastp", out_tsv="reciprocal_pairs.tsv"):
    """
    使用 Diamond 进行互惠最佳比对 (Reciprocal Best Hits)
    
    参数:
    seq_q : str, Query Fasta 文件的路径
    seq_r : str, Reference Fasta 文件的路径
    k     : int, Diamond 报告的最大比对数 (默认 1，RBH通常只需找 top 1)
    threads : int, 运行线程数
    program : str, 运行模式，一般蛋白质比对用 "blastp"，核酸比蛋白用 "blastx"
    out_tsv : str, 最终生成的互惠基因对列表保存路径
    
    返回:
    list of tuples: 包含 (Query_gene, Reference_gene, Q2R_Bitscore, R2Q_Bitscore) 的列表
    """
    
    print("building  Diamond database...")
    db_r = "temp_ref_db"
    db_q = "temp_query_db"
    
    # 1. 建库
    diamond_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),"diamond")
    print(f"running diamond")
    subprocess.run(
        [diamond_path, "makedb", "--in", os.path.abspath(seq_r), 
            "--db",f"{db_r}.dmnd", "--quiet"],
            check=True)
    subprocess.run(
        [diamond_path, "makedb", "--in", os.path.abspath(seq_q), 
            "--db",f"{db_q}.dmnd", "--quiet"],
            check=True)

    print("Commencing forward comparison (Query -> Reference)...")
    out_q2r = "temp_q2r.tsv"
    outfmt = "6 qseqid sseqid pident evalue bitscore"

    subprocess.run(
        [diamond_path, "blastp", "--db",f"{db_r}.dmnd", "--outfmt", "6","--quiet",
            "--query",  os.path.abspath(seq_q), "--threads", f"{threads}", "-o",out_q2r ,"-k", f"{k}", "--evalue","1e-5"],
            check=True)

    
    print("Commencing reverse comparison (Reference -> Query)...")
    out_r2q = "temp_r2q.tsv"
    subprocess.run(
        [diamond_path, "blastp", "--db",f"{db_q}.dmnd", "--outfmt", "6","--quiet",
            "--query",  os.path.abspath(seq_r), "--threads", f"{threads}" ,"-o", out_r2q ,"-k", f"{k}", "--evalue","1e-5"],
            check=True)
    
    print("Reciprocal gene pairs are being extracted....")
    
    # 内部解析函数：获取每个 Query 的 Best Hit (基于 Bitscore)
    def get_best_hits(tsv_file):
        best_hits = {} # q_id -> (s_id, bitscore)
        with open(tsv_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                q_id, s_id = parts[0], parts[1]
                bitscore = float(parts[11])
                
                # 如果这个基因还没被记录，或者遇到了更高的分数，则更新
                if q_id not in best_hits or bitscore > best_hits[q_id][1]:
                    best_hits[q_id] = (s_id, bitscore)
        return best_hits

    # 获取正反向的 Best Hit
    q2r_best = get_best_hits(out_q2r)
    r2q_best = get_best_hits(out_r2q)
    
    # 寻找互惠对 (Reciprocal)
    reciprocal_pairs = []
    for q_gene, (r_gene, q2r_score) in q2r_best.items():
        if r_gene in r2q_best:
            best_q_for_r, r2q_score = r2q_best[r_gene]
            # 如果 R 的最佳匹配也是 Q，则它们互为最佳比对 (RBH)
            if best_q_for_r == q_gene:
                reciprocal_pairs.append((q_gene, r_gene, q2r_score, r2q_score))
                
    # 写入结果文件
    with open(out_tsv, 'w') as f:
        f.write("Gene_Q\tGene_R\tBitscore_Q2R\tBitscore_R2Q\n")
        for pair in reciprocal_pairs:
            f.write(f"{pair[0]}\t{pair[1]}\t{pair[2]}\t{pair[3]}\n")
            
    print(f"Done! A total {len(reciprocal_pairs)} reciprocal gene pairs were identified.")
    print(f"The results have been saved to: {out_tsv}")
    
    # 清理临时文件
    for tmp_file in [f"{db_r}.dmnd", f"{db_q}.dmnd", out_q2r, out_r2q]:
        if os.path.exists(tmp_file):
            os.remove(tmp_file)
    
    df = pd.DataFrame(reciprocal_pairs, columns=["Gene_Q", "Gene_R", "Bitscore_Q2R", "Bitscore_R2Q"])

    return df

def csq_calculation2db(a_data,b_data,genelist,output,type):
    results = []
    ouput_withoutIntron=open(f"gene_withoutIntron.txt",'w+')
    #print(genelist)
    for i in range(len(genelist)):
        #print(genelist.iloc[i][:2])
        gene1,gene_2=genelist.iloc[i][:2]
        if gene1!="0":
            if gene1 not in a_data:
                ouput_withoutIntron.write(f'{gene1}\n')
                continue
            if gene_2 not in b_data:
                ouput_withoutIntron.write(f'{gene_2}\n')
                continue
            score= calc_csq(a_data[gene1],b_data[gene_2],type)
            results.append({
                'GeneA': gene1,
                'GeneB': gene_2,
                **score
                })

    ouputfile=f'{output}.csq.tsv'
    pd.DataFrame(results).to_csv(ouputfile, sep='\t', index=False)
    print(f"Results have save in {ouputfile}.")  
    return ouputfile

def read_table(file):
    busco_info=pd.read_csv(file,sep="\t",comment="#",names=["gene1","gene2"],usecols=[0, 1])
    genelist=busco_info[["gene1","gene2"]]
    genelist=genelist.fillna("0")
    return genelist
    
def drawdenstiy(file):
    if not os.path.isfile(file) or os.path.getsize(file) == 0:
        print(f"No CSQ records in {file}; skip density plot.")
        return
    dataset1=pd.read_csv(file,sep="\t",header=0)
    if dataset1.empty or dataset1['CSQ'].dropna().empty:
        print(f"No CSQ records in {file}; skip density plot.")
        return
    rs_95=dataset1['CSQ'].quantile(0.25)
    plt.figure(figsize=(6, 4)) 
    ax=sns.kdeplot(
        data=dataset1,
        x='CSQ',
        #fill=True,
        color='gray',
        alpha=0.5,)
    if not ax.lines:
        print(f"CSQ values in {file} are too few or degenerate to draw a density curve; skip.")
        plt.close()
        return
    line = ax.lines[0]
    x = line.get_xdata()
    y = line.get_ydata()
    max_density = y.max()

    max_density = y.max()
    plt.axvline(rs_95, color='red', linestyle=':', alpha=0.3)
    plt.text(x=rs_95-0.1,y=max_density*0.75,s=f'Q1:{rs_95:.2}')
    plt.savefig("CSQ-density.png")
    plt.savefig("CSQ-density.svg")


def main(args_list=None):
    parser = argparse.ArgumentParser(description='Process GTF and genome files to extract splice site information.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--query_fasta', required=True, help='Input query genome FASTA file (required)')
    parser.add_argument('-g', '--query_gtf', required=True, help='Input query GTF annotation file (required)')
    parser.add_argument('-F', '--ref_fasta', required=True, help='Input reference genome FASTA file (required)')
    parser.add_argument('-G', '--ref_gtf', required=True, help='Input reference GTF annotation file (required)')
    parser.add_argument('-w', '--weight_type', required=True, help='Based on the species input, select plant or animal weight coefficients to calculate CSQ.[plant or animal]]  (required)')
    parser.add_argument('-l', '--genepairs', help='Provide high-confidence gene pair information; \n\
                        If omitted, DIAMOND maps gene pairs using the selected alignment mode')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--reciprocal', dest='reciprocal', action='store_true',
                       help='Enable reciprocal best-hit DIAMOND alignment (default)')
    group.add_argument('--no-reciprocal', dest='reciprocal', action='store_false',
                       help='Use one-way DIAMOND alignment')
    parser.set_defaults(reciprocal=True)
    parser.add_argument('-k', '--max-target-seqs', dest='max_target_seqs', type=int, default=2, help='hit the number of Diamonds')
    parser.add_argument('-d', '--type', default="CDS", help='mode:CDS or exon.')

    parser.add_argument('-t', '--threads', default=8, help='cpu number')
    parser.add_argument('--report', action="store_true",  default=False, help='Whether to generate report(.md)' )
    parser.add_argument('-o', '--output', default='output', help='output prefix')
    parser.add_argument('--fig_format',default="svg", help='The selection of graphic format. such as:svg, png, pdf')
    parser.add_argument('--debug', action='store_true', help='Enable debug mode')

    args = parser.parse_args(args_list) if args_list else parser.parse_args()

    #读取数据
    print(args)
    db_query_file=build_junciton_lib(args.query_fasta,args.query_gtf,args.type,"query",args)
    db_ref_file=build_junciton_lib(args.ref_fasta,args.ref_gtf,args.type,"ref",args)
    for lib in (db_query_file, db_ref_file):
        if os.path.getsize(lib) == 0:
            sys.exit(f"No spliced (multi-exon) genes were found in {lib}; "
                     f"CSQ analysis is based on splice junctions and requires genes with introns.")
    db_query=load_junction_data(db_query_file)
    db_ref=load_junction_data(db_ref_file)
    print("The data has been loaded.")
    #形成对应列表
    
    if args.genepairs == None:
        queryfile=getProt(db_query,"query",False)
        reffile=getProt(db_ref,"ref",False)
        if args.reciprocal:
            genepairs=reciprocal_diamond_blast(queryfile,reffile,args.threads,k=args.max_target_seqs)
        else:
            genepairs=single_diamond_blast(queryfile,reffile,args.threads,k=args.max_target_seqs)
        pairsnum=len(genepairs)
        print(genepairs)
        qCSStsvfile=csq_calculation2db(db_query,db_ref,genepairs,args.output,args.weight_type)
        drawdenstiy(qCSStsvfile)
    else:
        genelist=read_table(args.genepairs)
        pairsnum=len(genelist)
        print(genelist)
        qCSStsvfile=csq_calculation2db(db_query,db_ref,genelist,args.output,args.weight_type)
        drawdenstiy(qCSStsvfile)

    #summarize(qCSStsvfile,threshold_table,args.output)
    dirpath = Path(__file__).parent.parent
    current_dir=os.getcwd()
    if args.report:
        templatePath=os.path.join(dirpath,"Template/Template.md")
        target_template=os.path.join(os.getcwd(),"Template.md")
        safe_copy_dir(templatePath, target_template)
        pattern = os.path.join(current_dir, f"*density.png")
        copypict(pattern,"CSQ2sp.png","Template.md")

        pattern = os.path.join(os.getcwd(), f"{args.output}.csq.tsv")
        CSQ2sptsv=glob.glob(pattern)[0]
        CSQ2spdf = pd.read_csv(CSQ2sptsv,header=0,sep="\t")
        CSQ2spdf = CSQ2spdf[["GeneA","GeneB","cdsA_number","cdsB_number","CSQ"]]
        reportfile=os.path.join(os.getcwd(),"Template.md/CSQ2sp-template.md")
        print(CSQ2spdf.head())
        generate_report_tables(CSQ2spdf,reportfile)


        replaceKeywordNumber(pairsnum,"Number of homologous gene pairs: ", reportfile)
        
        #除了CSQ-template.md其他报告模板都删了
        dir_path = Path("Template.md")
        exclude_file = "CSQ2sp-template.md"
        for filename in os.listdir(dir_path):
            filepath = os.path.join(dir_path, filename)
            if os.path.isfile(filepath) and filename != exclude_file:
                os.remove(filepath)

        #除了CSQ2sp.png其他图片都删了  
        dir_path = Path("Template.md/_v_images")
        exclude_file = "CSQ2sp.png"
        for filename in os.listdir(dir_path):
            filepath = os.path.join(dir_path, filename)
            if os.path.isfile(filepath) and filename != exclude_file:
                os.remove(filepath)       

        print("The report has been generated.")

if __name__ == '__main__':
    main()
