import sys
from pathlib import Path
import importlib.util
import argparse
import os 
import subprocess
import pandas as pd
import shutil
import glob
from collections import defaultdict
import itertools
import math
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, str(Path(__file__).parent))
try:
    from build_junction_lib_v5 import *
    from csq_oneStep import *
    from loaddata import *
    from csq_calculation import *
    from csq2sp_oneStep import csq_calculation2db,read_table
    from CSQ_statv4_forBUSCO import *
    from runThresholdv2 import *
except ImportError as e:
    print(f"Error: Required module overlap_resolver was not found. Ensure it is in the parent directory of the script directory.")
    print(f"Error details: {e}")
    sys.exit(1)

## 1. 文件目录结构：
#|____busco
#| |____Brap_full_table.tsv
#| |____Mper_full_table.tsv
#| |____TAIR12_full_table.tsv
#|____genome
#| |____sp3
#| | |____GCA_978657495.1_TAIR12_genomic.gff
#| | |____GCA_978657495.1_TAIR12_genomic_clean.fna
#| |____sp2
#| | |____GWHERGL00000000.genome.fasta
#| | |____GWHERGL00000000.gff
#| |____sp1
#| | |____Mperfoliatum_583_v2.0.fa
#| | |____Mperfoliatum_583_v2.1.gene.gff3

## 2 建库
##从输入文件中读取对应的基因组和注释
def build_lib_Sps(genome_dir,args):
    species_data = {}
    liblist=[]
    # 遍历所有子目录（跳过文件，只处理目录）
    for species_dir in [d for d in glob.glob(os.path.join(genome_dir, "*")) if os.path.isdir(d)]:
        species_name = os.path.basename(species_dir)  # 子目录名作为物种名
        print(f'loading directory: {species_name}')
        
        genome_files =(glob.glob(os.path.join(species_dir, "*.f*a")) + glob.glob(os.path.join(species_dir, "*.f*a.gz"))) # 匹配 .fa/.fna/.fasta
        if not genome_files:
            print(f"⚠️ FATAL: Genome file not found in {species_name}.")
            continue
        genome_path = genome_files[0]

        annot_files = (glob.glob(os.path.join(species_dir, "*.gff*"))+glob.glob(os.path.join(species_dir, "*.gtf")))
        if not annot_files:
            print(f"⚠️ FATAL: Annotation file not found in {species_name}.")
            continue
        annot_path=annot_files[0]
        print(genome_path)
        print(annot_path)
        output_lib=os.path.join(species_dir,species_name)
        output_lib_file=output_lib+"_lib.tsv"
        if not os.path.exists(output_lib_file):
            outputfile=build_junciton_lib(genome_path,annot_path,"CDS",output_lib,args)
        
        liblist.append(output_lib_file)
    return liblist

## 2 提取保守busco的id和信息
def find_conserved_busco_ids(busco_tsv_dir, output_file):
    file_list=[d for d in glob.glob(os.path.join(busco_tsv_dir, "*.tsv"))]
    # 1. 确定保守数量阈值
    num_files = len(file_list)
    
    if num_files < 5:
        threshold = num_files
    elif 5 <= num_files < 10:
        threshold = num_files - 1
    else:
        threshold = int(num_files * 0.9)  # 取90%并向下取整
    
    print(f"Analysis of {num_files} files, with the conservative threshold set to: {threshold}")
    
    # 2. 读取所有文件中的Complete/Duplicated BUSCO ID
    busco_counts = defaultdict(int)
    total_buscos = set()
    
    for file_path in file_list:          
        with open(file_path, 'r') as f:
            file_buscos = set()
            for line in f:
                #if line.startswith(('Complete', 'Duplicated')):
                parts = line.strip().split('\t')
                if len(parts) >= 3 and parts[1].lower() in ['complete', 'duplicated']:
                    busco_id = parts[0]
                    #print(busco_id)
                    file_buscos.add(busco_id)
            
            # 更新全局统计
            for busco_id in file_buscos:
                busco_counts[busco_id] += 1
            total_buscos.update(file_buscos)
    #print(busco_counts)
    # 3. 筛选满足阈值的BUSCO ID
    conserved_buscos = [
        busco_id for busco_id, count in busco_counts.items() 
        if count >= threshold
    ]
    
    # 4.输出到文件
    output=os.path.join(busco_tsv_dir,f'{output_file}.id')
    if output:
        with open(output, 'w') as f:
            f.write("\n".join(conserved_buscos))
        print(f"Conservative BUSCO IDs have been saved to:  {output_file}.id")
    
    print(f"Found  {len(conserved_buscos)} conserved BUSCO IDs")

    output=os.path.join(busco_tsv_dir,f"{output_file}.txt")
    with open(output, 'w') as o:
        for buscoid in conserved_buscos:
            for file_path in file_list:          
                with open(file_path, 'r') as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if  buscoid== parts[0]:
                            o.write(line)

    return conserved_buscos ,f'{output_file}.id', f'{output_file}.txt'

## 3 基因对信息计算
def interCombination(input_file, output_prefix):
    # 1. 读取基因ID并去重（直接使用Pandas Series）
    print(input_file)
    output_filename = f"{output_prefix}.txt"
    gene_ids = pd.read_csv(input_file, header=None, names=["id"])["id"].drop_duplicates()
    gene_ids = list(gene_ids)

    # 2. 生成组合迭代器（避免内存中的列表过大）
    pairs = itertools.combinations(gene_ids, 2)

        # 3. 批量写入（减少I/O次数）
        
    with open(output_filename, "w") as f:
            # 一次性写入所有组合（生成器表达式在写入时逐个求值，节约内存）
        f.writelines(f"{a}\t{b}\n" for a, b in pairs)

        # 计算数量
    num_genes = len(gene_ids)
    num_pairs = num_genes * (num_genes - 1) // 2

## 4 提取合适的保守基因
def merge_gene_data(*dicts):
    merged = {}
    for d in dicts:  # 遍历所有输入的字典
        for busco_id, genes in d.items():
            if busco_id not in merged:
                merged[busco_id] = {}
            merged[busco_id].update(genes)  # 合并 gene_id 层级（相同 gene_id 会被覆盖）
    return merged

##6 整理所有文件。
def create_combined_summary(buscodir,output='conservative_BUSCO_summary.tsv'):
    with open(output, 'w+') as outfile:
        for root, _, files in os.walk(buscodir):
            for file in files:
                if file.endswith('_overall_summary.tsv'):
                    filepath = os.path.join(root, file)
                    print(f'reading {filepath}')
                    with open(filepath, 'r') as infile:
                        for line in infile:
                            if "BUSCO_id" not in line:
                                outfile.write(line)

def filter_busco_ids(input='conservative_BUSCO_summary.tsv',output='filter_busco.id',
                    min_exon_num = 5,min_genepairs = 10,min_score = 0.9):
    filtered_ids = []
    
    with open(input, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            try:
                if (float(parts[2]) >= min_exon_num and 
                    float(parts[1]) > min_genepairs and 
                    float(parts[11]) > min_score):
                    # Extract ID from path (3rd component when split by '/')
                    busco_id = parts[0].split('/')[2]
                    filtered_ids.append(busco_id)
            except (IndexError, ValueError):
                continue  # Skip malformed lines
    
    with open(output, 'w') as f:
        for busco_id in filtered_ids:
            f.write(f"{busco_id}\n")


def process_threshold_files(buscodir,buscofile="filter_busco.id",output='filter_busco_threshold.txt'):
    # Prepare output file
    with open(output, 'w') as outfile:
        # Read filtered BUSCO IDs
        with open(buscofile, 'r') as f:
            for busco_id in f:
                busco_id = busco_id.strip()
                input_file = f'{buscodir}/{busco_id}/{busco_id}_threshold.txt'
                if os.path.exists(input_file):
                    with open(input_file, 'r') as infile:
                        # Skip header (1st line)
                        next(infile)  
                        for line in infile:
                            # Add BUSCO ID as first column
                            line=line.strip().split("\t")
                            info="\t".join(line[1:])
                            outfile.write(f"{busco_id}\t{info}\n")

def final_processing(output='output',input='filter_busco_threshold.txt', min_csq=0.8, min_length = 100):
    header = 'busco\tgenepairs_number\texon_number\tthreshold\tqscs\n'

    # Read and filter data
    filtered_lines = []
    with open(input, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            try:
                if (float(parts[3]) > min_csq and 
                    float(parts[1]) > min_length):
                    filtered_lines.append(line)
            except (IndexError, ValueError):
                continue
    
    # Write final output
    with open(output, 'w') as f:
        f.write(header)
        f.writelines(filtered_lines)
    
    print(f"Processing complete. Final file contains {len(filtered_lines)} records.")

def process_single_busco(busco, cbdf, merged_data, weight_type="plant"):
    """
    独立的工作函数：处理单个 busco 的完整逻辑
    """
    buscodir = os.path.join("Dataset", busco)
    os.makedirs(buscodir, exist_ok=True)
    
    tmpdf = cbdf.loc[cbdf["busco"] == busco]
    buscogene = tmpdf["geneid"]
    buscogene_file = os.path.join(buscodir, f'{busco}_gene.id')
    
    with open(buscogene_file, "w+") as output:
        for gid in buscogene:
            output.write(str(gid) + "\n")
            
    # 形成基因对文件：
    pair_file = os.path.join(buscodir, 'genepairs')
    interCombination(buscogene_file, pair_file)
    genelist = read_table(pair_file + ".txt")
    output_path = os.path.join(buscodir, f"{busco}")
    if os.path.exists(f"{output_path}.csq.tsv") and  os.path.exists(f"{output_path}_threshold.txt") and os.path.exists(f"{output_path}_gene_stat.tsv"):
        pass
    else:
        try:
            qCSStsvfile = csq_calculation2db(merged_data, merged_data, genelist, output_path, weight_type)
            
            # 计算阈值 (修复了原代码的 & 符号 Bug)
            genenum = len(tmpdf["geneid"])
            if genenum <= 5:
                top_n = genenum
            elif 5 < genenum <= 10:  # 修复了这里的逻辑
                top_n = genenum - 1
            else:
                top_n = round(genenum * 0.8)

            # 打印加上 [busco] 前缀，避免多进程打印混乱
            print(f'[{busco}] Achieve {top_n} gene as reference dataset.')
            print(f'[{busco}] {qCSStsvfile}')
            
            stats_df, stats_df_gene, summary_df, fig = process_gene_data(
                qCSStsvfile, output_path,
                target_col="CSQ",
                use_prefix=False,
                top_n=top_n,
                visualize=False,
                do_summary=True,
            )

            # 保存结果
            stats_df.to_csv(output_path + "_species_stat.tsv", index=False, sep="\t")
            stats_df_gene.to_csv(output_path + "_gene_stat.tsv", index=False, sep="\t")
            
            # 保存汇总文件
            if summary_df is not None and not summary_df.empty:
                summary_df.to_csv(output_path + "_overall_summary.tsv", index=False, sep="\t")

            df = pd.read_csv(qCSStsvfile, header=0, sep="\t")
            threshold, csq = plot_boxplot_and_get_threshold(df, busco, buscodir, column_name='CSQ')
            csq_str = [str(i) for i in csq]
            exon_m = round(np.array(df["cdsA_number"]).mean(), 3)
            csq_num = str(len(df["cdsB_number"]))

            outputfile = os.path.join(buscodir, f"{busco}_threshold.txt")
            with open(outputfile, 'w+') as f:
                f.write('busco\tgenepairs_number\texon_number\tthreshold\tcsq\n')
                f.write(f'{output_path}\t{csq_num}\t{str(exon_m)}\t{threshold}\t{",".join(csq_str)}\n')          
            print(f"[{busco}] --- Done ---")
            return busco  # 成功完成，返回 busco 名称
        except pd.errors.EmptyDataError:
            print(f'[{busco}] Genes in {busco} without splicing.')
            print(f"[{busco}] ------------------------------------------------------------")
            return None  # 出现异常，返回 None

def main(args_list=None):
    parser = argparse.ArgumentParser(description='Self-build conservative database.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--genome_dir', required=True, help='Input directory including genome(fa) and annotation(gxf) of different species')
    parser.add_argument('-l', '--busco_tsv_dir', required=True, help='Input directory including busco gene pairs list(tsv) of different species')
    parser.add_argument('-o', '--output', default='output', help='output prefix of database directory')
    parser.add_argument('-n', '--cpu', default=10, type=int, help='cpu number used in parallel jobs')
    parser.add_argument('--debug', action='store_true', help='Enable debug mode')
    parser.add_argument('-w', '--weight_type', default='plant', choices=['plant', 'animal'],
                help='Weight coefficients used when building the database')

    args = parser.parse_args(args_list) if args_list else parser.parse_args()
    #建库
    liblist=build_lib_Sps(args.genome_dir,args)
    conserved_buscos,cb_id,cb_file=find_conserved_busco_ids(args.busco_tsv_dir,"conserved_buscos")
    merged_data={}
    for lib in liblist:
        df=load_junction_data(lib)
        merged_data = merge_gene_data(merged_data, df)

    #创建文件：
    os.makedirs("Dataset",exist_ok=True)

    cbdf=pd.read_csv(os.path.join(args.busco_tsv_dir, cb_file),sep="\t",
                    names=["busco","status","geneid","score","length","OrthoDB","utr","desc."])
    conserved_buscos_ref=[]
    cores_to_use = args.cpu
    print(f"Starting parallel processing with {cores_to_use} cores...")
    with ProcessPoolExecutor(max_workers=cores_to_use) as executor:
        # 提交所有任务到进程池
        futures = [
            executor.submit(process_single_busco, busco, cbdf, merged_data, args.weight_type)
            for busco in conserved_buscos
        ]
        
        # 获取结果 (as_completed 会在任务完成时立刻返回，不用等待按顺序结束)
        for future in as_completed(futures):
            result = future.result()
            # 如果成功返回了 busco (不是 None)，则添加到列表中
            if result is not None:
                conserved_buscos_ref.append(result)
    print(f"Total conserved buscos refs saved: {len(conserved_buscos_ref)}")

    
    #合并阈值和参考集合
    min_length=0
    if len(liblist)<=5:
        min_length=math.comb(int(len(liblist)*0.9),2)
    else:
        min_length=math.comb(int(len(liblist)*0.8),2)
    
    create_combined_summary("Dataset",output='conservative_BUSCO_summary.tsv')
    filter_busco_ids(input='conservative_BUSCO_summary.tsv',output='filter_busco.id',min_genepairs = min_length)
    process_threshold_files("Dataset",buscofile="filter_busco.id",output='filter_busco_threshold.txt')
    final_processing(output=f'{args.output}_threshold.txt',input='filter_busco_threshold.txt', min_csq=0.8,min_length = min_length)

    ##读取所有
    output_file = f'{args.output}_lib.txt'
   # merged_data = pd.DataFrame(merged_data)
    threshold_df=pd.read_csv(f"{args.output}_threshold.txt",sep="\t",header=0)

    f=open(output_file,"w+")
    for busco in threshold_df["busco"]:
        gene_stat_file = Path(f"Dataset/{busco}/{busco}_gene_stat.tsv")
        try:
            df_stat = pd.read_csv(gene_stat_file, sep="\t", header=0)           
            top_genes = df_stat.loc[df_stat["Top_Flag"]=="Top"]["Gene"]
            #print(top_genes)
            for gene in top_genes:
                tmpinfo=merged_data[gene]

                formatted_data = [
                    f"({','.join(map(str, item))})"  # 将子列表转为字符串
                    for item in tmpinfo["cds"]
                ]
                cds_result = "|".join(formatted_data)
                seq_result = "|".join(tmpinfo["seq"])
                f.write(f'{busco}\t{gene}\t{tmpinfo["cds_count"]}\t{tmpinfo["cds_len"]}\t{cds_result}\t{seq_result}\n')
        except:
            print(f"⚠️ Unable to read or process {gene_stat_file}")
    f.close()
            

    ##删除文件
    files_to_remove = ["filter_busco_threshold.txt", "filter_busco.id","gene_withoutIntron.txt"]
    for file_path in files_to_remove:
        try:
            os.remove(file_path)  # 删除文件
            print(f"✅ Successfully deleted: {file_path}")
        except FileNotFoundError:
            print(f"⚠️ File does not exist: {file_path}")
        except Exception as e:
            print(f"❌ Failed to delete: {file_path} - Error: {e}")

if __name__ == "__main__":
    main()

