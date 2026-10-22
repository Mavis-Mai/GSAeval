#!/usr/bin/env python3
import importlib.util
from pathlib import Path
import sys
import argparse
import os
import warnings
import shutil
import glob
import subprocess
import shutil

sys.path.insert(0, str(Path(__file__).parent))
try:
    from bed_comparisonv6 import BedAnalyzer,validate_inputs,timeit
    from merge_bam import BAMMerger,get_bam_file
    from swtichRegtools_SJ import convert_bed
    from build_junction_lib_v5 import gff3_to_gtf
    from merge_gene_quality_info import merge_tables_by_geneid
    from AEDv5 import *
    from count_CDS_ratio import parse_cds_coords
    from generatereport import *
except ImportError as e:
    print(f"Error text: {e}")
    sys.exit(1)

def parse_args(args_list=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-g","--gxf", required=True, help="gene anotation (gff/gtf).")
    parser.add_argument("-b","--bam_floder_file", nargs='+',required=True, help="A directory with RNA-seq files(bam).")
    
    #parser.add_argument("-s",'--RNAseqSJ', required=True, help='Splice site information in the transcriptome (bed) ')
    parser.add_argument("-o", "--output", default="output", help="output prefiex")
    parser.add_argument("--force", default=False,action='store_true', help="Force regenerate the result (overwrite the existing file)")
    parser.add_argument("-f","--fasta",  help="gene assembly (fasta).")
    parser.add_argument("-m", "--min_overlap", type=float, default=0.5, help="Minimum overlap ratio (0-1)")
    parser.add_argument('-c', '--min_support', type=int, default=5,
                                help='The minimum number threshold of supported read segments (default: %(default)s)')
    parser.add_argument('-t', '--threads', type=int, default=8,
                                help=f'Number of parallel threads (default:8 CPUs)')
    parser.add_argument('--item', default="", help='extract exon information from item. with --gxfform')
    if args_list is None:
        args_list = sys.argv[1:]  # 默认使用命令行参数
    
    args, unknown = parser.parse_known_args(args_list)
    
    return args,unknown


class one_step:
    def __init__(self):
        self.version = "1.0.0"
        self.script_dir = os.path.dirname(os.path.abspath(__file__))

    def run_tool(self,module_name,args_list):
        script_path = os.path.join(self.script_dir, f"{module_name}.py")
        if not os.path.exists(script_path):
            raise FileNotFoundError(f"scirpt {module_name}.py cannot found in {self.script_dir} ")   
        # 加载模块
        spec = importlib.util.spec_from_file_location(module_name, script_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
                
        # 确保有main函数
        if not hasattr(module, 'main'):
            raise AttributeError(f"Tool module '{module_name}' has no 'main' function")
                
        # 执行工具的主函数并传递参数
        print(f"\nExecuting tool: {module_name}")
        import subprocess
        #print(args_list)
        module.main(args_list) 

    def merge_bam(self,bam_floder,gff,cpu,output_dir,force=False):
        #print(bam_floder)
        if os.path.isfile(bam_floder[0]):
            print(f"input file {bam_floder[0]}")
            return bam_floder[0]

        bam_file_paths=get_bam_file(bam_floder,1)
        print(bam_file_paths)
        if len(bam_file_paths)>1:
            merge=BAMMerger()
            os.makedirs(output_dir, exist_ok=True)
            print(f"{cpu} cpus are used")
            #检查是否已经合并了bam文件。
            if os.path.exists(os.path.join(output_dir,"merged_sorted.bam")) and not force:
                return  os.path.join(output_dir,"merged_sorted.bam")
            #先考虑进来的bam文件已经实现了sort
            #merge.parallel_sort_and_index_bamfiles(bam_file_paths, output_dir, cpu)
            #bam_file_paths = glob.glob(output_dir + "/*bam")           
            #print("sorting and indexing complete")

            chr_len=merge.get_chrom_lengths(bam_file_paths[0])
            gene_region_list,unsort_all_region_list=merge.split_gff_to_bed(gff,chr_len,5000,output_dir)
            all_region_list=sorted(unsort_all_region_list,key=lambda x: (x[0], x[1]))
            store_set=set()
            for tuple in gene_region_list:
                store_set.add(tuple[0])
            for chr_name in chr_len.keys():
                if chr_name not in store_set:
                    print(f"warning, {chr_name} cannot be found in the gff file,this chromosome will not be split region")
            print("split region succssfully")
            print(all_region_list[:10])

            big_dict=merge.store_bam_region_coverage_parallel(bam_file_paths,all_region_list)

            merge.merge_regions_to_bam(bam_file_paths,big_dict,os.path.join(output_dir,"merged_region.txt"),
                        all_region_list,output_dir, 100000)
            print("Start sorting and indexing the merged bam files")

            merge.sort_and_index_merged_bam(os.path.join(output_dir, "merged.bam"), 
            os.path.join(output_dir,"merged_sorted.bam"),
            os.path.join(output_dir,"merged_sorted.bam.bai"))
            
            print("Complete sorting and indexing ")


            bamfile=os.path.join(output_dir,"merged_sorted.bam")
        elif len(bam_file_paths)==1:
            bamfile=bam_file_paths[0]
        elif len(bam_file_paths)==0:
            sys.exit(1)
        return bamfile

    def run(self,args_list):
        print(args_list)
        
        args,unknown=parse_args(args_list)
        outputgtf=f"{args.gxf}.gtf"
        gff3_to_gtf(args.gxf,outputgtf)
        args.gxf=outputgtf

        try:
            gxf_index = args_list.index("--gxf")
        except ValueError:
            gxf_index = args_list.index("-g")
        args_list[gxf_index + 1] = outputgtf

        os.makedirs("Quality_info", exist_ok=True) 
        bamfile=self.merge_bam(args.bam_floder_file,args.gxf,args.threads,"Quality_info",args.force)
        print(f"bamfile: {bamfile}")
        
        if os.path.exists("Quality_info/merged_sorted_SJ_regtools.tsv") and not args.force:
            pass
        elif os.path.exists("merged_sorted_SJ_regtools.tsv") and not args.force:
            pass
        else:
            regtools_cmd=f'regtools junctions extract  {bamfile} -o Quality_info/merged_sorted_SJ_regtools.tsv -a 4  -m 20  -s XS'
            subprocess.run(regtools_cmd, shell=True, check=True)

        RNAseqSJ=f"{args.output}.SJfromBam.bed"
        try:
            convert_bed("Quality_info/merged_sorted_SJ_regtools.tsv", RNAseqSJ)
        except:
            convert_bed("merged_sorted_SJ_regtools.tsv", RNAseqSJ)
        args_list.append("--bam_file")
        args_list.append(bamfile)
        args_list.append("--RNAseqSJ")
        args_list.append(RNAseqSJ)

        genefromgff=f'{args.output}_genefromgff.bed'
        genefrombam=f'{args.output}_label.genefromRNA.bed'
        #genefrombam=f'{args.output}_r1.geneloci.bed'
            

        self.run_tool("extract_gene_from_GFFv2",args_list)
        print(f"generate file: {genefromgff}")
        self.run_tool("extract_gene_from_Bamv2",args_list)
        print(f"generate file: {genefrombam}")


        ##gene overlap
        validate_inputs(genefrombam, genefromgff)
        analyzer = BedAnalyzer(genefrombam, genefromgff, args.output, args.min_overlap)
        analyzer.run() 
        os.makedirs("Gene_Region", exist_ok=True) 
        for file_path in glob.glob("Gene*"):
            if os.path.isfile(file_path): 
                shutil.move(file_path, os.path.join("Gene_Region", os.path.basename(file_path)))
            
        ##SJ extarct
        self.run_tool("extract_SJ_from_GFF",args_list)
        sjfromgff=f'{args.output}_sjfromgff.bed'
        print(f"generate file: {sjfromgff}")
        
        args_list.append("--genebed")
        args_list.append(genefromgff)
        args_list.append('--geneSJ')
        args_list.append(sjfromgff)
        if os.path.exists(f"Splice_junctions/SJ_{args.output}_junction_match_stat.tsv") and not args.force:
            pass
        else:
            self.run_tool("checkGeneSJv5",args_list)
            os.makedirs("Splice_junctions", exist_ok=True) 
            for file_path in glob.glob("SJ*"):
                if os.path.isfile(file_path):  
                    shutil.move(file_path, os.path.join("Splice_junctions", os.path.basename(file_path)))

        

        for file_path in glob.glob(f"{args.output}*"):
            if os.path.isfile(file_path):
                shutil.move(file_path, os.path.join("Quality_info", os.path.basename(file_path)))

        ##计算基因的cds比例？
        df_cds_ratio=parse_cds_coords(args.gxf)
        df_cds_ratio.to_csv(f"Quality_info/{args.output}_cds_ratio.txt",sep="\t",index=False)
        ##计算基因间的aed
        compute_aed_for_overlapping_pairs(outputgtf, 
                    f"Quality_info/{args.output}_orf.gff", 
                    f"Quality_info/{args.output}.CDS_distance.txt",
                   )
        mergedf=merge_tables_by_geneid(f"Gene_Region/Gene_{args.output}_detailed_overlap.tsv",
                                       f"Splice_junctions/SJ_{args.output}_detail_comparing.tsv",
                                       f"Quality_info/{args.output}.CDS_distance.txt",
                                       f"Quality_info/{args.output}_cds_ratio.txt",
                                       args.gxf,
                                       f"Quality_info/{args.output}")
        print(mergedf.head())

        mergedf.to_csv("Total_gene_quality.tsv",sep="\t",index=False)
        dirpath = Path(__file__).parent.parent
        current_dir=os.getcwd()
        if "--report" in args_list:
            templatePath=os.path.join(dirpath,"Template/Template.md")
            target_template=os.path.join(os.getcwd(),"Template.md")
            safe_copy_dir(templatePath, target_template)
            #RNAseq
            pattern = os.path.join(current_dir,  "Gene_Region", "*_venn.png")
            copypict(pattern,"RNAseqvenn.png","Template.md")

            pattern = os.path.join(current_dir,  "Splice_junctions", "*_pieplot.png")
            copypict(pattern,"SJlevel.png","Template.md")

            pattern = os.path.join(current_dir,  "Quality_info", "*_riskscore_distribution.png")
            copypict(pattern,"Riskscore.png","Template.md")

            ### 读取RNA-seq generegion的部分信息
            template=os.path.join(os.getcwd(),"Template.md/RNAseq-template.md")

            pattern = os.path.join(current_dir, "Gene_Region", "*_stats.tsv")
            regionstat=glob.glob(pattern)[0]
            genedict={}
            numpattern = re.compile(r'(\d+)\((\d+\.\d+)%?\)')
            with open(regionstat,'r+')as f:
                for line in f:
                    if "\t" in line:
                        line=line.strip().split("\t")
                        #print(line)
                        text,num=line
                        genedict[text]=num
            #nornanum=re.find(genedict["without RNA-seq support"],numpattern)[0]
            grouplist=re.findall(numpattern,genedict["without RNA-seq support"]) 
            #print(genedict["without RNA-seq support"])
            overlapregions=int(grouplist[0][0])
            overlapregionspic=float(grouplist[0][1])
            text =f'{overlapregions} ({overlapregionspic}%)'
            replace_markdown_table(text,"genes without RNA-seq support",template)

            grouplist=re.findall(numpattern,genedict["Total overlap region:"])
            overlapregions=int(grouplist[0][0])
            overlapregionspic=float(grouplist[0][1])*0.01
            abtranscriptregions=int(genedict["total transcript region "])
            text=f'{abtranscriptregions-overlapregions} ({100*(1-overlapregionspic)}%)'
            replace_markdown_table(text,"RNA-seq support regions without gene models",template)

        ### 读取SJ的信息
            pattern = os.path.join(current_dir, "Splice_junctions", "*_gene_stat.tsv")
            regionstat=glob.glob(pattern)[0]
            sjleveldf=pd.read_csv(regionstat,sep="\t",header=0,index_col=0)
            levelb_num = int(sjleveldf.loc["LEVEL B", "Status"])
            levelb_pic = sjleveldf.loc["LEVEL B", "Status"] / sjleveldf["Status"].sum()
            text = f'{levelb_num} ({levelb_pic * 100:.2f}%)'
            replace_markdown_table(text,"genes with conflicting splicing junctions (SJ)",template)

        ###修改RNA模式输出的内容方便读取
            pattern = os.path.join(current_dir, "Quality_info", "*_quantiles_statistics.tsv")
            quantiles_statistics=glob.glob(pattern)[0]
            quantiledf=pd.read_csv(quantiles_statistics,sep="\t",header=0,index_col=0)
            
            CDS_distanceQ95 = float(quantiledf.loc["CDS_distance", "Quantile_0.95"])
            Risk_scoreQ95 = float(quantiledf.loc["Risk_score", "Quantile_0.95"])
            SJ_scoreQ95 = float(quantiledf.loc["SJ_score (1 - sj_score)", "Quantile_0.95"])
            Overlap_scoreQ95 = float(quantiledf.loc["Overlap_score", "Quantile_0.95"])

            replace_markdown_table(f'{CDS_distanceQ95:.2f}',"CDS distance (Q0.95)",template)
            replace_markdown_table(f'{Risk_scoreQ95:.2f}',"Risk score (Q0.95)",template)
            replace_markdown_table(f'{Overlap_scoreQ95:.2f}',"Overlap score (Q0.05)",template)
            replace_markdown_table(f'{SJ_scoreQ95:.2f}',"SJ score (Q0.05)",template)

            #除了CSQ-template.md其他报告模板都删了
            dir_path = Path("Template.md")
            exclude_file = "RNAseq-template.md"
            for filename in os.listdir(dir_path):
                filepath = os.path.join(dir_path, filename)
                if os.path.isfile(filepath) and filename != exclude_file:
                    os.remove(filepath)

            #除了RNAseqvenn.png其他图片都删了  
            dir_path = Path("Template.md/_v_images")
            exclude_file = ["RNAseqvenn.png","SJlevel.png","Riskscore.png"]
            for filename in os.listdir(dir_path):
                filepath = os.path.join(dir_path, filename)
                if os.path.isfile(filepath) and filename not in  exclude_file:
                    os.remove(filepath)       

            print("The report has been generated.")

def main(args_list=None):
    parser = argparse.ArgumentParser(description="Compare the gene regions, junction, TES/TSS between RNA-seq files and annotation", 
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-g","--gxf", required=True, help="gene anotation (gff/gtf).")
    parser.add_argument("-b","--bam_floder_file", nargs='+',required=True, help="A directory with RNA-seq files(bam).")
    parser.add_argument("-f","--fasta", required=True, help="gene assembly (fasta).")   
    parser.add_argument("-o", "--output", default="output", help="output prefiex")
    parser.add_argument("--force", default=False,action='store_true', help="Force regenerate the result (overwrite the existing file)")
    parser.add_argument("-m", "--min_overlap", type=float, default=0.5, help="Minimum overlap ratio (0-1)")
    parser.add_argument('-c', '--min_support', type=int, default=5,
                                help='The minimum number threshold of supported read segments (default: %(default)s)')
    parser.add_argument('-t', '--threads', type=int, default=8,
                                help=f'Number of parallel threads (default:8 CPUs)')
    parser.add_argument('-r','--recordTyep', default='exon', help='extract junction from feature type')
    parser.add_argument('--item', default="", help='extract exon information from item. with --gxfform')
    parser.add_argument('--idItem', default="", help='extract gene name from item. with --gxfform')
    parser.add_argument('--report', action="store_true",  default=False, help='Whether to generate report(.md)' )
    args = parser.parse_args(args_list) if args_list else parser.parse_args()
    analyzer = one_step()
    analyzer.run(args_list)

if __name__ == "__main__":
    main()
