#!/usr/bin/env python3
import os
import subprocess
import sys
import argparse
import pandas as pd
from intervaltree import Interval, IntervalTree
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
# 检查是否导入成功，如果没有则显示帮助信息
try:
    import overlap_resolver
    from build_junction_lib_v5 import gff3_to_gtf
    from gtf_gff_orf_predict import *
except ImportError as e:
    print(f"Error text: {e}")
    sys.exit(1)


SCRIPT_NAME = "transcript_processor"
VERSION = "1.0.0"
DESCRIPTION = "Process BAM files to extract transcript information using StringTie and generate BED files"

def check_required_tools():
    required_tools = ['stringtie','gffread']
    missing_tools = []
    
    for tool in required_tools:
        try:
            subprocess.run([tool, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except FileNotFoundError:
            missing_tools.append(tool)
    
    if missing_tools:
        print("Error: The following necessary tools are not installed or are not in the PATH:")
        for tool in missing_tools:
            print(f"- {tool}")
        print("\n Please install these tools first before running this script.")
        sys.exit(1)


def print_help_on_error(func):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            error_msg = f"""
            ERROR in [{func.__name__}] 
            → Exception: {type(e).__name__} 
            → Message: {str(e)}
            → Arguments:
               - Positional args: {args}
               - Keyword args: {kwargs}
            """
            print(f'{error_msg}')
    return wrapper

@print_help_on_error
def run_stringtie(bam_file, gff3, output_gtf, threads, genome_fasta, output_prefix, isforce):
    """Run StringTie2 to process BAM file and generate GTF"""
    cmd = [
        "stringtie",
        bam_file,
        "-G", gff3,
        "-o", output_gtf,
        "-l", "STRG",
        "-a",f"{output_prefix}_cov.txt",
        "-m", "500",
        "-p", str(threads),
        "-c", "5",
        "-s", "10",
        "-j" ,"5",
        "-f" , "0.2"
    ]

    print(f"Running StringTie: {' '.join(cmd)}")
    if os.path.exists(output_gtf) and not isforce:
        print(f"{output_gtf} have already generated, skip stringtie process")
        pass
    elif os.path.exists(f'Quality_info/{output_gtf}') and not isforce:
        print(f"Quality_info/{output_gtf} have already generated, skip stringtie process")
        output_gtf=f"Quality_info/{output_gtf}"
        pass
    else:
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error: StringTie failed with exit code {e.returncode}", file=sys.stderr)
            sys.exit(1)

        
    return output_gtf



@print_help_on_error
def process_gtf_to_bed(gtf_file, output_bed):
    """Process GTF file to generate BED format and DataFrame
    
    Returns:
        pd.DataFrame: DataFrame containing transcript information
    """
    print(f"Processing GTF file: {gtf_file}")
    
    # Initialize data storage
    data = {
        'Chr': [],
        'Start': [],
        'End': [],
        'Strand': [],
        'cov': [],
    }
    
    with open(gtf_file) as f_in, open(output_bed, "w") as f_out:
        # Write BED header
        f_out.write("Chr\tStart\tEnd\tStrand\tgene_id\ttranscript_id\tCoverage\tFPKM\tTPM\n")
        
        for line in f_in: 
            if line.startswith("#") or not line.strip():
                continue            
            fields = line.strip().split("\t")
            attributes = fields[8].split(";")
            if fields[2] != "transcript":
                #print(line)
                continue
            #print("check")
            # Parse attributes
            attr_dict = {}
            #print(attributes)
            for attr in attributes:
                if "\"" in attr:
                    try:
                        key, value = attr.strip().split()
                        key = key.replace("\"", "").strip()
                        value = value.replace("\"", "").strip()
                        attr_dict[key] = value
                    except ValueError:
                        pass
                        #print(attributes)
            
            # Store data
            data['Chr'].append(fields[0])
            data['Start'].append(int(fields[3]))
            data['End'].append(int(fields[4]))
            data['Strand'].append(fields[6])
            data['cov'].append(float(attr_dict.get("cov", ".")))
            #data['full_attributes'].append(fields[8])
            
            # Write to BED
            output_fields = [
                fields[0],  # Chr
                fields[3],  # Start
                fields[4],  # End
                fields[6],  # Strand
                attr_dict.get("gene_id", "."),
                attr_dict.get("transcript_id", "."),
                attr_dict.get("cov", "."),
                attr_dict.get("FPKM", "."),
                attr_dict.get("TPM", ".")
            ]
            f_out.write("\t".join(output_fields) + "\n")
    
    # Create DataFrame
    df = pd.DataFrame(data)
    df["Length"]=df["End"]-df["Start"]
    return df

@print_help_on_error
def process_overlap(df):
    # Split the dataframe based on column 4 values
    #过滤完全重叠的位置：
    #print(df)
    result=df.groupby(["Chr", "Start", "End"], group_keys=False).apply(lambda x: x.loc[(x["cov"]).idxmax()]) 
    #print(result)
    # Combine the results and sort by original index
    return result.sort_index().reset_index(drop=True)

@print_help_on_error
def main(args_list=None):
    parser = argparse.ArgumentParser(description=f"{SCRIPT_NAME} (v{VERSION}) - {DESCRIPTION}")
    parser.add_argument("-b","--bam_file", default="", help="Input BAM file containing aligned RNA-seq reads")
    parser.add_argument("-g", "--gxf",  help='Input GFF/GFF3 file')
    parser.add_argument("-o","--output", default="output",help="Output prefix for files (will add extensions)")
    parser.add_argument("-t","--threads", default=5,type=int, help="Number of threads to use for processing")
    parser.add_argument("-f","--fasta",  help="gene assembly (fasta).")
    parser.add_argument("--force", default=False,action='store_true', help="")
    #args = parser.parse_args(args_list) if args_list else parser.parse_args()
    #args, _ = parser.parse_known_args(args_list) if args_list else parser.parse_args()
    if args_list:
        args, _ = parser.parse_known_args(args_list)
    else:
        args = parser.parse_args()
    #print(args)
    # Define output files
    stringtie_gtf = f"{args.output}_stringtie.gtf"
    output_gff = f"{args.output}_orf.gff"
    output_bed = f"{args.output}.Totalgene.bed"
    output_geneloci_0 = f"{args.output}_r1.geneloci.bed"
    output_geneloci = f"{args.output}_label.genefromRNA.bed"
    
    # Step 1: Run StringTie
    stringtie_gtf=run_stringtie(args.bam_file, args.gxf, stringtie_gtf, args.threads,args.fasta,args.output,args.force)

    # Step 1: Run prediction
    if os.path.exists(output_gff) and not args.force:
        print(f"{output_gff} have already generated, skip stringtie process")
        pass
    elif os.path.exists(f'Quality_info/{output_gff}') and not args.force:
        print(f"{output_gff} have already generated, skip stringtie process")
        output_gff=f"Quality_info/{output_gff}"
        pass
    else:
        #tmpgtf= f"{args.output}_orf.gff" 
        predict_and_write_orfs(
        gtf_path=stringtie_gtf,
        genome_fasta=args.fasta,
        output_gff=output_gff,
        min_aa=50,
        require_atg=True
        )
        
        #gff3_to_gtf(tmpgtf,output_gtf)
       # os.remove(tmpgtf) 
  
    # Step 2.2: Process GTF to BED and get DataFrame
    transcript_df = process_gtf_to_bed(stringtie_gtf, output_bed)   
    filter_df=process_overlap(transcript_df)
    filter_df.to_csv(output_geneloci_0,sep="\t",index=False,header=["Chr","Start","End","Strand","Coverage","Length"])

    regions = []
    # 假设你的数据在一个名为data的DataFrame中
    for row in filter_df.itertuples():
        regions.append(
            overlap_resolver.Region(
                Chr=row.Chr,
                Start=row.Start,
                End=row.End,
                Strand=row.Strand,
                cov=row.cov,
                Length=row.Length
            )
        )
    
    processed_regions = overlap_resolver.process_regions(regions, gap=50)
    result_df = pd.DataFrame([(r.Chr, r.Start, r.End, r.Strand, r.cov, r.Length, r.flag) 
                            for r in processed_regions],
                            columns=['Chr', 'Start', 'End', 'Strand', 'Coverage', 'Length', 'Flag'])

    # Step 3: Generate gene loci BED
    result_df.to_csv(output_geneloci,sep="\t",index=False)

    print("\nProcessing completed successfully")
    print(f"Results saved to:")
    print(f"- GTF file: {stringtie_gtf}")
    print(f"- GTF file: {output_gff}")
    print(f"- Total BED file: {output_bed}")
    print(f"- Gene loci BED: {output_geneloci}")

if __name__ == "__main__":
    try:
        main() 
    except Exception as e:
        parser.print_help()
        print(f"\nError: An unexpected error occurred: {e}")
        print("\nComplete help information:")
        sys.exit(1)