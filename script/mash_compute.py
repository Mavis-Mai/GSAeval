#!/usr/bin/env python3
"""
主程序集成专用版 - 彻底解决spawn.Pool调用问题
"""
import os
import sys
import argparse
import subprocess
from pathlib import Path
from multiprocessing import get_context

class MashIntegrator:
    """专为主程序集成设计的处理器"""
    
    # 静态方法解决pickle问题
    @staticmethod
    def _run_mash(cmd):
        """安全执行Mash命令"""
        try:
            subprocess.run(cmd, check=True, 
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
            return True
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Mash error in : {e.stderr.decode().strip()}")

    @staticmethod
    def sketch_worker(args):
        """工作函数必须完全自包含"""
        file_path, out_dir, k, s, m = args
        try:
            output = Path(out_dir)/f"{Path(file_path).stem}.msh"
            MashIntegrator._run_mash([
                "mash", "sketch",
                "-k", str(k),
                "-s", str(s),
                "-o", str(output),
                "-m", str(m),
                str(file_path)
            ])
            return str(output)
        except Exception as e:
            print(f"[Worker Error] {file_path}: {str(e)}", file=sys.stderr)
            return None

def execute_mash_analysis(args):
    """主程序调用入口点"""
    try:
        # 初始化环境
        sketch_dir = Path("mash_sketches")
        sketch_dir.mkdir(exist_ok=True)
        
        # 获取绝对路径文件列表
        input_files = [
            str(f.resolve()) for f in Path(args.input).glob("*")
            if f.suffix.lower() in {'.fasta','.fa','.fastq','.fq','.gz'}
        ]
        
        print(f"Starting to switch the mash format of {len(input_files)} samples")
        
        # 关键修正：使用独立进程上下文
        ctx = get_context('spawn')
        with ctx.Pool(args.processes) as pool:
            tasks = [(f, "mash_sketches", args.kmer, args.s, args.cov) 
                    for f in input_files]
            
            # 修正：使用pool.map替代imap保证稳定性
            results = pool.map(MashIntegrator.sketch_worker, tasks)
            
            # 过滤有效结果
            msh_files = [f for f in results if f is not None]
            if len(msh_files) < 2:
                raise ValueError("not enough the mash files for calculating the distance")
            
            # 生成距离矩阵（必须在主进程执行）
            dist_file = f"{args.output}.dist"
            with open(dist_file, 'w') as f:
                subprocess.run(["mash", "triangle"] + msh_files,
                             check=True, stdout=f)
            
            print(f"✅ finish the analysis, and the results save in {dist_file}")
            return True
            
    except Exception as e:
        print(f"💥 Error: {str(e)}", file=sys.stderr)
        raise

if __name__ == "__main__":
    # 命令行接口
    parser = argparse.ArgumentParser(description='Mash calculation',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", required=True, help=f"input diretory with fastq")
    parser.add_argument("-k", "--kmer", type=int, default=21, help=f"k-mer size (bp)")
    parser.add_argument("-s", type=int, default=1000000, help=f"sketch size")
    parser.add_argument("-m", "--cov", type=int, default=5, help="minimum depth")
    parser.add_argument("-p", "--processes", type=int, default=4, help=f"CPU number")
    parser.add_argument("-o", "--output", default="mash_out", help=f"output prefix")
    args = parser.parse_args()
    try:
        execute_mash_analysis(args)
    except Exception:
        sys.exit(1)
