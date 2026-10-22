import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import sys
import os

def plot_boxplot_and_get_threshold(df, output, dir_path, column_name='CSQ'):
    """
    绘制指定列的箱线图并返回3/4分位数阈值
    
    参数:
    df -- 输入的数据框(pandas.DataFrame)
    column_name -- 要分析的列名(默认为'CSQ')
    
    返回:
    该列的3/4分位数(75百分位数)阈值
    """
    # 绘制箱线图
    #plt.style.use('seaborn-v0_8-pastel') 
    plt.figure(figsize=(2.5, 5))
    boxprops = dict(linestyle='-', linewidth=1.5, color='black')  # 箱子线条
    flierprops = dict(marker='o', markersize=6, markerfacecolor='#b8b8b8', alpha=0.5)  # 异常点样式
    whiskerprops = dict(linestyle='--', linewidth=1, color='#7f8c8d')  # 须线样式
    plt.boxplot(
        df[column_name].dropna(),  # 自动忽略NaN值
        widths=0.3,  # 控制宽度
        boxprops=boxprops,
        flierprops=flierprops,
        whiskerprops=whiskerprops,
        medianprops=dict(color='#2c3e50', linewidth=2)  # 中位数线
    )
    x_jitter = np.random.normal(1, 0.05, size=len(df[column_name]))
    plt.scatter(
        x_jitter,
        df[column_name],
        color='orange',
        alpha=0.5,
        s=20,  # 点大小
        edgecolor='black',
        linewidth=0.5,
        label='All Data Points'
    )

  #  plt.title(f'{column_name} distribution box-plot')
   
    plt.ylabel('CSQ value')
    
    
    # 计算并返回3/4分位数
    threshold = df[column_name].quantile(0.25)
    plt.axhline(y=threshold, color='r', linestyle='--', 
                label=f'25%: {threshold:.2f}')

    plt.text(x=0.3,y=threshold-.005,color='r',s=f' {threshold:.2f}')
    ouputfile=os.path.join(dir_path,f'{output}.png')
    plt.savefig(ouputfile)
    

    return threshold,list(df[column_name].dropna())


def main():
    parser = argparse.ArgumentParser(
        description='Take the CSQ threshold of BUSCO homologous set ',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--input', required=True, help='Input file path')
    parser.add_argument('-o', '--output', help='ouput prefix', default='output')
    args = parser.parse_args()
    
    dir_path = os.path.dirname(args.input)

    df=pd.read_csv(args.input,header=0,sep="\t")
    threshold,qscs=plot_boxplot_and_get_threshold(df, args.output,dir_path, column_name='CSQ',)
    qscs=[str(i) for i in qscs]
    exon_m=round(np.array(df["cdsA_number"]).mean(),3)
    qscs_num=str(len(df["cdsB_number"]))

    outputfile=os.path.join(dir_path,f"{args.output}_threshold.txt")
    with open(outputfile,'w+')as f:
        f.write('busco\tgenepairs_number\texon_number\tthreshold\tcsq\n')
        f.write(f'{args.output}\t{qscs_num}\t{str(exon_m)}\t{threshold}\t{",".join(qscs)}\n')

    sys.exit(0)

if __name__ == '__main__':
    main()
