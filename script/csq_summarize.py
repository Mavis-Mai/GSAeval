import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict
import matplotlib.transforms as transforms
from scipy import  stats


def plot_pie(df,output,format,lib):
    expected_labels = ['Good', 'Moderate', 'Poor']
    existing_data = df.set_index('label').to_dict('index')

    # 3. 构造新的数据：如果存在则读取，如果缺失则补充为 0
    full_data = []
    for label in expected_labels:
        if label in existing_data:
            full_data.append({
                'label': label,
                'count': existing_data[label].get('count', 0),
                'frequency(%)': existing_data[label].get('frequency(%)', 0.0)
            })
        else:
            # 缺失的类别，自动补全 count 为 0，frequency(%) 为 0.0
            full_data.append({
                'label': label,
                'count': 0,
                'frequency(%)': 0.0
            })

    colors = ['#EF888A','#b8b8b8','#45465e'] 
    plot_data = pd.DataFrame(full_data)

    # 绘制饼图
    plt.figure(figsize=(9, 6))
    wedges, texts, autotexts = plt.pie(
        plot_data['frequency(%)'],
        labels=plot_data['label'],
        explode=[0.02] * len(plot_data),
        colors=colors,
        autopct="%1.2f%%",
        startangle=90,
        textprops={'fontsize': 12},
        labeldistance=1.2,
        pctdistance=0.8 
    )
    #plt.setp(autotexts, size=12)
    #plt.title(title, fontsize=15, pad=20)
    # 使用 adjustText 调整标签

    legend_labels = [f"{l}: {c}" for l, c in zip(plot_data['label'], plot_data['count'])]
    missingcount=plot_data.loc[plot_data['label']=="Missing"]['count']
    legend_labels.append(f"Missing: {missingcount}")  
    plt.legend(
        title=f"Database: {lib}",
        loc="upper left",
        bbox_to_anchor=(1, 0, 0.5, 1),
        fontsize=12,  # 统一设置图例字体大小
        title_fontsize=14,  # 设置图例标题字体大小
        labels=legend_labels
    )
    plt.tight_layout()
    plt.subplots_adjust(left=0.1, right=0.7) 

    outputfile=f'{output}_pie.{format}'
    outputfile2=f'{output}_pie.png'
    plt.savefig(outputfile)
    plt.savefig(outputfile2)


def main(args_list=None):
    import argparse
    parser = argparse.ArgumentParser(description='Statistically analyze and plot the CSQ information')
    parser.add_argument('-i', '--CSQ_file', required=True, help='CSQ information of each BUSCO of the species')
    parser.add_argument('-t', '--threshold_table', required=True, help='threshold information of each BUSCO')
    parser.add_argument('--fig_format',default="svg", help='The selection of graphic format. such as:svg, png, pdf')
    parser.add_argument('--csq_threshold',default=0.90, help='a CSQ experience threshold')
    parser.add_argument('--miss', help='missing busco list')
    parser.add_argument('-b', '--library', help='The highly conserved BUSCO gene library (eudicots/monocot/embryophyta)', metavar='LIB' )
    parser.add_argument('-f', '--fasta', help='Input genome FASTA file (required for analysis mode)',  metavar='FILE')
    parser.add_argument('-o', '--output', default='output', 
                       help='output file prefiex (default: output)')

    #args = parser.parse_args()
    args = parser.parse_args(args_list) if args_list else parser.parse_args()
    database_name=args.library
    df_=pd.read_csv(args.CSQ_file, sep="\t",header=0,index_col=False) 
    df_=df_.fillna(0)
    result = df_.loc[
        df_.groupby(["BUSCO", "GeneA"])["CSQ"].idxmax()
    ].reset_index(drop=True)

    df_thresthold=pd.read_csv(args.threshold_table, sep="\t",header=0,index_col=False) 

    labellsit,q_value=[],[]

    for i in range(len(result)):
        busco_id = result.iloc[i]["BUSCO"]
        #print(busco_id)
        CSQ_tmp = result.iloc[i]["CSQ"]
        try:
            thresthold=df_thresthold.loc[df_thresthold["busco"]==busco_id]["threshold"].values[0]
        except IndexError:
            print(busco_id)
        #print(thresthold)
        if CSQ_tmp >= thresthold:
            label='Good'
        elif CSQ_tmp >= args.csq_threshold:
            label='Good'
        elif CSQ_tmp <0.6 :
            label= 'Poor'
        else:
            label= 'Moderate'   
        labellsit.append(label)

    result["label"]=labellsit

    # 保存结果
    result.to_csv(f'{args.output}_label_CSQ.txt', sep="\t", index=False)
    # 打印结果
    print(result)
    result_busco=result.loc[
                result.groupby(["BUSCO"])["CSQ"].idxmax()
                ].reset_index(drop=True)
    result_busco.to_csv(f'{args.output}_max_CSQ.txt', sep="\t", index=False)

    ##开始画图
    df_counts = result_busco["label"].value_counts().rename('count').reset_index()
    df_counts.columns = ['label', 'count']
    filename=f'{args.output}.CSQ_stat.txt'
    if args.miss:
        missdf=pd.read_csv(args.miss,names=["buscoid"])
        misslist=list(missdf["buscoid"])
        df_counts["frequency(%)"]=round(df_counts["count"]/df_counts["count"].sum()*100,2)
        df_counts.loc[len(df_counts)] = ['Missing', len(misslist),None]
        df_counts.to_csv(filename,sep="\t",index=False)
        print(df_counts)
        plot_pie(df_counts, args.output+"_busco",args.fig_format,database_name)
    else:        
        df_counts["frequency(%)"]=round(df_counts["count"]/df_counts["count"].sum()*100,2)
        df_counts.to_csv(filename,sep="\t",index=False)
        print(df_counts)

        #plot_frequency_distribution(freq_dist_busco, args.output+"_busco")
        plot_pie(df_counts, args.output+"_busco",args.fig_format,database_name)

    with open(filename,'a+') as f:
        f.write(f"The database used in CSQ mode:{database_name}\n")
        f.write(f"The fasta used in CSQ mode:{args.fasta}\n")
        

if __name__ == '__main__':
    main()

