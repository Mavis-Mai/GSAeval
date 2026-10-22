import sys
import numpy as np
from pathlib import Path
from Bio.Align import PairwiseAligner
from Bio.Seq import translate
from Bio.pairwise2 import format_alignment
sys.path.insert(0, str(Path(__file__).parent))
try:
    from dtw_match import find_unique_top_matches,len_matrix,phase_matrix
except ImportError as e:
    #print(f"Error: 必需的模块 overlap_resolver 未找到。请确保它在脚本所在目录的父目录中。")
    print(f"Error text: {e}")
    sys.exit(1)

##增加动态匹配的机制
def getMath_results(cds_a,cds_b):
    cds_a_,cds_b_=cds_a["cds"],cds_b["cds"]
    seq_a,seq_b=cds_a["seq"],cds_b["seq"]
    ## pahse_result
    results,seqlist=[],[]
    
    phase1=[key[0] for key in cds_a_]
    phase2=[key[0] for key in cds_b_]
    cost_matrix=phase_matrix(phase1,phase2)
    optimal_paths_1 = find_unique_top_matches(cost_matrix)
    
    for group in optimal_paths_1:
        pahse_result,seq_result=[],[]
        for pairs in group:
            #optimal_paths格式如下
            #([(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6)], 3.8175711346821175)
            idx_a, idx_b = pairs 
            pahse_result.append((cds_a_[idx_a],cds_b_[idx_b]))
            seq_result.append((seq_a[idx_a],seq_b[idx_b]))
            
        results.append(pahse_result)
        seqlist.append(seq_result)

    ##cds_result
    length1=[key[1] for key in cds_a_]
    length2=[key[1] for key in cds_b_]
    cost_matrix=len_matrix(length1,length2)

    optimal_paths_2 = find_unique_top_matches(cost_matrix)

    for group in optimal_paths_2:
        len_result,seq_result=[],[]
        for pairs in group:
            idx_a, idx_b = pairs 
            len_result.append((cds_a_[idx_a],cds_b_[idx_b]))
            seq_result.append((seq_a[idx_a],seq_b[idx_b]))
        results.append(len_result)
        seqlist.append(seq_result)
    
    return results,seqlist

           

def calc_csq_(cds_a,cds_b,results,seqlist,type):
    weightdict={
        'plant':{"cns":0.238912,
        "iss":0.214510,
        "pms":0.202859,
        "cls":0.183402,
        "tls":0.160316  
        },
        'animal':{"cns":0.046990,
        "iss":0.146535,
        "pms":0.317778,
        "cls":0.331030,
        "tls":0.157667  
        }
    }

    cds_a_,cds_b_=cds_a["cds"],cds_b["cds"]
    len_a,len_b=cds_a["cds_len"],cds_b["cds_len"]
    #seq_a,seq_b=cds_a["seq"],cds_b["seq"]
    #默认introna>intronb
    total_csq=[]
    total_info=[]
    total_AED=[]

    #遍历每一个reasult
    for index_ in range(len(results)):
        group=results[index_]
        if group !=[]:
            seq_group=seqlist[index_]
            length_bias=1-abs(len_a-len_b)/max((len_a+len_b)/2,1) 
            if length_bias<0:
                length_bias=0

            matched_pairs = []
            for pair in group:
                pa, la, donor_a, acceptor_a = pair[0]   
                pb, lb, donor_b, acceptor_b= pair[1]   
                matched_pairs.append(( la, lb, pa, pb, donor_a, donor_b, acceptor_a, acceptor_b))

            # print(len(matched_pairs))
            # 1. CDS_count 不一定大于0，为了避免负分出现，分拆原来的指标
            IA, IB = len(cds_a_), len(cds_b_)
            mathLen=len(group)
            count_match = 1 - abs(IA - IB) / max((IA+IB)/2, 1) 
            if count_match < 0:
                count_match= 0

            # 2. PhaseMatch 
            phase_match=[1 if pa == pb else 0 for _, _, pa,pb, _, _ , _, _,in matched_pairs]
            #考虑上罚分：  
            phase_match = 1- (max(IA,IB)-sum(phase_match))/(max(IA,IB)) if phase_match else 0
            
            # 3. LenSim
            len_diffs = []
            for la, lb, _, _,  _, _, _, _,in matched_pairs:
                la, lb = max(0, int(la)), max(0, int(lb))  # 强制非负整数
                if la + lb + 1 == 0:  # 防御性编程
                    continue
                norm_diff = abs(la - lb) / max((la + lb)/2, 1) 
                len_diffs.append(norm_diff)  # 强制 ≤1
                
            #还是考虑上罚分吧：  
            len_diffs=len_diffs+(max(IA,IB)-mathLen)*[1]
            len_sim = 1 - (sum(len_diffs)/len(len_diffs)) if len_diffs else 0
            len_sim = max(0.0, min(1.0, len_sim))  # 最终钳制到[0,1]
    
            # 4. intron_site
            donor_score=[1 if da == db else 0 for _, _, _, _, da, db, _,_ in matched_pairs]
            acceptor_socre=[1 if aa == ab else 0 for _, _, _, _,  _,_, aa, ab in matched_pairs]
            #还是考虑上罚分吧：  
            intronSite_match=1- (2*max(IA,IB)-(sum(donor_score)+sum(acceptor_socre)))/(2*max(IA,IB))

            # 计算CSQ
            
            csq = weightdict[type]['cns']*count_match + weightdict[type]['pms']*phase_match + weightdict[type]['cls']*len_sim + weightdict[type]['iss']*intronSite_match + weightdict[type]['tls']*length_bias 

            total_csq.append(csq)
            total_info.append([IA, IB ,count_match,  phase_match,len_sim, intronSite_match,length_bias]) #要再考虑一下权重问题
    return total_csq,total_info

def calc_csq(cdss_a, cdss_b,type):
    """计算两个基因之间的qSCS分数"""
    results,seqlist=getMath_results(cdss_a,cdss_b)
    total_csq,total_info=calc_csq_(cdss_a,cdss_b,results,seqlist,type)
    csq_max=max(total_csq)
    IA, IB ,count_match,phase_match,len_sim,intronSite_match,length_bias =total_info[total_csq.index(csq_max)]
    return {
            'cdsA_number':IA,
            'cdsB_number':IB,
            'Coding-exon_number_score': round(count_match, 3),
            'Phase_match_score': round(phase_match, 3),
            'Coding-exon_length_score': round(len_sim, 3),
            'Intron_site_score':round(intronSite_match, 3),
            'Total_length_score':round(length_bias,3),
            'CSQ': round(csq_max, 3)
        }


    