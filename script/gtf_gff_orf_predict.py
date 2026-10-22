from typing import Dict, List, Tuple, Optional
import gffutils
import os 
from collections import defaultdict
import sys
#这里可能预测出问题。可能还是用原来的方法比较合适。transcode

def predict_and_write_orfs(
    gtf_path: str,
    genome_fasta: Dict[str, str],
    output_gff: str,
    min_aa: int = 30,
    require_atg: bool = True
    ):
    """
    从GTF预测ORF并输出GFF3文件
    
    Args:
        gtf_path: 输入GTF文件路径
        genome_fasta: 基因组序列文件
        output_gff: 输出GFF3文件路径
        min_aa: 最小ORF长度(AA)
        require_atg: 是否要求起始密码子ATG
    """
    
    # 辅助函数
    def revcomp(s):
        return s.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]
    
    def parse_fasta_to_dict(fasta_path: str) -> Dict[str, str]:
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"FASTA file does not exist: {fasta_path}")

        genome_dict = {}
        current_chr = None
        current_seq = []

        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()  # 移除头尾空白字符
                if not line:  # 跳过空行
                    continue
                
                if line.startswith('>'):  # 染色体头行
                    if current_chr is not None:  # 保存上一个染色体的序列
                        genome_dict[current_chr] = ''.join(current_seq).upper()
                        current_seq = []
                    
                    # 提取染色体名（去掉 '>' 和描述部分）
                    current_chr = line[1:].split()[0]
                else:
                    # 合并多行序列（移除数字/空格）
                    cleaned_line = ''.join(c for c in line if c.isalpha())
                    current_seq.append(cleaned_line)
        
        # 添加最后一个染色体
        if current_chr is not None:
            genome_dict[current_chr] = ''.join(current_seq).upper()

        if not genome_dict:
            raise ValueError("No valid sequences found in the FASTA file")
        return genome_dict

    # ---------- GTF attrs ----------
    def parse_attrs_gtf(attr: str) -> Dict[str, str]:
        d = {}
        for field in attr.strip().split(";"):
            f = field.strip()
            if not f: continue
            parts = f.split(None, 1)
            if len(parts) != 2: continue
            k, v = parts
            d[k] = v.strip().strip('"')
        return d
    # ---------- Mapping ----------
    def map_transcript_span_to_genome_blocks(exons_tx_order: List[Tuple[int,int,int,int]],
                                            t_start0: int, t_end0: int, strand: str
                                            ) -> List[Tuple[int,int]]:
        blocks = []
        for g_s, g_e, t_s, t_e in exons_tx_order:
            ov_s = max(t_start0, t_s)
            ov_e = min(t_end0,   t_e)
            if ov_s >= ov_e: continue
            offset = ov_s - t_s
            length = ov_e - ov_s
            if strand == "+":
                bs = g_s + offset
                be = bs + length - 1
            else:
                be = g_e - offset
                bs = be - length + 1
            blocks.append((bs, be))
        return blocks

    # ---------- ORF ----------
    STOP = {"TAA","TAG","TGA"}
    def find_longest_orf(seq: str, require_atg=True, permissive=False) -> Optional[Tuple[int,int]]:
        best = None
        best_len = 0
        n = len(seq)
        for frame in (0,1,2):
            if require_atg:
                i = frame
                while i + 3 <= n:
                    if seq[i:i+3] == "ATG":
                        k = i + 3
                        while k + 3 <= n and seq[k:k+3] not in STOP:
                            k += 3
                        if k + 3 <= n and seq[k:k+3] in STOP:
                            L = (k + 3) - i
                            if L > best_len: best_len, best = L, (i, k + 3)
                        elif permissive:
                            L = n - i - ((n - i) % 3)
                            if L > best_len: best_len, best = L, (i, i + L)
                    i += 3
            else:
                last = frame
                j = frame
                seen_stop = False
                while j + 3 <= n:
                    if seq[j:j+3] in STOP:
                        seen_stop = True
                        L = j - last
                        L -= (L % 3)
                        if L > best_len: best_len, best = L, (last, last + L)
                        last = j + 3
                    j += 3
                if permissive and (seen_stop is False or last < n):
                    L = n - last - ((n - last) % 3)
                    if L > best_len: best_len, best = L, (last, last + L)
        return best

    def compute_phases(block_lengths_in_tx_order: List[int]) -> List[int]:
        phases, cum = [], 0
        for L in block_lengths_in_tx_order:
            phases.append(0 if cum % 3 == 0 else (3 - (cum % 3)) % 3)
            cum += L
        return phases

# ---------- Main ---------- 修改出来

    genome = parse_fasta_to_dict(genome_fasta)

    # Collect per-transcript structure & flags
    tx_exons: Dict[str, Dict] = {}                 # tid -> {chrom,strand_raw,exons[],gene}
    gene_to_tx: Dict[str, List[str]] = defaultdict(list)
    tx_has_evi: Dict[str, bool] = {}
    gene_has_evi: Dict[str, bool] = defaultdict(bool)

    with open(gtf_path, "r+") as fin:
        for raw in fin:
            if raw.startswith("#") or not raw.strip(): continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 9: continue
            chrom, source, feat, s, e, score, strand, frame, attr = cols
            try:
                s_i = int(s); e_i = int(e)
            except: 
                continue
            attrs = parse_attrs_gtf(attr)
            gid = attrs.get("gene_id")
            tid = attrs.get("transcript_id")

            if not tid: 
                continue
            if tid not in tx_exons:
                tx_exons[tid] = {"chrom":chrom, "strand_raw": strand, "exons":[], "gene":gid}
                if gid: gene_to_tx[gid].append(tid)
            else:
                if gid and tx_exons[tid].get("gene") is None:
                    tx_exons[tid]["gene"] = gid
                    gene_to_tx[gid].append(tid)

            if feat.lower() == "exon":
                tx_exons[tid]["exons"].append((min(s_i,e_i), max(s_i,e_i)))

    # Fill missing genes (rare)
    for tid, info in list(tx_exons.items()):
        if info.get("gene") is None:
            fake = f"{tid}.gene"; info["gene"] = fake; gene_to_tx[fake].append(tid)
            if tx_has_evi.get(tid, False): gene_has_evi[fake] = True

    # Build transcript structures, **robust to '.' strand**
    tx_struct = {}  # tid -> dict with chosen strand, seq, ORF, maps, etc.
    for tid, info in tx_exons.items():
        chrom = info["chrom"]; strand_raw = info["strand_raw"]
        exs = sorted(info["exons"], key=lambda t: t[0])
        if not exs: continue
        chrom_seq = genome.get(chrom)
        if chrom_seq is None:
            sys.exit(f"[error] Chromosome {chrom} not found in FASTA.")

        def make_seq_and_map(ex_order, rc: bool):
            tpos = 0
            ex_map = []
            parts = []
            for g_s, g_e in ex_order:
                frag = chrom_seq[g_s-1:g_e]
                if rc: frag = revcomp(frag)
                parts.append(frag)
                ex_map.append((g_s, g_e, tpos, tpos + len(frag)))
                tpos += len(frag)
            return "".join(parts), ex_map

        # Case 1: known '+' or '-' strand
        if strand_raw in {"+","-"}:
            ex_tx_order = exs if strand_raw == "+" else exs[::-1]
            seq, ex_map = make_seq_and_map(ex_tx_order, rc=(strand_raw == "-"))
            orf = find_longest_orf(seq)
            aa_len = 0 if not orf else (orf[1] - orf[0]) // 3
            strand = strand_raw
        else:
            # Case 2: unknown '.' → evaluate BOTH strands, pick the better ORF
            # plus
            ex_plus = exs
            seq_plus, map_plus = make_seq_and_map(ex_plus, rc=False)
            orf_plus = find_longest_orf(seq_plus)
            aa_plus = 0 if not orf_plus else (orf_plus[1] - orf_plus[0]) // 3
            # minus
            ex_minus = exs[::-1]
            seq_minus, map_minus = make_seq_and_map(ex_minus, rc=True)
            orf_minus = find_longest_orf(seq_minus)
            aa_minus = 0 if not orf_minus else (orf_minus[1] - orf_minus[0]) // 3
            # choose
            if aa_minus > aa_plus:
                strand, seq, ex_tx_order, ex_map, orf, aa_len = "-", seq_minus, ex_minus, map_minus, orf_minus, aa_minus
            else:
                strand, seq, ex_tx_order, ex_map, orf, aa_len = "+", seq_plus, ex_plus, map_plus, orf_plus, aa_plus

        tx_struct[tid] = {
            "chrom":chrom, "strand":strand, "exons":exs,
            "ex_tx_order":ex_tx_order, "ex_tx_map":ex_map,
            "tx_seq":seq, "orf":orf, "aa_len":aa_len
        }

    # Write GFF3
    out = open(output_gff, "w+")
    out.write("##gff-version 3\n")
    filter_min_aa=50
    for gid in sorted(gene_to_tx.keys()):
        tids_all = gene_to_tx[gid]
        # apply ORF length filter
        if filter_min_aa is not None:
            tids = [t for t in tids_all if t in tx_struct and tx_struct[t]["aa_len"] >= filter_min_aa]
        else:
            tids = [t for t in tids_all if t in tx_struct]
        if not tids: continue

        # gene span & strand (if consistent)
        chroms = set(); strands = set()
        g_min = 10**18; g_max = -1
        for t in tids:
            chroms.add(tx_struct[t]["chrom"])
            s = tx_struct[t]["strand"]
            if s in "+-": strands.add(s)
            for xs, xe in tx_struct[t]["exons"]:
                if xs < g_min: g_min = xs
                if xe > g_max: g_max = xe
        chrom = sorted(chroms)[0]
        g_strand = list(strands)[0] if len(strands) == 1 else "."

        # gene line (source=GeneCode-X, optional EVI=RNAseq)
        attrs = [f"ID={gid}"]
        out.write("\t".join([chrom, "GSA", "gene", str(g_min), str(g_max),
                             ".", g_strand, ".", "; ".join(attrs)]) + "\n")

        # transcripts
        for tid in sorted(tids):
            st = tx_struct[tid]
            chrom_t, strand_t = st["chrom"], st["strand"]
            exs = st["exons"]
            tx_start, tx_end = exs[0][0], exs[-1][1]
            out.write("\t".join([chrom_t, "GSA", "mRNA", str(tx_start), str(tx_end),
                                 ".", strand_t, ".", f"ID={tid};Parent={gid}"]) + "\n")

            # exons
            order_map = {ex: i+1 for i, ex in enumerate(st["ex_tx_order"])}
            for gs, ge in sorted(exs, key=lambda t: (t[0], t[1])):
                out.write("\t".join([chrom_t, "GSA", "exon", str(gs), str(ge),
                                     ".", strand_t, ".", f"ID=exon:{tid}:{order_map[(gs,ge)]};Parent={tid}"]) + "\n")

            # CDS/UTR if ORF passes --min-aa
            orf, aa_len = st["orf"], st["aa_len"]
            #print(st)
            if orf and aa_len >= min_aa:
                t0, t1 = orf
                cds_blocks = map_transcript_span_to_genome_blocks(st["ex_tx_map"], t0, t1, strand_t)
                blens = [e - s + 1 for (s, e) in cds_blocks]
                phases = compute_phases(blens)
                cds_out = sorted([(s,e,p) for (s,e),p in zip(cds_blocks, phases)], key=lambda x: (x[0], x[1]))
                for s_cds, e_cds, phase in cds_out:
                    out.write("\t".join([chrom_t, "GSA", "CDS", str(s_cds), str(e_cds),
                                         ".", strand_t, str(phase), f"ID=CDS:{tid};Parent={tid}"]) + "\n")

    out.write("##end\n")
    if out is not sys.stdout:
        out.close()


# 使用示例
if __name__ == '__main__':
    genome = {}  # 从FASTA加载 {'chr1': 'ATGCGT...', ...}
    
    predict_and_write_orfs(
        gtf_path="input.gtf",
        genome_fasta=genome,
        output_gff="output.gff3",
        min_aa=50,
        require_atg=True
    )
