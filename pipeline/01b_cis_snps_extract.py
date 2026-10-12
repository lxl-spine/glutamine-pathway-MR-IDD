# ============================================================
# 步骤 01b · 数据准备：cis±2Mb 区域全量 SNP 提取
# ------------------------------------------------------------
# 功能: 从 BESD 提取基因列表中每个基因的 cis±2Mb 区域全部 SNP
#       （b/se/p；coloc / SMR-HEIDI / 多 IV MR 的输入）
# 输入: data/GTEx_v8_westlake/{组织}.besd/.epi/.esi + gene_list.txt
# 输出: results/cis_snps_full_{组织}.csv
#       （列: gene,probeID,chr,probe_bp,SNP,SNP_bp,A1,A2,b,se,p）
# 用法: python 01b_cis_snps_extract.py <besd_prefix> <gene_list> <out>
#       （对 4 个组织各调用一次）
# ============================================================
# parse_besd_cis.py — 从 BESD 提取指定基因的 cis±2Mb 区域全部 SNP (b/se/p)
# 用法: python parse_besd_cis.py <besd_prefix> <gene_list.txt> <out_csv>
# gene_list.txt: 每行一个 gene symbol
import sys, struct, mmap, math, csv

def parse_esi(path):
    snp_chr, snp_rs, snp_bp, snp_a1, snp_a2 = [], [], [], [], []
    with open(path, 'rt') as f:
        for line in f:
            c = line.rstrip('\n').split('\t')
            if len(c) < 6:
                continue
            snp_chr.append(int(c[0]))
            snp_rs.append(c[1])
            snp_bp.append(int(c[3]))
            snp_a1.append(c[4])
            snp_a2.append(c[5])
    return snp_chr, snp_rs, snp_bp, snp_a1, snp_a2

def parse_epi(path):
    epi_chr, epi_prb, epi_bp, epi_gene = [], [], [], []
    with open(path, 'rt') as f:
        for line in f:
            c = line.rstrip('\n').split('\t')
            if len(c) < 6:
                continue
            epi_chr.append(int(c[0]))
            epi_prb.append(c[1])
            epi_bp.append(int(c[3]))
            epi_gene.append(c[4])
    return epi_chr, epi_prb, epi_bp, epi_gene

def main():
    prefix, genefile, out = sys.argv[1], sys.argv[2], sys.argv[3]
    want = set()
    with open(genefile) as f:
        for line in f:
            g = line.strip()
            if g:
                want.add(g)
    snp_chr, snp_rs, snp_bp, snp_a1, snp_a2 = parse_esi(prefix + '.esi')
    epi_chr, epi_prb, epi_bp, epi_gene = parse_epi(prefix + '.epi')
    nSnp = len(snp_rs)
    nProbe = len(epi_prb)
    # probe 索引: gene -> list of (probe_idx)
    from collections import defaultdict
    gene2idx = defaultdict(list)
    for i, g in enumerate(epi_gene):
        if g in want:
            gene2idx[g].append(i)
    print(f"esi {nSnp}, epi {nProbe}, target genes {len(gene2idx)} found", flush=True)

    with open(prefix + '.besd', 'rb') as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        gflag = struct.unpack_from('<i', mm, 0)[0]
        assert gflag == 3
        valNum = struct.unpack_from('<Q', mm, 0x40)[0]
        colNum = 2 * nProbe + 1
        cols = struct.unpack_from(f'<{colNum}Q', mm, 0x48)
        rowSTART = 0x40 + 8 + colNum * 8
        valSTART = rowSTART + valNum * 4
        import numpy as np
        rowids = np.frombuffer(mm, dtype=np.uint32, count=valNum, offset=rowSTART)
        vals = np.frombuffer(mm, dtype=np.float32, count=valNum, offset=valSTART)

        W = 2000  # kb
        with open(out, 'w', newline='') as o:
            w = csv.writer(o, delimiter='\t')
            w.writerow(['gene', 'probeID', 'chr', 'probe_bp', 'SNP', 'SNP_bp', 'A1', 'A2', 'b', 'se', 'p'])
            n_out = 0
            for g, idxs in gene2idx.items():
                for pi in idxs:
                    pchr = epi_chr[pi]
                    pbp = epi_bp[pi]
                    lo = pbp - W * 1000
                    hi = pbp + W * 1000
                    b_start = cols[2 * pi]
                    b_end = cols[2 * pi + 1]
                    s_start = cols[2 * pi + 1]
                    s_end = cols[2 * pi + 2] if 2 * pi + 2 < colNum else valNum
                    num = b_end - b_start
                    if num == 0:
                        continue
                    rid = rowids[b_start:b_end]
                    bv = vals[b_start:b_end]
                    sv = vals[s_start:s_end]
                    # 只取 cis 窗口内 SNP
                    keep = []
                    for j in range(num):
                        r = int(rid[j])
                        if r < 0 or r >= nSnp:
                            continue
                        if snp_chr[r] == pchr and lo <= snp_bp[r] <= hi:
                            keep.append(j)
                    if not keep:
                        continue
                    zz = np.abs(bv[keep]) / np.maximum(sv[keep], 1e-300)
                    with np.errstate(over='ignore'):
                        pp = 2.0 * (1.0 - 0.5 * (1.0 + np.vectorize(math.erf)(zz / np.sqrt(2.0))))
                    for k, j in enumerate(keep):
                        r = int(rid[j])
                        w.writerow([g, epi_prb[pi], pchr, pbp,
                                    snp_rs[r], snp_bp[r], snp_a1[r], snp_a2[r],
                                    bv[j], sv[j], pp[k]])
                        n_out += 1
                    if n_out % 200000 == 0:
                        print(f"  {n_out} rows", flush=True)
        del rowids, vals
        try:
            mm.close()
        except BufferError:
            pass
    print(f"DONE {n_out} cis-SNP rows -> {out}")

if __name__ == '__main__':
    main()
