# ============================================================
# 步骤 01 · 数据准备：GTEx v8 四组织 lead cis-eQTL 提取
# ------------------------------------------------------------
# 功能: 从 Westlake SMR BESD（与 GTEx v8 allpairs 同源）解析每组织
#       每基因 top cis-eQTL（cis±2Mb 窗口内 p 最小 SNP）→ 工具表
# 输入: data/GTEx_v8_westlake/{组织}.besd/.epi/.esi
# 输出: results/lead_eqtl_{Adipose,Muscle,Nerve,WholeBlood}.csv
#       （每基因一行 lead SNP；列: probeID,gene,chr,probe_bp,SNP,
#         SNP_bp,A1,A2,freq,b,se,p,z）
# 用法: python 01_data_prep_lead_eqtl.py <besd_prefix> <out_csv>
#       （对 4 个组织各调用一次）
# ============================================================
# parse_besd.py — 解析 SMR sparse BESD v3，提取每基因 top cis-eQTL (lead eQTL)
# 用法: python parse_besd.py <besd_prefix> <out_csv>
# 输出: probeID, gene, chr, probe_bp, SNP, SNP_bp, A1, A2, freq, b, se, p, z
import sys, struct, mmap, math

def parse_esi(path):
    snp_chr, snp_rs, snp_bp, snp_a1, snp_a2, snp_freq = [], [], [], [], [], []
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
            snp_freq.append(float(c[6]) if c[6] not in ('NA', '') else -9.0)
    return snp_chr, snp_rs, snp_bp, snp_a1, snp_a2, snp_freq

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
    prefix = sys.argv[1]
    out = sys.argv[2]
    esi = parse_esi(prefix + '.esi')
    epi = parse_epi(prefix + '.epi')
    snp_chr, snp_rs, snp_bp, snp_a1, snp_a2, snp_freq = esi
    epi_chr, epi_prb, epi_bp, epi_gene = epi
    nSnp = len(snp_rs)
    nProbe = len(epi_prb)
    print(f"esi {nSnp} SNPs, epi {nProbe} probes", flush=True)

    with open(prefix + '.besd', 'rb') as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        gflag = struct.unpack_from('<i', mm, 0)[0]
        assert gflag == 3, f"not sparse v3: {gflag}"
        valNum = struct.unpack_from('<Q', mm, 0x40)[0]
        colNum = 2 * nProbe + 1
        cols_off = 0x48
        cols = struct.unpack_from(f'<{colNum}Q', mm, cols_off)
        rowSTART = 0x40 + 8 + colNum * 8
        valSTART = rowSTART + valNum * 4
        print(f"gflag={gflag} valNum={valNum} rowSTART={rowSTART} valSTART={valSTART}", flush=True)
        # 校验文件大小
        size = mm.size()
        expect = valSTART + valNum * 4
        print(f"file {size} expect {expect} ok={size == expect}", flush=True)
        assert size == expect, "file size mismatch"

        import numpy as np
        rowids = np.frombuffer(mm, dtype=np.uint32, count=valNum, offset=rowSTART)
        vals = np.frombuffer(mm, dtype=np.float32, count=valNum, offset=valSTART)

        out_rows = []
        n_ok = 0
        for i in range(nProbe):
            b_start = cols[2 * i]
            b_end = cols[2 * i + 1]
            s_start = cols[2 * i + 1]
            s_end = cols[2 * i + 2] if 2 * i + 2 < colNum else valNum
            num = b_end - b_start
            if num == 0:
                continue
            rid = rowids[b_start:b_end]
            bv = vals[b_start:b_end]
            sv = vals[s_start:s_end]
            if len(sv) < num:
                continue
            # p 值
            z = np.abs(bv) / np.maximum(sv, 1e-300)
            with np.errstate(over='ignore'):
                p = 2.0 * (1.0 - 0.5 * (1.0 + np.vectorize(math.erf)(z / np.sqrt(2.0))))
            j = int(np.argmin(p))
            r = rid[j]
            if r < 0 or r >= nSnp:
                continue
            out_rows.append((epi_prb[i], epi_gene[i], epi_chr[i], epi_bp[i],
                             snp_rs[r], snp_bp[r], snp_a1[r], snp_a2[r],
                             snp_freq[r], bv[j], sv[j], p[j], z[j]))
            n_ok += 1
            if n_ok % 5000 == 0:
                print(f"  {n_ok} genes done", flush=True)
        # 先释放 numpy 视图再 close mmap
        del rowids, vals

    with open(out, 'w') as o:
        o.write("probeID\tgene\tchr\tprobe_bp\tSNP\tSNP_bp\tA1\tA2\tfreq\tb\tse\tp\tz\n")
        for r in out_rows:
            o.write('\t'.join(str(x) for x in r) + '\n')
    try:
        mm.close()
    except BufferError:
        pass
    print(f"DONE {n_ok} genes with lead eQTL -> {out}")

if __name__ == '__main__':
    main()
