# USAGE:
# Rscript simulate_sv_genome.R te_fa genome_fa n_ins n_grad
# Rscript simulate_sv_genome.R te_ccs.fa origion_genome.fa 1 2

library(Rsamtools)
library(regioneR)

# 从命令行获取参数 ---- te ccs & ref genome
args <- commandArgs(trailingOnly = TRUE)
n_ins = as.integer(args[3]) # 每种transposon生成的insertion数量
n_grad = as.integer(args[4]) # 长度划分梯度的数目

# 读取transposon consensus序列，以及genome
te_fa <- open(FaFile(args[1]))
te_idx <- scanFaIndex(te_fa)
te_idx <- narrow(te_idx, end = -2)
genome_fa <- open(FaFile(args[2]))
genome_idx <- scanFaIndex(genome_fa)
genome_idx <- narrow(genome_idx, end = -2)

# 首先生成tsd区间
tsd_idx <- createRandomRegions(nregions = n_ins*length(te_idx)*n_grad, length.mean = 8, length.sd = 1, genome = genome_idx)
tsd_idx = restrict(tsd_idx, start = 1)

# 对于每种transposon，生成插入片段，与对应的tsd拼接起来，并且输出到文件当中
k=0
grad = 0.9/n_grad
for (i in 1:length(te_idx)) {
  for (j in 1:n_grad) {
    l = k+j
    start_idx = (l-1)*n_ins + 1
    end_idx = l*n_ins
    
    # 随机生成TE insert片段
    sd = 0.1*grad*width(te_idx[i])
    tmp_te_idx = createRandomRegions(nregions = n_ins, length.mean = round(j*grad*width(te_idx[i])), length.sd = sd, genome = te_idx[i], non.overlapping = FALSE)
    tmp_te_idx = restrict(tmp_te_idx, start = 1)
    tmp_tsd_idx = tsd_idx[start_idx:end_idx]
    
    # 随机设定TE片段是forward inserted还是reverse inserted
    tmp_tsd_df = as.data.frame(tmp_tsd_idx)
    tmp_te_df = as.data.frame(tmp_te_idx)
    tmp_te_df$strand <- sample(c("-", "+"), size = nrow(tmp_te_df), replace = TRUE)
    
    # 获取TE片段以及TSD片段，编辑每个insertion的name
    tmp_te_seq = as.character(getSeq(te_fa, makeGRangesFromDataFrame(tmp_te_df)))
    tmp_tsd_seq = as.character(getSeq(genome_fa, tmp_tsd_idx))
    tmp_names = paste(tmp_tsd_df$seqnames, tmp_tsd_df$end-1, tmp_tsd_df$end, tmp_te_df$seqnames, tmp_te_df$strand, sep = "_")
    
    # 将一轮循环种产生的insertion信息输出到文件当中
    tmp_ins_df = data.frame(chrom = tmp_tsd_df$seqnames, start = tmp_tsd_df$end-1, end = tmp_tsd_df$end,
                            type = rep("insertion", n_ins), ins_seq = paste(tmp_te_seq, tmp_tsd_seq, sep = ""), add_len = rep(0, n_ins),
                            name = tmp_names, strand = tmp_te_df$strand, TE = tmp_te_df$seqnames,
                            te_start = tmp_te_df$start-1, te_len = tmp_te_df$width, te_seq = tmp_te_seq,
                            tsd_start = tmp_tsd_df$start-1, tsd_len = tmp_tsd_df$width, tsd_seq = tmp_tsd_seq)
    
    # write.table(tmp_ins_df, file = "./simulated_sv.bed", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    if (l==1) {write.table(tmp_ins_df, file = "./simulated_sv.summary", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE)}
    else {write.table(tmp_ins_df, file = "./simulated_sv.summary", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)}
  }
  k = k + n_grad
}