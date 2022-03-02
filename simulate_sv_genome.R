# USAGE:
# Rscript simulate_sv_genome.R te_fa genome_fa n t
# Rscript simulate_sv_genome.R te_ccs.fa origion_genome.fa 12 4
library(Rsamtools)
library(regioneR)


seq_mutate <- function(idx, te_seq, te_idx, te_fa, genome_idx, genome_fa) {
  
  # 随机生成一些突变的位点以及类型
  n_var = sample(1:3, size = 1, prob = c(0.8, 0.1, 0.1))
  var_type = sample(1:4, size = n_var, prob = rep(0.25, 4))
  if (width(idx) <= 1000) {
    var_idx = createRandomRegions(nregions=n_var, length.mean=0.1*width(idx), length.sd=1, genome=idx, non.overlapping=TRUE)
  } else {
    var_idx = createRandomRegions(nregions=n_var, length.mean=100, length.sd=10, genome=idx, non.overlapping=TRUE)
  }
  var_idx = data.frame(var_idx)[order(data.frame(var_idx)$end, decreasing=TRUE),]
  
  # 将te_idx, genome_idx汇总为ref_idx并且对ref_idx进行一些限制
  ref_idx = makeGRangesFromDataFrame(rbind(as.data.frame(te_idx), as.data.frame(genome_idx)))
  ref_idx = narrow(ref_idx, end = -2)
  
  # 从var_idx中筛除一些可能导致bug的突变位点，若全部被筛除，则返回原本的te_seq
  var_idx = var_idx[var_idx$start>0 & var_idx$width>0 & var_idx$end<width(idx),]
  n_var = nrow(var_idx)
  if (nrow(var_idx) == 0) {return(te_seq)}
  
  for (i in 1:nrow(var_idx)) {
    if (var_type[i] == 1) {         # deletion
      del_seq = subseq(te_seq, var_idx[i,]$start, var_idx[i,]$end)
      te_seq = paste(subseq(te_seq, 1, var_idx[i,]$start-1), subseq(te_seq, var_idx[i,]$end+1, nchar(te_seq)), sep = "")
      
      tmp_df <- data.frame(chr=var_idx[i,]$seqnames, start=var_idx[i,]$start-1, end=var_idx[i,]$end,
                           type="deletion", sequence=del_seq, strand=var_idx[i,]$strand)
      write.table(tmp_df, file = "./interal_mutation", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
    
    if (var_type[i] == 2) {         # invertion
      rc_seq = as.character(reverseComplement(DNAString(subseq(te_seq, var_idx[i,]$start, var_idx[i,]$end))))
      te_seq = paste(subseq(te_seq, 1, var_idx[i,]$start-1), rc_seq, subseq(te_seq, var_idx[i,]$end+1, nchar(te_seq)), sep = "")
      
      tmp_df <- data.frame(chr=var_idx[i,]$seqnames, start=var_idx[i,]$start-1, end=var_idx[i,]$end,
                           type="invertion", sequence=rc_seq, strand=var_idx[i,]$strand)
      write.table(tmp_df, file = "./interal_mutation", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
    
    if (var_type[i] == 3) {         # duplication
      dup_seq = paste(subseq(te_seq, var_idx[i,]$start, var_idx[i,]$end), subseq(te_seq, var_idx[i,]$start, var_idx[i,]$end), sep = "")
      te_seq = paste(subseq(te_seq, 1, var_idx[i,]$start-1), dup_seq, subseq(te_seq, var_idx[i,]$end+1, nchar(te_seq)), sep = "")
      
      tmp_df <- data.frame(chr=var_idx[i,]$seqnames, start=var_idx[i,]$start-1, end=var_idx[i,]$end,
                           type="duplication", sequence=dup_seq, strand=var_idx[i,]$strand)
      write.table(tmp_df, file = "./interal_mutation", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
    
    if (var_type[i] == 4) {         # insertion
      ins_idx = createRandomRegions(nregions=1, length.mean=0.1*width(idx), length.sd=1, genome=ref_idx)
      ins_idx = ins_idx[start(ins_idx)>=1]
      if (length(ins_idx)==0) {next}
      if (as.character(seqnames(ins_idx)) %in% seqnames(te_idx)) {
        ins_seq = as.character(getSeq(te_fa, ins_idx))
      } else {
        ins_seq = as.character(getSeq(genome_fa, ins_idx))
      }
      te_seq = paste(subseq(te_seq, 1, var_idx[i,]$start), ins_seq, subseq(te_seq, var_idx[i,]$start+1, nchar(te_seq)), sep = "")
      
      tmp_df <- data.frame(chr=var_idx[i,]$seqnames, start=var_idx[i,]$start-1, end=var_idx[i,]$start,
                           type="insertion", sequence=ins_seq, strand=var_idx[i,]$strand)
      write.table(tmp_df, file = "./interal_mutation", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
  }
  
  return(te_seq)
}


simulate_ins <- function(te_fa, te_idx, genome_fa, genome_idx, tsd_idx, n, t) {
  
  # 初始化中间变量k，计算长度梯度grad，计算每种transposon在每个梯度应产生的insertion number -- n_ins
  k=0
  grad = 1/t
  n_ins = round(0.5*n/(t-1))
  
  # 对于每种transposon，生成长度随机的插入片段，与对应的tsd拼接起来，并且输出到文件当中
  for (i in 1:length(te_idx)) {
    for (j in 1:(t-1)) {
      l = k+j
      start_idx = (l-1)*n_ins + 1
      end_idx = l*n_ins
      
      # 随机生成TE insert片段
      sd = 0.1*grad*width(te_idx[i])
      tmp_te_idx = createRandomRegions(nregions = n_ins, length.mean = round(j*grad*width(te_idx[i])), length.sd = sd, genome = te_idx[i], non.overlapping = FALSE)
      tmp_te_idx = restrict(tmp_te_idx, start = 1)
      tmp_tsd_idx = tsd_idx[start_idx:end_idx]
      
      # 随机设定TE片段是forward inserted还是reverse inserted
      tmp_tsd_df = as.data.frame(tmp_tsd_idx, row.names = NULL)
      tmp_te_df = as.data.frame(tmp_te_idx, row.names = NULL)
      tmp_te_df$strand <- sample(c("-", "+"), size = nrow(tmp_te_df), replace = TRUE)
      
      # 获取TE片段以及TSD片段，编辑每个insertion的name
      tmp_tsd_seq = as.character(getSeq(genome_fa, tmp_tsd_idx))
      tmp_te_seq = as.character(getSeq(te_fa, makeGRangesFromDataFrame(tmp_te_df)))
      tmp_names = paste(tmp_tsd_df$seqnames, tmp_tsd_df$end-1, tmp_tsd_df$end, tmp_te_df$seqnames, tmp_te_df$strand, sep = "_")
      names(tmp_te_seq) <- tmp_names
      
      # 向TE片段中人为添加一些突变
      if (l==1) {writeLines(paste("chr", "start", "end", "type", "sequence", "strand", sep="\t"), "./interal_mutation", sep="\n")}
      origin_seq = getSeq(te_fa, makeGRangesFromDataFrame(tmp_te_df))
      names(origin_seq) = tmp_names
      origin_idx = makeGRangesFromDataFrame(data.frame(chr=origin_seq@ranges@NAMES, start=rep(1,length(origin_seq)), end=width(origin_seq)))
      origin_seq = as.character(origin_seq)
      for (x in 1:length(tmp_te_seq)) {
        if (sample(c(TRUE, FALSE), size = 1, prob = c(0.2, 0.8))) {
          tmp_te_seq[x] = seq_mutate(idx=origin_idx[x], te_seq=tmp_te_seq[x], te_idx, te_fa, genome_idx, genome_fa)
        }
      }
      
      # 将一轮循环种产生的insertion信息输出到文件当中
      tmp_ins_df = data.frame(chrom = tmp_tsd_df$seqnames, start = tmp_tsd_df$end-1, end = tmp_tsd_df$end,
                              type = rep("insertion", n_ins), ins_seq = paste(tmp_te_seq, tmp_tsd_seq, sep = ""),
                              add_len = sample(c(0,1), n_ins, prob=c(0.8, 0.2), replace=TRUE), name = tmp_names, strand = tmp_te_df$strand, TE = tmp_te_df$seqnames,
                              te_start = tmp_te_df$start-1, te_end = tmp_te_df$end, te_seq = tmp_te_seq,
                              tsd_start = tmp_tsd_df$start-1, tsd_len = tmp_tsd_df$end, tsd_seq = tmp_tsd_seq,
                              origin_seq = origin_seq, FullLength = rep("False", n_ins))
      
      if (l==1) {write.table(tmp_ins_df, file = "./simulated_sv.summary", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE)}
      else {write.table(tmp_ins_df, file = "./simulated_sv.summary", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)}
    }
    k = k + (t-1)
  }
}


simulate_ins_full <- function(te_fa, te_idx, genome_fa, genome_idx, tsd_idx, n) {
  
  # 设置每种transposon在此梯度应产生的insertion number -- n_ins
  n_ins = round(0.5*n)
  
  # 对于每种transposon，生成全长插入片段，与对应的tsd拼接起来，并且输出到文件当中
  for (i in 1:length(te_idx)) {
    start_idx = (i-1)*n_ins + 1
    end_idx = i*n_ins
    
    # 随机生成全长片段
    tmp_te_idx = createRandomRegions(nregions = n_ins, length.mean = width(te_idx[i]), length.sd = 0, genome = te_idx[i], non.overlapping = FALSE)
    tmp_te_idx = restrict(tmp_te_idx, start = 1)
    tmp_tsd_idx = tsd_idx[start_idx:end_idx]
    
    # 随机设定TE片段是forward inserted还是reverse inserted
    tmp_tsd_df = as.data.frame(tmp_tsd_idx, row.names = NULL)
    tmp_te_df = as.data.frame(tmp_te_idx, row.names = NULL)
    tmp_te_df$strand <- sample(c("-", "+"), size = nrow(tmp_te_df), replace = TRUE)
    
    # 获取TE片段以及TSD片段，编辑每个insertion的name
    tmp_tsd_seq = as.character(getSeq(genome_fa, tmp_tsd_idx))
    tmp_te_seq = as.character(getSeq(te_fa, makeGRangesFromDataFrame(tmp_te_df)))
    tmp_names = paste(tmp_tsd_df$seqnames, tmp_tsd_df$end-1, tmp_tsd_df$end, tmp_te_df$seqnames, tmp_te_df$strand, sep = "_")
    names(tmp_te_seq) <- tmp_names
    
    # 向TE片段中人为添加一些突变
    origin_seq = getSeq(te_fa, makeGRangesFromDataFrame(tmp_te_df))
    names(origin_seq) = tmp_names
    origin_idx = makeGRangesFromDataFrame(data.frame(chr=origin_seq@ranges@NAMES, start=rep(1,length(origin_seq)), end=width(origin_seq)))
    origin_seq = as.character(origin_seq)
    for (x in 1:length(tmp_te_seq)) {
      if (sample(c(TRUE, FALSE), size = 1, prob = c(0.4, 0.6))) {
        tmp_te_seq[x] = seq_mutate(idx=origin_idx[x], te_seq=tmp_te_seq[x], te_idx, te_fa, genome_idx, genome_fa)
      }
    }
    
    # 将一轮循环种产生的insertion信息输出到文件当中
    tmp_ins_df = data.frame(chrom = tmp_tsd_df$seqnames, start = tmp_tsd_df$end-1, end = tmp_tsd_df$end,
                            type = rep("insertion", n_ins), ins_seq = paste(tmp_te_seq, tmp_tsd_seq, sep = ""),
                            add_len = sample(c(0,1), n_ins, prob=c(0.8, 0.2), replace=TRUE), name = tmp_names, strand = tmp_te_df$strand, TE = tmp_te_df$seqnames,
                            te_start = tmp_te_df$start-1, te_end = tmp_te_df$end, te_seq = tmp_te_seq,
                            tsd_start = tmp_tsd_df$start-1, tsd_len = tmp_tsd_df$end, tsd_seq = tmp_tsd_seq,
                            origin_seq = origin_seq, FullLength = rep("False", n_ins))
    
    write.table(tmp_ins_df, file = "./simulated_sv.summary", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
}




# 从命令行获取参数 ---- te ccs & ref genome
args <- commandArgs(trailingOnly = TRUE)
n = as.integer(args[3]) # 每种transposon生成的insertion总数
t = as.integer(args[4]) # 长度梯度划分数目

# 读取transposon consensus序列，以及genome
te_fa <- open(FaFile(args[1]))
te_idx <- scanFaIndex(te_fa)
te_idx <- narrow(te_idx, end = -2)
genome_fa <- open(FaFile(args[2]))
genome_idx <- scanFaIndex(genome_fa)
genome_idx <- narrow(genome_idx, end = -2)

# 首先生成tsd区间，nregions略多于n*length(te_idx)
tsd_idx <- createRandomRegions(nregions = n*length(te_idx)+50 , length.mean = 8, length.sd = 1, genome = genome_idx)
tsd_idx = restrict(tsd_idx, start = 1)

simulate_ins(te_fa, te_idx, genome_fa, genome_idx, tsd_idx[1:round(0.5*length(tsd_idx))], n, t)
simulate_ins_full(te_fa, te_idx, genome_fa, genome_idx, tsd_idx[(round(0.5*length(tsd_idx))+1):length(tsd_idx)], n)