files <- list.files(path = "~/task_bioinf/RNA-SEQ/counts_sra/", pattern = "counts", full.names = TRUE)
df    <- do.call(cbind,lapply(files,function(fn)read.table(fn,header=FALSE, sep="\t")))
df    <- cbind(df)
a <- which(colnames(df)=="NumReads" )
a <- c(1,a)
df <- df[,c(1,2,4,6,8,10,12)]
colnames(df) <- c('ENSEMBL', 'SRR3414629', 'SRR3414630', 'SRR3414631', 'SRR3414635', 'SRR3414636', 'SRR3414637')

samples <- c('SRR3414629', 'SRR3414630', 'SRR3414631', 'SRR3414635', 'SRR3414636', 'SRR3414637')
condition <- c(rep('treatment', 3), rep('control', 3))
metadata <- data.frame(samples, condition)

library(DESeq2)

row.names(df) <- df$ENSEMBL
df[1] <- NULL

ddsMat <- DESeqDataSetFromMatrix(countData = df,
                                 colData = metadata,
                                 design = ~ condition)
dds <- estimateSizeFactors(ddsMat)
dds <- dds[rowMeans(counts(dds, normalized=TRUE)) > 0 ]

dds <- DESeq(dds)
resultsNames(dds) 
res <- results(dds)
resultsNames(dds)
res05 <- results(dds, alpha=0.05)
summary(res05)
resOrdered <- res05[order(res05$pvalue),]
write.csv(resOrdered, '~/task_bioinf/RNA-SEQ/results_deseq2.csv')

