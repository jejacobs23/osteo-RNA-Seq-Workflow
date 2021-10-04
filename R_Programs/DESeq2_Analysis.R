#Take the RNA count data and meta data and run a differential gene expression analysis using DESeq2.
#
#Based off of: https://lashlock.github.io/compbio/R_presentation.html

library(DESeq2)
library(ggplot2)

setwd("C:/Users/Jamie/Work/Projects & Publications/Immunodeficiency and Cancer Project/Data/RNA-Seq")

#Load the RNA-Seq raw count data
countData <- read.table("countData.txt", header = TRUE)
head(countData)

#Load the meta data
metaData <- read.table("metaData_U.txt", header = TRUE)
metaData

#Create the DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData=countData,
                              colData=metaData,
                              design=~Group, tidy = TRUE)

dds

dds <- DESeq(dds)

res <- results(dds)
head(results(dds, tidy=TRUE))
summary(res) #summary of results

#Sort the results by the adjusted p-value (fdr)
res <- res[order(res$padj),]
head(res)
outfile_name="DiffExpressGenes_U.txt"
write.table(res, file = outfile_name, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

#Make a volcano plot
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
