library(biomaRt)
library(ggplot2)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
getwd()

#E14.5 forebrain lung limb for rep2: expression analysis
#and lungE14.5 vs lungP0

#read annotation file
annotation <- read.csv(file="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/annotation.table.csv",row.names=1,header=TRUE)

#put raw data in data.frame
counts_forebrain_raw <- as.data.frame(read.delim(file="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/mm_E14.5_forebrain_RNAseq_rep2ReadsPerGene.out.edit.tab",header=FALSE,row.names=1))
counts_limb_raw <- as.data.frame(read.delim(file="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/mm_RNA_E14.5_limb_cat_rep2ReadsPerGene.out.edit.tab",header=FALSE,row.names=1))
counts_lung_raw <- as.data.frame(read.delim(file="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/mm_RNA_E14.5_lung_cat_rep2ReadsPerGene.out.edit.tab",header=FALSE,row.names=1))
counts_lung_raw_P0 <- as.data.frame(read.delim(file="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/mm_RNA_P0_lung_cat_rep2ReadsPerGene.out.edit.tab",header=FALSE,row.names=1))
#For replicates 1
counts_forebrain_raw1 <- as.data.frame(read.delim(file="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/mm_E14.5_forebrain_RNAseq_rep1ReadsPerGene.out.edit.tab",header=FALSE,row.names=1))
counts_limb_raw1 <- as.data.frame(read.delim(file="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/mm_RNA_E14.5_limb_cat_rep1ReadsPerGene.out.edit.tab",header=FALSE,row.names=1))
counts_lung_raw1 <- as.data.frame(read.delim(file="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/mm_RNA_E14.5_lung_cat_rep1ReadsPerGene.out.edit.tab",header=FALSE,row.names=1))
counts_lung_raw_P01 <- as.data.frame(read.delim(file="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/mm_RNA_P0_lung_cat_rep1ReadsPerGene.out.edit.tab",header=FALSE,row.names=1))

#real counts per tissue
counts_forebrain <- counts_forebrain_raw[,1]
counts_limb <- counts_limb_raw[,1]
counts_lung <- counts_lung_raw[,1]
counts_lung_P0 <- counts_lung_raw_P0[,1]
#for reps 1
counts_forebrain1 <- counts_forebrain_raw1[,1]
counts_limb1 <- counts_limb_raw1[,1]
counts_lung1 <- counts_lung_raw1[,1]
counts_lung_P01 <- counts_lung_raw_P01[,1]

#get comparisons
forebrain_vs_limb <- as.matrix( cbind( (cbind(counts_forebrain1,counts_forebrain )), (cbind(counts_limb1, counts_limb))))
forebrain_vs_lung <- as.matrix( cbind( (cbind(counts_forebrain1,counts_forebrain )), (cbind(counts_lung1, counts_lung))))
limb_vs_lung <- as.matrix(cbind( (cbind(counts_limb1, counts_limb)) , (cbind(counts_lung1, counts_lung)))) 
lungE14.5_vs_lungP0 <-as.matrix(cbind( (cbind(counts_lung1, counts_lung)) , (cbind(counts_lung_P01, counts_lung_P0))))

#get respective annotations: colnames(counts) == rownames(annotation)
annotation.forebrain_vs_limb <- rbind(rbind(annotation[1,], annotation[2,]) ,rbind(annotation[9,], annotation[10,]))
annotation.forebrain_vs_lung <- rbind(rbind(annotation[1,], annotation[2,]),rbind(annotation[5,], annotation[6,]))
annotation.limb_vs_lung <- rbind(rbind(annotation[9,], annotation[10,]),rbind(annotation[5,], annotation[6,]))
annotation.lungE14.5_vs_lungP0 <- rbind(rbind(annotation[5,], annotation[6,]),rbind(annotation[7,], annotation[8,]))

# define DESeq2 object
dds_forebrain_vs_limb <- DESeqDataSetFromMatrix(countData = forebrain_vs_limb,colData = annotation.forebrain_vs_limb, design = ~1)
dds_forebrain_vs_lung <- DESeqDataSetFromMatrix(countData = forebrain_vs_lung,colData = annotation.forebrain_vs_lung, design = ~ 1)
dds_limb_vs_lung <- DESeqDataSetFromMatrix(countData = limb_vs_lung,colData = annotation.limb_vs_lung, design = ~ 1)
dds_lungE14.5_vs_lungP0 <- DESeqDataSetFromMatrix(countData = lungE14.5_vs_lungP0,colData = annotation.lungE14.5_vs_lungP0, design = ~ 1)

# keep genes with at least 10 reads
keep <- rowSums(counts(dds_forebrain_vs_limb)) >= 10
dds_forebrain_vs_limb <- dds_forebrain_vs_limb[keep,]

keep <- rowSums(counts(dds_forebrain_vs_lung)) >= 10
dds_forebrain_vs_lung <- dds_forebrain_vs_lung[keep,]

keep <- rowSums(counts(dds_limb_vs_lung)) >= 10
dds_limb_vs_lung <- dds_limb_vs_lung[keep,]

keep <- rowSums(counts(dds_lungE14.5_vs_lungP0)) >= 10
dds_lungE14.5_vs_lungP0 <- dds_lungE14.5_vs_lungP0[keep,]

# differentially expressed genes (and other calculations)
dseq_forebrain_vs_limb <- DESeq(dds_forebrain_vs_limb)
res_forebrain_vs_limb <- results(dseq_forebrain_vs_limb)

dseq_forebrain_vs_lung <- DESeq(dds_forebrain_vs_lung)
res_forebrain_vs_lung <- results(dseq_forebrain_vs_lung)

dseq_limb_vs_lung <- DESeq(dds_limb_vs_lung)
res_limb_vs_lung <- results(dseq_limb_vs_lung)

dseq_lungE14.5_vs_lungP0 <- DESeq(dds_lungE14.5_vs_lungP0)
res_lungE14.5_vs_lungP0 <- results(dseq_lungE14.5_vs_lungP0)

#####get differentially expressed genes filtered by  |log2FC| > 0.5 and  padj < 0.01
sum_forebrain_vs_limb <- sum(res_forebrain_vs_limb$padj < 0.01 & res_forebrain_vs_limb$log2FoldChange > 0.5, na.rm=TRUE)
cat("Number of differentially expressed genes between forebrain and limb in E14: ", sum_forebrain_vs_limb)

sum_forebrain_vs_lung <- sum(res_forebrain_vs_lung$padj < 0.01 & res_forebrain_vs_lung$log2FoldChange > 0.5, na.rm=TRUE)
cat("Number of differentially expressed genes between forebrain and lung in E14: ", sum_forebrain_vs_lung)

sum_limb_vs_lung <- sum(res_limb_vs_lung$padj < 0.01 & res_limb_vs_lung$log2FoldChange > 0.5, na.rm=TRUE)
cat("Number of differentially expressed genes between limb and lung in E14: ", sum_limb_vs_lung)

sum_lungE14.5_vs_lungP0 <- sum(res_lungE14.5_vs_lungP0$padj < 0.01 & res_lungE14.5_vs_lungP0$log2FoldChange > 0.5, na.rm=TRUE)
cat("Number of differentially expressed genes between lung in E14.5 and lung in P0: ", sum_lungE14.5_vs_lungP0)



####shrinkage for vizualisation
resultsNames(dds_forebrain_vs_limb)
resLFCn_forebrain_vs_limb <- lfcShrink(dseq_forebrain_vs_limb, coef = resultsNames(dseq_forebrain_vs_limb), type="normal")
summary(resLFCn_forebrain_vs_limb)

# creates a display with three windows
par(mfrow=c(1,3))
#### plot MA-Plot
plotMA(resLFCn_forebrain_vs_limb, ylim=c(-3,3),main="forebrain vs limb (normal)")
par(mfrow=c(1,3))
plotMA(resLFCn_forebrain_vs_lung, ylim=c(-3,3),main="forebrain vs lung (normal)")
par(mfrow=c(1,3))
plotMA(resLFCn_limb_vs_lung, ylim=c(-3,3),main="limb vs lung (normal)")
par(mfrow=c(1,3))
plotMA(resLFCn_lungE14.5_vs_lungP0, ylim=c(-3,3),main="lung in E14.5 vs lung in P0 (normal)")


#volcano plot
####################################

threshold <- abs(res_forebrain_vs_limb$log2FoldChange) > 0.5 & res_forebrain_vs_limb$adj < 0.1
res_forebrain_vs_limb$threshold <- threshold
ggplot(as.data.frame(res_forebrain_vs_limb)) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  ggtitle("forebrain vs limb") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  

####################################


threshold <- abs(res_forebrain_vs_lung$log2FoldChange) > 0.5 & res_forebrain_vs_lung$padj < 0.01
res_forebrain_vs_lung$threshold <- threshold
ggplot(as.data.frame(res_forebrain_vs_lung)) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  ggtitle("forebrain vs lung") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  
##################################

threshold <- abs(res_limb_vs_lung$log2FoldChange) > 0.5 & res_limb_vs_lung$padj < 0.01
res_limb_vs_lung$threshold <- threshold
ggplot(as.data.frame(res_limb_vs_lung)) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  ggtitle("limb vs lung") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  


###################################################

threshold <- res_lungE14.5_vs_lungP0$padj < 0.01 & abs(res_lungE14.5_vs_lungP0$log2FoldChange) > 0.5
res_lungE14.5_vs_lungP0$threshold <- threshold
ggplot(as.data.frame(res_lungE14.5_vs_lungP0)) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  ggtitle("lung in E14.5 vs lung in P0") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  




######.all f�r forebrain-vs-limb
dds_forebrain_vs_limb.all <- DESeqDataSetFromMatrix(countData = forebrain_vs_limb, colData = annotation.forebrain_vs_limb, design = ~1)
keep <- rowSums(counts(dds_forebrain_vs_limb.all)) >=10
dds_forebrain_vs_limb.all <- dds_forebrain_vs_limb.all[keep,]

dds_forebrain_vs_limb.all <- DESeq(dds_forebrain_vs_limb.all)
res.dds_forebrain_vs_limb.all <- results(dds_forebrain_vs_limb.all)
res.dds_forebrain_vs_limb.all
summary(res.dds_forebrain_vs_limb.all)

ntd <- normTransform(dds_forebrain_vs_limb.all)
vsd <- vst(dds_forebrain_vs_limb.all,blind=FALSE)
tr.forebrain_vs_limb.counts <- assay(vsd)
nr.forebrain_vs_limb.counts <- assay(ntd)

select <- order(abs(res.dds_forebrain_vs_limb.all$log2FoldChange), decreasing=TRUE)[1:100]
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, show_colnames=TRUE, fontsize_col=8, filename="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/heatmap1_forebrain_limb.pdf")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, show_colnames=TRUE, fontsize_col=8, filename="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/heatmap2_forebrain_limb.pdf")
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE, show_colnames=TRUE, fontsize_col=8, annotation_col=annotation.forebrain_vs_limb, filename="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/heatmap3_forebrain_limb.pdf")
pheatmap(assay(vsd)[select,], show_rownames=FALSE, fontsize_col=8, annotation_col=annotation.forebrain_vs_limb, cutree_cols=3, cutree_rows=3, main="Heatmap with annotation", filename="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/heatmapCl_forebrain_limb.pdf")


#.all f�r forebrain-vs-lung
dds_forebrain_vs_lung.all <- DESeqDataSetFromMatrix(countData = forebrain_vs_lung, colData = annotation.forebrain_vs_lung, design = ~1)
keep <- rowSums(counts(dds_forebrain_vs_lung.all)) >=10
dds_forebrain_vs_lung.all <- dds_forebrain_vs_lung.all[keep,]

dds_forebrain_vs_lung.all <- DESeq(dds_forebrain_vs_lung.all)
res.dds_forebrain_vs_lung.all <- results(dds_forebrain_vs_lung.all)
res.dds_forebrain_vs_lung.all
summary(res.dds_forebrain_vs_lung.all)

ntd <- normTransform(dds_forebrain_vs_lung.all)
vsd <- vst(dds_forebrain_vs_lung.all,blind=FALSE)
tr.forebrain_vs_lung.counts <- assay(vsd)
nr.forebrain_vs_lung.counts <- assay(ntd)

select <- order(abs(res.dds_forebrain_vs_lung.all$log2FoldChange), decreasing=TRUE)[1:100]
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, show_colnames=TRUE, fontsize_col=8, filename="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/heatmap1_forebrain_lung.pdf")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, show_colnames=TRUE, fontsize_col=8, filename="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/heatmap2_forebrain_lung.pdf")
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE, show_colnames=TRUE, fontsize_col=8, annotation_col=annotation.forebrain_vs_lung, filename="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/heatmap3_forebrain_lung.pdf")
pheatmap(assay(vsd)[select,], show_rownames=FALSE, fontsize_col=8, annotation_col=annotation.forebrain_vs_lung, cutree_cols=3, cutree_rows=3, main="Heatmap with annotation", filename="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/heatmapCl_forebrain_lung.pdf")


#.all f�r limb-vs-lung
dds_limb_vs_lung.all <- DESeqDataSetFromMatrix(countData = limb_vs_lung, colData = annotation.limb_vs_lung, design = ~1)
keep <- rowSums(counts(dds_limb_vs_lung.all)) >=10
dds_limb_vs_lung.all <- dds_limb_vs_lung.all[keep,]

dds_limb_vs_lung.all <- DESeq(dds_limb_vs_lung.all)
res.dds_limb_vs_lung.all <- results(dds_limb_vs_lung.all)
res.dds_limb_vs_lung.all
summary(res.dds_limb_vs_lung.all)

ntd <- normTransform(dds_limb_vs_lung.all)
vsd <- vst(dds_limb_vs_lung.all,blind=FALSE)
tr.limb_vs_lung.counts <- assay(vsd)
nr.limb_vs_lung.counts <- assay(ntd)

select <- order(abs(res.dds_limb_vs_lung.all$log2FoldChange), decreasing=TRUE)[1:100]
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, show_colnames=TRUE, fontsize_col=8, filename="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/heatmap1_limb_vs_lung.pdf")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, show_colnames=TRUE, fontsize_col=8, filename="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/heatmap2_limb_vs_lung.pdf")
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE, show_colnames=TRUE, fontsize_col=8, annotation_col=annotation.limb_vs_lung, filename="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/heatmap3_limb_vs_lung.pdf")
pheatmap(assay(vsd)[select,], show_rownames=FALSE, fontsize_col=8, annotation_col=annotation.limb_vs_lung, cutree_cols=3, cutree_rows=3, main="Heatmap with annotation", filename="C:/Users/feige/Desktop/functionalGenomics/praesentation5/data/heatmapCl_limb_vs_lung.pdf")

#.all f�r lung14-vs-lungP0
dds_lungE14.5_vs_lungP0.all <- DESeqDataSetFromMatrix(countData = lungE14.5_vs_lungP0, colData = annotation.lungE14.5_vs_lungP0, design = ~1)
keep <- rowSums(counts(dds_lungE14.5_vs_lungP0.all)) >=10
dds_lungE14.5_vs_lungP0.all <- dds_lungE14.5_vs_lungP0.all[keep,]

dds_lungE14.5_vs_lungP0.all <- DESeq(dds_lungE14.5_vs_lungP0.all)
res.dds_lungE14.5_vs_lungP0.all <- results(dds_lungE14.5_vs_lungP0.all)
res.dds_lungE14.5_vs_lungP0.all
summary(res.dds_lungE14.5_vs_lungP0.all)

ntd <- normTransform(dds_lungE14.5_vs_lungP0.all)
vsd <- vst(dds_lungE14.5_vs_lungP0.all,blind=FALSE)
tr.lungE14.5_vs_lungP0.counts <- assay(vsd)
nr.lungE14.5_vs_lungP0.counts <- assay(ntd)

select <- order(abs(res.dds_lungE14.5_vs_lungP0.all$log2FoldChange), decreasing=TRUE)[1:100]
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, show_colnames=TRUE, fontsize_col=8, filename="C:/Users/feige/Desktop/functionalGenomics/prasentation5/data/heatmap1_lungE14.5_vs_lungP0.pdf")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, show_colnames=TRUE, fontsize_col=8, filename="C:/Users/feige/Desktop/functionalGenomics/prasentation5/data/heatmap2_lungE14.5_vs_lungP0.pdf")
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE, show_colnames=TRUE, fontsize_col=8, annotation_col=annotation.lungE14.5_vs_lungP0, filename="C:/Users/feige/Desktop/functionalGenomics/prasentation5/data/heatmap3_lungE14.5_vs_lungP0.pdf")
pheatmap(assay(vsd)[select,], show_rownames=FALSE, fontsize_col=8, annotation_col=annotation.lungE14.5_vs_lungP0, cutree_cols=3, cutree_rows=3, main="Heatmap with annotation", filename="C:/Users/feige/Desktop/functionalGenomics/prasentation5/data/heatmapCl_lungE14.5_vs_lungP0.pdf")






#################################
#test
dim(counts_forebrain_raw)
head(counts_forebrain_raw)
head(annotation)

head(counts_forebrain)
head(counts_limb)
head(counts_lung)

head(forebrain_vs_limb)
head(forebrain_vs_lung)
head(limb_vs_lung)

annotation.forebrain_vs_limb 
annotation.forebrain_vs_lung
annotation.limb_vs_lung

dds_forebrain_vs_limb

summary(res_forebrain_vs_limb)

sum_forebrain_vs_limb


