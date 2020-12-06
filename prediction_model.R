library(biomaRt)
library(ggplot2)
library(DESeq2)
library(Rsubread)
library(RColorBrewer)
library(pheatmap)
library(stringr)
#library(tidyverse)
library(roperators)
library(GenomicRanges)
library(Rsamtools)
library(bamsignals)
library(GenomicFeatures)
library(TriMatch)
library(reshape2)
library(dplyr)
library(Momocs)

############
#Test-Area:########################################################################################################################################
############
#...auskommeniert
#Dies sollte automatisch die Files in enem Ordner erkennen und Variablennamen dafür deklariere nund die Funktionen ausführen.
#lediglich die letzte For-Loop funktioniert hier noch nicht ganz, sollte aber vom Prinzip ger das richtige tun.. 
if(FASLE){

#save count_basenames
files <- list.files(path="C:/Users/feige/Desktop/functionalGenomics/praesentation6/data/", pattern="*.bam$", full.names=TRUE, recursive=FALSE)
count_files <- vector()
for( i in files){
   count_files <- c(count_files, paste0("counts_",basename(i)))
}

######apply bamCount() to all files, saved in corresponding variable with name= count_basename(x)
lapply(files, function(x) {
   assign(paste0("counts_",basename(x)), bamCount(x, proms, verbose=FALSE), envir=globalenv())
})

######normalization:
#usage: works only with a complete set of data
#exception: forebrain-files have to be treated seperatly, 'caused of differet naming convention.
   # --> rename the forebrain_files! ContChIP-Data should be 7th in order per block (just like our own datasets are ordered)
j=2 #here it is 2, for local testing ! change j to 7 for server run 
for(i in length(count_files)){     
   j = j-1
   if(i %% 2 == 0){
      j = 2 #change it also here to 7 for a server run
   }else{
      assign(paste0("m_",count_files[i]), median(((get(count_files[i]))+1)/((get(count_files[(i+j)]))+1)))
      assign(paste0("Snorm_",count_files[i]), (get(count_files[i])+1)/(get(count_files[(i+j)])+1) * (1/m), envir=globalenv())
      get(tail(ls(),n=1)) <- log(get(tail(ls(),n=1)))
      assign(paste0("Snormscaled_",count_files[i]), scale(get(tail(ls(),n=1)),center = 1, scale = TRUE), envir=globalenv())
   }
}

} #Kommentarblock Ende !
###################################################################################################################################################


#p0 Lung
bampathLungP0Cont <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ContChIP_P0_lung_rep1.filt.nodup.srt.bam"#
bampathLungP0_H3K27me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_P0_lung_H3K27me3_rep1.filt.nodup.srt.bam"
bampathLungP0_H3K4me1  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_P0_lung_H3K4me1_rep1.filt.nodup.srt.bam"
bampathLungP0_H3K4me3  <-"C:/Users/lacky/Desktop/Softwareprojek/tmm_ChIP_P0_lung_H3K4me3_rep1.filt.nodup.srt.bam"
bampathLungP0_H3K9me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_P0_lung_H3K9me3_rep1.filt.nodup.srt.bam"
bampathLungP0_H3k27ac  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_P0_lung_H3k27ac_rep1.filt.nodup.srt.bam"
bampathLungP0_H3k36me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_P0_lung_H3k36me3_rep1.filt.nodup.srt.bam"

#P0 Limb
bampathlimbP0Cont <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ContChIP_P0_limb_rep1.filt.nodup.srt.bam"#
bampathlimbP0_H3K27me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_P0_limb_H3K27me3_rep1.filt.nodup.srt.bam"
bampathlimbP0_H3K4me1  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_P0_limb_H3K4me1_rep1.filt.nodup.srt.bam"
bampathlimbP0_H3K4me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_P0_limb_H3K4me3_rep1.filt.nodup.srt.bam"
bampathlimbP0_H3K9me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_P0_limb_H3K9me3_rep1.filt.nodup.srt.bam"
bampathlimbP0_H3k27ac  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_P0_limb_H3k27ac_rep1.filt.nodup.srt.bam"
bampathlimbP0_H3k36me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_P0_limb_H3k36me3_rep1.filt.nodup.srt.bam"

#P0 forebrain
bampathforebrainP0Cont <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_P0_forebrain_Control_rep1.bam"#
bampathforebrainP0_H3K27me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_P0_forebrain_H3K27me3_rep1.bam"
bampathforebrainP0_H3K4me1  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_P0_forebrain_H3K4me1_rep1.bam"
bampathforebrainP0_H3K4me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_P0_forebrain_H3K4me3_rep1.bam"
bampathforebrainP0_H3K9me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_P0_forebrain_H3K9me3_rep1.bam"
bampathforebrainP0_H3k27ac  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_P0_forebrain_H3K27ac_rep1.bam"
bampathforebrainP0_H3k36me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_P0_forebrain_H3k36me3_rep1.bam"

#E14 Lung
bampathLungE14.5Cont <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ContChIP_E14.5_lung_rep1.filt.nodup.srt.bam"#
bampathLungE14.5_H3K27me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_E14.5_lung_H3K27me3_rep1.filt.nodup.srt.bam"
bampathLungE14.5_H3K4me1  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_E14.5_lung_H3K4me1_rep1.filt.nodup.srt.bam"
bampathLungE14.5_H3K4me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_E14.5_lung_H3K4me3_rep1.filt.nodup.srt.bam"
bampathLungE14.5_H3K9me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_E14.5_lung_H3K9me3_rep1.filt.nodup.srt.bam"
bampathLungE14.5_H3k27ac  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_E14.5_lung_H3k27ac_rep1.filt.nodup.srt.bam"
bampathLungE14.5_H3k36me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_E14.5_lung_H3k36me3_rep1.filt.nodup.srt.bam"

#E14 Limb
bampathlimbE14.5Cont <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ContChIP_E14.5_limb_rep1.filt.nodup.srt.bam"#
bampathlimbE14.5_H3K27me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_E14.5_limb_H3K27me3_rep1.filt.nodup.srt.bam"
bampathlimbE14.5_H3K4me1  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_E14.5_limb_H3K4me1_rep1.filt.nodup.srt.bam"
bampathlimbE14.5_H3K4me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_E14.5_limb_H3K4me3_rep1.filt.nodup.srt.bam"
bampathlimbE14.5_H3K9me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_E14.5_limb_H3K9me3_rep1.filt.nodup.srt.bam"
bampathlimbE14.5_H3k27ac  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_E14.5_limb_H3k27ac_rep1.filt.nodup.srt.bam"
bampathlimbE14.5_H3k36me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_ChIP_E14.5_limb_H3k36me3_rep1.filt.nodup.srt.bam"

#E14 forebrain
bampathforebrainE14.5Cont <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_E14.5_forebrain_Control_rep1.bam"#
bampathforebrainE14.5_H3K27me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_E14.5_forebrain_H3K27me3_rep1.bam"
bampathforebrainE14.5_H3K4me1  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_E14.5_forebrain_H3K4me1_rep1.bam"
bampathforebrainE14.5_H3K4me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_E14.5_forebrain_H3K4me3_rep1.bam"
bampathforebrainE14.5_H3K9me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_E14.5_forebrain_H3K9me3_rep1.bam"
bampathforebrainE14.5_H3k27ac  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_E14.5_forebrain_H3K27ac_rep1.bam"
bampathforebrainE14.5_H3k36me3  <-"C:/Users/lacky/Desktop/Softwareprojekt/mm_E14.5_forebrain_H3k36me3_rep1.bam"



ensembl.mouse <- useMart(biomart="ENSEMBL_MART_ENSEMBL",  dataset="mmusculus_gene_ensembl") 
prot.gene <- getBM(attributes=c("ensembl_gene_id","gene_biotype","strand","chromosome_name","start_position","end_position"), mart = ensembl.mouse, filters = "biotype", values = "protein_coding")
prot.gene <- prot.gene[prot.gene$chromosome_name %in% 1:19, ] #filter autosomal chromosomes (thx Schoepflin for "%in%")
prot.gene$strand <- gsub(pattern='-1', replacement='-', prot.gene$strand) # Zur anpassung der Strand-Notation
prot.gene$strand <- gsub(pattern='1', replacement='+', prot.gene$strand)
#make GRanges-Object:
#online-tutorial sagt, es braucht nur seqnames, ranges und strand -Informationen und kann die werte aus getBM-file-output$whatever verwenden
grange <- GRanges(seqnames = Rle (paste("chr", prot.gene$chromosome_name, sep="")), ranges = IRanges(start=prot.gene$start_position, end=prot.gene$end_position, names=prot.gene$ensembl_gene_id), strand = Rle(strand(prot.gene$strand)))
proms <- GenomicRanges::promoters(grange, upstream=2000, downstream=500)
head(proms)

#P0 Counts Lung
# countsLungP0Cont      <- bamCount(bampathLungP0Cont, proms, verbose=FALSE)
# countsLungP0_H3K27me3 <- bamCount(bampathLungP0_H3K27me3, proms, verbose=FALSE)
# countsLungP0_H3K4me1  <- bamCount(bampathLungP0_H3K4me1, proms, verbose=FALSE)
# countsLungP0_H3K4me3  <- bamCount(bampathLungP0_H3K4me3, proms, verbose=FALSE)
# countsLungP0_H3K9me3  <- bamCount(bampathLungP0_H3K9me3, proms, verbose=FALSE)
# countsLungP0_H3k27ac  <- bamCount(bampathLungP0_H3k27ac, proms, verbose=FALSE)
# countsLungP0_H3k36me3 <- bamCount(bampathLungP0_H3k36me3, proms, verbose=FALSE)
# 
# #P0 counts limb
# countslimbP0Cont      <- bamCount(bampathlimbP0Cont, proms, verbose=FALSE)
# countslimbP0_H3K27me3 <- bamCount(bampathlimbP0_H3K27me3, proms, verbose=FALSE)
# countslimbP0_H3K4me1  <- bamCount(bampathlimbP0_H3K4me1, proms, verbose=FALSE)
# countslimbP0_H3K4me3  <- bamCount(bampathlimbP0_H3K4me3, proms, verbose=FALSE)
# countslimbP0_H3K9me3  <- bamCount(bampathlimbP0_H3K9me3, proms, verbose=FALSE)
# countslimbP0_H3k27ac  <- bamCount(bampathlimbP0_H3k27ac, proms, verbose=FALSE)
# countslimbP0_H3k36me3 <- bamCount(bampathlimbP0_H3k36me3, proms, verbose=FALSE)
# 
# #P0 counts forebrain
# countsforebrainP0Cont      <- bamCount(bampathforebrainP0Cont, proms, verbose=FALSE)
# countsforebrainP0_H3K27me3 <- bamCount(bampathforebrainP0_H3K27me3, proms, verbose=FALSE)
# countsforebrainP0_H3K4me1  <- bamCount(bampathforebrainP0_H3K4me1, proms, verbose=FALSE)
# countsforebrainP0_H3K4me3  <- bamCount(bampathforebrainP0_H3K4me3, proms, verbose=FALSE)
# countsforebrainP0_H3K9me3  <- bamCount(bampathforebrainP0_H3K9me3, proms, verbose=FALSE)
# countsforebrainP0_H3k27ac  <- bamCount(bampathforebrainP0_H3k27ac, proms, verbose=FALSE)
# countsforebrainP0_H3k36me3 <- bamCount(bampathforebrainP0_H3k36me3, proms, verbose=FALSE)

#E14 counts lung
countsLungE14.5Cont      <- bamCount(bampathLungE14.5Cont, proms, verbose=FALSE)
countsLungE14.5_H3K27me3 <- bamCount(bampathLungE14.5_H3K27me3, proms, verbose=FALSE)
countsLungE14.5_H3K4me1  <- bamCount(bampathLungE14.5_H3K4me1, proms, verbose=FALSE)
countsLungE14.5_H3K4me3  <- bamCount(bampathLungE14.5_H3K4me3, proms, verbose=FALSE)
countsLungE14.5_H3K9me3  <- bamCount(bampathLungE14.5_H3K9me3, proms, verbose=FALSE)
countsLungE14.5_H3k27ac  <- bamCount(bampathLungE14.5_H3k27ac, proms, verbose=FALSE)
countsLungE14.5_H3k36me3 <- bamCount(bampathLungE14.5_H3k36me3, proms, verbose=FALSE)

names(countsLungE14.5Cont)<- names(ranges(proms))
names(countsLungE14.5_H3K27me3)<- names(ranges(proms))
names(countsLungE14.5_H3K4me1)<- names(ranges(proms))
names(countsLungE14.5_H3K4me3)<- names(ranges(proms))
names(countsLungE14.5_H3K9me3)<- names(ranges(proms))
names(countsLungE14.5_H3k27ac)<- names(ranges(proms))
names(countsLungE14.5_H3k36me3)<- names(ranges(proms))

#E14 counts limb
countslimbE14.5Cont      <- bamCount(bampathlimbE14.5Cont, proms, verbose=FALSE)
countslimbE14.5_H3K27me3 <- bamCount(bampathlimbE14.5_H3K27me3, proms, verbose=FALSE)
countslimbE14.5_H3K4me1  <- bamCount(bampathlimbE14.5_H3K4me1, proms, verbose=FALSE)
countslimbE14.5_H3K4me3  <- bamCount(bampathlimbE14.5_H3K4me3, proms, verbose=FALSE)
countslimbE14.5_H3K9me3  <- bamCount(bampathlimbE14.5_H3K9me3, proms, verbose=FALSE)
countslimbE14.5_H3k27ac  <- bamCount(bampathlimbE14.5_H3k27ac, proms, verbose=FALSE)
countslimbE14.5_H3k36me3 <- bamCount(bampathlimbE14.5_H3k36me3, proms, verbose=FALSE)

# #E14 counts forebrain
# countsforebrainE14.5Cont      <- bamCount(bampathforebrainE14.5Cont, proms, verbose=FALSE)
# countsforebrainE14.5_H3K27me3 <- bamCount(bampathforebrainE14.5_H3K27me3, proms, verbose=FALSE)
# countsforebrainE14.5_H3K4me1  <- bamCount(bampathforebrainE14.5_H3K4me1, proms, verbose=FALSE)
# countsforebrainE14.5_H3K4me3  <- bamCount(bampathforebrainE14.5_H3K4me3, proms, verbose=FALSE)
# countsforebrainE14.5_H3K9me3  <- bamCount(bampathforebrainE14.5_H3K9me3, proms, verbose=FALSE)
# countsforebrainE14.5_H3k27ac  <- bamCount(bampathforebrainE14.5_H3k27ac, proms, verbose=FALSE)
# countsforebrainE14.5_H3k36me3 <- bamCount(bampathforebrainE14.5_H3k36me3, proms, verbose=FALSE)
# 
# 
# head(countsLungE14.5_H3K27me3)
# countsLungP0_H3K27me3
#########################################################################Normalisierung################################################# 
######################################################################################################Lunge1

#H3K27me3
m=median((countsLungE14.5_H3K27me3+1)/(countsLungE14.5Cont+1))
SnormLungE14.5_H3K27me3= (countsLungE14.5_H3K27me3+1)/(countsLungE14.5Cont+1) * 1/m
#SnormLungE14.5_H3K27me3
SnormLungE14.5_H3K27me3 <- log(SnormLungE14.5_H3K27me3)
#SnormLungE14.5_H3K27me3
SnormLungE14.5_H3K27me3_scaled <- scale(SnormLungE14.5_H3K27me3,center = 1, scale = TRUE)
#SnormLungE14.5_H3K27me3_scaled

#H3K4me1
m=median((countsLungE14.5_H3K4me1+1)/(countsLungE14.5Cont+1))
SnormLungE14.5_H3K4me1= (countsLungE14.5_H3K4me1+1)/(countsLungE14.5Cont+1) * 1/m
#SnormLungE14.5_H3K4me1
SnormLungE14.5_H3K4me1 <- log(SnormLungE14.5_H3K4me1)
#SnormLungE14.5_H3K4me1
SnormLungE14.5_H3K4me1_scaled <- scale(SnormLungE14.5_H3K4me1,center = 1, scale = TRUE)
#SnormLungE14.5_H3K4me1_scaled

#H3K4me3
m=median((countsLungE14.5_H3K4me3+1)/(countsLungE14.5Cont+1))
SnormLungE14.5_H3K4me3= (countsLungE14.5_H3K4me3+1)/(countsLungE14.5Cont+1) * 1/m
#SnormLungE14.5_H3K4me3
SnormLungE14.5_H3K4me3 <- log(SnormLungE14.5_H3K4me3)
#SnormLungE14.5_H3K4me3
SnormLungE14.5_H3K4me3_scaled <- scale(SnormLungE14.5_H3K4me3,center = 1, scale = TRUE)
#SnormLungE14.5_H3K4me3_scaled

#H3K9me3
m=median((countsLungE14.5_H3K9me3+1)/(countsLungE14.5Cont+1))
SnormLungE14.5_H3K9me3= (countsLungE14.5_H3K9me3+1)/(countsLungE14.5Cont+1) * 1/m
#SnormLungE14.5_H3K9me3
SnormLungE14.5_H3K9me3 <- log(SnormLungE14.5_H3K9me3)
#SnormLungE14.5_H3K9me3
SnormLungE14.5_H3K9me3_scaled <- scale(SnormLungE14.5_H3K9me3,center = 1, scale = TRUE)
#SnormLungE14.5_H3K9me3_scaled

#H3k27ac
m=median((countsLungE14.5_H3k27ac+1)/(countsLungE14.5Cont+1))
SnormLungE14.5_H3k27ac= (countsLungE14.5_H3k27ac+1)/(countsLungE14.5Cont+1) * 1/m
#SnormLungE14.5_H3k27ac
SnormLungE14.5_H3k27ac <- log(SnormLungE14.5_H3k27ac)
#SnormLungE14.5_H3k27ac
SnormLungE14.5_H3k27ac_scaled <- scale(SnormLungE14.5_H3k27ac,center = 1, scale = TRUE)
#SnormLungE14.5_H3k27ac_scaled

#H3k36me3
m=median((countsLungE14.5_H3k36me3+1)/(countsLungE14.5Cont+1))
SnormLungE14.5_H3k36me3= (countsLungE14.5_H3k36me3+1)/(countsLungE14.5Cont+1) * 1/m
#SnormLungE14.5_H3k36me3
SnormLungE14.5_H3k36me3 <- log(SnormLungE14.5_H3k36me3)
#SnormLungE14.5_H3k36me3
SnormLungE14.5_H3k36me3_scaled <- scale(SnormLungE14.5_H3k36me3,center = 1, scale = TRUE)
#SnormLungE14.5_H3k36me3_scaled

rownames(SnormLungE14.5_H3K27me3_scaled)<- names(ranges(proms))
rownames(SnormLungE14.5_H3K4me1_scaled)<- names(ranges(proms))
rownames(SnormLungE14.5_H3K4me3_scaled)<- names(ranges(proms))
rownames(SnormLungE14.5_H3K9me3_scaled)<- names(ranges(proms))
rownames(SnormLungE14.5_H3k27ac_scaled)<- names(ranges(proms))
rownames(SnormLungE14.5_H3k36me3_scaled)<- names(ranges(proms))

#head(SnormLungE14.5_H3K27me3_scaled)

#################################################################################################Limb
#H3K27me3
m=median((countslimbE14.5_H3K27me3+1)/(countslimbE14.5Cont+1))
SnormlimbE14.5_H3K27me3= (countslimbE14.5_H3K27me3+1)/(countslimbE14.5Cont+1) * 1/m
#SnormlimbE14.5_H3K27me3
SnormlimbE14.5_H3K27me3 <- log(SnormlimbE14.5_H3K27me3)
#SnormlimbE14.5_H3K27me3
SnormlimbE14.5_H3K27me3_scaled <- scale(SnormlimbE14.5_H3K27me3,center = 1, scale = TRUE)
#SnormlimbE14.5_H3K27me3_scaled

#H3K4me1
m=median((countslimbE14.5_H3K4me1+1)/(countslimbE14.5Cont+1))
SnormlimbE14.5_H3K4me1= (countslimbE14.5_H3K4me1+1)/(countslimbE14.5Cont+1) * 1/m
#SnormlimbE14.5_H3K4me1
SnormlimbE14.5_H3K4me1 <- log(SnormlimbE14.5_H3K4me1)
#SnormlimbE14.5_H3K4me1
SnormlimbE14.5_H3K4me1_scaled <- scale(SnormlimbE14.5_H3K4me1,center = 1, scale = TRUE)
#SnormlimbE14.5_H3K4me1_scaled

#H3K4me3
m=median((countslimbE14.5_H3K4me3+1)/(countslimbE14.5Cont+1))
SnormlimbE14.5_H3K4me3= (countslimbE14.5_H3K4me3+1)/(countslimbE14.5Cont+1) * 1/m
#SnormlimbE14.5_H3K4me3
SnormlimbE14.5_H3K4me3 <- log(SnormlimbE14.5_H3K4me3)
#SnormlimbE14.5_H3K4me3
SnormlimbE14.5_H3K4me3_scaled <- scale(SnormlimbE14.5_H3K4me3,center = 1, scale = TRUE)
#SnormlimbE14.5_H3K4me3_scaled

#H3K9me3
m=median((countslimbE14.5_H3K9me3+1)/(countslimbE14.5Cont+1))
SnormlimbE14.5_H3K9me3= (countslimbE14.5_H3K9me3+1)/(countslimbE14.5Cont+1) * 1/m
#SnormlimbE14.5_H3K9me3
SnormlimbE14.5_H3K9me3 <- log(SnormlimbE14.5_H3K9me3)
#SnormlimbE14.5_H3K9me3
SnormlimbE14.5_H3K9me3_scaled <- scale(SnormlimbE14.5_H3K9me3,center = 1, scale = TRUE)
#SnormlimbE14.5_H3K9me3_scaled

#H3k27ac
m=median((countslimbE14.5_H3k27ac+1)/(countslimbE14.5Cont+1))
SnormlimbE14.5_H3k27ac= (countslimbE14.5_H3k27ac+1)/(countslimbE14.5Cont+1) * 1/m
#SnormlimbE14.5_H3k27ac
SnormlimbE14.5_H3k27ac <- log(SnormlimbE14.5_H3k27ac)
#SnormlimbE14.5_H3k27ac
SnormlimbE14.5_H3k27ac_scaled <- scale(SnormlimbE14.5_H3k27ac,center = 1, scale = TRUE)
#SnormlimbE14.5_H3k27ac_scaled

#H3k36me3
m=median((countslimbE14.5_H3k36me3+1)/(countslimbE14.5Cont+1))
SnormlimbE14.5_H3k36me3= (countslimbE14.5_H3k36me3+1)/(countslimbE14.5Cont+1) * 1/m
#SnormlimbE14.5_H3k36me3
SnormlimbE14.5_H3k36me3 <- log(SnormlimbE14.5_H3k36me3)
#SnormlimbE14.5_H3k36me3
SnormlimbE14.5_H3k36me3_scaled <- scale(SnormlimbE14.5_H3k36me3,center = 1, scale = TRUE)
#SnormlimbE14.5_H3k36me3_scaled

rownames(SnormlimbE14.5_H3K27me3_scaled)<- names(ranges(proms))
rownames(SnormlimbE14.5_H3K4me1_scaled)<- names(ranges(proms))
rownames(SnormlimbE14.5_H3K4me3_scaled)<- names(ranges(proms))
rownames(SnormlimbE14.5_H3k27ac_scaled)<- names(ranges(proms))
rownames(SnormlimbE14.5_H3k36me3_scaled)<- names(ranges(proms))
rownames(SnormLungE14.5_H3K9me3_scaled)<- names(ranges(proms))

# ########################################################################################################Forebrain
# #H3K27me3
# m=median((countsforebrainE14.5_H3K27me3+1)/(countsforebrainE14.5Cont+1))
# SnormforebrainE14.5_H3K27me3= (countsforebrainE14.5_H3K27me3+1)/(countsforebrainE14.5Cont+1) * 1/m
# #SnormforebrainE14.5_H3K27me3
# SnormforebrainE14.5_H3K27me3 <- log(SnormforebrainE14.5_H3K27me3)
# #SnormforebrainE14.5_H3K27me3
# SnormforebrainE14.5_H3K27me3_scaled <- scale(SnormforebrainE14.5_H3K27me3,center = 1, scale = TRUE)
# #SnormforebrainE14.5_H3K27me3_scaled
# 
# #H3K4me1
# m=median((countsforebrainE14.5_H3K4me1+1)/(countsforebrainE14.5Cont+1))
# SnormforebrainE14.5_H3K4me1= (countsforebrainE14.5_H3K4me1+1)/(countsforebrainE14.5Cont+1) * 1/m
# #SnormforebrainE14.5_H3K4me1
# SnormforebrainE14.5_H3K4me1 <- log(SnormforebrainE14.5_H3K4me1)
# #SnormforebrainE14.5_H3K4me1
# SnormforebrainE14.5_H3K4me1_scaled <- scale(SnormforebrainE14.5_H3K4me1,center = 1, scale = TRUE)
# #SnormforebrainE14.5_H3K4me1_scaled
# 
# #H3K4me3
# m=median((countsforebrainE14.5_H3K4me3+1)/(countsforebrainE14.5Cont+1))
# SnormforebrainE14.5_H3K4me3= (countsforebrainE14.5_H3K4me3+1)/(countsforebrainE14.5Cont+1) * 1/m
# #SnormforebrainE14.5_H3K4me3
# SnormforebrainE14.5_H3K4me3 <- log(SnormforebrainE14.5_H3K4me3)
# #SnormforebrainE14.5_H3K4me3
# SnormforebrainE14.5_H3K4me3_scaled <- scale(SnormforebrainE14.5_H3K4me3,center = 1, scale = TRUE)
# #SnormforebrainE14.5_H3K4me3_scaled
# 
# #H3K9me3
# m=median((countsforebrainE14.5_H3K9me3+1)/(countsforebrainE14.5Cont+1))
# SnormforebrainE14.5_H3K9me3= (countsforebrainE14.5_H3K9me3+1)/(countsforebrainE14.5Cont+1) * 1/m
# #SnormforebrainE14.5_H3K9me3
# SnormforebrainE14.5_H3K9me3 <- log(SnormforebrainE14.5_H3K9me3)
# #SnormforebrainE14.5_H3K9me3
# SnormforebrainE14.5_H3K9me3_scaled <- scale(SnormforebrainE14.5_H3K9me3,center = 1, scale = TRUE)
# #SnormforebrainE14.5_H3K9me3_scaled
# 
# #H3k27ac
# m=median((countsforebrainE14.5_H3k27ac+1)/(countsforebrainE14.5Cont+1))
# SnormforebrainE14.5_H3k27ac= (countsforebrainE14.5_H3k27ac+1)/(countsforebrainE14.5Cont+1) * 1/m
# #SnormforebrainE14.5_H3k27ac
# SnormforebrainE14.5_H3k27ac <- log(SnormforebrainE14.5_H3k27ac)
# #SnormforebrainE14.5_H3k27ac
# SnormforebrainE14.5_H3k27ac_scaled <- scale(SnormforebrainE14.5_H3k27ac,center = 1, scale = TRUE)
# #SnormforebrainE14.5_H3k27ac_scaled
# 
# #H3k36me3
# m=median((countsforebrainE14.5_H3k36me3+1)/(countsforebrainE14.5Cont+1))
# SnormforebrainE14.5_H3k36me3= (countsforebrainE14.5_H3k36me3+1)/(countsforebrainE14.5Cont+1) * 1/m
# #SnormforebrainE14.5_H3k36me3
# SnormforebrainE14.5_H3k36me3 <- log(SnormforebrainE14.5_H3k36me3)
# #SnormforebrainE14.5_H3k36me3
# SnormforebrainE14.5_H3k36me3_scaled <- scale(SnormforebrainE14.5_H3k36me3,center = 1, scale = TRUE)
# #SnormforebrainE14.5_H3k36me3_scaled
# 
# rownames(SnormforebrainE14.5_H3K27me3_scaled)<- names(ranges(proms))
# rownames(SnormforebrainE14.5_H3K4me1_scaled)<- names(ranges(proms))
# rownames(SnormforebrainE14.5_H3K4me3_scaled)<- names(ranges(proms))
# rownames(SnormforebrainE14.5_H3k27ac_scaled)<- names(ranges(proms))
# rownames(SnormforebrainE14.5_H3k36me3_scaled)<- names(ranges(proms))
# rownames(SnormLungE14.5_H3K9me3_scaled)<- names(ranges(proms))

################################################################################Matrix bilden########################################################

LungE14Matrix <- as.matrix(cbind(cbind(cbind(cbind(SnormLungE14.5_H3K27me3_scaled,SnormLungE14.5_H3K4me1_scaled), cbind(SnormLungE14.5_H3K4me3_scaled,SnormLungE14.5_H3K9me3_scaled)),SnormLungE14.5_H3k27ac_scaled),SnormLungE14.5_H3k36me3_scaled ))
colnames(LungE14Matrix) <- c("H3K27me3","H3K4me1","H3K4me3","H3K9me3","H3k27ac","H3k36me3")
#rownames(LungE14Matrix) <- rownames(SnormforebrainE14.5_H3K27me3_scaled)
head(LungE14Matrix)

limbE14Matrix <- as.matrix(cbind(cbind(cbind(cbind(SnormlimbE14.5_H3K27me3_scaled,SnormlimbE14.5_H3K4me1_scaled), cbind(SnormlimbE14.5_H3K4me3_scaled,SnormlimbE14.5_H3K9me3_scaled)),SnormlimbE14.5_H3k27ac_scaled),SnormlimbE14.5_H3k36me3_scaled ))
colnames(limbE14Matrix) <- c("H3K27me3","H3K4me1","H3K4me3","H3K9me3","H3k27ac","H3k36me3")

# forebrainE14Matrix <- as.matrix(cbind(cbind(cbind(cbind(SnormforebrainE14.5_H3K27me3_scaled,SnormforebrainE14.5_H3K4me1_scaled), cbind(SnormforebrainE14.5_H3K4me3_scaled,SnormforebrainE14.5_H3K9me3_scaled)),SnormforebrainE14.5_H3k27ac_scaled),SnormforebrainE14.5_H3k36me3_scaled ))
# colnames(forebrainE14Matrix) <- c("H3K27me3","H3K4me1","H3K4me3","H3K9me3","H3k27ac","H3k36me3")


#########################################################################EXONLENGTH BERECHNUNG ########################################################

txdb<-makeTxDbFromGFF(file="C:/Users/lacky/Desktop/Softwareprojekt/gencode.vM24.annotation.gtf",format="gtf")
# Return a list of exons grouped by gene
exons <-exonsBy(txdb,by='gene')
exon_length<-lapply(exons, function(e) { 
               flattened <-reduce(e)
               s2 <-sum(width(flattened))
               return(s2) })
#exon_length
###################################################Splicing der GenIds

#Testlist speichert die Namen der exons zwischen nur zur sicherheit da exon_length sehr lange berechnet hat
testlist <- exon_length
desired_length <- length(testlist) 
empty_list <- vector(mode = "list", length = desired_length) # empty_list wird die beschnittenen Exonnamen entgegennehmen

for(i in (1:length(testlist))){
  empty_list[i] <- unlist(strsplit(names(testlist[i]),"\\."))[1]
}
names(exon_length) <- empty_list
head(exon_length)


#########################################################################RNA-Seq load#############################################
rnaLung1 <- read.table(file="C:/Users/lacky/Desktop/Softwareprojekt/NeueRNA/mm_RNA_E14.5_lung_cat_rep1ReadsPerGene.out.edit.tab",header=FALSE,sep ='\t', row.names = 1 )
#rnaLung2 <- read.table(file="C:/Users/lacky/Desktop/Softwareprojekt/NeueRNA/mm_RNA_E14.5_lung_cat_rep2ReadsPerGene.out.edit.tab",header=FALSE,sep ='\t', row.names = 1 )
rnaLimb1 <- read.table(file="C:/Users/lacky/Desktop/Softwareprojekt/NeueRNA/mm_RNA_E14.5_limb_cat_rep1ReadsPerGene.out.edit.tab",header=FALSE,sep ='\t', row.names = 1 )
#rnaLimb2 <- read.table(file="C:/Users/lacky/Desktop/Softwareprojekt/NeueRNA/mm_RNA_E14.5_limb_cat_rep2ReadsPerGene.out.edit.tab",header=FALSE,sep ='\t', row.names = 1 )

#L�schen der nicht ben�tigten Spalten
rnaLung1[,4] <- NULL
rnaLung1[,3] <- NULL
rnaLung1[,2] <- NULL

# rnaLung2[,4] <- NULL
# rnaLung2[,3] <- NULL
# rnaLung2[,2] <- NULL

rnaLimb1[,4] <- NULL
rnaLimb1[,3] <- NULL
rnaLimb1[,2] <- NULL

# rnaLimb2[,4] <- NULL
# rnaLimb2[,3] <- NULL
# rnaLimb2[,2] <- NULL

########################Aufbereiten der Liste an benannten Vectoren Exon_length. Wird so zum DataFrame um cbind nutzen zu k�nnen

tmp <- exon_length
exon_length_df <- melt(tmp)
rownames(exon_length_df) <- exon_length_df[,2]
exon_length_df[,2] <- NULL
write.table(exon_length_df, file = "exon_length_df.txt",sep = "\t",row.names= TRUE, col.names = TRUE)
########################################################## exonlegth merge#######################################################FOREBRAIN FEHLT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
head(rnaLung1)
head(exon_length_df)

###### Sortieren von rnaLung1 (Anhand der Ids)
ord.rnaLung1<-cbind(rownames(rnaLung1)[order(rownames(rnaLung1))], rnaLung1[order(rownames(rnaLung1)),])
#Da ich eine Spalte zu viel hatte nachdem sortieren habe ich die gew�nschten spalten in einen neuen Dataframe gepackt.
ord2.rnaLung1 <- data.frame(rna = ord.rnaLung1[,2],exon_length=exon_length_df[,1])
rownames(ord2.rnaLung1) <- ord.rnaLung1[,1]
head(ord2.rnaLung1) 

# ####RNALung2
# ord.rnaLung2<-cbind(rownames(rnaLung2)[order(rownames(rnaLung2))], rnaLung2[order(rownames(rnaLung2)),])
# #Da ich eine Spalte zu viel hatte nachdem sortieren habe ich die gew�nschten spalten in einen neuen Dataframe gepackt.
# ord2.rnaLung2 <- data.frame(rna = ord.rnaLung2[,2],exon_length=exon_length_df[,1])
# rownames(ord2.rnaLung2) <- ord.rnaLung2[,1]
# head(ord2.rnaLung2)

##RnaLimb1
ord.rnaLimb1<-cbind(rownames(rnaLimb1)[order(rownames(rnaLimb1))], rnaLimb1[order(rownames(rnaLimb1)),])
#Da ich eine Spalte zu viel hatte nachdem sortieren habe ich die gew�nschten spalten in einen neuen Dataframe gepackt.
ord2.rnaLimb1 <- data.frame(rna = ord.rnaLimb1[,2],exon_length=exon_length_df[,1])
rownames(ord2.rnaLimb1) <- ord.rnaLimb1[,1]
head(ord2.rnaLimb1)

# #### RnaLimb2
# ord.rnaLimb2<-cbind(rownames(rnaLimb2)[order(rownames(rnaLimb2))], rnaLimb2[order(rownames(rnaLimb2)),])
# #Da ich eine Spalte zu viel hatte nachdem sortieren habe ich die gew�nschten spalten in einen neuen Dataframe gepackt.
# ord2.rnaLimb2 <- data.frame(rna = ord.rnaLimb2[,2],exon_length=exon_length_df[,1])
# rownames(ord2.rnaLimb2) <- ord.rnaLimb2[,1]
# head(ord2.rnaLimb2)

########################### RPKM BERECHNUNG!!!###########################################################

#LUNG1
#Diese Umwandlung in Numerische Vectoren ist leider n�tig weil eine meiner Funktionen die ich oben genutzt habe die zahlen werte in Factors umwandelt, auf denen man keine Operationen ausf�hren kann
rnaVectorLung1 <- ord2.rnaLung1[,1]
rnaVectorLung1 <- as.numeric((as.character(rnaVectorLung1)))

exonVectorLung1 <- ord2.rnaLung1[,2]
exonVectorLung1 <- as.numeric((as.character(exonVectorLung1)))

RPKMLung1 <- rnaVectorLung1 / (exonVectorLung1/1000 * (sum(rnaVectorLung1))/1000000)
head(RPKMLung1)

rownames(ord2.rnaLung1)
names(RPKMLung1) <- rownames(ord2.rnaLung1)
head(RPKMLung1)
rpkmMatrixLunge1 <- as.matrix(RPKMLung1,row.names = names(RPKMLung1))############################################################################
colnames(rpkmMatrixLunge1) <- c("values")
#write.table(rpkmMatrixLunge1, file = "rpkmMatrixLunge1.txt",sep = "\t",row.names= TRUE, col.names = TRUE)

# #Lung2
# rnaVectorLung2 <- ord2.rnaLung2[,1]
# rnaVectorLung2 <- as.numeric((as.character(rnaVectorLung2)))
# 
# exonVectorLung2 <- ord2.rnaLung2[,2]
# exonVectorLung2 <- as.numeric((as.character(exonVectorLung2)))
# 
# RPKMLung2 <- rnaVectorLung2 / (exonVectorLung2/1000 * (sum(rnaVectorLung2))/1000000)
# head(RPKMLung2)
# 
# rpkmMatrixLunge2 <- as.matrix(RPKMLung2,row.names = names(RPKMLung2))############################################################################
# colnames(rpkmMatrixLunge2) <- c("values")

#Limb1
rnaVectorLimb1 <- ord2.rnaLimb1[,1]
rnaVectorLimb1 <- as.numeric((as.character(rnaVectorLimb1)))

exonVectorLimb1 <- ord2.rnaLimb1[,2]
exonVectorLimb1 <- as.numeric((as.character(exonVectorLimb1)))

RPKMLimb1 <- rnaVectorLimb1 / (exonVectorLimb1/1000 * (sum(rnaVectorLimb1))/1000000)
head(RPKMLimb1)

names(RPKMLimb1) <- rownames(ord2.rnaLimb1)
rpkmMatrixlimb1 <- as.matrix(RPKMLimb1,row.names = names(RPKMLimb1))############################################################################
colnames(rpkmMatrixlimb1) <- c("values")

# #Limb2
# rnaVectorLimb2 <- ord2.rnaLimb2[,1]
# rnaVectorLimb2 <- as.numeric((as.character(rnaVectorLimb2)))
# 
# exonVectorLimb2 <- ord2.rnaLimb2[,2]
# exonVectorLimb2 <- as.numeric((as.character(exonVectorLimb2)))
# 
# RPKMLimb2 <- rnaVectorLimb2 / (exonVectorLimb2/1000 * (sum(rnaVectorLimb2))/1000000)
# head(RPKMLimb2)
# rpkmMatrixlimb2 <- as.matrix(RPKMlimb2,row.names = names(RPKMlimb2))############################################################################
# colnames(rpkmMatrixlimb2) <- c("values")


#########################
#RPKM FOREBRAIN FEHLT!!!!!!!!!!!!!! -------> DADURCH AUCH RPKMFOREBRAIN UND DIE RPKMMATRIX DAZU!!!!!

###############################################################################SPLITTEN DER DATENSETS
LungMatrix <- cbind(LungE14Matrix, rpkmMatrixLunge1[, "values"][match(rownames(LungE14Matrix), rownames(rpkmMatrixLunge1))])
colnames(LungMatrix) <- c("H3K27me3","H3K4me1","H3K4me3","H3K9me3","H3k27ac","H3k36me3","RPKM")
head(LungMatrix)
LimbMatrix  <- cbind(limbE14Matrix, rpkmMatrixlimb1[, "values"][match(rownames(limbE14Matrix), rownames(rpkmMatrixlimb1))])
colnames(LimbMatrix) <- c("H3K27me3","H3K4me1","H3K4me3","H3K9me3","H3k27ac","H3k36me3","RPKM")
#forebrainMatrix  <- cbind(forebrainE14Matrix, rpkmMatrixforebrain1[, "values"][match(rownames(forebrainE14Matrix), rownames(rpkmMatrixforebrain1))])
nrow(LungMatrix)


###################################### Lm() berechnung Lung #################################################
#trainLung <- as.data.frame(LungMatrix[1:ceiling(((nrow(LungMatrix)/100)*75)),])
#testLung <- as.data.frame(LungMatrix[ceiling(((nrow(LungMatrix)/100)*25)):nrow(LungMatrix),])

set.seed(101) # Set Seed so that same sample can be reproduced in future also
# Now Selecting 75% of data as sample from total 'n' rows of the data
sample <- sample.int(n = nrow(LungMatrix), size = floor(.75*nrow(LungMatrix)), replace = F)
trainLung <- as.data.frame(LungMatrix[sample, ])
testLung  <- as.data.frame(LungMatrix[-sample,])
sample

fmLung <- lm( log(RPKM+1) ~ . ,data = trainLung)

coefficients(fmLung)
summary(fmLung)

fittedLungTrain <- fitted.values(fmLung)
fittedLungTest <- predict(fmLung,testLung)

summary(fittedLungTest)


RPKMLUNGTMP <- trainLung[,"RPKM"]
CorrelationLungTrain <- cor(log(trainLung[,"RPKM"]+1),fittedLungTrain)
CorrelationLungTest <- cor(log(testLung[,"RPKM"]+1),fittedLungTest)

CorrelationLungTrain
CorrelationLungTest

RSSLungTrain <- sum((RPKMLUNGTMP - fittedLungTrain)^2)
TSSLungTrain <- sum((RPKMLUNGTMP - mean(RPKMLUNGTMP))^2)
RSquareLungTrain <- 1 - (RSSLungTrain/TSSLungTrain)
RSquareLungTrain

RSSLungTest <- sum((testLung[,"RPKM"] - fittedLungTest)^2)
TSSLungTest <- sum((testLung[,"RPKM"] - mean(testLung[,"RPKM"]))^2)
RSquareLungTest <- 1 - (RSSLungTest/TSSLungTest)
RSquareLungTest

smoothScatter(fittedLungTest,log(testLung[,"RPKM"]+1),  xlab = "Predicted Values" , ylab = "True Values",nbin = 500)
###################################### LM() LIMB###########################################################################


#trainLimb <- as.data.frame(LimbMatrix[1:ceiling(((nrow(LimbMatrix)/100)*75)),])
#testLimb <- as.data.frame(LimbMatrix[ceiling(((nrow(LimbMatrix)/100)*25)):nrow(LimbMatrix),])

set.seed(102) # Set Seed so that same sample can be reproduced in future also
# Now Selecting 75% of data as sample from total 'n' rows of the data
sample <- sample.int(n = nrow(LimbMatrix), size = floor(.75*nrow(LimbMatrix)), replace = F)
trainLimb <- as.data.frame(LimbMatrix[sample, ])
testLimb  <- as.data.frame(LimbMatrix[-sample,])


fmLimb <- lm( log(RPKM+1) ~ . ,data = trainLimb)

coefficients(fmLimb)
summary(fmLimb)

fittedLimbTrain <- fitted.values(fmLimb)
fittedLimbTest <- predict(fmLimb,testLimb)
summary(fittedLimbTest)


RPKMLimbTMP <- trainLimb[,"RPKM"]
CorrelationLimbTrain <- cor(log(trainLimb[,"RPKM"]+1),fittedLimbTrain)
CorrelationLimbTest <- cor(log(testLimb[,"RPKM"]+1),fittedLimbTest)

CorrelationLimbTrain
CorrelationLimbTest

RSSLimbTrain <- sum((RPKMLimbTMP - fittedLimbTrain)^2)
TSSLimbTrain <- sum((RPKMLimbTMP - mean(RPKMLimbTMP))^2)
RSquareLimbTrain <- 1 - (RSSLimbTrain/TSSLimbTrain)
RSquareLimbTrain

RSSLimbTest <- sum((testLimb[,"RPKM"] - fittedLimbTest)^2)
TSSLimbTest <- sum((testLimb[,"RPKM"] - mean(testLimb[,"RPKM"]))^2)
RSquareLimbTest <- 1 - (RSSLimbTest/TSSLimbTest)
RSquareLimbTest

smoothScatter(fittedLimbTest,log(testLimb[,"RPKM"]+1) , xlab = "Predicted Values" , ylab = "True Values",nbin = 500)

############################################################### Experimentelles Validieren auf Fremdtissues

fremdPredict <- predict(fmLung,trainLimb)
cor(fremdPredict,log(trainLimb[,"RPKM"]+1) )
summary(fremdPredict)
######How to compare???
RSSFremd <-  sum((RPKMLimbTMP - fremdPredict)^2)
TSSFremd <- sum((RPKMLimbTMP - mean(RPKMLimbTMP))^2)
RSquareFremd <- 1 - (RSSLimbTrain/TSSLimbTrain)
RSquareFremd
