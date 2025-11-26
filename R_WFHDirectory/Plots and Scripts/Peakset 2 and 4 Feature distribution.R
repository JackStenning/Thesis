#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
#BiocManager::install("ChIPseeker")
#BiocManager::install("EnsDb.Hsapiens.v86")
#BiocManager::install("clusterProfiler")
#BiocManager::install("AnnotationDbi")
#BiocManager::install("org.Hs.eg.db")


# Load libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)

#Set working directory
setwd("//userfs/jps558/w2k/Desktop/R_WorkingDirectory") 

#Load data
samplefiles <- list.files("peaks/", pattern='.bed', full.names = T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("Peakset2", "Peakset4")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#Get Annotations
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE)
peakAnnoList

#Plot data barchart
plotAnnoBar(peakAnnoList, title = "Feature distribution of ESR1-HyPB binding loci in peakset 2 and 4")

#Distance to TSS
plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci \n relative to TSS")

#Annotations for each peakset
Set2_annot <- data.frame(peakAnnoList[["Peakset2"]]@anno)
Set4_annot <- data.frame(peakAnnoList[["Peakset4"]]@anno)

#Map EntrezID symbols
#Get IDs
entrez2 <- Set2_annot$geneId
entrez4 <- Set4_annot$geneId

#Return gene symbol

annotations_ebd2 <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                          keys = entrez2,
                                          columns = c("GENENAME"),
                                          keytype = "ENTREZID")
annotations_ebd4 <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                          keys = entrez4,
                                          columns = c("GENENAME"),
                                          keytype = "ENTREZID")


#Change IDs to Character type to merge
annotations_ebd2$ENTREZID <- as.character(annotations_ebd2$ENTREZID)


annotations_ebd4$ENTREZID <- as.character(annotations_ebd4$ENTREZID)

#Write to file
Set2_annot %>%
  left_join(annotations_ebd2, by=c("geneId"="ENTREZID"), relationship = "many-to-many") %>% 
  write.table(file="output/Set2_PeakAnnotation.txt", sep="\t", quote=F, row.names=F)
Set4_annot %>%
  left_join(annotations_ebd4, by=c("geneId"="ENTREZID"), relationship = "many-to-many") %>% 
  write.table(file="output/Set4_PeakAnnotation.txt", sep="\t", quote=F, row.names=F)

#Overrepresentation of Gene ontology terms
ego2 <- enrichGO(gene=entrez2,
                 keyType = "ENTREZID",
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 readable = TRUE)

ego4 <- enrichGO(gene=entrez4,
                 keyType = "ENTREZID",
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 readable = TRUE)

#Output results to table
cluster_summary2 <- data.frame(ego2)
write.csv(cluster_summary2, "output/ClusterProfiler_set2.csv")

cluster_summary4 <- data.frame(ego4)
write.csv(cluster_summary4, "output/ClusterProfiler_set4.csv")

#Dot Plot visualisation
dotplot(ego4, showCategory=50)


#KEGG
ekegg4 <- enrichKEGG(gene = entrez4,
                     organism = 'hsa',
                     qvalueCutoff = 0.05)
dotplot(ekegg4)
