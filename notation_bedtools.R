library(ChIPseeker)
library(org.Hs.eg.db)
library(GenomicRanges)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

ICM_special_peaks <- readPeakFile("../pval_12_data/ICM_special_peaks_pval_12.bed")
TE_special_peaks <- readPeakFile('../pval_12_data/TE_special_peaks_pval_12.bed')


# ICM
peak_list <- list(ICM_special_peaks=ICM_special_peaks, TE_special_peaks=TE_special_peaks)
# covplot(ICM_ATAC_special_peaks, weightCol="V5")

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- lapply(peak_list, getTagMatrix, windows=promoter)

plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region",
            ylab = "Read Count Frequency", facet="row")

# Annotation 1
ICM_special_peaks_Anno <- annotatePeak(ICM_special_peaks, tssRegion=c(-3000, 3000),
                                 TxDb=txdb, annoDb="org.Hs.eg.db", verbose=FALSE)
plotAnnoPie(ICM_special_peaks_Anno)
write.table(as.data.frame(ICM_special_peaks_Anno), "../pval_12_data/ICM_peaks_Anno.xls", quote=F, row.names=F, sep="\t")

# 2
TE_special_peaks_Anno <- annotatePeak(TE_special_peaks, tssRegion=c(-3000, 3000),
                                  TxDb=txdb, annoDb="org.Hs.eg.db", verbose=FALSE)
plotAnnoPie(TE_special_peaks_Anno)
write.table(as.data.frame(TE_special_peaks_Anno), "../pval_12_data/TE_peaks_Anno.xls", quote=F, row.names=F, sep="\t")

# KEGG
anno_list <- lapply(peak_list, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE, annoDb="org.Hs.eg.db")
genes = lapply(anno_list, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster=genes, fun="enrichKEGG", pvalueCutoff=0.05, pAdjustMethod="BH")
dotplot(compKEGG, showCategory=15, title="KEGG Pathway Enrichment Analysis")

genes = lapply(anno_list, function(i) as.data.frame(i)$geneId)
vennplot(genes)
