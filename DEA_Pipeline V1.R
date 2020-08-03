library(TCGAbiolinks)
library(SummarizedExperiment)
library(recount)
library(TCGAutils)
library(limma)
library(biomaRt)
library(maftools)
library(dplyr)
library(enhancedvolcano)

convert.ENSG.Symbol <- function(genes){
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
  #merge(DEG.ucs,G_list,by.x="gene",by.y="ensembl_gene_id")
  return(G_list)
}

!!!.recount.gtex <- TCGAquery_recount2(project = "GTEX", tissue = "###")
!!!.recount.tcga <- TCGAquery_recount2(project = "TCGA", tissue = "#####")

SE.!!!.recount.gtex <- __.recount.gtex$GTEX_###
SE.!!!.recount.tcga <- __.recount.tcga$TCGA_#####
query.!!! <- GDCquery(project = "TCGA-!!!",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts")

samplesDown.!!! <- getResults(query.!!!, cols = c("cases"))
dataSmTP.!!! <- TCGAquery_SampleTypes(barcode = samplesDown.!!!,
                                       typesample = "TP")   ### Or "TB" if blood-derived cancer ###
dataSmNT.!!! <- TCGAquery_SampleTypes(barcode = samplesDown.!!!,
                                       typesample = "NT")
                                       
!!!.eset.gtex <- assays(scale_counts(!!!.recount.gtex$GTEX_###, round = TRUE))$counts
!!!.eset.tcga <- assays(scale_counts(!!!.recount.tcga$TCGA_#####, round = TRUE))$counts

rse_scaled <- scale_counts(!!!.recount.gtex$GTEX_###, round = TRUE)
summary(colSums(assays(rse_scaled)$counts)) / 1e6

colnames(!!!.eset.tcga) <- colData(!!!.recount.tcga$TCGA_#####)$gdc_cases.samples.portions.analytes.aliquots.submitter_id

rownames(!!!.eset.gtex) <- gsub("\\..*", "", rownames(!!!.eset.gtex))
rownames(!!!.eset.tcga) <- gsub("\\..*", "", rownames(!!!.eset.tcga))

!!!.eset.tcga.cancer <- !!!.eset.tcga[,which(colData(!!!.recount.tcga$TCGA_#####)$gdc_cases.samples.sample_type=="Primary Tumor")]
!!!.eset.tcga.normal <- !!!.eset.tcga[,which(colData(!!!.recount.tcga$TCGA_#####)$gdc_cases.samples.sample_type=="Solid Tissue Normal")]

dataPrep.!!! <- merge(as.data.frame(!!!.eset.gtex), as.data.frame(!!!.eset.tcga.cancer), by = 0, all = TRUE)
rownames(dataPrep.!!!) <- dataPrep.!!!$Row.names
dataPrep.!!!$Row.names <- NULL

dataNorm.!!! <- TCGAanalyze_Normalization(tabDF = dataPrep.!!!,
                                           geneInfo = geneInfoHT,
                                           method = "gcContent")

dataFilt.!!! <- TCGAanalyze_Filtering(tabDF = dataNorm.!!!,
                                       method = "quantile",
                                       qnt.cut = 0.25)

DEG.!!! <- TCGAanalyze_DEA(mat1 = dataFilt.!!![,colnames(!!!.eset.gtex)],
                           mat2 = dataFilt.!!![,colnames(!!!.eset.tcga.cancer)],
                           metadata = FALSE,
                           pipeline = "limma",
                           voom = TRUE,
                           Cond1type = "Normal",
                           Cond2type = "Tumor",
                           method = "glmLRT")

!!!.conversion.table <- convert.ENSG.Symbol(rownames(DEG.!!!))
!!!.conversion.inter.DEG <- intersect(!!!.conversion.table[-which(!!!.conversion.table$hgnc_symbol==""),]$ensembl_gene_id, rownames(DEG.!!!))
!!!.conversion.table2 <- !!!.conversion.table[which(!!!.conversion.table$ensembl_gene_id %in% !!!.conversion.inter.DEG),]
rownames(!!!.conversion.table2) <- !!!.conversion.table2$ensembl_gene_id
!!!.conversion.table[-which(!!!.conversion.table$hgnc_symbol==""),]

## If recive error: "non-unique value when setting 'row.names': ‘ENSGXXXXXXXXXX’" then, and repeat above ##
DEG.!!! <- DEG.!!![-which(rownames(DEG.!!!) %in% c("ENSGXXXXXXXXXX")), ]

DEG.!!!.hgnc <- DEG.!!![!!!.conversion.inter.DEG,]
DEGs.!!!.hgnc <- merge(DEG.!!!.hgnc, !!!.conversion.table2, by = 0)
rownames(DEGs.!!!.hgnc) <- DEGs.!!!.hgnc$Row.names
DEGs.!!!.hgnc$Row.names <- NULL
rownames(DEGs.!!!.hgnc) <- DEGs.!!!.hgnc$hgnc_symbol

PID.DEGs.!!!.hgnc <- subset(DEGs.!!!.hgnc, hgnc_symbol %in% Panel_PIDGenes$X1)

write.csv(PID.DEGs.!!!.hgnc, "PID.DEGs!!!.csv")

## Plot DEA after filtering for PID-related genes
EnhancedVolcano(PID.DEGs.!!!.hgnc,
                lab = rownames(PID.DEGs.!!!.hgnc),
                x = 'logFC',
                y = 'adj.P.Val',
                ylab = bquote(~-Log[10]~ "(FDR corrected-P values)")
                title = "?????? (TCGA-!!!) vs. Normal (GTEX-###)",
                subtitle = "Differential expression using limma-voom; showing PID-related genes",
                caption = "FC cutoff, 2; p-value cutoff, 10e-16",
                pCutoff = 10e-16,
                FCcutoff = 2,
                pointSize = 5.0,
                labSize = 5.0,
                legendPosition = 'none')               
