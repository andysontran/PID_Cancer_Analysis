library(TCGAbiolinks)
library(SummarizedExperiment)
library(recount)
library(TCGAutils)
library(limma)
library(biomaRt)
library(maftools)
library(dplyr)

convert.ENSG.Symbol <- function(genes){
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
  #merge(DEG.ucs,G_list,by.x="gene",by.y="ensembl_gene_id")
  return(G_list)
}

lung.recount.gtex <- TCGAquery_recount2(project = "GTEX", tissue = "lung")
lung.recount.tcga <- TCGAquery_recount2(project = "TCGA", tissue = "lung")

SE.lung.recount.gtex <- lung.recount.gtex$GTEX_blood
SE.lung.recount.tcga <- lung.recount.tcga$TCGA_lymph_nodes
query.dlbc <- GDCquery(project = "TCGA-DLBC",
+                        data.category = "Transcriptome Profiling",
+                        data.type = "Gene Expression Quantification",
+                        workflow.type = "HTSeq - Counts")

samplesDown.dblc <- getResults(query.dlbc, cols = c("cases"))
dataSmTP.dlbc <- TCGAquery_SampleTypes(barcode = samplesDown.dblc,
+                                        typesample = "TP")
dataSmNT.dlbc <- TCGAquery_SampleTypes(barcode = samplesDown.dblc,
+                                        typesample = "NT")

eset.gtex <- assays(scale_counts(bld.recount.gtex$GTEX_blood, round = TRUE))$counts
eset.tcga <- assays(scale_counts(dlbc.recount.tcga$TCGA_lymph_nodes, round = TRUE))$counts
rse_scaled <- scale_counts(bld.recount.gtex$GTEX_blood, round = TRUE)
summary(colSums(assays(rse_scaled)$counts)) / 1e6

colnames(eset.tcga) <- colData(dlbc.recount.tcga$TCGA_lymph_nodes)$gdc_cases.samples.portions.analytes.aliquots.submitter_id

rownames(eset.gtex) <- gsub("\\..*", "", rownames(eset.gtex))
rownames(eset.tcga) <- gsub("\\..*", "", rownames(eset.tcga))

eset.tcga.cancer <- eset.tcga[,which(colData(dlbc.recount.tcga$TCGA_lymph_nodes)$gdc_cases.samples.sample_type=="Primary Tumor")]
eset.tcga.normal <- eset.tcga[,which(colData(dlbc.recount.tcga$TCGA_lymph_nodes)$gdc_cases.samples.sample_type=="Solid Tissue Normal")]

dataPrep.dlbc <- merge(as.data.frame(eset.gtex), as.data.frame(eset.tcga.cancer), by = 0, all = TRUE)
rownames(dataPrep.dlbc) <- dataPrep.dlbc$Row.names
dataPrep.dlbc$Row.names <- NULL

dataNorm.dlbc <- TCGAanalyze_Normalization(tabDF = dataPrep.dlbc,
+                                            geneInfo = geneInfoHT,
+                                            method = "gcContent")

dataFilt.dlbc <- TCGAanalyze_Filtering(tabDF = dataNorm.dlbc,
+                                        method = "quantile",
+                                        qnt.cut = 0.25)

DEG.dlbc <- TCGAanalyze_DEA(mat1 = dataFilt.dlbc[,colnames(eset.gtex)],
+                             mat2 = dataFilt.dlbc[,colnames(eset.tcga.cancer)],
+                             metadata = FALSE,
+                             pipeline = "limma",
+                             voom = TRUE,
+                             Cond1type = "Normal",
+                             Cond2type = "Tumor",
+                             fdr.cut = 0.01,
+                             logFC.cut = 1,
+                             method = "glmLRT")

conversion.table <- convert.ENSG.Symbol(rownames(DEG.dlbc))
conversion.inter.DEG <- intersect(conversion.table[-which(conversion.table$hgnc_symbol==""),]$ensembl_gene_id, rownames(DEG.dlbc))
conversion.table2 <- conversion.table[which(conversion.table$ensembl_gene_id %in% conversion.inter.DEG),]
conversion.table2 <- conversion.table[which(conversion.table$ensembl_gene_id %in% conversion.inter.DEG),]
rownames(conversion.table2) <- conversion.table2$ensembl_gene_id
non-unique value when setting 'row.names': ‘ENSG00000187510’ 
DEG.dlbc <- DEG.dlbc[-which(rownames(DEG.dlbc) %in% "ENSG00000187510"), ]

conversion.table <- convert.ENSG.Symbol(rownames(DEG.dlbc))
conversion.inter.DEG <- intersect(conversion.table[-which(conversion.table$hgnc_symbol==""),]$ensembl_gene_id, rownames(DEG.dlbc))
conversion.table2 <- conversion.table[which(conversion.table$ensembl_gene_id %in% conversion.inter.DEG),]
rownames(conversion.table2) <- conversion.table2$ensembl_gene_id
conversion.table[-which(conversion.table$hgnc_symbol==""),]
  
DEG.dlbc.hgnc <- DEG.dlbc[conversion.inter.DEG,]
DEGs.dlbc.hgnc <- merge(DEG.dlbc.hgnc, conversion.table2, by = 0)
rownames(DEGs.dlbc.hgnc) <- DEGs.dlbc.hgnc$Row.names
DEGs.dlbc.hgnc$Row.names <- NULL
rownames(DEGs.dlbc.hgnc) <- DEGs.dlbc.hgnc$hgnc_symbol
write.csv(DEGs.dlbc.hgnc, "DEGsDLBC.csv")

TCGAVisualize_volcano(DEGs.dlbc.hgnc$logFC, DEGs.dlbc.hgnc$adj.P.Val,
+                       filename = "DLBC_+5.pdf", xlab = "logFC",
+                       names = DEGs.dlbc.hgnc$hgnc_symbol, show.names = "highlighted",
+                       x.cut = 1, y.cut = 0.01,
+                       highlight = DEGs.dlbc.hgnc$hgnc_symbol[which(DEGs.dlbc.hgnc$logFC >= 5)][1:50],
+                       highlight.color = "orange")
TCGAVisualize_volcano(DEGs.dlbc.hgnc$logFC, DEGs.dlbc.hgnc$adj.P.Val,
+                       filename = "DLBC_-5.pdf", xlab = "logFC",
+                       names = DEGs.dlbc.hgnc$hgnc_symbol, show.names = "highlighted",
+                       x.cut = 1, y.cut = 0.01,
+                       highlight = DEGs.dlbc.hgnc$hgnc_symbol[which(DEGs.dlbc.hgnc$logFC <= -5)][1:50],
+                       highlight.color = "blue")


PID.DEGs.dlbc.hgnc <- subset(DEGs.dlbc.hgnc, hgnc_symbol %in% Panel_PIDGenes$X1)

TCGAVisualize_volcano(PID.DEGs.dlbc.hgnc$logFC, PID.DEGs.dlbc.hgnc$adj.P.Val,
+                       filename = "PID.DLBC_+2.pdf", xlab = "logFC",
+                       names = PID.DEGs.dlbc.hgnc$hgnc_symbol, show.names = "highlighted",
+                       x.cut = 1, y.cut = 0.01,
+                       highlight = PID.DEGs.dlbc.hgnc$hgnc_symbol[which(PID.DEGs.dlbc.hgnc$logFC >= 2)][1:50],
+                       highlight.color = "orange")
TCGAVisualize_volcano(PID.DEGs.dlbc.hgnc$logFC, PID.DEGs.dlbc.hgnc$adj.P.Val,
+                       filename = "PID.DLBC_-2.pdf", xlab = "logFC",
+                       names = PID.DEGs.dlbc.hgnc$hgnc_symbol, show.names = "highlighted",
+                       x.cut = 1, y.cut = 0.01,
+                       highlight = PID.DEGs.dlbc.hgnc$hgnc_symbol[which(PID.DEGs.dlbc.hgnc$logFC <= -2)][1:50],
+                       highlight.color = "blue")

EA.Genelist.dlbc <- rownames(DEGs.dlbc.hgnc)
system.time(ansEA.dlbc <- TCGAanalyze_EAcomplete(TFname = "DEA genes Normal (GTEX-Blood) Vs Tumour (TCGA-DLBC)", EA.Genelist.dlbc))
TCGAvisualize_EAbarplot(tf = rownames(ansEA.dlbc$ResBP),
+                         GOBPTab = ansEA.dlbc$ResBP,
+                         GOCCTab = ansEA.dlbc$ResCC,
+                         GOMFTab = ansEA.dlbc$ResMF,
+                         PathTab = ansEA.dlbc$ResPat,
+                         nRGTab = EA.Genelist.dlbc,
+                         nBar = 10,
+                         filename = "ALL.DLBC.EA.pdf")
        
PID.EAGenelist.dlbc <- rownames(PID.DEGs.dlbc.hgnc)
system.time(PID.ansEA.dlbc <- TCGAanalyze_EAcomplete(TFname = "DEA PID genes Normal (GTEX-Blood) Vs Tumour (TCGA-DLBC)", PID.EAGenelist.dlbc))
TCGAvisualize_EAbarplot(tf = rownames(PID.ansEA.dlbc$ResBP),
+                         GOBPTab = PID.ansEA.dlbc$ResBP,
+                         GOCCTab = PID.ansEA.dlbc$ResCC,
+                         GOMFTab = PID.ansEA.dlbc$ResMF,
+                         PathTab = PID.ansEA.dlbc$ResPat,
+                         nRGTab = PID.EAGenelist.dlbc,
+                         nBar = 10,
+                         filename = "PID.DLBC.EA.pdf")
