### Add RNAseq annotation and calculate the minP across the genes ####

library(data.table)

## 1. Read Annotation
anno <- fread("/restricted/projectnb/bioinfo_lin/data/ROSMAP/ROSMAP_RNAseq_FPKM_gene_overlapping_with_Gencode_v26.txt",
              select = c("Chr","ID", "GID","Gene_type"))
## restricted to autosome and protein coding genes
anno<-anno[!(anno$Chr %in% c("X","Y")),]
anno<-anno[anno$Gene_type=="protein_coding",]
anno$Gene<-anno$GID

anno$Gene[anno$ID=="ENSG00000108387.10"]<-"SEPT4"
anno$Gene[anno$ID=="ENSG00000140623.9"]<-"SEPTIN12" 
anno$Gene[anno$ID=="ENSG00000144583.4"]<-"MARCHF4"
anno$Gene[anno$ID=="ENSG00000145416.9"]<-"MARCH1"
anno$Gene[anno$ID=="ENSG00000154997.8"]<-"SEPTIN14"
anno$Gene[anno$ID=="ENSG00000164402.9"]<-"SEPTIN8"
anno$Gene[anno$ID=="ENSG00000165406.9"]<-"MARCH8"
anno$Gene[anno$ID=="ENSG00000173838.7"]<-"MARCHF10"
anno$Gene[anno$ID=="ENSG00000173926.5"]<-"MARCHF3"
anno$Gene[anno$ID=="ENSG00000180096.7"]<-"SEPT1"
anno$Gene[anno$ID=="ENSG00000183291.11"]<-"SELENOF"
anno$Gene[anno$ID=="ENSG00000183654.7"]<-"MARCHF11"
anno$Gene[anno$ID=="ENSG00000186522.10"]<-"SEPTIN10"
anno$Gene[anno$ID=="ENSG00000099785.5"]<-"MARCH2"
anno$Gene[anno$ID=="ENSG00000100167.15"]<-"SEPTIN3"
anno$Gene[anno$ID=="ENSG00000117791.11"]<-"MARC2"
anno$Gene[anno$ID=="ENSG00000122545.12"]<-"SEPTIN7"
anno$Gene[anno$ID=="ENSG00000136536.9"]<-"MARCH7"
anno$Gene[anno$ID=="ENSG00000138758.7"]<-"SEPTIN11"
anno$Gene[anno$ID=="ENSG00000139266.5"]<-"MARCH9"
anno$Gene[anno$ID=="ENSG00000145495.9"]<-"MARCH6"
anno$Gene[anno$ID=="ENSG00000168385.13"]<-"SEPT2"
anno$Gene[anno$ID=="ENSG00000173077.10"]<-"DEC1"
anno$Gene[anno$ID=="ENSG00000184640.12"]<-"SEPT9"
anno$Gene[anno$ID=="ENSG00000184702.13"]<-"SEPTIN5"
anno$Gene[anno$ID=="ENSG00000186205.8"]<-"MARC1"
anno$Gene[anno$ID=="ENSG00000198060.4"]<-"MARCHF5"

anno<-anno[,c("ID","Gene")]
anno <- anno[!is.na(Gene) & Gene != ""]

setkey(anno, ID)


## 2. Read and combine all CpG assoc result files
files <- list.files(
  "/restricted/projectnb/bioinfo_lin/group/biqiwang/ROSMAP_omics_cogtraj/revision/internal_validation/expression_res/",
  pattern = "\\.txt$",
  full.names = TRUE
)


all_rna <- rbindlist(
  lapply(files, fread),
  use.names = TRUE,
  fill = TRUE
)

names(all_rna)[4]<-"ID"

## 3. Add gene annotation to all the result file (some CpGs might be removed because they did not annotate to a gene)
all_rna2 <- anno[
  all_rna,
  on = "ID",
  nomatch = 0
]


## 4. Minimum P by Gene × fold
gene_fold <- all_rna2[
  !is.na(pval),
  .(
    exp_minP = min(pval),
    exp_minP_rna = ID[which.min(pval)],
    n_rna = uniqueN(ID)
  ),
  by = .(Gene, folds)
][
  , adj_minP := 1 - (1 - exp_minP)^n_rna
]


setorder(gene_fold, folds, adj_minP)

write.table(gene_fold,"RNAseq_Genelevel_globalcog_slope_minPval_adjusted_10folds_980_sort.csv",quote=F,col=T,row=F,sep=",")

q('no')
