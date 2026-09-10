### Add methylation annotation and calculate the minP across the genes ####

library(data.table)

## 1. Read Annotation
anno <- fread("/restricted/projectnb/afgen/biqiwang/AF_multi_omics/FHS_mQTL/autosome_methy_anno_all.txt",
              select = c("V2", "V4"))
names(anno)<-c("cpgs","Gene")
anno <- anno[!is.na(Gene) & Gene != ""]

setkey(anno, cpgs)


## 2. Read and combine all CpG assoc result files
files <- list.files(
  "/restricted/projectnb/bioinfo_lin/group/biqiwang/ROSMAP_omics_cogtraj/revision/internal_validation/methy_res/",
  pattern = "\\.txt$",
  full.names = TRUE
)


all_cpg <- rbindlist(
  lapply(files, fread),
  use.names = TRUE,
  fill = TRUE
)


## 3. Add gene annotation to all the result file (some CpGs might be removed because they did not annotate to a gene)
all_cpg2 <- anno[
  all_cpg,
  on = "cpgs",
  nomatch = 0
]


## 4. Minimum P by Gene × fold
gene_fold <- all_cpg2[
  !is.na(pval),
  .(
    mth_minP = min(pval),
    mth_minP_cpgs = cpgs[which.min(pval)],
    n_cpgs = uniqueN(cpgs)
  ),
  by = .(Gene, folds)
][
  , adj_minP := 1 - (1 - mth_minP)^n_cpgs
]


setorder(gene_fold, folds, adj_minP)

write.table(gene_fold,"Methy_Genelevel_globalcog_slope_minPval_adjusted_10folds_980_sort.csv",quote=F,col=T,row=F,sep=",")

q('no')
