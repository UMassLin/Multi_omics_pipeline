### Add proteomics annotation and calculate the minP across the genes ####

library(data.table)


## Read proteomics results
  globres <- fread(
    paste0(
      "/restricted/projectnb/bioinfo_lin/group/biqiwang/",
      "ROSMAP_omics_cogtraj/revision/internal_validation/proteomic_res/",
      "Brain_proteomics_on_globalcog_slope_10folds.txt"
    ),
    data.table = FALSE
  )

  # Extract gene symbol
  globres$gene <- sub("\\|.*", "", globres$prot_id)
  
  # Remove unwanted genes
  globres2 <- globres[
    globres$gene != "0" &
    globres$gene != "uncharacterized",
  ]



## Find Minimum P by Gene × fold
globres2<-data.table(globres2)

gene_fold <- globres2[
  !is.na(pval),
  .(
    prot_minP = min(pval),
    prot_minP_prot = prot_id[which.min(pval)],
    n_prot = uniqueN(prot_id)
  ),
  by = .(gene, folds)
][
  , adj_minP := 1 - (1 - prot_minP)^n_prot
]


setorder(gene_fold, folds, adj_minP)

write.table(gene_fold,"Proteom_Genelevel_globalcog_slope_minPval_adjusted_10folds_980_sort.csv",quote=F,col=T,row=F,sep=",")

q('no')
