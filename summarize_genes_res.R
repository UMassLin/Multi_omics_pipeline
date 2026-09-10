### Summarize Gene-level results ####

library(data.table)
#library(VennDiagram)

#outcomes <- c(
#  "globalcog",
#  "Episodic_memory",
#  "Working_memory",
#  "Semantic_memory",
#  "Perceptual_speed",
#  "Visuospatial_ability"
#)

library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75


#for(outcome in outcomes){

  meth<-fread(paste0("/restricted/projectnb/bioinfo_lin/group/biqiwang/ROSMAP_omics_cogtraj/revision/internal_validation/",
  "Methy_Genelevel_globalcog_slope_minPval_adjusted_10folds_980_sort.csv"),data.table=F)
  
  rna<-fread(paste0("/restricted/projectnb/bioinfo_lin/group/biqiwang/ROSMAP_omics_cogtraj/revision/internal_validation/",
  "RNAseq_Genelevel_globalcog_slope_minPval_adjusted_10folds_980_sort.csv"),data.table=F)
  
  prot<-fread(paste0("/restricted/projectnb/bioinfo_lin/group/biqiwang/ROSMAP_omics_cogtraj/revision/internal_validation/",
  "Proteom_Genelevel_globalcog_slope_minPval_adjusted_10folds_980_sort.csv"),data.table=F)

names(prot)[1]<-"Gene"

### Merge the three omics data
## Convert to data.table
setDT(meth)
setDT(rna)
setDT(prot)

## Rename P-value columns before merging
setnames(meth, "adj_minP", "methylation_P")
setnames(rna,  "adj_minP", "expression_P")
setnames(prot,  "adj_minP", "proteomics_P")

## Merge by Gene + fold
multiomics <- Reduce(
  function(x, y)
    merge(
      x, y,
      by = c("Gene", "folds"),
      all = FALSE,
      sort = FALSE
    ),
  list(meth, rna, prot)
)


## fixed weights
den <- 645 + 567 + 539

w1 <- 645 / den
w2 <- 567 / den
w3 <- 539 / den

## denominator of weighted Z
zden <- sqrt(w1^2 + w2^2 + w3^2)

## calculate directly
multiomics[
  ,
  `:=`(
    z = (
      w1 * qnorm(methylation_P) +
      w2 * qnorm(expression_P) +
      w3 * qnorm(proteomics_P)
    ) / zden
  )
]

multiomics[, p := pnorm(z)]

## fast in-place sorting
setorder(multiomics, folds, p)

for (f in unique(multiomics$folds)) {
  
  fwrite(
    multiomics[folds == f],
    file = paste0("multiomics_fold", f, ".txt"),
    sep = "\t"
  )
}

## NETWAS format
res <- genes(
  edb,
  filter = GeneNameFilter(unique(multiomics$Gene)),
  return.type = "data.frame"
)

setDT(res)

out2 <- res[
  multiomics,
  on = .(gene_name = Gene)
]

out2$chr<-as.numeric(out2$seq_name)

out3<-out2[!is.na(out2$chr),]
out3$rep<-1000

for (f in unique(multiomics$folds)){
out33<-out3[out3$folds==f,c("chr","gene_name","n_cpgs","rep","gene_seq_start","gene_seq_end","z","p","gene_id","proteomics_P")]
write.table(out33,paste0("Three_omics_metal_ssizewts_fold",f,"_forNetwas_globalcog_slope.txt"),quote=T,col=F,row=F,sep="\t")
}


q('no')
