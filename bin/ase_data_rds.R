#!/usr/bin/env Rscript

args=commandArgs(TRUE)

ase_hetsnp_fn = args[1] # 'ASE_on_hetSNPs.gz'
exons_bed_fn = args[2] #'../ref_data/gtf.exonMergedByGenes.bed'
out_fn = args[3]

##### filter SNP data
ase_data_full = read.table(gzfile(ase_hetsnp_fn),header=T)

exons_all = read.table(exons_bed_fn)
exons_all = as.data.frame(unique(matrix(unlist(strsplit(as.character(exons_all[,4]),":")),nrow(exons_all),6,byrow=T)))
exons_all[,2] = as.numeric(exons_all[,2])
exons_all[,3] = as.numeric(exons_all[,3])

ase_data = list(ase_data_full,exons_all)
names(ase_data) = c('ase_df','union_exons_per_gene')
saveRDS(ase_data, file = out_fn)

