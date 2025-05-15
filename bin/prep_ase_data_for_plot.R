#!/usr/bin/env Rscript

library(dplyr)

args=commandArgs(TRUE)

ase_hetsnp_fn = args[1] # 'ASE_on_hetSNPs.gz'
exons_bed_fn = args[2] #'../ref_data/gtf.exonMergedByGenes.bed'
ase_dp_min = as.numeric(args[3]) # 10
out_fn = args[4]

##### filter SNP data
ase_data_full = read.table(gzfile(ase_hetsnp_fn),header=T)

ase_selc = ase_data_full %>% filter( 
    (totalCount >= ase_dp_min) & 
    (!is.na(gene_contaminated_0.05)) & 
    (gene_contaminated_0.05 == 0) & 
    (!is.na(mergedExons_all)) & 
    (!grepl(':', mergedExons_all)) &
    (nonAltFreq_perRNAid < 0.05)) %>% select(
    RNAid,contig,position,totalCount,strand,rawASE,mergedExons_all,PatFreq,PatDepth,MatDepth,gene_from_exons) %>% rename(gene = gene_from_exons)

ase_selc = cbind(ase_selc, t(matrix(unlist(lapply(ase_selc$gene, function(x){strsplit(x,';')[[1]]})),2,)) %>% `colnames<-`(c('gene_id','gene_name')))

ase_selc_phased = subset(ase_selc, !is.na(PatFreq))

patfreq_bygene = aggregate(ase_selc_phased[,which(colnames(ase_selc_phased)%in%c('PatDepth','MatDepth','totalCount'))],by=list(ase_selc_phased$RNAid,ase_selc_phased$gene_id,ase_selc_phased$gene_name),sum)

colnames(patfreq_bygene)[1:3] = c('RNAid','gene_id','gene_name')

patfreq_bygene$PatFreq_PerGenePerRNAid = patfreq_bygene$PatDepth/patfreq_bygene$totalCount

exons_all = read.table(exons_bed_fn)
exons_all = as.data.frame(unique(matrix(unlist(strsplit(as.character(exons_all[,4]),";")),nrow(exons_all),6,byrow=T)))
exons_all[,2] = as.numeric(exons_all[,2])
exons_all[,3] = as.numeric(exons_all[,3])

ase_data = list(ase_selc,ase_selc_phased,patfreq_bygene,exons_all)
names(ase_data) = c('ase_filt','ase_filt_phased','patfreq_bygene','exons_all')
saveRDS(ase_data, file = out_fn)

