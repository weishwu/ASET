#!/usr/bin/env python

import pandas as pd
import sys

hetsnp_ase = sys.argv[1]
ase_x_exons = sys.argv[2]
ase_x_genes = sys.argv[3]
ase_exon_label = sys.argv[4]

ase=pd.read_csv(hetsnp_ase,sep='\t',low_memory=False,compression='gzip')
exon_all=pd.read_csv(ase_x_exons,header=None,sep='\t',low_memory=False)
genes=pd.read_csv(ase_x_genes,header=None,sep='\t',low_memory=False)

# all exons
exon_all.columns=['contig','position','strand','exons_merged']
exon_all.loc[exon_all['strand']=="+",'strand']="plus"
exon_all.loc[exon_all['strand']=="-",'strand']="minus"
exon_all_merged=exon_all.groupby(['contig','position','strand'])['exons_merged'].apply(";".join).reset_index()

# genes
genes.columns=['contig','position','strand','genes']
genes.loc[genes['strand']=="+",'strand']="plus"
genes.loc[genes['strand']=="-",'strand']="minus"
genes_merged=genes.groupby(['contig','position','strand'])['genes'].apply(";".join).reset_index()

aseExon=ase.merge(exon_all_merged,on=['contig','position','strand'],how="left").merge(genes_merged,on=['contig','position','strand'],how="left")

def gene_from_exons(x):
   if pd.isnull(x):
      y='NA'
   else:
      y=';'.join([i.split(':')[4]+':'+i.split(':')[5] for i in x.split(';')])
   return(y)

def gene_id_from_exons(x):
   if pd.isnull(x):
      y='NA'
   else:
      y=':'.join([i.split(';')[4] for i in x.split(':')])
   return(y)

def gene_name_from_exons(x):
   if pd.isnull(x):
      y='NA'
   else:
      y=':'.join([i.split(';')[5] for i in x.split(':')])
   return(y)

#aseExon['gene_id_from_exons'] = aseExon['exons_merged'].apply(gene_id_from_exons)
#aseExon['gene_name_from_exons'] = aseExon['exons_merged'].apply(gene_name_from_exons)
aseExon['genes_exonic'] = aseExon['exons_merged'].apply(gene_from_exons)

aseExon.fillna("NA").to_csv(ase_exon_label,sep="\t",index=False,compression='gzip')

