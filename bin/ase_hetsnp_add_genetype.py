#!/usr/bin/env python

import pandas as pd
import sys

ase_hetsnp_fn = sys.argv[1]
genes_bed_fn = sys.argv[2]
out_fn = sys.argv[3]

d=pd.read_csv(ase_hetsnp_fn,compression='gzip',sep='\t',low_memory=False)
gd=d[['contig','position','strand','RNAid','exons_merged']].copy()

g=gd.loc[pd.isnull(gd['exons_merged'])==False]

k=pd.DataFrame(g['exons_merged'].str.split(';').apply(pd.Series,1).stack())
k.index=k.index.droplevel(-1)
k.columns=["exons_merged_split"]
gx=g.merge(k,left_index=True,right_index=True).reset_index().drop(['index'],axis=1)
gx['gene_id']=gx['exons_merged_split'].str.split(':',expand=True)[4]

t=pd.read_csv(genes_bed_fn,sep='\t',low_memory=False,header=None)
t.columns=['contig','start','end','gene_id_name','score','strand','gene_type','gene_id','gene_name']

gxt=gx.merge(t,on='gene_id')
gxt['exons_merged_split_genetype']=gxt['exons_merged_split'].astype(str) + ':' + gxt['gene_type'].astype(str)

gxtm=gxt.groupby(['contig_x','position','strand_x','RNAid'])['exons_merged_split_genetype'].apply(';'.join)
#gxtm=gxt.groupby(['contig_x','position','strand_x','RNAid'])['gene_type'].apply(':'.join)

dm0=d.drop(['exons_merged'],axis=1).merge(gxtm.reset_index().rename(columns={'contig_x':'contig','strand_x':'strand','exons_merged_split_genetype':'exons_merged'}),on=['contig','position','strand','RNAid'],how='outer')

dm1=dm0.loc[pd.isnull(dm0['exons_merged']),].copy()
dm=dm0.loc[pd.isnull(dm0['exons_merged'])==False,].copy()

dm['gene_type_exonic']=dm['exons_merged'].apply(lambda x: ';'.join(sorted(list(set([y.split(':')[-1] for y in x.split(';')])))))

dm1['gene_type_exonic']='.'
#dm1['exons_merged_geneType_code']='.'
#dm['exons_merged_geneType_code']='.'
#dm.loc[((pd.isnull(dm['gene_type_exonic'])==False) & (dm['gene_type_exonic'].isin(['protein_coding']))),['exons_merged_geneType_code']]='pc'
#dm.loc[((pd.isnull(dm['gene_type_exonic'])==False) & (dm['gene_type_exonic'].isin(['lncRNA']))),['exons_merged_geneType_code']]='lnc'
#dm.loc[((pd.isnull(dm['gene_type_exonic'])==False) & (dm['gene_type_exonic'].isin(['snoRNA']))),['exons_merged_geneType_code']]='sno'
#dm.loc[((pd.isnull(dm['gene_type_exonic'])==False) & (dm['gene_type_exonic'].isin(['snRNA']))),['exons_merged_geneType_code']]='sn'
#dm.loc[((pd.isnull(dm['gene_type_exonic'])==False) & (dm['gene_type_exonic'].isin(['miRNA']))),['exons_merged_geneType_code']]='mi'
#dm.loc[((pd.isnull(dm['gene_type_exonic'])==False) & (dm['gene_type_exonic'].isin(['lncRNA:snoRNA']))),['exons_merged_geneType_code']]='sno_lnc'
#dm.loc[((pd.isnull(dm['gene_type_exonic'])==False) & (dm['gene_type_exonic'].isin(['misc_RNA']))),['exons_merged_geneType_code']]='misc'
#dm.loc[((pd.isnull(dm['gene_type_exonic'])==False) & (dm['gene_type_exonic'].isin(['transcribed_unitary_pseudogene']))),['exons_merged_geneType_code']]='tutp'
#dm.loc[((pd.isnull(dm['gene_type_exonic'])==False) & (dm['gene_type_exonic'].isin(['transcribed_unprocessed_pseudogene']))),['exons_merged_geneType_code']]='tup'
#dm.loc[((pd.isnull(dm['gene_type_exonic'])==False) & (dm['gene_type_exonic'].isin(['transcribed_processed_pseudogene']))),['exons_merged_geneType_code']]='tpp'
#dm.loc[((pd.isnull(dm['gene_type_exonic'])==False) & (dm['gene_type_exonic'].isin(['processed_pseudogene']))),['exons_merged_geneType_code']]='pp'
#dm.loc[((pd.isnull(dm['gene_type_exonic'])==False) & (dm['gene_type_exonic'].isin(['unprocessed_pseudogene']))),['exons_merged_geneType_code']]='up'
#dm.loc[((pd.isnull(dm['gene_type_exonic'])==False) & (dm['gene_type_exonic'].isin(['TEC']))),['exons_merged_geneType_code']]='tec'
#dm.loc[((pd.isnull(dm['gene_type_exonic'])==False) & (dm['gene_type_exonic'].isin(['TEC:lncRNA']))),['exons_merged_geneType_code']]='lnc_tec'
#dm.loc[((pd.isnull(dm['gene_type_exonic'])==False) & (dm['gene_type_exonic'].isin(['lncRNA:protein_coding']))),['exons_merged_geneType_code']]='pc_lnc'
#dm.loc[((pd.isnull(dm['gene_type_exonic'])==False) & (dm['gene_type_exonic'].isin(['protein_coding:transcribed_processed_pseudogene']))),['exons_merged_geneType_code']]='pc_tpp'
#dm.loc[((pd.isnull(dm['gene_type_exonic'])==False) & (dm['gene_type_exonic'].isin(['processed_pseudogene:protein_coding']))),['exons_merged_geneType_code']]='pc_pp'

pd.concat([dm,dm1],axis=0).sort_values(['contig','position','strand','RNAid']).fillna('NA').to_csv(out_fn,index=False,sep='\t',compression='gzip')

