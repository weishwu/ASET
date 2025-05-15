#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys

hetsnp_fn = sys.argv[1]
out_fn = sys.argv[2]

g=pd.read_csv(hetsnp_fn,sep='\t',compression='gzip',low_memory=False)
g['RNAid']=g['RNAid'].astype('str')
k=g.loc[g['totalCount']>0,].copy()
k['rawASE']=abs(0.5-k['refCount']/k['totalCount'])

d=k.loc[((k['totalCount']>9) & (k['nonAltFreq_perRNAid']<0.05) & (pd.isnull(k['genes_all'])==False) & (k['genes_all'].str.contains(':')==False)),].copy()

gexonic=d.loc[pd.isnull(d['mergedExons_all'])==False,].copy()
gintronic=d.loc[pd.isnull(d['mergedExons_all'])==True,].copy()

t=pd.concat([gexonic,gintronic])[['contig','position','strand','RNAid','rawASE','genes_all']].copy()
h=t.groupby(['RNAid','genes_all'])['rawASE'].agg(['mean','std','count']).reset_index().rename(columns={'mean':'geneASE_mean','std':'geneASE_stdev','count':'geneSNP_count'})
h['geneASE_mean_plus2std']=h['geneASE_mean']+2*h['geneASE_stdev']
h['geneASE_mean_minus2std']=h['geneASE_mean']-2*h['geneASE_stdev']

dx=d.merge(h,on=['RNAid','genes_all'])
dx['geneASE_outlier']=dx.apply(lambda x: 1 if ((x['rawASE']>x['geneASE_mean_plus2std']) | (x['rawASE']<x['geneASE_mean_minus2std'])) else 0,axis=1)
dx=dx[['contig','position','strand','RNAid','geneASE_mean','geneASE_stdev','geneSNP_count','geneASE_mean_plus2std','geneASE_mean_minus2std','geneASE_outlier']]

kx=k.merge(dx,on=['contig','position','strand','RNAid'],how='left')
kx.fillna('NA').to_csv(out_fn,sep='\t',index=False,compression='gzip')

