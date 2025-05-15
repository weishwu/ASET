#!/usr/bin/env python

import pandas as pd
import sys

ase_hetsnp_fn = sys.argv[1]
ase_homref_exonic_fn = sys.argv[2]
ase_dp_min = float(sys.argv[3])
out_fn = sys.argv[4]

homRef_exonic = pd.read_csv(ase_homref_exonic_fn,sep='\t',compression='gzip',low_memory=False)
hetSNP = pd.read_csv(ase_hetsnp_fn,sep='\t',compression='gzip',low_memory=False)
hetSNP['RNAid'] = hetSNP['RNAid'].astype('str')
homRef_exonic['RNAid'] = homRef_exonic['RNAid'].astype('str')

# add per gene per sample nonRefFreq mean
if len(homRef_exonic) > 0:
  homRef_exonic['RNAid'] = homRef_exonic['RNAid'].astype('str')
  homRef_exonic = homRef_exonic.loc[pd.isnull(homRef_exonic['MatGT']) == False].copy()
  homRef_exonic['MatGT_alt'] = homRef_exonic['MatGT'].apply(lambda x: len([i for i in x.split('/') if not i in ['0','.']]))
  homRef_exonic_perSampleGene = homRef_exonic.groupby(['gene','RNAid'])['nonRefFreq'].agg('mean').reset_index().rename(columns = {'nonRefFreq':'homRef_nonRefFreq_mean_perGene_perRNAid'})
  homRefMat_exonic_perSampleGene = homRef_exonic.loc[homRef_exonic['MatGT_alt'] > 0,].groupby(['gene','RNAid'])['nonRefFreq'].agg('mean').reset_index().rename(columns = {'nonRefFreq':'homRef_nonRefFreq_atMatAlt_mean_perGene_perRNAid'})
  
  hetSNP_genic = hetSNP.loc[((pd.isnull(hetSNP['genes_exonic']) == False) & (hetSNP['genes_exonic'].str.contains(';') == False))].copy()
  hetSNP_other = hetSNP.loc[((pd.isnull(hetSNP['genes_exonic']) == True) | (hetSNP['genes_exonic'].str.contains(';') == True))].copy()

  hetSNP_genic['gene'] = hetSNP_genic['genes_exonic'].str.split(':',expand=True)[0]  
  hetSNP_genic = hetSNP_genic.merge(homRef_exonic_perSampleGene, on=['RNAid','gene'], how='left').merge(homRefMat_exonic_perSampleGene, on=['RNAid','gene'], how='left').drop(['gene'],axis=1).copy()
  
  hetSNP_other['homRef_nonRefFreq_mean_perGene_perRNAid'] = float('nan')
  hetSNP_other['homRef_nonRefFreq_atMatAlt_mean_perGene_perRNAid'] = float('nan')
  
  hetSNP = pd.concat([hetSNP_genic, hetSNP_other], axis=0).sort_values(by=['RNAid','contig','position','strand']).copy()
else:
  hetSNP['homRef_nonRefFreq_mean_perGene_perRNAid'] = float('nan')
  hetSNP['homRef_nonRefFreq_atMatAlt_mean_perGene_perRNAid'] = float('nan')

hetSNP.drop(['homRef_nonRefFreq_mean_perGene_perRNAid'],axis=1).fillna('NA').to_csv(out_fn,sep='\t',index=False,compression='gzip')

