#!/usr/bin/env python

import pandas as pd
import scipy
import sys

ase_hetsnp_fn = sys.argv[1]
ase_homsnp_fn = sys.argv[2]
ase_homref_fn = sys.argv[3]
ase_dp_min = float(sys.argv[4])
out_fn = sys.argv[5]

# load ASE of hetSNP, homSNP and homRef
hetSNP = pd.read_csv(ase_hetsnp_fn,sep='\t',low_memory=False)
homSNP = pd.read_csv(ase_homsnp_fn,sep='\t',low_memory=False)
homRef = pd.read_csv(ase_homref_fn,sep='\t',low_memory=False)
hetSNP['RNAid'] = hetSNP['RNAid'].astype('str')
homSNP['RNAid'] = homSNP['RNAid'].astype('str')
homRef['RNAid'] = homRef['RNAid'].astype('str')

# assign higher and lower allele read count
hetSNP['hiCount'] = hetSNP[['refCount','altCount']].max(axis=1)
hetSNP['loCount'] = hetSNP[['refCount','altCount']].min(axis=1)

# add ASE and p-value
hetSNP = hetSNP.loc[hetSNP['totalCount'] > 0].copy()
hetSNP['rawASE'] = abs(0.5 - hetSNP['refCount']/hetSNP['totalCount'])
hetSNP['rawASE_pvalue'] = hetSNP.apply(lambda x: scipy.stats.binomtest(x['refCount'],x['totalCount'],p=0.5, alternative='two-sided').pvalue, axis=1)

# add per-sample non-ref freq means
if len(homRef) > 0:
  homRef_dpPass = homRef.loc[(homRef['totalCount'] + homRef['otherBases']) >= ase_dp_min,].copy()
  homRef_dpPass['nonRefFreq_perRNAid'] = 1 - homRef_dpPass['refCount']/(homRef_dpPass['totalCount'] + homRef_dpPass['otherBases'])
  homRef_byRNAid = homRef_dpPass.groupby(['RNAid'])['nonRefFreq_perRNAid'].agg('mean').reset_index()
else:
  homRef_byRNAid = pd.DataFrame(['.','.']).T.set_axis(['RNAid','nonRefFreq_perRNAid'],axis=1)

# add per-sample non-alt freq means
if len(homSNP) > 0:
  homSNP_dpPass = homSNP.loc[(homSNP['totalCount']+homSNP['otherBases']) >= ase_dp_min,].copy()
  homSNP_dpPass['nonAltFreq_perRNAid'] = 1 - homSNP_dpPass['altCount']/(homSNP_dpPass['totalCount'] + homSNP_dpPass['otherBases'])
  homSNP_byRNAid = homSNP_dpPass.groupby(['RNAid'])['nonAltFreq_perRNAid'].agg('mean').reset_index()
else:
  homSNP_byRNAid = pd.DataFrame(['.','.']).T.set_axis(['RNAid','nonAltFreq_perRNAid'],axis=1)

hetSNP = hetSNP.merge(homSNP_byRNAid,on='RNAid',how='left').merge(homRef_byRNAid,on='RNAid',how='left')

hetSNP.fillna('NA').to_csv(out_fn,sep='\t',index=False,compression='gzip')


