#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys

ase_data_fn = sys.argv[1]
sample_phase_fn = sys.argv[2]
out_fn = sys.argv[3]

smp_phase = pd.read_csv(sample_phase_fn,dtype=str).set_index('sample').to_dict(orient='index')

phasingData={}
for i in smp_phase.keys():
   fn = list(smp_phase[i].values())[0]
   if not (pd.isnull(fn) | (fn == '.')):
      phasingData[i] = pd.read_csv(fn,sep='\t',low_memory=False).assign(RNAid = i)

phasingDataMerged = pd.concat(phasingData.values(),axis=0)

aseDataMerged = pd.read_csv(ase_data_fn,sep='\t',compression='gzip',low_memory=False).merge(phasingDataMerged,on=['contig','position','refAllele','altAllele','RNAid'],how='left')
ase = aseDataMerged.loc[pd.isna(aseDataMerged['phasedGT'])==False,].copy()
aseNophased = aseDataMerged.loc[pd.isna(aseDataMerged['phasedGT']),].copy()

ase['PatAllele']=ase['refAllele']
ase['MatAllele']=ase['altAllele']
ase.loc[ase['phasedGT']=="1|0",'PatAllele']=ase.loc[ase['phasedGT']=="1|0",'altAllele']
ase.loc[ase['phasedGT']=="1|0",'MatAllele']=ase.loc[ase['phasedGT']=="1|0",'refAllele']

ase['PatDepth']=ase['refCount']
ase['MatDepth']=ase['altCount']
ase.loc[ase['phasedGT']=="1|0",'PatDepth']=ase.loc[ase['phasedGT']=="1|0",'altCount']
ase.loc[ase['phasedGT']=="1|0",'MatDepth']=ase.loc[ase['phasedGT']=="1|0",'refCount']

ase_part1=ase.loc[ase['totalCount']>0].copy()
ase_part2=ase.loc[ase['totalCount']==0].copy()

ase_part1['PatFreq']=ase_part1['PatDepth']/ase_part1['totalCount']
ase_part1['MatFreq']=ase_part1['MatDepth']/ase_part1['totalCount']
ase_part2['PatFreq']=np.NaN
ase_part2['MatFreq']=np.NaN


for f in ['PatAllele','MatAllele','PatDepth','MatDepth','PatFreq','MatFreq']:
    aseNophased[f]=np.NaN

colsIN=list(ase_part1.columns)
alldata=pd.concat([ase_part1,ase_part2,aseNophased],axis=0).sort_values(by=["RNAid","contig","position","strand"])[colsIN]

alldata.fillna("NA").to_csv(out_fn,sep="\t",index=False,compression='gzip')

