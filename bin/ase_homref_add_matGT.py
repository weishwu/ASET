#!/usr/bin/env python

import sys
import pandas as pd

ase_homref_fn = sys.argv[1]
smp_mat_fn = sys.argv[2]
ase_homref_out = sys.argv[3]

smp_mat = pd.read_csv(smp_mat_fn,dtype=str).set_index('sample').to_dict(orient='index')

MatGTdata=[]

for i in smp_mat.keys():
  fn = list(smp_mat[i].values())[0]
  if not (pd.isnull(fn) | (fn == '.')):
    snps = pd.read_csv(fn,sep='\t',dtype=str,comment='#',header=None,compression='gzip').iloc[:,[0,1,-1]].set_axis(['contig','position','MatGTinfo'],axis=1).assign(RNAid = i)
    snps['anyMatGT'] = snps.iloc[:,-2].str.split(':',expand=True).iloc[:,0]
    MatGTdata.append(snps.drop('MatGTinfo',axis=1))

MatGTall=pd.concat(MatGTdata,axis=0)

if len(MatGTall) == 0:
   pd.DataFrame(['.','.','.','.']).T.set_axis(['contig','position','RNAid','anyMatGT'],axis=1)

pd.read_csv(ase_homref_fn,sep="\t",dtype=str).merge(MatGTall,on=['contig','position','RNAid'],how='left').fillna('NA').to_csv(ase_homref_out,sep='\t',index=False, compression='gzip')

