#!/usr/bin/env python

import sys
import pandas as pd

ase_homref = sys.argv[1]
mat_snp_vcf = sys.argv[2]
homref_matGT = sys.argv[3]

ase = pd.read_csv(ase_homref,sep="\t",dtype=str)

mat_snp = pd.read_csv(mat_snp_vcf,sep='\t',dtype=str,comment='#',header=None,compression='gzip').iloc[:,[0,1,-1]]
mat_snp.columns=['contig','position','MatGTinfo']
mat_snp['anyMatGT'] = mat_snp.iloc[:,-1].str.split(':',expand=True).iloc[:,0]
ase.merge(mat_snp.drop('MatGTinfo',axis=1),on=['contig','position'],how='left').fillna('NA').to_csv(homref_matGT,sep='\t',index=False,compression='zip')

