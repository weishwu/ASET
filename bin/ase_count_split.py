#!/usr/bin/env python

import pandas as pd
import scipy
import sys

sample_id = sys.argv[1]
vcf_fn = sys.argv[2]
ase_sense = sys.argv[3]
ase_antisense = sys.argv[4]
ase_hetsnp_out_fn = sys.argv[5]
ase_homsnp_out_fn = sys.argv[6]
ase_homref_out_fn = sys.argv[7]

vcfdata=pd.read_csv(vcf_fn,header=None,comment="#",sep="\t",dtype=str,compression='gzip')
vcfdata.columns=['contig','position','variantID','refAllele','altAllele','QUAL','FILTER','INFO','FORMAT','GT']

hetsnp=vcfdata.loc[(vcfdata['GT'].str.startswith('0/1') | vcfdata['GT'].str.startswith('0|1') | vcfdata['GT'].str.startswith('1|0')),].iloc[:,0:5]
homsnp=vcfdata.loc[(vcfdata['GT'].str.startswith('1/1') | vcfdata['GT'].str.startswith('1|1')),].iloc[:,0:5]
homref=vcfdata.loc[(vcfdata['GT'].str.startswith('0/0') | vcfdata['GT'].str.startswith('0|0')),].iloc[:,0:5]

senseASE=pd.read_csv(ase_sense,sep="\t",dtype=str)
antisenseASE=pd.read_csv(ase_antisense,sep="\t",dtype=str)

senseASE['strand']="plus"
antisenseASE['strand']='minus'

ASE=pd.concat([senseASE,antisenseASE],axis=0).assign(RNAid = sample_id)

hetsnpASE=hetsnp.merge(ASE,on=list(hetsnp.columns))
homsnpASE=homsnp.merge(ASE,on=list(homsnp.columns))
homrefASE=homref.merge(ASE,on=list(homref.columns))

hetsnpASE.to_csv(ase_hetsnp_out_fn,sep='\t',index=False)
homsnpASE.to_csv(ase_homsnp_out_fn,sep='\t',index=False)
homrefASE.to_csv(ase_homref_out_fn,sep='\t',index=False)

