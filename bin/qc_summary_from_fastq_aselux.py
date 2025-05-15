#!/usr/bin/env python

import pandas as pd
import glob
import sys
import os
from multiprocessing import Pool

ase_hetsnp_fn = sys.argv[1]
sample_fn = sys.argv[2]
out_fn = sys.argv[3]
ase_dp_min = float(sys.argv[4])
n_cores = int(sys.argv[5]) if len(sys.argv) > 4 else 1

ase_hetSNP = pd.read_csv(ase_hetsnp_fn,sep='\t',compression='gzip',low_memory=False)
ase_hetSNP['RNAid'] = ase_hetSNP['RNAid'].astype('str')
rnaIDs = list(set(ase_hetSNP['RNAid'].to_list()))

sample_vcfs = pd.read_csv(sample_fn,dtype=str).iloc[:,[0,3]].set_index('sample').to_dict(orient='index')

def getQC(rnaID):
  qc_data = {}
  qc_data['RNAid'] = rnaID
  trimsummary = 'Sample_'+rnaID+'_trimmomatic.log'
  hetSNPvcf = os.path.basename(list(sample_vcfs[rnaID].values())[0])

  qc_data['rawReadPairs'] = int([x for x in open(trimsummary) if x.startswith('Input Read Pairs:')][0].split(':')[1].split()[0])
  qc_data['trimmedReadPairs'] = int([x for x in open(trimsummary) if x.startswith('Input Read Pairs:')][0].split(':')[2].split()[0])
  qc_data['hetSNP_total'] = sum(pd.read_csv(hetSNPvcf,sep='\t',comment='#',header=None,compression='gzip',low_memory=False).iloc[:,9].apply(lambda x: x.split(':')[0] in ['0/1','0|1','1|0']))
  qc_data['hetSNP_DPpass'] = len(ase_hetSNP.loc[((ase_hetSNP['totalCount'] > ase_dp_min) & (ase_hetSNP['RNAid']==rnaID)),['contig','position']].drop_duplicates())
  qc_data['hetSNP_DPpass_exonic'] = len(ase_hetSNP.loc[((ase_hetSNP['totalCount'] >= ase_dp_min) & (ase_hetSNP['RNAid']==rnaID) & (pd.isnull(ase_hetSNP['exons_merged'])==False)),['contig','position']].drop_duplicates())
  qc_data['nonAltFreq_atHomSNP'] = 'NA'
  qc_data['nonRefFreq_atHomRef'] = 'NA'
  qc_data['hetSNP_DPpass_over_total'] = qc_data['hetSNP_DPpass']/qc_data['hetSNP_total']
  return(qc_data)

with Pool(n_cores) as pool:
    qcAll = pool.map(getQC, rnaIDs)

qcAll = pd.DataFrame(qcAll)[['RNAid','rawReadPairs','trimmedReadPairs','hetSNP_total','hetSNP_DPpass','hetSNP_DPpass_over_total','hetSNP_DPpass_exonic','nonAltFreq_atHomSNP','nonRefFreq_atHomRef']]

qcAll.to_csv(out_fn,index=False,sep='\t')


