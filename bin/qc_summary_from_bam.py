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
rnaIDs = os.path.basename(list(set(ase_hetSNP['RNAid'].to_list())))

sample_vcfs = pd.read_csv(sample_fn,dtype=str)[['sample','snps']].set_index('sample').to_dict(orient='index')

qc_allsamples={}

def getQC(rnaID):
  qc_data = {}
  qc_data['RNAid'] = rnaID
  allhcFlgst = rnaID + '.coord_sorted.hc.bam.flagstat'
  deduppedhcFlagst = rnaID + '.aln.hc.dedup.bam.flagstat'
  picardMetrics = rnaID + '.hc.bam.RNAseqMetrics.txt'
  hetSNPvcf = list(sample_vcfs[rnaID].values())[0]

  qc_data['hcAln']=int([x for x in open(allhcFlgst) if 'mapped (' in x][0].split()[0])
  qc_data['dedupped_hcAlan']=int([x for x in open(deduppedhcFlagst) if 'mapped (' in x][0].split()[0])
  qc_data['hetSNP_total'] = sum(pd.read_csv(hetSNPvcf,sep='\t',comment='#',header=None,compression='gzip',low_memory=False).iloc[:,9].apply(lambda x: x.split(':')[0] in ['0/1','0|1','1|0']))
  qc_data['hetSNP_DPpass'] = len(ase_hetSNP.loc[((ase_hetSNP['totalCount'] > ase_dp_min) & (ase_hetSNP['RNAid']==rnaID)),['contig','position']].drop_duplicates())
  qc_data['hetSNP_DPpass_exonic'] = len(ase_hetSNP.loc[((ase_hetSNP['totalCount'] >= ase_dp_min) & (ase_hetSNP['RNAid']==rnaID) & (pd.isnull(ase_hetSNP['exons_merged'])==False)),['contig','position']].drop_duplicates())
  qc_data['nonAltFreq_atHomSNP'] = ase_hetSNP.loc[ase_hetSNP['RNAid']==rnaID,'nonAltFreq_perRNAid'].drop_duplicates().iloc[0]
  qc_data['nonRefFreq_atHomRef'] = ase_hetSNP.loc[ase_hetSNP['RNAid']==rnaID,'nonRefFreq_perRNAid'].drop_duplicates().iloc[0]
  qc_data['dedupped_hcAlan_over_hcAln']=qc_data['dedupped_hcAlan']/qc_data['hcAln']
  qc_data['hetSNP_DPpass_over_total']=qc_data['hetSNP_DPpass']/qc_data['hetSNP_total']

  metrics = [x for x in open(picardMetrics) if len(x)>0]
  pfline = [x for x in range(0,len(metrics)) if metrics[x].startswith('PF_BASES')][0]
  for i,j in zip(metrics[pfline].split(),metrics[pfline+1].split()):
     qc_data[i]=j

  return(qc_data)

with Pool(n_core) as pool:
   qcAll = pool.map(getQC,rnaIDs)

qcAll = pd.DataFrame(qcAll)[['RNAid','hcAln','dedupped_hcAlan','dedupped_hcAlan_over_hcAln','hetSNP_total','hetSNP_DPpass','hetSNP_DPpass_over_total','hetSNP_DPpass_exonic','nonAltFreq_atHomSNP','nonRefFreq_atHomRef','PCT_RIBOSOMAL_BASES','PCT_CODING_BASES','PCT_UTR_BASES','PCT_INTRONIC_BASES','PCT_INTERGENIC_BASES','PCT_MRNA_BASES','MEDIAN_5PRIME_BIAS','MEDIAN_3PRIME_BIAS','MEDIAN_5PRIME_TO_3PRIME_BIAS','MEDIAN_CV_COVERAGE','PCT_CORRECT_STRAND_READS']]

qcAll.to_csv(out_fn,index=False,sep='\t')

