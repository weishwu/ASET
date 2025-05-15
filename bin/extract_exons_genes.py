#!/usr/bin/env python

import pandas as pd
import sys

gtf_df = sys.argv[1]
exon_out = sys.argv[2]
gene_out = sys.argv[3]

gtf=pd.read_csv(gtf_df,comment='#',sep='\t', low_memory=False)

gtf=gtf.loc[gtf['feature'] == 'exon'].copy()
gtf['gene_id_name']=gtf['gene_id']+':'+gtf['gene_name']
gtf['start_1']=gtf['start'].astype('int')-1
exons=gtf[['chr','start_1','end','transcript_id','score1','strand','gene_id_name','gene_id','gene_name']].copy()
exons.to_csv(exon_out, sep='\t',index=False, header=None)

genes = gtf.groupby(['chr','gene_id_name','score1','strand','gene_type','gene_id','gene_name'])['start_1'].agg('min').reset_index().merge(gtf.groupby(['chr','gene_id_name','score1','strand','gene_type','gene_id','gene_name'])['end'].agg('max').reset_index(), on=['chr','gene_id_name','score1','strand','gene_type','gene_id','gene_name'])[['chr','start_1','end','gene_id_name','score1','strand','gene_type','gene_id','gene_name']]
genes.to_csv(gene_out, sep='\t',index=False, header=None)

