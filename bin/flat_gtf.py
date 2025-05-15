#!/usr/bin/env python

import pandas as pd
import sys
from multiprocessing import Pool

gtf_fn = sys.argv[1]
gtf_df_out = sys.argv[2]
n_cores = int(sys.argv[3]) if len(sys.argv) > 3 else 1

gtf = pd.read_csv(gtf_fn, comment='#', sep='\t', header=None)
gtf.columns = ['chr', 'group', 'feature', 'start', 'end', 'score1', 'strand', 'score2', 'info']

def get_fields(info_line):
    field_dict = {}
    tokens = info_line.split()
    for i, token in enumerate(tokens):
        if token in all_fields:
            val = tokens[i+1].replace('"','').replace(';','')
            field_dict[token] = val
    return field_dict

all_fields = ['gene_id', 'transcript_id', 'gene_type', 'gene_name', 'transcript_type', 'transcript_name',
              'exon_number', 'exon_id', 'level', 'transcript_support_level', 'hgnc_id', 'tag',
              'havana_gene', 'havana_transcript']

with Pool(n_cores) as pool:
    records = pool.map(get_fields, gtf['info'].tolist())

fields_df = pd.DataFrame(records)
# Ensure all expected columns are present
fields_df = fields_df.reindex(columns=all_fields)

_pd_all = pd.concat([gtf[['chr', 'group', 'feature', 'start', 'end', 'score1', 'strand', 'score2']], fields_df], axis=1).fillna('.')
_pd_all.to_csv(gtf_df_out, sep='\t', index=False)

