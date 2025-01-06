#!/usr/bin/env python
# coding: utf-8

import pyranges as pr
import pandas as pd
import config

OUTPUT_RANGES_PATH = '../outputs/ATAC.tsv'
TARGET_GENE_NAME = 'Tmx3'
#TARGET_GENE_NAME = 'Cul2'

gencode_data_path = config.GENCODE_DATA_PATH
chr = config.CHR

ranges = pd.read_csv(OUTPUT_RANGES_PATH, sep='\t')
input_ranges = pr.PyRanges(ranges)
genes = pr.read_gtf(gencode_data_path)
genes = genes[genes.Feature == 'gene']

genes_df = genes.df
genes_df = genes_df[genes_df["gene_name"] == TARGET_GENE_NAME]
gene_ranges = pr.PyRanges(genes_df)
gene_ranges.Chromosome = gene_ranges.Chromosome.str.lower().str.removeprefix('chr')

overlap_dict = {'gene':gene_ranges, 'ranges':input_ranges}
overlaps = pr.count_overlaps(overlap_dict)
overlaps_df = overlaps.df
overlapping_genes = pr.PyRanges(overlaps_df[(overlaps_df['gene'] != 0)])

print("\n\nTARGET_GENE_NAME:", TARGET_GENE_NAME)
print("OUTPUT_RANGES_PATH:", OUTPUT_RANGES_PATH)
print(overlapping_genes)

