#!/usr/bin/env python
# coding: utf-8

import pyranges as pr
import pandas as pd
import config

OUTPUT_RANGES_PATH = '../outputs/ATAC.tsv'
TARGET_GENE_NAME = 'Tmx3'
GENCODE_DATA_PATH = config.GENCODE_DATA_PATH

STAGES = {'S0' : 'ESC', 'S1' : 'NPC', 'S2' : 'NC', 'AC' : 'TEST'}

def get_open_gene(target_gene_name = TARGET_GENE_NAME, output_ranges_path = OUTPUT_RANGES_PATH, gencode_data_path = GENCODE_DATA_PATH):
    ranges = pd.read_csv(output_ranges_path, sep='\t')
    input_ranges = pr.PyRanges(ranges)
    genes = pr.read_gtf(gencode_data_path)
    genes = genes[genes.Feature == 'gene']

    genes_df = genes.df
    genes_df = genes_df[genes_df["gene_name"] == target_gene_name]
    if len(genes_df) == 0:
        raise ValueError(f'No genes found for: {target_gene_name}')
    gene_ranges = pr.PyRanges(genes_df)
    gene_ranges.Chromosome = gene_ranges.Chromosome.str.lower().str.removeprefix('chr')
    gene_neighbourhood = get_neighbourhood(gene_ranges, target_gene_name)

    print(input_ranges)

    overlap_dict = {'gene_neighb':gene_neighbourhood, 'ranges':input_ranges}
    overlaps = pr.count_overlaps(overlap_dict)
    overlaps_df = overlaps.df
    overlapping_genes = pr.PyRanges(overlaps_df[(overlaps_df['gene_neighb'] != 0)])

    print("\n\nTARGET_GENE_NAME:", target_gene_name)
    print(gene_ranges)
    print("\n\nGENE_NEIGHBOURHOOD:", target_gene_name)
    print(gene_neighbourhood)
    print("OUTPUT_RANGES_PATH:", output_ranges_path)
    print(overlapping_genes)

    overlapping_genes = pr.PyRanges(overlaps_df[(overlaps_df['gene_neighb'] != 0)])
    open_length = overlapping_genes[overlapping_genes.ranges == 1].length
    total_length = overlapping_genes.length

    output = {'gene' : target_gene_name, 'open': (open_length / total_length), 'stage' : STAGES[output_ranges_path[-6:-4]]}
    return pd.DataFrame([output])

def get_neighbourhood(gene_ranges, target_gene_name):
    if len(gene_ranges) != 1:
        raise ValueError(f'Multiple gene ranges identified: {target_gene_name}')

    neighbour_ranges = gene_ranges.copy()

    if gene_ranges.Strand[0] == '+':
        neighbour_ranges.Start = gene_ranges.Start - 1000
        neighbour_ranges.End = gene_ranges.Start + 1000
    elif gene_ranges.Strand[0] == '-':
        neighbour_ranges.Start = gene_ranges.End - 1000
        neighbour_ranges.End = gene_ranges.End + 1000
    else:
        raise ValueError('Invalid strand type.')

    return neighbour_ranges

def get_expression_data(gene_name, chr, stage):
    output_path ='../outputs/expression/ranges2_chr' + chr + '_kmeans_TT_'+ stage +'.tsv'
    gencode_path = '../inputs/gencode_chr' + chr + '_M25.gtf'
    return get_open_gene(gene_name, output_path, gencode_path)

if __name__ == '__main__':
    import sys

    print(get_open_gene(sys.argv[1]))
