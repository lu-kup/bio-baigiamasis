#!/usr/bin/env python
# coding: utf-8

import pyranges as pr
import pandas as pd
import config

OUTPUT_RANGES_PATH = '../outputs/ATAC.tsv'
TARGET_GENE_NAME = 'Tmx3'
GENES_PATH = '../outputs/genes_rna_valid.tsv'
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

def get_open_vector(chr, output_ranges_path = '../outputs/expression/ranges_chr3_kmeans_TT_S0.tsv', genes_path = GENES_PATH):
    ranges = pd.read_csv(output_ranges_path, sep='\t')
    genes_rna_valid = pd.read_csv(genes_path, sep='\t')
    genes_rna_chr = genes_rna_valid[genes_rna_valid['Chromosome'] == int(chr)]
    open_fraction_title = 'open_' + output_ranges_path[-6:-4]

    gr = pr.PyRanges(ranges)
    gr_rna = pr.PyRanges(genes_rna_chr[['Chromosome', 'Start', 'End', 'gene_name', 'gene_id', 'Strand']])

    gr_neighbour = get_neighbourhood_vector(gr_rna)
    intersection = gr_neighbour.intersect(gr)

    intersection.Length = intersection.End - intersection.Start
    df_summed = intersection.df[['gene_name', 'Length']].groupby("gene_name", as_index=False).sum()
    df_summed[open_fraction_title] = df_summed['Length'] / 2000
    basic_data_chr = genes_rna_chr[['Chromosome', 'gene_name', 'gene_id']].copy().set_index('gene_name')

    gene_ocrs = basic_data_chr.join(df_summed.set_index('gene_name'))
    gene_ocrs[open_fraction_title] = gene_ocrs[open_fraction_title].fillna(0)
    gene_ocrs['Length'] = gene_ocrs['Length'].fillna(0)
    gene_ocrs['Stage'] = STAGES[output_ranges_path[-6:-4]]
    gene_ocrs['gene_name'] = gene_ocrs.index
    gene_ocrs = gene_ocrs.reset_index(drop=True)

    return gene_ocrs

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

def get_neighbourhood_vector(gene_ranges):
    positive_ranges = gene_ranges[gene_ranges.Strand == '+'].copy()
    negative_ranges = gene_ranges[gene_ranges.Strand == '-'].copy()

    p_ranges = positive_ranges.copy()
    n_ranges = negative_ranges.copy()

    p_ranges.Start = positive_ranges.Start - 1000
    p_ranges.End = positive_ranges.Start + 1000

    n_ranges.Start = negative_ranges.End - 1000
    n_ranges.End = negative_ranges.End + 1000

    if len(p_ranges) + len(n_ranges) != len(gene_ranges):
        raise ValueError('Missing strand data.')

    neighbour_df = pd.concat([p_ranges.df, n_ranges.df], axis=0)
    neighbour_ranges = pr.PyRanges(neighbour_df)

    return neighbour_ranges

def get_expression_data(gene_name, chr, stage):
    output_path ='../outputs/expression/ranges2_chr' + chr + '_kmeans_TT_'+ stage +'.tsv'
    gencode_path = '../inputs/gencode_chr' + chr + '_M25.gtf'
    return get_open_gene(gene_name, output_path, gencode_path)

def get_expression_vector(chr, stage, genes_path):
    output_path ='../outputs/expression/ranges2_chr' + chr + '_kmeans_TT_'+ stage +'.tsv'
    return get_open_vector(chr, output_path, genes_path)

if __name__ == '__main__':
    import sys

    print(get_open_gene(sys.argv[1]))
