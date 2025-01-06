#!/usr/bin/env python
# coding: utf-8

import pyreadr
import pandas as pd
import pyranges as pr
import config
from matplotlib import cbook

chr = config.CHR
chr_size = config.CHR_SIZE
gencode_data_path = config.GENCODE_DATA_PATH

def evaluate(merged_ranges, model_name = '', chr_subset_size = chr_size):
    chr_subset_size_adj = chr_subset_size - 3000000

    df2 = pyreadr.read_r('../inputs/ATAC_dt.RDS')[None]
    df2.rename(columns={"start": "Start", "end": "End", "seqnames": "Chromosome"}, inplace=True)
    atac = pr.PyRanges(df2)

    df3 = pyreadr.read_r('../inputs/DNAseq_dt.RDS')[None]
    df3.rename(columns={"start": "Start", "end": "End", "seqnames": "Chromosome"}, inplace=True)
    dnaseq = pr.PyRanges(df3)

    if len(merged_ranges) == 0:
        print("The model failed to identify open ranges.")
        data = [model_name]
        row_labels = ['model_name']
        return pd.DataFrame(data, index=row_labels)

    merged_ranges.Length = merged_ranges.End - merged_ranges.Start
    evaluation = merged_ranges.Length.describe()
    evaluation['mode'] = merged_ranges.Length.mode().iloc[0].item()

    atac_intersection = atac.intersect(merged_ranges)
    atac_int_length = atac_intersection.length
    # Edited for subset
    atac_length = atac[(atac.Chromosome == chr) & (atac.Start < chr_subset_size)].length
    percentage_overlap = atac_int_length / atac_length
    evaluation['ATACseq_overlap'] = percentage_overlap

    percentage_overlap_ranges = atac_int_length / merged_ranges.length
    evaluation['ATACseq_ranges_overlap'] = percentage_overlap_ranges

    # Jaccard similarity
    atac_subset = atac[(atac.Chromosome == chr) & (atac.Start < chr_subset_size)]
    atac_union = atac_subset.set_union(merged_ranges)

    if atac_union.length == 0:
        evaluation['ATACseq_Jaccard'] = 0.0
    else:
        evaluation['ATACseq_Jaccard'] = atac_int_length / atac_union.length


    dnaseq_intersection = merged_ranges.intersect(dnaseq)
    dnaseq_int_length = dnaseq_intersection.length
    # Edited for subset
    dnaseq_length = dnaseq[(dnaseq.Chromosome == chr) & (dnaseq.Start < chr_subset_size)].length
    percentage_overlap2 = dnaseq_int_length / dnaseq_length
    evaluation['DNAseq_overlap'] = percentage_overlap2

    percentage_overlap_ranges2 = dnaseq_int_length / merged_ranges.length
    evaluation['DNAseq_ranges_overlap'] = percentage_overlap_ranges2

    # Jaccard similarity
    dnaseq_subset = dnaseq[(dnaseq.Chromosome == chr) & (dnaseq.Start < chr_subset_size)]
    dnaseq_union = dnaseq_subset.set_union(merged_ranges)

    if dnaseq_union.length == 0:
        evaluation['DNAseq_Jaccard'] = 0.0
    else:
        evaluation['DNAseq_Jaccard'] = dnaseq_int_length / dnaseq_union.length


    gencode = pr.read_gtf(gencode_data_path)

    # {gene,transcript,exon,CDS,UTR,start_codon,stop_codon,Selenocysteine}
    genes_chr = gencode[gencode.Feature == 'gene']
    genes_chr.Chromosome = genes_chr.Chromosome.str.lower().str.removeprefix('chr')
    merged_genes_chr = genes_chr.merge()

    genes_intersection = merged_ranges.intersect(merged_genes_chr)
    genes_int_length = genes_intersection.length

    output_length = merged_ranges.length
    percentage_genes = genes_int_length / output_length

    evaluation['genes_overlap'] = percentage_genes
    # Edited for subset
    evaluation['percentage_open'] = output_length / chr_subset_size
    evaluation['percentage_open_adj'] = output_length / chr_subset_size_adj
    evaluation['model_name'] = model_name
    evaluation['boxplot_json'] = get_boxplot_stats(merged_ranges.Length)

    return evaluation

def process_ranges(merged_ranges):
    if len(merged_ranges) == 0:
        print("Failed to identify open ranges.")
        return merged_ranges

    merged_ranges.Length = merged_ranges.End - merged_ranges.Start
    processed_ranges = merged_ranges[merged_ranges.Length >= 75]
    processed_ranges = processed_ranges.merge(slack = 50)

    return processed_ranges

def evaluate_processed(merged_ranges, model_name = '', chr_subset_size = chr_size):
    return evaluate(process_ranges(merged_ranges), model_name, chr_subset_size)

def count_clusters(filename = "../outputs/output_prototypes.csv", dataframe = None, model_name = ''):
    if dataframe is None:
        df = pd.read_csv(filename, sep = '\t', index_col=0)
    else:
        df = dataframe

    return pd.Series([max(df['labels']) + 1], name='no_clusters')

def get_boxplot_stats(series):
    stats = cbook.boxplot_stats(series)
    stats_df = pd.DataFrame(stats)
    stats_json = stats_df.loc[0].to_json()

    return stats_json

if __name__ == "__main__":
    print(evaluate())