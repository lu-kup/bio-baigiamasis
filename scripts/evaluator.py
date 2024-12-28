#!/usr/bin/env python
# coding: utf-8

import pyreadr
import pandas as pd
import pyranges as pr

CHR18_SIZE = 90720763
CHR18_SUBSET_SIZE = 9717594

def evaluate(merged_ranges, model_name = '', chr_subset_size = CHR18_SUBSET_SIZE):
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

    atac_intersection = atac.intersect(merged_ranges)
    atac_int_length = atac_intersection.length
    # Edited for subset
    atac_length = atac[(atac.Chromosome == '18') & (atac.Start < chr_subset_size)].length
    percentage_overlap = atac_int_length / atac_length
    evaluation['ATACseq_overlap'] = percentage_overlap

    percentage_overlap_ranges = atac_int_length / merged_ranges.length
    evaluation['ATACseq_ranges_overlap'] = percentage_overlap_ranges

    dnaseq_intersection = merged_ranges.intersect(dnaseq)
    dnaseq_int_length = dnaseq_intersection.length
    # Edited for subset
    dnaseq_length = dnaseq[(dnaseq.Chromosome == '18') & (dnaseq.Start < chr_subset_size)].length
    percentage_overlap2 = dnaseq_int_length / dnaseq_length
    evaluation['DNAseq_overlap'] = percentage_overlap2

    percentage_overlap_ranges2 = dnaseq_int_length / merged_ranges.length
    evaluation['DNAseq_ranges_overlap'] = percentage_overlap_ranges2

    gencode = pr.read_gtf("../inputs/gencode_chr18.gtf")

    # {gene,transcript,exon,CDS,UTR,start_codon,stop_codon,Selenocysteine}
    genes18 = gencode[gencode.Feature == 'gene']
    genes18.Chromosome = '18'
    merged_genes18 = genes18.merge()

    genes_intersection = merged_ranges.intersect(merged_genes18)
    genes_int_length = genes_intersection.length

    output_length = merged_ranges.length
    percentage_genes = genes_int_length / output_length


    evaluation['genes_overlap'] = percentage_genes
    # Edited for subset
    evaluation['percentage_open'] = output_length / chr_subset_size
    evaluation['percentage_open_adj'] = output_length / chr_subset_size_adj
    evaluation['model_name'] = model_name

    return evaluation

def process_ranges(merged_ranges):
    if len(merged_ranges) == 0:
        print("Failed to identify open ranges.")
        return merged_ranges

    merged_ranges.Length = merged_ranges.End - merged_ranges.Start
    processed_ranges = merged_ranges[merged_ranges.Length >= 75]
    processed_ranges = processed_ranges.merge(slack = 50)

    return processed_ranges

def evaluate_processed(merged_ranges, model_name = '', chr_subset_size = CHR18_SUBSET_SIZE):
    return evaluate(process_ranges(merged_ranges), model_name, chr_subset_size)

def count_clusters(filename = "../outputs/output_prototypes.csv", dataframe = None, model_name = ''):
    if dataframe is None:
        df = pd.read_csv(filename, sep = '\t', index_col=0)
    else:
        df = dataframe

    return pd.Series([max(df['labels']) + 1], name='no_clusters')

if __name__ == "__main__":
    print(evaluate())