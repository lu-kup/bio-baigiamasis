#!/usr/bin/env python
# coding: utf-8

import pyreadr
import pandas as pd
import pyranges as pr

CHR18_SIZE = 90720763
CHR18_SUBSET_SIZE = 9717594
CHR18_SUBSET_SIZE_ADJ = CHR18_SUBSET_SIZE - 3000000
SIGNAL_LABEL = 'TT_S0'

def evaluate(filename = "../outputs/output_prototypes.csv", dataframe = None, model_name = ''):
    if dataframe is None:
        df = pd.read_csv(filename, sep = '\t', index_col=0)
    else:
        df = dataframe

    input = pr.PyRanges(df)

    count0 = len(input[input.labels == 0])
    count1 = len(input[input.labels == 1])

    if count0 > count1:
        open_label = 1
        closed_label = 0
    else:
        open_label = 0
        closed_label = 1

    open_ranges = input[input.labels == open_label]
    merged_ranges = open_ranges.merge()

    df2 = pyreadr.read_r('../inputs/ATAC_dt.RDS')[None]
    df2.rename(columns={"start": "Start", "end": "End", "seqnames": "Chromosome"}, inplace=True)
    atac = pr.PyRanges(df2)

    df3 = pyreadr.read_r('../inputs/DNAseq_dt.RDS')[None]
    df3.rename(columns={"start": "Start", "end": "End", "seqnames": "Chromosome"}, inplace=True)
    dnaseq = pr.PyRanges(df3)

    m_lengths = merged_ranges.lengths()
    merged_ranges.Length = merged_ranges.End - merged_ranges.Start

    evaluation = merged_ranges.Length.describe()

    atac_intersection = atac.intersect(merged_ranges)
    atac_int_length = atac_intersection.length
    # Edited for subset
    atac_length = atac[(atac.Chromosome == '18') & (atac.Start < CHR18_SUBSET_SIZE)].length
    percentage_overlap = atac_int_length / atac_length
    evaluation['ATACseq_overlap'] = percentage_overlap


    dnaseq_intersection = merged_ranges.intersect(dnaseq)
    dnaseq_int_length = dnaseq_intersection.length
    # Edited for subset
    dnaseq_length = dnaseq[(dnaseq.Chromosome == '18') & (dnaseq.Start < CHR18_SUBSET_SIZE)].length
    percentage_overlap2 = dnaseq_int_length / dnaseq_length
    evaluation['DNAseq_overlap'] = percentage_overlap2

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
    evaluation['percentage_open'] = output_length / CHR18_SUBSET_SIZE
    evaluation['percentage_open_adj'] = output_length / CHR18_SUBSET_SIZE_ADJ
    evaluation['model_name'] = model_name

    evaluation['min_open_signal'] = df.loc[df["labels"] == open_label, SIGNAL_LABEL].min()
    evaluation['max_open_signal'] = df.loc[df["labels"] == open_label, SIGNAL_LABEL].max()
    evaluation['min_closed_signal'] = df.loc[df["labels"] == closed_label, SIGNAL_LABEL].min()
    evaluation['max_closed_signal'] = df.loc[df["labels"] == closed_label, SIGNAL_LABEL].max()

    return evaluation

if __name__ == "__main__":
    print(evaluate())