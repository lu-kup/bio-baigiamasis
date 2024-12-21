#!/usr/bin/env python
# coding: utf-8

import pyreadr
import pandas as pd
import pyranges as pr

CHR18_SIZE = 90720763
CHR18_SUBSET_SIZE = 25556940
CHR18_SUBSET_SIZE_ADJ = 25556940 - 10000018

def evaluate(filename = "output_prototypes.csv", dataframe = None):
    if dataframe is None:
        df = pd.read_csv(filename, sep = '\t', index_col=0)
    else:
        df = dataframe

    input = pr.PyRanges(df)

    count0 = len(input[input.labels == 0])
    count1 = len(input[input.labels == 1])

    if count0 > count1:
        target_label = 1
    else:
        target_label = 0


    open_ranges = input[input.labels == target_label]
    merged_ranges = open_ranges.merge()

    df2 = pyreadr.read_r('ATAC_dt.RDS')[None]
    df2.rename(columns={"start": "Start", "end": "End", "seqnames": "Chromosome"}, inplace=True)
    atac = pr.PyRanges(df2)

    df3 = pyreadr.read_r('DNAseq_dt.RDS')[None]
    df3.rename(columns={"start": "Start", "end": "End", "seqnames": "Chromosome"}, inplace=True)
    dnaseq = pr.PyRanges(df3)

    m_lengths = merged_ranges.lengths()
    merged_ranges = merged_ranges.assign('Length', lambda x : m_lengths)

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

    gencode = pr.read_gtf("gencode_chr18.gtf")

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

    return evaluation

if __name__ == "__main__":
    print(evaluate())