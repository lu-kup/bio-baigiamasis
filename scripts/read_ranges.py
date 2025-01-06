import pyranges as pr
import pandas as pd
import numpy as np
import config

BIN_SIZE = 100
CROMOSOME_NO = "18"

gencode_data_path = config.GENCODE_DATA_PATH

#pd.set_option('display.min_rows', 500)
#pd.set_option('display.max_rows', 500)

def get_ranges(from_array, to_array):
    gr = pr.read_gtf(gencode_data_path)

    # {gene,transcript,exon,CDS,UTR,start_codon,stop_codon,Selenocysteine}
    genes = gr[gr.Feature == 'gene']
    transcripts = gr[gr.Feature == 'transcript']
    exons = gr[gr.Feature == 'exon']
    CDSs = gr[gr.Feature == 'CDS']

    grs = {'genes': genes, 'transcripts': transcripts, 'exons': exons, 'CDSs': CDSs}
    bins = pr.PyRanges(chromosomes=CROMOSOME_NO, starts=from_array, ends=to_array)

    overlaps_result = pr.count_overlaps(grs, bins)
    overlaps_df = overlaps_result.df

    # print(overlaps_df)
    return overlaps_df

if __name__ == "__main__":
    get_ranges()