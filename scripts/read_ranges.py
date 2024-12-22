import pyranges as pr
import pandas as pd
import numpy as np

#pd.set_option('display.min_rows', 500)
#pd.set_option('display.max_rows', 500)

def get_ranges():
    gr = pr.read_gtf("../inputs/gencode_chr18.gtf")

    # {gene,transcript,exon,CDS,UTR,start_codon,stop_codon,Selenocysteine}
    genes = gr[gr.Feature == 'gene']
    transcripts = gr[gr.Feature == 'transcript']
    exons = gr[gr.Feature == 'exon']
    CDSs = gr[gr.Feature == 'CDS']

    grs = {'genes': genes, 'transcripts': transcripts, 'exons': exons, 'CDSs': CDSs}

    # HARDCODED
    chromo_start = 3000016
    chromo_end = 9717593
    offset = 0

    breaks = list(np.arange(chromo_start + offset, chromo_end, 100))
    from_array = [chromo_start] + breaks
    to_array = breaks + [chromo_end]

    bins = pr.PyRanges(chromosomes="18", starts=from_array, ends=to_array)

    overlaps_result = pr.count_overlaps(grs, bins)
    overlaps_df = overlaps_result.df

    # print(overlaps_df)
    return overlaps_df

if __name__ == "__main__":
    get_ranges()