#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pyranges as pr
import evaluator
import chromatin_algos

def generate_consensus_5x20(algorithm, params):
    open_ranges_dict = {}

    for i in range(20):
        current_offset = i * 5
        output_ranges = algorithm(*params, bin_offset = current_offset)
        open_ranges_dict['offset_' + str(current_offset)] = output_ranges

    overlaps = pr.count_overlaps(open_ranges_dict)
    print("OVERLAPS", overlaps)

    overlaps_df = overlaps.df
    filtered = overlaps_df.filter(regex = "offset_")
    overlaps_df['predict_state_SUM'] = filtered.sum(axis = 1)
    overlaps_df['labels'] = overlaps_df["predict_state_SUM"] > 10
    overlaps_df['labels'] = overlaps_df['labels'].astype(int)
    overlaps_filtered = pr.PyRanges(overlaps_df[overlaps_df['labels'] == 1]).merge()
    print("OVERLAPS FILTERED", overlaps_filtered)

    overlaps_filtered.df.to_csv('../outputs/output_consensus.csv', sep = '\t')
    return overlaps_filtered

if __name__ == "__main__":
    out = generate_consensus_5x20(chromatin_algos.algo5d_aggregated, [1, 2, 1.5])
    print(evaluator.evaluate(merged_ranges = out))
