import pyreadr
import pandas as pd
import pyranges as pr
import numpy as np
import config

SIGNAL_COLUMN = config.SIGNAL_COLUMN

def calculate_norm_signal(ranges, targets_df):
    ranges = ranges.sort()
    targets_df['target_count'] = 1
    

    from_array_open = list(ranges.df['Start'])
    to_array_open = list(ranges.df['End'])

    bins = pd.IntervalIndex.from_arrays(from_array_open, to_array_open, closed="left")

    targets_df['ranges_bin'] = pd.cut(targets_df['start'], bins=bins)
    summed_ranges = targets_df.groupby('ranges_bin')[[SIGNAL_COLUMN, 'target_count']].sum()
    summed_ranges["norm_signal"] = np.divide(summed_ranges[SIGNAL_COLUMN], summed_ranges['target_count'], where=summed_ranges['target_count'] != 0, out=np.zeros_like(summed_ranges[SIGNAL_COLUMN], dtype=float))

    """     targets_df['target_count'] = 1
    targets = pr.PyRanges(targets_df.rename(columns={"start": "Start", "end": "End", "seqnames": "Chromosome"}))
    ranges_targets = ranges.join(targets, how='left')
    print(ranges_targets)
    summed_ranges = ranges_targets.df.groupby("Start")[[SIGNAL_COLUMN, 'target_count']].sum() """

    return summed_ranges

def get_pre_regions(ranges):
    pre_regions = ranges.copy()
    pre_regions.Start = ranges.Start - 100
    pre_regions.End = ranges.Start

    return pre_regions

def get_post_regions(ranges):
    post_regions = ranges.copy()
    post_regions.Start = ranges.End
    post_regions.End = ranges.End + 100

    return post_regions

def calculate_enrichment_factors(ranges, targets_df):
    ranges_signal = calculate_norm_signal(ranges, targets_df)
    
    pre_signal = calculate_norm_signal(get_pre_regions(ranges), targets_df)
    pre_signal = pre_signal.rename(columns={'norm_signal': 'pre_signal'})
    
    post_signal = calculate_norm_signal(get_post_regions(ranges), targets_df)
    post_signal = post_signal.rename(columns={'norm_signal': 'post_signal'})

    ranges_signal = pd.concat([ranges_signal.reset_index(), pre_signal.reset_index(drop=True)], axis=1)
    ranges_signal = pd.concat([ranges_signal.reset_index(), post_signal.reset_index(drop=True)], axis=1)
    ranges_signal = ranges_signal.drop(columns=[SIGNAL_COLUMN, 'target_count', 'index'], axis=1)

    ranges_signal['surround_signal'] = ranges_signal['pre_signal'] + ranges_signal['post_signal']
    ranges_signal['enrichment_factor1'] = np.divide(ranges_signal['norm_signal'], ranges_signal['surround_signal'], where=ranges_signal['surround_signal'] != 0, out=np.zeros_like(ranges_signal['norm_signal'], dtype=float))

    general_norm_signal = caculate_dataframe_norm_signal(targets_df)
    print("Dataframe normalized signal:", general_norm_signal)

    if general_norm_signal == 0:
        ranges_signal['enrichment_factor2'] = 0
    else:
        ranges_signal['enrichment_factor2'] = ranges_signal['norm_signal'] / general_norm_signal

    return ranges_signal

def caculate_dataframe_norm_signal(targets_df):
    return targets_df[SIGNAL_COLUMN].sum() / len(targets_df)

def calculate_enrichment_factors_file(ranges_filepath, targets_filepath):
    targets_df = pyreadr.read_r(targets_filepath)[None]
    ranges = pr.PyRanges(pd.read_csv(ranges_filepath, sep='\t'))

    return calculate_enrichment_factors(ranges, targets_df)