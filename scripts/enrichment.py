import pyreadr
import pandas as pd
import pyranges as pr
import numpy as np
import config

sample_filepath = config.SAMPLE_FILEPATH

def calculate_norm_signal(ranges, targets_df, signal_column):
    ranges = ranges.sort()
    targets_df['target_count'] = 1

    from_array_open = list(ranges.df['Start'])
    to_array_open = list(ranges.df['End'])

    bins = pd.IntervalIndex.from_arrays(from_array_open, to_array_open, closed="left")

    targets_df['ranges_bin'] = pd.cut(targets_df['start'], bins=bins)
    summed_ranges = targets_df.groupby('ranges_bin')[[signal_column, 'target_count']].sum()
    summed_ranges["norm_signal"] = np.divide(summed_ranges[signal_column], summed_ranges['target_count'], where=summed_ranges['target_count'] != 0, out=np.zeros_like(summed_ranges[signal_column], dtype=float))

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

def calculate_enrichment_factors(ranges, signal_column):
    targets_df = pyreadr.read_r(sample_filepath)[None]
    ranges_signal = calculate_norm_signal(ranges, targets_df, signal_column)
    
    pre_signal = calculate_norm_signal(get_pre_regions(ranges), targets_df, signal_column)
    pre_signal = pre_signal.rename(columns={'norm_signal': 'pre_signal'})
    
    post_signal = calculate_norm_signal(get_post_regions(ranges), targets_df, signal_column)
    post_signal = post_signal.rename(columns={'norm_signal': 'post_signal'})

    ranges_signal = pd.concat([ranges_signal.reset_index(), pre_signal.reset_index(drop=True)], axis=1)
    ranges_signal = pd.concat([ranges_signal.reset_index(), post_signal.reset_index(drop=True)], axis=1)
    ranges_signal = ranges_signal.drop(columns=[signal_column, 'target_count', 'index'], axis=1)

    ranges_signal['surround_signal'] = ranges_signal['pre_signal'] + ranges_signal['post_signal']
    ranges_signal['enrichment_factor1'] = np.divide(ranges_signal['norm_signal'], ranges_signal['surround_signal'], where=ranges_signal['surround_signal'] != 0, out=np.zeros_like(ranges_signal['norm_signal'], dtype=float))

    general_norm_signal = caculate_dataframe_norm_signal(targets_df, signal_column)
    print("Dataframe normalized signal:", general_norm_signal)

    if general_norm_signal == 0:
        ranges_signal['enrichment_factor2'] = 0
    else:
        ranges_signal['enrichment_factor2'] = ranges_signal['norm_signal'] / general_norm_signal

    return ranges_signal

def caculate_dataframe_norm_signal(targets_df, signal_column):
    return targets_df[signal_column].sum() / len(targets_df)

def calculate_enrichment_factors_file(ranges_filepath, targets_filepath):
    targets_df = pyreadr.read_r(targets_filepath)[None]
    ranges = pr.PyRanges(pd.read_csv(ranges_filepath, sep='\t'))

    return calculate_enrichment_factors(ranges, targets_df)