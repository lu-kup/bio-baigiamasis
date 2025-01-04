import numpy as np
import pyreadr
import pandas as pd
import pyranges as pr
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from kmodes.kprototypes import KPrototypes
import gower
import read_ranges
import math

BIN_SIZE = 100
SAMPLE_FILEPATH = '../inputs/subset2.rds'

SIGNAL_COLUMN = 'TT_S1'
OTHER_SIGNALS = ['TT_S0', 'TT_S2']

def add_bins(offset, chromosome):
    chromo_start = chromosome["start"].min() - 1
    chromo_end = chromosome["start"].max()

    breaks = list(np.arange(chromo_start + offset, chromo_end, BIN_SIZE))
    from_array = [chromo_start] + breaks
    to_array = breaks + [chromo_end]

    bins = pd.IntervalIndex.from_arrays(from_array, to_array)

    print("BINS ", bins)

    chromosome['bin_offset_' + str(offset)] = pd.cut(chromosome['start'], bins=bins)

    return from_array, to_array

def algo1d(bin_offset = 0):
    df = pyreadr.read_r(SAMPLE_FILEPATH)[None]

    print("HEAD")
    print(df.head())

    add_bins(bin_offset, df)

    print("HEAD with BINS")
    print(df.head())
    print("TAIL with BINS")
    print(df.tail())
    print()

    df.drop(['seqnames', 'start', 'end', 'CG_ID'], axis=1, inplace=True)
    sumos = df.groupby('bin_offset_' + str(bin_offset)).sum()

    sumos.drop(columns=OTHER_SIGNALS, axis=1, inplace=True)

    print("SUMOS")
    print(sumos[:50])
    print(sumos.shape)
    print(sumos.dtypes)

    # Skaičiuojam
    kmeans = KMeans(n_clusters=2)
    kmeans.fit(sumos)
    inertia = kmeans.inertia_
    labels = list(kmeans.labels_)
    sumos['labels'] = labels
    open_ratio = labels.count(1)/len(labels)

    print(open_ratio)
    print(kmeans.labels_[:1000])

    np.savetxt('../outputs/labels.dat', kmeans.labels_, fmt='%d', delimiter=',')

    return sumos

def algo2d(scale = 1, bin_offset = 0):
    df = pyreadr.read_r(SAMPLE_FILEPATH)[None]

    print("HEAD")
    print(df.head())

    add_bins(bin_offset, df)

    print("HEAD with BINS")
    print(df.head())
    print("TAIL with BINS")
    print(df.tail())
    print()

    df.drop(['seqnames', 'start', 'end', 'CG_ID'], axis=1, inplace=True)
    sumos = df.groupby('bin_offset_' + str(bin_offset)).sum().reset_index()
    sumos['starting_nt'] = sumos['bin_offset_' + str(bin_offset)].apply(lambda x: x.left)

    sumos.drop(['bin_offset_' + str(bin_offset)], axis=1, inplace=True)
    sumos.drop(columns=OTHER_SIGNALS, axis=1, inplace=True)

    print("SUMOS")
    print(sumos[:50])
    print(sumos.shape)
    print(sumos.dtypes)

    sumos['signal'] = sumos[SIGNAL_COLUMN] * scale
    sumos.drop([SIGNAL_COLUMN], axis=1, inplace=True)

    print("SUMOS NEW")
    print(sumos[:50])

    # Skaičiuojam
    kmeans = KMeans(n_clusters=2)
    kmeans.fit(sumos)
    inertia = kmeans.inertia_
    labels = list(kmeans.labels_)
    sumos['labels'] = labels
    open_ratio = labels.count(1)/len(labels)

    print(open_ratio)
    print(kmeans.labels_[:1000])

    sumos.to_csv('../outputs/output_2d.csv', sep = '\t')

    return sumos

def algo5d(scale = 1, n_clusters = 2, bin_offset = 0):
    df = pyreadr.read_r(SAMPLE_FILEPATH)[None]

    print("HEAD")
    print(df.head())

    from_array, to_array = add_bins(bin_offset, df)
    df['target_count'] = 1

    print("HEAD with BINS")
    print(df.head())
    print("TAIL with BINS")
    print(df.tail())
    print()

    df.drop(['seqnames', 'start', 'end', 'CG_ID'], axis=1, inplace=True)
    sumos = df.groupby('bin_offset_' + str(bin_offset)).sum().reset_index()
    sumos['starting_nt'] = sumos['bin_offset_' + str(bin_offset)].apply(lambda x: x.left)

    sumos.drop(columns=OTHER_SIGNALS, axis=1, inplace=True)

    print("SUMOS")
    print(sumos[:50])
    print(sumos.shape)
    print(sumos.dtypes)

    ranges = read_ranges.get_ranges(from_array, to_array)
    print("sumos ilgis:", len(sumos))
    print("ranges ilgis:", len(ranges))

    sumos_features = pd.concat([sumos.reset_index(drop=True), ranges.reset_index(drop=True)], axis=1)
    model_input = sumos_features.copy()
    model_input.drop(['target_count', 'starting_nt', 'Chromosome', 'Start', 'End', 'bin_offset_' + str(bin_offset)], axis=1, inplace=True)
    print("Model input - NO SCALE")
    print(model_input)

    model_input['genes_scaled'] = model_input['genes'] * scale
    model_input.drop(['genes'], axis=1, inplace=True)

    model_input['transcripts_scaled'] = model_input['transcripts'] * scale
    model_input.drop(['transcripts'], axis=1, inplace=True)

    model_input['exons_scaled'] = model_input['exons'] * scale
    model_input.drop(['exons'], axis=1, inplace=True)

    model_input['CDSs_scaled'] = model_input['CDSs'] * scale
    model_input.drop(['CDSs'], axis=1, inplace=True)

    print("Model input - SCALED")
    print(model_input)

    # Skaičiuojam
    kmeans = KMeans(n_clusters=n_clusters)
    kmeans.fit(model_input)
    inertia = kmeans.inertia_
    labels = list(kmeans.labels_)
    sumos_features['labels'] = kmeans.labels_
    open_ratio = labels.count(1)/len(labels)

    print(open_ratio)
    print(kmeans.labels_[:1000])

    sumos_features.to_csv('../outputs/output_5d.csv', sep = '\t')

    return sumos_features

def algo_prototypes(gamma = 1, n_clusters = 2, bin_offset = 0):
    df = pyreadr.read_r(SAMPLE_FILEPATH)[None]

    print("HEAD")
    print(df.head())

    from_array, to_array = add_bins(bin_offset, df)
    df['target_count'] = 1

    print("HEAD with BINS")
    print(df.head())
    print("TAIL with BINS")
    print(df.tail())
    print()

    df.drop(['seqnames', 'start', 'end', 'CG_ID'], axis=1, inplace=True)
    sumos = df.groupby('bin_offset_' + str(bin_offset)).sum().reset_index()
    sumos['starting_nt'] = sumos['bin_offset_' + str(bin_offset)].apply(lambda x: x.left)

    sumos.drop(columns=OTHER_SIGNALS, axis=1, inplace=True)

    print("SUMOS")
    print(sumos[:50])
    print(sumos.shape)
    print(sumos.dtypes)

    ranges = read_ranges.get_ranges(from_array, to_array)
    print("sumos ilgis:", len(sumos))
    print("ranges ilgis:", len(ranges))

    sumos_features = pd.concat([sumos.reset_index(drop=True), ranges.reset_index(drop=True)], axis=1)
    model_input = sumos_features.copy()
    model_input.drop(['target_count', 'starting_nt', 'Chromosome', 'Start', 'End', 'bin_offset_' + str(bin_offset)], axis=1, inplace=True)
    print(model_input)

    # Skaičiuojam
    kprotos = KPrototypes(n_clusters=n_clusters, gamma=gamma, verbose=1)
    kprotos.fit(model_input, categorical=[1, 2, 3, 4])
    labels = list(kprotos.labels_)
    sumos_features['labels'] = kprotos.labels_
    open_ratio = labels.count(1)/len(labels)

    print(open_ratio)
    print(kprotos.labels_[:1000])

    sumos_features.to_csv('../outputs/output_prototypes.csv', sep = '\t')

    return sumos_features

def algo_dbscan(eps = 0.01, min_samples = 6, bin_offset = 0):
    BATCH_SIZE = 10000

    df = pyreadr.read_r(SAMPLE_FILEPATH)[None]

    print("HEAD")
    print(df.head())

    from_array, to_array = add_bins(bin_offset, df)
    df['target_count'] = 1

    print("HEAD with BINS")
    print(df.head())
    print("TAIL with BINS")
    print(df.tail())
    print()

    df.drop(['seqnames', 'start', 'end', 'CG_ID'], axis=1, inplace=True)
    sumos = df.groupby('bin_offset_' + str(bin_offset)).sum().reset_index()
    sumos['starting_nt'] = sumos['bin_offset_' + str(bin_offset)].apply(lambda x: x.left)

    sumos.drop(columns=OTHER_SIGNALS, axis=1, inplace=True)

    print("SUMOS")
    print(sumos[:50])
    print(sumos.shape)
    print(sumos.dtypes)

    ranges = read_ranges.get_ranges(from_array, to_array)
    print("sumos ilgis:", len(sumos))
    print("ranges ilgis:", len(ranges))

    sumos_features = pd.concat([sumos.reset_index(drop=True), ranges.reset_index(drop=True)], axis=1)

    start_index = 0
    end_index = BATCH_SIZE
    total_values = len(sumos_features)
    batch_no = math.ceil(total_values / BATCH_SIZE)
    labels = pd.DataFrame(index = range(total_values), columns = ['labels'])

    for i in range(batch_no):
        model_input = sumos_features.copy()[start_index:end_index]
        model_input.drop(['target_count', 'starting_nt', 'Chromosome', 'Start', 'End', 'bin_offset_' + str(bin_offset)], axis=1, inplace=True)
        print(model_input)

        # Skaičiuojam
        print("Length of input:", len(model_input))
        print("Calculating distances...")
        distance_matrix = gower.gower_matrix(model_input)
        clustering = DBSCAN(eps=eps, min_samples=min_samples, metric = "precomputed")

        print("Running DBSCAN...")
        clustering.fit(distance_matrix)

        labels.loc[start_index:end_index - 1, 'labels'] = clustering.labels_

        start_index += BATCH_SIZE
        end_index += BATCH_SIZE

        if (end_index >= total_values):
            end_index = total_values - 1

    sumos_features['labels'] = labels
    final_labels = list(sumos_features['labels'])
    open_ratio = final_labels.count(1)/total_values

    print(open_ratio)
    print("Number of clusters:", max(sumos_features['labels']) + 1)

    sumos_features.to_csv('../outputs/output_dbscan.csv', sep = '\t')

    return sumos_features

def algo_dbscan_aggregated(eps = 0.01, min_samples = 6, threshold = 1.5, bin_offset = 0):
    BATCH_SIZE = 10000

    df = pyreadr.read_r(SAMPLE_FILEPATH)[None]

    print("HEAD")
    print(df.head())

    from_array, to_array = add_bins(bin_offset, df)
    df['target_count'] = 1

    print("HEAD with BINS")
    print(df.head())
    print("TAIL with BINS")
    print(df.tail())
    print()

    df.drop(['seqnames', 'start', 'end', 'CG_ID'], axis=1, inplace=True)
    sumos = df.groupby('bin_offset_' + str(bin_offset)).sum().reset_index()
    sumos['starting_nt'] = sumos['bin_offset_' + str(bin_offset)].apply(lambda x: x.left)

    sumos.drop(columns=OTHER_SIGNALS, axis=1, inplace=True)

    print("SUMOS")
    print(sumos[:50])
    print(sumos.shape)
    print(sumos.dtypes)

    ranges = read_ranges.get_ranges(from_array, to_array)
    print("sumos ilgis:", len(sumos))
    print("ranges ilgis:", len(ranges))

    sumos_features = pd.concat([sumos.reset_index(drop=True), ranges.reset_index(drop=True)], axis=1)

    start_index = 0
    end_index = BATCH_SIZE
    total_values = len(sumos_features)
    batch_no = math.ceil(total_values / BATCH_SIZE)
    labels = pd.DataFrame(index = range(total_values), columns = ['labels'])
    signal_density = get_binned_signal_density(sumos_features)

    for i in range(batch_no):
        batch = sumos_features.copy()[start_index:end_index]
        model_input = batch.copy()
        model_input.drop(['target_count', 'starting_nt', 'Chromosome', 'Start', 'End', 'bin_offset_' + str(bin_offset)], axis=1, inplace=True)
        print(model_input)

        # Skaičiuojam
        print("Length of input:", len(model_input))
        print("Calculating distances...")
        distance_matrix = gower.gower_matrix(model_input, cat_features=[False, True, True, True, True])
        clustering = DBSCAN(eps=eps, min_samples=min_samples, metric = "precomputed")

        print("Running DBSCAN...")
        clustering.fit(distance_matrix)

        batch['labels'] = clustering.labels_
        update_labels(batch, threshold, signal_density)
        labels.loc[start_index:end_index - 1, 'labels'] = batch['labels']

        start_index += BATCH_SIZE
        end_index += BATCH_SIZE

        if (end_index >= total_values):
            end_index = total_values - 1

    sumos_features['labels'] = labels
    final_labels = list(sumos_features['labels'])
    open_ratio = final_labels.count(1)/total_values

    print(open_ratio)
    print("Number of clusters:", max(sumos_features['labels']) + 1)

    sumos_features.to_csv('../outputs/output_dbscan_aggregated.csv', sep = '\t')

    return bins_to_ranges(sumos_features, 'DBSCAN')

def algo5d_aggregated(scale = 1, n_clusters = 2, threshold = 1.5, bin_offset = 0):
    model_output = algo5d(scale, n_clusters, bin_offset)
    update_labels(model_output, threshold)
    model_output.to_csv('../outputs/output_algo5d_aggregated.csv', sep = '\t')
    return bins_to_ranges(model_output, 'k-means 5D')

def algo_prototypes_aggregated(gamma = 1, n_clusters = 2, threshold = 1.5, bin_offset = 0):
    model_output = algo_prototypes(gamma, n_clusters, bin_offset)
    update_labels(model_output, threshold)
    model_output.to_csv('../outputs/output_prototypes_aggregated.csv', sep = '\t')
    return bins_to_ranges(model_output, 'k-prototypes')

def map_clusters(dataframe, threshold, signal_density = None):
    cluster_labels = dataframe['labels'].unique()
    if signal_density is None:
        signal_density = get_binned_signal_density(dataframe)
    threshold_density = threshold * signal_density
    print("threshold density", threshold_density)
    mapping = {}

    for cluster_label in cluster_labels:
        cluster = dataframe[dataframe['labels'] == cluster_label]
        if len(cluster) == 0:
            continue
        cluster_signal = cluster[SIGNAL_COLUMN].sum()
        cluster_target_count = cluster['target_count'].sum()
        print("cluster signal", cluster_signal)
        print("cluster target count", cluster_target_count)

        if cluster_target_count != 0:
            cluster_density = cluster_signal / cluster_target_count
        else:
            cluster_density = 0

        print("!!cluster density", cluster_density)
        if cluster_density > threshold_density:
            mapping[cluster_label] = 1
        else:
            mapping[cluster_label] = 0

    return mapping 

def update_labels(dataframe, threshold = 1.5, signal_density = None):
    mapping = map_clusters(dataframe, threshold, signal_density)
    print("\nLabels mapping\n", mapping)
    dataframe['labels'] = dataframe['labels'].map(mapping)

def evaluate_bins(df, model_name):
    evaluation = {}
    evaluation['model_name'] = model_name
    evaluation['min_open_signal'] = df.loc[df["labels"] == open_label, SIGNAL_COLUMN].min()
    evaluation['max_open_signal'] = df.loc[df["labels"] == open_label, SIGNAL_COLUMN].max()
    evaluation['min_closed_signal'] = df.loc[df["labels"] != open_label, SIGNAL_COLUMN].min()
    evaluation['max_closed_signal'] = df.loc[df["labels"] != open_label, SIGNAL_COLUMN].max()
    evaluation['no_clusters'] = max(df['labels']) + 1

    return evaluation

def bins_to_ranges(df, model_name):
    input = pr.PyRanges(df)

    count0 = len(input[input.labels == 0])
    count1 = len(input[input.labels == 1])

    if count0 > count1:
        open_label = 1
    else:
        open_label = 0

    print("DEBUG:", count0, count1, "\n\n", input[input.labels == 0])

    open_ranges = input[input.labels == open_label]
    merged_ranges = open_ranges.merge()

    return merged_ranges

def get_binned_signal_density(dataframe):
    if (dataframe['target_count'].sum() != 0):
        return dataframe[SIGNAL_COLUMN].sum() / dataframe['target_count'].sum()
    else:
        print("ERROR: Zero signal density in dataframe")
        return float('inf')

if __name__ == "__main__":
    import evaluator
    out = algo_dbscan_aggregated(eps = 0.001, threshold = 1.5)
    #out = pd.read_csv('../outputs/output_dbscan.csv', sep = '\t', index_col=0)
    print(evaluator.evaluate(dataframe = out))
