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
SIGNAL_COLUMN = 'TT_S0'

def add_bins(offset, chromosome):
    chromo_start = chromosome["start"].min() - 1
    chromo_end = chromosome["start"].max()

    breaks = list(np.arange(chromo_start + offset, chromo_end, BIN_SIZE))
    from_array = [chromo_start] + breaks
    to_array = breaks + [chromo_end]

    bins = pd.IntervalIndex.from_arrays(from_array, to_array)

    print("BINS ", bins)

    chromosome['bin_offset_' + str(offset)] = pd.cut(chromosome['start'], bins=bins)

    return bins

def algo1d(bin_offset = 0):
    result = pyreadr.read_r('../inputs/subset1.rds')
    df = result[None]    

    print("HEAD")
    print(df.head())

    add_bins(bin_offset, df)

    print("HEAD w BINS")
    print(df.head())
    print("TAIL w BINS")
    print(df.tail())
    print()

    df.drop(['seqnames', 'start', 'end', 'CG_ID'], axis=1, inplace=True)
    sumos = df.groupby('bin_offset_' + str(bin_offset)).sum()

    sumos.drop(['TT_S1', 'TT_S2'], axis=1, inplace=True)

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

    a = 10
    counter = 0
    for labl in labels:
        b = a
        a = labl
        if a != b:
            counter += 1
    print(counter - 1)

    print(open_ratio)
    print(kmeans.labels_[:1000])

    np.savetxt('../outputs/labels.dat', kmeans.labels_, fmt='%d', delimiter=',')

    return sumos

def algo2d(scale = 1, bin_offset = 0):
    df = pyreadr.read_r('../inputs/subset1.rds')[None]

    print("HEAD")
    print(df.head())

    add_bins(bin_offset, df)

    print("HEAD w BINS")
    print(df.head())
    print("TAIL w BINS")
    print(df.tail())
    print()

    df.drop(['seqnames', 'start', 'end', 'CG_ID'], axis=1, inplace=True)
    sumos = df.groupby('bin_offset_' + str(bin_offset)).sum().reset_index()
    sumos['starting_nt'] = sumos['bin_offset_' + str(bin_offset)].apply(lambda x: x.left)

    sumos.drop(['TT_S1', 'TT_S2', 'bin_offset_' + str(bin_offset)], axis=1, inplace=True)

    print("SUMOS")
    print(sumos[:50])
    print(sumos.shape)
    print(sumos.dtypes)

    sumos['signal'] = sumos['TT_S0'] * scale
    sumos.drop(['TT_S0'], axis=1, inplace=True)

    print("SUMOS NEW")
    print(sumos[:50])

    # Skaičiuojam
    kmeans = KMeans(n_clusters=2)
    kmeans.fit(sumos)
    inertia = kmeans.inertia_
    labels = list(kmeans.labels_)
    sumos['labels'] = labels
    open_ratio = labels.count(1)/len(labels)

    a = 10
    counter = 0
    for labl in labels:
        b = a
        a = labl
        if a != b:
            counter += 1
    print(counter - 1)

    print(open_ratio)
    print(kmeans.labels_[:1000])

    sumos.to_csv('../outputs/output_2d.csv', sep = '\t')

    return sumos

def algo5d(scale = 1, n_clusters = 2, bin_offset = 0):
    result = pyreadr.read_r('../inputs/subset1.rds')
    df = result[None]

    print("HEAD")
    print(df.head())

    add_bins(bin_offset, df)

    print("HEAD w BINS")
    print(df.head())
    print("TAIL w BINS")
    print(df.tail())
    print()

    df.drop(['seqnames', 'start', 'end', 'CG_ID'], axis=1, inplace=True)
    sumos = df.groupby('bin_offset_' + str(bin_offset)).sum().reset_index()
    sumos['starting_nt'] = sumos['bin_offset_' + str(bin_offset)].apply(lambda x: x.left)

    sumos.drop(['TT_S1', 'TT_S2'], axis=1, inplace=True)

    print("SUMOS")
    print(sumos[:50])
    print(sumos.shape)
    print(sumos.dtypes)

    ranges = read_ranges.get_ranges(bin_offset)
    print("sumos ilgis:", len(sumos))
    print("ranges ilgis:", len(ranges))

    sumos_features = pd.concat([sumos.reset_index(drop=True), ranges.reset_index(drop=True)], axis=1)
    model_input = sumos_features.copy()
    model_input.drop(['starting_nt', 'Chromosome', 'Start', 'End', 'bin_offset_' + str(bin_offset)], axis=1, inplace=True)
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
    result = pyreadr.read_r('../inputs/subset1.rds')
    df = result[None]

    print("HEAD")
    print(df.head())

    add_bins(bin_offset, df)

    print("HEAD w BINS")
    print(df.head())
    print("TAIL w BINS")
    print(df.tail())
    print()

    df.drop(['seqnames', 'start', 'end', 'CG_ID'], axis=1, inplace=True)
    sumos = df.groupby('bin_offset_' + str(bin_offset)).sum().reset_index()
    sumos['starting_nt'] = sumos['bin_offset_' + str(bin_offset)].apply(lambda x: x.left)

    sumos.drop(['TT_S1', 'TT_S2'], axis=1, inplace=True)

    print("SUMOS")
    print(sumos[:50])
    print(sumos.shape)
    print(sumos.dtypes)

    ranges = read_ranges.get_ranges(bin_offset)
    print("sumos ilgis:", len(sumos))
    print("ranges ilgis:", len(ranges))

    sumos_features = pd.concat([sumos.reset_index(drop=True), ranges.reset_index(drop=True)], axis=1)
    model_input = sumos_features.copy()
    model_input.drop(['starting_nt', 'Chromosome', 'Start', 'End', 'bin_offset_' + str(bin_offset)], axis=1, inplace=True)
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

    result = pyreadr.read_r('../inputs/subset1.rds')
    df = result[None]

    print("HEAD")
    print(df.head())

    add_bins(bin_offset, df)

    print("HEAD w BINS")
    print(df.head())
    print("TAIL w BINS")
    print(df.tail())
    print()

    df.drop(['seqnames', 'start', 'end', 'CG_ID'], axis=1, inplace=True)
    sumos = df.groupby('bin_offset_' + str(bin_offset)).sum().reset_index()
    sumos['starting_nt'] = sumos['bin_offset_' + str(bin_offset)].apply(lambda x: x.left)

    sumos.drop(['TT_S1', 'TT_S2'], axis=1, inplace=True)

    print("SUMOS")
    print(sumos[:50])
    print(sumos.shape)
    print(sumos.dtypes)

    ranges = read_ranges.get_ranges(bin_offset)
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
        model_input.drop(['starting_nt', 'Chromosome', 'Start', 'End', 'bin_offset_' + str(bin_offset)], axis=1, inplace=True)
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

    result = pyreadr.read_r('../inputs/subset1.rds')
    df = result[None]

    print("HEAD")
    print(df.head())

    add_bins(bin_offset, df)

    print("HEAD w BINS")
    print(df.head())
    print("TAIL w BINS")
    print(df.tail())
    print()

    df.drop(['seqnames', 'start', 'end', 'CG_ID'], axis=1, inplace=True)
    sumos = df.groupby('bin_offset_' + str(bin_offset)).sum().reset_index()
    sumos['starting_nt'] = sumos['bin_offset_' + str(bin_offset)].apply(lambda x: x.left)

    sumos.drop(['TT_S1', 'TT_S2'], axis=1, inplace=True)

    print("SUMOS")
    print(sumos[:50])
    print(sumos.shape)
    print(sumos.dtypes)

    ranges = read_ranges.get_ranges(bin_offset)
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
        model_input = sumos_features.copy()[start_index:end_index]
        model_input.drop(['starting_nt', 'Chromosome', 'Start', 'End', 'bin_offset_' + str(bin_offset)], axis=1, inplace=True)
        print(model_input)

        # Skaičiuojam
        print("Length of input:", len(model_input))
        print("Calculating distances...")
        distance_matrix = gower.gower_matrix(model_input)
        clustering = DBSCAN(eps=eps, min_samples=min_samples, metric = "precomputed")

        print("Running DBSCAN...")
        clustering.fit(distance_matrix)

        model_input['labels'] = clustering.labels_
        update_labels(model_input, threshold, signal_density)
        labels.loc[start_index:end_index - 1, 'labels'] = model_input['labels']

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

    return sumos_features

def algo5d_aggregated(scale = 1, n_clusters = 2, threshold = 1.5, bin_offset = 0):
    model_output = algo5d(scale, n_clusters, bin_offset)
    update_labels(model_output, threshold)
    model_output.to_csv('../outputs/output_algo5d_aggregated.csv', sep = '\t')
    return model_output

def algo_prototypes_aggregated(gamma = 1, n_clusters = 2, threshold = 1.5, bin_offset = 0):
    model_output = algo_prototypes(gamma, n_clusters, bin_offset)
    update_labels(model_output, threshold)
    model_output.to_csv('../outputs/output_prototypes_aggregated.csv', sep = '\t')
    return model_output

def map_clusters(dataframe, threshold, signal_density = None):
    cluster_labels = dataframe['labels'].unique()
    if signal_density is None:
        signal_density = get_binned_signal_density(dataframe)
    threshold_density = threshold * signal_density
    print("threshold density", threshold_density)
    mapping = {}

    for cluster_label in cluster_labels:
        cluster = dataframe[dataframe['labels'] == cluster_label]
        cluster_length = len(cluster) * BIN_SIZE
        if cluster_length == 0:
            continue
        cluster_signal = cluster[SIGNAL_COLUMN].sum()
        print("cluster signal", cluster_signal)
        print("cluster length", cluster_length)
        cluster_density = cluster_signal / cluster_length
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

def get_binned_signal_density(dataframe):
    return dataframe[SIGNAL_COLUMN].sum() / (dataframe['starting_nt'].max() + BIN_SIZE - dataframe['starting_nt'].min())

def generate_consensus_5x20(algorithm, params):
    dataframe = pyreadr.read_r('../inputs/subset1.rds')[None]

    for i in range(20):
        current_offset = i * 5
        output_dataframe = algorithm(*params, bin_offset = current_offset)
        output_dataframe = output_dataframe.rename(columns = {'labels' : 'predicted_state_' + str(current_offset)})
        states_df = output_dataframe[['bin_offset_' + str(current_offset), 'predicted_state_' + str(current_offset)]]
        dataframe = dataframe.join(states_df, on = 'bin_offset_' + str(current_offset))

    filtered = dataframe.filter(regex="predicted_state_")
    dataframe['predict_state_SUM'] = filtered.sum(axis=1)
    dataframe['labels'] = dataframe["predict_state_SUM"] > 10
    dataframe['labels'] = dataframe['labels'].astype(int)

    merged_output = pr.PyRanges(dataframe).merge(by = 'labels')
    merged_output.df.to_csv('../outputs/output_consensus.csv', sep = '\t')

if __name__ == "__main__":
    import evaluator
    out = algo_dbscan_aggregated(eps = 0.001, threshold = 1.5)
    #out = pd.read_csv('../outputs/output_dbscan.csv', sep = '\t', index_col=0)
    print(evaluator.evaluate(dataframe = out))
