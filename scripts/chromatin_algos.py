import numpy as np
import pyreadr
import pandas as pd
import pyranges as pr
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from kmodes.kprototypes import KPrototypes
import gower
import read_ranges

def add_bins(offset, chromosome):
    chromo_start = chromosome["start"].min() - 1
    chromo_end = chromosome["start"].max()

    breaks = list(np.arange(chromo_start + offset, chromo_end, 100))
    from_array = [chromo_start] + breaks
    to_array = breaks + [chromo_end]

    bins = pd.IntervalIndex.from_arrays(from_array, to_array)

    print("BINS ", bins)

    chromosome['bin_offset_' + str(offset)] = pd.cut(chromosome['start'], bins=bins)

    return bins

def algo1d():
    result = pyreadr.read_r('../inputs/subset1.rds')
    df = result[None]
    bin_offset = 0

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

def algo2d(scale = 1):
    df = pyreadr.read_r('../inputs/subset1.rds')[None]
    bin_offset = 0

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
    sumos['starting_nt'] = sumos['bin_offset_0'].apply(lambda x: x.left)

    sumos.drop(['TT_S1', 'TT_S2', 'bin_offset_0'], axis=1, inplace=True)

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

def algo5d(scale = 1):
    result = pyreadr.read_r('../inputs/subset1.rds')
    df = result[None]
    bin_offset = 0

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
    sumos['starting_nt'] = sumos['bin_offset_0'].apply(lambda x: x.left)

    sumos.drop(['TT_S1', 'TT_S2', 'bin_offset_0'], axis=1, inplace=True)

    print("SUMOS")
    print(sumos[:50])
    print(sumos.shape)
    print(sumos.dtypes)

    ranges = read_ranges.get_ranges()
    print("sumos ilgis:", len(sumos))
    print("ranges ilgis:", len(ranges))

    sumos_features = pd.concat([sumos.reset_index(drop=True), ranges.reset_index(drop=True)], axis=1)
    model_input = sumos_features.copy()
    model_input.drop(['starting_nt', 'Chromosome', 'Start', 'End'], axis=1, inplace=True)
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
    kmeans = KMeans(n_clusters=2)
    kmeans.fit(model_input)
    inertia = kmeans.inertia_
    labels = list(kmeans.labels_)
    sumos_features['labels'] = kmeans.labels_
    open_ratio = labels.count(1)/len(labels)

    print(open_ratio)
    print(kmeans.labels_[:1000])

    sumos_features.to_csv('../outputs/output_5d.csv', sep = '\t')

    return sumos_features

def algo_prototypes(gamma = 1):
    result = pyreadr.read_r('../inputs/subset1.rds')
    df = result[None]
    bin_offset = 0

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
    sumos['starting_nt'] = sumos['bin_offset_0'].apply(lambda x: x.left)

    sumos.drop(['TT_S1', 'TT_S2', 'bin_offset_0'], axis=1, inplace=True)

    print("SUMOS")
    print(sumos[:50])
    print(sumos.shape)
    print(sumos.dtypes)

    ranges = read_ranges.get_ranges()
    print("sumos ilgis:", len(sumos))
    print("ranges ilgis:", len(ranges))

    sumos_features = pd.concat([sumos.reset_index(drop=True), ranges.reset_index(drop=True)], axis=1)
    model_input = sumos_features.copy()
    model_input.drop(['starting_nt', 'Chromosome', 'Start', 'End'], axis=1, inplace=True)
    print(model_input)

    # Skaičiuojam
    kprotos = KPrototypes(n_clusters=2, gamma=gamma, verbose=1)
    kprotos.fit(model_input, categorical=[1, 2, 3, 4])
    labels = list(kprotos.labels_)
    sumos_features['labels'] = kprotos.labels_
    open_ratio = labels.count(1)/len(labels)

    print(open_ratio)
    print(kprotos.labels_[:1000])

    sumos_features.to_csv('../outputs/output_prototypes.csv', sep = '\t')

    return sumos_features

def algo_dbscan(eps=0.05, min_samples=6):
    result = pyreadr.read_r('../inputs/subset1.rds')
    df = result[None]
    bin_offset = 0

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
    sumos['starting_nt'] = sumos['bin_offset_0'].apply(lambda x: x.left)

    sumos.drop(['TT_S1', 'TT_S2', 'bin_offset_0'], axis=1, inplace=True)

    print("SUMOS")
    print(sumos[:50])
    print(sumos.shape)
    print(sumos.dtypes)

    ranges = read_ranges.get_ranges()
    print("sumos ilgis:", len(sumos))
    print("ranges ilgis:", len(ranges))

    BATCH_SIZE = 10000
    sumos_features = pd.concat([sumos.reset_index(drop=True), ranges.reset_index(drop=True)], axis=1)
    model_input = sumos_features.copy()[:BATCH_SIZE]
    model_input.drop(['starting_nt', 'Chromosome', 'Start', 'End'], axis=1, inplace=True)
    print(model_input)

    # Skaičiuojam
    print("Length of input:", len(model_input))
    print("Calculating distances...")
    distance_matrix = gower.gower_matrix(model_input)
    clustering = DBSCAN(eps=eps, min_samples=min_samples, metric = "precomputed")

    print("Running DBSCAN...")
    clustering.fit(distance_matrix)

    labels = list(clustering.labels_)
    sumos_features.loc[:BATCH_SIZE - 1, 'labels'] = clustering.labels_
    open_ratio = labels.count(1)/len(labels)

    print(open_ratio)
    print(clustering.labels_)
    print("Number of clusters:", max(labels))

    sumos_features.to_csv('../outputs/output_dbscan.csv', sep = '\t')

    return sumos_features

if __name__ == "__main__":
    algo5d()
