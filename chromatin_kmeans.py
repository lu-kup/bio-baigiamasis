import numpy as np
import pyreadr
import pandas as pd
import pyranges as pr
from sklearn.cluster import KMeans
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
    result = pyreadr.read_r('subset1.rds')
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

    np.savetxt('labels.dat', kmeans.labels_, fmt='%d', delimiter=',')

def algo2d():
    result = pyreadr.read_r('subset1.rds')
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

    max_mean = sumos['TT_S0'].nlargest(10).mean()
    print("Max mean", max_mean)
    scale = sumos['starting_nt'].max() / max_mean
    print("Scale", scale)
    sumos['signal'] = sumos['TT_S0'] * scale * 10
    sumos.drop(['TT_S0'], axis=1, inplace=True)

    print("SUMOS NEW")
    print(sumos[:50])

    # Skaičiuojam
    kmeans = KMeans(n_clusters=2)
    kmeans.fit(sumos)
    inertia = kmeans.inertia_
    labels = list(kmeans.labels_)
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

    np.savetxt('labels.dat', kmeans.labels_, fmt='%d', delimiter=',')

def algo5d():
    result = pyreadr.read_r('subset1.rds')
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

    """     max_mean = sumos['TT_S0'].nlargest(10).mean()
    print("Max mean", max_mean)
    scale = sumos['starting_nt'].max() / max_mean
    print("Scale", scale)
    sumos['signal'] = sumos['TT_S0'] * scale * 10
    sumos.drop(['TT_S0'], axis=1, inplace=True)

    print("SUMOS NEW")
    print(sumos[:50]) """

    ranges = read_ranges.get_ranges()
    print("sumos ilgis:", len(sumos))
    print("ranges ilgis:", len(ranges))

    sumos_features = pd.concat([sumos.reset_index(drop=True), ranges.reset_index(drop=True)], axis=1)
    print(sumos_features)

    return
    # Skaičiuojam
    kmeans = KMeans(n_clusters=2)
    kmeans.fit(sumos)
    inertia = kmeans.inertia_
    labels = list(kmeans.labels_)
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

    np.savetxt('labels.dat', kmeans.labels_, fmt='%d', delimiter=',')

if __name__ == "__main__":
    algo5d()
