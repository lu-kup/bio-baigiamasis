import numpy as np
import pyreadr
import pandas as pd
import pyranges as pr
from sklearn.cluster import KMeans

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

def main():
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

    print("SUMOS")
    print(sumos.head())
    print(sumos.shape)

    data = list(zip(x, y))
    kmeans = KMeans(n_clusters=2)
    kmeans.fit(data)
    inertia = kmeans.inertia_
    print(kmeans.labels_[:1000])

if __name__ == "__main__":
    main()
