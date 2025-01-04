import chromatin_algos
import evaluator
import pandas as pd
import generate_consensus

out = generate_consensus.generate_consensus_5x20(chromatin_algos.algo5d_aggregated, [16, 10, 3])
out.to_csv('../outputs/final/ranges_kmeans_s2.csv', sep = '\t')

processed_out = evaluator.process_ranges(out)
processed_out.to_csv('../outputs/final/ranges2_kmeans_s2.csv', sep = '\t')
