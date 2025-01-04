import chromatin_algos
import evaluator
import pandas as pd
import generate_consensus

out = generate_consensus.generate_consensus_5x20(chromatin_algos.algo_dbscan_aggregated, [0.001, 6, 3])
out.to_csv('../outputs/final/ranges_dbscan_s2.csv', sep = '\t')

processed_out = evaluator.process_ranges(out)
processed_out.to_csv('../outputs/final/ranges2_dbscan_s2.csv', sep = '\t')
