import chromatin_algos
import evaluator
import pandas as pd
import generate_consensus

out = generate_consensus.generate_consensus_5x20(chromatin_algos.algo_prototypes_aggregated, [1, 4, 2.5])
out.to_csv('../outputs/final/ranges_kprototypes_s1.csv', sep = '\t')

processed_out = evaluator.process_ranges(out)
processed_out.to_csv('../outputs/final/ranges2_kprototypes_s1.csv', sep = '\t')
