import chromatin_algos
import evaluator
import pandas as pd
import generate_consensus
import config
import sys

signal_column = sys.argv[1]
chr = config.CHR

out = generate_consensus.generate_consensus_5x20(chromatin_algos.algo5d_aggregated, [256, 8, 3, signal_column])
out.to_csv('../outputs/expression/ranges_chr' + chr + '_kmeans_' + signal_column + '.tsv', sep = '\t')

processed_out = evaluator.process_ranges(out)
processed_out.to_csv('../outputs/expression/ranges2_chr' + chr + '_kmeans_' + signal_column + '.tsv', sep = '\t')
