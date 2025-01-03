import chromatin_algos
import evaluator
import enrichment
import pandas as pd
import generate_consensus
import pyreadr

out = generate_consensus.generate_consensus_5x20(chromatin_algos.algo_dbscan_aggregated, [0.001, 6, 3])
out.to_csv('../outputs/final/ranges_dbscan.csv', sep = '\t')

evaluation = evaluator.evaluate(merged_ranges = out, model_name="dbscan")
evaluation.to_csv('../outputs/final/evaluate_dbscan.csv', sep = '\t')

processed_out = evaluator.process_ranges(out)
processed_out.to_csv('../outputs/final/ranges2_dbscan.csv', sep = '\t')

evaluation2 = evaluator.evaluate(merged_ranges = processed_out, model_name="dbscan processed")
evaluation2.to_csv('../outputs/final/evaluate2_dbscan.csv', sep = '\t')

targets_df = pyreadr.read_r('../inputs/subset2.rds')[None]
enrichment_coeffs = enrichment.calculate_enrichment_factors(processed_out, targets_df)
enrichment_coeffs.to_csv('../outputs/final/enrichment_dbscan.csv', sep = '\t')
