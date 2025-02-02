import chromatin_algos
import evaluator
import enrichment
import pandas as pd
import generate_consensus
import pyreadr

out = generate_consensus.generate_consensus_5x20(chromatin_algos.algo_prototypes_aggregated, [1, 4, 1.5])
out.to_csv('../outputs/final4/ranges_kprototypes.csv', sep = '\t')

evaluation = evaluator.evaluate(merged_ranges = out, model_name="k-prototypes")
evaluation.to_csv('../outputs/final4/evaluate_kprototypes.csv', sep = '\t')

processed_out = evaluator.process_ranges(out)
processed_out.to_csv('../outputs/final4/ranges2_kprototypes.csv', sep = '\t')

evaluation2 = evaluator.evaluate(merged_ranges = processed_out, model_name="k-prototypes processed")
evaluation2.to_csv('../outputs/final4/evaluate2_kprototypes.csv', sep = '\t')

enrichment_coeffs = enrichment.calculate_enrichment_factors(processed_out)
enrichment_coeffs.to_csv('../outputs/final4/enrichment_kprototypes.csv', sep = '\t')
