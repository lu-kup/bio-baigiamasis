import chromatin_algos
import evaluator
import pandas as pd
import generate_consensus

ITERATIONS = 10
FILEPATH = '../outputs/zresults_dbscan.tsv'

thresholds = [1, 1.5, 2, 2.5, 3]

result = evaluator.evaluate(merged_ranges = chromatin_algos.algo5d_aggregated(), model_name="test run")

for i in range(ITERATIONS):
    for threshold in thresholds:
        eps = 0.00001 * (2 ** i)
        print("Runnning DBSCAN, epsilon", eps)
        out = generate_consensus.generate_consensus_5x20(chromatin_algos.algo_dbscan_aggregated, [eps, 6, threshold])
        evaluation = evaluator.evaluate(merged_ranges = out, model_name=f"DBSCAN eps={eps}, threshold={threshold}")
        result = pd.concat([result, evaluation], axis=1)
        result.transpose().to_csv(FILEPATH, sep = '\t')

result.transpose().to_csv(FILEPATH, sep = '\t')
