import chromatin_algos
import evaluator
import pandas as pd
import generate_consensus

ITERATIONS = 10
FILEPATH = '../outputs/zresults_kmeans_scale.tsv'

thresholds = [1, 1.5, 2, 2.5, 3]

result = evaluator.evaluate(merged_ranges = chromatin_algos.algo5d_aggregated(), model_name="test run")

for i in range(ITERATIONS):
    for threshold in thresholds:
        scale = 2 ** i
        print("Runnning k-means, 10 clusters, scale=", scale)
        out = generate_consensus.generate_consensus_5x20(chromatin_algos.algo5d_aggregated, [scale, 10, threshold])
        evaluation = evaluator.evaluate(merged_ranges = out, model_name=f"k-means 5 features n_clusters=10, scale={scale}, threshold={threshold}")
        result = pd.concat([result, evaluation], axis=1)
        result.transpose().to_csv(FILEPATH, sep = '\t')

result.transpose().to_csv(FILEPATH, sep = '\t')
