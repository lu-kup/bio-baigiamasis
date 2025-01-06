import chromatin_algos
import evaluator
import pandas as pd
import generate_consensus

ITERATIONS = 8
FILEPATH = '../outputs/sresults_kprototypes.csv'

thresholds = [1, 1.5, 2, 2.5, 3]
n_clusters_list = [2, 4]

result = evaluator.evaluate(merged_ranges = chromatin_algos.algo5d_aggregated(), model_name="test run")

for n_clusters in n_clusters_list:
    for threshold in thresholds:
        print("Runnning k-prototypes, n-clusters", n_clusters)
        out = generate_consensus.generate_consensus_5x20(chromatin_algos.algo_prototypes_aggregated, [1, n_clusters, threshold])
        evaluation = evaluator.evaluate(merged_ranges = out, model_name=f"k-prototypes 5 features n_clusters={n_clusters}, threshold={threshold}")
        result = pd.concat([result, evaluation], axis=1)
        result.transpose().to_csv(FILEPATH, sep = '\t')

result.transpose().to_csv(FILEPATH, sep = '\t')
