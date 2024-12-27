import chromatin_algos
import evaluator
import pandas as pd

ITERATIONS = 10
FILEPATH = '../outputs/results_regression_test.csv'

thresholds = [1, 1.5, 2, 2.5, 3]

result = evaluator.evaluate(dataframe = chromatin_algos.algo5d(), model_name="test run")

for i in range(ITERATIONS):
    for threshold in thresholds:
        eps = 0.0001 * (2 ** i)
        print("Runnning DBSCAN, epsilon", eps)
        evaluation = evaluator.evaluate(dataframe = chromatin_algos.algo_dbscan_aggregated(eps=eps, threshold=threshold), model_name=f"DBSCAN eps={eps}, threshold={threshold}")
        result = pd.concat([result, evaluation], axis=1)
        result.to_csv(FILEPATH, sep = '\t')

for i in range(ITERATIONS):
    for threshold in thresholds:
        n_clusters = 2 ** i
        print("Runnning k-means, n-clusters", n_clusters)
        evaluation = evaluator.evaluate(dataframe = chromatin_algos.algo5d_aggregated(n_clusters=n_clusters, threshold=threshold), model_name=f"k-means 5 features n_clusters={n_clusters}, threshold={threshold}")
        result = pd.concat([result, evaluation], axis=1)
        result.to_csv(FILEPATH, sep = '\t')

for i in range(ITERATIONS):
    for threshold in thresholds:
        n_clusters = 2 ** i
        print("Runnning k-prototypes, n-clusters", n_clusters)
        evaluation = evaluator.evaluate(dataframe = chromatin_algos.algo_prototypes_aggregated(n_clusters=n_clusters, threshold=threshold), model_name=f"k-prototypes 5 features n_clusters={n_clusters}, threshold={threshold}")
        result = pd.concat([result, evaluation], axis=1)
        result.to_csv(FILEPATH, sep = '\t')

result.to_csv(FILEPATH, sep = '\t')
