import chromatin_algos
import evaluator
import pandas as pd

ITERATIONS = 10

result = evaluator.evaluate(dataframe = chromatin_algos.algo1d(), model_name="k-means 1 dimension")

for i in range(ITERATIONS):
    scale = 2 ** i
    print("Runnning k-means, scale", scale)
    evaluation = evaluator.evaluate(dataframe = chromatin_algos.algo5d(scale), model_name=f"k-means 5 features scale={scale}")
    result = pd.concat([result, evaluation], axis=1)
    result.to_csv('../outputs/results2.csv', sep = '\t')

for i in range(ITERATIONS):
    gamma = 2 ** i
    print("Runnning k-prototypes, gamma", gamma)
    evaluation = evaluator.evaluate(dataframe = chromatin_algos.algo_prototypes(gamma), model_name=f"k-prototypes 5 features gamma={gamma}")
    result = pd.concat([result, evaluation], axis=1)
    result.to_csv('../outputs/results2.csv', sep = '\t')

for i in range(ITERATIONS):
    eps = 0,1 * (2 ** i)
    print("Runnning DBSCAN, epsilon", eps)
    evaluation = evaluator.evaluate(dataframe = chromatin_algos.algo_dbscan(eps), model_name=f"DBSCAN eps={eps}")
    result = pd.concat([result, evaluation], axis=1)
    result.to_csv('../outputs/results2.csv', sep = '\t')

result.to_csv('../outputs/results2.csv', sep = '\t')