import chromatin_algos
import evaluator
import pandas as pd
import generate_consensus

ITERATIONS = 10
FILEPATH = '../outputs/zresults_kprototypes_gamma.tsv'

thresholds = [2, 3]

result = evaluator.evaluate(merged_ranges = chromatin_algos.algo5d_aggregated(), model_name="test run")

for threshold in thresholds:
    for i in range(ITERATIONS):
        gamma = 2 ** i
        print("Runnning k-prototypes, 2 clusters, gamma=", gamma)
        out = generate_consensus.generate_consensus_5x20(chromatin_algos.algo_prototypes_aggregated, [gamma, 2, threshold])
        evaluation = evaluator.evaluate(merged_ranges = out, model_name=f"k-prototypes 5 features 2 clusters, gamma={gamma}, threshold={threshold}")
        result = pd.concat([result, evaluation], axis=1)
        result.transpose().to_csv(FILEPATH, sep = '\t')

result.transpose().to_csv(FILEPATH, sep = '\t')
