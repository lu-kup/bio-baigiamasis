import pyranges as pr

#path = pr.get_example_path("gencode_chr18.gtf")
gr = pr.read_gtf("gencode_chr18.gtf")

genes = gr[gr.Feature == 'gene']

print(genes)