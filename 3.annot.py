import os
import pandas as pd
import numpy as np


work_dir = "/home/minzhang/workspace/RNA-Seq"
gtf_file = f"{work_dir}/gencode.v47.basic.annotation.gtf"

dir_suffix = {
    'expr': '.genes.results',
    'deg': '.csv',
}


# process GTF file to get gene/transcrpt informations
genes = dict()
transcripts = dict()
with open(gtf_file, 'r') as gtf:
    for line in gtf:
        # print(line)
        if line.startswith("#"):
            continue
        columns = line.strip().split("\t")
        attributes = columns[-1]
        #print(attributes)
        attributes = {i.split(maxsplit=1)[0]: i.split(maxsplit=1)[1].strip('"')
                      for i in attributes.split(";")
                      if len(i) > 1}
        if columns[2] == "gene":
            genes[attributes["gene_id"]] = attributes
        elif columns[2] == "transcript":
            transcripts[attributes["transcript_id"]] = attributes
gene_attr_list = [set(genes[gene_id].keys()) for gene_id in genes]
gene_attr_list = list(set.union(*gene_attr_list))
transcript_attr_list = [set(transcripts[transcript_id].keys()) for transcript_id in transcripts]
transcript_attr_list = list(set.union(*transcript_attr_list))
gene_data = []
transcript_data = []

for gene_id in genes:
    gene_info = genes[gene_id]
    gene_row = {attr: gene_info.get(attr, np.nan) for attr in gene_attr_list}
    gene_data.append(gene_row)
gene_df = pd.DataFrame(gene_data, columns=gene_attr_list)
gene_df.set_index("gene_id", inplace=True)
for transcript_id in transcripts:
    transcript_info = transcripts[transcript_id]
    transcript_row = {attr: transcript_info.get(attr, np.nan) for attr in transcript_attr_list}
    transcript_data.append(transcript_row)
transcript_df = pd.DataFrame(transcript_data, columns=transcript_attr_list)
transcript_df.set_index("transcript_id", inplace=True)
gene_df.to_csv(f"{work_dir}/gene_info.csv")
transcript_df.to_csv(f"{work_dir}/transcript_info.csv")


# Add gene symbols to indicated tables
gene_df = pd.read_csv(f"{work_dir}/gene_info.csv", index_col="gene_id")
gene_df = gene_df[['gene_name', 'gene_type']]
for dir in dir_suffix:
    suffix = dir_suffix[dir]
    dir = f"{work_dir}/{dir}"
    for expr_file in os.listdir(dir):
        if expr_file.endswith(suffix):
            expr_path = os.path.join(dir, expr_file)
            if suffix.endswith("csv"):
                expr_df = pd.read_csv(expr_path, sep=",")
            else:
                expr_df = pd.read_csv(expr_path, sep="\t")
            expr_df = pd.merge(expr_df, gene_df, on="gene_id", how="left")
            columns_order = ['gene_name', 'gene_type'] + [col for col in expr_df.columns if col not in ['gene_name', 'gene_type']]
            expr_df = expr_df[columns_order]
            output_file = os.path.join(dir, expr_file.replace(suffix, f"{suffix}.annoted.csv"))
            expr_df.to_csv(output_file, index=False)
