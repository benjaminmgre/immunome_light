import pandas as pd
import plotly.express as px

from data_processing import expression_matrix as em
from pathways import gene_set
from plotting import barchart


def z_normalization(x, mean, std):
    return (x - mean) / std

def score_pathways():
    sample = em.ExpressionDF('data/GSE205161/GSE205161_20220525-geo_raw_counts.csv', 'GSE205161')
    control = em.ExpressionDF('data/GTEx_data/gene_reads_2017-06-05_v8_whole_blood.gct', 'GTEx')

    gene_sets = gene_set.load_gene_sets('data/Gene_sets/GeneSets_June13_2023.xlsx')
    gene_set_names = list(gene_sets.keys())

    # Pathway scores dataframe will have samples for rows and gene sets for columns
    sample_ids = sample.get_headers()
    pathway_scores = pd.DataFrame(index=sample_ids, columns=gene_set_names)


    for name, g_set in gene_sets.items():
        # For each gene set, create a dataframe where genes are rows and sample #s are columns
        gene_scores = pd.DataFrame(columns=sample_ids)
        for gene in g_set:
            sample_counts = sample.get_gene_counts(gene)
            control_counts = control.get_gene_counts(gene)

            if sample_counts.empty:
                continue
            #TODO What if control counts is empty?

            control_std = control_counts.std()
            control_mean = control_counts.mean()

            # Calculate the z-score for each sample
            z_scores = sample_counts.apply(z_normalization, args=(control_mean, control_std,))

            # Next: create a dataframe initialized with column headings with the sample names and iteratively concatenate
            # the z-scores. Then, combine all the z-scores (find mean or something).

            # https://stackoverflow.com/questions/44156051/add-a-series-to-existing-dataframe
            gene_scores.loc[gene] = z_scores

        # Calculate the combined z score
        # https://sparkbyexamples.com/pandas/pandas-get-column-average-mean/#:~:text=To%20calculate%20the%20mean%20of,wise%20mean%20of%20the%20DataFrame.
        gene_mean = gene_scores.mean(axis=0)
        print(gene_mean)

        # Immunogram score is mean of z-scores + 3
        im_scores = gene_mean + 3

        # Score the pathway just using the first sample
        pathway_scores[name] = im_scores

    # pathway_df = pd.DataFrame(dict(
    #     r=pathway_scores.tolist(),
    #     theta=pathway_scores.keys(),
    # ))

    fig = px.line_polar(pathway_scores.T, r="CFB2071", theta=gene_set_names, range_r=[0, 6], line_close=True)
    fig.show()

    # fig = make_subplots(rows=4, cols=8)
    #
    # fig.add_trace(go.Scatterpolar(
    #     r=pathway_scores.iloc[0],
    #     theta=gene_set_names,
    #     mode='lines',
    #     name=pathway_scores.iloc[0].name,
    # ), 1, 1)
    #
    # fig.add_trace(go.Scatterpolar(
    #     r=pathway_scores.iloc[1],
    #     theta=gene_set_names,
    #     mode='lines',
    #     name=pathway_scores.iloc[1].name,
    # ), 1, 2)
    #
    # fig.add_trace(go.Scatterpolar(
    #     r=pathway_scores.iloc[2],
    #     theta=gene_set_names,
    #     mode='lines',
    #     name=pathway_scores.iloc[2].name,
    # ), 1, 3)
    #
    # fig.add_trace(go.Scatterpolar(
    #     r=pathway_scores.iloc[3],
    #     theta=gene_set_names,
    #     mode='lines',
    #     name=pathway_scores.iloc[3].name,
    # ), 1, 4)
    #
    # fig.show()

    barchart.score_save_barchart(pathway_scores)

