import numpy as np
import pandas as pd


from data_processing import gene_id


def load_gene_sets(gene_file):
    """
    Loads all gene sets from excel sheet
    :param gene_file:
    :return:
    """
    gene_sets_df = pd.read_excel(gene_file, header=0, engine='openpyxl')
    print(gene_sets_df.to_string())

    hgnc_genes = []
    ensembl_genes = []
    for (name, column) in gene_sets_df.items():
        for (i, gene_name) in column.items():
            if gene_name not in hgnc_genes:
                hgnc_genes.append(gene_name)

    # Convert hgnc set to ensembl set, and then make ensembl gene sets
    gid_df = gene_id.GeneID(hgnc_genes, "HGNC_symbol")
    ensembl_genes = gid_df.id_conversion(hgnc_genes, "HGNC_symbol", "Ensembl")

    # Generate an HGNC symbol to Ensembl dict
    hgnc2ens = dict(zip(hgnc_genes, ensembl_genes))

    # Generate dict of gene sets
    gene_sets = {}
    for (name, column) in gene_sets_df.items():
        gene_sets[name] = [hgnc2ens[hgnc] for hgnc in column.tolist()]

    # Remove all 'NA' values
    for (name, g_list) in gene_sets.items():
        gene_sets[name] = [gene for gene in g_list if gene != 'NA']

    return gene_sets
