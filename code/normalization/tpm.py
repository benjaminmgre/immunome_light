import pandas as pd
import numpy as np


def normalize_tpm(counts: pd.DataFrame):
    # Use https://github.com/lucynwosu/TPM-Transcripts-Per-Million-Normalization-Python/blob/main/TPM-Transcripts-Per-Million-Normalization.ipynb
    gene_lengths = pd.read_csv('data/Gencode/gencode.v43lift37.annotation.gtf')
    print(gene_lengths.head())

    counts.merge(gene_lengths, left_on='ensemble_id_version', )
