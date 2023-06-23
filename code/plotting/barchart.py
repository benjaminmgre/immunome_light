import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math


def score_barchart(data: pd.DataFrame):
    # n_rows is n samples and n cols is n gene sets
    n_rows = data.shape[0]
    n_cols = data.shape[1]

    n_fig_cols = math.ceil(n_rows / 3)
    fig, axes = plt.subplots(3, n_fig_cols)

    i = 0
    for (sample_index, row) in data.iterrows():
        axs = sns.barplot(x=data.iloc[i].index, y=data.iloc[i].values, ax=axes[i // n_fig_cols][i % n_fig_cols],
                          label=sample_index)
        axs.set_ylim(bottom=0, top=6)
        i += 1

    plt.show()
    return


def score_save_barchart(data: pd.DataFrame):
    # n_rows is n samples and n cols is n gene sets
    n_rows = data.shape[0]
    n_cols = data.shape[1]

    i = 0
    for (sample_index, row) in data.iterrows():
        plt.figure()
        axs = sns.barplot(x=data.iloc[i].index, y=data.iloc[i].values,
                          label=sample_index)
        axs.set_ylim(bottom=0, top=6)
        plt.savefig(f'generated_files/plots/{sample_index}')
        plt.close()
        i += 1
    return
