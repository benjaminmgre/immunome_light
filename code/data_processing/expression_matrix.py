import numpy as np
import pandas as pd


class ExpressionDF:
    def __init__(self, expression_data, data_format):
        """
        Generates a new expression dataframe from the specified format
        :param expression_data: the data file path
        :param data_format: format type of the data. Currently supported: 'GSE205161', 'GTEx'
        """
        with open(expression_data) as data_fh:
            if data_format == 'GSE205161':
                self.expression_df = pd.read_table(data_fh, sep=',', header=0, index_col=0)
            elif data_format == 'GTEx':
                self.expression_df = pd.read_table(data_fh, sep='\t', header=2, index_col=1)
                # Remove the first two unnecessary columns
                # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.drop.html#pandas-dataframe-drop
                self.expression_df.drop(['id', 'Description'], inplace=True, axis=1)
        return

    def get_gene_counts(self, ensembl_name):
        """
        Gets the sample gene counts from the expression dataframe
        :param ensembl_name:
        :return:
        """
        # Query the database for any gene that begins with ensembl name
        # (some ensembl ids have version suffix that needs to be ignored)
        query_df = self.expression_df[self.expression_df.index.str.startswith(ensembl_name)]
        if len(query_df.index) > 1:
            print(f'Error! There are multiple rows for {ensembl_name} in dataframe!')
            return query_df.iloc[0, :]
        elif query_df.empty:
            print(f'No records exist for {ensembl_name} in sample.')
            return query_df
        else:
            # Extract the first and only row and only the sample expression counts and return as pandas series
            # self.ensembl_col: + 1 only selects columns past the gene id (i.e. it selects all the sample count columns)
            return query_df.iloc[0, :]
        # TODO: RETURN EMPTY SERIES


    def get_headers(self):
        """
        Returns the headers of the columns
        :return:
        """
        return list(self.expression_df.columns.values)