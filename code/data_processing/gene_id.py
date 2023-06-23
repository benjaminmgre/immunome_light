import pandas


from config import config


def id_tab(id_format):
    """
    Converts the given ID format string into the string which can be used to index the gene
    identity dataframes.

    Parameters:
    ----------
    id_format: str
        A string representing the identification forms of the initial gene list.

    Returns:
    --------
    id_t: str
        A string that can be used to index the gene identification dataframes. Matches id_format
    """
    if id_format == 'HGNC_symbol':
        id_t = 'Approved symbol'
    elif id_format == 'HGNC_id':
        id_t = 'HGNC ID'
    elif id_format == 'RefSeq':
        id_t = 'RefSeq IDs'
    elif id_format == 'NCBI':
        id_t = 'NCBI Gene ID'
    elif id_format == 'Ensembl':
        id_t = 'Ensembl gene ID'
    else:
        id_t = ''
        print("""
            Error! id_type was not entered correctly. Valid str: "HGNC_symbol", "HGNC_id", "RefSeq", "NCBI",
            "Ensembl".
              """)
    return id_t


class GeneID:
    """
    A dataframe for referencing various forms of a gene identifications.
    Supports HGNC IDs, HGNC Approved Symbols, RefSeq IDs, NCBI Gene IDs, Ensembl gene IDs

    Attributes:
    -----------
    self.id_format: str
        A dataframe column header which represents the given identification form from the gene list.
    self.gid_df: pandas dataframe
        Gene Identifier Dataframe: Generated from the hgnc identification file
    self.gene_log: pandas dataframe
        A filtered self.gid_df: Only contains genes given in the gene list
    """

    def __init__(self, gene_list, id_form):
        """
        Initiates a gene identity dataframe for all genes in the gene list.

        Parameters:
        -----------
        gene_list: list
            A list of genes in the specified id format
        id_form: str
            String representing the ID form of the gene list.
            'HGNC_symbol': Approved HGNC symbols
            'HGNC_id': HGNC ID
            'RefSeq': RefSeq ID
            'NCBI': NCBI Gene ID
            'Ensembl': Ensembl gene ID
        """
        self.id_format = id_tab(id_form)
        with open(config.hgnc_filename) as gid:
            self.gid_df = pandas.read_table(gid, dtype='str')
        self.gene_searched_list = []
        self.unknown_gene_list = []
        self.gene_log = pandas.DataFrame()
        for gene in gene_list:
            # DEBUG: print(gene)
            # Check for nan genes
            if pandas.isnull(gene):
                continue
            if gene in self.gene_searched_list:
                continue
            else:
                gene_search = self.gid_df[self.gid_df[self.id_format] == gene]
                # print(gene_search)
                if gene_search.empty:
                    # print(f'DEBUG: {gene} is not in GID dictionary.')
                    self.unknown_gene_list.append(gene)
                else:
                    if len(gene_search) > 1:
                        gene_search = gene_search[gene_search["Status"] == "Approved"]
                    # Only add the gene to the log if there are results
                    self.gene_log = pandas.concat([self.gene_log, gene_search])
                self.gene_searched_list.append(gene)

        # Search for unknown genes
        # Unknown genes can be either alias, previous symbols, or just unknown!
        # Alias gene symbols will be prioritized over previous symbols
        # print("DEBUG: Starting search for unknown genes")
        self.alias_dict = {}
        self.previous_dict ={}
        for i in range(0, len(self.gid_df)):
            # If there are no more unknown genes, stop iterating
            if len(self.unknown_gene_list) == 0:
                break
            entry = self.gid_df.loc[i]
            alias = entry.at["Alias symbols"]
            previous = entry.at["Previous symbols"]
            # Previous or alias will be NaN (null) if no entries
            if not pandas.isnull(alias):
                alias_list = alias.split(", ")
                for a in alias_list:
                    if a in self.unknown_gene_list:
                        if a in self.alias_dict.keys():
                            print(f'ERROR! {a} has multiple alias!')
                        self.alias_dict[a] = entry
            if not pandas.isnull(previous):
                previous_list = previous.split(", ")
                for p in previous_list:
                    if p in self.unknown_gene_list:
                        if p in self.previous_dict.keys():
                            print(f'ERROR! {p} has multiple previous symbols!')
                        self.previous_dict[p] = entry
        # print('DEBUG: Done iterating through genes for alias and previous')
        # print(f'DEBUG: alias dict: {self.alias_dict} previous: {self.previous_dict}')
        return

    def id_conversion(self, gene_list, from_id, to_id):
        """
        Converts a list of genes from one identification form into another.

        Parameters:
        -----------
        gene_list: list
            Genes to be returned in a different form
        from_id: str
            The ID format the gene list is in
        to_id: str
            The ID format the gene list will be converted to

        See __init__ for ID forms.

        Returns:
        --------
        converted_list: list
            Genes listed in the to_id format.
        """
        from_form = id_tab(from_id)
        to_form = id_tab(to_id)
        converted_list = []
        for gene in gene_list:
            # Check that there is a valid gene
            # pandas.isnull(gene) checks for nan values
            if pandas.isnull(gene) or gene == 'None' or gene == 'NA':
                converted_list.append('NA')
                continue
            gene_search = self.gene_log[self.gene_log[from_form] == gene]
            if gene_search.empty:
                if gene in self.alias_dict.keys():
                    converted_form = self.alias_dict[gene].at[to_form]
                elif gene in self.previous_dict.keys():
                    converted_form = self.previous_dict[gene].at[to_form]
                else:
                    print(f'ERROR: {gene} can not be converted from {from_id} to {to_id}.')
            elif len(gene_search) > 1:
                print(f'ERROR: Multiple {to_id} IDs available for {from_id} {gene}.')
                print(gene_search.to_string())
                converted_form = gene_search.iloc[0][to_form]
            else:
                converted_form = gene_search.iloc[0][to_form]
            converted_list.append(converted_form)
        return converted_list

    def single_conversion(self, gene, from_id, to_id):
        """
        Duplicate of id_conversion method, but only takes a single gene and returns a single gene.
        """
        from_form = id_tab(from_id)
        to_form = id_tab(to_id)
        gene_search = self.gene_log[self.gene_log[from_form] == gene]
        if gene_search.empty:
            converted_form = gene
            print(f'DEBUG: {gene} cannot be converted from {from_id} to {to_id}.')
        elif len(gene_search) > 1:
            print(f'DEBUG: Multiple {to_id} IDs available for {from_id} {gene}.')
            print(gene_search)
            converted_form = gene_search.iloc[0][to_form]
        else:
            converted_form = gene_search.iloc[0][to_form]
        return converted_form
