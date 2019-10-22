import pandas as pd

##################################
##
## Immunology data
##
#################################

immuno_gene_counts_files = "/home/team9_hs19/immunoData/dice/"
ctypes = ["B_CELL_NAIVE", "CD4_NAIVE", "CD4_STIM", "CD8_NAIVE", "CD8_STIM", "M2", "MONOCYTES", "NK", "TFH", "TH17", "TH1", "TH2", "THSTAR", "TREG_MEM", "TREG_NAIVE"]

# reading in celltype sepcific data
data_list = []
for ctype in ctypes:
    data_list.append( pd.read_csv(immuno_gene_counts_files + ctype + "_TPM.csv", index_col=0)
            .drop(["Transcript_Length(bp)", "Additional_annotations"], axis=1)
            .mean(axis=1))

# aggregating datat to one dataframe
imm_aggregate = pd.concat(data_list, axis=1, sort=True).max(axis=1)

# have to reindex because brain data does not have sub index in the Ensemb gene names
imm_aggregate.index = [ x.strip().split(".")[0] for x in imm_aggregate.index]

# Have to sum up the counts if the gene name collapse due to reindexing
imm_aggregate = imm_aggregate.groupby(by = imm_aggregate.index).sum()

imm_genes = set(imm_aggregate.index)

##################################
##
## Brain data
##
#################################

# This Assumes we have filtered the protein coding genes right
brain_gene_counts_file = "/home/team9_hs19/immunoData/avi/data/brain/brain_mean_counts.csv"
brain_aggregate = pd.read_csv(brain_gene_counts_file, index_col=0)["max"]
brain_genes = set(brain_aggregate.index)

##################################
##
## Generating Protein Interaction 
##
#################################

ppi_file = "/home/team9_hs19/immunoData/avi/data/brain/human-PPI-pairs.csv"
ppi_ntwk = pd.read_csv(ppi_file, skipinitialspace=True, usecols=["ENSEMBL_A", "ENSEMBL_B"]).dropna()
ppi_ntwk = ppi_ntwk[ ppi_ntwk["ENSEMBL_A"] != ppi_ntwk["ENSEMBL_B" ] ]

ppi_ntwk = ppi_ntwk[ ppi_ntwk["ENSEMBL_A"].isin(brain_genes) ]
ppi_ntwk = ppi_ntwk[ ppi_ntwk["ENSEMBL_B"].isin(brain_genes) ]
ppi_ntwk = ppi_ntwk[ ppi_ntwk["ENSEMBL_A"].isin(imm_genes) ]
ppi_ntwk = ppi_ntwk[ ppi_ntwk["ENSEMBL_B"].isin(imm_genes) ]

##################################
##
## Extracting ppi for imm & Brain
##
#################################

ppi_ntwk["brain_a"] = brain_aggregate.loc[ ppi_ntwk["ENSEMBL_A"] ].values
ppi_ntwk["brain_b"] = brain_aggregate.loc[ ppi_ntwk["ENSEMBL_B"] ].values

ppi_ntwk["immo_a"] = imm_aggregate.loc[ ppi_ntwk["ENSEMBL_A"] ].values
ppi_ntwk["immo_b"] = imm_aggregate.loc[ ppi_ntwk["ENSEMBL_B"] ].values

######################################
##
## Extracting plotting cooredinates
##
#####################################


ppi_ntwk["X"] = ppi_ntwk[["brain_a", "brain_b"]].min(axis=1)
ppi_ntwk["Y"] = ppi_ntwk[["immo_a", "immo_b"]].min(axis=1)

col_names_out = ["ENSEMBL_A", "ENSEMBL_B", "X", "Y"]
ppi_ntwk[col_names_out].to_csv("plotting.txt", index=False)
