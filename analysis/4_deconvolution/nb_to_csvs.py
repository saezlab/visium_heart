#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
# coding: utf-8
import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
import gc
data_type = 'float32'


# In[29]:


# Here we will define the path location of all nb models

# Set paths to data and results used through the document:
reg_data_folder = '/Users/ricardoramirez/Dropbox/PhD/Research/mi_atlas/results/nb_estimates/'
reg_data_folder_out = '/Users/ricardoramirez/Dropbox/PhD/Research/mi_atlas/results/nb_estimates_csv/'
nb_folders = [f for f in os.listdir(reg_data_folder) if os.path.isdir(reg_data_folder + f)]
nb_out = [reg_data_folder + f + "/sc.h5ad" for f in nb_folders]


# In[15]:


def get_infer_average(nb_file):
    # First generate the reference data frame
    adata_snrna_raw = sc.read(nb_file)
    # Column name containing cell type annotations
    covariate_col_names = 'cell_type'
    # Extract a pd.DataFrame with signatures from anndata object
    inf_aver = adata_snrna_raw.raw.var.copy()
    inf_aver = inf_aver.loc[:, [f'means_cov_effect_{covariate_col_names}_{i}' for i in adata_snrna_raw.obs[covariate_col_names].unique()]]
    from re import sub
    inf_aver.columns = [sub(f'means_cov_effect_{covariate_col_names}_{i}', '', i) for i in adata_snrna_raw.obs[covariate_col_names].unique()]
    inf_aver = inf_aver.iloc[:, inf_aver.columns.argsort()]
    # normalise by average experiment scaling factor (corrects for sequencing depth)
    inf_aver = inf_aver * adata_snrna_raw.uns['regression_mod']['post_sample_means']['sample_scaling'].mean()
    # eliminate duplicates in inf aver too just in case
    x = pd.Series([i for i in inf_aver.index.values]) 
    inf_aver = inf_aver.iloc[~x.duplicated(keep = "first").values,]
    return(inf_aver)


# In[23]:


inf_aver_list = {i: get_infer_average(nb_out[i]) for i in range(0,5)}


# In[30]:


for i in range(0,5):
    inf_aver_list[i].to_csv(reg_data_folder_out + nb_folders[i] + ".csv")


# In[19]:





# In[28]:





# In[ ]:




