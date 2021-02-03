#!/usr/bin/env python
# coding: utf-8

# In[16]:


import numpy as np
import pandas as pd
import scanpy as sc
import os
import os.path


# In[18]:


meta_folder = "/net/data.isilon/ag-saez/bq_shared/scellMI/processed_RNA/meta_data/"
meta_onlyfiles = [f for f in os.listdir(meta_folder) if os.path.isfile(os.path.join(meta_folder, f))]


# In[55]:


raw_folder = "/net/data.isilon/ag-saez/bq_shared/scellMI/raw_RNA/"
raw_onlyfiles = [f for f in os.listdir(raw_folder)]


# In[27]:


meta_index = [f.split(r".")[0] for f in meta_onlyfiles]
raw_index = [f.split(r".")[0] for f in raw_onlyfiles]


# In[ ]:


out_folder = "/net/data.isilon/ag-saez/bq_shared/scellMI/processed_RNA/pymirrors/"


# In[39]:


raw_df = pd.DataFrame({'data_file': [raw_folder + f + "/filtered_feature_bc_matrix/" for f in raw_onlyfiles]}, 
                  index = raw_index)


# In[50]:


path_df = pd.DataFrame({'meta_file': [meta_folder + f for f in meta_onlyfiles]}, 
                  index = meta_index)


# In[40]:


path_df["out_file"] = [out_folder + f + r".h5ad" for f in path_df.index.values]
path_df["data_file"] = raw_df.loc[path_df.index.values, "data_file"]


# In[49]:


def generate_mirror(row):
    mfile = row[0]
    ofile = row[1]
    dfile = row[2]
    adata = sc.read_10x_mtx(dfile,  # the directory with the `.mtx` file
                            var_names = 'gene_symbols', # use gene symbols for the variable names (variables-axis index)
                            cache = False)
    mdata = pd.read_table(mfile)
    mdata["id"] = [i.split("-")[0] for i in mdata["id"]]
    adata.var_names_make_unique()
    
    # Just adjust the labels -----------------------------------------
    adata.obs["raw_cell_id"] = adata.obs.index.values
    adata.obs["id"] = [i.split("-")[0] for i in adata.obs.index.values]
    # Merge cell annotations -----------------------------------------
    cell_anns = pd.merge(adata.obs, mdata, on='id', how='left')
    adata.obs["cell_type"] = pd.Series.tolist(cell_anns["cell_type"])
    adata.obs["major_cell_type"] = pd.Series.tolist(cell_anns["major_cell_type"])
    print(adata.shape)
    adata = adata[[not x for x in pd.isna(adata.obs.cell_type)], :]
    print(adata.shape)
    adata.write(ofile)
    return(adata.shape)


# In[6]:


out = path_df.apply(generate_mirror, axis=1, result_type = 'expand')

