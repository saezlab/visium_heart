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

# this line forces theano to use the GPU and should go before importing cell2location
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=' + data_type + ',force_device=True'

import cell2location

import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns

# Set paths to data and results used through the document:
sp_data_folder = '/net/data.isilon/ag-saez/bq_shared/scellMI/raw_visium/'
results_folder = '/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/deconvolution/c2l_statescollapsed/'
regression_model_output = 'RegressionGeneBackgroundCoverageTorch_27covariates_51906cells_12443genes_statescollapsed'
reg_path = f'{results_folder}{regression_model_output}/'
sample_names = ["157771", "157772", "157775", "157777", "157779", "157781", "157782", "157785"]

# Function to import a single slide
def read_and_qc(sample_name, path=sp_data_folder, force_filter = True):
    r""" This function reads the data for one 10X spatial experiment into the anndata object.
    It also calculates QC metrics. Modify this function if required by your workflow.

    :param sample_name: Name of the sample
    :param path: path to data
    """

    adata = sc.read_visium(path + str(sample_name) + '/outs',
                           count_file='filtered_feature_bc_matrix.h5', 
                           load_images=True)
    
    adata.obs['sample'] = sample_name
    adata.var['SYMBOL'] = adata.var_names

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.var['mt'] = [gene.startswith('MT-') for gene in adata.var['SYMBOL']]
    adata.var['rps'] = [gene.startswith('RPS') for gene in adata.var['SYMBOL']]
    adata.var['mrp'] = [gene.startswith('MRP') for gene in adata.var['SYMBOL']]
    adata.var['rpl'] = [gene.startswith('RPL') for gene in adata.var['SYMBOL']]
    adata.obs['mt_frac'] = adata[:,adata.var['mt'].tolist()].X.sum(1).A.squeeze()/adata.obs['total_counts']

    # add sample name to obs names
    adata.obs["sample"] = [str(i) for i in adata.obs['sample']]
    adata.obs_names = adata.obs["sample"] + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'
    adata.var["duplicated"] = adata.var['SYMBOL'].duplicated(keep = "first")
    adata = adata[:, ~adata.var['duplicated'].values]
    
    if force_filter:
        # First filter: mt and rb genes
        # mitochondria-encoded (MT) genes should be removed for spatial mapping
        adata.obsm['mt'] = adata[:,   adata.var['mt'].values | 
                              adata.var['rps'].values |
                              adata.var['mrp'].values |
                              adata.var['rpl'].values].X.toarray() 
        
        adata = adata[:, ~ (adata.var['mt'].values | 
                              adata.var['rps'].values |
                              adata.var['mrp'].values |
                              adata.var['rpl'].values)]
        
        # Second filter
        # Genes expressed in less than 10 spots
        adata = adata[:, adata.var['n_cells_by_counts'].values > 10]
        
        # Third filter
        # spots with no information (less than 300 genes and 500 UMIs)
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        adata = adata[(adata.obs['n_genes_by_counts'].values > 300) & 
              (adata.obs['total_counts'].values > 500), :] 
        

    return adata


# First generate the reference data frame
adata_snrna_raw = sc.read(f'{reg_path}sc.h5ad')
# Column name containing cell type annotations
covariate_col_names = 'deconv_col'
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

del adata_snrna_raw
gc.collect()

# Now do the whole pipeline per slide

for s in sample_names:
    
    print(s)
    
    adata = read_and_qc(s, path = sp_data_folder)
    adata_vis_s = adata.copy()
    adata_vis_s.raw = adata_vis_s
    
    r = cell2location.run_cell2location(

      # Single cell reference signatures as pd.DataFrame
      # (could also be data as anndata object for estimating signatures
      #  as cluster average expression - `sc_data=adata_snrna_raw`)
      sc_data=inf_aver,
      # Spatial data as anndata object
      sp_data=adata_vis_s,

      # the column in sc_data.obs that gives cluster idenitity of each cell
      summ_sc_data_args={'cluster_col': "deconv_col",
                         # select marker genes of cell types by specificity of their expression signatures
                         'selection': "cluster_specificity",
                         # specificity cutoff (1 = max, 0 = min)
                         'selection_specificity': 0.3
                        },

      train_args={'use_raw': True, # By default uses raw slots in both of the input datasets.
                  'n_iter': 40000, # Increase the number of iterations if needed (see QC below)

                  # When analysing the data that contains multiple experiments,
                  # cell2location automatically enters the mode which pools information across experiments
                  'sample_name_col': 'sample'}, # Column in sp_data.obs with experiment ID (see above)


      export_args={'path': results_folder, # path where to save results
                   'run_name_suffix': "_" + s # optinal suffix to modify the name the run
                  },

      model_kwargs={ # Prior on the number of cells, cell types and co-located groups

                    'cell_number_prior': {
                        # - N - the expected number of cells per location:
                        'cells_per_spot': 8,
                        # - A - the expected number of cell types per location:
                        'factors_per_spot': 9,
                        # - Y - the expected number of co-located cell type groups per location
                        'combs_per_spot': 5
                    },

                     # Prior beliefs on the sensitivity of spatial technology:
                    'gene_level_prior':{
                        # Prior on the mean
                        'mean': 1/2,
                        # Prior on standard deviation,
                        # a good choice of this value should be at least 2 times lower that the mean
                        'sd': 1/4
                    }
      }
    )

