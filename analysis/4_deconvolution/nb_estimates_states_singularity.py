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
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=' + data_type + ',force_device=True'

import cell2location
import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns

# silence scanpy that prints a lot of warnings
import warnings
warnings.filterwarnings('ignore')

# read data

sc_data = sc.read("/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/integrated_data_fordeconv.h5ad")

# remove cells and genes with 0 counts everywhere
sc.pp.filter_cells(sc_data, min_genes=1)
sc.pp.filter_genes(sc_data, min_cells=1)

# calculate the mean of each gene across non-zero cells
sc_data.var['n_cells'] = (sc_data.X.toarray() > 0).sum(0)
sc_data.var['nonz_mean'] = sc_data.X.toarray().sum(0) / sc_data.var['n_cells']
sc_data.var['SYMBOL'] = sc_data.var.index.tolist()

nonz_mean_cutoff = 0.05
cell_count_cutoff = np.log10(sc_data.shape[0] * 0.0005)
cell_count_cutoff2 = np.log10(sc_data.shape[0] * 0.03)

sc_data[:,(np.array(np.log10(sc_data.var['nonz_mean']) > nonz_mean_cutoff)
         | np.array(np.log10(sc_data.var['n_cells']) > cell_count_cutoff2))
      & np.array(np.log10(sc_data.var['n_cells']) > cell_count_cutoff)].shape

# select genes based on mean expression in non-zero cells
sc_data = sc_data[:,(np.array(np.log10(sc_data.var['nonz_mean']) > nonz_mean_cutoff)
         | np.array(np.log10(sc_data.var['n_cells']) > cell_count_cutoff2))
      & np.array(np.log10(sc_data.var['n_cells']) > cell_count_cutoff)]

sc_data.raw = sc_data


from cell2location import run_regression

r, adata_sc_data = run_regression(sc_data, # input data object]

                   verbose=True, return_all=True,

                   train_args={
                    'covariate_col_names': ['deconv_col'], # column listing cell type annotation
                    'sample_name_col': 'orig_ident', # column listing sample ID for each cell
                    'tech_name_col': None,
                    'stratify_cv': 'deconv_col', # stratify cross-validation by cell type annotation
                    'n_epochs': 100, 
                    'minibatch_size': 1024,
                    'learning_rate': 0.01,
                    'use_cuda': True, # use GPU?
                    'train_proportion': 0.9, # proportion of cells in the training set (for cross-validation)
                    'l2_weight': True,  # uses defaults for the model

                    'readable_var_name_col': 'SYMBOL', 'use_raw': True},

                   model_kwargs={}, # keep defaults
                   posterior_args={}, # keep defaults

                   export_args={'path': '/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/deconvolution/c2l_statescollapsed/', # where to save results
                                'save_model': True, # save pytorch model?
                                'run_name_suffix': '_statescollapsed'})

reg_mod = r['mod']

