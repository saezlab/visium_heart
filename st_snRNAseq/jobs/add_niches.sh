#!/bin/bash

cd /Users/ricardoramirez/Dropbox/PhD/Research/mi_atlas;

Rscript ./analysis/utils/add_niche_info.R \
        --visium_folder "./processed_visium/objects/" \
        --pseudobulk_file "./processed_visium/integration/ps_integrated_slides_niches.rds";
