# Analysis of smFISH data for human myocardial infarction samples

Docker containers used to run tools for this analysis can be found under the following links:

* Mesmer: https://hub.docker.com/r/vanvalenlab/deepcell-applications
* RS-FISH: https://hub.docker.com/repository/docker/wuennemannflorian/rs_fish
* Jupyter-Scipy: https://hub.docker.com/r/jupyter/scipy-notebook

## Scripts

 * [segment_nuclei_smfish_images.sh](./segment_nuclei_smfish_images.sh)             : Used to segment nuclei from smFISH images and compute centroid positions for nuclear masks.
 * [count_spots.run_rsfish.sh](./count_spots.run_rsfish.sh)                           : Run RS-FISH on all .tif images in folder and count RNA spots.
 * [count_spots_for_CMs.Rmd](./count_spots_for_CMs.Rmd)                               : Count spots per image for quantification of NPPB and ANKDR1 signal relative to TNNT2.
 * [assign_spots_to_nuclei.Macrophages.Rmd](./assign_spots_to_nuclei.Macrophages.Rmd) : Assign spots from RS-FISH to closest nuclei positions to assign positive cell counts for markers.

