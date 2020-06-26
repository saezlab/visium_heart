
0. run_all.py ---> run step1 to step12 automatically for all samples 
1. peaks_2_putative_gene.R ---->generate putativie peak-to-gene links
2. findNN_rna.R  ----> find atac nearest rna cell based on coembeding distance
3. findNN_atac.R  ----> find atac top 50 neareast neighbour cells and aggregate these 50 cells for each cell
4. fast_cor_patients.R  ---> calculate the peak-to-gene correlation based on aggregated ATAC matrix and the neighbour RNA matrix
5. chromosome_range.py ---> based on ordered correlation file generate each chromosome range and peak related genes
6. chromosome_statistic.R ---> generate a null hypothesis peak-to-gene correlation set for each chromosome
7. test_hypothesis.R ----> test each peak-to-gene link generate p-value for the null hypothesis.
8. filter_pval.R ----> remove link padjust-value >=0.5, distance <=2000
9. 1_p2g_heatmap_celltype.R ---> heatmap for ATAC matrix only, based on significant peak-to-gene links.  Column order: put same celltype column(cell) together. Row order: pam clustering based on correlation distance.
10. 1_p2g_heatmap_celltype_reorder_tiff.R ----> heatmap for ATAC and RNA after rearrange the row cluster order manually.
11. detail_corr.py ---> add detail info onto the peak-to-gene link csv files.
12. table_add_clusters_each.R ---> assign celltype annotation for each peak-to-gene link.
