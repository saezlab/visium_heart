import os
from datetime import datetime
from multiprocessing import Process 

sample_list = ["CK166", "CK167", "CK168", #"CK169", 
            "CK170", "CK171", "CK173", "CK174"]


print("run putative", datetime.now())
lst = []
for sample in sample_list:
    f = os.system
    paramenter = "Rscript peaks_2_putative_gene.R -s  %s" % sample 
    p = Process(target=f, args=(paramenter, ))
    lst.append(p) 
    p.start()

for p in lst:
    p.join()  


## findNN 
print("run all findNN", datetime.now())
lst = []
for sample in sample_list:
    f = os.system
    paramenter = "Rscript findNN_rna.R -s %s" % sample 
    p = Process(target=f, args=(paramenter, ))
    lst.append(p) 
    p.start()

for sample in sample_list:
    f = os.system
    paramenter = "Rscript findNN_atac.R -s %s" % sample 
    p = Process(target=f, args=(paramenter, ))
    lst.append(p) 
    p.start()

for p in lst:
    p.join()  




### fast_cor_patients
print("run all fast_cor_patients", datetime.now())
lst = []
for sample in sample_list:
    f = os.system
    paramenter = "Rscript fast_cor_patients.R -s %s" % sample 
    p = Process(target=f, args=(paramenter, ))
    lst.append(p) 
    p.start()

for p in lst:
    p.join()  



print("run chromosome_range", datetime.now())
os.system("python chromosome_range.py")



print("run chromosome_statistic", datetime.now())
lst = []
for sample in sample_list:
    f = os.system
    paramenter = "Rscript chromosome_statistic.R -s %s" % sample 
    p = Process(target=f, args=(paramenter, ))
    lst.append(p) 
    p.start()

for p in lst:
    p.join()  



print("run test_hypothesis", datetime.now())
os.system("Rscript test_hypothesis.R")


print("run filtering", datetime.now())
lst = []
for sample in sample_list:
    f = os.system
    paramenter = "Rscript filter_pval.R -f 0.01 -s %s" % sample 
    p = Process(target=f, args=(paramenter, ))
    lst.append(p) 
    p.start()

for p in lst:
    p.join()  


print("run heatmap 0.05", datetime.now())
lst = []
for sample in sample_list:
    f = os.system
    paramenter = "Rscript 1_p2g_heatmap_celltype.R -g YES -d YES -r YES -s %s" % sample 
    p = Process(target=f, args=(paramenter, ))
    lst.append(p) 
    p.start()

for p in lst:
    p.join()  

print("run heatmap all 0.05", datetime.now())
lst = []
for sample in sample_list:
    f = os.system
    paramenter = "Rscript 1_p2g_heatmap_celltype_reorder_tiff.R -s %s" % sample 
    p = Process(target=f, args=(paramenter, ))
    lst.append(p) 
    p.start()

for p in lst:
    p.join()  

os.system("Rscript table_add_clusters_each.R")

