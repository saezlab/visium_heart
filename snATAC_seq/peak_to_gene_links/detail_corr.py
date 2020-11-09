from collections import defaultdict 
from pathlib import Path
import getopt
import sys

"""
python detail_corr.py
python detail_corr.py -f save/corr_pval_final_0.1.tsv
python detail_corr.py -f save/corr_pval_final_0.2.tsv
python detail_corr.py -f save/corr_pval_final_topn_10000.tsv
python detail_corr.py -f save/corr_pval_final.tsv
python detail_corr.py -f save/corr_pval_filtered.tsv
python detail_corr.py -f save/corr_pval.tsv
"""


fname = "save/corr_pval_final_0.3.tsv"

try:
    options,args = getopt.getopt(sys.argv[1:],"f:s:")
except getopt.GetoptError:
    print("Erorr Parametes")
    sys.exit()
for name,value in options:
    if name in "-f":
        fname = value
        print("processing file:", fname) 
    if name in "-s":
        sample = value
        print("On:", sample) 




p = Path(fname)
if not p.exists():
    print("No this file")
    sys.exit()


out_file = p.parent / ("detail_%s" % p.name)

fname_sample = "save/peak2gene_putative_%s.tsv" % sample

f_sample = open(fname_sample)
sample_first_line = f_sample.readline()


f = open(fname)
first_line = f.readline()
other_header = "\t".join(first_line.strip("\n").split("\t")[3:]) 

fw = open(out_file, "w")
fw.write("%s\t%s\n" % (sample_first_line.strip("\n"), other_header))


d = defaultdict(str)
for line in f_sample:
    items_sample = line.strip("\n").split("\t")
    peak_gene_sample = items_sample[0:2]
    d[tuple(peak_gene_sample)] = "\t".join(items_sample)



for line in f:
    items = line.strip("\n").split("\t")
    other_info = "\t".join(items[3:])
    peak_gene = items[0:2]
    
    out = d.get(tuple(peak_gene), "")
    if not out:
        continue            
    fw.write("%s\t%s\n" % (out, other_info))
fw.close()        

