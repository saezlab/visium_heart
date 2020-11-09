import sys
import pandas as pd
from pysam import Samfile
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.Util import GenomeData

tf_file = sys.argv[1]
bam_file = sys.argv[2]
output_file = sys.argv[3]

gr_tfs = GenomicRegionSet(name="TFs")
gr_tfs.read(filename=tf_file)
gr_genes = gr_tfs.gene_association(organism="hg38")

# Fetching chromosome sizes
genome_data = GenomeData("hg38")
chrom_sizes_file_name = genome_data.get_chromosome_sizes()
chrom_sizes_file = open(chrom_sizes_file_name, "r")
chrom_sizes_dict = dict()
for chrom_sizes_entry_line in chrom_sizes_file:
	chrom_sizes_entry_vec = chrom_sizes_entry_line.strip().split("\t")
	chrom_sizes_dict[chrom_sizes_entry_vec[0]] = int(chrom_sizes_entry_vec[1])
chrom_sizes_file.close()

bam = Samfile(bam_file, "rb")

tf_list = list()
gene_list = list()
tc_list = list()

for i, r in enumerate(gr_tfs):
	tf = r.name.split(".")[-1]
	gene = gr_genes[i].name
	if gene == "." or "+" in gene or "-" in gene or ":" in gene:
		continue

	mid = (r.initial + r.final) / 2
	p1 = max(mid - 100, 0)
	p2 = min(mid + 100, chrom_sizes_dict[r.chrom])

	iter = bam.fetch(reference=r.chrom, start=p1, end=p2)
	tc = 0
	for alignment in iter:
		tc += 1
		
	tf_list.append(tf)
	gene_list.append(gene)
	tc_list.append(tc)


df = pd.DataFrame([tf_list, gene_list, tc_list])
df = df.transpose()
df.rename(columns={0: 'TF', 1: 'Gene', 2: 'TC'}, inplace=True)
df = df.groupby(['TF', 'Gene']).sum().reset_index()
df.to_csv(output_file, header=False, index=False, sep='\t')
