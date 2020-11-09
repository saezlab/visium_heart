from __future__ import print_function
# Split the bam file by cluster ID.
# Credited to https://divingintogeneticsandgenomics.rbind.io/post/split-a-10xscatac-bam-file-by-cluster/

import os
import sys
import pysam
import csv

cluster_file = sys.argv[1]
bam_file = sys.argv[2]
output_location = sys.argv[3]
output_prefix = sys.argv[4]

cluster_dict = {}
with open(cluster_file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    # skip header
    header = next(csv_reader)
    for row in csv_reader:
        cluster_dict[row[0]] = row[1]

clusters = set(x for x in cluster_dict.values())

fin = pysam.AlignmentFile(bam_file, "rb")

# open the number of bam files as the same number of clusters, and map the out file handler to the cluster id,
# write to a bam with wb
fouts_dict = {}
for cluster in clusters:
    output_filename = os.path.join(output_location, "{}_{}.bam".format(output_prefix, cluster))
    fout = pysam.AlignmentFile(output_filename, "wb", template=fin)
    fouts_dict[cluster] = fout

for read in fin:
    if not read.is_proper_pair:
        continue

    tags = read.tags
    CB_list = [x for x in tags if x[0] == "CB"]
    if CB_list:
        cell_barcode = CB_list[0][1]
    # the bam files may contain reads not in the final clustered barcodes
    # will be None if the barcode is not in the clusters.csv file
    else:
        continue
    cluster_id = cluster_dict.get(cell_barcode)
    if cluster_id:
        fouts_dict[cluster_id].write(read)

# do not forget to close the files
fin.close()
for fout in fouts_dict.values():
    fout.close()
