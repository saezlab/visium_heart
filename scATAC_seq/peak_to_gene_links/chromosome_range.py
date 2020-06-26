from collections import defaultdict

#f = open("data/peak_annotation_expand_newformat.tsv")

sample_list = ["CK166", "CK167", "CK168",  
               "CK170", "CK171", "CK173", "CK174"]

#sample_list = ["CK166",]


#sample_list = ["CK170"]
for sample in sample_list:
    f = open("save/corr_%s.tsv" %  sample)
    f.readline()
    last_chromosome = "xxx"

    d = defaultdict(list)
    genes_dict = defaultdict(set)

    for idx, line in enumerate(f):
        items = line.strip("\n").split("\t")
        peak = items[1]
        gene = items[2]
        chromosome = peak.split(":")[0]
        genes_dict[chromosome].add(gene)
        if last_chromosome == "xxx":
            d[chromosome].append(idx+1)
            last_chromosome = chromosome
            continue
            
        if last_chromosome != chromosome and last_chromosome != "xxx":
            d[last_chromosome].append(idx)
            d[chromosome].append(idx+1)
            last_chromosome = chromosome
        

    d[chromosome].append(idx)


    fw = open("save/chromosome_range_%s.txt" % sample, "w")
    fw.write("chromosome start end genes\n")
    for k, v in d.items():
        genes = genes_dict[k]
        fw.write("%s %d %d %s\n" %(k, v[0], v[1], ",".join(genes)))
    fw.close()
