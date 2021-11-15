from __future__ import print_function

import os
import glob


for filename in glob.glob("./CK*.txt"):
    base=os.path.basename(filename)
    sample = os.path.splitext(base)[0]
    job_name = "split_{}".format(sample)
    command = "sbatch -J " + job_name + " -o " + "./cluster_out/" + job_name + "_out.txt -e " + \
              "./cluster_err/" + job_name + "_err.txt -t 120:00:00 --mem 180G ./run.zsh"
    os.system(command + " " + sample)
