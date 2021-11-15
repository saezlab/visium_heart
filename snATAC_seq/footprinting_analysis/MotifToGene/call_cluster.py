import os
import subprocess


############################
# test
# sample_list = ["CK166"]
###########################
celltype_list = ['CM', 'Endo', 'Fib', 'Lymphoid', 'Myeloid', 'Pericyte', 'vSMCs']

for celltype in celltype_list:
    job_name = f"{celltype}"
    subprocess.run(["sbatch", "-J", job_name,
                    "-o", f"./cluster_out/{job_name}.txt",
                    "-e", f"./cluster_err/{job_name}.txt",
                    "--time", "120:00:00",
                    "--mem", "180G",
                    "run.zsh", celltype])
