import os
import subprocess

sra_list = ['SRR6426203'] # Adult failing heart
sra_list = ['SRR6426182', 'SRR6426203'] # Adult non-failing heart

for sra in sra_list:
    job_name = sra
    subprocess.run(["sbatch", "-J", job_name,
                    "-o", f"./cluster_out/{job_name}.txt",
                    "-e", f"./cluster_err/{job_name}.txt",
                    "--time", "5:00:00",
                    "--mem", "180G",
                    "-c", "48",
                    "-A", "rwth0429",
                    "run.zsh", sra])
