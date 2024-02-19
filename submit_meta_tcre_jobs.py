import os
from pathlib import Path
import sys

wd = os.path.abspath("/osc-fs_home/hon-chun/analysis/ramzan/ALS/scafe/count/merged_CRE_per_lib/")
sd = os.path.dirname(os.path.realpath(__file__))
Path(sd + "/temp").mkdir(parents=True, exist_ok=True)
i = 0
for entry in os.scandir(wd):
    if entry.is_dir():
        directory_name = entry.name
        directory_path = entry.path
        if not directory_name.startswith("TFHS"):
            print(f"Skipping {directory_name}")
            continue
        with open(f"{sd}/temp/job{i}", 'w+') as f:
            f.write("#!/bin/bash\n")
            f.write(f"#SBATCH --job-name=meta_sum{i}\n")
            f.write(f"#SBATCH --output={sd}/temp/job{i}.out\n")
            f.write(f"#SBATCH --error={sd}/temp/job{i}.out\n")
            f.write(f"#SBATCH --mem=64G\n")
            # f.write(f"#SBATCH --cpus-per-task=4\n")
            f.write(f"#SBATCH --partition=batch\n")
            f.write(f"python run_meta_tcre_job.py {directory_name}")
        os.system(f"sbatch {sd}/temp/job{i}")
        i += 1
