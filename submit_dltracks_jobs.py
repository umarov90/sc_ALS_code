import os
from pathlib import Path
import sys
import subprocess
import time

# Function to count the number of jobs currently running
def count_running_jobs():
    # Execute the squeue command and get the output
    output = subprocess.check_output(["squeue", "-u", os.getlogin()], encoding='utf-8')
    # Count the number of lines in the output (excluding the header)
    return len(output.strip().split('\n')) - 1

wd = os.path.abspath("/osc-fs_home/hon-chun/analysis/ramzan/ALS/scafe/solo/")
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
        while count_running_jobs() >= 10:
            print("Maximum number of running jobs reached. Waiting for 30 seconds...")
            time.sleep(30)  # Wait for 30 seconds
        with open(f"{sd}/temp/job{i}", 'w+') as f:
            f.write("#!/bin/bash\n")
            f.write(f"#SBATCH --job-name=dl_track{i}\n")
            f.write(f"#SBATCH --output={sd}/temp/job{i}.out\n")
            f.write(f"#SBATCH --error={sd}/temp/job{i}.out\n")
            f.write(f"#SBATCH --mem=64G\n")
            f.write(f"#SBATCH --cpus-per-task=4\n")
            f.write(f"#SBATCH --partition=batch\n")
            f.write(f"python -u run_dltrack.py {directory_name}")
        os.system(f"sbatch {sd}/temp/job{i}")
        i += 1
