#!/bin/bash
#SBATCH -J batch_download # Job name
#SBATCH --mail-user=npsingh@umd.edu # Email for job info
#SBATCH --mail-type=fail,end# Get email for begin, end, and fail
#SBATCH --time=0-18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=32gb

source /fs/cbcb-scratch/npsingh/virtualenv/env/bin/activate

python read_fastq.py -f $1 -d $2