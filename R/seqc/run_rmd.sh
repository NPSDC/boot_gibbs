#!/bin/bash
#SBATCH -J batch_download # Job name
#SBATCH --mail-user=npsingh@umd.edu # Email for job info
#SBATCH --mail-type=fail,end# Get email for begin, end, and fail
#SBATCH --time=0-18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=32gb

module load pandoc/1.17.0.3
module load R/4.0.2/3.11
R_LIBS_USER=/cbcb/sw/RedHat-7-x86_64/users/npsingh
echo $1
Rscript create_rmd.R $1
