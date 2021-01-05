#!/bin/bash
#SBATCH -J run_terminus # Job name
#SBATCH --mail-user=npsingh@umd.edu
#SBATCH --mail-type=end,end# Get email for begin, end, and fail
#SBATCH --time=1-18:00:00
#SBATCH --qos=large
#SBATCH --mem=36gb
#SBATCH --cpus-per-task=10

TERM_PATH=/fs/cbcb-lab/rob/students/noor/Uncertainity/terminus/target/release/terminus
SAL_DIR=/fs/cbcb-lab/rob/students/noor/Uncertainity/boot_gibbs/quant_output/Drosophilla
TERM_OUT=/fs/cbcb-lab/rob/students/noor/Uncertainity/boot_gibbs/terminus/Drosophilla

### run terminus group
# for f in {sim_1e6,sim_5e6,sim_1e7,sim_5e7}
# #for f in sim_1e6
# do
#     SF=${SAL_DIR}/$f
#     TERM_OUT_F=${TERM_OUT}/$f
#     for i in {1..2}
#     do
#         OUT=${SF}/sample_0${i}_GS
#         $TERM_PATH group -m 0.05 -d $OUT -o $TERM_OUT_F
#     done
# done

for f in {sim_1e6,sim_5e6,sim_1e7,sim_5e7}
do
    TERM_OUT_F=${TERM_OUT}/$f
    SF=${SAL_DIR}/$f
    $TERM_PATH collapse -d ${SF}/sample_01_GS ${SF}/sample_02_GS -c 0.5 -o ${TERM_OUT_F}
done