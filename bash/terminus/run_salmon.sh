#!/bin/bash
#SBATCH -J run_salmon # Job name
#SBATCH --mail-user=npsingh@umd.edu
#SBATCH --mail-type=end,end# Get email for begin, end, and fail
#SBATCH --time=1-18:00:00
#SBATCH --qos=large
#SBATCH --mem=36gb
#SBATCH --cpus-per-task=10
TRANS_FASTA=/fs/cbcb-lab/rob/students/noor/Uncertainity/ase-sim/transcripts.fa
IND=/fs/cbcb-lab/rob/students/noor/Uncertainity/boot_gibbs/index/Droph_ens_100
SAMP_DIR=/fs/cbcb-lab/rob/students/noor/Uncertainity/ase-sim
OUT_DIR=/fs/cbcb-lab/rob/students/noor/Uncertainity/boot_gibbs/quant_output/Drosophilla
TYPE=$1
source ~/.bashrc
if [ $TYPE != "B" ] && [ $TYPE != "GS" ]
	then
		echo "Invalid Sampling parameters"
		exit
fi
###Creating index
#salmon index -t $TRANS_FASTA -i $IND -k 31 --keepDuplicates -p 10

#for f in {sim_1e6,sim_5e6}
for f in {sim_1e7,sim_5e7}
do
    COMB=${SAMP_DIR}/$f
    DIR=${OUT_DIR}/$f
    for i in {1..2}
    do
        fasta1=${COMB}/sample_0${i}_1.shuffled.fa
        fasta2=${COMB}/sample_0${i}_2.shuffled.fa
        out=${DIR}/sample_0${i}_${TYPE}
        if [ $TYPE == "GS" ]
            then
            salmon quant -p 10 -i $IND -l A --gcBias -o $out -1 $fasta1 -2 $fasta2 --thinningFactor 100 --numGibbsSamples 100 -d
        else
            salmon quant -p 10 -i $IND -l A --gcBias -o $out -1 $fasta1 -2 $fasta2 --numBootstraps 100 -d
        fi
    done
done
#salmon quant -p 10 -i $IND -l A --gcBias -o $dir -1 $f1 -2 $f2 --thinningFactor 100 --numGibbsSamples $N -d