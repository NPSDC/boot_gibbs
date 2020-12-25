#!/bin/bash
#SBATCH -J run_bootstrap # Job name
#SBATCH --mail-user=npsingh@umd.edu
#SBATCH --mail-type=fail,end# Get email for begin, end, and fail
#SBATCH --time=1-18:00:00
#SBATCH --qos=large
#SBATCH --mem=36gb
#SBATCH --cpus-per-task=10

unset PYTHONPATH
FASTQ=$1
IND=$2
TYPE=$3
N=$4
INPDIR=$5
OUTDIR=$6

f1=${INPDIR}/${FASTQ}_1.fastq.gz
f2=${INPDIR}/${FASTQ}_2.fastq.gz

if [ ! -e $f1 ]
        then
                echo "File ${f1} does not exist"
                exit
fi

if [ ! -e $f2 ]
        then
                echo "File ${f2} does not exist"
                exit
fi

if [ ! -d $IND ]
        then
                echo "Invalid index directory, ${IND}"
                exit
fi

if [ $TYPE != "B" ] && [ $TYPE != "GS" ]
	then
		echo "Invalid Sampling parameters"
		exit
fi

if [ ! -d $OUTDIR ]
        then
                echo "Invalid output directory, ${OUTDIR}"
                exit
fi

dir=${OUTDIR}/${FASTQ}_${TYPE}_${N}
echo $dir
if [ ! -d $dir ]
        then
                mkdir $dir
                echo "Created directory, ${dir}"
fi

source ~/.bashrc
conda activate salmon

if [ $TYPE == "B" ]
	then
		salmon quant -p 10 -i $IND -l A --gcBias -o $dir -1 $f1 -2 $f2 --numBootstraps $N -d

else
		salmon quant -p 10 -i $IND -l A --gcBias -o $dir -1 $f1 -2 $f2 --thinningFactor 100 --numGibbsSamples $N -d
fi

