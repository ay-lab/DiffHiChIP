#!/bin/bash
#SBATCH -t 80:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem-per-cpu=150G
# source ~/.bashrc
source ~/.bash_profile

conda activate BASIC_R_PYTHON

CodeDir='/mnt/bioadhoc-temp/Groups/vd-ay/sourya/Projects/2022_DiffHiChIP/src'
cd $CodeDir

CodeExec=$CodeDir'/DiffLoop_Code_FitHiChIP.R'

InpTableFile='/mnt/bioadhoc-temp/Groups/vd-ay/sourya/Projects/2022_DiffHiChIP/scripts/CD4N_vs_CD8N_InputTable.txt'

ChrSizeFile='/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/chrsize/hg38_mod.chrom.sizes'

OutDir='/mnt/bioadhoc-temp/Groups/vd-ay/sourya/Projects/2022_DiffHiChIP/Results/DiffLoop_CD4N_vs_CD8N_NEW'
mkdir -p $OutDir

for PreFilt in 0 1; do
    # for DistStrat in 0 1; do
    for DistStrat in 1; do
        for Model in 0 1 2 3; do
            Rscript $CodeExec \
                --InpTable $InpTableFile \
                --OutDir ${OutDir} \
                --ChrSizeFile ${ChrSizeFile} \
                --FDRThr 0.01 \
                --Model $Model \
                --PreFilt $PreFilt \
                --DistStrat $DistStrat
        done
    done
done


