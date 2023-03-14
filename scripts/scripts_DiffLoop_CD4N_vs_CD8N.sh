#!/bin/bash -ex
#PBS -l nodes=1:ppn=1
#PBS -l mem=200GB
#PBS -l walltime=80:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
#source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR

CodeDir='/home/sourya/Projects/DiffHiChIP/src'

CodeExec=$CodeDir'/DiffLoop_Code_FitHiChIP.R'

InpTableFile='/home/sourya/Projects/DiffHiChIP/scripts/CD4N_vs_CD8N_InputTable.txt'

ChrSizeFile='/home/sourya/genomes/chrsize/hg38.chrom.sizes'

OutDir='/home/sourya/results/DiffHiChIP/DiffLoop_CD4N_vs_CD8N'
mkdir -p $OutDir

##========================
## using All the loops as background - basically contacts present (significant or not) in at least one sample 
##========================
cd $CodeDir

##===================
## DESeq2 on the complete set of loops (model 0)
## recommended - lists both traditional FDR and FDR with IHW correction
##===================
Rscript $CodeExec --InpTable $InpTableFile --ChrSizeFile ${ChrSizeFile} --OutDir ${OutDir} --FDRThr 0.01 --DiffModel 0 --UseDESeq2 1

##===================
## DESeq2 on the complete set of loops (model 0) with pre-filtered background
## lists both traditional FDR and FDR with IHW correction
##===================
 Rscript $CodeExec --InpTable $InpTableFile --ChrSizeFile ${ChrSizeFile} --OutDir ${OutDir} --FDRThr 0.01 --PreFilt 1 --DiffModel 0 --UseDESeq2 1 

##===================
## DESeq2 using distance stratification
##===================
Rscript $CodeExec --InpTable $InpTableFile --ChrSizeFile ${ChrSizeFile} --OutDir ${OutDir} --FDRThr 0.01 --DiffModel 1 --UseDESeq2 1 

##===================
## DESeq2 using distance stratification with pre-filtered background
##===================
Rscript $CodeExec --InpTable $InpTableFile --ChrSizeFile ${ChrSizeFile} --OutDir ${OutDir} --FDRThr 0.01 --PreFilt 1 --DiffModel 1 --UseDESeq2 1

##===================
## edgeR on the complete set of loops (model 0)
## recommended - lists both traditional FDR and FDR with IHW correction
##===================
for edgermodel in 0 1 2 3 4; do
	Rscript $CodeExec --InpTable $InpTableFile --ChrSizeFile ${ChrSizeFile} --OutDir ${OutDir} --FDRThr 0.01 --DiffModel 0 --UseDESeq2 0 --EdgeRModel $edgermodel
done

##===================
## edgeR on the complete set of loops (model 0) with pre-filtered background
## lists both traditional FDR and FDR with IHW correction
##===================
for edgermodel in 0 1 2 3 4; do
	Rscript $CodeExec --InpTable $InpTableFile --ChrSizeFile ${ChrSizeFile} --OutDir ${OutDir} --FDRThr 0.01 --PreFilt 1 --DiffModel 0 --UseDESeq2 0 --EdgeRModel $edgermodel
done

##===================
## edgeR using distance stratification
##===================
for edgermodel in 0 1 2 3 4; do 
	Rscript $CodeExec --InpTable $InpTableFile --ChrSizeFile ${ChrSizeFile} --OutDir ${OutDir} --FDRThr 0.01 --DiffModel 1 --UseDESeq2 0 --EdgeRModel $edgermodel
done

##===================
## edgeR using distance stratification with pre-filtered background
##===================
for edgermodel in 0 1 2 3 4; do
	Rscript $CodeExec --InpTable $InpTableFile --ChrSizeFile ${ChrSizeFile} --OutDir ${OutDir} --FDRThr 0.01 --PreFilt 1 --DiffModel 1 --UseDESeq2 0 --EdgeRModel $edgermodel
done



