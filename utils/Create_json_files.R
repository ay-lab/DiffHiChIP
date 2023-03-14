#!/usr/bin/env Rscript

##========================
## script to automatically create the .json files
##========================

##=========================
## parameters - to edit
## follow the "scripts_DiffLoop*.sh" file 
##=========================

## FitHiChIP Significant loops in WashU browser compatible files
## separate lists for two categories
## edit the file names and paths

## FitHiChIP significant WashU files for the first category
Cat1_FitHiChIP_SigFiles <- c('FitHiChIP_Cat1_R1_SigLoop.WashU.bed.gz', 'FitHiChIP_Cat1_R2_SigLoop.WashU.bed.gz')

## FitHiChIP significant WashU files for the second category
Cat2_FitHiChIP_SigFiles <- c('FitHiChIP_Cat2_R1_SigLoop.WashU.bed.gz', 'FitHiChIP_Cat2_R2_SigLoop.WashU.bed.gz')

## base output directory containing all differential HiChIP analysis results
BASEDIR <- '/home/sourya/results/DiffHiChIP/DiffLoop_CD4N_vs_CD8N'

## if ChIP-seq data is present, this is the ChIP coverage threshold to determine the non-differential 1D bins
## specified by the parameter (--CovThr)
## do not need to edit
CovThr <- 25

## two categories 
## use the same spelling as used in the "scripts_DiffLoop_Exec.sh" file 
Cat1 <- 'CD4N'
Cat2 <- 'CD8N'

## two categories - short notation names
Cat1_Short <- 'CD4N'
Cat2_Short <- 'CD8N'

## FDR and logFC thresholds for differential loops
## do not need to edit
FDR_Thr <- 0.05
logFC_Thr <- 1

## distance stratification bin size
## do not need to edit
DistStrat_BinSize <- 10000

## colors for two different categories
## do not need to edit
Category_ColorList <- c('#FF00FF', '#00FFFF')

## colors for differential ChIP bins - HD, LD or ND
## do not need to edit
ChIPBinType_ColorList <- c('#FF8000', '#00CCCC', '#009900')

## colors for exclusive loops
## Excl_OneCat, Excl_OneCat_Sig_AllRepl
## do not need to edit
Excl_Loops_ColorList <- c('#0D526F', '#006633')

## colors for different differential loop calling models
## All Loop FDR, All loop FDR with IHW, Distance stratification, per chromosome FDR, per chromosome FDR with IHW
## do not need to edit
DiffLoops_ColorList <- c('#0000FF', '#3A19C0', '#663300', '#51731F', '#3399FF')

## parameters indicating which categories of loops will be displayed

## if TRUE, plots loops exclusive to one category, and significant in at least one replicate
## default: FALSE
Plot_ExclOnecat_Sign_OneRepl <- FALSE

## if TRUE, plots loops exclusive to one category, and significant in all replicates of that category
## default: TRUE
Plot_ExclOnecat_Sign_ALLRepl <- TRUE

## if TRUE, plots differential loops derived by conventional FDR correction
## default: FALSE
Plot_DiffLoops_FDR <- FALSE

## if TRUE, plots differential loops derived by FDR + IHW correction
## default: TRUE
Plot_DiffLoops_FDR_IHW <- TRUE

## if TRUE, plots differential loops derived by distance stratification
## default: FALSE
Plot_DiffLoops_DistStat <- TRUE

## if TRUE, plots differential loops derived by per chromosome analysis
## default: FALSE
Plot_DiffLoops_PerChrAnalysis <- FALSE

##=========================
## end parameters - to edit
##=========================

##==================
## functions to add entries in the target .json file
##==================
Add_bedgraph_Info <- function(inpfilepath, label, colorval, textfile) {
	if (file.exists(inpfilepath)) {
		outtext <- paste0("\n\t {\n\t type:\"bedgraph\",\n\t url:\"",inpfilepath, "\",\n\t name:\"", label, "\",\n\t mode:\"Density\",\n\t height:40,\n\t fixedscale:{min:0,max:40},\n\t showOnHubLoad:true,\n\t options: {\n\t\t color:\"", colorval, "\",\n\t\t },\n\t }, ")
		cat(outtext, file=textfile, append=TRUE, sep="\n")
	}
}

Add_longrange_Info <- function(inpfilepath, label, colorval, textfile) {
	if (file.exists(inpfilepath)) {
		outtext <- paste0("\n\t {\n\t type:\"longrange\",\n\t url:\"",inpfilepath, "\",\n\t name:\"", label, "\",\n\t showOnHubLoad:true,\n\t options: {\n\t\t height:300,\n\t\t color:\"", colorval, "\",\n\t\t displayMode:\"arc\",\n\t\t fetchViewWindowOnly:true,\n\t\t bothAnchorsInView:true,\n\t\t scoreScale:\"fixed\",\n\t\t scoreMax:10,\n\t\t scoreMin:0 \n\t\t },\n\t }, ")
		cat(outtext, file=textfile, append=TRUE, sep="\n")
	}
}

##====================
## process the differential loops with non-differential coverage
## and using the complete set of contacts as the background
##====================
for (bckgrndtype in c('PreFilt_Bckground', 'AllLoop_Bckground')) {

	TargetBaseDir <- paste0(BASEDIR, '/NonDiff_1D_Cov_pct_', CovThr)
	TargetBaseDirLoops <- paste0(TargetBaseDir, '/', bckgrndtype)
	if (file.exists(TargetBaseDir)) {
		OutJSONFile <- paste0(BASEDIR, '/NonDiff_1D_Cov_pct_', CovThr, '_', bckgrndtype, '_WashU_TrackList.json')
		outtext <- paste0("[")
		cat(outtext, file=OutJSONFile, append=FALSE, sep="\n")

		## scaled ChIP coverages for categories 1 and 2
		Add_bedgraph_Info(paste0(TargetBaseDir, '/scaled_ChIP_Coverage_', Cat1, '.bedgraph.gz'), paste0('ChIPCov_', Cat1_Short), Category_ColorList[1], OutJSONFile)
		Add_bedgraph_Info(paste0(TargetBaseDir, '/scaled_ChIP_Coverage_', Cat2, '.bedgraph.gz'), paste0('ChIPCov_', Cat2_Short), Category_ColorList[2], OutJSONFile)

		## add bin types (HD, LD or ND)
		Add_bedgraph_Info(paste0(TargetBaseDir, '/HD_Bins.bedgraph.gz'), 'HD_Bins', ChIPBinType_ColorList[1], OutJSONFile)
		Add_bedgraph_Info(paste0(TargetBaseDir, '/LD_Bins.bedgraph.gz'), 'LD_Bins', ChIPBinType_ColorList[2], OutJSONFile)
		Add_bedgraph_Info(paste0(TargetBaseDir, '/ND_Bins.bedgraph.gz'), 'ND_Bins', ChIPBinType_ColorList[3], OutJSONFile)

		## add FitHiChIP significant loops for the first category
		for (i in 1:length(Cat1_FitHiChIP_SigFiles)) {
			Add_longrange_Info(Cat1_FitHiChIP_SigFiles[i], paste0(Cat1_Short, '_R', i), Category_ColorList[1], OutJSONFile)
		}
		## add FitHiChIP significant loops for the second category
		for (i in 1:length(Cat2_FitHiChIP_SigFiles)) {
			Add_longrange_Info(Cat2_FitHiChIP_SigFiles[i], paste0(Cat2_Short, '_R', i), Category_ColorList[2], OutJSONFile)
		}

		## add exclusive FitHiChIP loops in both categories
		## exclusive in one category, significant in at least one sample
		if (Plot_ExclOnecat_Sign_OneRepl == TRUE) {
			for (methodname in c('DESeq2', 'EdgeR_exactTest', 'EdgeR_glmLRT', 'EdgeR_glmQLFTest', 'EdgeR_glmTreat', 'EdgeR_glmTreat_glmFit')) {
				targetfile <- paste0(TargetBaseDirLoops, '/DiffModel_0_AllLoops/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Exclusive_OneCategory/DiffLoops_ALL_WashU.bed.gz')
				if (file.exists(targetfile)) {
					Add_longrange_Info(targetfile, 'Excl_OneCat', Excl_Loops_ColorList[1], OutJSONFile)
					break
				}
			}		
		}

		## add exclusive FitHiChIP loops in both categories
		## exclusive in one category, significant in all replicates
		if (Plot_ExclOnecat_Sign_ALLRepl == TRUE) {
			for (methodname in c('DESeq2', 'EdgeR_exactTest', 'EdgeR_glmLRT', 'EdgeR_glmQLFTest', 'EdgeR_glmTreat', 'EdgeR_glmTreat_glmFit')) {
				targetfile1 <- paste0(TargetBaseDirLoops, '/DiffModel_0_AllLoops/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Exclusive_OneCategory_Sig_AllRepl/DiffLoops_ALL_WashU.bed.gz')
				targetfile2 <- paste0(TargetBaseDirLoops, '/DiffModel_0_AllLoops/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Exclusive_OneCategory_Sig_AllRepl/DiffLoops_ALL_Upreg_', Cat1, '_WashU.bed.gz')
				targetfile3 <- paste0(TargetBaseDirLoops, '/DiffModel_0_AllLoops/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Exclusive_OneCategory_Sig_AllRepl/DiffLoops_ALL_Upreg_', Cat2, '_WashU.bed.gz')
				if (file.exists(targetfile1) & file.exists(targetfile2) & file.exists(targetfile3)) {
					Add_longrange_Info(targetfile1, 'Excl_OneCat_Sig_AllRepl', Excl_Loops_ColorList[2], OutJSONFile)
					## loops upregulated in category 1
					Add_longrange_Info(targetfile2, paste0('Excl_OneCat_Sig_AllRepl_UP_', Cat1_Short), Excl_Loops_ColorList[2], OutJSONFile)
					## loops upregulated in category 2
					Add_longrange_Info(targetfile3, paste0('Excl_OneCat_Sig_AllRepl_UP_', Cat2_Short), Excl_Loops_ColorList[2], OutJSONFile)
					break
				}
			}
		}

		## add conventional FDR based loops for the model 0 (all loops) 
		## and for all different computational methods (DESeq2 or EdgeR)
		if (Plot_DiffLoops_FDR == TRUE) {
			for (methodname in c('DESeq2', 'EdgeR_exactTest', 'EdgeR_glmLRT', 'EdgeR_glmQLFTest', 'EdgeR_glmTreat', 'EdgeR_glmTreat_glmFit')) {
				## differential loops
				Add_longrange_Info(paste0(TargetBaseDirLoops, '/DiffModel_0_AllLoops/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Significant_FDR/SigLoop_AtLeastOneRepl/DiffLoops_ALL_WashU.bed.gz'), paste0(methodname, '_AllLoop_FDR'), DiffLoops_ColorList[1], OutJSONFile)
				## differential loops upregulated in category 1
				Add_longrange_Info(paste0(TargetBaseDirLoops, '/DiffModel_0_AllLoops/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Significant_FDR/SigLoop_AtLeastOneRepl/DiffLoops_ALL_Upreg_', Cat1, '_WashU.bed.gz'), paste0(methodname, '_AllLoop_FDR_UP_', Cat1_Short), DiffLoops_ColorList[1], OutJSONFile)
				## differential loops upregulated in category 2
				Add_longrange_Info(paste0(TargetBaseDirLoops, '/DiffModel_0_AllLoops/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Significant_FDR/SigLoop_AtLeastOneRepl/DiffLoops_ALL_Upreg_', Cat2, '_WashU.bed.gz'), paste0(methodname, '_AllLoop_FDR_UP_', Cat2_Short), DiffLoops_ColorList[1], OutJSONFile)
			}
		}

		if (Plot_DiffLoops_FDR_IHW == TRUE) {
			## add FDR + IHW based loops for the model 0 (all loops) and for all different computational methods (DESeq2 or EdgeR)
			for (methodname in c('DESeq2', 'EdgeR_exactTest', 'EdgeR_glmLRT', 'EdgeR_glmQLFTest', 'EdgeR_glmTreat', 'EdgeR_glmTreat_glmFit')) {
				## differential loops
				Add_longrange_Info(paste0(TargetBaseDirLoops, '/DiffModel_0_AllLoops/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Significant_FDR_IHW/SigLoop_AtLeastOneRepl/DiffLoops_ALL_WashU.bed.gz'), paste0(methodname, '_AllLoop_FDR_IHW'), DiffLoops_ColorList[2], OutJSONFile)
				## differential loops upregulated in category 1
				Add_longrange_Info(paste0(TargetBaseDirLoops, '/DiffModel_0_AllLoops/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Significant_FDR_IHW/SigLoop_AtLeastOneRepl/DiffLoops_ALL_Upreg_', Cat1, '_WashU.bed.gz'), paste0(methodname, '_AllLoop_FDR_IHW_UP_', Cat1_Short), DiffLoops_ColorList[2], OutJSONFile)
				## differential loops upregulated in category 2
				Add_longrange_Info(paste0(TargetBaseDirLoops, '/DiffModel_0_AllLoops/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Significant_FDR_IHW/SigLoop_AtLeastOneRepl/DiffLoops_ALL_Upreg_', Cat2, '_WashU.bed.gz'), paste0(methodname, '_AllLoop_FDR_IHW_UP_', Cat2_Short), DiffLoops_ColorList[2], OutJSONFile)
			}			
		}

		## add loops for the model 1 (distance stratification) and for all different computational methods (DESeq2 or EdgeR)
		if (Plot_DiffLoops_DistStat == TRUE) {
			for (methodname in c('DESeq2', 'EdgeR_exactTest', 'EdgeR_glmLRT', 'EdgeR_glmQLFTest', 'EdgeR_glmTreat', 'EdgeR_glmTreat_glmFit')) {
				## differential loops
				Add_longrange_Info(paste0(TargetBaseDirLoops, '/DiffModel_1_DistStrat/Binsize_', DistStrat_BinSize, '/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Significant_FDR/SigLoop_AtLeastOneRepl/DiffLoops_ALL_WashU.bed.gz'), paste0(methodname, '_DistStat_', (DistStrat_BinSize/1000), 'kb'), DiffLoops_ColorList[3], OutJSONFile)
				## differential loops upregulated in category 1
				Add_longrange_Info(paste0(TargetBaseDirLoops, '/DiffModel_1_DistStrat/Binsize_', DistStrat_BinSize, '/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Significant_FDR/SigLoop_AtLeastOneRepl/DiffLoops_ALL_Upreg_', Cat1, '_WashU.bed.gz'), paste0(methodname, '_DistStat_', (DistStrat_BinSize/1000), 'kb_UP_', Cat1_Short), DiffLoops_ColorList[3], OutJSONFile)
				## differential loops upregulated in category 2
				Add_longrange_Info(paste0(TargetBaseDirLoops, '/DiffModel_1_DistStrat/Binsize_', DistStrat_BinSize, '/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Significant_FDR/SigLoop_AtLeastOneRepl/DiffLoops_ALL_Upreg_', Cat2, '_WashU.bed.gz'), paste0(methodname, '_DistStat_', (DistStrat_BinSize/1000), 'kb_UP_', Cat2_Short), DiffLoops_ColorList[3], OutJSONFile)
			}
		}

		## add conventional FDR based loops for the model 2 (per chromosome analysis) and for all different computational methods (DESeq2 or EdgeR)
		if (Plot_DiffLoops_PerChrAnalysis == TRUE) {
			if (Plot_DiffLoops_FDR == TRUE) {		
				for (methodname in c('DESeq2', 'EdgeR_exactTest', 'EdgeR_glmLRT', 'EdgeR_glmQLFTest', 'EdgeR_glmTreat', 'EdgeR_glmTreat_glmFit')) {
					## differential loops
					Add_longrange_Info(paste0(TargetBaseDirLoops, '/DiffModel_2_PerChr/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Significant_FDR/SigLoop_AtLeastOneRepl/DiffLoops_ALL_WashU.bed.gz'), paste0(methodname, '_Per_chr_FDR'), DiffLoops_ColorList[4], OutJSONFile)
					## differential loops upregulated in category 1
					Add_longrange_Info(paste0(TargetBaseDirLoops, '/DiffModel_2_PerChr/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Significant_FDR/SigLoop_AtLeastOneRepl/DiffLoops_ALL_Upreg_', Cat1, '_WashU.bed.gz'), paste0(methodname, '_Per_chr_FDR_UP_', Cat1_Short), DiffLoops_ColorList[4], OutJSONFile)
					## differential loops upregulated in category 2
					Add_longrange_Info(paste0(TargetBaseDirLoops, '/DiffModel_2_PerChr/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Significant_FDR/SigLoop_AtLeastOneRepl/DiffLoops_ALL_Upreg_', Cat2, '_WashU.bed.gz'), paste0(methodname, '_Per_chr_FDR_UP_', Cat2_Short), DiffLoops_ColorList[4], OutJSONFile)
				}
			}

			## add FDR + IHW based loops for the model 2 (per chromosome analysis) and for all different computational methods (DESeq2 or EdgeR)
			if (Plot_DiffLoops_FDR_IHW == TRUE) {
				for (methodname in c('DESeq2', 'EdgeR_exactTest', 'EdgeR_glmLRT', 'EdgeR_glmQLFTest', 'EdgeR_glmTreat', 'EdgeR_glmTreat_glmFit')) {
					## differential loops
					Add_longrange_Info(paste0(TargetBaseDirLoops, '/DiffModel_2_PerChr/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Significant_FDR_IHW/SigLoop_AtLeastOneRepl/DiffLoops_ALL_WashU.bed.gz'), paste0(methodname, '_Per_chr_FDR_IHW'), DiffLoops_ColorList[5], OutJSONFile)
					## differential loops upregulated in category 1
					Add_longrange_Info(paste0(TargetBaseDirLoops, '/DiffModel_2_PerChr/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Significant_FDR_IHW/SigLoop_AtLeastOneRepl/DiffLoops_ALL_Upreg_', Cat1, '_WashU.bed.gz'), paste0(methodname, '_Per_chr_FDR_IHW_UP_', Cat1_Short), DiffLoops_ColorList[5], OutJSONFile)
					## differential loops upregulated in category 2
					Add_longrange_Info(paste0(TargetBaseDirLoops, '/DiffModel_2_PerChr/FDR_', FDR_Thr, '_LOG2FC_', logFC_Thr, '/', methodname, '/Significant_FDR_IHW/SigLoop_AtLeastOneRepl/DiffLoops_ALL_Upreg_', Cat2, '_WashU.bed.gz'), paste0(methodname, '_Per_chr_FDR_IHW_UP_', Cat2_Short), DiffLoops_ColorList[5], OutJSONFile)
				}
			}
		}

		outtext <- paste0("]")
		cat(outtext, file=OutJSONFile, append=TRUE, sep="\n")

		## once the .json file is written
		## replace the entry "/mnt/BioAdHoc"
		## with the entry "https://informaticsdata.liai.org/BioAdHoc"
		system(paste0("sed -i \'s/\\/mnt\\/BioAdHoc/https:\\/\\/informaticsdata.liai.org\\/BioAdHoc/g\' ", OutJSONFile))
	}
}


