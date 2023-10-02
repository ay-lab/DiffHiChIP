#!/usr/bin/env Rscript

#===========================================================
## differential analysis on HiChIP and Hi-C loops
# Author: Sourya Bhattacharyya
# Vijay-Ay lab, LJI
#===========================================================

source("Header.R")
source("Functions.R")

option_list <- list(

	make_option("--InpTable",
				type = "character", 
				default = NULL, 
				help = "Input table file with loop file, covariate information. Default = NULL. Mandatory parameter."),

	make_option("--OutDir", 
				type = "character", 
				default = NULL, 
				help = "Base Output directory. Mandatory parameter."),

	make_option("--ChrSizeFile", 
				type = "character", 
				default = NULL, 
				help = "Chromosome size file corresponding to reference genome. Mandatory parameter."),

	make_option("--MidList", 
				type="character", 
				default=NULL, 
				help="Comma or colon separated list of boolean values (0/1) corresponding to the input loop files. 
						The list should have equal length as the --AllLoopList parameter.
						The order should be same as the files mentioned in the --AllLoopList parameter.
						1: loop file has midpoints of each bin 
						0: loop file has bin interval (chr, start, end). 
						Instead of a list, a single value (0 or 1) can also be provided. 
						In such a case, all files will be assumed to have the same configuration. 
						Default: list of all zeroes (means all loop files have bin interval information)."),

	make_option("--BinSize", 
				type="integer", 
				action="store", 
				default=5000, 
				help="Bin size / resolution of the input loops. 
					Used only if all the entries in --MidList option is 1, thus all loop files contain only midpoints of the bin intervals. Mandatory parameter in such a case. Default = 5000 means 5 Kb."),	

	make_option("--CCColList", 
				type="character", 
				default=NULL, 
				help="Comma or colon separated list of integers depicting the column numbers of the input loop files storing the raw contact count information. 
						The list should have equal length as the --AllLoopList parameter.
						The order should be same as the files mentioned in the --AllLoopList parameter.
						If all the input files have similar settings, user can specify only one value.
						Default = 7, same as FitHiChIP output format."),

	make_option("--QColList", 
				type="character", 
				default=NULL, 
				help="Comma or colon separated list of integers depicting the column numbers of the input loop files storing the q-value (FDR) information. 
						The list should have equal length as the --AllLoopList parameter.
						The order should be same as the files mentioned in the --AllLoopList parameter.
						If all the input files have similar settings, user can specify only one value.
						Default is the last column, same as FitHiChIP output format."),

	make_option("--FDRThr", 
				type="numeric", 
				action="store", 
				default=0.01, 
				help="FDR (q-value) threshold used for determining Hi-C or HiChIP significant loops. Default = 0.01."),	

	#==================
	# differential analysis and package related parameters

	make_option("--Model", 
				type="integer", 
				action="store", 
				default=4, 
				help="Differential analysis model.
						0: DESeq2
						1: edgeR exact Test
						2: edgeR GLM with likelihood ratio tests (glmFit + glmLRT)
						3: edgeR GLM with Quasi-likelihood (QL) F-Test (glmQLFit + glmQLFTest)
						4: Wilcoxon rank sum test on edgeR TMM normalized count
						Default = 1
						** Note: Models 0, 2, 3, are applicable if both categories have > 1 replicates."),

	make_option("--PreFilt", 
				type="integer", 
				action="store", 
				default=0, 
				help="Binary value. 
						0: all FitHiChIP loops from all samples are used as background.
						1: Union of loops having FitHiChIP FDR < 0.1 in at least one input sample are used as background. 
						Default = 0, means no filtering is done."),

  	make_option("--DiffFDRThr", 
				type="numeric", 
				action="store", 
				default=0.05, 
				help="FDR threshold for differential loops. Default = 0.05"),

	make_option("--LFCThr", 
				type="numeric", 
				action="store", 
				default=1, 
				help="log2 fold change threshold for differential loops. Default = 1."),

	make_option("--DistStrat", 
				type="integer", 
				action="store", 
				default=0, 
				help="Binary value. Possible values are: 
				[0]: No distance stratification. Use FDR and IHW corrected FDR for differential analysis.
				[1]: Use distance stratification (equal occupancy binning) and FDR. 
				Default = 0, means no distance stratification is employed."),

	make_option("--DSBinSize", 
				type="integer", 
				action="store", 
				default=10000, 
				help="Applicable if the parameter --DistStrat is 1. 
						Specifies the bin size (in bytes) employed for stratification. 
						Default = 10000 means 10 Kb is the equal occupancy bin size."),

  	make_option("--bcv", 
				type="numeric", 
				action="store", 
				default=0.4, 
				help="If edgeR is used with one replicate in at least one input category, this value is the square-root-dispersion. 
					For datasets arising from well-controlled experiments are 0.4 for human data, 
					0.1 for data on genetically identical model organisms, 
					or 0.01 for technical replicates. 
					For details, see the edgeR manual. 
					By default, the value is set as 0.4."),

	#==================
	# HiChIP related parameters
	
	make_option("--ChIPAlignFileList", 
				type="character", 
				default=NULL, 
				help="Applicable only for HiChIP data processing. 
					Comma or colon separated list of ChIP-seq alignment files. 
					Either BAM or bedgraph formatted file is supported. 
					Default = NULL. 
					User can: a) either provide two files, one for each category, or 
					b) provide ChIP seq alignment files for each sample. 
					If this parameter is empty, differential loops would be analyzed alright 
					but their underlying 1D (ChIP) differences would not."),

	make_option("--CovThr", 
				type="integer", 
				action="store", 
				default=25, 
				help="Applicable only for HiChIP data processing. 
				Threshold signifying the maximum allowed deviation of ChIP coverage between two categories for which the bins would be considered as 1D invariant (ND). 
				Default (or if user specifies value <=0 or > 100) = 25, 
				means 25% deviation is set as maximum. 
				If user chooses 50, 50% maximum ChIP seq coverage deviation would be allowed."),

	#==================
	# Miscellaneous parameters

  	make_option("--SuffixStr", 
				type="character", 
				default=NULL, 
				help="This custom string can be used to denote specific input conditions. 
				Such as DESeq2_Covariate. 
				If --Overwrite is 0, using different suffix strings would help to store outputs in different folders without re-computing majority of the differential analysis. 
				Default = NULL."),

	make_option("--Overwrite", 
				type="integer", 
				action="store", 
				default=0, 
				help="Binary variable. If 1, overwrites existing results. 
				Default = 0.")

); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#					
# ============= process the input arguments =============
# ============= check the parameters and exit condition if required ==================
# 

if (is.null(opt$ChrSizeFile)) {
	stop("ERROR !!!!!!! File containing chromosome sizes corresponding to the reference genome are not provided - check the option --ChrSizeFile \n", call.=FALSE)
	print_help(opt_parser)
}

if (is.null(opt$OutDir)) {
	print_help(opt_parser)
	stop("ERROR !!!!!!! Output directory is not provided - check the option --OutDir \n", call.=FALSE)
}

if (is.null(opt$InpTable)) {
	print_help(opt_parser)
	stop("ERROR !!!!!!! Input table file containing FitHiChIP loops, replicate information and covariates is not provided - check the option --InpTable \n", call.=FALSE) # nolint
}
InpTableData <- data.table::fread(opt$InpTable, header=TRUE)

## list of FitHiChIP loop files
AllLoopList <- as.vector(InpTableData[,1])

if (is.null(opt$MidList)) {
	## default initialization	
	MidList <- rep(DEFAULT_MIDLIST_ARG, length(AllLoopList))
} else {
	CurrList <- as.integer(unlist(strsplit(opt$MidList, "[,:]")))
	if ((length(CurrList) != length(AllLoopList)) & (length(CurrList) != 1)) {
		print_help(opt_parser)
		stop("ERROR !!!!!!! Number of entries in --MidList option should be either 1 (which will be applied to all input files) or the same as the number of input loop files. Check the option --MidList \n", call.=FALSE)
	}
	if (length(CurrList) == 1) {
		MidList <- rep(CurrList[1], length(AllLoopList))
	} else {
		MidList <- CurrList
	}
}

if (is.null(opt$CCColList)) {
	# default initialize as DEFAULT_CC_COL_ARG
	CCColList <- rep(DEFAULT_CC_COL_ARG, length(AllLoopList))
} else {
	CurrList <- as.integer(unlist(strsplit(opt$CCColList, "[,:]")))
	if ((length(CurrList) != length(AllLoopList)) & (length(CurrList) != 1)) {
		print_help(opt_parser)
		stop("ERROR !!!!!!! Number of entries in --CCColList option should be either 1 (which will be applied to all input files) or the same as the number of input loop files. Check the option --CCColList \n", call.=FALSE)
	}
	if (length(CurrList) == 1) {
		CCColList <- rep(CurrList[1], length(AllLoopList))
	} else {
		CCColList <- CurrList
	}
}

if (is.null(opt$QColList)) {
	# default initialize as DEFAULT_FDR_COL_ARG, means the last column - stores q-values in FitHiChIP
	QColList <- rep(DEFAULT_FDR_COL_ARG, length(AllLoopList))
} else {
	CurrList <- as.integer(unlist(strsplit(opt$QColList, "[,:]")))
	if ((length(CurrList) != length(AllLoopList)) & (length(CurrList) != 1)) {
		print_help(opt_parser)
		stop("ERROR !!!!!!! Number of entries in --QColList option should be either 1 (which will be applied to all input files) or the same as the number of input loop files. Check the option --QColList \n", call.=FALSE)
	}
	if (length(CurrList) == 1) {
		QColList <- rep(CurrList[1], length(AllLoopList))
	} else {
		QColList <- CurrList
	}
}

# FitHiChIP q-value threshold
if (is.null(opt$FDRThr)) {
	FDR_Th_FitHiChIP <- DEFAULT_FDR_LOOP
} else {
	FDR_Th_FitHiChIP <- as.numeric(opt$FDRThr)	
}

# DESeq2 / EdgeR FDR threshold
if (is.null(opt$DiffFDRThr)) {
	FDR_Th_DESeq <- DEFAULT_DIFFLOOP_FDR_THR
} else {
	FDR_Th_DESeq <- as.numeric(opt$DiffFDRThr)	
}

# DESeq2 / EdgeR log fold change threshold
if (is.null(opt$LFCThr)) {
	FOLD_Change_Thr <- DEFAULT_LFC_THR
} else {
	FOLD_Change_Thr <- as.numeric(opt$LFCThr)	
}

# EdgeR bcv threshold
if (is.null(opt$bcv)) {
	bcv_Thr <- DEFAULT_bcv_THR
} else {
	bcv_Thr <- as.numeric(opt$bcv)
}

# max allowed deviation for ChIP seq coverage
if (is.null(opt$CovThr)) {
	ChIP_Cov_Thr <- DEFAULT_CHIP_COVERAGE_THR
} else {
	ChIP_Cov_Thr <- (as.numeric(opt$CovThr) * 1.0) / 100
	if ((ChIP_Cov_Thr <= 0) | (ChIP_Cov_Thr > 100)) {
		ChIP_Cov_Thr <- DEFAULT_CHIP_COVERAGE_THR
	}
}

## list of loop categories
CategoryList <- unique(as.vector(InpTableData[,2]))

## replica counts for two categories
ReplicaCount <- c(length(which(InpTableData[,2] == CategoryList[1])), length(which(InpTableData[,2] == CategoryList[2])))

## replica labels for two categories
ReplicaLabels1 <- paste0('R', seq(1, ReplicaCount[1]))
ReplicaLabels2 <- paste0('R', seq(1, ReplicaCount[2]))

## ChIP-seq files - for HiChIP differential loops
if (!is.null(opt$ChIPAlignFileList)) {
	ChIPAlignFileList <- as.character(unlist(strsplit(opt$ChIPAlignFileList,"[,:]")))
	if ((length(ChIPAlignFileList) != 2) & (length(ChIPAlignFileList) != length(AllLoopList))) {
		stop("ERROR !!!!!!! Number of ChIP-seq alignment files in the --ChIPAlignFileList parameter should be either 2 (one for each category) or should be equal to the number of samples !! quit !!! \n", call.=FALSE)
	}
} else {
	ChIPAlignFileList <- c()
}

cat(sprintf("\n *** Distance stratification option (parameter --DistStrat) : %s ", opt$DistStrat))
if ((opt$DistStrat < 0) | (opt$DistStrat > 1)) {
	stop("ERROR !!!!!!! The parameter --DistStrat should have a value either 0 (no stratification) or 1 (distance stratification). Check the parameter. \n", call.=FALSE)
}

# ## DESeq2 or edgeR GLM + no replicates for at least one category
# if (((ReplicaCount[1] == 1) | (ReplicaCount[2] == 1)) & (opt$Model != 1) & (opt$Model != 4)) {
# 	stop("ERROR !!!!! At least one of the input categories have only one sample (i.e. no replicates). Only edgeR exact test model (with default test or wilcoxon rank sum test) is applicable - bur the parameter --Model is not 1 or 4 - exit !!! \n", call.=FALSE)
# }

# check the distance stratification condition
if ((opt$DistStrat == 1) & (opt$DSBinSize <= 0)) {
	stop("ERROR !!!!! distance stratification is employed (parameter --DistStrat is 1) but the corresponding bin size (parameter --DSBinSize) is not positive - Exit !!! \n", call.=FALSE)
}

# combining all the replica labels
AllRepLabels <- c(paste0(CategoryList[1], '_', ReplicaLabels1), paste0(CategoryList[2], '_', ReplicaLabels2))

# count the number of ChIP-seq alignment files
if (length(ChIPAlignFileList) == 2) {
	ChIPAlignFileCountVec <- c(1,1)
} else if (length(ChIPAlignFileList) > 0) {
	ChIPAlignFileCountVec <- ReplicaCount
} else {
	ChIPAlignFileCountVec <- c()
}

if (length(ChIPAlignFileCountVec) == 0) {
	HiC_Data_Process <- TRUE
} else {
	HiC_Data_Process <- FALSE
}

# create the base output directory
system(paste("mkdir -p", opt$OutDir))

##====================
## various output directories
##====================

## output log directory
OutLogfileDir <- paste0(opt$OutDir, '/logs')
system(paste("mkdir -p", OutLogfileDir))

## directory for storing 1D ChIP
if (HiC_Data_Process == FALSE) {
	ChIP_OutDir <- paste0(opt$OutDir, '/Diff_ChIP_1D')
	system(paste("mkdir -p", ChIP_OutDir))
}

if (opt$PreFilt == 0) {
	## directory storing differential loops using all contacts as background
	MainDir <- paste0(opt$OutDir, '/Background_AllContacts')
} else {
	## directory storing differential loops using filtered contacts as background
	MainDir <- paste0(opt$OutDir, '/Background_FilteredContacts')
}
system(paste("mkdir -p", MainDir))

## directory storing the exclusive loops for a category
## which are significant in all the replicates of that category
## and not significant in any replicate of the other category
ExclLoopDir <- paste0(opt$OutDir, '/Exclusive_OneCategory_Sig_AllRepl')
system(paste("mkdir -p", ExclLoopDir))

#
# ==========  output text file storing the parameters
#
textfile <- paste0(OutLogfileDir, '/out_', gsub("-","_",gsub(" ","_",format(Sys.time(), "%F %H-%M"))), '.log')
sink(textfile)

#
# ========== create multi-processing environment
#
ncore <- detectCores()
cat(sprintf("\n\n ==>>>> Number of cores in the system: %s ", ncore))

## register upto 4 cores
register(MulticoreParam(min(ncore, 4)))
# register(MulticoreParam(ncore))
# register(SerialParam())

# =============== get the bin size of the input loops ===============
idxList <- which(MidList == 0)
if (length(idxList) == 0) {
	# every input loop file has midpoints specified; use the input parameter bin size
	BINSIZE <- as.integer(opt$BinSize)
	if (BINSIZE == 0) {
		stop("ERROR !!!!! All input files have midpoints specified as the bin interval. But the parameter --BinSize is not specified. User should specify the correct bin size - exit !!! \n", call.=FALSE)
	}
} else {
	if (tools::file_ext(AllLoopList[idxList[1]]) == "gz") {
		BINSIZE <- as.integer(system(paste("zcat", AllLoopList[1], " | awk \'{if (NR==2) {print ($3-$2)}}\' - "), intern = TRUE))
	} else {
		BINSIZE <- as.integer(system(paste("awk \'{if (NR==2) {print ($3-$2)}}\'", AllLoopList[idxList[1]]), intern = TRUE))	
	}	
}
cat(sprintf("\n **** Bin size of Input loops : %s \n", BINSIZE))

# =============== get the complete list of chromosomes ===============
chrDF <- data.table::fread(opt$ChrSizeFile, header=F)
CHRLIST_NAMENUM <- as.vector(chrDF[,1])
cat(sprintf("\n\n Complete list of chromosomes : %s ", paste(as.vector(CHRLIST_NAMENUM), collapse=" ")))

## raw contact count and q-values
NUMFEAT <- 2

if (HiC_Data_Process == FALSE) {
	NUMCOL_SEGMENT_UNIONLOOP <- 12	
} else {
	NUMCOL_SEGMENT_UNIONLOOP <- 6
}

#
# =========== paste the configuration options within this execution ===========
#

cat(sprintf("\n ********* Listing all the input parameters ****** "))
cat(sprintf("\n ===>>> Input table contaning loop files and metadata information : %s \n ", opt$InpTable))
cat(sprintf("\n ===>>> Input loop files : %s \n ", paste(AllLoopList, collapse="\n")))
cat(sprintf("\n ===>>> Input chromosome size file : %s ", opt$ChrSizeFile))
cat(sprintf("\n ===>>> bin size of loops : %s ", BINSIZE))
cat(sprintf("\n ===>>> Midpoint of bin interval parameter : %s ", paste(MidList, collapse=" ")))
cat(sprintf("\n ===>>> Contact count columns : %s ", paste(CCColList, collapse=" ")))
cat(sprintf("\n ===>>> q-value columns : %s ", paste(QColList, collapse=" ")))
cat(sprintf("\n ===>>> FDR threshold of HiC or HiChIP loops: %s ", FDR_Th_FitHiChIP))
cat(sprintf("\n ===>>> Output directory: %s ", opt$OutDir))
cat(sprintf("\n ===>>> Input categories : \t %s ", paste(CategoryList, collapse="\t")))
cat(sprintf("\n ===>>> Count of replicates for first category : %s and for the second category : %s ", ReplicaCount[1], ReplicaCount[2]))
cat(sprintf("\n ===>>> Labels associated with replicates of category 1 : \t %s ", paste(ReplicaLabels1, collapse="\t")))
cat(sprintf("\n ===>>> Labels associated with replicates of category 2 : \t %s ", paste(ReplicaLabels2, collapse="\t")))
if (HiC_Data_Process == FALSE) {
	cat(sprintf("\n ===>>> ChIP-seq alignment files : \n %s ", paste(ChIPAlignFileList, collapse="\n")))
	cat(sprintf("\n ===>>> ChIP-seq coverage difference (fraction): %s ", ChIP_Cov_Thr))
}

cat(sprintf("\n ===>>> Differential loop analysis model employed: (0: All loops, 1: Distance stratification) : %s ", opt$DistStrat))
cat(sprintf("\n ===>>> Distance stratification bin size ? (20000 means 20 Kb) : %s ", opt$DSBinSize))
cat(sprintf("\n ===>>> Differential analysis model (0: DESeq2, 1: edgeR exact Test, 2: edgeR GLM + LRT, 3: edgeR GLM + QL F-test, 4: Wilcoxon rank sum test on edgeR TMM normalized count : %s ", opt$Model))

cat(sprintf("\n ===>>> FDR threshold for differential analysis (DESeq / EdgeR) : %s ", FDR_Th_DESeq))
cat(sprintf("\n ===>>> Log Fold change threshold for differential analysis : %s ", FOLD_Change_Thr))
cat(sprintf("\n ===>>> bcv option (constant) for differential analysis (specific to EdgeR) : %s ", bcv_Thr))

cat(sprintf("\n ===>>> Choosing background for differential analysis ====== only using loops with FitHiChIP q-value < 0.1 in at least one sample (1) or using all possible contacts (0)  : %s ", opt$PreFilt))
cat(sprintf("\n ===>>> overwrite existing output ? (1: yes, 0: no) : %s ", opt$Overwrite))

#
# =============== (if ChIP-seq data is provided)
# =============== merge ChIP-seq coverage + scale + combine ===============
#
if (HiC_Data_Process == FALSE) {		
	##
	## merge input ChIP-seq coverage values
	##
	
	MergedChIPCovFile <- paste0(ChIP_OutDir, '/ChIP_Coverage_ALL.bed')
	if ((file.exists(MergedChIPCovFile) == FALSE) | (opt$Overwrite == 1)) {
		Create_Merged_ChIP_Coverage_File(MergedChIPCovFile, ChIP_OutDir, 
										opt$ChrSizeFile, BINSIZE, ChIPAlignFileList, 
										ChIPAlignFileCountVec, CategoryList)
	}
	
	##
	## apply edgeR to get the non-differential 1D bins
	##	
	
	# this file stores the non-differential 1D bins (derived by edgeR default model)
	ChIP_1D_NonDiff_BinFile <- paste0(ChIP_OutDir, '/ChIP_Coverage_EdgeR_Default_NonSIG.bed')
	
	## this file stores the differential 1D bins (derived by edgeR default model)
	SignificantBinFile <- paste0(ChIP_OutDir, '/ChIP_Coverage_EdgeR_Default_SIG.bed')

	if ((file.exists(ChIP_1D_NonDiff_BinFile) == FALSE) | 
		(file.exists(SignificantBinFile) == FALSE) | (opt$Overwrite == 1)) {
		
		# merged ChIP-seq coverage for all input samples
		Merged_ChIPCovData <- data.table::fread(MergedChIPCovFile, header=T)
	
		# count matrix for EdgeR
		# 1. bin indices, 2. midpoints 3. their coverages (scaled) for both categories
		CountData_1D <- cbind.data.frame(
							seq(1,nrow(Merged_ChIPCovData)), 
							Merged_ChIPCovData[,1], 
							(Merged_ChIPCovData[,2] + Merged_ChIPCovData[,3])/2, 
							Merged_ChIPCovData[,4:ncol(Merged_ChIPCovData)])
		colnames(CountData_1D) <- c("Idx", "chr", "mid", colnames(Merged_ChIPCovData)[4:ncol(Merged_ChIPCovData)])
		CountDataColNames_1D <- colnames(CountData_1D)		
		if (0) {
			CountDataFile_1D <- paste0(ChIP_OutDir, '/count_matrix_1D.bed')
			write.table(CountData_1D, CountDataFile_1D, row.names=F, col.names=T, quote=F, sep="\t")
		}
	
		# appy EdgeR on ChIP coverage 
		# consider only FDR criterion - no fold change criterion is used
		ApplyEdgeR_ChIP(ChIP_OutDir, CountData_1D[, 4:ncol(CountData_1D)], 
						CategoryList, ChIPAlignFileCountVec, 'ChIP_Coverage', bcv_Thr)
		cat(sprintf("\n\n *** Performed EdgeR on ChIP coverage *** \n\n"))

		# write the input ChIP coverage along with EdgeR output for all 1D bins
		EdgeR_ResFile <- paste0(ChIP_OutDir, '/ChIP_Coverage_EdgeR_Res_temp.bed')
		AllBinFile <- paste0(ChIP_OutDir, '/ChIP_Coverage_EdgeR_Default.bed')
		system(paste("paste", MergedChIPCovFile, EdgeR_ResFile, ">", AllBinFile))

		# now extract the significant differential 1D bins
		# consider only FDR information		
		system(paste0("awk \'((NR==1) || ($NF <= ", FDR_Th_DESeq, "))\' ", AllBinFile, " > ", SignificantBinFile))

		# also extract non-significant differential 1D bins
		# since differential loops involving non-differential 1D bins are the final target	
		system(paste0("awk \'((NR==1) || ($NF >= (2 * ", FDR_Th_DESeq, ")))\' ", AllBinFile, " > ", ChIP_1D_NonDiff_BinFile))

		# delete temporary variables
		if (exists("CountData_1D")) {
			rm(CountData_1D)
		}
		if (exists("Merged_ChIPCovData")) {
			rm(Merged_ChIPCovData)
		}
		gc()	# call garbage collector routine

		## delete temporary files
		if (file.exists(EdgeR_ResFile)) {
			system(paste("rm", EdgeR_ResFile))
		}

	}	# end file exist condition

	##
	## derive the  1D bins with non-differential ChIP-seq coverage, w.r.t the specified coverage percentage threshold
	## label the 1D bins as ND (no difference), LD (low difference) and HD (high difference)
	##	
	ChIP_OutDir_CovThr <- paste0(ChIP_OutDir, '/NonDiff_1D_Cov_pct_', opt$CovThr)
	system(paste("mkdir -p", ChIP_OutDir_CovThr))

	##
	## scale the ChIP coverages according to the sequencing depth
	##
	scaled_ChIPCovFile1 <- paste0(ChIP_OutDir_CovThr, '/scaled_ChIP_Coverage_', CategoryList[1], '.bed')
	scaled_ChIPCovFile2 <- paste0(ChIP_OutDir_CovThr, '/scaled_ChIP_Coverage_', CategoryList[2], '.bed')

	MergedChIPCovFile_Scaled <- paste0(ChIP_OutDir_CovThr, '/ChIP_Coverage_ALL_Scaled_Labeled.bed')
	system(paste("cp", MergedChIPCovFile, MergedChIPCovFile_Scaled))
	
	if ((file.exists(scaled_ChIPCovFile1) == FALSE) | (file.exists(scaled_ChIPCovFile2) == FALSE) | (opt$Overwrite == 1)) {
		Create_Scaled_ChIP_Coverage_with_Label(
			MergedChIPCovFile_Scaled, ChIPAlignFileCountVec, 
			scaled_ChIPCovFile1, scaled_ChIPCovFile2, ChIP_1D_NonDiff_BinFile)
	}

	## create WashU browser tracks
	## scaled ChIP-seq coverage files
	for (scaledfilename in c(scaled_ChIPCovFile1, scaled_ChIPCovFile2)) {
		scaledfilename_washu <- gsub(".bed", ".bedgraph", scaledfilename)
		if ((file.exists(paste0(scaledfilename_washu, ".gz")) == FALSE) | (file.exists(paste0(scaledfilename_washu, ".gz.tbi")) == FALSE)) {
			system(paste0("awk \'(NR>1)\' ", scaledfilename, " | cut -f1-4 > ", scaledfilename_washu))
			if (file.exists(paste0(scaledfilename_washu, ".gz")) == FALSE) {
				system(paste("bgzip", scaledfilename_washu))
			}
			if (file.exists(paste0(scaledfilename_washu, ".gz.tbi")) == FALSE) {
				system(paste0("tabix -p bed ", scaledfilename_washu, ".gz"))
			}
		}
	}

	## separate tracks for ND, LD, and HD bins
	for (BinCat in c('ND', 'LD', 'HD')) {
		targetfile <- paste0(dirname(MergedChIPCovFile_Scaled), '/', BinCat, '_Bins.bedgraph')
		if ((file.exists(paste0(targetfile, ".gz")) == FALSE) | (file.exists(paste0(targetfile, ".gz.tbi")) == FALSE)) {
			system(paste0("awk \'{if ((NR>1) && ($NF==\"", BinCat, "\")) {print $1\"\t\"$2\"\t\"$3\"\t10\"}}\' ", MergedChIPCovFile_Scaled, " > ", targetfile))
			if (file.exists(paste0(targetfile, ".gz")) == FALSE) {
				system(paste("bgzip", targetfile))
			}		
			if (file.exists(paste0(targetfile, ".gz.tbi")) == FALSE) {
				system(paste0("tabix -p bed ", targetfile, ".gz"))	
			}		
		}
	}

}	# end condition 1D bin differential analysis

##***************
## now compute the differential loops
##***************

#
# ====================== union of loops from all input samples ======================
#
UnionLoopFile <- paste0(MainDir, '/MasterSheet_', CategoryList[1], '_', CategoryList[2], '_Loops.bed')
if ((file.exists(UnionLoopFile) == FALSE) | (opt$Overwrite == 1)) {
	if (opt$PreFilt == 1) {
		## previously we used q-value threshold of 0.01 (i.e. FitHiChIP significant in at least one input sample)
		## now we are using q-value threshold of 0.1 (slightly lenient)
		MergeLoops(AllLoopList, MidList, QColList, BINSIZE, UnionLoopFile, CHRLIST_NAMENUM, 
					FDRFilt=TRUE, FDRThr=DEFAULT_FDR_LOOP_FILT_BCK)
	} else {
		MergeLoops(AllLoopList, MidList, QColList, BINSIZE, UnionLoopFile, CHRLIST_NAMENUM, FDRFilt=FALSE)
	}
}
cat(sprintf("\n\n *** Created union of FitHiChIP loops *** \n\n"))

#
# ====================== adding features to the union of loops ======================
#
UnionLoopFeatureFile <- paste0(MainDir, '/MasterSheet_', CategoryList[1], '_', CategoryList[2], '_Loops_Features.bed')
if ((file.exists(UnionLoopFeatureFile) == FALSE) | (opt$Overwrite == 1)) {
	if (HiC_Data_Process == FALSE) {
		FillFeatureValues(UnionLoopFile, UnionLoopFeatureFile, 
							AllLoopList, BINSIZE, 
							c(scaled_ChIPCovFile1, scaled_ChIPCovFile2), 
							AllRepLabels, CategoryList, CCColList, QColList, MidList)
	} else {
		FillFeatureValues(UnionLoopFile, UnionLoopFeatureFile, 
							AllLoopList, BINSIZE, c(), 
							AllRepLabels, CategoryList, CCColList, QColList, MidList)
	}
}
cat(sprintf("\n\n *** Created master sheet of loops (with feature values) *** \n\n"))

##==============================
## read the master sheet for all loops  and their features
MasterSheetData <- data.table::fread(UnionLoopFeatureFile, header=T, sep="\t", stringsAsFactors=F)

# find the minimum and maximum loop distance values
MinLoopDist <- min(abs(MasterSheetData[, 5] - MasterSheetData[, 2]))
MaxLoopDist <- max(abs(MasterSheetData[, 5] - MasterSheetData[, 2]))
cat(sprintf("\n\n ********* minimum loop distance computed from MasterSheetData : %s \n maximum loop distance computed from MasterSheetData : %s ******* \n", MinLoopDist, MaxLoopDist))

#
# ============ get the feature columns for this master sheet ============
#

# list of columns in the master sheet data storing the raw contact counts
RawCC_ColList <- NUMCOL_SEGMENT_UNIONLOOP + ((seq(1, length(AllRepLabels)) - 1) * NUMFEAT) + 1	

# list of columns in the master sheet data storing the q-value
QVal_ColList <- RawCC_ColList + 1

#
# =========== sequence depth equalized raw contact count per loop
# =========== applicable when both input categories have exactly one replicate
# =========== currently not used - sourya
#
if (0) {
	if ((ReplicaCount[1] == 1) & (ReplicaCount[2] == 1)) {	
		# sum of contact counts for all loops in category 1
		sumCC_Cat1 <- sum(MasterSheetData[, RawCC_ColList[1]])
		# sum of contact counts for all loops in category 2
		sumCC_Cat2 <- sum(MasterSheetData[, RawCC_ColList[2]])	
		# normalization factor
		norm_factor <- ((max(sumCC_Cat1, sumCC_Cat2) * 1.0) / min(sumCC_Cat1, sumCC_Cat2))
		if (sumCC_Cat1 > sumCC_Cat2) {
			# up-scale the raw contact counts of second category
			MasterSheetData[, RawCC_ColList[2]] <- as.integer(MasterSheetData[, RawCC_ColList[2]] * norm_factor)
		} else {
			# up-scale the raw contact counts of first category
			MasterSheetData[, RawCC_ColList[1]] <- as.integer(MasterSheetData[, RawCC_ColList[1]] * norm_factor)
		}
		# write down the modified master sheet
		write.table(MasterSheetData, paste0(MainDir, '/MasterSheet_', CategoryList[1], '_', CategoryList[2], '_Loops_Scaled_CC.bed'), row.names=F, col.names=T, sep="\t", quote=F, append=F)
	}
}

#
# ================== find the list of chromosomes in the master sheet ==================
#
ChrList_MasterSheet <- unique(MasterSheetData[,1])
cat(sprintf("\n **** Chromsome list in master sheet : %s \n", paste(as.vector(ChrList_MasterSheet), collapse=" ")))

#
# ============ get the exclusive loops in either categories ============
#

## number of significant replicates
SigVec_Cat1 <- GetCntNumRepVec(MasterSheetData, 
					RawCC_ColList[1:ReplicaCount[1]], 
					QVal_ColList[1:ReplicaCount[1]], 
					FDR_Th_FitHiChIP)
SigVec_Cat2 <- GetCntNumRepVec(MasterSheetData, 
					RawCC_ColList[(ReplicaCount[1]+1):(ReplicaCount[1]+ReplicaCount[2])], 
					QVal_ColList[(ReplicaCount[1]+1):(ReplicaCount[1]+ReplicaCount[2])], 
					FDR_Th_FitHiChIP)
SigReplDF <- cbind.data.frame(MasterSheetData, data.frame(f1=SigVec_Cat1, f2=SigVec_Cat2))
colnames(SigReplDF) <- c(colnames(MasterSheetData), c(paste0(CategoryList[1], '_SigRepl'), paste0(CategoryList[2], '_SigRepl')))

## exclusive loops in either categories
idx <- which(((SigReplDF[, (ncol(SigReplDF)-1)] == ReplicaCount[1]) & (SigReplDF[, ncol(SigReplDF)] == 0)) | ((SigReplDF[, (ncol(SigReplDF)-1)] == 0) & (SigReplDF[, ncol(SigReplDF)] == ReplicaCount[2])))
if (length(idx) > 0) {
	CompleteLoopFile <- paste0(ExclLoopDir, '/DiffLoops_ALL.bed')
	write.table(SigReplDF[idx, ], CompleteLoopFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
	Categorize_ExclusiveLoops(ExclLoopDir, CompleteLoopFile, CategoryList, QVal_ColList, HiC_Data_Process, opt$CovThr)
}

#
# ========== create output directory for the current differential model
#
if (opt$DistStrat == 0) {
	DiffLoopDir <- paste0(MainDir, '/No_DistStrat')
} else {
	DiffLoopDir <- paste0(MainDir, '/With_DistStrat_', opt$DSBinSize)
}

##==========
## output folder names
##==========
if (opt$Model == 0) {
	DiffLoopDir <- paste0(DiffLoopDir, '/DESeq2')
} else if (opt$Model == 1) {
	DiffLoopDir <- paste0(DiffLoopDir, '/EdgeR_exactTest')
} else if (opt$Model == 2) {
	DiffLoopDir <- paste0(DiffLoopDir, '/EdgeR_glmLRT')
} else if (opt$Model == 3) {
	DiffLoopDir <- paste0(DiffLoopDir, '/EdgeR_glmQLFTest')
} else if (opt$Model == 4) {
	DiffLoopDir <- paste0(DiffLoopDir, '/Wilcox_rank_sum')
}

if (!is.null(opt$SuffixStr)) {
	DiffLoopDir <- paste0(DiffLoopDir, '_', opt$SuffixStr)
}
system(paste("mkdir -p", DiffLoopDir))

#
# ================== apply DESeq for differential analysis ==================
#
if (opt$Model == 0) {

	## count matrix file
	CountDataFile <- paste0(DiffLoopDir, '/count_matrix.bed')

	## DESeq2 sample condition file
	SampleInfoFile <- paste0(DiffLoopDir, '/condition_file.bed')

	## file to contain DESeq2 output
	DESeq_ResFile <- paste0(DiffLoopDir, '/DESEQ_Results.bed')

	## file to contain complete DESeq2 output - final results
	CompleteLoopFile <- paste0(DiffLoopDir, '/FINAL_DESEQ_Results.bed')

	## stores count data file for the current distance value
	## applicable for the opt$DistStrat = 1
	CountDataFile_CurrDist <- paste0(DiffLoopDir, '/count_matrix_CurrDist.bed')			

	## check if the DESeq2 results are already computed
	if (file.exists(CompleteLoopFile) == FALSE) {

		if (opt$DistStrat == 0) {
			
			##=============
			# no distance stratification - process all loops at once 
			##=============
			# create count matrix
			# Create_DESeq2_Compatible_Count_Matrix(MasterSheetData, RawCC_ColList, CountDataFile)
			write.table(MasterSheetData[, RawCC_ColList], CountDataFile, row.names=F, col.names=T, quote=F, sep="\t")
					
			# add pseudo count to all counts of 0 in the count matrix
			# this is because, DESEQ2 requires at least one gene whose all elements are non zero
			# check https://help.galaxyproject.org/t/error-with-deseq2-every-gene-contains-at-least-one-zero/564			
			CountData <- data.table::fread(CountDataFile, header=T)						
			CountData[CountData == 0] <- 1

			## used to label the plot files
			if (is.null(opt$SuffixStr)) {
				suffixStr <- 'AllLoops'
			} else {				
				suffixStr <- paste0('AllLoops_', opt$SuffixStr)
			}

			## call DESEQ2 for this data
			Perform_DESeq2(DiffLoopDir, opt$InpTable, CountData, 
							SampleInfoFile, DESeq_ResFile, opt$Model, opt$DistStrat, 
							suffixStr, FDR_Th_DESeq)

			## read the DESEQ2 results
			DESEQRes <- data.table::fread(DESeq_ResFile, header=T)
			
			## append the mastersheet data and the significance of individual samples
			SigVec_Cat1 <- GetCntNumRepVec(MasterSheetData, 
											RawCC_ColList[1:ReplicaCount[1]], 
											QVal_ColList[1:ReplicaCount[1]], 
											FDR_Th_FitHiChIP)
			SigVec_Cat2 <- GetCntNumRepVec(MasterSheetData, 
								RawCC_ColList[(ReplicaCount[1]+1):(ReplicaCount[1]+ReplicaCount[2])], QVal_ColList[(ReplicaCount[1]+1):(ReplicaCount[1]+ReplicaCount[2])], FDR_Th_FitHiChIP)
			SigReplDF <- data.frame(f1=SigVec_Cat1, f2=SigVec_Cat2)
			colnames(SigReplDF) <- c(paste0(CategoryList[1], '_SigRepl'), paste0(CategoryList[2], '_SigRepl'))
			OutDESeqAllRes <- cbind.data.frame(MasterSheetData, DESEQRes, SigReplDF)
			write.table(OutDESeqAllRes, CompleteLoopFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
		
		} else {			

			##=============
			# distance stratification
			##=============
			bool_DESEQ_Out_Dump <- FALSE

			NumDistBins <- ceiling((MaxLoopDist - MinLoopDist) / opt$DSBinSize)
			cat(sprintf("\n\n ********* MinLoopDist : %s MaxLoopDist : %s  DSBinSize : %s NumDistBins : %s ****** \n\n", MinLoopDist, MaxLoopDist, opt$DSBinSize, NumDistBins))

			#======================
			# for distance stratification,
			# equal ocupancy binning is employed				
			#======================

			# average number of locus pairs per bin
			avgloci <- (nrow(MasterSheetData) / NumDistBins)
			cat(sprintf("\n total number of elements in master sheet : %s average number of elements per distance stratify bins : %s ", nrow(MasterSheetData), avgloci))

			carry_fwd_idx <- c()
			for (binidx in 1:NumDistBins) {
				startDist <- MinLoopDist + (binidx - 1) * opt$DSBinSize
				endDist <- startDist + opt$DSBinSize
				cat(sprintf("\n distance stratification - processing bin idx : %s startDist : %s endDist : %s ", binidx, startDist, endDist))
				# extract master sheet data for the current distance
				loopidx <- which((abs(MasterSheetData[,5] - MasterSheetData[,2]) >= startDist) & (abs(MasterSheetData[,5] - MasterSheetData[,2]) < endDist))
				cat(sprintf("\n --->> number of elements for DESeq2 processing in the current bin : %s number of carry forward samples: %s ", length(loopidx), length(carry_fwd_idx)))
				loopidx <- c(loopidx, carry_fwd_idx)
				# if the number of loops < avgloci, continue	
				if ((length(loopidx) < avgloci) & (binidx <= NumDistBins)) {
					carry_fwd_idx <- loopidx
					cat(sprintf("\n *********** total sample count < %s --- using it as a carry forward ", avgloci))
					next
				}
				# extract master sheet and count data for the current distance
				MasterSheetData_CurrDist <- MasterSheetData[loopidx, ]
				
				# create count matrix
				# Create_DESeq2_Compatible_Count_Matrix(MasterSheetData_CurrDist, RawCC_ColList, CountDataFile_CurrDist)
				write.table(MasterSheetData_CurrDist[, RawCC_ColList], CountDataFile_CurrDist, row.names=F, col.names=T, quote=F, sep="\t")

				CountData <- data.table::fread(CountDataFile_CurrDist, header=T)
				# add pseudo count to all counts of 0 in the count matrix
				# this is because, DESEQ2 requires at least one gene whose all elements are non zero
				# check https://help.galaxyproject.org/t/error-with-deseq2-every-gene-contains-at-least-one-zero/564
				CountData[CountData == 0] <- 1					

				## used to label the plot files
				mindist <- min(MasterSheetData_CurrDist[, 5] - MasterSheetData_CurrDist[, 2])
				maxdist <- max(MasterSheetData_CurrDist[, 5] - MasterSheetData_CurrDist[, 2])
				if (maxdist == mindist) {
					suffixStr <- paste0('dist_', mindist)
				} else {
					suffixStr <- paste0('dist_', mindist, '_to_', maxdist)
				}

				# function to call DESEQ2 for this distance value
				Perform_DESeq2(DiffLoopDir, opt$InpTable, CountData, 
								SampleInfoFile, DESeq_ResFile, opt$Model, opt$DistStrat, 
								suffixStr, FDR_Th_DESeq)

				if (file.exists(DESeq_ResFile)) {
					## read the DESEQ2 results
					DESEQRes <- data.table::fread(DESeq_ResFile, header=T, sep="\t", stringsAsFactors=F)

					## append the mastersheet data at the beginning
					OutDESeqAllRes <- cbind.data.frame(MasterSheetData_CurrDist, DESEQRes)
					write.table(OutDESeqAllRes, CompleteLoopFile, row.names=F, col.names=!bool_DESEQ_Out_Dump, sep="\t", quote=F, append=bool_DESEQ_Out_Dump)

					# update the DESEQ output dump condition after processing one distance value
					bool_DESEQ_Out_Dump <- TRUE

					# # reset the carry_fwd_idx
					# carry_fwd_idx <- c()

				} else {
					cat(sprintf("\n\n DESeq2 model could not converge for this input !!! \n\n "))

					# ## update the carry forward samples
					# carry_fwd_idx <- loopidx
					# cat(sprintf("\n\n Updated number of carry forward samples : %s ", length(carry_fwd_idx)))
				}

				# reset the carry_fwd_idx
				carry_fwd_idx <- c()				

			}	# end bin loop

			##==== once all bins are processed, employ global FDR correction
			##==== for the complete results
			OutDESeqAllRes <- data.table::fread(CompleteLoopFile, header=T)
			OutDESeqAllRes$FDR_BH <- p.adjust(OutDESeqAllRes$pvalue, method = "BH")

			##=== also append the number of significant replicates
			SigVec_Cat1 <- GetCntNumRepVec(OutDESeqAllRes, 
											RawCC_ColList[1:ReplicaCount[1]], 
											QVal_ColList[1:ReplicaCount[1]], 
											FDR_Th_FitHiChIP)
			SigVec_Cat2 <- GetCntNumRepVec(OutDESeqAllRes, 
						RawCC_ColList[(ReplicaCount[1]+1):(ReplicaCount[1]+ReplicaCount[2])], QVal_ColList[(ReplicaCount[1]+1):(ReplicaCount[1]+ReplicaCount[2])], FDR_Th_FitHiChIP)
			SigReplDF <- data.frame(f1=SigVec_Cat1, f2=SigVec_Cat2)
			colnames(SigReplDF) <- c(paste0(CategoryList[1], '_SigRepl'), paste0(CategoryList[2], '_SigRepl'))			

			##===== write the updated results
			OutDESeqAllRes <- cbind.data.frame(OutDESeqAllRes, SigReplDF)
			write.table(OutDESeqAllRes, CompleteLoopFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)

		}	# end if - distance stratification condition

		##====== remove temporary files and objects
		if (exists("CountData")) {
			rm(CountData)
		}
		if (exists("MasterSheetData")) {
			rm(MasterSheetData)
		}
		if (exists("MasterSheetData_CurrDist")) {
			rm(MasterSheetData_CurrDist)
		}
		if (exists("SigVec_Cat1")) {
			rm(SigVec_Cat1)
		}
		if (exists("SigVec_Cat2")) {
			rm(SigVec_Cat2)
		}	
		if (exists("SigReplDF")) {
			rm(SigReplDF)
		}		
		if (exists("DESEQRes")) {
			rm(DESEQRes)
		}	
		if (exists("OutDESeqAllRes")) {
			rm(OutDESeqAllRes)
		}
		if (file.exists(CountDataFile))	{
			system(paste("rm", CountDataFile))
		}
		# if (file.exists(SampleInfoFile))	{
		# 	system(paste("rm", SampleInfoFile))
		# }
		if (file.exists(DESeq_ResFile))	{
			system(paste("rm", DESeq_ResFile))
		}
		if (file.exists(CountDataFile_CurrDist))	{
			system(paste("rm", CountDataFile_CurrDist))
		}
		gc()	# call garbage collector routine

	}	# end if DESeq2 results already exist

	##================
	## extract significant loops
	##================
	if (opt$DistStrat == 0) {
		
		## No distance stratification
		## two sets of significant loops
		## one for the conventional FDR, one for the IHW corrected FDR
		## loops significant in at least one input sample are reported
		
		##
		##========= set 1: conventional FDR
		##
		# DESeq2_SigDir <- paste0(DiffLoopDir, '/DiffLoop_Sig_FDR_BH')
		DESeq2_SigDir <- paste0(DiffLoopDir, '/DiffLoop_BH_FDR_', FDR_Th_DESeq, '_LOG2FC_', FOLD_Change_Thr)
		system(paste("mkdir -p", DESeq2_SigDir))
		OutLoopFile <- paste0(DESeq2_SigDir, '/DiffLoops_ALL.bed')
		system(paste0("awk -F\'[\t]\' \'function abs(v) {return v < 0 ? -v : v} {if ((NR==1) || (($(NF-3) <= ", FDR_Th_DESeq, ") && (abs($(NF-6)) >= ", FOLD_Change_Thr, "))) {print $0}}\' ", CompleteLoopFile, " | awk -F\'[\t]\' \'((NR==1) || ($(NF-1)>0) || ($NF>0))\' - > ", OutLoopFile))
		
		## categorize these DESeq2 significant loops 
		## according to differential / non-differential ChIP coverage		
		Categorize_DiffLoops(DESeq2_SigDir, OutLoopFile, CategoryList, QVal_ColList, opt$Model, HiC_Data_Process, opt$CovThr)	
		
		##
		##========= set 2: FDR computed by IHW
		##		
		# DESeq2_SigDir <- paste0(DiffLoopDir, '/DiffLoop_Sig_FDR_IHW')
		DESeq2_SigDir <- paste0(DiffLoopDir, '/DiffLoop_IHW_FDR_', FDR_Th_DESeq, '_LOG2FC_', FOLD_Change_Thr)
		system(paste("mkdir -p", DESeq2_SigDir))
		OutLoopFile <- paste0(DESeq2_SigDir, '/DiffLoops_ALL.bed')
		system(paste0("awk -F\'[\t]\' \'function abs(v) {return v < 0 ? -v : v} {if ((NR==1) || (($(NF-2) <= ", FDR_Th_DESeq, ") && (abs($(NF-6)) >= ", FOLD_Change_Thr, "))) {print $0}}\' ", CompleteLoopFile, " | awk -F\'[\t]\' \'((NR==1) || ($(NF-1)>0) || ($NF>0))\' - > ", OutLoopFile))

		## categorize these DESeq2 significant loops 
		## according to differential / non-differential ChIP coverage
		Categorize_DiffLoops(DESeq2_SigDir, OutLoopFile, CategoryList, QVal_ColList, opt$Model, HiC_Data_Process, opt$CovThr)	
		
	} else {

		##
		##========= set 1: computed by using FDR_BH column (global FDR)
		##
		# DESeq2_SigDir <- paste0(DiffLoopDir, '/DiffLoop_Sig_FDR_BH')
		DESeq2_SigDir <- paste0(DiffLoopDir, '/DiffLoop_BH_FDR_', FDR_Th_DESeq, '_LOG2FC_', FOLD_Change_Thr)
		system(paste("mkdir -p", DESeq2_SigDir))
		OutLoopFile <- paste0(DESeq2_SigDir, '/DiffLoops_ALL.bed')		
		system(paste0("awk -F\'[\t]\' \'function abs(v) {return v < 0 ? -v : v} {if ((NR==1) || (($(NF-2) <= ", FDR_Th_DESeq, ") && (abs($(NF-6)) >= ", FOLD_Change_Thr, "))) {print $0}}\' ", CompleteLoopFile, " | awk -F\'[\t]\' \'((NR==1) || ($(NF-1)>0) || ($NF>0))\' - > ", OutLoopFile))

		## categorize these DESeq2 significant loops 
		## according to differential / non-differential ChIP coverage
		Categorize_DiffLoops(DESeq2_SigDir, OutLoopFile, CategoryList, QVal_ColList, opt$Model, HiC_Data_Process, opt$CovThr)
		
	}				

	cat(sprintf("\n\n *** Applied DESeq2 for differential loop finding *** \n\n"))
	
}	# end applying DESeq2 condition

#
# ================== apply EdgeR for differential analysis ==================
#
if ((opt$Model >= 1) & (opt$Model <= 4)) {

	## count matrix file
	CountDataFile <- paste0(DiffLoopDir, '/count_matrix.bed')

	## edgeR sample condition file
	SampleInfoFile <- paste0(DiffLoopDir, '/condition_file.bed')	

	## file to contain EdgeR output
	EdgeR_ResFile <- paste0(DiffLoopDir, '/edgeR_Results.bed')

	## file to contain complete DESeq2 output - final results
	CompleteLoopFile <- paste0(DiffLoopDir, '/FINAL_edgeR_Results.bed')

	## stores count data file for the current distance value
	## applicable for the Distance Stratification
	CountDataFile_CurrDist <- paste0(DiffLoopDir, '/count_matrix_CurrDist.bed')	

	## check if the EdgeR results are already computed
	if (file.exists(CompleteLoopFile) == FALSE) {	

		if (opt$DistStrat == 0) {
			
			##=============
			# No Distance stratification
			##=============
			# create count matrix
			Create_EdgeR_Compatible_Count_Matrix(MasterSheetData, RawCC_ColList, DiffLoopDir, CountDataFile)			

			## read the count data
			CountData <- data.table::fread(CountDataFile, header=T)

			## used to label the plot files
			if (is.null(opt$SuffixStr)) {
				suffixStr <- 'AllLoops'
			} else {
				suffixStr <- paste0('AllLoops_', opt$SuffixStr)
			}

			## call EdgeR for this data
			Perform_EdgeR(CountData, DiffLoopDir, opt$InpTable, 
							SampleInfoFile, EdgeR_ResFile, 
							opt$Model, opt$DistStrat, suffixStr, FDR_Th_DESeq, bcv_Thr)
			
			## read the edgeR results
			EdgeRRes <- data.table::fread(EdgeR_ResFile, header=T)
			
			## append the mastersheet data and the significance of individual samples
			SigVec_Cat1 <- GetCntNumRepVec(MasterSheetData, 
											RawCC_ColList[1:ReplicaCount[1]], 
											QVal_ColList[1:ReplicaCount[1]], 
											FDR_Th_FitHiChIP)
			SigVec_Cat2 <- GetCntNumRepVec(MasterSheetData, 
							RawCC_ColList[(ReplicaCount[1]+1):(ReplicaCount[1]+ReplicaCount[2])], QVal_ColList[(ReplicaCount[1]+1):(ReplicaCount[1]+ReplicaCount[2])], FDR_Th_FitHiChIP)
			SigReplDF <- data.frame(f1=SigVec_Cat1, f2=SigVec_Cat2)
			colnames(SigReplDF) <- c(paste0(CategoryList[1], '_SigRepl'), paste0(CategoryList[2], '_SigRepl'))
			
			## final edgeR results
			OutEdgeRAllRes <- cbind.data.frame(MasterSheetData, EdgeRRes, SigReplDF)
			write.table(OutEdgeRAllRes, CompleteLoopFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
		
		} else {

			##=============
			# distance stratification
			##=============
			bool_EdgeR_Out_Dump <- FALSE

			# distance stratification, w.r.t "DSBinSize"
			NumDistBins <- ceiling((MaxLoopDist - MinLoopDist) / opt$DSBinSize)
			cat(sprintf("\n\n ********* MinLoopDist : %s MaxLoopDist : %s  DSBinSize : %s NumDistBins : %s ****** \n\n", MinLoopDist, MaxLoopDist, opt$DSBinSize, NumDistBins))

			#======================
			# approach 2 - recommended
			# equal ocupancy binning is employed				
			#======================
				
			# average number of locus pairs per bin
			avgloci <- (nrow(MasterSheetData) / NumDistBins)
			cat(sprintf("\n total number of elements in master sheet : %s average number of elements per distance stratify bins : %s ", nrow(MasterSheetData), avgloci))

			carry_fwd_idx <- c()
			for (binidx in 1:NumDistBins) {
				startDist <- MinLoopDist + (binidx - 1) * opt$DSBinSize
				endDist <- startDist + opt$DSBinSize
				cat(sprintf("\n distance stratification - processing bin idx : %s startDist : %s endDist : %s ", binidx, startDist, endDist))
				# extract master sheet data for the current distance
				loopidx <- which((abs(MasterSheetData[,5] - MasterSheetData[,2]) >= startDist) & (abs(MasterSheetData[,5] - MasterSheetData[,2]) < endDist))
				cat(sprintf("\n --->> number of elements for EdgeR processing in the current bin : %s number of carry forward samples: %s ", length(loopidx), length(carry_fwd_idx)))
				loopidx <- c(loopidx, carry_fwd_idx)

				# if the number of loops < avgloci, continue	
				if ((length(loopidx) < avgloci) & (binidx <= NumDistBins)) {
					carry_fwd_idx <- loopidx
					cat(sprintf("\n *********** total sample count < %s --- using it as a carry forward ", avgloci))
					next
				}

				# extract master sheet and count data for the current distance
				MasterSheetData_CurrDist <- MasterSheetData[loopidx, ]

				## create count file compatible to this mastersheet
				Create_EdgeR_Compatible_Count_Matrix(MasterSheetData_CurrDist, RawCC_ColList, DiffLoopDir, CountDataFile_CurrDist)
				
				# read the count data matrix for EdgeR based processing		
				CountData <- data.table::fread(CountDataFile_CurrDist, header=T)	

				## used to label the plot files
				mindist <- min(MasterSheetData_CurrDist[, 5] - MasterSheetData_CurrDist[, 2])
				maxdist <- max(MasterSheetData_CurrDist[, 5] - MasterSheetData_CurrDist[, 2])
				if (maxdist == mindist) {
					suffixStr <- paste0('dist_', mindist)
				} else {
					suffixStr <- paste0('dist_', mindist, '_to_', maxdist)
				}

				## call EdgeR for this data
				Perform_EdgeR(CountData, DiffLoopDir, opt$InpTable, 
								SampleInfoFile, EdgeR_ResFile, 
								opt$Model, opt$DistStrat, suffixStr, FDR_Th_DESeq, bcv_Thr)

				if (file.exists(EdgeR_ResFile)) {
					
					## read the edgeR results
					EdgeRRes <- data.table::fread(EdgeR_ResFile, header=T)

					## append the mastersheet data at the beginning
					OutEdgeRAllRes <- cbind.data.frame(MasterSheetData_CurrDist, EdgeRRes)
					write.table(OutEdgeRAllRes, CompleteLoopFile, row.names=F, col.names=!bool_EdgeR_Out_Dump, sep="\t", quote=F, append=bool_EdgeR_Out_Dump)

					# update the EdgeR output dump condition after processing one distance value
					bool_EdgeR_Out_Dump <- TRUE

					# # reset the carry_fwd_idx
					# carry_fwd_idx <- c()

				} else {
					cat(sprintf("\n\n edgeR model could not converge for this input !!! \n\n "))

					# ## update the carry forward samples
					# carry_fwd_idx <- loopidx
					# cat(sprintf("\n\n Updated number of carry forward samples : %s ", length(carry_fwd_idx)))
				}

				# reset the carry_fwd_idx
				carry_fwd_idx <- c()				

			}	# end bin loop

			##==== once all bins are processed, employ global FDR correction
			##==== for the complete results
			OutEdgeRAllRes <- data.table::fread(CompleteLoopFile, header=T)
			OutEdgeRAllRes$FDR_BH <- p.adjust(OutEdgeRAllRes$pvalue, method = "BH")

			##=== also append the number of significant replicates
			SigVec_Cat1 <- GetCntNumRepVec(OutEdgeRAllRes, 
											RawCC_ColList[1:ReplicaCount[1]], 
											QVal_ColList[1:ReplicaCount[1]], 
											FDR_Th_FitHiChIP)
			SigVec_Cat2 <- GetCntNumRepVec(OutEdgeRAllRes, 
								RawCC_ColList[(ReplicaCount[1]+1):(ReplicaCount[1]+ReplicaCount[2])], QVal_ColList[(ReplicaCount[1]+1):(ReplicaCount[1]+ReplicaCount[2])], FDR_Th_FitHiChIP)
			SigReplDF <- data.frame(f1=SigVec_Cat1, f2=SigVec_Cat2)
			colnames(SigReplDF) <- c(paste0(CategoryList[1], '_SigRepl'), paste0(CategoryList[2], '_SigRepl'))			

			##===== write the updated results
			OutEdgeRAllRes <- cbind.data.frame(OutEdgeRAllRes, SigReplDF)
			write.table(OutEdgeRAllRes, CompleteLoopFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)		

		}	# end if - distance stratification

		##====== remove temporary files and objects		
		if (exists("MasterSheetData")) {
			rm(MasterSheetData)
		}
		if (exists("MasterSheetData_CurrDist")) {
			rm(MasterSheetData_CurrDist)
		}
		if (exists("SigVec_Cat1")) {
			rm(SigVec_Cat1)
		}
		if (exists("SigVec_Cat2")) {
			rm(SigVec_Cat2)
		}	
		if (exists("SigReplDF")) {
			rm(SigReplDF)
		}		
		if (exists("EdgeRRes")) {
			rm(EdgeRRes)
		}	
		if (exists("OutEdgeRAllRes")) {
			rm(OutEdgeRAllRes)
		}
		if (file.exists(CountDataFile))	{
			system(paste("rm", CountDataFile))
		}
		if (file.exists(EdgeR_ResFile))	{
			system(paste("rm", EdgeR_ResFile))
		}
		if (file.exists(CountDataFile_CurrDist))	{
			system(paste("rm", CountDataFile_CurrDist))
		}
		gc()	# call garbage collector routine

	}	# end edgeR results already exist
	
	##================
	## extract significant loops
	##================
	if (opt$DistStrat == 0) {

		## two sets of significant loops
		## one for the conventional FDR, one for the IHW corrected FDR
		## loops significant in at least one input sample are reported
		
		##
		##========= set 1: conventional FDR
		##
		# edgeR_SigDir <- paste0(DiffLoopDir, '/DiffLoop_Sig_FDR_BH')
		edgeR_SigDir <- paste0(DiffLoopDir, '/DiffLoop_BH_FDR_', FDR_Th_DESeq, '_LOG2FC_', FOLD_Change_Thr)
		system(paste("mkdir -p", edgeR_SigDir))
		OutLoopFile <- paste0(edgeR_SigDir, '/DiffLoops_ALL.bed')		
		system(paste0("awk -F\'[\t]\' \'function abs(v) {return v < 0 ? -v : v} {if ((NR==1) || (($(NF-3) <= ", FDR_Th_DESeq, ") && (abs($(NF-6)) >= ", FOLD_Change_Thr, "))) {print $0}}\' ", CompleteLoopFile, " | awk -F\'[\t]\' \'((NR==1) || ($(NF-1)>0) || ($NF>0))\' - > ", OutLoopFile))
		## categorize these edgeR significant loops 
		## according to differential / non-differential ChIP coverage		
		Categorize_DiffLoops(edgeR_SigDir, OutLoopFile, CategoryList, QVal_ColList, opt$Model, HiC_Data_Process, opt$CovThr)	


		##
		##========= set 2: FDR computed by IHW
		##		
		# edgeR_SigDir <- paste0(DiffLoopDir, '/DiffLoop_Sig_FDR_IHW')
		edgeR_SigDir <- paste0(DiffLoopDir, '/DiffLoop_IHW_FDR_', FDR_Th_DESeq, '_LOG2FC_', FOLD_Change_Thr)
		system(paste("mkdir -p", edgeR_SigDir))
		OutLoopFile <- paste0(edgeR_SigDir, '/DiffLoops_ALL.bed')
		system(paste0("awk -F\'[\t]\' \'function abs(v) {return v < 0 ? -v : v} {if ((NR==1) || (($(NF-2) <= ", FDR_Th_DESeq, ") && (abs($(NF-6)) >= ", FOLD_Change_Thr, "))) {print $0}}\' ", CompleteLoopFile, " | awk -F\'[\t]\' \'((NR==1) || ($(NF-1)>0) || ($NF>0))\' - > ", OutLoopFile))
		## categorize these edgeR significant loops 
		## according to differential / non-differential ChIP coverage
		Categorize_DiffLoops(edgeR_SigDir, OutLoopFile, CategoryList, QVal_ColList, opt$Model, HiC_Data_Process, opt$CovThr)	
		
	} else {

		##
		##========= set 1: computed by using FDR_BH column (global FDR)
		##
		# edgeR_SigDir <- paste0(DiffLoopDir, '/DiffLoop_Sig_FDR_BH')
		edgeR_SigDir <- paste0(DiffLoopDir, '/DiffLoop_BH_FDR_', FDR_Th_DESeq, '_LOG2FC_', FOLD_Change_Thr)
		system(paste("mkdir -p", edgeR_SigDir))
		OutLoopFile <- paste0(edgeR_SigDir, '/DiffLoops_ALL.bed')		
		system(paste0("awk -F\'[\t]\' \'function abs(v) {return v < 0 ? -v : v} {if ((NR==1) || (($(NF-2) <= ", FDR_Th_DESeq, ") && (abs($(NF-6)) >= ", FOLD_Change_Thr, "))) {print $0}}\' ", CompleteLoopFile, " | awk -F\'[\t]\' \'((NR==1) || ($(NF-1)>0) || ($NF>0))\' - > ", OutLoopFile))
		## categorize these edgeR significant loops 
		## according to differential / non-differential ChIP coverage		
		Categorize_DiffLoops(edgeR_SigDir, OutLoopFile, CategoryList, QVal_ColList, opt$Model, HiC_Data_Process, opt$CovThr)
	}

	cat(sprintf("\n\n *** Applied EdgeR for differential loop finding *** \n\n"))

}	# end EdgeR condition

