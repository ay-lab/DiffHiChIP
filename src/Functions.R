#!/usr/bin/env Rscript

#===========================================================
## differential analysis on HiChIP and Hi-C loops
# Author: Sourya Bhattacharyya
# Vijay-Ay lab, LJI
#===========================================================

source("Header.R")

##==========
## number of lines in a file
##==========
GetNumLines <- function(inpfile) {
	nline <- as.integer(system(paste("cat", inpfile, "| wc -l"), intern = TRUE))
	return(nline)
}

##==========
## overlap of 1D segments
##==========
Overlap1D <- function(Inpdata1, Inpdata2, boundary=1, offset=0, uniqov=TRUE, chrWise=FALSE) {

	if (chrWise == FALSE) {
		ov1 <- as.data.frame(findOverlaps(GRanges(Inpdata1[,1], IRanges(Inpdata1[,2]+boundary-offset, Inpdata1[,3]-boundary+offset)),GRanges(Inpdata2[,1], IRanges(Inpdata2[,2]+boundary-offset, Inpdata2[,3]-boundary+offset))))
		if (uniqov == TRUE) {
			ov_idx_file1 <- unique(ov1[,1])
			ov_idx_file2 <- unique(ov1[,2])		
		} else {
			ov_idx_file1 <- ov1[,1]
			ov_idx_file2 <- ov1[,2]
		}
		nonov_idx_file1 <- setdiff(seq(1, nrow(Inpdata1)), ov_idx_file1)
		nonov_idx_file2 <- setdiff(seq(1, nrow(Inpdata2)), ov_idx_file2)

	} else {
		# get the chromosome list from the given peaks
		ChrList <- union(unique(Inpdata1[,1]), unique(Inpdata2[,1]))

		ov_idx_file1 <- c()
		ov_idx_file2 <- c()
		nonov_idx_file1 <- c()
		nonov_idx_file2 <- c()
		
		# process individual chromosomes in a loop
		for (idx in (1:length(ChrList))) {
			CurrChr <- ChrList[idx]
			# peak indices of first data corresponding to the current chromosome
			Inpdata1_CurrChr_IdxSet <- which(Inpdata1[,1] == CurrChr)
			# peak indices of second data corresponding to the current chromosome
			Inpdata2_CurrChr_IdxSet <- which(Inpdata2[,1] == CurrChr)
			# proceed if both data have some entries for the current chromosome
			if ((length(Inpdata1_CurrChr_IdxSet) == 0) | (length(Inpdata2_CurrChr_IdxSet) == 0)) {
	 			next
			}
			# extract the input datasets limited to the current chromosome
			Inpdata1_CurrChr <- Inpdata1[Inpdata1_CurrChr_IdxSet, ]
			Inpdata2_CurrChr <- Inpdata2[Inpdata2_CurrChr_IdxSet, ]
			
			# compute the overlap
			ov1 <- as.data.frame(findOverlaps(GRanges(Inpdata1_CurrChr[,1], IRanges(Inpdata1_CurrChr[,2]+boundary-offset, Inpdata1_CurrChr[,3]-boundary+offset)),GRanges(Inpdata2_CurrChr[,1], IRanges(Inpdata2_CurrChr[,2]+boundary-offset, Inpdata2_CurrChr[,3]-boundary+offset))))
			
			# overlapping / non-overlapping indices for both categories
			if (uniqov == TRUE) {
				temp_ov_idx_file1 <- unique(ov1[,1])
				temp_ov_idx_file2 <- unique(ov1[,2])		
			} else {
				temp_ov_idx_file1 <- ov1[,1]
				temp_ov_idx_file2 <- ov1[,2]
			}
			temp_nonov_idx_file1 <- setdiff(seq(1, nrow(Inpdata1_CurrChr)), temp_ov_idx_file1)
			temp_nonov_idx_file2 <- setdiff(seq(1, nrow(Inpdata2_CurrChr)), temp_ov_idx_file2)

			# indices with respect to original loop indices
			# (i.e. with respect to Inpdata1 and Inpdata2)
			if (length(temp_ov_idx_file1) > 0) {
				ov_idx_file1_currchr <- Inpdata1_CurrChr_IdxSet[temp_ov_idx_file1]	
				# append in the final list
				ov_idx_file1 <- c(ov_idx_file1, ov_idx_file1_currchr)
			}
			
			if (length(temp_ov_idx_file2) > 0) {
				ov_idx_file2_currchr <- Inpdata2_CurrChr_IdxSet[temp_ov_idx_file2]
				# append in the final list
				ov_idx_file2 <- c(ov_idx_file2, ov_idx_file2_currchr)
			}
			
			if (length(temp_nonov_idx_file1) > 0) {
				nonov_idx_file1_currchr <- Inpdata1_CurrChr_IdxSet[temp_nonov_idx_file1]
				# append in the final list
				nonov_idx_file1 <- c(nonov_idx_file1, nonov_idx_file1_currchr)
			}
			
			if (length(temp_nonov_idx_file2) > 0) {
				nonov_idx_file2_currchr <- Inpdata2_CurrChr_IdxSet[temp_nonov_idx_file2]
				# append in the final list
				nonov_idx_file2 <- c(nonov_idx_file2, nonov_idx_file2_currchr)
			}
		}	# end chromosomes loop
	}

	# return the overlapping and non-overlapping set of indices
	newList <- list(A_AND_B = ov_idx_file1, B_AND_A = ov_idx_file2, A_MINUS_B = nonov_idx_file1, B_MINUS_A = nonov_idx_file2)
	return(newList)

}

##==========
## extract chromosome wise 3D data
##==========
ExtractChrData <- function(InpFile, chrName, OutFile=NULL, header=TRUE, dist=c(-1,-1), mid=FALSE) {
	if (is.null(OutFile)) {
		OutFile <- paste0(dirname(InpFile), "/temp_Chr_data.bed")
	}

	# process the distance thresholds
	# and insert in two variables
	if ((dist[1] > 0) & (dist[2] > 0) & (dist[2] > dist[1])) {
		distthrlow <- dist[1]
		distthrhigh <- dist[2]		
	} else {
		distthrlow <- -1
		distthrhigh <- -1
	}

	# condition based on using gzipped input file
	# or plain text file

	if (file_ext(InpFile) == "gz") {
		if (header == TRUE) {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		} else {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))						
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		}
	} else {
		if (header == TRUE) {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		} else {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))						
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		}
	}

}	# end function

#=======================
# This function creates two new files for storing the scaled ChIP coverage values
# from the merged ChIP coverage data 
# we use the mean of the ChIP-seq coverages
## also assigns labels to the (ND, LD, HD) 
## and assigns labels to the Merged_ChIPCovData
#=======================
Create_Scaled_ChIP_Coverage_with_Label <- function(MergedChIPCovFile, ChIPAlignFileCountVec, scaled_ChIPCovFile1, scaled_ChIPCovFile2, ChIP_1D_NonDiff_BinFile) {

	# merged ChIP-seq coverage for all input samples
	Merged_ChIPCovData <- data.table::fread(MergedChIPCovFile, header=T)	

	# read non-differential 1D bins
	NonDiff_ChIP_BinData <- data.table::fread(ChIP_1D_NonDiff_BinFile, header=T)

	##==============
	## scaling of ChIP-seq coverages 
	##==============

	# ChIP coverage of the first category - row means operation
	if (ChIPAlignFileCountVec[1] > 1) {
		ChIPCovCat1 <- rowMeans(Merged_ChIPCovData[, 4:(4+ChIPAlignFileCountVec[1]-1)])
	} else {
		ChIPCovCat1 <- Merged_ChIPCovData[, 4]
	}

	# ChIP coverage of the second category - row means operation
	if (ChIPAlignFileCountVec[2] > 1) {
		ChIPCovCat2 <- rowMeans(Merged_ChIPCovData[, (4+ChIPAlignFileCountVec[1]):(4+ChIPAlignFileCountVec[1]+ChIPAlignFileCountVec[2]-1)])
	} else {
		ChIPCovCat2 <- Merged_ChIPCovData[, (4+ChIPAlignFileCountVec[1])]
	}

	# scaling factor to be applied on the coverage values of first category
	scaling_factor <- (sum(ChIPCovCat2) * 1.0) / sum(ChIPCovCat1)
	cat(sprintf("\n ===>> Scaling ChIP coverage - scaling_factor : %s ====>> \n", scaling_factor))

	# scale the ChIP-seq coverage of the first category
	for (i in (1:ChIPAlignFileCountVec[1])) { 
		Merged_ChIPCovData[, (4+i-1)] <- as.integer(round(Merged_ChIPCovData[, (4+i-1)] * scaling_factor))
	}	

	##==============
	# mean scaled ChIP coverage of the first category
	##==============	
	if (ChIPAlignFileCountVec[1] > 1) {
		ChIPCovCat1 <- rowMeans(Merged_ChIPCovData[, 4:(4+ChIPAlignFileCountVec[1]-1)])
	} else {
		ChIPCovCat1 <- Merged_ChIPCovData[, 4]
	}
	# mean scaled ChIP coverage of the second category
	if (ChIPAlignFileCountVec[2] > 1) {
		ChIPCovCat2 <- rowMeans(Merged_ChIPCovData[, (4+ChIPAlignFileCountVec[1]):(4+ChIPAlignFileCountVec[1]+ChIPAlignFileCountVec[2]-1)])
	} else {
		ChIPCovCat2 <- Merged_ChIPCovData[, (4+ChIPAlignFileCountVec[1])]
	}

	cat(sprintf("\n\n *** Performed scaling of ChIP coverage values for uniform ChIP coverage *** \n\n"))

	##==============
	# label the ChIP bins as HD, LD or ND
	# depending on their EdgeR output, and ChIP coverage comparison
	##==============

	# 1D bins with very low deviation in their ChIP-seq coverage (ND)
	MaxCovVec <- pmax(ChIPCovCat1, ChIPCovCat2)
	MinCovVec <- pmin(ChIPCovCat1, ChIPCovCat2)
	BinIdxLowCovDiff <- which(MinCovVec > ((1 - ChIP_Cov_Thr) * MaxCovVec))

	# label the ChIP bins as HD, LD or ND
	# depending on their EdgeR output, and ChIP coverage comparison
	ov <- Overlap1D(Merged_ChIPCovData[, 1:3], NonDiff_ChIP_BinData[, 1:3])

	## ND: non-differential + very low coverage deviation
	ND_Label_Idx <- intersect(ov$A_AND_B, BinIdxLowCovDiff)

	## LD: non-differential + non-ND
	LD_Label_Idx <- setdiff(ov$A_AND_B, ND_Label_Idx)

	## HD: differential
	HD_Label_Idx <- setdiff(setdiff(seq(1, nrow(Merged_ChIPCovData)), LD_Label_Idx), ND_Label_Idx)

	cat(sprintf("\n Number of 1D bins: %s \n number of 1D bins classified as HD: %s \n number of 1D bins classified as LD: %s \n number of 1D bins classified as ND: %s ", nrow(Merged_ChIPCovData), length(HD_Label_Idx), length(LD_Label_Idx), length(ND_Label_Idx)))

	# assign the 1D label information
	ChIPLabelVec <- rep('HD', nrow(Merged_ChIPCovData))
	ChIPLabelVec[ND_Label_Idx] <- 'ND'
	ChIPLabelVec[LD_Label_Idx] <- 'LD'

	# append the 1D label in the merged ChIP coverage data
	# and write back to the merged ChIP coverage file
	CN <- colnames(Merged_ChIPCovData)
	Merged_ChIPCovData <- cbind.data.frame(Merged_ChIPCovData, ChIPLabelVec)
	colnames(Merged_ChIPCovData) <- c(CN, 'Label')
	write.table(Merged_ChIPCovData, MergedChIPCovFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)

	# create two data frames using the scaled mean of ChIP coverage values
	covdf1 <- cbind.data.frame(Merged_ChIPCovData[, 1:3], ChIPCovCat1, Merged_ChIPCovData[, ncol(Merged_ChIPCovData)])
	colnames(covdf1) <- c(colnames(Merged_ChIPCovData)[1:3], paste0(CategoryList[1], '_meanChIPCov'), 'Label')

	covdf2 <- cbind.data.frame(Merged_ChIPCovData[, 1:3], ChIPCovCat2, Merged_ChIPCovData[, ncol(Merged_ChIPCovData)])
	colnames(covdf2) <- c(colnames(Merged_ChIPCovData)[1:3], paste0(CategoryList[2], '_meanChIPCov'), 'Label')

	write.table(covdf1, scaled_ChIPCovFile1, row.names=F, col.names=T, sep="\t", quote=F, append=F)
	write.table(covdf2, scaled_ChIPCovFile2, row.names=F, col.names=T, sep="\t", quote=F, append=F)

	## remve temporary objects
	if (exists("ChIPCovCat1")) {
		rm("ChIPCovCat1")
	}
	if (exists("ChIPCovCat2")) {
		rm("ChIPCovCat2")
	}
	if (exists("covdf1")) {
		rm("covdf1")
	}
	if (exists("covdf2")) {
		rm("covdf2")
	}
	if (exists("Merged_ChIPCovData")) {
		rm("Merged_ChIPCovData")
	}
	if (exists("ov")) {
		rm(ov)
	}
	if (exists("NonDiff_ChIP_BinData")) {
		rm(NonDiff_ChIP_BinData)
	}
	gc()

}	# end function

#======================
# this function creates merged ChIP-seq coverage file, from the input ChIP-seq files of both categories
#======================
Create_Merged_ChIP_Coverage_File <- function(MergedChIPCovFile, OutDir, 
											ChrSizeFile, BinSize, ChIPAlignFileList, ChIPAlignFileCountVec, CategoryList) {

	tempfile1 <- paste0(OutDir, '/ChIP_Coverage_temp1.bed')
	tempfile2 <- paste0(OutDir, '/ChIP_Coverage_temp2.bed')

	# first get the binned distribution of specified chromosome size file
	TargetBinnedChrFile <- paste0(OutDir, '/Binned_Chromosome_Size_RefGenome.bed')
	system(paste("bedtools makewindows -g", ChrSizeFile, "-w", BinSize, "| sort -k1,1 -k2,2n >", TargetBinnedChrFile))

	# process individual alignment files and first get the ChIP seq coverage for those files
	# using the TargetBinnedChrFile
	for (i in (1:length(ChIPAlignFileList))) {
		if (tools::file_ext(ChIPAlignFileList[i]) == "bam") {
			cat(sprintf("\n ===>> Merging ChIP coverage - processing ChIP-seq alignment file in BAM format : %s ====>> \n", ChIPAlignFileList[i]))
		} else if (tools::file_ext(ChIPAlignFileList[i]) == "gz") {
			cat(sprintf("\n ===>> Merging ChIP coverage - processing ChIP-seq alignment file in gzipped bedgraph (4 column) format : %s ====>> \n", ChIPAlignFileList[i]))
		} else {
			cat(sprintf("\n ===>> Merging ChIP coverage - processing ChIP-seq alignment file in bedgraph (4 column) format : %s ====>> \n", ChIPAlignFileList[i]))
		}
		if (tools::file_ext(ChIPAlignFileList[i]) == "gz") {
			# gzipped bedgraph format
			# bedgraph files are assumed to be sorted
			if (i == 1) {		
				system(paste("zcat", ChIPAlignFileList[i], "| bedtools coverage -sorted -a", TargetBinnedChrFile, "-b stdin -counts >", MergedChIPCovFile))
			} else {
				system(paste("zcat", ChIPAlignFileList[i], "| bedtools coverage -sorted -a", TargetBinnedChrFile, "-b stdin -counts | cut -f4 >", tempfile1))
				system(paste("paste", MergedChIPCovFile, tempfile1, ">", tempfile2))
				system(paste("mv", tempfile2, MergedChIPCovFile))			
			}
		} else {
			# either BAM file or plain bedgraph format
			# bedgraph files are assumed to be sorted
			if (i == 1) {
				if (tools::file_ext(ChIPAlignFileList[i]) == "bam") {	
					system(paste("bedtools coverage -a", TargetBinnedChrFile, "-b", ChIPAlignFileList[i], "-counts >", MergedChIPCovFile))
				} else {
					system(paste("bedtools coverage -sorted -a", TargetBinnedChrFile, "-b", ChIPAlignFileList[i], "-counts >", MergedChIPCovFile))
				}
			} else {	
				if (tools::file_ext(ChIPAlignFileList[i]) == "bam") {	
					system(paste("bedtools coverage -a", TargetBinnedChrFile, "-b", ChIPAlignFileList[i], "-counts | cut -f4 >", tempfile1))
				} else {					
					system(paste("bedtools coverage -sorted -a", TargetBinnedChrFile, "-b", ChIPAlignFileList[i], "-counts | cut -f4 >", tempfile1))
				}
				system(paste("paste", MergedChIPCovFile, tempfile1, ">", tempfile2))
				system(paste("mv", tempfile2, MergedChIPCovFile))
			}		
		}
	}	# end input file loop

	# now scale the ChIP coverages according to different categories
	Merged_ChIPCovData <- data.table::fread(MergedChIPCovFile, header=F)

	# now write the scaled coverage values
	colnames(Merged_ChIPCovData) <- c('chr', 'start', 'end', paste0(CategoryList[1], '_ChIPCov_', seq(1,ChIPAlignFileCountVec[1])), paste0(CategoryList[2], '_ChIPCov_', seq(1,ChIPAlignFileCountVec[2])))
	write.table(Merged_ChIPCovData, MergedChIPCovFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)

	# now remove the temporary files
	if (file.exists(tempfile1)) {
		system(paste("rm", tempfile1))
	}
	if (file.exists(tempfile2)) {
		system(paste("rm", tempfile2))
	}	
	if (file.exists(TargetBinnedChrFile)) {
		system(paste("rm", TargetBinnedChrFile))
	}

	## remove temporary objects
	if (exists("Merged_ChIPCovData")) {
		rm("Merged_ChIPCovData")
	}
	if (exists("ChIPCovCat1")) {
		rm("ChIPCovCat1")
	}
	if (exists("ChIPCovCat2")) {
		rm("ChIPCovCat2")	
	}
	gc()	

}	# end function

#======================================================
# applies EdgeR to the input count matrix for 1D ChIP coverage data
# returns a set of differential bins 
# parameters:
# MainDir: main base directory for EdgeR outputs
# CountData: only the count matrix employed for EdgeR design
# CategoryList: two input categories
# ReplicaCount: list of two elements, corresponding to the count of replicates in each categories
# prefixstr: if 'Loops', computes differential loops (3D), else computes differential (Hi)ChIP (1D)
# bcv: default threshold for single replicates data
#======================================================
ApplyEdgeR_ChIP <- function(MainDir, CountData, CategoryList, ReplicaCount, prefixstr, bcv=0.4) {

	# directory to store the plots related to EdgeR
	PlotDir <- paste0(MainDir, '/Plots_EdgeR_ChIP')
	system(paste("mkdir -p", PlotDir))

	# category and replicate information for EdgeR structure
	GroupDistrVec <- c(rep(CategoryList[1], ReplicaCount[1]), rep(CategoryList[2], ReplicaCount[2]))
	cat(sprintf("\n ===>>> GroupDistrVec : %s ", paste0(GroupDistrVec, sep="\t")))
	
	# create the EdgeR count data structure (7th column onwards)
	y <- DGEList(counts=CountData, group=GroupDistrVec)

	# generate the design (according to the group distribution)
	design <- model.matrix(~GroupDistrVec) 

	if ((ReplicaCount[1] > 1) & (ReplicaCount[2] > 1)) {
		# for data with multiple replicates, better to normalize the DGEList object
		y <- calcNormFactors(y)			
	}

	# Estimate Common, Trended and Tagwise Negative Binomial dispersions 
	# by weighted likelihood empirical Bayes
	y <- estimateDisp(y, design)

	if ((ReplicaCount[1] > 1) & (ReplicaCount[2] > 1)) {
		# approach 1 - identifying DE tags - qCML method
		# exact test (applicable for experiments with single factor) 
		# to generate matrix of pseudo counts
		et <- exactTest(y, dispersion="trended", pair=c(CategoryList[1], CategoryList[2]))	
		# get the DE genes - useful for plotting MA 
		et_tags <- topTags(et)			

	} else {
		# default exact test when no replicates are available
		et <- exactTest(y, dispersion=bcv^2)
		# get the DE genes - useful for plotting MA 
		et_tags <- topTags(et)			
	}

	# # get log-counts-per-million (cpm)
	# logcpm <- cpm(y, prior.count=2, log=TRUE)

	#===============
	if ((ReplicaCount[1] > 1) & (ReplicaCount[2] > 1)) {
	
		# statistics
		# using plotMDS function to generate a plot in which distances between samples 
		# correspond to leading biological coefficient of variation (BCV) between those samples:
		MDSPlotFile <- paste0(PlotDir, '/EdgeR_MDS_Plot.pdf')
		pdf(MDSPlotFile, width=8, height=6)
		plotMDS(y)
		dev.off()

		# statistics
		# dispersion estimates can be viewed in a BCV plot
		BCVPlotFile <- paste0(PlotDir, '/EdgeR_BCV_Plot.pdf')
		pdf(BCVPlotFile, width=8, height=6)
		plotBCV(y)
		dev.off()

		# Plot log-fold change against log-counts per million, 
		# with DE genes highlighted
		# using qlf (output of glmQLFTest)
		et_MDPLotFile <- paste0(PlotDir, '/EdgeR_et_MD_Plot.pdf')
		pdf(et_MDPLotFile, width=8, height=6)
		plotMD(et, hl.cex = 0.1)
		abline(h=c(-1, 1), col="blue")
		dev.off()		

		# smear plot - MA for groups of samples
		# approach 1 - exact test
		SmearPlotFile <- paste0(PlotDir, '/EdgeR_Smear_Plot_et.pdf')
		pdf(SmearPlotFile, width=8, height=6)
		plotSmear(et, pair=c(CategoryList[1], CategoryList[2]), de.tags=et_tags)
		dev.off()

	}	# end replica count condition

	#===============
	# for each of the DE tag finding techniques (qCML or CR methods)
	# convert the p-values to the corresponding q-values, using BH correction
	#===============
	Qval_et <- p.adjust(et$table$PValue, method = "BH")		# approach 1 - exact test

	# just print the EdgeR significance results,
	# without main loop / bin information
	EdgeRRes_et <- cbind.data.frame(et$table$logFC, et$table$logCPM, et$table$PValue, Qval_et)
	colnames(EdgeRRes_et) <- c('logFC', 'logCPM', 'PValue', 'FDR')

	ResOutFile <- paste0(MainDir, '/', prefixstr, '_EdgeR_Res_temp.bed')
	write.table(EdgeRRes_et, ResOutFile, row.names=F, col.names=T, quote=F, sep="\t", append=F)	

	## remove temporary objects
	if (exists("y")) {
		rm("y")
	}
	if (exists("design")) {
		rm("design")
	}
	if (exists("et")) {
		rm("et")
	}
	if (exists("et_tags")) {
		rm("et_tags")
	}
	gc()

}	# end ApplyEdgeR_ChIP function


#======================================================
# this function merges input FitHiChIP loops (basically clubs first 6 fields)
# parameters:
# LoopList: list of all the input files 
# MidList: vector of 0/1 specifying whether the midpoints of bin intervals are specified
# QColList: q-value column list of all input files
# BinSize: bin size of loops
# UnionLoopFile: where output is written
# ChrNames: list of chromosome names
# FDRFilt: if TRUE, last field is checked and only the loops with FDR value < FDRThr are considered for merging
# FDRThr: default 0.1 unless provided
#======================================================
MergeLoops <- function(LoopList, MidList, QColList, BinSize, 
						UnionLoopFile, ChrNames, FDRFilt=FALSE, FDRThr=0.1) {
	
	cat(sprintf("\n\n **** within function - MergeLoops *** \n\n"))
	
	## parallel processing for individual chromosomes
	lapply( 1:length(ChrNames), function(chrIdx) {	
	# for (chrIdx in (1:length(ChrNames))) {

		currChr <- ChrNames[chrIdx]
		cat(sprintf("\n processing chromosome : %s ", currChr))

		UnionLoopTempFile <- paste0(dirname(UnionLoopFile), '/temp_merge_union_', currChr, '.bed')
		UnionLoopTempFile_1 <- paste0(dirname(UnionLoopFile), '/temp_merge_union1_', currChr, '.bed')

		if (file.exists(UnionLoopTempFile) == TRUE) {
			system(paste("rm", UnionLoopTempFile))
		}
		if (file.exists(UnionLoopTempFile_1) == TRUE) {
			system(paste("rm", UnionLoopTempFile_1))
		}		

		for (i in (1:length(LoopList))) {
			# boolean value depicting the midpoint of the bin for the current loop file
			midval <- as.integer(MidList[i])
			# column storing the q-value for the current loop file			
			qcol <- as.integer(QColList[i])

			if (tools::file_ext(LoopList[i]) == "gz") {
				# gzipped input file
				if (FDRFilt == FALSE) {
					if (midval == 1) {
						system(paste0("zcat ", LoopList[i], " | awk -v b=", BinSize, " \'{if ((NR>1) && ($1==\"", currChr, "\")) {print $1\"\t\"($2-(b/2))\"\t\"($2+(b/2))\"\t\"$3\"\t\"($4-(b/2))\"\t\"($4+(b/2))}}\' - >> ", UnionLoopTempFile_1))
					} else {
						system(paste0("zcat ", LoopList[i], " | awk \'{if ((NR>1) && ($1==\"", currChr, "\")) {print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}}\' - >> ", UnionLoopTempFile_1))
					}
				} else {
					if (midval == 1) {
						if (qcol == -1) {
							system(paste0("zcat ", LoopList[i], " | awk -v b=", BinSize, " \'{if ((NR>1) && ($1==\"", currChr, "\") && ($NF != \"NA\") && (sprintf(\"%0.400f\",$NF) <", FDRThr, ")) {print $1\"\t\"($2-(b/2))\"\t\"($2+(b/2))\"\t\"$3\"\t\"($4-(b/2))\"\t\"($4+(b/2))}}\' - >> ", UnionLoopTempFile_1))
						} else {
							system(paste0("zcat ", LoopList[i], " | awk -v b=", BinSize, " \'{if ((NR>1) && ($1==\"", currChr, "\") && ($", qcol, " != \"NA\") && (sprintf(\"%0.400f\",$", qcol, ") <", FDRThr, ")) {print $1\"\t\"($2-(b/2))\"\t\"($2+(b/2))\"\t\"$3\"\t\"($4-(b/2))\"\t\"($4+(b/2))}}\' - >> ", UnionLoopTempFile_1))
						}
					} else {
						if (qcol == -1) {
							system(paste0("zcat ", LoopList[i], " | awk \'{if ((NR>1) && ($1==\"", currChr, "\") && ($NF != \"NA\") && (sprintf(\"%0.400f\",$NF) <", FDRThr, ")) {print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}}\' - >> ", UnionLoopTempFile_1))
						} else {
							system(paste0("zcat ", LoopList[i], " | awk \'{if ((NR>1) && ($1==\"", currChr, "\") && ($", qcol, " != \"NA\") && (sprintf(\"%0.400f\",$", qcol, ") <", FDRThr, ")) {print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}}\' - >> ", UnionLoopTempFile_1))
						}
					}
				}
			} else {
				# non-gzipped input file
				if (FDRFilt == FALSE) {
					if (midval == 1) {
						system(paste0("awk -v b=", BinSize, " \'{if ((NR>1) && ($1==\"", currChr, "\")) {print $1\"\t\"($2-(b/2))\"\t\"($2+(b/2))\"\t\"$3\"\t\"($4-(b/2))\"\t\"($4+(b/2))}}\' ", LoopList[i], " >> ", UnionLoopTempFile_1))
					} else {
						system(paste0("awk \'{if ((NR>1) && ($1==\"", currChr, "\")) {print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}}\' ", LoopList[i], " >> ", UnionLoopTempFile_1))
					}
				} else {
					if (midval == 1) {
						if (qcol == -1) {
							system(paste0("awk -v b=", BinSize, " \'{if ((NR>1) && ($1==\"", currChr, "\") && ($NF != \"NA\") && (sprintf(\"%0.400f\",$NF) <", FDRThr, ")) {print $1\"\t\"($2-(b/2))\"\t\"($2+(b/2))\"\t\"$3\"\t\"($4-(b/2))\"\t\"($4+(b/2))}}\' ", LoopList[i], " >> ", UnionLoopTempFile_1))
						} else {
							system(paste0("awk -v b=", BinSize, " \'{if ((NR>1) && ($1==\"", currChr, "\") && ($", qcol, " != \"NA\") && (sprintf(\"%0.400f\",$", qcol, ") <", FDRThr, ")) {print $1\"\t\"($2-(b/2))\"\t\"($2+(b/2))\"\t\"$3\"\t\"($4-(b/2))\"\t\"($4+(b/2))}}\' ", LoopList[i], " >> ", UnionLoopTempFile_1))
						}
					} else {
						if (qcol == -1) {
							system(paste0("awk \'{if ((NR>1) && ($1==\"", currChr, "\") && ($NF != \"NA\") && (sprintf(\"%0.400f\",$NF) <", FDRThr, ")) {print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}}\' ", LoopList[i], " >> ", UnionLoopTempFile_1))
						} else {
							system(paste0("awk \'{if ((NR>1) && ($1==\"", currChr, "\") && ($", qcol, " != \"NA\") && (sprintf(\"%0.400f\",$", qcol, ") <", FDRThr, ")) {print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}}\' ", LoopList[i], " >> ", UnionLoopTempFile_1))
						}
					}
				}
			}
		}	# end loop input files

		## combine all the loops
		system(paste("sort -k2,2n -k5,5n", UnionLoopTempFile_1, " | uniq > ", UnionLoopTempFile))

		## remove the temporary files
		if (file.exists(UnionLoopTempFile_1)) {
			system(paste("rm", UnionLoopTempFile_1))
		}

	})	# end processing chromosomes - parallel processing

	## create the final output file
	system(paste0("cat ", dirname(UnionLoopFile), "/temp_merge_union_*.bed > ", UnionLoopFile))

	## remove the temporary file
	system(paste0("rm ", dirname(UnionLoopFile), "/temp_merge_union_*.bed"))

	cat(sprintf("\n\n **** exit from function - MergeLoops *** \n\n"))

}	# end function

#================================
# this function annotates individual loops in the merged union set of loops
# according to the feature vectors of individual loops 
# parameters:
# UnionLoopFile: file containing merged set of loops from all input replicates and all categories
# UnionLoopFeatureFile: file which will contain all features corresponding to all loops in UnionLoopFile
# AllLoopList: Loops with FitHiChIP significance for different input replicates and categories
# CHRLIST_NAMENUM: list of chromosome names
# BinSize: bin size of loops
# ChIPCovFileList: ChIP-seq coverage file list (one per input category) for the bin size of loops
# AllRepLabels: labels of individual replicates of individual categories
# CategoryList: two categories experimented
# CCColList: contact count column list for all input files
# QColList: q-value column list for all input files
# MidList: sequence of 1s or 0s where 1 means midpoint of a bin is mentioned
#================================
FillFeatureValues <- function(UnionLoopFile, UnionLoopFeatureFile, 
								AllLoopList, BinSize, 
								ChIPCovFileList, AllRepLabels, 
								CategoryList, CCColList, QColList, MidList) {

	## get the chromosome list containing loops
	tempchrfile <- paste0(dirname(UnionLoopFile), '/temp_chrlist.bed')
	system(paste("cut -f1", UnionLoopFile, "|sort -k1,1 | uniq > ", tempchrfile))
	CHRLIST_DF <- data.table::fread(tempchrfile, header=F)
	CHRLIST <- as.vector(CHRLIST_DF[,1])

	lapply( 1:length(CHRLIST), function(chr_idx) {	
	# for (chr_idx in (1:length(CHRLIST_NAMENUM))) {
		chrName <- CHRLIST[chr_idx]
		# cat(sprintf("\n ===>>> Within function FillFeatureValues --- processing chromosome : %s ", chrName))

		## union loops for the current chromosome
		UnionLoopTempFile1 <- paste0(dirname(UnionLoopFile), '/temp_Merged_Loops_', chrName, '.bed')
		## temporary loop + feature containing file used per iteration
		temp_final_UnionLoopFile <- paste0(dirname(UnionLoopFile), '/temp_mastersheet_', chrName, '.bed')
		## stores temporary FitHiChIP loops for the current chromosome
		InpTempFitHiChIPLoopFile <- paste0(dirname(UnionLoopFile), '/temp_InpFitHiChIPLoopFile_', chrName, '.bed')
		## stores temporary ChIP coverage for the current chromosome
		InpTempChIPCoverageFile <- paste0(dirname(UnionLoopFile), '/temp_InpChIPCoverageFile_', chrName, '.bed')

		# first extract the loops involving current chromosome		
		ExtractChrData(UnionLoopFile, chrName, UnionLoopTempFile1, header=FALSE)
		MergedIntTempData <- data.table::fread(UnionLoopTempFile1, header=F)

		# also get the interacting bins (start position divided by the bin size)
		AllLoop_BinDF <- cbind.data.frame((MergedIntTempData[,2] / BinSize), (MergedIntTempData[,5] / BinSize))
		colnames(AllLoop_BinDF) <- c('B1', 'B2')

		#
		# ========== process loop files of all input samples for the current chromosome ============
		# 
		
		# vectors to store feature values of individual loops
		RawCC_Categ <- list()
		QVal_Categ <- list()

		for (i in (1:length(AllLoopList))) {
			# default initialization for the current chromosome and current input FitHiChIP loop file
			rawccvec <- rep(0, nrow(MergedIntTempData))
			qvec <- rep(1, nrow(MergedIntTempData))

			# input file containing FitHiChIP significance for all loops
			inpfile <- AllLoopList[i]
			cccol <- as.integer(CCColList[i])
			qcol <- as.integer(QColList[i])

			# extract the interacting bins
			# and the raw contact count and the q-values
			if (tools::file_ext(inpfile) == "gz") {
				if (QColList[i] == -1) {
					if (MidList[i] == 1) {
						system(paste0("zcat ", inpfile, " | awk -v b=", BinSize, " \'{if ((NR>1) && ($1==\"", chrName, "\")) {print (($2-(b/2))/b)\"\t\"(($4-(b/2))/b)\"\t\"$", cccol, "\"\t\"$NF}}\' - > ", InpTempFitHiChIPLoopFile))
					} else {
						system(paste0("zcat ", inpfile, " | awk -v b=", BinSize, " \'{if ((NR>1) && ($1==\"", chrName, "\")) {print ($2/b)\"\t\"($5/b)\"\t\"$", cccol, "\"\t\"$NF}}\' - > ", InpTempFitHiChIPLoopFile))
					}
				} else {
					if (MidList[i] == 1) {
						system(paste0("zcat ", inpfile, " | awk -v b=", BinSize, " \'{if ((NR>1) && ($1==\"", chrName, "\")) {print (($2-(b/2))/b)\"\t\"(($4-(b/2))/b)\"\t\"$", cccol, "\"\t\"$", qcol, "}}\' - > ", InpTempFitHiChIPLoopFile))
					} else {
						system(paste0("zcat ", inpfile, " | awk -v b=", BinSize, " \'{if ((NR>1) && ($1==\"", chrName, "\")) {print ($2/b)\"\t\"($5/b)\"\t\"$", cccol, "\"\t\"$", qcol, "}}\' - > ", InpTempFitHiChIPLoopFile))
					}
				}				
			} else {
				if (QColList[i] == -1) {
					if (MidList[i] == 1) {
						system(paste0("awk -v b=", BinSize, " \'{if ((NR>1) && ($1==\"", chrName, "\")) {print (($2-(b/2))/b)\"\t\"(($4-(b/2))/b)\"\t\"$", cccol, "\"\t\"$NF}}\' ", inpfile, " > ", InpTempFitHiChIPLoopFile))
					} else {
						system(paste0("awk -v b=", BinSize, " \'{if ((NR>1) && ($1==\"", chrName, "\")) {print ($2/b)\"\t\"($5/b)\"\t\"$", cccol, "\"\t\"$NF}}\' ", inpfile, " > ", InpTempFitHiChIPLoopFile))
					}
				} else {
					if (MidList[i] == 1) {
						system(paste0("awk -v b=", BinSize, " \'{if ((NR>1) && ($1==\"", chrName, "\")) {print (($2-(b/2))/b)\"\t\"(($4-(b/2))/b)\"\t\"$", cccol, "\"\t\"$", qcol, "}}\' ", inpfile, " > ", InpTempFitHiChIPLoopFile))
					} else {
						system(paste0("awk -v b=", BinSize, " \'{if ((NR>1) && ($1==\"", chrName, "\")) {print ($2/b)\"\t\"($5/b)\"\t\"$", cccol, "\"\t\"$", qcol, "}}\' ", inpfile, " > ", InpTempFitHiChIPLoopFile))
					}
				}
			}			

			# number of loops for the current chromosome and for the current sample
			nreadInp <- GetNumLines(InpTempFitHiChIPLoopFile)
			if (nreadInp > 0) {
				InpTempData <- data.table::fread(InpTempFitHiChIPLoopFile, header=F)
				colnames(InpTempData) <- c('B1', 'B2', 'RawCC', 'qval')
				cat(sprintf("\n ***** Computing overlap of merged loops with the FitHiChIP loop file index: %s  name: %s for the chromosome : %s ***** \n", i, inpfile, chrName))

				## merge the data frames - keeping all rows in "AllLoop_BinDF" intact
				mergeDF <- dplyr::left_join(AllLoop_BinDF, InpTempData)

				## store the feature values, and replace the NA entries
				rawccvec <- mergeDF$RawCC
				qvec <- mergeDF$qval
				rawccvec[is.na(rawccvec)] <- 0
				qvec[is.na(qvec)] <- 1
				
			}	# end number of reads condition

			# assign the feature vectors for the current input file
			# into the final list of feature arrays
			RawCC_Categ[[i]] <- rawccvec
			QVal_Categ[[i]] <- qvec

			cat(sprintf("\n ***** Assigned contact count and q-values for the FitHiChIP loop file index: %s for the chromosome : %s ***** \n", i, chrName))

		}	# end loop FitHiChIP significance files

		#
		# ============= process input ChIP-seq coverage files for two input categories ==========
		#
		if (length(ChIPCovFileList) == 2) {

  			# feature vectors for ChIP information (loop wise)
	  		Seg1_ChIP_Coverage <- list()
			Seg2_ChIP_Coverage <- list()
			Seg1_ChIP_Label <- rep('HD', nrow(MergedIntTempData))
			Seg2_ChIP_Label <- rep('HD', nrow(MergedIntTempData))
		
			for (i in (1:length(ChIPCovFileList))) {
				seg1_chip_coverage_vec <- rep(0, nrow(MergedIntTempData))
				seg2_chip_coverage_vec <- rep(0, nrow(MergedIntTempData))
				currcovfile <- ChIPCovFileList[i]
				cat(sprintf("\n Merging reference ChIP-seq coverage values: file number : %s  file name : %s ", i, currcovfile))

				# extract ChIP-seq coverage for the current chromosome
				# field 1: chromosome name, fields 2 and 3: bin 
				# field 4: coverage, field 5: bin label - LD, HD or LD
				# here we extract bin index, coverage and bin label
				system(paste0("awk -v b=", BinSize, " \'{if (($1==\"", chrName, "\") && ($2 ~ /^[0-9]+$/) && ($3 ~ /^[0-9]+$/)) {print ($2/b)\"\t\"$4\"\t\"$5}}\' ", currcovfile, " > ", InpTempChIPCoverageFile))

				# read the ChIP coverage for the current sample and for the current chromosome
				ChIPCoverageData <- data.table::fread(InpTempChIPCoverageFile, header=F)

				# *** Note: the below mentioned MATCH operation is only possible 
				# since we perform exact overlap of the segments
				
				# get the overlap of the first interacting bin (B1) of AllLoop_BinDF
				# with the bins in ChIPCoverageData
				m <- match(AllLoop_BinDF[,1], ChIPCoverageData[,1])
				idx_Loop <- which(!is.na(m))
				idx_Cov <- m[idx_Loop]
				# copy in the segment 1 ChIP coverage vector
				seg1_chip_coverage_vec[idx_Loop] <- ChIPCoverageData[idx_Cov, 2]
				if (i == 1) {
					Seg1_ChIP_Label[idx_Loop] <- ChIPCoverageData[idx_Cov, ncol(ChIPCoverageData)]
				}

				# get the overlap of the second interacting bin (B2) of AllLoop_BinDF
				# with the bins in ChIPCoverageData
				m <- match(AllLoop_BinDF[,2], ChIPCoverageData[,1])
				idx_Loop <- which(!is.na(m))
				idx_Cov <- m[idx_Loop]
				# copy in the segment 2 ChIP coverage vector
				seg2_chip_coverage_vec[idx_Loop] <- ChIPCoverageData[idx_Cov, 2]
				if (i == 1) {
					Seg2_ChIP_Label[idx_Loop] <- ChIPCoverageData[idx_Cov, ncol(ChIPCoverageData)]
				}

				# assign the feature vectors for the current input file
				# into the final list of feature arrays
				Seg1_ChIP_Coverage[[i]] <- seg1_chip_coverage_vec
				Seg2_ChIP_Coverage[[i]] <- seg2_chip_coverage_vec

				cat(sprintf("\n ***** Assigned reference ChIP-seq coverage information for the file index: %s for the chromosome : %s ***** \n", i, chrName))

			}	# end ChIP coverage file loop

		}	# end ChIP coverage file presence condition

		#
		# ============= for the current chromosome, merge the interacting regions =============
		# ============= along with the features accumulated	=============
		#
		namesvec <- c("chr1", "start1", "end1", "chr2", "start2", "end2")

		# if ChIP-seq coverage information is computed, insert coverage and label (differential)
		if (length(ChIPCovFileList) == 2) {
			# ChIP coverage 
			for (i in (1:2)) {
				MergedIntTempData <- cbind.data.frame(MergedIntTempData, Seg1_ChIP_Coverage[[i]], Seg2_ChIP_Coverage[[i]])
			}
			# ChIP label
			MergedIntTempData <- cbind.data.frame(MergedIntTempData, Seg1_ChIP_Label, Seg2_ChIP_Label)
			# update the column labels as well
			namesvec <- c(namesvec, paste0(CategoryList[1], c('_ChIPCov1', '_ChIPCov2')), paste0(CategoryList[2], c('_ChIPCov1', '_ChIPCov2')), 'Bin1_Label', 'Bin2_Label')
		}	# end ChIP coverage file presence condition

		# insert for individual samples (replicates within categories), different feature values
		for (i in (1:length(AllLoopList))) { 
			MergedIntTempData <- cbind.data.frame(MergedIntTempData, RawCC_Categ[[i]], QVal_Categ[[i]])
		}
		# update column names
		appendnamevec <- c()
		for (i in (1:length(AllRepLabels))) {
			for (v in c('_RawCC', '_QVal')) {
				appendnamevec <- c(appendnamevec, paste0(AllRepLabels[i], v))
			}
		}
		namesvec <- c(namesvec, appendnamevec)
		colnames(MergedIntTempData) <- namesvec

		write.table(MergedIntTempData, temp_final_UnionLoopFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)

		# delete the temporary file
		if (file.exists(UnionLoopTempFile1) == TRUE) {
			system(paste("rm", UnionLoopTempFile1))
		}

		# delete the temporary file
		if (file.exists(InpTempChIPCoverageFile) == TRUE) {
			system(paste("rm", InpTempChIPCoverageFile))
		}

		# delete the temporary file
		if (file.exists(InpTempFitHiChIPLoopFile) == TRUE) {
			system(paste("rm", InpTempFitHiChIPLoopFile))
		}

	}) 	# end chromosome index loop

	## finally merge the results to the final output file
	for (chr_idx in (1:length(CHRLIST))) {
		chrName <- CHRLIST[chr_idx]
		temp_final_UnionLoopFile <- paste0(dirname(UnionLoopFile), '/temp_mastersheet_', chrName, '.bed')
		if (chr_idx == 1) {
			system(paste("cat", temp_final_UnionLoopFile, " > ", UnionLoopFeatureFile))
		} else {
			system(paste0("awk \'(NR>1)\' ", temp_final_UnionLoopFile, " >> ", UnionLoopFeatureFile))
		}
		## remove the temporary file
		system(paste("rm", temp_final_UnionLoopFile))
	}

	if (file.exists(tempchrfile)) {
		system(paste("rm", tempchrfile))
	}

	# ## remve temporary objects
	# if (exists("MergedIntTempData")) {
	# 	rm("MergedIntTempData")
	# }	
	# if (exists("AllLoop_BinDF")) {
	# 	rm("AllLoop_BinDF")
	# }
	# if (exists("RawCC_Categ")) {
	# 	rm("RawCC_Categ")
	# }
	# if (exists("QVal_Categ")) {
	# 	rm("QVal_Categ")
	# }
	# if (exists("InpTempData")) {
	# 	rm("InpTempData")
	# }
	# if (exists("mergeDF")) {
	# 	rm("mergeDF")
	# }
	# if (exists("rawccvec")) {
	# 	rm("rawccvec")
	# }
	# if (exists("qvec")) {
	# 	rm("qvec")
	# }
	# if (exists("Seg1_ChIP_Coverage")) {
	# 	rm("Seg1_ChIP_Coverage")
	# }
	# if (exists("Seg2_ChIP_Coverage")) {
	# 	rm("Seg2_ChIP_Coverage")
	# }
	# if (exists("ChIPCoverageData")) {
	# 	rm("ChIPCoverageData")
	# }
	gc()	

}	# end function

#=========================
# this function creates DESeq2 compatible count matrix file
# parameter:
# MasterSheetData: union of all loops with feature information
# RawCC_ColList: set of columns in MasterSheetData with contact count information
# CountDataFile: output file to store the count data
#=========================
Create_DESeq2_Compatible_Count_Matrix <- function(MasterSheetData, RawCC_ColList, CountDataFile) {

	# dump the count matrix 
	write.table(MasterSheetData[, RawCC_ColList], CountDataFile, row.names=F, col.names=T, quote=F, sep="\t")

}	# end function

#======================================================
# this function returns a vector (equal to the size of number of loops)
# which counts for each loop the number of replicates showing significant interaction
# parameters:
# IntDF: loops from all replicates - DESEQ input structure like 
# CCCols: columns storing raw contact count
# QValCols: columns storing q-values
# FDR_thr_Loop: FDR threshold for FitHiChIP loop significance (default = 0.01)
#======================================================
GetCntNumRepVec <- function(IntDF, CCCols, QValCols, FDR_thr_Loop=0.01) {	
	CntVec <- matrix(0, nrow=nrow(IntDF), ncol=1)
	for (i in (1:length(CCCols))) {
		# get the loop indices such that 
		# current replicate is significant
		curr_cc_col <- CCCols[i]
		curr_q_col <- QValCols[i]
		Idx <- which((IntDF[, curr_cc_col] > 0) & (IntDF[, curr_q_col] <= FDR_thr_Loop))
		# increment the counter of those loops
		CntVec[Idx] <- CntVec[Idx] + 1
	}
	return(CntVec)
}

#===============
## this function performs DESeq2 
## parameters:
## DiffLoopDir: output directory
## InpTableFile: Grouping and covariate information table
## CountData: DESeq2 compatible count data 
## SampleInfoFile: file to contain the grouping information for DESeq2 execution
## DESeq_ResFile: DESeq2 results output file
## Inp_DiffModel: 0: use all loops, 1: use distance stratification
## suffixStr: suffix string at the plot files
#===============
Perform_DESeq2 <- function(DiffLoopDir, InpTableFile, CountData, SampleInfoFile, 
							DESeq_ResFile, Inp_DiffLoopModel, DistStrat, suffixStr, FDR_Th_DESeq) {

	CountDataColNames <- colnames(CountData)

	## remove the existing DESeq_ResFile
	if (file.exists(DESeq_ResFile)) {
		system(paste("rm", DESeq_ResFile))
	}

	## read the InpTableFile, containing all metadata information
	## write the DESeq2 copatible condition file
	InpTableData <- data.table::fread(InpTableFile, header=T)
	if (ncol(InpTableData) > 2) {
		coldata <- cbind.data.frame(data.frame(Name=CountDataColNames), 
									InpTableData[, c(2:ncol(InpTableData))])
		colnames(coldata) <- c("Name", "Condition", 
								paste0("Covariate_", seq(1, (ncol(coldata)-2))))
	} else {
		coldata <- cbind.data.frame(data.frame(Name=CountDataColNames), 
									InpTableData[, 2])
		colnames(coldata) <- c("Name", "Condition")
	}
	write.table(coldata, SampleInfoFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)

	## define the category list
	CategoryList <- as.vector(unique(InpTableData[, 2]))
	cat(sprintf("\n\n ******** CategoryList : %s **** \n\n", paste(CategoryList, collapse=" ")))

	## read the condition file according to the DESeq2 copatible format
	# coldata <- coldata[,c("Condition")]	#, "Type")]
	coldata <- read.table(SampleInfoFile, header=T, row.names=1)
	x <- colnames(coldata)

	# construct a DESeqDataSet using the given count matrix and the generated condition variable 
	## Note: we use dynamic formula
	## check https://support.bioconductor.org/p/86455/
	if (ncol(coldata) == 1) {
		formula_str <- paste("~", x)
	} else {
		formula_str <- paste("~", x[1], paste(paste("+", x[2:length(x)]), collapse=" "), collapse=" ")
	}
	cat(sprintf("\n\n ******** DESeq2 formula string : %s **** \n\n", formula_str))

	# Note: only one factor, namely the category information, is used to generate design matrix
	#design= ~ Condition + Type	
	dds <- DESeqDataSetFromMatrix(countData=CountData, 
									colData=coldata, 
									design=formula(formula_str))

	cat(sprintf("\n ********** performing DESeq on dds ********** \n"))

	## old code 
	## commented - sourya
	# ddsMF <- DESeq2::DESeq(dds, parallel=T)

	##========= new code - error handling
	## using try-catch module of R		
	out <- tryCatch({    			
    	ddsMF <- DESeq2::DESeq(dds, parallel=T)
	}, error=function(cond) {		
        cat(sprintf("\n\n ****** DESeq2 model error **** \n"))
    })

	## condition to check if error occurred
	if ((any(class(out) == "error") == "FALSE") & (exists("ddsMF"))) {

		PlotDir <- paste0(DiffLoopDir, '/Plots_DESeq2')
		system(paste("mkdir -p", PlotDir))

		## commented - sourya
		if (0) {
			# also plot dispersion estimates 
			DispEstPlotFile <- paste0(PlotDir, '/Dispersion_Estimate_Plot_', suffixStr, '.pdf')
			pdf(DispEstPlotFile, width=8, height=6)
			DESeq2::plotDispEsts(ddsMF)
			dev.off()
		}	# end if

		## default DESEQ results
		cat(sprintf("\n ********** DESeq results - contrast"))
		## check https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt
		DESEQRes <- DESeq2::results(ddsMF, parallel=T, 
									independentFiltering=FALSE, 
									contrast=c('Condition', CategoryList[1], CategoryList[2]))

		## commented - sourya
		if (0) {
			## MA plot using DEseq results
			MAPlotFile <- paste0(PlotDir, '/MA_Plot_', suffixStr, '.pdf')
			pdf(MAPlotFile, width=8, height=6)
			# DESeq2::plotMA(DESEQRes, ylim=c((-1) * FOLD_Change_Thr, FOLD_Change_Thr), xlim=c(1,10000), main=methodstr)
			DESeq2::plotMA(DESEQRes, main='MA_Plot')
			dev.off()
		}	# end if 

		## selective fields - so that DESeq2 and EdgeR outputs have similar format
		FinalDESEQRes <- cbind.data.frame(DESEQRes$log2FoldChange, 
											DESEQRes$baseMean, 
											DESEQRes$pvalue, 
											DESEQRes$padj)
		colnames(FinalDESEQRes) <- c('logFC', 'baseMean', 'pvalue', 'FDR')

		##==================
		## perform IHW if distance stratification is not employed
		##==================
		if (DistStrat == 0) {

			## diagnostic plot showing the utility of "baseMean" as a covariate
			Covariate_Plotfile <- paste0(PlotDir, '/pValue_vs_baseMean_', suffixStr, '.png')
			dataDF <- data.frame(pvalue = DESEQRes$pvalue, 
									covariate = rank(DESEQRes$baseMean)/nrow(DESEQRes), covariate_type="baseMean") 
			p <- ggplot(dataDF, aes(x = covariate, y = -log10(pvalue))) + geom_hex(bins = 100) + facet_grid( . ~ covariate_type) + ylab(expression(-log[10]~p))
			ggsave(Covariate_Plotfile, plot = p, width=6, height=6)

			## compute weighted p-values using covariate "baseMean" - use IHW 
			ihwRes <- as.data.frame(
						IHW::ihw(pvalues=FinalDESEQRes$pvalue, 
						covariates=FinalDESEQRes$baseMean, 
						alpha=FDR_Th_DESeq))	#0.1
			## append the adjusted p-value from IHW in the DESeq2 results
			FinalDESEQRes$FDR_IHW <- ihwRes$adj_pvalue
		}		

		##==================
		## write the DESeq2 results
		##==================
		write.table(FinalDESEQRes, DESeq_ResFile, row.names=F, col.names=T, quote=F, sep="\t", append=F)

	}	# end check error condition

	## remve temporary objects
	if (exists("dds")) {
		rm("dds")
	}
	if (exists("ddsMF")) {
		rm("ddsMF")
	}
	if (exists("DESEQRes")) {
		rm("DESEQRes")
	}
	if (exists("FinalDESEQRes")) {
		rm("FinalDESEQRes")
	}
	gc()	

}	# end function

#================================
## converts loops to WashU browser format (new)
#================================
ConvertLoopsWashU <- function(inpfile, QVal_ColList) {
	n <- GetNumLines(inpfile)
	if (n > 1) {
		outfile <- gsub(".bed", "_WashU.bed", inpfile)		
		# system(paste0("awk \'{if (NR>1) {print $1\"\t\"(($2+$3)/2-1)\"\t\"(($2+$3)/2+1)\"\t\"$4\":\"(($5+$6)/2-1)\"-\"(($5+$6)/2+1)\",10\t\"(NR-1)\"\t.\"}}\' ", inpfile, " | sort -k1,1 -k2,2n > ", outfile))
		## for each loop, use the minimum FDR
		system(paste0("cut -f1-6,", paste(QVal_ColList, collapse=","), " ", inpfile, " | awk \'{if (NR>1) {v=1; for (i=7;i<=NF;i++){v=(v<$i)?v:$i}; if (v==0) {q=20} else {q=(-log(v)/log(10))}; print $1\"\t\"(($2+$3)/2-1)\"\t\"(($2+$3)/2+1)\"\t\"$4\":\"(($5+$6)/2-1)\"-\"(($5+$6)/2+1)\",\"q\"\t\"(NR-1)\"\t.\"}}\' - | sort -k1,1 -k2,2n > ", outfile))
		if (file.exists(paste0(outfile, '.gz'))) {
			system(paste0("rm ", outfile, ".gz"))	
		}
		system(paste("bgzip", outfile))
		if (file.exists(paste0(outfile, '.gz.tbi'))) {
			system(paste0("rm ", outfile, ".gz.tbi"))	
		}
		system(paste("tabix -p bed", paste0(outfile, ".gz")))
	}
}

#===============================
# function to categorize the exclusive loops (significant in one category) according to the 
# ChIP-seq differential coverage information
#===============================
Categorize_ExclusiveLoops <- function(MainOutDir, InpLoopFile, CategoryList, QVal_ColList, HiC_Data_Process, CovThr) {

	ConvertLoopsWashU(InpLoopFile, QVal_ColList)

	## create two files - one containing the upegulated interactions in the first category
	## and one containing the upegulated interactions in the second category 		
	UpReg_Cat1_File <- gsub(".bed", paste0("_Upreg_", CategoryList[1], '.bed'), InpLoopFile)
	UpReg_Cat2_File <- gsub(".bed", paste0("_Upreg_", CategoryList[2], '.bed'), InpLoopFile)

	system(paste0("awk \'((NR==1) || (($(NF-1)>0) && ($NF==0)))\' ", InpLoopFile, " > ", UpReg_Cat1_File))
	system(paste0("awk \'((NR==1) || (($(NF-1)==0) && ($NF>0)))\' ", InpLoopFile, " > ", UpReg_Cat2_File))

	ConvertLoopsWashU(UpReg_Cat1_File, QVal_ColList)
	ConvertLoopsWashU(UpReg_Cat2_File, QVal_ColList)

	## if the ChIP-seq data is provided, we also classify these loops according to the 1D invariance 
	## of the underlying loops
	if (HiC_Data_Process == FALSE) {

		categoryDir <- paste0(MainOutDir, '/Differential_ChIP_Cov_pct_', CovThr, '/ND_ND')
		system(paste("mkdir -p", categoryDir))
		allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($11==\"ND\") && ($12==\"ND\"))\' ", InpLoopFile, ">", allLoop_file))
		ConvertLoopsWashU(allLoop_file, QVal_ColList)	
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)>0) && ($NF==0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[1], '.bed')))
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)==0) && ($NF>0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[2], '.bed')))

		categoryDir <- paste0(MainOutDir, '/Differential_ChIP_Cov_pct_', CovThr, '/LD_ND')
		system(paste("mkdir -p", categoryDir))
		allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
		system(paste("awk -F\'[\t]\' \'(NR==1) || ((($11==\"LD\") && ($12==\"ND\")) || (($11==\"ND\") && ($12==\"LD\")))\' ", InpLoopFile, ">", allLoop_file))
		ConvertLoopsWashU(allLoop_file, QVal_ColList)	
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)>0) && ($NF==0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[1], '.bed')))
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)==0) && ($NF>0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[2], '.bed')))

		categoryDir <- paste0(MainOutDir, '/Differential_ChIP_Cov_pct_', CovThr, '/LD_LD')
		system(paste("mkdir -p", categoryDir))
		allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($11==\"LD\") && ($12==\"LD\"))\' ", InpLoopFile, ">", allLoop_file))
		ConvertLoopsWashU(allLoop_file, QVal_ColList)	
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)>0) && ($NF==0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[1], '.bed')))
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)==0) && ($NF>0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[2], '.bed')))

		categoryDir <- paste0(MainOutDir, '/Differential_ChIP_Cov_pct_', CovThr, '/HD_HD')
		system(paste("mkdir -p", categoryDir))
		allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($11==\"HD\") && ($12==\"HD\"))\' ", InpLoopFile, ">", allLoop_file))
		ConvertLoopsWashU(allLoop_file, QVal_ColList)	
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)>0) && ($NF==0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[1], '.bed')))
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)==0) && ($NF>0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[2], '.bed')))

		categoryDir <- paste0(MainOutDir, '/Differential_ChIP_Cov_pct_', CovThr, '/HD_LD_or_ND')
		system(paste("mkdir -p", categoryDir))
		allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
		system(paste("awk -F\'[\t]\' \'(NR==1) || ((($11==\"HD\") && ($12!=\"HD\")) || (($11!=\"HD\") && ($12==\"HD\")))\' ", InpLoopFile, ">", allLoop_file))
		ConvertLoopsWashU(allLoop_file, QVal_ColList)	
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)>0) && ($NF==0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[1], '.bed')))
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)==0) && ($NF>0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[2], '.bed')))
	}
}

#===============================
# function to categorize the differential loops according to the 
# ChIP-seq differential coverage information
#===============================
Categorize_DiffLoops <- function(MainOutDir, InpLoopFile, CategoryList, QVal_ColList, DiffLoopModel, HiC_Data_Process, CovThr, UpDownLoops=1) {

	ConvertLoopsWashU(InpLoopFile, QVal_ColList)
	
	if (UpDownLoops == 1) {
		## create two files - one containing the upegulated interactions in the first category
		## and one containing the upegulated interactions in the second category 		
		UpReg_Cat1_File <- gsub(".bed", paste0("_Upreg_", CategoryList[1], '.bed'), InpLoopFile)
		UpReg_Cat2_File <- gsub(".bed", paste0("_Upreg_", CategoryList[2], '.bed'), InpLoopFile)

		## for DESeq2, fold change > 0 means first category is upregulated 
		## for edgeR, fold change > 0 means second category is upregulated 
		if (DiffLoopModel == 0) {
			## DESeq2
			system(paste0("awk \'((NR==1) || ($(NF-6)>0))\' ", InpLoopFile, " > ", UpReg_Cat1_File))
			system(paste0("awk \'((NR==1) || ($(NF-6)<0))\' ", InpLoopFile, " > ", UpReg_Cat2_File))		
		} else {
			## edgeR
			system(paste0("awk \'((NR==1) || ($(NF-6)>0))\' ", InpLoopFile, " > ", UpReg_Cat2_File))
			system(paste0("awk \'((NR==1) || ($(NF-6)<0))\' ", InpLoopFile, " > ", UpReg_Cat1_File))
		}
		ConvertLoopsWashU(UpReg_Cat1_File, QVal_ColList)
		ConvertLoopsWashU(UpReg_Cat2_File, QVal_ColList)
	}

	## if the ChIP-seq data is provided, we also classify these loops according to the 1D invariance 
	## of the underlying loops
	if (HiC_Data_Process == FALSE) {

		categoryDir <- paste0(MainOutDir, '/Differential_ChIP_Cov_pct_', CovThr, '/ND_ND')
		system(paste("mkdir -p", categoryDir))
		allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($11==\"ND\") && ($12==\"ND\"))\' ", InpLoopFile, ">", allLoop_file))
		ConvertLoopsWashU(allLoop_file, QVal_ColList)	
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)>0) && ($NF==0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[1], '.bed')))
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)==0) && ($NF>0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[2], '.bed')))

		categoryDir <- paste0(MainOutDir, '/Differential_ChIP_Cov_pct_', CovThr, '/LD_ND')
		system(paste("mkdir -p", categoryDir))
		allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
		system(paste("awk -F\'[\t]\' \'(NR==1) || ((($11==\"LD\") && ($12==\"ND\")) || (($11==\"ND\") && ($12==\"LD\")))\' ", InpLoopFile, ">", allLoop_file))
		ConvertLoopsWashU(allLoop_file, QVal_ColList)	
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)>0) && ($NF==0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[1], '.bed')))
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)==0) && ($NF>0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[2], '.bed')))

		categoryDir <- paste0(MainOutDir, '/Differential_ChIP_Cov_pct_', CovThr, '/LD_LD')
		system(paste("mkdir -p", categoryDir))
		allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($11==\"LD\") && ($12==\"LD\"))\' ", InpLoopFile, ">", allLoop_file))
		ConvertLoopsWashU(allLoop_file, QVal_ColList)	
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)>0) && ($NF==0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[1], '.bed')))
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)==0) && ($NF>0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[2], '.bed')))

		categoryDir <- paste0(MainOutDir, '/Differential_ChIP_Cov_pct_', CovThr, '/HD_HD')
		system(paste("mkdir -p", categoryDir))
		allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($11==\"HD\") && ($12==\"HD\"))\' ", InpLoopFile, ">", allLoop_file))
		ConvertLoopsWashU(allLoop_file, QVal_ColList)	
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)>0) && ($NF==0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[1], '.bed')))
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)==0) && ($NF>0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[2], '.bed')))

		categoryDir <- paste0(MainOutDir, '/Differential_ChIP_Cov_pct_', CovThr, '/HD_LD_or_ND')
		system(paste("mkdir -p", categoryDir))
		allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
		system(paste("awk -F\'[\t]\' \'(NR==1) || ((($11==\"HD\") && ($12!=\"HD\")) || (($11!=\"HD\") && ($12==\"HD\")))\' ", InpLoopFile, ">", allLoop_file))
		ConvertLoopsWashU(allLoop_file, QVal_ColList)	
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)>0) && ($NF==0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[1], '.bed')))
		system(paste("awk -F\'[\t]\' \'(NR==1) || (($(NF-1)==0) && ($NF>0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[2], '.bed')))
	}

}	# end function

#=========================
# this function creates EdgeR compatible count matrix file
# parameters:
# MasterSheetData: union of all loops with feature information
# RawCC_ColList: set of columns in MasterSheetData with contact count information
# DiffLoopDir: output directory to contain the count matrix
# CountDataFile: output count data file
#=========================
Create_EdgeR_Compatible_Count_Matrix <- function(MasterSheetData, RawCC_ColList, DiffLoopDir, CountDataFile) {

	## list of chromosomes for the current input mastersheet file
	CHRLIST_NAMENUM <- as.vector(unique(MasterSheetData[,1]))

	# number of interacting bins considered so far
	numBin <- 0	

	for (chrIdx in (1:length(CHRLIST_NAMENUM))) {
		
		## extract master sheet data for the current chromosome
		CurrChr <- CHRLIST_NAMENUM[chrIdx]
		MasterSheetData_CurrChr <- MasterSheetData[which(MasterSheetData[, 1] == CurrChr), ]
		
		# count matrix to be used for EdgeR - using raw count data (RawCC)
		CountData <- MasterSheetData_CurrChr[, RawCC_ColList]
		CountDataColNames <- colnames(CountData)
		
		# extract unique interacting bins
		Bin1DF <- unique(MasterSheetData_CurrChr[,1:3])
		colnames(Bin1DF) <- c("chr1", "start1", "end1")
		Bin2DF <- unique(MasterSheetData_CurrChr[,4:6])
		colnames(Bin2DF) <- c("chr1", "start1", "end1")
		IntervalData <- unique(rbind.data.frame(Bin1DF, Bin2DF))

		## also mention indices of the interacting bins - they should be unique
		IntervalData$Idx1 <- seq((numBin+1), (numBin+nrow(IntervalData)))

		# prepare a copy, for joining operation
		IntervalData2 <- data.frame(IntervalData)
		colnames(IntervalData2) <- c("chr2", "start2", "end2", "Idx2")

		# increment the variable "numBin"
		numBin <- numBin + nrow(IntervalData)

		# join the set of loops together with these bin index data frames
		finaldf <- dplyr::inner_join(MasterSheetData_CurrChr[, 1:6], IntervalData) %>% dplyr::inner_join(IntervalData2)
		cat(sprintf("\n after merging -  dimension of finaldf : %s X %s ", nrow(finaldf), ncol(finaldf)))

		# prepare the count data (for EdgeR structure) containing the indices as well
		CountData <- cbind.data.frame(finaldf$Idx1, finaldf$Idx2, finaldf$chr1, ((finaldf$start1 + finaldf$end1) / 2), finaldf$chr2, ((finaldf$start2 + finaldf$end2) / 2), CountData)
		colnames(CountData) <- c("Idx1", "Idx2", "chr1", "mid1", "chr2", "mid2", CountDataColNames)

		if (chrIdx == 1) {
			write.table(CountData, CountDataFile, row.names=F, col.names=T, quote=F, sep="\t", append=F)
		} else {
			write.table(CountData, CountDataFile, row.names=F, col.names=F, quote=F, sep="\t", append=T)
		}		

		## remove temporary objects
		if (exists("Bin1DF")) {
			rm("Bin1DF")
		}
		if (exists("Bin2DF")) {
			rm("Bin2DF")
		}
		if (exists("IntervalData")) {
			rm("IntervalData")
		}
		if (exists("IntervalData2")) {
			rm("IntervalData2")
		}
		if (exists("finaldf")) {
			rm("finaldf")
		}
		if (exists("CountData")) {
			rm("CountData")
		}
		gc()

	}	# end chromosome index loop

	# final summary of the total number of loops
	N <- GetNumLines(CountDataFile)
	cat(sprintf("\n\n ========> FINAL number of loops - EdgeR count matrix: %s \n\n", N))

}	# end function

#===================
# function to apply EdgeR on loops
## parameters:
## CountData_Orig: input count data for edgeR processing
## DiffLoopDir: output directory
## InpTableFile: Input categpry and covariate information
## SampleInfoFile: Dumped category information used in edgeR
## EdgeRModel: EdgeR model employed
## EdgeR_ResFile: output file to contain edgeR results
## Inp_DiffModel: 0: use all loops, 1: use distance stratification
## inp_suffixStr: suffix string
## bcv: by default 0.4 
#===================
Perform_EdgeR <- function(CountData_Orig, DiffLoopDir, InpTableFile, SampleInfoFile, EdgeR_ResFile, Inp_DiffLoopModel, DistStrat, suffixStr, FDR_Th_DESeq, bcv=0.4) {

	## reset
	if (file.exists(EdgeR_ResFile)) {
		system(paste("rm", EdgeR_ResFile))
	}

	## extract the counts from the input "CountData_Orig"
	CountData <- CountData_Orig[, 7:ncol(CountData_Orig)]
	CountDataColNames <- colnames(CountData)

	# directory to store the plots related to EdgeR
	PlotDir <- paste0(DiffLoopDir, '/Plots_EdgeR')
	system(paste("mkdir -p", PlotDir))

	## read the sample and group information (and optional covariate information)
	InpTableData <- data.table::fread(InpTableFile, header=T)

	CategoryList <- as.vector(unique(InpTableData[, 2]))
	cat(sprintf("\n\n ******** CategoryList : %s **** \n\n", paste(CategoryList, collapse=" ")))

	## replica counts for two categories
	ReplicaCount <- c(length(which(InpTableData[,2] == CategoryList[1])), length(which(InpTableData[,2] == CategoryList[2])))
	cat(sprintf("\n\n ******** ReplicaCount : %s **** \n\n", paste(ReplicaCount, collapse=" ")))
	
	coldata <- cbind.data.frame(data.frame(Name=CountDataColNames), InpTableData[, 2])
	colnames(coldata) <- c("Name", "Group")
	write.table(coldata, SampleInfoFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)

	## create the design matrix
	## the 0+ in the model formula eliminate the intercept column and just includes the covariate columns
	Group <- factor(coldata$Group)
	design <- model.matrix(~0+Group) 
	colnames(design) <- levels(Group)

	# create the EdgeR count data structure (7th column onwards)
	y <- DGEList(counts=CountData, group=Group)

	# ## Remove rows conssitently have zero or very low counts
	# keep <- filterByExpr(y)
	# y <- y[keep, keep.lib.sizes = FALSE]

	# for data with multiple replicates, better to normalize the DGEList object
	# if ((ReplicaCount[1] > 1) & (ReplicaCount[2] > 1)) {
	if (1) {
		y <- calcNormFactors(y, method="TMM")
	}
	
	if (Inp_DiffLoopModel >= 1 & Inp_DiffLoopModel <= 3) {
		##=======
		## for edgeR exactTest or GLM models, estimate the dispersions
		##=======
		tryCatch({    			
			# Estimate Common, Trended and Tagwise Negative Binomial dispersions 
			# by weighted likelihood empirical Bayes
			y <- estimateDisp(y, design)		
		}, error=function(cond) {		
			cat(sprintf("\n\n ****** edgeR model - estimateDisp - error **** \n"))
		})

		## according to the input "Inp_DiffLoopModel" (i.e. exactTest or GLM model specification)
		## estimate dispersions, tags, and significance
		if (Inp_DiffLoopModel == 1) {	
			## exact test - experiments with single factor - qCML method
			## check the number of replicates in both categories for estimating dispersions
			## comparison between two categories (Category 1 - Category 2)	
			if (exists("y")) {
				tryCatch({			
					if ((ReplicaCount[1] > 1) & (ReplicaCount[2] > 1)) {
						edgeR_model <- exactTest(y, dispersion="trended", pair=c(CategoryList[1], CategoryList[2]))
					} else {		
						## if at least one of the input categories do not have multiple replicates
						## we manually supply the dispersion
						edgeR_model <- exactTest(y, dispersion=bcv^2, pair=c(CategoryList[1], CategoryList[2]))
					}
				}, error=function(cond) {		
					cat(sprintf("\n\n ****** edgeR model - exactTest - error **** \n"))
				})				
			}
		} else if (Inp_DiffLoopModel == 3) {
			tryCatch({
				## quasi-likelihood F-tests using GLM
				if ((ReplicaCount[1] > 1) & (ReplicaCount[2] > 1)) {
					fit_model <- glmQLFit(y, design)
				} else {
					## if at least one of the input categories do not have multiple replicates
					## we manually supply the dispersion
					fit_model <- glmQLFit(y, design, dispersion=bcv^2)
				}
			}, error=function(cond) {		
				cat(sprintf("\n\n ****** edgeR model - glmQLFit - error **** \n"))
			})			
			## explicitly specify the contrast argument
			## comparison between two categories (Category 1 - Category 2)	
			if (exists("fit_model")) {
				tryCatch({		
					edgeR_model <- glmQLFTest(fit_model, contrast=c(1, -1))
				}, error=function(cond) {		
					cat(sprintf("\n\n ****** edgeR model - glmQLFTest - error **** \n"))
				})
			}
		} else if (Inp_DiffLoopModel == 2) {
			tryCatch({
				# likelihood ratio tests (LRT) using GLM
				if ((ReplicaCount[1] > 1) & (ReplicaCount[2] > 1)) {
					fit_model <- glmFit(y, design)
				} else {
					## if at least one of the input categories do not have multiple replicates
					## we manually supply the dispersion
					fit_model <- glmFit(y, design, dispersion=bcv^2)
				}
			}, error=function(cond) {		
				cat(sprintf("\n\n ****** edgeR model - glmFit - error **** \n"))
			})			
			## explicitly specify the contrast argument
			## comparison between two categories (Category 1 - Category 2)			
			if (exists("fit_model")) {				
				tryCatch({	
					edgeR_model <- glmLRT(fit_model, contrast=c(1, -1))	
				}, error=function(cond) {		
					cat(sprintf("\n\n ****** edgeR model - glmLRT - error **** \n"))
				})
			}
		} 

		## get the DE genes 
		## will be used for MA plot
		## and also the final results
		if (exists("edgeR_model")) {	
			edgeR_model_tags <- topTags(edgeR_model)
		}

		##===== plot statistics
		if (exists("y") & (ReplicaCount[1] > 1) & (ReplicaCount[2] > 1)) {

			## commented - sourya
			if (0) {

				## 1. plotMDS function
				## generates a plot in which distances between samples 
				## correspond to leading biological coefficient of variation (BCV) between those samples
				if (Inp_DiffLoopModel == 1) {		
					MDSPlotFile <- paste0(PlotDir, '/EdgeR_MDS_Plot_', suffixStr, '.pdf')
					pdf(MDSPlotFile, width=8, height=6)
					plotMDS(y)
					dev.off()
				}

				## 2. plotBCV function
				## dispersion estimates can be viewed in a BCV plot
				if (Inp_DiffLoopModel == 1) {		
					BCVPlotFile <- paste0(PlotDir, '/EdgeR_BCV_Plot_', suffixStr, '.pdf')
					pdf(BCVPlotFile, width=8, height=6)
					plotBCV(y)
					dev.off()
				}

				## 3. plotMD function
				## applicable when both cateopries have multiple replicates
				## Plot log-fold change against log-counts per million, with DE genes highlighted	
				MDPLotFile <- paste0(PlotDir, '/EdgeR_MD_Plot_', suffixStr, '.pdf')
				pdf(MDPLotFile, width=8, height=6)
				plotMD(edgeR_model, hl.cex = 0.1)
				abline(h=c(-1, 1), col="blue")
				dev.off()

				## 4. plotSmear function
				## smear plot - MA for groups of samples
				SmearPlotFile <- paste0(PlotDir, '/EdgeR_Smear_Plot_', suffixStr, '.pdf')
				pdf(SmearPlotFile, width=8, height=6)
				plotSmear(edgeR_model, pair=c(CategoryList[1], CategoryList[2]), de.tags=edgeR_model_tags)
				dev.off()

				## 5. plotQLDisp function
				## plot the QL dispersions
				## using output from glmQLFit function
				if (Inp_DiffLoopModel == 3) {
					QLDispersionPlotFile <- paste0(PlotDir, '/EdgeR_QL_Dispersion_Plot_', suffixStr, '.pdf')
					pdf(QLDispersionPlotFile, width=8, height=6)
					plotQLDisp(fit_model)
					dev.off()
				}

			} 	# end if 	
		}

		#===============
		# for each of the DE tag finding techniques (qCML or CR methods)
		# convert the p-values to the corresponding q-values, using BH correction
		#===============
		if (exists("edgeR_model")) {	
			QVal <- p.adjust(edgeR_model$table$PValue, method = "BH")
			EdgeRRes <- cbind.data.frame(edgeR_model$table$logFC, edgeR_model$table$logCPM, edgeR_model$table$PValue, QVal)		
			colnames(EdgeRRes) <- c('logFC', 'logCPM', 'pvalue', 'FDR')
		}

	} else if (Inp_DiffLoopModel == 4) {
		##=======
		## for Wilcoxon rank sum test on edgeR TMM normalized data
		##=======
		time_val1 = Sys.time()

		## get counts-per-million (cpm)
		count_norm <- as.data.frame(cpm(y))
		
		## get average log CPM for each row
		avglogCPMVal <- edgeR::aveLogCPM(y)

		time_val2 = Sys.time()
		cat(sprintf("\n Start of Wilcoxon test p-value computation - time elapsed : %s \n", (time_val2 - time_val1)))

		ncore <- detectCores()

		## Run the Wilcoxon rank-sum test for each entry
		pvalues_list <- parallel:::mclapply( 1:nrow(count_norm), mc.cores = ncore, function(i) {
		# pvalues <- sapply( 1:nrow(count_norm), function(i) {
		
			data <- cbind.data.frame(loop = as.numeric(t(count_norm[i,])), Group)
			p <- wilcox.test(loop~Group, data)$p.value
			return(p)
		})

		## as mclapply returns a list, we need to unlist to get the vector of p-values
		pvalues = unlist(pvalues_list)

		time_val3 = Sys.time()
		cat(sprintf("\n End of Wilcoxon test p-value computation - time elapsed : %s \n", (time_val3 - time_val2)))

		## calculate FDR (using BH correction)
		fdrval <- p.adjust(pvalues, method = "BH")

		time_val4 = Sys.time()
		cat(sprintf("\n FDR of Wilcoxon test p-value computation - time elapsed : %s \n", (time_val4 - time_val3)))

		## Calculate the fold-change for each loop
		conditionsLevel <- levels(Group)
		dataCon1 <- count_norm[,c(which(Group==conditionsLevel[1]))]
		dataCon2 <- count_norm[,c(which(Group==conditionsLevel[2]))]
		logFCVal <- log2((rowMeans(dataCon2) + 1) / (rowMeans(dataCon1) + 1))

		EdgeRRes <- data.frame(logFC=logFCVal, logCPM=avglogCPMVal, pvalue=pvalues, FDR=fdrval)
	}

	if (exists("EdgeRRes")) {

		## perform IHW if distance stratification is not employed
		if (DistStrat == 0) {
			## diagnostic plot showing the utility of "logCPM" as a covariate
			Covariate_Plotfile <- paste0(PlotDir, '/pValue_vs_logCPM_', suffixStr, '.png')
			dataDF <- data.frame(pvalue = EdgeRRes$pvalue, 
								covariate = rank(EdgeRRes$logCPM)/nrow(EdgeRRes), 
								covariate_type="logCPM") 
			p <- ggplot(dataDF, aes(x = covariate, y = -log10(pvalue))) + geom_hex(bins = 100) + facet_grid( . ~ covariate_type) + ylab(expression(-log[10]~p))
			ggsave(Covariate_Plotfile, plot = p, width=6, height=6)		

			## compute weighted p-values using covariate "logCPM" - use IHW
			ihwRes <- as.data.frame(IHW::ihw(pvalues=EdgeRRes$pvalue, 
				covariates=EdgeRRes$logCPM, alpha=FDR_Th_DESeq))
			## append the adjusted p-value from IHW in the edgeR results
			EdgeRRes$FDR_IHW <- ihwRes$adj_pvalue
		}

		# ## As the "EdgeRRes" contains only the filtered entries from "filterByExpr"
		# ## generate the complete results
		# df1 <- data.frame(idx=seq(1,nrow(CountData)))
		# EdgeRRes$idx <- which(keep == TRUE)
		# cat(sprintf("\n length df1 idx : %s length EdgeRRes idx : %s ", length(df1$idx), length(EdgeRRes$idx)))
		# FinalEdgeRRes <- dplyr::full_join(df1, EdgeRRes)
		
		# ## drop the idx column
		# FinalEdgeRRes <- subset(FinalEdgeRRes, select = -c(idx))
		
		## write the edgeR results
		write.table(EdgeRRes, EdgeR_ResFile, row.names=F, col.names=T, quote=F, sep="\t", append=F)		

	}

	## remove temporary objects
	if (exists("CountData")) {
		rm("CountData")
	}
	if (exists("y")) {
		rm("y")
	}
	if (exists("design")) {
		rm("design")
	}
	if (exists("edgeR_model")) {
		rm("edgeR_model")
	}
	if (exists("edgeR_model_tags")) {
		rm("edgeR_model_tags")
	}
	if (exists("EdgeRRes")) {
		rm("EdgeRRes")
	}
	if (exists("FinalEdgeRRes")) {
		rm("FinalEdgeRRes")
	}
	gc()
	
}	# end function


