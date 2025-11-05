DiffHiChIP
=============

Framework to detect statistically significant differential loops between two conditions / categories for HiChIP (and other compatible 3C protocols like PLAC-seq, ChIA-PET, or Hi-C) data.

DiffHiChIP is now published at Cell Reports Methods (https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(25)00250-4)

If you are using DiffHiChIP, please cite:

Bhattacharyya et al., DiffHiChIP: Identifying differential chromatin contacts from HiChIP data, Cell Reports Methods (2025), https://doi.org/10.1016/j.crmeth.2025.101214

Authors (and Contacts)
=======================

Sourya Bhattacharyya (sourya@lji.org)

Daniela Salgado Figueroa (dsfigueroa@lji.org)

Ferhat Ay (ferhatay@lji.org)


Description (and Novelty)
==========================

Existing differential HiChIP loop callers like [FitHiChIP](https://www.nature.com/articles/s41467-019-11950-y), [diffloop](https://academic.oup.com/bioinformatics/article/34/4/672/4282662) and [HiC-DC+](https://www.nature.com/articles/s41467-021-23749-x) employ count-based tests from [edgeR](https://academic.oup.com/bioinformatics/article/26/1/139/182458) or [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) which are mainly used for RNA-seq analysis.

These methods, by default, cannot model the exponential distance decay of HiChIP contacts and fail to detect differences in long-range loops (>500Kb). To counter the distance decay, stratification approaches perform equal-distance binning of contacts by their genomic distance and process individual bins separately, but still exhibit lower statistical power in detecting long-range loops.

DiffHiChIP is the first comprehensive framework to assess these differential HiChIP loop calling models. 

1. It incorporates both DESeq2 and edgeR with false discovery rate (FDR) and [independent hypothesis weighting (IHW)](https://bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html) corrected p-values. 

2. Implements stratification by equal-occupancy binning (similar contacts per bin) to assign higher statistical power to long-range loops.

3. DiffHiChIP further incorporates edgeR with generalized linear model (GLM), and defines four additional models employing either of quasi-likelihood F-test, likelihood ratio test, or fold-change specific thresholds (TREAT).

	3a. The GLM-based models show higher precision in recovering short and long-range differential contacts, which are highly supported by the respective Hi-C backgrounds.

4. DiffHiChIP further categorizes differential loops by the underlying ChIP-seq differences of the interacting bins, and the significance of loops in individual samples.


Installation
==============

To execute DiffHiChIP, user needs to install the following environments and packages:

1. R (we recommend using R version 4.3 and the latest Bioconductor so to use the latest versions of DESeq2 and edgeR)

2. Following R libraries need to be installed:

	optparse, data.table, GenomicRanges, DESeq2, edgeR, dplyr, ggplot2, IHW, ashr, BiocParallel, parallel

3. If user needs to perform loop calling on HiChIP data first, we recommed using FitHiChIP for HiChIP loop calling. Please check the [manuscript](https://www.nature.com/articles/s41467-019-11950-y) and [GitHub page](https://github.com/ay-lab/FitHiChIP) and a step by step [documentation](https://ay-lab.github.io/FitHiChIP/html/index.html) of FitHiChIP to execute it.


Execution
==========

1. User first needs to edit the text file *scripts/CD4N_vs_CD8N_InputTable.txt* (user can edit the file name as well) which contains the paths and categories of the input samples (HiChIP loop files).

**Note** The FitHiChIP loops should be the set of all possible contacts along with their statistical significance values, not only the significant contacts. In FitHiChIP, the file names are *PREFIX*.interactions_FitHiC.bed (check the [documentation](https://ay-lab.github.io/FitHiChIP/html/usage/output.html) for details).

2. Then check the file *scripts_DiffLoop_CD4N_vs_CD8N.sh* (user can edit the file name as well) and edit the file paths and other configuration parameters. Then run *qsub scripts_DiffLoop_CD4N_vs_CD8N.sh* to execute the job.

**Note** The chromosome size file needs to be provided according to the reference genome employed.


Configuration parameters
=========================

The above mentioned script executes various options of DiffHiChIP all at once. However, not all options / configuration options may be needed. Below are the list of configuration options that users need to look into before running the above mentioned script (and also to interpret the outputs for different settings).

The script basically invokes the main source code *src/DiffLoop_Code_FitHiChIP.r* along with different options, which are summarized below.

Usage: DiffLoop_Code_FitHiChIP.r [options]

Options:
	
	--InpTable=INPTABLE
		Input table file with loop file, covariate information. Default = NULL. Mandatory parameter.

	--OutDir=OUTDIR
		Base Output directory. Mandatory parameter.

	--ChrSizeFile=CHRSIZEFILE
		Chromosome size file corresponding to reference genome. Mandatory parameter.

	--MidList=MIDLIST
		Comma or colon separated list of boolean values (0/1) corresponding to the input loop files. 
						The list should have equal length as the --AllLoopList parameter.
						The order should be same as the files mentioned in the --AllLoopList parameter.
						1: loop file has midpoints of each bin 
						0: loop file has bin interval (chr, start, end). 
						Instead of a list, a single value (0 or 1) can also be provided. 
						In such a case, all files will be assumed to have the same configuration. 
						Default: list of all zeroes (means all loop files have bin interval information).

	--BinSize=BINSIZE
		Bin size / resolution of the input loops. 
					Used only if all the entries in --MidList option is 1, thus all loop files contain only midpoints of the bin intervals. Mandatory parameter in such a case. Default = 5000 means 5 Kb.

	--CCColList=CCCOLLIST
		Comma or colon separated list of integers depicting the column numbers of the input loop files storing the raw contact count. 
						The list should have equal length as the --AllLoopList parameter.
						The order should be same as the files mentioned in the --AllLoopList parameter.
						If all the input files have similar settings, user can specify only one value.
						Default = 7, same as FitHiChIP output format.

	--QColList=QCOLLIST
		Comma or colon separated list of integers depicting the column numbers of the input loop files storing the q-value (FDR). 
						The list should have equal length as the --AllLoopList parameter.
						The order should be same as the files mentioned in the --AllLoopList parameter.
						If all the input files have similar settings, user can specify only one value.
						Default is the last column, same as FitHiChIP output format.

	--FDRThr=FDRTHR
		FDR (q-value) threshold used for determining Hi-C or HiChIP significant loops. Default = 0.01.

	--Model=MODEL
		Differential analysis model.
						0: DESeq2
						1: edgeR exact Test
						2: edgeR GLM with likelihood ratio tests (glmFit + glmLRT)
						3: edgeR GLM with Quasi-likelihood (QL) F-Test (glmQLFit + glmQLFTest)
						4: Wilcoxon rank sum test on edgeR TMM normalized count
						Default = 1
						** Note: Models 0, 2, 3, are applicable if both categories have > 1 replicates.

	--PreFilt=PREFILT
		Binary value. 
						0: all FitHiChIP loops from all samples are used as background.
						1: Union of loops having FitHiChIP FDR < 0.1 in at least one input sample are used as background. 
						Default = 0, means no filtering is done.

	--DiffFDRThr=DIFFFDRTHR
		FDR threshold for differential loops. Default = 0.05

	--LFCThr=LFCTHR
		log2 fold change threshold for differential loops. Default = 1.

	--DistStrat=DISTSTRAT
		Binary value. Possible values are: 
				[0]: No distance stratification. Use FDR and IHW corrected FDR for differential analysis.
				[1]: Use distance stratification (equal occupancy binning) and FDR. 
				Default = 0, means no distance stratification is employed.

	--DSBinSize=DSBINSIZE
		Applicable if the parameter --DistStrat is 1. 
						Specifies the bin size (in bytes) employed for stratification. 
						Default = 10000 means 10 Kb is the equal occupancy bin size.

	--bcv=BCV
		If edgeR is used with one replicate in at least one input category, this value is the square-root-dispersion. 
					For datasets arising from well-controlled experiments are 0.4 for human data, 
					0.1 for data on genetically identical model organisms, 
					or 0.01 for technical replicates. 
					For details, see the edgeR manual. 
					By default, the value is set as 0.4.

	--ChIPAlignFileList=CHIPALIGNFILELIST
		Applicable only for HiChIP data processing. 
					Comma or colon separated list of ChIP-seq alignment files. 
					Either BAM or bedgraph formatted file is supported. 
					Default = NULL. 
					User can: (a) either provide two files, one for each category, or 
					(b) provide ChIP seq alignment files for each sample. 
					If this parameter is empty, differential loops would be analyzed alright 
					but their underlying 1D (ChIP) differences would not.

	--CovThr=COVTHR
		Applicable only for HiChIP data processing. 
				Threshold signifying the maximum allowed deviation of ChIP coverage between two categories for which the bins would be considered as 1D invariant (ND). 
				Default (or if user specifies value <=0 or > 100) = 25, 
				means 25% deviation is set as maximum. 
				If user chooses 50, 50% maximum ChIP seq coverage deviation would be allowed.

	--SuffixStr=SUFFIXSTR
		This custom string can be used to denote specific input conditions. 
				Such as DESeq2_Covariate. 
				If --Overwrite is 0, using different suffix strings would help to store outputs in different folders 
				without re-computing majority of the differential analysis. 
				Default = NULL.

	--Overwrite=OVERWRITE
		Binary variable. If 1, overwrites existing results. 
				Default = 0.

	-h, --help
		Show this help message and exit


Output description
===================

Within the mentioned output directory **OUTDIR**, the output differential loops are stored in the following directory structure, depending on the input parameters:

$${\color{blue}BACKGROUND}/{\color{blue}STRAT}/{\color{green}MODEL}/{\color{brown}FDRTYPE}/**DiffLoops_ALL.bed**$$

Details of the directory structure:

$${\color{blue}BACKGROUND}$$: 

Denotes the background set of contacts employed for differential loop calling.

	1. If --PreFilt=0 (default), the directory name is *Background_AllContacts*. 

	2. Otherwise, the folder name is *Background_FilteredContacts*.

$${\color{blue}STRAT}$$: 

Denotes whether distance stratification is employed or not.

	1. If --DistStrat=0, no distance stratification is employed. The folder name is *No_DistStrat*. 

	2. Otherwise, for distance stratification (--DistStrat=1), the directory name is *With_DistStrat_{DSBINSIZE}* where *DSBINSIZE* is the bin size employed for stratification (see the parameter --DSBinSize)

$${\color{green}MODEL}$$: 

Specifies the differential analysis model.

1. If the parameter --Model=0, DESeq2 is employed. The directory name is *DESeq2*

2. If the parameter --Model=1, edgeR exact test is employed. The directory name is *EdgeR_exactTest*

3. If the parameter --Model=2, edgeR GLM with likelihood ratio test (LRT) is employed. The directory name is *EdgeR_glmLRT*

4. If the parameter --Model=3, edgeR GLM with quasi-likelihood F test is employed. The directory name is *EdgeR_glmQLFTest*

5. If the parameter --Model=4, Wilcoxon rank sum test with edgeR TMM normalized counts is employed. The directory name is *Wilcox_rank_sum*



$${\color{brown}FDRTYPE}$$: 

Denotes the FDR (q-value) derivation method.

1. DiffLoop_BH_FDR_{DIFFFDRTHR}_LOG2FC_{LFCTHR}: Directory storing outputs when the conventional BH correction is used to determine FDR. Here {DIFFFDRTHR} is the FDR threshold for differential loops (see the parameter --DiffFDRThr) and {LFCTHR} is the corresponding log fold change threshold (see the parameter --LFCThr).

2. DiffLoop_IHW_FDR_{DIFFFDRTHR}_LOG2FC_{LFCTHR}: Directory storing outputs when the [IHW  correction](https://bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html) is used to determine FDR.

**Note:**

1. The **DiffLoops_ALL.bed** file within this directory structure stores the significant differential loops.

2. The WashU browser compatible files for these differential loops are stored in the same directory, by the names **DiffLoops_ALL_WashU.bed.gz** and **DiffLoops_ALL_WashU.bed.gz.tbi**.


Utility scripts
================

**1. utils/Creates_json_files.R**

Execution of the above mentioned differential analysis script creates lots of differential loop files, which need to be visualized in a genome browser (like [WashU epigenome browser](https://epigenomegateway.wustl.edu/)) for visualization and comparison between different settings.

User can check this script and edit the parameters (lines 10 - 90) according to the differential analysis scripts and execution paths.

The output file will be a .json file which can be easily loaded in WashU epigenome browser for visualization.

Data repository
================
A stable version of the DiffHiChIP code (as of October 2025) is available on Zenodo: https://zenodo.org/records/17410330

Differential loop results for all datasets and DiffHiChIP settings can be accessed via our web server (https://ay-lab-tools.lji.org/DiffHiChIP/) and through the same Zenodo repository.

Queries
=========

For any queries, please e-mail to the authors.

Sourya Bhattacharyya (sourya@lji.org)

Daniela Salgado Figueroa (dsfigueroa@lji.org)

Ferhat Ay (ferhatay@lji.org)

