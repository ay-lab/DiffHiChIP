DiffHiChIP
=============

Statistically significant differential loops between two conditions / categories for HiChIP (and other compatible 3C protocols like PLAC-seq, ChIA-PET, or Hi-C) data.

Authors (and Contacts)
=======================

Sourya Bhattacharyya (sourya@lji.org)

Daniela Salgado Figueroa (dsfigueroa@lji.org)


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

1. R (we have used R version 4.1.0)

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
	
	--AllLoopList=ALLLOOPLIST
		Comma or colon separated list of loop files of all categories and all replicates. Loop files without any FDR thresholding (i.e. having all loops) should be provided. For example, considering FitHiChIP, files $PREFIX$*.interactions.FitHiC.bed should be provided. Loop files for the first category are to be mentioned first, followed by those of second category. For example, suppose there are two replicates for each category. Then, the file list would be like Cat1_Repl1_File:Cat1_Repl2_File:Cat2_Repl1_File:Cat2_Repl2_File where, Cat1 denotes input category 1, and Repl1 means replicate 1, and so on. Loop files can be gzipped as well. All the loops in all files should have equal resolution (bin size). Mandatory parameter.

	--ChrSizeFile=CHRSIZEFILE
		File containing size of chromosomes for reference genome. Mandatory parameter.

	--MidList=MIDLIST
		Boolean values (0/1) where 1 means midpoints of each bin are specified, rather than the bin interval. If a single value is provided, all of the input files are assumed to have the same property. Otherwise, if different input files have different configuration, user needs to specify a Comma or colon separated list of boolean values. The list should have equal length as the --AllLoopList parameter, and order of input files should be similar as --AllLoopList parameter. Default: list of all zeroes (means bin interval is provided for all the files).

	--BinSize=BINSIZE
		Applicable if all the entries in --MidList option is 1, that is all input files contain only midpoints. In such a case, user needs to specify the bin size of input loops. Mandatory parameter in such a case. Default = 0

	--CCColList=CCCOLLIST
		Comma or colon separated list of integers depicting the column numbers of individual input files, which contain the raw contact counts. If all the input files have similar settings, user may just specify one number, which will be applied to all the files. Default = 7, same as FitHiChIP output format.

	--QColList=QCOLLIST
		Comma or colon separated list of integers depicting the column numbers of individual input files, which contain the q-value (FDR). If all the input files have similar settings, user may just specify one number, which will be applied to all the files. Default is the last column for all the input loops, same as FitHiChIP output format.

	--FDRThr=FDRTHR
		FDR (q-value) threshold used for determining HiC or HiChIP significant loops. Default = 0.01.

	--OutDir=OUTDIR
		Base Output directory. Mandatory parameter.

	--CategoryList=CATEGORYLIST
		Comma or colon separated list of strings, corresponding to the two main categories (whose replicates are present). Default: Category1, Category2.

	--ReplicaCount=REPLICACOUNT
		Comma or colon separated list of the count of replicates for individual categories. Default: 1,1 (means that we are considering one replicate per sample).

	--ReplicaLabels1=REPLICALABELS1
		Comma or colon separated list of the label of replicates for the first category. Default: R1,R2, etc.

	--ReplicaLabels2=REPLICALABELS2
		Comma or colon separated list of the label of replicates for the second category. Default: R1,R2, etc.

	--ChIPAlignFileList=CHIPALIGNFILELIST
		Applicable only for HiChIP data processing. Comma or colon separated list of ChIP-seq alignment files. Either BAM or bedgraph formatted file is supported. Default is NULL. User can 1) either provide two files, one for each category, 2) or provide ChIP seq alignment files one for each sample. If this parameter is empty, HiC like processing is performed.

	--CovThr=COVTHR
		Applicable only for HiChIP data processing. Threshold signifying the max allowed deviation of ChIP coverage between two categories, to consider those bins as 1D invariant (ND). Default (or if user specifies value <=0 or > 100) is 25, means 25% deviation is set as maximum. If user chooses 50, 50% maximum ChIP seq coverage deviation would be allowed.

	--DiffModel=DIFFMODEL
		Differential loop model. Possible values are: [0]: Use complete set of loops, and report FDR and IHW corrected FDR. [1]: Use distance stratification (equal occupancy binning, based on the bin size specified). [2]: Differential loops per chromosome, useful for very large number of input loops. Default = 0.

	--DSBinSize=DSBINSIZE
		Applicable if the parameter --DiffModel is 1. Specifies the bin size (in bytes) employed for stratification. Default = 10000 means 10 Kb is the equal occupancy bin size.

	--UseDESeq2=USEDESEQ2
		If specified 1, uses DESeq2 to find the differential loops, provided that the input data has replicates in both categories. Default = 0, means edgeR is used to find the differential loops.

	--EdgeRModel=EDGERMODEL
		The exact statistical model to be employed for EdgeR (if used). Values are: [0] - exact Test, [1] - Quasi-likelihood (QL) F-Test (glmQLFit + glmQLFTest), [2] - (glmQLFit + glmTreat), [3] - likelihood ratio tests (glmFit + glmLRT), [4] - (glmFit + glmTreat). Default = 0. *** Note: values 1 to 4 are applicable only when both input categories have multiple replicates.

	--PreFilt=PREFILT
		If 1, pre-filtered loops which are significant in at least one sample, are used as background. Otherwise (if 0) no such pre-filtering is done. Default = 0.

	--DiffFDRThr=DIFFFDRTHR
		FDR threshold for DESeq / EdgeR. Default is 0.05, means that loops with FDR < 0.05, and fold change >= log2(FoldChangeThr) would be considered as differential.

	--LFCThr=LFCTHR
		DESeq / EdgeR log2 fold change threshold. Default = 1.

	--bcv=BCV
		If EdgeR is used with single samples (replica count = 1 for any of the categories), this value is the square-root-dispersion. For datasets arising from well-controlled experiments are 0.4 for human data, 0.1 for data on genetically identical model organisms, or 0.01 for technical replicates. For details, see the edgeR manual. By default, the value is set as 0.4.

	--Overwrite=OVERWRITE
		Binary variable. If 1, overwrites existing results. Default = 0.

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

	1. If --DiffModel=0, no distance stratification is employed. The folder name is *DiffModel_0_AllLoops*. 

	2. Otherwise, for distance stratification (--DiffModel=1), the directory name is *DiffModel_1_DistStrat_{DSBINSIZE}* where *DSBINSIZE* is the bin size employed for stratification (see the parameter --DSBinSize)

$${\color{green}MODEL}$$: 

Specifies the differential analysis model.

1. If --UseDESeq2=1, DESeq2 is employed. The directory name is *DESeq2*

2. Otherwise, edgeR is employed for differential loop finding.

	2a. --EdgeRModel=0 - directory: *EdgeR_exactTest*

	2b. --EdgeRModel=1 - directory: *EdgeR_glmQLFTest*

	2c. --EdgeRModel=2 - directory: *EdgeR_glmTreatQL*

	2d. --EdgeRModel=3 - directory: *EdgeR_glmLRT*

	2e. --EdgeRModel=4 - directory: *EdgeR_glmTreat*

	**Note** The --EdgeRModel=0 setting is the default. --EdgeRModel 1 to 4 values are only applicable when both input categories have multiple replicates.

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


Queries
=========

For any queries, please e-mail to the authors.

Sourya Bhattacharyya (sourya@lji.org)

Daniela Salgado Figueroa (dsfigueroa@lji.org)

