#!/usr/bin/env Rscript

#===========================================================
## differential analysis on HiChIP and Hi-C loops
# Author: Sourya Bhattacharyya
# Vijay-Ay lab, LJI
#===========================================================

suppressMessages(library(GenomicRanges))
library(optparse)
library(DESeq2)
library(edgeR)
library(tools)
library(dplyr)
library(ggplot2)
library(data.table)
# Independent Hypothesis Weighting (Ignatiadis et al. 2016)
# P-value weighting
library(IHW)
library(ashr)
library(BiocParallel)
library(parallel)	# library for parallel processing

# disable converting scientific notations
options(scipen = 10)

# export data tables as data frames
options(datatable.fread.datatable=FALSE)

## default midlist parameter
## midpoints are not specified
DEFAULT_MIDLIST_ARG <- 0

## default CCColList parameter
## contact count column = 7
DEFAULT_CC_COL_ARG <- 7

## default QColList parameter
## FDR column = -1 (last column)
DEFAULT_FDR_COL_ARG <- -1

## default FDR threshold for FitHiChIP loops
DEFAULT_FDR_LOOP <- 0.01

## Applicable only for HiChIP data processing. 
## Threshold signifying the max allowed deviation of ChIP coverage between two categories, 
## to consider those bins as 1D invariant (ND).
DEFAULT_CHIP_COVERAGE_THR <- 0.25

## FDR threshold for DESeq / EdgeR
DEFAULT_DIFFLOOP_FDR_THR <- 0.05

## log fold change threshold for DESeq / EdgeR
DEFAULT_LFC_THR <- 1

## bcv threshold for EdgeR
DEFAULT_bcv_THR <- 0.4



