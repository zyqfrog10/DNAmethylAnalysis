
## ---------------------------
##
## Script name: BSseq_analysis.R
##
## Purpose of script: This script to run analyses on WGBS/RRBS/EPICseq with mulitple options packages
## 
## Note: BSseq (K.D. Hansen et al) will read in data, DSS/DMRseq (K. Korthauer et al) will smooth and call DMRs
##
## Author: Yingqian Zhan
## 
## Affiliation: Center for Epigenetics Research, MSKCC 
##
## Date Created: 2019-11-25
## Current update: 2019-12-07
## Copyright (c) Yingqian Zhan, 2019
## Email: yzhan@mskcc.org
##
## ---------------------------

source("BSseq_functions.R")
message("\n***BSseq_analysis.R\n\n")

## Parsing the command line ###
library(argparser)
if (packageVersion("argparser") < 0.3) {
  stop("argparser version (", packageVersion("argparser"), ") is out of date - 0.3 or later is required. Please open R and run install.packages('argparser') to update.")
}

message("\n")

args = commandArgs(trailingOnly=T)

p <- arg_parser("Run analysis on WGBS/RRBS/EPICseq", name="Rscript BSseq_analysis.R psa_condition.txt")
p <- add_argument(p, "--sampleSheet", help="A text file containing information of sample_name, dir, file and condition. One file per sample. Fileformat should be bismark coverage or genome wide cytosine.")
p <- add_argument(p, "--contrastFile", help="Comma separated text file with the samples to compare. Condition names must be the same as
                        condition in the sampleSheet file.")
p <- add_argument(p, "--BSseqData", help="A BSseq object in .rds format.")
p <- add_argument(p, "--outDir",help="Output directory. Default: The same direcotry of sampleSheet or bSseq file.")
p <- add_argument(p, "--analysis",help="Type of analyses: <explore|dmrseq|DSS>, one or multiple with comma separated.")
p <- add_argument(p, "--RNAseqList", help="A text file containing a list of file paths to RNAseq DESeq2 result table. Comparison must match with the contrastFile.")
p <- add_argument(p, "--genome",help="Type of genome: mm9,mm10,hg19, or hg38.")
opts <- parse_args(p, args)

## argparser read-in ###
if(packageVersion("argparser") < 0.4)
{
  names(opts) <- gsub("-", "_", names(opts))
}
message("Read in input ...\n")
sample_sheet = opts[["sampleSheet"]]
contrast = opts[["contrastFile"]]
BSseqData = opts[["BSseqData"]]
outDir = opts[["outDir"]]
analysis_list = strsplit(opts[["analysis"]], "\\,")[[1]]
RNAseq_sheet = opts[["RNAseqList"]]
genome = opts[["genome"]]

## Check the input data
#if(is.na(sample_sheet) & is.na(BSseqData)){
#  stop("Please provide either a sampleSheet or/and a BSseq class file.\n")
#}else if(!file.exists(sample_sheet) & !file.exists(BSseqData)){
#  stop("Please provide a sampleSheet or/and a BSseq class file.\n")
#}
#if(!file.exists(sample_sheet) & !file.exists(BSseqData) ){
#  stop("Please provide a sampleSheet or/and a BSseq class file.\n")
#}

if(is.na(sample_sheet) || rlang::is_empty(sample_sheet)){
  message("Sample sheet is not provided in this analysis.")
  sample_sheet <- as.character(sample_sheet)
  if(is.na(BSseqData) || rlang::is_empty(BSseqData)){
    stop("Please provide either a sampleSheet or/and a BSseq class file.\n")
  }else if(!file.exists(BSseqData)){
    stop("BSseqData cannot be accessed.")
  }else{
    message("BSseqData is: ", BSseqData, "\n")
  }
}else{
  if(file.exists(sample_sheet)){
    message("sample sheet is: ", sample_sheet, "\n")
  }else{
    message("Sample sheet cannot be accessed.")
  }

  if(is.na(BSseqData) || rlang::is_empty(BSseqData)){
    BSseqData <- as.character(BSseqData)
    message("BSseqData is not provided")
  }else if(!file.exists(BSseqData)){
    message("BSseqData cannot be accessed.","\n")
    if(!file.exists(sample_sheet)){
      stop("Both sample sheet and BSseqData cannot be accessed.")
    }
  }else{
    message("BSseqData is: ", BSseqData, "\n")
  }
}

if((is.na(outDir) || rlang::is_empty(outDir)) & file.exists(sample_sheet)){
  outDir <- dirname(sample_sheet)
  outDir <- paste(outDir,"Results",sep="/")
  message("Output directory: ",outDir)
}else if((is.na(outDir) || rlang::is_empty(outDir)) & file.exists(BSseqData)){
  outDir <- dirname(BSseqData)
  outDir <- paste(outDir,"Results",sep="/")
  message("Output directory: ",outDir)
}

if("dmrseq" %in% analysis_list & (is.na(contrast) || rlang::is_empty(contrast))){
  stop("DMR analysis requires a contrast file.")
}

if("DSS" %in% analysis_list & (is.na(contrast) || rlang::is_empty(contrast))){
  stop("DMR analysis requires a contrast file.")
}

if(is.na(contrast)){
  contrast <- as.character(contrast)
}

if(is.na(RNAseq_sheet)){
  RNAseq_sheet <- as.character(RNAseq_sheet)
}

if(!(genome %in% c("mm9","mm10","hg19","hg38"))){
  stop("genome must be one of: mm9, mm10, hg19, hg38")
}

if(!dir.exists(outDir)){
  message('Creating ', outDir, '\n')
  dir.create(outDir)
}
setwd(outDir)


## Set environment ###
message("Loading packages and dependencies...\n")
library(data.table)
library(BiocParallel)
library(bsseq)
library(aod)
library(DelayedMatrixStats)
library(ggplot2)
library(dendextend)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggcorrplot)
library(rafalib)
library(factoextra)
library(ggpubr)
library(httr)
set_config(config(ssl_verifypeer = 0L)) # solve the biomart connection problem that may happen sometimes
library(annotatr)
numCores <- 4
#register(MulticoreParam(workers = numCores), default = TRUE)

## maint code ###
if(file.exists(sample_sheet) & !file.exists(BSseqData)){
  min_cols <- c("sample_name","dir","file","condition")
  ### load in bismark files from sampleSheet
  sampleSheet <- fread(sample_sheet)
  if(sum(min_cols %in% colnames(sampleSheet)) != 4){
    stop("A sampleSheet contains at least four columns if BSseq file is not provided, <sample_name>,<dir>,<file>, and <condition>.\n")
  }else{
    sampleSheet$bismark.file <- file.path(sampleSheet$dir,sampleSheet$file)
    #temp_dir <- file.path(outDir,'tmp')
  #dir.create(temp_dir)
  #temp_file <- tempfile("BSseq",tmpdir=temp_dir)
    bismarkBSseq <- read.bismark(sampleSheet$bismark.file,
                            colData = DataFrame(row.names = sampleSheet$sample_name),
                            rmZeroCov = TRUE,
                            BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
                           # BACKEND = "HDF5Array",
                           # dir = temp_file,
                            nThread = numCores,
                            verbose = 2)
  ### add metadata for easy accessing the condition information
    bfile <- file.path(outDir, "bsseq.rds")
    if(ncol(pData(bismarkBSseq))==0){
      pData(bismarkBSseq) <- data.frame(condition = sampleSheet$condition,
                                    row.names = sampleSheet$sample_name,
                                    stringsAsFactors = FALSE)
      saveRDS(bismarkBSseq, bfile)}
  }
}else if(!file.exists(sample_sheet) & file.exists(BSseqData)){
  bismarkBSseq <- readRDS(BSseqData)
  message("Samples in BSseq object: ", paste(sampleNames(bismarkBSseq),collapse=","),"\n")
}else if(file.exists(sample_sheet) & file.exists(BSseqData)){
  bismarkBSseq <- readRDS(BSseqData)
  sampleSheet <- fread(sample_sheet)
  len <- length(sampleSheet$sample_name)
  min_cols <- c("sample_name","condition")
  if(sum(min_cols %in% colnames(sampleSheet)) != 2){
    stop("A sampleSheet contains at least two columns if BSseq file is provided, <sample_name> and <condition>.\n")
  }else{
    if(sum(sampleSheet$sample_name %in% sampleNames(bismarkBSseq)) == 0){
        stop("None of sample_name in sampleSheet matches with sample names in BSseq file.\n")
    }else{
      idx <- which(sampleSheet$sample_name %in% sampleNames(bismarkBSseq))
      message("Samples in the original BSseq object: ", paste(sampleNames(bismarkBSseq),collapse=","),"\n")
      message("Samples in the BSseq object adapted to sampleSheet: ", paste(sampleSheet$sample_name[idx],collapse=","),"\n")
      bismarkBSseq <- bismarkBSseq[,sampleSheet$sample_name[idx]]
      pData(bismarkBSseq)$condition <- sampleSheet$condition[idx]
    } 
  }
}
print(bismarkBSseq)
head(getCoverage(bismarkBSseq), n=4)
pData(bismarkBSseq)

### smooth data --- memory insufficient
#bismarkBSseq.fit <- BSmooth(BSseq = bismarkBSseq,
 #                           BPPARAM = MulticoreParam(workers = numCores,progressbar = TRUE),
  #                          verbose = TRUE)

### common variables
conditions <- unique(pData(bismarkBSseq)$condition)
sample_number <- length(sampleNames(bismarkBSseq))
fig_height <- ceiling(sample_number/2)

### Filter CpGs
seqn <- as.character(unique(seqnames(bismarkBSseq))) # remove unwanted random contigs
seqpreserve <- seqn[grep("random",seqn,invert=TRUE)]
seqpreserve <- seqpreserve[grep("alt",seqpreserve,invert=TRUE)]
seqpreserve <- seqpreserve[grep("chrUn",seqpreserve,invert=TRUE)]
seqpreserve <- seqpreserve[grep("fix",seqpreserve,invert=TRUE)]
seqdrop <- seqn[grep("random|alt|chrUn|fix",seqn)]
bismarkBSseq <- bismarkBSseq[which(seqnames(bismarkBSseq) %in% seqpreserve)]
bismarkBSseq <- dropSeqlevels(bismarkBSseq, seqdrop, pruning.mode="coarse")
### remove all CpGs that have coverage less than 10 in each sample
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov")>=10) >= sample_number)
#loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov")>=5) >= sample_number)
bismarkBSseq.filt <- bismarkBSseq[loci.idx,]
### discards the bases that have more than 99.9th percentile of coverage in each sample
BS.cov <- getCoverage(bismarkBSseq.filt,type="Cov")
upper <- quantile(BS.cov,probs =  0.999) 
keepLoci.exhi <- which(DelayedMatrixStats::rowSums2(BS.cov <= upper) >= sample_number)
bismarkBSseq.filt <- bismarkBSseq.filt[keepLoci.exhi,]

### analysis
#### exploratory plots
if("explore" %in% analysis_list){
  ### Downsizing to 50,000 CpGs by randomization
  set.seed(399)
  name <- "allSample"
  loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov")>=1) >= sample_number)
  sample_bs <- bismarkBSseq[loci.idx,]
  if(nrow(sample_bs) > 100000){
   sample_bs <- sample_bs[sample(1:nrow(sample_bs),50000),]}
  sample_bsCov <- getCoverage(sample_bs, type="Cov")
  sample_bsMethCov <- getCoverage(sample_bs, type="M")
  sample_bsBeta <- sample_bsMethCov/(sample_bsCov + 1) 
  pdf(paste("Hist_Cov_",name,"_raw.pdf",sep=""),width=10,height=10)
  print(covPlot(sample_bsCov,sample_bs,type="hist",var="Coverage"))
  dev.off()
  pdf(paste("Hist_MethCov_",name,"_raw.pdf",sep=""),width=10,height=10)
  print(covPlot(sample_bsMethCov,sample_bs,type="hist",var="Methyl Coverage"))
  dev.off()
  pdf(paste("Hist_beta_",name,"_raw.pdf",sep=""),width=10,height=10)
  print(covPlot(sample_bsBeta,sample_bs,type="hist",var="Methyl level"))
  dev.off()
  pdf(paste("Boxplot_Cov_",name,"_raw.pdf",sep=""),width=7,height=10)
  print(covPlot(sample_bsCov,sample_bs,type="boxplot",var="Coverage"))
  dev.off()
  pdf(paste("Boxplot_MethCov_",name,"_raw.pdf",sep=""),width=7,height=10)
  print(covPlot(sample_bsMethCov,sample_bs,type="boxplot",var="Methyl Coverage"))
  dev.off()
  pdf(paste("Boxplot_beta_",name,"_raw.pdf",sep=""),width=7,height=10)
  print(covPlot(sample_bsBeta,sample_bs,type="boxplot",var="Methyl level"))
  dev.off()
  rm(sample_bs)
 # bismarkBSseq.reduced <- bismarkBSseq.filt[sample(1:nrow(bismarkBSseq.filt),50000),]
  bismarkBSseq.reduced <- bismarkBSseq.filt
  if(nrow(bismarkBSseq.reduced) > 100000){
   bismarkBSseq.reduced <- bismarkBSseq.reduced[sample(1:nrow(bismarkBSseq.reduced),50000),]}
  ### methylation and coverage
  meth.cov <- getCoverage(bismarkBSseq.reduced, type="M")
  BS.cov <- getCoverage(bismarkBSseq.reduced, type="Cov")
  umeth.cov <- BS.cov-meth.cov
  M.value <- log2((meth.cov+1)/(BS.cov-meth.cov+1))
  beta.value <- meth.cov/(BS.cov + 1)
  ### hist and box plot
  pdf(paste("Hist_Cov_",name,"_filt.pdf",sep=""),width=10,height=10)
  print(covPlot(BS.cov,bismarkBSseq.reduced,type="hist",var="Coverage"))
  dev.off()
  pdf(paste("Hist_MethCov_",name,"_filt.pdf",sep=""),width=10,height=10)
  print(covPlot(meth.cov,bismarkBSseq.reduced,type="hist",var="Methyl Coverage"))
  dev.off()
  pdf(paste("Hist_beta_",name,"_filt.pdf",sep=""),width=10,height=10)
  print(covPlot(beta.value,bismarkBSseq.reduced,type="hist",var="Methyl level"))
  dev.off()
  pdf(paste("Boxplot_Cov_",name,"_filt.pdf",sep=""),width=7,height=10)
  print(covPlot(BS.cov,bismarkBSseq.reduced,type="boxplot",var="Coverage"))
  dev.off()
  pdf(paste("Boxplot_MethCov_",name,"_filt.pdf",sep=""),width=7,height=10)
  print(covPlot(meth.cov,bismarkBSseq.reduced,type="boxplot",var="Methyl Coverage"))
  dev.off()
  pdf(paste("Boxplot_beta_",name,"_filt.pdf",sep=""),width=7,height=10)
  print(covPlot(beta.value,bismarkBSseq.reduced,type="boxplot",var="Methyl level"))
  dev.off()

  ### correlation Heatmap
  name <- "methyl_overall"
  #png(paste("Heatmap_cor_excldiag_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
  pdf(paste("Heatmap_cor_excldiag_",name,".pdf",sep=""),width=10,height=10)
  print(correlationPlot(meth.cov))
  dev.off()
  #correlationPlot(umeth.cov,"umet_overall",15)
  name <- "M_overall"
  #png(paste("Heatmap_cor_excldiag_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
  pdf(paste("Heatmap_cor_excldiag_",name,".pdf",sep=""),width=10,height=10)
  print(correlationPlot(M.value))
  dev.off()

  name <- "beta_overall"
  #png(paste("Heatmap_cor_excldiag_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
  pdf(paste("Heatmap_cor_excldiag_",name,".pdf",sep=""),width=10,height=10)
  print(correlationPlot(beta.value))
  dev.off()
  ### dendrogram
  name <- "methyl_overall"
  #png(paste("Dendrogam_",name,"_mqc.png",sep=""),res=300,width=sample_number,height=6*1.5,units="cm",type="cairo")
  pdf(paste("Dendrogam_",name,".pdf",sep=""),width=sample_number/2,height=6)
  par(mar = c(12, 2, 2, 2))
  plot(dendPlot(meth.cov,bismarkBSseq.reduced,"complete"))
  dev.off()
  #dendPlot(umeth.cov,bismarkBSseq,"complete","umet_overall",sample_number,sample_number*1.5)
  name <- "M_wardD2_overall"
  #png(paste("Dendrogam_",name,"_mqc.png",sep=""),res=300,width=sample_number,height=6*1.5,units="cm",type="cairo")
  pdf(paste("Dendrogam_",name,".pdf",sep=""),width=sample_number/2,height=6)
  par(mar = c(12, 2, 2, 2))
  plot(dendPlot(M.value,bismarkBSseq.reduced,"ward.D2"))
  dev.off()

  name <- "beta_wardD2_overall"
  #png(paste("Dendrogam_",name,"_mqc.png",sep=""),res=300,width=sample_number,height=6*1.5,units="cm",type="cairo")
  pdf(paste("Dendrogam_",name,".pdf",sep=""),width=sample_number/2,height=6)
  par(mar = c(12, 2, 2, 2))
  plot(dendPlot(beta.value,bismarkBSseq.reduced,"ward.D2"))
  dev.off()
  ### PCA and Umap
  name <- "methyl_overall"
  #png(paste("PCA_",name,"_mqc.png",sep=""),res=300,width=30,height=15,units="cm",type="cairo")
  pdf(paste("PCA_",name,".pdf",sep=""),width=30/2,height=15/2.54)
  print(PCAPlot(meth.cov,bismarkBSseq.reduced))
  dev.off()

  if(sample_number > 20){
  pdf(paste("umap_",name,".pdf",sep=""),width=15/2,height=15/2.54)
  print(umapPlot(meth.cov,bismarkBSseq.reduced))
  dev.off()
  }

  pdf(paste("MDS_",name,".pdf",sep=""),width=15/2,height=15/2.54)
  print(mdsPlot(meth.cov,bismarkBSseq.reduced))
  dev.off()

  #PCAPlot(umeth.cov,bismarkBSseq,"umet_overall",15,7.5)
  name <- "M_overall"
  #png(paste("PCA_",name,"_mqc.png",sep=""),res=300,width=30,height=15,units="cm",type="cairo")
  pdf(paste("PCA_",name,".pdf",sep=""),width=30/2,height=15/2.54)
  print(PCAPlot(M.value,bismarkBSseq.reduced))
  dev.off()

  if(sample_number > 20){
  pdf(paste("umap_",name,".pdf",sep=""),width=15/2,height=15/2.54)
  print(umapPlot(M.value,bismarkBSseq.reduced))
  dev.off()
  }
  pdf(paste("MDS_",name,".pdf",sep=""),width=15/2,height=15/2.54)
  print(mdsPlot(M.value,bismarkBSseq.reduced))
  dev.off()

  name <- "beta_overall"
  #png(paste("PCA_",name,"_mqc.png",sep=""),res=300,width=30,height=15,units="cm",type="cairo")
  pdf(paste("PCA_",name,".pdf",sep=""),width=30/2,height=15/2.54)
  print(PCAPlot(beta.value,bismarkBSseq.reduced))
  dev.off()
  if(sample_number > 20){
  pdf(paste("umap_",name,".pdf",sep=""),width=15/2,height=15/2.54)
  print(umapPlot(beta.value,bismarkBSseq.reduced))
  dev.off()
  }
  pdf(paste("MDS_",name,".pdf",sep=""),width=15/2,height=15/2.54)
  print(mdsPlot(beta.value,bismarkBSseq.reduced))
  dev.off()  

  ### distribution of coverage and beta by group and sample
  #png("stat_Distribution_overall_mqc.png",res=300,width=15,height=7.5,units="cm",type="cairo")
  pdf(paste0("stat_Distribution_overall_",name,".pdf"),width=8,height=3.5)
  print(distPlot(bismarkBSseq.reduced))
  dev.off()
  ### clean up
  rm(BS.cov,meth.cov,umeth.cov,M.value,beta.value,bismarkBSseq.reduced)
}

### DMR analysis
if(file.exists(contrast)){
  comparisonsGrid <- t(as.matrix(read.csv(contrast, header=FALSE))) 
  check_cond <- comparisonsGrid %in% pData(bismarkBSseq)$condition
  if(sum(check_cond) < length(check_cond)){
    message("Conditions in BSseq object: ",paste(unique(pData(bismarkBSseq)$condition),collapse=","), "\n")
    message("Conditions in contrast File: ",paste(unique(c(comparisonsGrid[1,],comparisonsGrid[2,])),collapse=","),"\n")
    stop("Unknown condition(s) in contrast File.\n")
  }else if(length(unique(c(comparisonsGrid[1,],comparisonsGrid[2,]))) == 1){
    stop("Invalid comparison. Only one condition specified.")
  }
}else{
  comparisonsGrid <- NULL
}

if(!is.null(comparisonsGrid)){
  for (i in 1:ncol(comparisonsGrid)) {
    tryCatch({
    # Determine which conditions we're comparing in this iteration.
    compPair <- comparisonsGrid[, i]
    cond1 <- compPair[1]
    cond2 <- compPair[2]

    if(cond1 == cond2){
      stop("Invalid comparison on the same condition.\n")
    }else{
      cond1_idx <- which(pData(bismarkBSseq)$condition %in% cond1)
      cond2_idx <- which(pData(bismarkBSseq)$condition %in% cond2)
      samples.sub <- rownames(pData(bismarkBSseq))[c(cond1_idx,cond2_idx)]
      bismarkBSseq.sub <- bismarkBSseq.filt[,samples.sub]
      compPairName <- sprintf("%s_vs_%s", cond1, cond2)
      message(sprintf("Comparing conditions: %s vs. %s", cond1, cond2))

      # Ensure output subfolder that's specific to this pair of conditions.
      compPairOutfolder <- file.path(outDir, compPairName)
      R.utils::mkdirs(compPairOutfolder, mustWork=TRUE)
      setwd(compPairOutfolder)
      
      # dmrseq
      if("dmrseq" %in% analysis_list){
        message("Running dmrseq ... \n")
        name_sort <- sort(c(cond1,cond2))
        compPairName <- sprintf("%s_vs_%s", name_sort[2], name_sort[1])
        message(sprintf("Comparing conditions: %s vs. %s", name_sort[2], name_sort[1]))
        resfile <- file.path(compPairOutfolder, paste("dmrseq",genome,compPairName,"regions.rds",sep="_"))

        if(!file.exists(resfile)){
          library(dmrseq)
          run_dmrseq(bismarkBSseq.sub,compPairOutfolder,genome,compPairName)
        }
        if(file.exists(RNAseq_sheet)){
          regions <- readRDS(resfile)
          rnaseq_dir <- read.table(RNAseq_sheet,skip=i-1,nrows=1,stringsAsFactors=FALSE)
          rnaseq <- read.table(rnaseq_dir$V1)
          if(cond1 == name_sort[1]){
            rnaseq$log2FoldChange <- -rnaseq$log2FoldChange
          }
          methy_TSSann <- annotTSS(regions,genome,OutputType="df")
          methy_rna(rnaseq,methy_TSSann,compPairName,"dmrseq",genome,compPairOutfolder)
        } 
      }

      # DSS
      if("DSS" %in% analysis_list){
        message("Running DSS ... \n")
        compPairName <- sprintf("%s_vs_%s", cond1, cond2)
        message(sprintf("Comparing conditions: %s vs. %s", cond1, cond2))
        pval <- 0.1
        diffM <- 0.05
        dmls.file <- file.path(compPairOutfolder,paste("DSS",genome,compPairName,"DMCtest.rds",sep="_"))
        dmrs.file <- file.path(compPairOutfolder,paste0("DSS_DMR_",genome,"_",compPairName,"_delta0.1_p",pval,".rds"))        
        
        if(!file.exists(dmls.file) | !file.exists(dmrs.file)){
          library(DSS)
          run_DSS(bismarkBSseq.sub,compPairOutfolder,genome,condition1=cond1,condition2=cond2,compPairName,pval,diffM)
        }
        if(file.exists(RNAseq_sheet)){
          dmls <- readRDS(dmls.file)
          dmrs <- readRDS(dmrs.file)
          rnaseq_dir <- read.table(RNAseq_sheet,skip=i-1,nrows=1,stringsAsFactors=FALSE)
          rnaseq <- read.table(rnaseq_dir$V1)
          ### DMCs
          gr_DM <- with(dmls,GRanges(chr,IRanges(pos,pos),mu1=mu1,mu2=mu2,diffMethy=diff,pval=pval,fdr=fdr))
          dmls_TSSann <- annotTSS(gr_DM,genome,OutputType="df")
          methy_rna(rnaseq,dmls_TSSann,compPairName,"DSS",genome,compPairOutfolder)
          ### DMRs
          gr_DM <- with(dmrs,GRanges(chr,IRanges(start,end),nCG=nCG,meanMethy1=meanMethy1,meanMethy2=meanMethy2,diffMethy=diff.Methy,areaStat=areaStat))
          dmrs_TSSann <- annotTSS(gr_DM,genome,OutputType="df")
          methy_rna(rnaseq,dmrs_TSSann,compPairName,"DSS",genome,compPairOutfolder)
        }
      }
    }
   },error=function(e){cat("Error",conditionMessage(e), "\n")})
  }
}
  

### multivariate: PLSDA
if("PLSDA" %in% analysis_list){ 
  ## Apply PLS-DA 
  ### methylation and coverage
  meth.cov <- getCoverage(bismarkBSseq.filt, type="M")
  BS.cov <- getCoverage(bismarkBSseq.filt, type="Cov")
  #umeth.cov <- BS.cov-meth.cov
  M.value <- log2((meth.cov+1)/(BS.cov-meth.cov+1))
  beta.value <- meth.cov/(BS.cov + 1)
  
  library("mixOmics")
  pldaPlot(meth.cov,bismarkBSseq.filt,"methy")
 # pldaPlot(umeth.cov,bismarkBSseq.filt,"umet",sample_number)
  pldaPlot(M.value,bismarkBSseq.filt,"M")
  pldaPlot(beta.value,bismarkBSseq.filt,"beta")
}
