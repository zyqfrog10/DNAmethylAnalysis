
## ---------------------------
##
## Script name: BSseq_functions.R
##
## Purpose of script: This script contains functions used by BSseq_analysis.R
## 
## Author: Yingqian Zhan
## 
## Affiliation: Center for Epigenetics Research, MSKCC 
##
## Date Created: 2019-11-25
## Current update: 2020-08-24
## Copyright (c) Yingqian Zhan, 2019
## Email: yzhan@mskcc.org
##
## ---------------------------

## Functions ###
### correlation Plot
#correlationPlot <- function(matrix,name,size){
correlationPlot <- function(matrix){
    cormat <- round(cor(matrix,
                      use = "pairwise.complete.obs", 
                      method = "spearman"),2)

    diag(cormat) <- 0
    max <- max(cormat)
    diag(cormat) <- 1
    min <- min(cormat)

    ggplotHeatmap <- ggcorrplot(cormat, hc.order = TRUE,hc.method = "ward.D2", 
                              show.diag = FALSE,
                              outline.col = "white") + 
                    scale_fill_gradient2(limit = c(min,max), 
                                      high = "#FC4E07", low =  "#00AFBB", 
                                      mid = "#E7B800", midpoint = median(cormat))
   # options(bitmapType='cairo',device="png")
   # ggsave(paste("Heatmap_cor_excldiag_",name,"_mqc.png",sep=""),
   #       width=size,height=size,units="cm",
   #       plot=ggplotHeatmap,device="png")
   return(ggplotHeatmap)
   # rm(cormat)
}

### dendrogram
#dendPlot <- function(matrix,bs,hcmethod,name,wd,ht){
dendPlot <- function(matrix,bs,hcmethod){
    num.color <- pData(bs) %>% unique() %>% nrow()
    value.color <- get_palette(palette="npg",num.color)
    names(value.color) <- pData(bs)$condition %>% unique() %>% sort()
    c.order <- value.color[pData(bs)$condition]
    d <- hclust(dist(t(matrix),
                  method="euclidean"),method=hcmethod)
    d$labels <- rownames(pData(bs))
    dend <- as.dendrogram(d)
    #labels_colors(dend) <- as.fumeric(pData(bs)$condition)[order.dendrogram(dend)]
    labels_colors(dend) <- c.order[order.dendrogram(dend)]
    return(dend)
    #png(paste("Dendrogam_",name,"_mqc.png",sep=""),res=300,width=wd,height=ht,units="cm",type="cairo")
    #par(mar = c(10, 2, 2, 2))
    ##plot(dend)
    #print(plot(dend))
    #dev.off()
    #rm(d)
    #rm(dend)
}

### PCA 
#PCAPlot <- function(matrix,bs,name,wd,ht){
PCAPlot <- function(matrix,bs){
    pca.mat <- prcomp(t(matrix),center=T, scale.=F, retx=T, tol=NULL)
    num.color <- pData(bs) %>% unique() %>% nrow()
    #png(paste("PCA_scree_ind_",name,"_mqc.png",sep=""),res=300,width=wd,height=ht,units="cm",type="cairo")
    #print(fviz_eig(pca.mat))
    #dev.off()
    #png(paste("PCA_ind_",name,"_mqc.png",sep=""),res=300,width=wd,height=wd,units="cm",type="cairo")
    #print(fviz_pca_ind(pca.mat,
    p1 <- fviz_pca_ind(pca.mat, geom=c("point", "text"),
                  #habillage = pData(bs)$condition,
                  col.ind=pData(bs)$condition,
                  legend.title = "Groups",
                  repel = TRUE,
                  show.legend.text = F)+
                  scale_color_manual(values=get_palette(palette="npg",num.color))+
                  scale_shape_manual(values=rep(19,num.color))#+ ggsave("pca_test.pdf")
              #)
    #dev.off()
    #rm(pca.mat)
    p2 <- fviz_eig(pca.mat)
    p <- ggarrange(p1,p2,ncol=2,nrow=1)
    return(p)
}

### umap 
#umapPlot <- function(matrix,bs,name,wd,ht){
umapPlot <- function(matrix,bs){
    umap.obj <- umap::umap(t(matrix))
    umap.label <- rownames(t(matrix))
    umap.label <- pData(bs)[umap.label,"condition"]
    umap.obj <- data.frame(umap.obj$layout)
    umap.obj[,"Groups"] <- umap.label
    umap.obj[,"name"] <- rownames(umap.obj)
    num.color <- pData(bs) %>% unique() %>% nrow()
    p <- ggscatter(umap.obj, x="X1", y="X2",color = "Groups", 
                  label = "name", repel = TRUE,
                  show.legend.text = F)+
                 # ellipse = TRUE) +
            scale_color_manual(values=get_palette(palette="npg",num.color))+
            labs(title="UMAP",x="X1",y="X2")+
            theme(legend.position="right")+
            theme_minimal() #+ ggsave("umap_test.pdf")

    return(p)
    #png(paste("umap_",name,"_mqc.png",sep=""),res=300,width=wd,height=wd,units="cm",type="cairo")
    #print(p)
    #dev.off()
}

### MDS
#mdsPlot <- function(matrix,bs,name,wd,ht){
mdsPlot <- function(matrix,bs){
    df <- data.frame(t(matrix))
    df <- scale(df)
    df.dist <- dist(df)
    mdsout <- data.frame(cmdscale(df.dist))
    mdsout[,"Groups"] <- pData(bs)[rownames(mdsout),"condition"]
    mdsout[,"name"] <- rownames(mdsout)
    num.color <- pData(bs) %>% unique() %>% nrow()
    p <- ggscatter(mdsout, x="X1", y="X2",color = "Groups", 
                  label = "name", repel = TRUE,
                  show.legend.text = F) +
                  # ellipse = TRUE) +
            scale_color_manual(values=get_palette(palette="npg",num.color))+
            labs(title="MDS",x="Dim1",y="Dim2")+
            theme(legend.position="right")+
            theme_minimal()
    return(p)
    #png(paste("MDS_",name,"_mqc.png",sep=""),res=300,width=wd,height=wd,units="cm",type="cairo")
    #print(p)
    #dev.off()
}

### PLSDA plot 
#pldaPlot <- function(matrix,bs,name,wd){
pldaPlot <- function(matrix,bs,name){
    X <- t(matrix)
    Y <- as.character(pData(bs)$condition)
    result.plda <- plsda(X,Y)

    background <- background.predict(result.plda, comp.predicted=2,
                                dist = "max.dist")
   num.color <- pData(bs) %>% unique() %>% nrow() 
  # pdf(paste("PLSDA_ind_",name,".pdf",sep=""),width=10,height=10)
   # print(plotIndiv(result.plda, comp = 1:2, group = Y,
  # plotIndiv(result.plda, comp = 1:2, group = Y,
  #        ind.names = FALSE, title = "PLSDA: Maximum distance", ellipse = TRUE,
  #        legend = TRUE,  background = background) #)
  #  dev.off()

    pdf(paste("PLSDA_ind_",name,".pdf",sep=""),width=15/2,height=15/2.54)
   # print(plotIndiv(result.plda, comp = 1:2, group = Y,
    plotIndiv(result.plda,
              title = "PLSDA",legend=TRUE,
              col.per.group=get_palette(palette="npg",num.color),
              size.title=12,
              size.legend=9,
              size.lengend.title=10,
              legend.title="Groups")#)
    dev.off()
   #pdf(paste("PLSDA_var_",name,"_mqc.png",sep=""),width=15/2.54,height=15/2.54)
   # print(plotVar(result.plda, var.names=FALSE))
   #p2 <- plotVar(result.plda, var.names=FALSE)
   #dev.off()
   pdf(paste("PLSDA_loading1_",name,".pdf",sep=""),width=15/2.54,height=10)
   # print(
   plotLoadings(result.plda, comp = 1, title = 'Loadings on comp 1', contrib = 'max', 
             method = 'median', ndisplay = 10, size.name = 0.6)#)
   dev.off()

   pdf(paste("PLSDA_loading2_",name,".pdf",sep=""),width=15/2.54,height=10)
   # print(
   plotLoadings(result.plda, comp = 2, title = 'Loadings on comp 2', contrib = 'max', 
             method = 'median', ndisplay = 10, size.name = 0.6)#)
   dev.off()
   #p <- ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
  # return(p)
  }

### heatmap
#heatmapPlot <- function(matrix,clustRow,labelRow,clustCol=2, name,wd,ht){
heatmapPlot <- function(matrix,bs,clustRow,labelRow){
  # set the custom distance and clustering functions, per your example
  hclustfunc <- function(x) hclust(x, method="complete")
  distfunc <- function(x) dist(x, method="euclidean")

  # samples and conditions
  condition <- data.frame(pData(bs))
  clustCol <- length(unique(condition$condition))
  condition_group <- condition %>% mutate(group_no = as.integer(as.factor(condition)))

  # perform clustering on rows and columns
  cl.row <- hclustfunc(distfunc(matrix))
  cl.col <- hclustfunc(distfunc(t(matrix)))

  # extract cluster assignments; i.e. k=clustRow (rows) k=2 (columns)
  gr.row <- cutree(cl.row, clustRow)
  gr.col <- cutree(cl.col, clustCol)

  # assign color by condition
  #gr.col <- condition_group$group_no
  #names(gr.col) <- rownames(condition)  

  # set color
  col1 <- colorspace::qualitative_hcl(clustRow,palette = "dark3") # by clusters for rows

  ## by sample conditions for columns
  #col2 <- colorspace::qualitative_hcl(clustCol,palette = "warm") 
  col2 <- get_palette(palette="npg",clustCol)
  names(col2) <- pData(bs)$condition %>% unique() %>% sort()

  # plot the heatmap
  #if(labelRow!=NULL | labelRow!=""){
   # labelRow <- rownames(matrix)
  #}
  library(gplots)
  #png(paste("Heatmap_",name,"_mqc.png",sep=""),res=300,width=wd,height=ht,units="cm",type="cairo")
  #print(
  p <-  heatmap.2(matrix,scale="row",trace="none",labRow = labelRow,
			            hclustfun=hclustfunc, distfun=distfunc,
			            key=TRUE, symkey=FALSE,key.title="",key.ylab="",cexRow=1,cexCol=1,
			            margins=c(12,8),
                  srtCol=45,
			            ColSideColors = col2[pData(bs)$condition],
			            RowSideColors=col1[gr.row],#lhei = c(1.2,8),#lwid = c(0.5,4),
			            col = colorspace::diverging_hcl(11,"Blue-Red"),
			            denscol = "grey")#)
  return(p)
  #dev.off()
}

### distribution of coverage and beta by group and sample
#distPlot <- function(bs,name,wd,ht){ 
distPlot <- function(bs){
  num.color <- pData(bs) %>% unique() %>% nrow()
  grp.color <- get_palette(palette="npg",num.color)
  names(grp.color) <- pData(bs)$condition %>% unique() %>% sort()
  condition_group <- data.frame(pData(bs)) %>% mutate(group_no = as.integer(as.factor(condition)))
  cov.dist <- dmrseq::plotEmpiricalDistribution(bs, 
                           bySample = TRUE,
                            testCovariate = "condition",
                            type = "Cov") + theme(legend.position = "none") +
                            guides(linetype=FALSE) +
                            scale_color_manual(values=grp.color[pData(bs)$condition])+
                            scale_linetype_manual(values=rep(1,nrow(pData(bs))))
                           # scale_linetype_manual(values=condition_group$group_no)#+
                           # ylim(0,0.5)
                            

  beta.dist <- dmrseq::plotEmpiricalDistribution(bs, 
                            bySample = TRUE,
                            testCovariate = "condition",
                            adj = 3) +
                            guides(linetype=FALSE) +
                            scale_color_manual(values=grp.color[pData(bs)$condition])+
                            scale_linetype_manual(values=rep(1,nrow(pData(bs))))
                           # scale_linetype_manual(values=condition_group$group_no)#+
                          #  ylim(0,2)
                           # scale_linetype_manual(values = sample)
  #png(paste("stat_Distribution_",name,"_mqc.png",sep=""),res=300,width=wd,height=ht,units="in",type="cairo")
  #print(
  p <- ggarrange(cov.dist,beta.dist, 
            ncol = 2, nrow = 1,
            # widths = c(0.6, 1))
            widths = c(0.8, 1))#)
  return(p)
  #dev.off()
  #rm(cov.dist, beta.dist)
}

### Coverage histogram and boxplot
covPlot <- function(matrix,bs,type="hist",var="Coverage"){
  #BS.cov <- getCoverage(bs, type="Cov")
  groups <- pData(bs) %>% data.frame()
  num.color <- nrow(unique(groups))
  df_long <- matrix %>% data.frame() %>%
                melt(variable.name="sample_name",
                        value.name="Coverage",
                        id.vars = NULL)
  df_long$Groups <- groups[df_long$sample_name,1]
  
  if(type=="hist"){
    p <- df_long %>% na.omit() %>%
          ggplot(aes(x=Coverage+0.001,fill=Groups))+
           geom_histogram(color="black")+
           xlab(var)+
           scale_fill_manual(values=get_palette(palette="npg",num.color))+
          # scale_fill_brewer(palette="jco")+
           scale_x_log10()+
           theme_minimal() + facet_wrap(~sample_name) 
  }else{
    p <- df_long %>% na.omit() %>%
          ggplot(aes(x=sample_name,y=Coverage+0.001,fill=Groups))+
           geom_boxplot(color="black")+
           ylab(var)+
           scale_fill_manual(values=get_palette(palette="npg",num.color))+
          # scale_fill_brewer(palette="jco")+
           scale_y_log10()+
           xlab("")+coord_flip()+
           theme_minimal()
  }
  return(p)
}

### gene annotation and CpG annotation plot
annotat_dm <- function(grDM,genome,name,resultDir){
  annots <- c(paste(genome,"_cpgs",sep=""), 
             paste(genome,"_basicgenes",sep=""), 
             paste(genome,"_genes_intergenic",sep=""),
             paste(genome,"_genes_intronexonboundaries",sep=""))
  # annots = c('hg38_custom_TSS')     
  annotations <- build_annotations(genome = genome, annotations = annots)
      
  dm_annotated <- annotate_regions(regions = grDM,
                                  annotations = annotations,
                                  ignore.strand = TRUE,
                                  quiet = FALSE)
  #if(reDMann){return(dm_annotated)}
  df_dm_annotated <- data.frame(dm_annotated)
  #write.table(df_dm_annotated,"DSS_DMR_delta0.1_p0.05_geneAnn.txt",sep="\t",quote=FALSE,row.names=FALSE)
  csvfile <- file.path(resultDir, paste(name,"_geneAnn.csv",sep=""))
  write.csv(df_dm_annotated,csvfile,row.names=FALSE)
  resfile <- file.path(resultDir, paste(name,"_geneAnn.rds",sep=""))
  saveRDS(dm_annotated, file=resfile)
  return(dm_annotated)
}

annotatPlot <- function(grDM,genome,name,resultDir){
  resfile <- file.path(resultDir, paste(name,"_geneAnn.rds",sep=""))
  if (!file.exists(resfile)){
    dm_annotated <- annotat_dm(grDM,genome,name,resultDir)
  }else{
    dm_annotated <- readRDS(resfile)
  }
  annots_order <- c(paste(genome,"_genes_1to5kb",sep=""),
                    paste(genome,"_genes_promoters",sep=""),
                    paste(genome,"_genes_5UTRs",sep=""),
                    paste(genome,"_genes_exons",sep=""),
                    paste(genome,"_genes_intronexonboundaries",sep=""),
                    paste(genome,"_genes_introns",sep=""),
                    paste(genome,"_genes_3UTRs",sep=""),
                    paste(genome,"_genes_intergenic",sep=""))    
    
 # png(paste("stat_",name,"_gene_mqc.png",sep=""),res=300,width=10,height=10,units="cm",type="cairo")
 # print(
  p1 <- plot_annotation(annotated_regions = dm_annotated,
                          annotation_order = annots_order,
                          plot_title = '',
                          x_label = 'Known Gene Annotations',
                          y_label = 'Count')#)
  #dev.off()

  annots_order <- c(paste(genome,"_cpg_islands",sep=""),
                      paste(genome,"_cpg_shores",sep=""),
                      paste(genome,"_cpg_shelves",sep=""),
                      paste(genome,"_cpg_inter",sep=""))
    
  #png(paste("stat_",name,"_CpG_mqc.png",sep=""),res=300,width=10,height=10,units="cm",type="cairo")
  #print(
   p2 <- plot_annotation(annotated_regions = dm_annotated,
                          annotation_order = annots_order,
                          plot_title = '',
                          x_label = 'Known CpG Annotations',
                          y_label = 'Count')#)
  #dev.off()
  p <- ggarrange(p1,p2,ncol=2,nrow=1)
  return(p)
  #png(paste("stat_",name,"_scatter_mqc.png",sep=""),res=300,width=10,height=10,units="cm",type="cairo")
  #plot_numerical(annotated_regions = dm_annotated,
   #              x = 'meanMethy1',
    #             y = 'meanMethy2',
     #            facet = 'annot.type',
      #           facet_order = c(paste(genome,"_genes_1to5kb",sep=""),
       #                          paste(genome,"_genes_promoters",sep=""),
        #                         paste(genome,"_genes_5UTRs",sep=""),
         #                        paste(genome,"_genes_3UTRs",sep=""), 
          #                       paste(genome,"_genes_intergenic",sep=""), 
           #                      paste(genome,"_cpg_islands",sep=""), 
            #                     paste(genome,"_cpg_shores",sep="")),
             #   plot_title = 'Methylation comparison',
              #  x_label = 'condition1',
               # y_label = 'condition2')
  #dev.off()
}

### filter the loci by coverage (remove those with zero coverage in at least one condition).
filterOne <- function(bs){
    conditions <- unique(pData(bs)$condition)
    BS.cov <- getCoverage(bs, type="Cov")
    filter <- pmax( 1*(rowSums2(BS.cov[,pData(bs)$condition == conditions[1]]) == 0),
                  1*(rowSums2(BS.cov[,pData(bs)$condition == conditions[2]]) == 0))
    filter <- which(filter > 0)
    bs.filtered <- bs[-filter,]
    rm(BS.cov)
    return(bs.filtered)
}

### filter out .1% high coverage and low coverage less than 2 in each condition
filterTwo <- function(bs){
    conditions <- unique(pData(bs)$condition)
    BS.cov <- getCoverage(bs,type="Cov")
    upper <- quantile(BS.cov,probs =  0.999)
    n <- min(table(pData(bs)$condition)[1],table(pData(bs)$condition)[2])
    ##n <- ceiling(n/2)
    keepLoci.ex <- which(rowSums(BS.cov[, pData(bs)$condition == conditions[1]] >= 2) >= n &
                      rowSums(BS.cov[, pData(bs)$condition == conditions[2]] >= 2) >= n)
    keepLoci.exhi <- which(BS.cov <= upper)
    keepLoci <- intersect(keepLoci.ex,keepLoci.exhi)
    ##bismarkBSseq.filt <- bismarkBSseq[keepLoci.ex,]
    bs.filtered <- bs[keepLoci,]
    rm(BS.cov)
    return(bs.filtered)
}

# annotate DMR with custom TSS
annotTSS <- function(grDM, genome,uplim=1000,dwlim=1000,OutputType="df"){
    if(genome == "hg19"){
        require(TxDb.Hsapiens.UCSC.hg19.knownGene)
        require(org.Hs.eg.db)
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
        species <- "Hs"
    }else if(genome == "hg38"){
        require(TxDb.Hsapiens.UCSC.hg38.knownGene)
        require(org.Hs.eg.db)
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
        species <- "Hs"
    }else if(genome == "mm10"){
        require(TxDb.Mmusculus.UCSC.mm10.knownGene)
        require(org.Mm.eg.db)
        txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
        species <- "Mm"
    }else if(genome == "mm9"){
        require(TxDb.Mmusculus.UCSC.mm9.knownGene)
        require(org.Mm.eg.db)
        txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
        species <- "Mm"
    }else{
        stop("Currently only support hg19, hg38, mm10 and mm9.")
    }

    # build custom TSS region
    TSS_txdb <- promoters(genes(txdb), upstream = uplim, downstream = dwlim) # get TSS for genes not transcript

    # convert id to symbol 
    x = get(sprintf('org.%s.egSYMBOL', species))
    mapped_genes = mappedkeys(x)
    eg2symbol = as.data.frame(x[mapped_genes]) # gene_id and symbol table

    ## add symbol column matched to gene_id
    mcols(TSS_txdb)$symbol = eg2symbol[match(mcols(TSS_txdb)$gene_id,eg2symbol$gene_id),'symbol']

    TSS_txdb_df <- data.frame(TSS_txdb)
    TSS_txdb_df$width <- paste(TSS_txdb_df$seqnames,TSS_txdb_df$start,sep=":")
    TSS_txdb_df$placeholder5 <- '.'
    TSS_txdb_df <- TSS_txdb_df[,c(1,2,3,4,8,5,6,7)]

    filename <- paste("txdb",genome,"TSS.bed",sep="_")
    write.table(TSS_txdb_df,filename,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
    TSS_path <- paste(getwd(),filename,sep="/")
    annotatr::read_annotations(con=TSS_path,genome=genome,name='TSS',format='bed',
                    extraCols=c(gene_id='character',symbol='character')) # build custom TSS annotation

    name <- paste(genome,"custom_TSS",sep="_")
    annotations <- annotatr::annotatr_cache$get(name)

    dm_annotated <- annotatr::annotate_regions(regions = grDM,
                                    annotations = annotations,
                                    ignore.strand = TRUE,
                                    quiet = FALSE)

    if(file.exists(TSS_path)){unlink(TSS_path)}
    if(OutputType=="df" || OutputType=="dataframe"){
      df_dm_annotated <- data.frame(dm_annotated)
      return(df_dm_annotated)
    }else if(OutputType=="gr" || OutputType=="GRanges"){
      return(dm_annotated)
    }else{
      stop("Please specify the output data type: dataframe (df) or GRanges (gr)")
    }   
}

### run dmrseq and plot
run_dmrseq <- function(bs,resultDir,genomeName,compName,ncores=4) {
  sample_number <- length(sampleNames(bs))
  #### run DMRseq
  register(MulticoreParam(ncores*2), default = TRUE) 
  set.seed(399)
  resfile <- file.path(resultDir, paste("dmrseq",genomeName,compName,"regions.rds",sep="_"))
  if (!file.exists(resfile)){
    if (is.unsorted(bs)) {
      bs <- sort(bs)
    }
    regions <- dmrseq(bs, 
                      testCovariate = "condition",
                      bpSpan = 500,
                      maxGap = 500, maxPerms = 10,
                      verbose = TRUE)
    #regions$meanDiff <- meanDiff(bismarkBSseq.filt, dmrs=regions, testCovariate="condition")
    saveRDS(regions, file=resfile)
  }else{
    regions <- readRDS(resfile)
  }

  ##### Plot top 3 regions
  annoTrack <- getAnnot(genomeName)
  for(index in 1:3){
    #png(paste("dmrseq_DMRregions",genomeName,compName,"region",index,"mqc.png",sep="_"),res=300,width=15,height=7.5,units="cm",type="cairo")
    pdf(paste("dmrseq_DMRregions_",genomeName,"_",compName,"_region_",index,".pdf",sep=""),width=15/2.54,height=7.5/2.54)
    tryCatch({
	    plotDMRs(bs, regions=regions[index,], testCovariate="condition",
				annoTrack=annoTrack)
    },error=function(e){cat("Error",conditionMessage(e), "\n")})
  dev.off()}
  ### parameters needs to be optimized for calling blocks
  #blockfile <- file.path(outDir, "blocks.rds")
  #if (!file.exists(blockfile)){
  #  blocks <- dmrseq(bismarkBSseq.filt, 
   #                   testCovariate = "condition",
    #                  block = TRUE,
     #                 minInSpan = 500,
      #                bpSpan = 5e4,
       #               maxGapSmooth = 1e6,
        #              maxGap = 5e3,
         #             verbose = TRUE)
    #saveRDS(blocks, file=blockfile)
  #}else{
   # blocks <- readRDS(blockfile)
  #}

  ## Plots for DMRs
  if(is.null(regions$mean.methy)){
     DM.beta <- getCoverage(bs,regions=regions,what="perRegionTotal",type="M")/getCoverage(bs,regions=regions,what="perRegionTotal",type="Cov")
     regions$mean.methy <- rowMeans(DM.beta,na.rm=TRUE) 
     rm(DM.beta)
  }
  if(is.null(regions$log2FC) | is.null(regions$diffMethy)){
    name_sort <- sort(unique(pData(bs)$condition))
    grp1.samples <- rownames(pData(bs))[pData(bs)$condition == name_sort[1]]
    grp2.samples <- rownames(pData(bs))[pData(bs)$condition == name_sort[2]]
    DM.beta.grp1 <- getCoverage(bs[,grp1.samples],regions=regions,what="perRegionTotal",type="M")/getCoverage(bs[,grp1.samples],regions=regions,what="perRegionTotal",type="Cov")
    DM.beta.grp2 <- getCoverage(bs[,grp2.samples],regions=regions,what="perRegionTotal",type="M")/getCoverage(bs[,grp2.samples],regions=regions,what="perRegionTotal",type="Cov")
    DM.beta.mean.grp1 <- rowMeans(DM.beta.grp1,na.rm=TRUE)
    DM.beta.mean.grp2 <- rowMeans(DM.beta.grp2,na.rm=TRUE)
    regions$log2FC <- log2((DM.beta.mean.grp2+1e-9)/(DM.beta.mean.grp1+1e-9))
    regions$diffMethy <- DM.beta.mean.grp2-DM.beta.mean.grp1
    rm(DM.beta.grp1,DM.beta.grp2)
    saveRDS(regions, file=resfile)
  }
  df_regions <- data.frame(regions)
  df_regions$significance <- "NO"
  if(dim(df_regions[which(df_regions$qval <= 0.05 & abs(df_regions$log2FC) >= log2(1.5)),])[1] > 100){
    df_regions[which(df_regions$qval <= 0.05 & abs(df_regions$log2FC) >= log2(1.5)),]$significance <- "YES"
  }else{
    df_regions[which(df_regions$pval <= 0.05 & abs(df_regions$log2FC) >= log2(1.5)),]$significance <- "YES"
  }
  
  ### Volcano plot
  ggplot(df_regions, aes(x = beta, y = -log10(qval),colour=significance))+
    geom_point(size=1) +
    scale_colour_manual(values=c("black", "red")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"), 
      panel.border=element_blank(),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),  
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14),
      plot.title = element_text(size=16)) +
    ggtitle("Volcano plot (dmrseq)")+
    xlab(paste("beta coefficient", compName,sep="\n")) +
    ylab("-log10(adjusted p-value)") #+ 
    ggsave(paste("dmrseqDMR_betaCoeff_",genomeName,"_",compName,".pdf",sep=""),
          width=10,height=10,units="cm")
 # png(paste("dmrseqDMR_betaCoeff",genomeName,compName,"_mqc.png",sep="_"),res=300,width=10,height=10,units="cm",type="cairo")
 # print(p)
 # dev.off()

  ggplot(df_regions, aes(x = log2FC, y = -log10(qval),colour=significance))+
    geom_point(size=1) +
    scale_colour_manual(values=c("black", "red")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"), 
      panel.border=element_blank(),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),  
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14),
      plot.title = element_text(size=16)) +
    ggtitle("Volcano plot (dmrseq)")+
    xlab(paste("log2 fold change", compName,sep="\n")) +
    ylab("-log10(adjusted p-value)") #+ 
    ggsave(paste("dmrseqDMR_volcano_",genomeName,"_",compName,".pdf",sep=""),
          width=10,height=10,units="cm")
  #png(paste("dmrseqDMR_volcano",genomeName,compName,"_mqc.png",sep="_"),res=300,width=10,height=10,units="cm",type="cairo")
  #print(p1)
  #dev.off()
  
  ### MA plot
  ggplot(df_regions, aes(x = mean.methy, y = log2FC,colour=significance))+
    geom_point(size=1) +
    geom_hline(yintercept=0,colour="dark red",alpha=0.2)+
    geom_hline(yintercept=log2(1.5),colour="dark blue",alpha=0.4)+
    geom_hline(yintercept=-log2(1.5),colour="dark blue",alpha=0.4)+
    scale_colour_manual(values=c("black", "red")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"), 
      panel.border=element_blank(),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),  
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14),
      plot.title = element_text(size=16)) +
    ggtitle("MA plot (dmrseq)")+
    xlab(paste("Mean methylation level", compName,sep="\n")) +
    ylab("log2 fold change") #+
    ggsave(paste("dmrseqDMR_MA_",genomeName,"_",compName,".pdf",sep=""),
          width=10,height=10,units="cm")
  #png(paste("dmrseqDMR_MA",genomeName,compName,"_mqc.png",sep="_"),res=300,width=10,height=10,units="cm",type="cairo")
  #print(p2)
  #dev.off()

  if(length(regions[which(mcols(regions)$qval <= 0.05)]) > 100){
    regions_sig <- regions[which(mcols(regions)$qval <= 0.05)]
  }else{
    regions_sig <- regions[which(mcols(regions)$pval <= 0.05)]
  }
  DM.cov <- getCoverage(bs,regions=regions_sig,what="perRegionTotal",type="Cov")
  DM.meth <- getCoverage(bs,regions=regions_sig,what="perRegionTotal",type="M")
  zero.idx <- which(rowSums(DM.meth)==0)
  if(length(zero.idx) > 0){
    DM.meth <- DM.meth[-zero.idx,]
    DM.cov <- DM.cov[-zero.idx,]
  }
 # DM.umeth <- DM.cov-DM.meth
  DM.M <- na.omit(log2((DM.meth+1)/(DM.cov-DM.meth+1)))
  DM.beta <- na.omit(DM.meth/(DM.cov + 1))
  
  ### stats
  bismarkBSseq.DM <- subsetByOverlaps(bs,regions_sig)
  name <- paste("dmrseqDMC",genomeName,compName,sep="_")
  #png(paste("stat_Distribution_",name,"_mqc.png",sep=""),res=300,width=8,height=4,units="in",type="cairo")
  pdf(paste("stat_Distribution_",name,".pdf",sep=""),width=8,height=3.5)
  print(distPlot(bismarkBSseq.DM))
  dev.off()
  
  #png(paste("stat_",name,"_gene_mqc.png",sep=""),res=300,width=8,height=4,units="in",type="cairo")
  pdf(paste("stat_",name,"_gene.pdf",sep=""),width=8,height=4)
  print(annotatPlot(granges(bismarkBSseq.DM),genomeName,name,resultDir))
  dev.off()
  
  name <- paste("dmrseqDMR",genomeName,compName,"q0.05",sep="_")
 # png(paste("stat_",name,"_gene_mqc.png",sep=""),res=300,width=8,height=4,units="in",type="cairo")
  pdf(paste("stat_",name,"_gene.pdf",sep=""),width=8,height=4)
  print(annotatPlot(regions_sig,genomeName,name,resultDir))
  dev.off()
  ### Plot top 3 regions
  regions_sub <- data.frame(regions[1:3])
  regions_sub <- regions_sub[,c(1,2,3,4,6,9)]
  colnames(regions_sub) <- c("chr","start","end","length","nCG","areaStat")
  for(index in 1:3){
    #png(paste("dmrseq_split","DMRregions",genomeName,compName,"region",index,"mqc.png",sep="_"),res=300,width=7.5,height=sample_number,units="in",type="cairo")
    pdf(paste("dmrseq_split","_DMRregions_",genomeName,"_",compName,"_region_",index,".pdf",sep=""),width=7.5,height=sample_number)
    tryCatch({
	    DSS::showOneDMR(data.frame(regions_sub[index,]), bs)
    },error=function(e){cat("Error",conditionMessage(e), "\n")})
  dev.off()}
  ### PCA
  name <- paste("methyl_dmrseqDMR",genomeName,compName,"q0.05",sep="_")
  #png(paste("PCA_",name,"_mqc.png",sep=""),res=300,width=30,height=15,units="cm",type="cairo")
  pdf(paste("PCA_",name,".pdf",sep=""),width=30/2.54,height=15/2.54)
  print(PCAPlot(DM.meth,bs))
  dev.off()
  #PCAPlot(DM.umeth,bs,paste("umet_dmrseqDMR_q0.05",genomeName,compName,sep="_"),15,7.5)
  name <- paste("M_dmrseqDMR",genomeName,compName,"q0.05",sep="_")
  #png(paste("PCA_",name,"_mqc.png",sep=""),res=300,width=30,height=15,units="cm",type="cairo")
  pdf(paste("PCA_",name,".pdf",sep=""),width=30/2.54,height=15/2.54)
  print(PCAPlot(DM.M,bs))
  dev.off()
  name <- paste("beta_dmrseqDMR",genomeName,compName,"q0.05",sep="_")
  #png(paste("PCA_",name,"_mqc.png",sep=""),res=300,width=30,height=15,units="cm",type="cairo")
  pdf(paste("PCA_",name,".pdf",sep=""),width=30/2.54,height=15/2.54)
  print(PCAPlot(DM.beta,bs))
  dev.off()

  ### Dendrogram
  name <- paste("methyl_dmrseqDMR",genomeName,compName,"q0.05",sep="_")
  #png(paste("Dendrogam_",name,"_mqc.png",sep=""),res=300,width=sample_number,height=6*1.5,units="cm",type="cairo")
  pdf(paste("Dendrogam_",name,".pdf",sep=""),width=sample_number/2.54,height=6*1.5/2.54)
  plot(dendPlot(DM.meth,bs,"complete"))
  dev.off()
 # dendPlot(DM.umeth,bs,"complete",paste("umet_dmrseqDMR_q0.05",genomeName,compName,sep="_"),sample_number,6*1.5)
  name <- paste("M_wardD2_dmrseqDMR",genomeName,compName,"q0.05",sep="_")
  #png(paste("Dendrogam_",name,"_mqc.png",sep=""),res=300,width=sample_number,height=6*1.5,units="cm",type="cairo")
  pdf(paste("Dendrogam_",name,".pdf",sep=""),width=sample_number/2.54,height=6*1.5/2.54)
  plot(dendPlot(DM.M,bs,"ward.D2"))
  dev.off()
  name <- paste("beta_wardD2_dmrseqDMR",genomeName,compName,"q0.05",sep="_")
  #png(paste("Dendrogam_",name,"_mqc.png",sep=""),res=300,width=sample_number,height=6*1.5,units="cm",type="cairo")
  pdf(paste("Dendrogam_",name,".pdf",sep=""),width=sample_number/2.54,height=6*1.5/2.54)
  plot(dendPlot(DM.beta,bs,"ward.D2"))
  dev.off()
  ## Heatmap for DMR
  if(length(regions_sig) <= 10000){
    name <- paste("beta_dmrseqDMR",genomeName,compName,"q0.05",sep="_")
    #png(paste("Heatmap_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
    pdf(paste("Heatmap_",name,".pdf",sep=""),width=5.5,height=5.5)
    print(heatmapPlot(DM.beta,bs,5,""))
    dev.off()
    name <- paste("M_dmrseqDMR",genomeName,compName,"q0.05",sep="_")
    #png(paste("Heatmap_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
    pdf(paste("Heatmap_",name,".pdf",sep=""),width=5.5,height=5.5)
    print(heatmapPlot(DM.M,bs,5,""))
    dev.off()
  }
  if(length(regions_sig) > 10000){
    DM.cov <- getCoverage(bs,regions=regions_sig[1:10000],what="perRegionTotal",type="Cov")
    DM.meth <- getCoverage(bs,regions=regions_sig[1:10000],what="perRegionTotal",type="M")
    zero.idx <- which(rowSums(DM.meth)==0)
    if(length(zero.idx > 0)){
      DM.meth <- DM.meth[-zero.idx,]
      DM.cov <- DM.cov[-zero.idx,]
    }
    #DM.umeth <- DM.cov-DM.meth
    DM.M <- na.omit(log2((DM.meth+1)/(DM.cov-DM.meth+1)))
    DM.beta <- na.omit(DM.meth/(DM.cov + 1))
    
    name <- paste("beta_dmrseqDMR",genomeName,compName,"top10000",sep="_")
    #png(paste("Heatmap_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
    pdf(paste("Heatmap_",name,".pdf",sep=""),width=5.5,height=5.5)
    print(heatmapPlot(DM.beta,bs,5,""))
    dev.off()

    name <- paste("M_dmrseqDMR",genomeName,compName,"top10000",sep="_")
    #png(paste("Heatmap_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
    pdf(paste("Heatmap_",name,".pdf",sep=""),width=5.5,height=5.5)
    print(heatmapPlot(DM.M,bs,5,""))
    dev.off()
  }
  ## Heatmap Plot for DMRs with TSS annotation 
  resfile <- file.path(resultDir, paste("dmrseq",genomeName,compName,"regions_TSSann.rds",sep="_"))
  if (!file.exists(resfile)){
    regions_TSSann <- annotTSS(regions,genomeName,OutputType="gr")
    saveRDS(regions_TSSann, file=resfile)
  }else{
    regions_TSSann <- readRDS(resfile)
  }
  
  DM.cov <- getCoverage(bs,regions=regions_TSSann,what="perRegionTotal",type="Cov")
  DM.meth <- getCoverage(bs,regions=regions_TSSann,what="perRegionTotal",type="M")
  
  zero.idx <- which(rowSums(DM.meth)==0)
  if(length(zero.idx) > 0){
      DM.meth <- DM.meth[-zero.idx,]
      DM.cov <- DM.cov[-zero.idx,]
  }
  #DM.umeth <- DM.cov-DM.meth
  DM.M <- na.omit(log2((DM.meth+1)/(DM.cov-DM.meth+1)))
  DM.beta <- na.omit(DM.meth/(DM.cov + 1))

  ## Heatmap for top 50 DMRs with TSS annotation
  if(length(regions_TSSann) <= 50){
    name <- paste("beta_dmrseqDMR",genomeName,compName,"TSSann",sep="_")
    #png(paste("Heatmap_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
    pdf(paste("Heatmap_",name,".pdf",sep=""),width=5.5,height=5.5)
    print(heatmapPlot(DM.beta,bs,5,regions_TSSann$annot$symbol))
    dev.off()

    name <- paste("M_dmrseqDMR",genomeName,compName,"TSSann",sep="_")
    #png(paste("Heatmap_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
    pdf(paste("Heatmap_",name,".pdf",sep=""),width=5.5,height=5.5)
    print(heatmapPlot(DM.M,bs,5,regions_TSSann$annot$symbol))
    dev.off()
  }else{
    tot_num <- length(regions_TSSann)
    name <- paste("beta_dmrseqDMR",genomeName,compName,"TSSann_50",tot_num,sep="_")
    #png(paste("Heatmap_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
    pdf(paste("Heatmap_",name,".pdf",sep=""),width=5.5,height=5.5)
    print(heatmapPlot(DM.beta[1:50,],bs,5,regions_TSSann$annot$symbol[1:50]))
    dev.off()

    name <- paste("M_dmrseqDMR",genomeName,compName,"TSSann_50",tot_num,sep="_")
    #png(paste("Heatmap_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
    pdf(paste("Heatmap_",name,".pdf",sep=""),width=5.5,height=5.5)
    print(heatmapPlot(DM.M[1:50,],bs,5,regions_TSSann$annot$symbol[1:50]))
    dev.off()
  }
}

### run DSS and plot
run_DSS <- function(bs,resultDir,genomeName,condition1=NULL,condition2=NULL,compName,pval,diffM,ncores=4){
  sample_number <- length(sampleNames(bs))
  register(MulticoreParam(ncores*2), default = TRUE) 
  conds <- unique(pData(bs)$condition)
  if(is.null(condition1) | is.null(condition2)){
    condition1 <- conds[1]
    condition2 <- conds[2]
  }
  set.seed(399)
  ## DML test
  c1_index <- which(pData(bs)$condition %in% condition1)
  c1 <- rownames(pData(bs))[c1_index]
  c2_index <- which(pData(bs)$condition %in% condition2)
  c2 <- rownames(pData(bs))[c2_index]
  dmlTest.file <- file.path(resultDir,paste("DSS",genomeName,compName,"DMCtest.rds",sep="_"))
  if(!file.exists(dmlTest.file)){
    dmlTest.sm <- DMLtest(bs,group1=c1,group2=c2,
                          smoothing=TRUE,BPPARAM = bpparam())
    saveRDS(dmlTest.sm, file=dmlTest.file)
  }else{
    dmlTest.sm <- readRDS(dmlTest.file)
  }

  ### plot volcano and MA plot on a subset of CpG sites
  if(nrow(dmlTest.sm)>100000){
    dmlTest.sub <- dmlTest.sm[sample(1:nrow(dmlTest.sm),50000),]
  }else{dmlTest.sub <- dmlTest.sm}
  dmlTest.sub$significance <- "NO"
  dmlTest.sub[which(dmlTest.sub$fdr <= pval & abs(log2(dmlTest.sub$mu1/(dmlTest.sub$mu2+0.001))) >= log2(1.5)),]$significance <- "YES"
  #dmlTest.sub[which(dmlTest.sub$fdr <= 0.05 & abs(dmlTest.sub$diff) >= 0.1),]$significance <- "YES"
  #### Volcano plot DMC
  h.scale <- max(abs(log2(dmlTest.sub$mu1/(dmlTest.sub$mu2+0.001))))
  h.scale <- round(h.scale,digits=1)
  ggplot(dmlTest.sub, aes(x = log2((mu1+1e-9)/(mu2+1e-9)), y = -log10(fdr),colour=significance))+
    geom_point(size=0.5) + xlim(-h.scale-0.1,h.scale+0.1)+ 
    scale_colour_manual(values=c("black", "red")) +
    theme_bw()+
    theme(panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"), 
      panel.border=element_blank(),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),  
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14),
      plot.title = element_text(size=12)) +
    ggtitle("Volcano plot (DSS DMC)")+
    xlab(paste("log2 fold change", compName,sep="\n")) #+
    ggsave(paste("DSSDMC_volcano_",genomeName,"_",compName,".pdf",sep=""),
        width=10,height=10,units="cm")
  #png(paste("DSSDMC_volcano",genomeName,compName,"_mqc.png",sep="_"),res=300,width=10,height=10,units="cm",type="cairo")
  #print(p1)
  #dev.off()
  #### MA plot DMC
  ngrp1 <- table(pData(bs)$condition)[condition1]
  ngrp2 <- table(pData(bs)$condition)[condition2]
  ggplot(dmlTest.sub, aes(x = (mu1*ngrp1+mu2*ngrp2)/(ngrp1+ngrp2), y = log2((mu1+1e-9)/(mu2+1e-9)),colour=significance))+
    geom_point(size=0.5) + ylim(-h.scale-0.1,h.scale+0.1)+
    geom_hline(yintercept=0,colour="dark red",alpha=0.2)+
    geom_hline(yintercept=log2(1.5),colour="dark blue",alpha=0.4)+
    geom_hline(yintercept=-log2(1.5),colour="dark blue",alpha=0.4)+
    scale_colour_manual(values=c("black", "red")) +
    theme_bw()+
    theme(panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"), 
      panel.border=element_blank(),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),  
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14),
      plot.title = element_text(size=12)) +
    ggtitle("MA plot (DSS DMC)")+
    xlab(paste("Mean methylation level", compName,sep="\n"))+
    ylab(paste("log2 fold change")) #+
    ggsave(paste("DSSDMC_MA_",genomeName,"_",compName,".pdf",sep=""),
          width=10,height=10,units="cm")
  #png(paste("DSSDMC_MA",genomeName,compName,"_mqc.png",sep="_"),res=300,width=10,height=10,units="cm",type="cairo")
  #print(p2)
  #dev.off()

  ### volcano and MA plot with diff methylation
  dmlTest.sub$significance <- "NO"
  dmlTest.sub[which(dmlTest.sub$fdr <= pval & abs(dmlTest.sub$diff) >= diffM),]$significance <- "YES"
  
  #### Volcano plot DMC
  h.scale <- max(abs(dmlTest.sub$diff))
  ggplot(dmlTest.sub, aes(x = diff, y = -log10(fdr),colour=significance))+
    geom_point(size=0.5) + xlim(-h.scale-0.1,h.scale+0.1)+
    scale_colour_manual(values=c("black", "red")) +
    theme_bw()+
    theme(panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"), 
      panel.border=element_blank(),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),  
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14),
      plot.title = element_text(size=12)) +
    ggtitle("Volcano plot (DSS DMC)")+
    xlab(paste("diff Methy", compName,sep="\n")) #+
    ggsave(paste("DSSDMC_volcano_",genomeName,"_",compName,"_diffmethyX.pdf",sep=""),
        width=10,height=10,units="cm")

ggplot(dmlTest.sub, aes(x = (mu1*ngrp1+mu2*ngrp2)/(ngrp1+ngrp2), y = diff,colour=significance))+
    geom_point(size=0.5) + ylim(-h.scale-0.1,h.scale+0.1)+
    geom_hline(yintercept=0,colour="dark red",alpha=0.2)+
    geom_hline(yintercept=0.1,colour="dark blue",alpha=0.4)+
    geom_hline(yintercept=-0.1,colour="dark blue",alpha=0.4)+
    scale_colour_manual(values=c("black", "red")) +
    theme_bw()+
    theme(panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"), 
      panel.border=element_blank(),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),  
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14),
      plot.title = element_text(size=12)) +
    ggtitle("MA plot (DSS DMC)")+
    xlab(paste("Mean methylation level", compName,sep="\n"))+
    ylab(paste("diff Methy")) #+
    ggsave(paste("DSSDMC_MA_",genomeName,"_",compName,"_diffmethyY.pdf",sep=""),
          width=10,height=10,units="cm")

  ## call DML
  dmls.file <- file.path(resultDir,paste0("DSS_DMC_",genomeName,"_",compName,"_delta",diffM,"_p",pval,".rds"))
  if(!file.exists(dmls.file)){
    dmls <- callDML(dmlTest.sm, delta=diffM, p.threshold=pval)
    saveRDS(dmls, file=dmls.file)
  }else{
    dmls <- readRDS(dmls.file)
  }
  ## call DMR
  tryCatch({
  dmrs.file <- file.path(resultDir,paste0("DSS_DMR_",genomeName,"_",compName,"delta",diffM,"_p",pval,".rds"))
  if(!file.exists(dmrs.file)){
    dmrs <- callDMR(dmlTest.sm, delta=diffM,p.threshold=pval)
    saveRDS(dmrs,file=dmrs.file)
  }else{
    dmrs <- readRDS(dmrs.file)
  }

  ### Volcano plot DMR
  #ggplot(dmrs, aes(x = log2((meanMethy1+1e-9)/(meanMethy2+1e-9)), y = areaStat))+
  h.scale <- max(abs(dmrs$diff.Methy))
  ggplot(dmrs, aes(x = diff.Methy, y = areaStat))+
    geom_point(size=1) + xlim(-h.scale-0.1,h.scale+0.1)+
    theme_bw()+
    theme(panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"), 
      panel.border=element_blank(),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),  
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14),
      plot.title = element_text(size=16)) +
    ggtitle("areaStat_vs_diff.Methy (DSS DMR)")+
    xlab(paste("Differential methylation level", compName,sep="\n")) #+
    ggsave(paste("DSSDMR_areaStat_vs_diff.Methy_",genomeName,"_",compName,".pdf",sep=""),
          width=10,height=10,units="cm")
  #png(paste("DSSDMR_areaStat_vs_diff.Methy",genomeName,compName,"_mqc.png",sep="_"),res=300,width=10,height=10,units="cm",type="cairo")
  #print(p1)
 # dev.off()
  ### MA plot DMR
  ext <- max(abs(log2((dmrs$meanMethy1+1e-9)/(dmrs$meanMethy2+1e-9))))
  ext <- round(ext,digits=1)
  ggplot(dmrs, aes(x = (meanMethy1*ngrp1+meanMethy2*ngrp2)/(ngrp1+ngrp2), y = log2((meanMethy1+1e-9)/(meanMethy2+1e-9))))+
    geom_point(size=1) + ylim(-ext-0.1,ext+0.1)+
    geom_hline(yintercept=0,colour="dark red",alpha=0.2)+
    geom_hline(yintercept=log2(1.5),colour="dark blue",alpha=0.4)+
    geom_hline(yintercept=-log2(1.5),colour="dark blue",alpha=0.4)+
    theme_bw() +
    theme(panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"), 
      panel.border=element_blank(),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),  
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14),
      plot.title = element_text(size=12)) +
    ggtitle("MA plot (DSS DMR)")+
    xlab(paste("Mean methylation level", compName,sep="\n")) +
    ylab("log2 fold change") #+
    ggsave(paste("DSSDMR_MA_",genomeName,"_",compName,".pdf",sep=""),
          width=10,height=10,units="cm")
  #png(paste("DSSDMR_MA",genomeName,compName,"_mqc.png",sep="_"),res=300,width=10,height=10,units="cm",type="cairo")
 # print(p2)
  #dev.off()

  ggplot(dmrs, aes(x = (meanMethy1*ngrp1+meanMethy2*ngrp2)/(ngrp1+ngrp2), y = diff.Methy))+
    geom_point(size=1) + ylim(-h.scale-0.1,h.scale+0.1)+
    geom_hline(yintercept=0,colour="dark red",alpha=0.2)+
   # geom_hline(yintercept=log2(1.5),colour="dark blue",alpha=0.4)+
   # geom_hline(yintercept=-log2(1.5),colour="dark blue",alpha=0.4)+
    theme_bw() +
    theme(panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"), 
      panel.border=element_blank(),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),  
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14),
      plot.title = element_text(size=12)) +
    ggtitle("MA plot (DSS DMR)")+
    xlab(paste("Mean methylation level", compName,sep="\n")) +
    ylab("diff Methy") #+
    ggsave(paste("DSSDMR_MA_",genomeName,"_",compName,"_diffmethyY.pdf",sep=""),
          width=10,height=10,units="cm")

  ### Plot top 3 regions
  for(index in 1:3){
    #png(paste("DSS_DMRregions",genomeName,compName,"region",index,"mqc.png",sep="_"),res=300,width=7.5,height=sample_number,units="in",type="cairo")
    pdf(paste("DSS_DMRregions_",genomeName,"_",compName,"_region_",index,".pdf",sep=""),width=7.5,height=sample_number)
    tryCatch({
	    DSS::showOneDMR(dmrs[index,], bs)
    },error=function(e){cat("Error",conditionMessage(e), "\n")})
  dev.off()}
  },error=function(e){cat("Error",conditionMessage(e), "\n")})
  ### Plots for DMC and DMR
  for(dm in c("DMC","DMR")){
    if(dm=="DMC"){
      #gr_DM <- with(dmls,GRanges(chr,IRanges(pos,pos)))
      gr_DM <- makeGRangesFromDataFrame(dmls,keep.extra.columns=TRUE,start.field="pos",end.field="pos")
      bismarkBSseq.DM <- subsetByOverlaps(bs,gr_DM)
      DM.cov <- getCoverage(bismarkBSseq.DM,type="Cov")
      DM.meth <- getCoverage(bismarkBSseq.DM,type="M")
      row_num <- nrow(dmls)

      name <- paste("DSSDMC",genomeName,compName,sep="_")
      #png(paste("stat_Distribution_",name,"_mqc.png",sep=""),res=300,width=8,height=4,units="in",type="cairo")
      pdf(paste("stat_Distribution_",name,".pdf",sep=""),width=8,height=3.5)
      print(distPlot(bismarkBSseq.DM))
      dev.off()
      name <- paste("DSSDMC",genomeName,compName,sep="_")
      #png(paste("stat_",name,"_gene_mqc.png",sep=""),res=300,width=8,height=4,units="in",type="cairo")
      pdf(paste("stat_",name,"_gene.pdf",sep=""),width=8,height=4)
      print(annotatPlot(gr_DM,genomeName,name,resultDir))
      dev.off()
      rm(bismarkBSseq.DM)
    }else{
      dmrs$abs.areaStat <- abs(dmrs$areaStat)
      dmrs_sorted <- dmrs[order(-dmrs$abs.areaStat),]
      gr_DM <- with(dmrs_sorted,GRanges(chr,IRanges(start,end),L=nCG,stat=areaStat,diff.Methy=diff.Methy))
      DM.cov <- getCoverage(bs,regions=gr_DM,what="perRegionTotal",type="Cov")
      DM.meth <- getCoverage(bs,regions=gr_DM,what="perRegionTotal",type="M")
      row_num <- nrow(dmrs)

      name <- paste("DSSDMR",genomeName,compName,sep="_")
      #png(paste("stat_",name,"_gene_mqc.png",sep=""),res=300,width=8,height=4,units="in",type="cairo")
      pdf(paste("stat_",name,"_gene.pdf",sep=""),width=8,height=4)
      print(annotatPlot(gr_DM,genomeName,name,resultDir))
      dev.off()
      annoTrack <- dmrseq::getAnnot(genomeName)
      for(index in 1:3){
        #png(paste("DSS_DMRregions",genomeName,compName,"AnnRegion",index,"mqc.png",sep="_"),res=300,width=15,height=7.5,units="cm",type="cairo")
        pdf(paste("DSS_DMRregions_",genomeName,"_",compName,"_AnnRegion_",index,".pdf",sep=""),width=15/2.54,height=7.5/2.54)
        tryCatch({
	       dmrseq::plotDMRs(bs, regions=gr_DM[index,],qval=FALSE,
                          testCovariate="condition",annoTrack=annoTrack)
        },error=function(e){cat("Error",conditionMessage(e), "\n")})
      dev.off()}
      rm(dmrs_sorted)
    }
    
    zero.idx <- which(rowSums(DM.meth)==0)
    if(length(zero.idx > 0)){
      DM.meth <- DM.meth[-zero.idx,]
      DM.cov <- DM.cov[-zero.idx,]
    }
    #DM.umeth <- DM.cov-DM.meth
    DM.M <- na.omit(log2((DM.meth+1)/(DM.cov-DM.meth+1)))
    DM.beta <- na.omit(DM.meth/(DM.cov + 1))

    ### PCA
    name <- paste("methyl_DSS",dm,genomeName,compName,sep="_")
    #png(paste("PCA_",name,"_mqc.png",sep=""),res=300,width=30,height=15,units="cm",type="cairo")
    pdf(paste("PCA_",name,".pdf",sep=""),width=11.8,height=5.9)
    print(PCAPlot(DM.meth,bs))
    dev.off()
    #PCAPlot(DM.umeth,bs,paste("umet_DSS",dm,sep=""),15,7.5)
    name <- paste("M_DSS",dm,genomeName,compName,sep="_")
    #png(paste("PCA_",name,"_mqc.png",sep=""),res=300,width=30,height=15,units="cm",type="cairo")
    pdf(paste("PCA_",name,".pdf",sep=""),width=11.8,height=5.9)
    print(PCAPlot(DM.M,bs))
    dev.off()

    name <- paste("beta_DSS",dm,genomeName,compName,sep="_")
    #png(paste("PCA_",name,"_mqc.png",sep=""),res=300,width=30,height=15,units="cm",type="cairo")
    pdf(paste("PCA_",name,".pdf",sep=""),width=11.8,height=5.9)
    print(PCAPlot(DM.beta,bs))
    dev.off()
    ### Dendrogram
    name <- paste("methyl_DSS",dm,genomeName,compName,sep="_")
    #png(paste("Dendrogam_",name,"_mqc.png",sep=""),res=300,width=sample_number,height=6*1.5,units="cm",type="cairo")
    pdf(paste("Dendrogam_",name,".pdf",sep=""),width=sample_number/2.54,height=6*1.5/2.54)
    plot(dendPlot(DM.meth,bs,"complete"))
    dev.off()
    #dendPlot(DM.umeth,bs,"complete",paste("umet_DSS",dm,sep=""),sample_number,6*1.5)
    name <- paste("M_wardD2_DSS",dm,genomeName,compName,sep="_")
   # png(paste("Dendrogam_",name,"_mqc.png",sep=""),res=300,width=sample_number,height=6*1.5,units="cm",type="cairo")
    pdf(paste("Dendrogam_",name,".pdf",sep=""),width=sample_number/2.54,height=6*1.5/2.54)
    plot(dendPlot(DM.M,bs,"ward.D2"))
    dev.off()
    
    name <- paste("beta_wardD2_DSS",dm,genomeName,compName,sep="_")
    #png(paste("Dendrogam_",name,"_mqc.png",sep=""),res=300,width=sample_number,height=6*1.5,units="cm",type="cairo")
    pdf(paste("Dendrogam_",name,".pdf",sep=""),width=sample_number/2.54,height=6*1.5/2.54)
    plot(dendPlot(DM.beta,bs,"ward.D2"))
    dev.off()
    ## Heatmap for DMC/DMR
    if(row_num <= 10000){
      name <- paste("beta_DSS",dm,genomeName,compName,sep="_")
      #png(paste("Heatmap_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
      pdf(paste("Heatmap_",name,".pdf",sep=""),width=5.5,height=5.5)
      heatmapPlot(DM.beta,bs,5,"")
      dev.off()

      name <- paste("M_DSS",dm,genomeName,compName,sep="_")
      #png(paste("Heatmap_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
      pdf(paste("Heatmap_",name,".pdf",sep=""),width=5.5,height=5.5)
      heatmapPlot(DM.M,bs,5,"")
      dev.off()
    }
    ## Plot for top 1000 DMC/DMR 
    if(row_num > 10000){
      if(dm=="DMC"){
        gr_DM <- with(dmls[1:10000,],GRanges(chr,IRanges(pos,pos)))
        bismarkBSseq.DM <- subsetByOverlaps(bs,gr_DM)
        DM.cov <- getCoverage(bismarkBSseq.DM,type="Cov")
        DM.meth <- getCoverage(bismarkBSseq.DM,type="M")
      }else{
        dmrs$abs.areaStat <- abs(dmrs$areaStat)
        dmrs_sorted <- dmrs[order(-dmrs$abs.areaStat),]
        gr_DM <- with(dmrs_sorted[1:10000,],GRanges(chr,IRanges(start,end)))
        DM.cov <- getCoverage(bs,regions=gr_DM,what="perRegionTotal",type="Cov")
        DM.meth <- getCoverage(bs,regions=gr_DM,what="perRegionTotal",type="M")
      }
    
      zero.idx <- which(rowSums(DM.meth)==0)
      if(length(zero.idx > 0)){
        DM.meth <- DM.meth[-zero.idx,]
        DM.cov <- DM.cov[-zero.idx,]
      }
     # DM.umeth <- DM.cov-DM.meth
      DM.M <- na.omit(log2((DM.meth+1)/(DM.cov-DM.meth+1)))
      DM.beta <- na.omit(DM.meth/(DM.cov + 1))

      ### Heatmap for DMC/DMR
      name <- paste("beta_DSS",dm,genomeName,compName,"top10000",sep="_")
      #png(paste("Heatmap_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
      pdf(paste("Heatmap_",name,".pdf",sep=""),width=5.5,height=5.5)
      heatmapPlot(DM.beta,bs,5,"")
      dev.off()

      name <- paste("M_DSS",dm,genomeName,compName,"top10000",sep="_")
     # png(paste("Heatmap_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
      pdf(paste("Heatmap_",name,".pdf",sep=""),width=5.5,height=5.5)
      heatmapPlot(DM.M,bs,5,"")
      dev.off()
    }
  
  ## Heatmap Plot for DMC/DMRs with TSS annotation 
  resfile <- file.path(resultDir, paste("DSS_DMR",genomeName,compName,dm,"TSSann.rds",sep="_"))
  if (!file.exists(resfile)){
    regions_TSSann <- annotTSS(gr_DM,genomeName,OutputType="gr")
    saveRDS(regions_TSSann, file=resfile)
  }else{
    regions_TSSann <- readRDS(resfile)
  }
  
  DM.cov <- getCoverage(bs,regions=regions_TSSann,what="perRegionTotal",type="Cov")
  DM.meth <- getCoverage(bs,regions=regions_TSSann,what="perRegionTotal",type="M")
  
  zero.idx <- which(rowSums(DM.meth)==0)
  if(length(zero.idx) > 0){
      DM.meth <- DM.meth[-zero.idx,]
      DM.cov <- DM.cov[-zero.idx,]
  }
  #DM.umeth <- DM.cov-DM.meth
  DM.M <- na.omit(log2((DM.meth+1)/(DM.cov-DM.meth+1)))
  DM.beta <- na.omit(DM.meth/(DM.cov + 1))

  ## Heatmap for top 50 DMRs with TSS annotation
  if(length(regions_TSSann) <= 50){
    name <- paste("beta_DSS",dm,genomeName,compName,"TSSann",sep="_")
    #png(paste("Heatmap_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
    pdf(paste("Heatmap_",name,".pdf",sep=""),width=5.5,height=5.5)
    print(heatmapPlot(DM.beta,bs,5,regions_TSSann$annot$symbol))
    dev.off()

    name <- paste("M_DSS",dm,genomeName,compName,"TSSann",sep="_")
    #png(paste("Heatmap_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
    pdf(paste("Heatmap_",name,".pdf",sep=""),width=5.5,height=5.5)
    print(heatmapPlot(DM.M,bs,5,regions_TSSann$annot$symbol))
    dev.off()
  }else{
    tot_num <- length(regions_TSSann)
    name <- paste("beta_DSS",dm,genomeName,compName,"TSSann_50",tot_num,sep="_")
    #png(paste("Heatmap_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
    pdf(paste("Heatmap_",name,".pdf",sep=""),width=5.5,height=5.5)
    print(heatmapPlot(DM.beta[1:50,],bs,5,regions_TSSann$annot$symbol[1:50]))
    dev.off()

    name <- paste("M_DSS",dm,genomeName,compName,"TSSann_50",tot_num,sep="_")
    #png(paste("Heatmap_",name,"_mqc.png",sep=""),res=300,width=15,height=15,units="cm",type="cairo")
    pdf(paste("Heatmap_",name,".pdf",sep=""),width=5.5,height=5.5)
    print(heatmapPlot(DM.M[1:50,],bs,5,regions_TSSann$annot$symbol[1:50]))
    dev.off()
  }

  }
}  

# methylation and RNAseq integration
methy_rna <- function(RNAdeseq2,methy_geneAnn,compName,method,genome,outDir){
	tablename <- paste(method,genome,compName,"DMR","TSSAnn","RNAexp.csv",sep="_")
    tablepath <- file.path(outDir,tablename)
    if(!file.exists(tablepath)){
        # read in RNAseq DESeq2 result table with gene name in the first column
        # RNAseq comparison should match with methylation comparison
        rna_seq <- RNAdeseq2
        rna_seq$annot.symbol <- rownames(rna_seq)
        
        # merge methylation annotation and RNAseq annotation and output the table
        methy_rna <- merge(methy_geneAnn,rna_seq,by.x="annot.symbol",by.y="annot.symbol")
        
        # save the table
        if(method=="dmrseq"){
          if(!is.null(methy_rna$qval)){
            methy_rna <- methy_rna[order(methy_rna$qval*methy_rna$padj),]
            write.csv(methy_rna,tablepath,row.names=FALSE)
          }else{stop("Unrecognized annotation table from dmrseq.")}
        }else if(method=="DSS"){ 
            if(!is.null(methy_rna$fdr)){
              tablename <- paste(method,compName,"DMC","TSSAnn","RNAexp.rds",sep="_")
              tablepath <- file.path(outDir,tablename)
              methy_rna <- methy_rna[order(methy_rna$fdr*methy_rna$padj),]
              saveRDS(methy_rna,file.path(outDir,tablename))
            }else if(!is.null(methy_rna$areaStat)){ 
                methy_rna <- methy_rna[order(-abs(methy_rna$areaStat),methy_rna$padj),]
                write.csv(methy_rna,tablepath,row.names=FALSE)
            }else{stop("Unrecognized annotation table from DSS.")}
        }else{
            stop("Method can only be either dmrseq or DSS.")
        }       
    }else{
        methy_rna <- read.csv(tablepath)
    }
	
	# quadrant plot
   # figname <- paste("Quadrant",method,genome,compName,"TSS_DMR","RNAexp_mqc.png",sep="_")
   figname <- paste("Quadrant",method,genome,compName,"TSS_DMR","RNAexp.pdf",sep="_")
   figpath <- file.path(outDir,figname)
	library(ggrepel)
    methy_rna$significance <- "NO"

    if(method=="dmrseq"){
      dim_n <- dim(methy_rna[which(methy_rna$qval <= 0.05 & methy_rna$padj <= 0.05 & abs(methy_rna$log2FoldChange) >= log2(1.5) & abs(methy_rna$diffMethy) >= 0.1),])[1]
      if( dim_n > 100){
        methy_rna[which(methy_rna$qval <= 0.05 & methy_rna$padj <= 0.05 & abs(methy_rna$log2FoldChange) >= log2(1.5) & abs(methy_rna$diffMethy) >= 0.1),]$significance <- "YES"
      }else{
        methy_rna[which(methy_rna$pval <= 0.05 & methy_rna$padj <= 0.05 & abs(methy_rna$log2FoldChange) >= log2(1.5) & abs(methy_rna$diffMethy) >= 0.1),]$significance <- "YES"
        dim_n <- dim(methy_rna[which(methy_rna$pval <= 0.05 & methy_rna$padj <= 0.05 & abs(methy_rna$log2FoldChange) >= log2(1.5) & abs(methy_rna$diffMethy) >= 0.1),])[1]
      }
      if(dim_n >=1){
        if(dim_n >= 20){
            n <- 20
          }else{n <- dim_n}
        # methy_rna[which(methy_rna$qval <= 0.05 & methy_rna$padj <= 0.05 & abs(methy_rna$log2FoldChange) >= log2(1.5) & abs(methy_rna$diffMethy) >= 0.1),]$significance <- "YES"
        ggplotMethyRna <- ggplot(methy_rna, aes(x=log2FoldChange, y=diffMethy,color = significance)) +
            geom_point(shape = 16, size = 1, alpha = .4) +
            scale_colour_manual(values=c("#0091ff", "#f0650e")) +
            geom_hline(yintercept=0, linetype="dashed", color = "grey")+
            geom_vline(xintercept=0, linetype="dashed", color = "grey")+
            geom_text_repel(data=methy_rna[which(methy_rna$significance == "YES"),][1:n,],
                aes(log2FoldChange,diffMethy,label=annot.symbol))+
            theme_minimal()+
            ggtitle(paste(compName,"(dmrseq)",sep=" "))+
            xlab(paste("log2FC, RNAseq", compName,sep="\n")) +
            ylab("Diff Methy at TSS") #+
            #ggsave("test_methy_rna_quadrantPlot.pdf")
      }else{
        # methy_rna[which(methy_rna$qval <= 0.05 & methy_rna$padj <= 0.05 & abs(methy_rna$log2FoldChange) >= log2(1.5) & abs(methy_rna$diffMethy) >= 0.1),]$significance <- "YES"
        ggplotMethyRna <- ggplot(methy_rna, aes(x=log2FoldChange, y=diffMethy,color = significance)) +
            geom_point(shape = 16, size = 1, alpha = .4) +
            scale_colour_manual(values=c("#0091ff", "#f0650e")) +
            geom_hline(yintercept=0, linetype="dashed", color = "grey")+
            geom_vline(xintercept=0, linetype="dashed", color = "grey")+
            theme_minimal()+
            ggtitle(paste(compName,"(dmrseq)",sep=" "))+
            xlab(paste("log2FC, RNAseq", compName,sep="\n")) +
            ylab("Diff Methy at TSS") #+
            #ggsave("test_methy_rna_quadrantPlot.pdf")
      } 
    }else if(method=="DSS"){ 
      if(!is.null(methy_rna$fdr)){
        #figname <- paste(method,compName,"TSS_DMC","RNAexp_mqc.png",sep="_")
        figname <- paste(method,compName,"TSS_DMC","RNAexp.pdf",sep="_")
        figpath <- file.path(outDir,figname)
        dim_n <- dim(methy_rna[which(methy_rna$fdr <= 0.1 & methy_rna$padj <= 0.1 & abs(methy_rna$log2FoldChange) >= log2(1.5) & abs(methy_rna$diffMethy) >= 0.1),])[1]
        if(dim_n >= 1){
          if(dim_n >= 20){
            n <- 20
          }else{n <- dim_n}
          methy_rna[which(methy_rna$fdr <= 0.1 & methy_rna$padj <= 0.1 & abs(methy_rna$log2FoldChange) >= log2(1.5) & abs(methy_rna$diffMethy) >= 0.1),]$significance <- "YES"
          ggplotMethyRna <- ggplot(methy_rna, aes(x=log2FoldChange, y=diffMethy,color = significance)) +
              geom_point(shape = 16, size = 1, alpha = .4) +
              scale_colour_manual(values=c("#0091ff", "#f0650e")) +
              geom_hline(yintercept=0, linetype="dashed", color = "grey")+
              geom_vline(xintercept=0, linetype="dashed", color = "grey")+
              geom_text_repel(data=methy_rna[which(methy_rna$significance == "YES"),][1:n,],
                  aes(log2FoldChange,diffMethy,label=annot.symbol))+
              theme_minimal()+
              ggtitle(paste(compName,"(DSS DMC)",sep=" "))+
              xlab(paste("log2FC, RNAseq", compName,sep="\n")) +
              ylab("Diff Methy at TSS")
          }else{
            ggplotMethyRna <- ggplot(methy_rna, aes(x=log2FoldChange, y=diffMethy,color = significance)) +
              geom_point(shape = 16, size = 1, alpha = .4) +
              scale_colour_manual(values=c("#0091ff", "#f0650e")) +
              geom_hline(yintercept=0, linetype="dashed", color = "grey")+
              geom_vline(xintercept=0, linetype="dashed", color = "grey")+
              theme_minimal()+
              ggtitle(paste(compName,"(DSS DMC)",sep=" "))+
              xlab(paste("log2FC, RNAseq", compName,sep="\n")) +
              ylab("Diff Methy at TSS")
          } #+
            #ggsave("test_DSS_DMC_methy_rna_quadrantPlot.pdf")
      }else if(!is.null(methy_rna$areaStat)){ 
        dim_n <- dim(methy_rna[which(abs(methy_rna$diffMethy)>0.1 & abs(methy_rna$log2FoldChange) > log2(1.5) & methy_rna$padj <= 0.05),])[1]
        if(dim_n >=1 ){
          if(dim_n >= 20){
              n <- 20
            }else{n <- dim_n}
          methy_rna[which(abs(methy_rna$diffMethy)>0.1 & abs(methy_rna$log2FoldChange) > log2(1.5) & methy_rna$padj <= 0.05),]$significance <- "YES"
          ggplotMethyRna <- ggplot(methy_rna, aes(x=log2FoldChange, y=diffMethy,color = significance)) +
              geom_point(shape = 16, size = 1, alpha = .4) +
              scale_colour_manual(values=c("#0091ff", "#f0650e")) +
              geom_hline(yintercept=0, linetype="dashed", color = "grey")+
              geom_vline(xintercept=0, linetype="dashed", color = "grey")+
              #geom_text_repel(data=subset(methy_rna, abs.areaStat > 570),
              geom_text_repel(data=methy_rna[which(methy_rna$significance == "YES"),][1:n,],
                  aes(log2FoldChange,diffMethy,label=annot.symbol))+
              theme_minimal()+
              ggtitle(paste(compName,"(DSS DMR)",sep=" "))+
              xlab(paste("log2FC, RNAseq", compName,sep="\n")) +
              ylab("Diff Methy at TSS")
          }else{
            ggplotMethyRna <- ggplot(methy_rna, aes(x=log2FoldChange, y=diffMethy,color = significance)) +
              geom_point(shape = 16, size = 1, alpha = .4) +
              scale_colour_manual(values=c("#0091ff", "#f0650e")) +
              geom_hline(yintercept=0, linetype="dashed", color = "grey")+
              geom_vline(xintercept=0, linetype="dashed", color = "grey")+
              theme_minimal()+
              ggtitle(paste(compName,"(DSS DMR)",sep=" "))+
              xlab(paste("log2FC, RNAseq", compName,sep="\n")) +
              ylab("Diff Methy at TSS")
          } #+
            #scale_color_gradient(low = "#0091ff", high = "#f0650e") +
            #ggsave("test_DSS_DMR_methy_rna_quadrantPlot.pdf")
        }
    }
	
	#options(bitmapType='cairo',device="png")
	ggsave(figpath,
          width=14,height=12,units="cm",
          plot=ggplotMethyRna) #,device="png") 
}
  
