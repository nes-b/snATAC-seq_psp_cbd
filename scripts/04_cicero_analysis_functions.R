###########################
## TITLE: Cicero analysis functions
## Author: Nils Briel, Feodor-Lynen-Strasse 23, 81377 Munich, Bavaria
## Description: helper functions to process the scATAC data sets with cicero, tradeSeq, and custom vizualization functions
###########################
library(SnapATAC)
library(monocle3)
library(cicero)
library(Biobase)
library(tidyverse)
library(BiocGenerics)
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
setwd("..")
############################# some helper functions ######################################

### function to process and plot cicero-connections ---> 2 Co-Accessibility
conns_plot <- function(cicero_cds, conns = NULL, compare = NULL, genes = c('MAPT', 'STX6', 'MOBP'), gene_anno = gene_anno, CCANs = F){
  
  if(is.null(conns)){
    ## 2.1 Compute Primary Cicero Values
    conns <- run_cicero(cicero_cds, 
                        human.hg19.genome, 
                        sample_num = 100)
    head(dplyr::arrange(conns, desc(coaccess)))  
    
    saveRDS(conns, paste0('conns_',cicero_cds@metadata$dis,'.Rds'))
  }else{conns <- conns}
  
  temp_list <- list()
  if(is.null(compare)){
    temp_list <- mclapply(genes, mc.cores = 6, function(i){
      # define genomic region:
      gene_anno_sub <- subset(gene_anno, gene_name %in% i)
      start <- gene_anno_sub$start[1] - 10e5
      end <- gene_anno_sub$end[1] + 10e5
      chr <- gene_anno_sub$chromosome[1]
      
      # check for regulation events around the MAPT-locus (chr17:43,971,748-44,105,700(GRCh37/hg19))
      cairo_pdf(file =paste0('output/3_2.2_cicero_conns_plot_',cicero_cds@metadata$dis,'_',i,'.pdf'), width = 18, height = 18) 
      cicero::plot_connections(conns, chr, start, end, #1300000, 1400000, 
                               gene_model = gene_anno, 
                               gene_model_shape = 'box',
                               coaccess_cutoff = .5, 
                               connection_width = .5, 
                               collapseTranscripts = "longest", 
                               include_axis_track = TRUE,
                               alpha_by_coaccess = F)
      dev.off()
    })
  }else{
    for(i in genes){
      # define genomic region:
      gene_anno_sub <- subset(gene_anno, gene_name %in% i)
      start <- gene_anno_sub$start[1] - 5e4 #5e5
      end <- gene_anno_sub$end[1] + 5e4 #5e5
      chr <- gene_anno_sub$chromosome[1]
      
      # check for regulation events around the MAPT-locus (chr17:43,971,748-44,105,700(GRCh37/hg19))
      cairo_pdf(file =paste0('output/3_2.2_cicero_conns_plot_',cicero_cds@metadata$dis,'_',i,'.pdf'), width = 18, height = 18) 
      cicero::plot_connections(conns, chr, start, end, #1300000, 1400000, 
                               gene_model = gene_anno, 
                               gene_model_shape = 'box',
                               comparison_track = compare,
                               coaccess_cutoff = .5, 
                               connection_width = .5, 
                               collapseTranscripts = "longest", 
                               include_axis_track = TRUE,
                               alpha_by_coaccess = T)
      dev.off()
    }
  }
  
  if(CCANs){
    ## 2.4 Compute Cis-regulatory Co-Accessibility Networks
    register(MulticoreParam(8))
    CCAN_assigns <- generate_ccans(conns)
    cairo_pdf(paste0('output/3_2.2_cicero_ccans_qplot_',cicero_cds,'.tiff'),  width = 18, height = 10) 
    qplot(CCAN_assigns$CCAN)
    dev.off()
    saveRDS(CCAN_assigns, paste0('CCAN_assigns_',cicero_cds@metadata$dis,'.Rds'))
  }
}

### a wrapper function to read TF enrichment from scATAC data in snap-object
runChromVAR2 <- function(obj, x.sp = x.sp, genome = BSgenome.Hsapiens.UCSC.hg19, min.count = 10, species = "Homo sapiens"){
  data.use = t(obj@assays@data@listData$counts) 
  data.use <- data.use %>% as.matrix() %>% as.data.frame() %>% dplyr::select(!contains("s37d5")) %>% Seurat::as.sparse()
  peak.use = x.sp@peak %>% .[!seqnames(.) == "s37d5"]
  ncell = nrow(data.use)
  min.count = max(min.count, 0)
  idy = which(Matrix::colSums(data.use) >= min.count)
  
  data.use = data.use[, idy, dropping = TRUE]
  peak.use = peak.use[idy]
  
  rse <- SummarizedExperiment(assays = list(counts = t(data.use)),
                              rowRanges = peak.use, 
                              colData = DataFrame(Cell_Type = 1:nrow(data.use), 
                                                  depth = Matrix::rowSums(data.use)))
  
  rse <- addGCBias(rse, genome = genome)
  
  motifs <- getJasparMotifs(collection = "CORE", species = species)
  
  motif_mm <- matchMotifs(motifs, rse, genome = genome)
  
  dev <- computeDeviations(object = rse, annotations = motif_mm)
  dev_mat = t(assay(dev))
  
  return(dev_mat)
}

#### a helper function to identify the root principal points in trajectory analysis with cicero/monocle3
get_earliest_principal_node <- function(cds, mmat = NULL, TF = TF, genes_cut_off=500, TF_cut_off = 0.5){
  
  if(is.null(mmat)){
    cell_ids <- which(colData(cds)[, "num_genes_expressed"] > genes_cut_off)
  }else{
    motif_i <- mmat[1:2,] %>% as.data.frame() %>% dplyr::select(contains(TF)) %>% colnames()
    cds@colData$TF <- mmat[,motif_i]
    cell_ids <- which(colData(cds)[,"TF"] > TF_cut_off)
    #cicero_cds@colData$TF <<- cds@colData$TF
    #return(cds)
    print(paste0('Cells were ordered according to ', TF,' enrichment.'))
  }
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}


### convert Snap object to cicero-compatible cell data set
snapToCicero <- function(x.sp, preprocces_cic = F, preprocess_traj = F, prep = 'PCA', mat = 'pmat', umap.n_neighbors = 18, umap.min_dist = 1e-6, cores = 8, cell_bin = 5, sel=rownames(colData(cicero_cds))){    
  print('Translating snap-object to cicero.')
  ## SnapATAC input
  if(mat == 'pmat'){
    print('...using the pmat.')
      indata <- t(x.sp@pmat)
      # binarize the matrix
      # not needed, if bmat! : 
      indata@x[indata@x > 0] <- 1
      # format peak info
      peakinfo <- x.sp@peak %>% as.data.frame() %>% 
        dplyr::select(chr = 'seqnames', bp1=start, bp2=end)
      peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
      row.names(peakinfo) <- peakinfo$site_name
    }else{
      print('...using the gmat.')
      indata <- t(x.sp@gmat[,unique(colnames(x.sp@gmat))])
      genes = read.table("refs/hsapiens.hg19.genes.bed");
      genes.gr = GRanges(genes[,1], 
                         IRanges(genes[,2], genes[,3]), name=genes[,7])
      peakinfo <- genes.gr[which(genes.gr$name %in% colnames(x.sp@gmat))]
      peakinfo <- data.frame(site_name = as.vector(unique(peakinfo$name)))
      rownames(peakinfo) <- peakinfo$site_name
    }
  
  # format cell info
  print('Retrieving cell-barcode metadata.')
  cellinfo <- data.frame(x.sp@barcode) 
  row.names(cellinfo) <- cellinfo$V1
  names(cellinfo) <- "cells"
  rownames(indata) <- rownames(peakinfo)
  colnames(indata) <- rownames(cellinfo)
  # make CDS
  input_cds <-  suppressWarnings(new_cell_data_set(indata,cell_metadata = cellinfo,gene_metadata = peakinfo))
  input_cds@colData <- cbind(input_cds@colData, x.sp@metaData)
  input_cds@colData$case <- as.factor(input_cds@colData$case)

  if(preprocces_cic){
    print("Processing with cell aggregation: 'make_cicero_cds'.")
    # in case you only have peakinfo from above:
    set.seed(2017)
    
    input_cds <- estimate_size_factors(input_cds)
    input_cds <- preprocess_cds(input_cds, num_dim = 50, method = prep)
    input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', preprocess_method = prep, cores = cores,
                                  max_components = 50, umap.n_neighbors = umap.n_neighbors,umap.min_dist = umap.min_dist)
    input_cds <- cicero::make_cicero_cds(input_cds, k = cell_bin, reduced_coordinates = reducedDims(input_cds)$UMAP)
    
  }
  
  if(preprocess_traj){
    print(paste0('Processing with ', prep,', cosine and UMAP for trajectory inference.'))
    input_cds <- preprocess_cds(input_cds, method = prep, verbose = T)
    input_cds <- align_cds(input_cds, alignment_group = 'case', preprocess_method = prep)
    input_cds <- reduce_dimension(input_cds,umap.metric = "cosine",umap.min_dist = umap.min_dist, umap.n_neighbors = umap.n_neighbors, cores = cores) 
    
    # Running the clustering method. This is necessary to the construct the graph
    input_cds <- cluster_cells(input_cds, reduction_method = "UMAP", num_iter = 5,
                                random_seed = 42)
    input_cds@colData$case <- as.factor(input_cds@colData$case)
    
    # Construct the graph
    set.seed(22)
    learn_graph_control_params <- list(
      euclidean_distance_ratio = 20, 
      geodesic_distance_ratio = 10, 
      minimal_branch_len = 30,
      orthogonal_proj_tip = F,prune_graph = T, rann.k = 25)
    input_cds <- learn_graph(input_cds, learn_graph_control = learn_graph_control_params, use_partition = F)
    
    # We find all the cells that are close to the starting point
    input_cds <- order_cells(input_cds, root_pr_nodes=get_earliest_principal_node(input_cds,mmat = x.sp@mmat,TF = 'EMX2', TF_cut_off = 0.1))
    input_cds@colData$Pseudotime <- as.vector(pseudotime(input_cds))
  }
  return(input_cds)
}

### convert snap object to Single-Cell-Experiment
snapToSce <- function(x.sp, # snap object
                      counts # one of 'peaks', 'motifs', 'genes'
                      ){
  
  if(counts == 'peaks'){
    indata <- t(x.sp@pmat)
    indata@x[indata@x > 0] <- 1
    # format cell info
    cellinfo <- data.frame(x.sp@barcode) 
    row.names(cellinfo) <- cellinfo$V1
    names(cellinfo) <- "cells"
    
    # format peak info
    peakinfo <- x.sp@peak %>% as.data.frame() %>% 
      select(chr = 'seqnames', bp1=start, bp2=end)
    peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
    row.names(peakinfo) <- peakinfo$site_name
    rownames(indata) <- rownames(peakinfo)
    colnames(indata) <- rownames(cellinfo)
    sce <- SingleCellExperiment(list(counts=as.matrix(indata)),
                                colData=DataFrame(label=cellinfo),
                                rowData=DataFrame(peakinfo),
                                reducedDims = list(UMAP = x.sp@umap),
                                metadata=DataFrame(x.sp@metaData))
      return(sce)
  }else if(counts == 'motifs'){
    indata <- t(x.sp@mmat)
   
    # format cell info
    cellinfo <- data.frame(x.sp@barcode) 
    row.names(cellinfo) <- cellinfo$V1
    names(cellinfo) <- "cells"
    
    sce <- SingleCellExperiment(list(counts=as.matrix(indata)),
                                colData=DataFrame(label=cellinfo),
                                metadata=DataFrame(x.sp@metaData))
    return(sce)
  }else if(counts == 'genes'){
    indata <- t(x.sp@gmat)
    
    # format cell info
    cellinfo <- data.frame(x.sp@barcode) 
    row.names(cellinfo) <- paste0(1:nrow(cellinfo),'_',cellinfo[,1])
    names(cellinfo) <- "cells"
    
    # format peak info
    colnames(indata) <- cellinfo$x.sp.barcode
    sce <- SingleCellExperiment(list(counts=as.matrix(indata)),
                                colData=DataFrame(label=cellinfo),
                                metadata=DataFrame(x.sp@metaData))
    return(sce)
  }
}


# function to compute preprocessing dimension reduction and pseudotime from a congruent snap-cds object pair
comp_traj <- function(cicero_obj, 
                      sel = rownames(x.sp@umap),
                      snap = x.sp,
                      prep = 'PCA',
                      umap.min_dist = 0.33, 
                      umap.n_neighbors = 25,
                      root = NULL,
                      TF = 'EMX2'){
  
  cicero_obj <- estimate_size_factors(cicero_obj)
  cicero_obj <- preprocess_cds(cicero_obj,method = prep, verbose = T)
  cicero_obj <- align_cds(cicero_obj, alignment_group = 'case', preprocess_method = prep)
  cicero_obj <- reduce_dimension(cicero_obj, preprocess_method = prep, cores = 8, 
                                 umap.min_dist = umap.min_dist, umap.n_neighbors = umap.n_neighbors)
  genes = read.table("refs/hsapiens.hg19.genes.bed");
  colnames(genes) <- c('chr', 'bp1', 'bp2', 'ensemble_id', 'v', 'strand', 'gene', 'sec', 'func')
  cicero_obj <- suppressWarnings(annotate_cds_by_site(cicero_obj, genes, all=TRUE))
  head(fData(cicero_obj))
  
  cicero_obj <- cluster_cells(cicero_obj, reduction_method = 'UMAP')
  rown <- rownames(cicero_obj@int_colData@listData$reducedDims@listData$UMAP)
  cicero_obj@int_colData@listData$reducedDims@listData$UMAP <- snap@umap[sel,]
  rownames(cicero_obj@int_colData@listData$reducedDims@listData$UMAP) <- rown
  cicero_obj@clusters$UMAP$clusters <- snap@metaData$cluster
  learn_graph_control_params <- list(
    euclidean_distance_ratio = 20, #20
    geodesic_distance_ratio = 1, #0.5
    minimal_branch_len = 25,#25
    orthogonal_proj_tip = F,
    prune_graph = T,
    rann.k = 25)
  
  cicero_obj <- learn_graph(cicero_obj,learn_graph_control = learn_graph_control_params,use_partition = F, close_loop = F)
  if(is.null(root)){
    root <- get_earliest_principal_node(cicero_obj, 
                                        mmat = astclust_mmat[,], 
                                        TF = TF,
                                        TF_cut_off = 0.1)
  }else{root <- root}
  
  astclust_mmat <- snap@mmat[sel,]
  cicero_obj <- order_cells(cicero_obj, 
                            root_pr_nodes=root)
  
  cicero_obj@colData$TF <- astclust_mmat[,colnames(astclust_mmat) %like% TF]
  colnames(cicero_obj@colData)[ncol(cicero_obj@colData)] <- TF
  cicero_obj@colData$pt <- pseudotime(cicero_obj)
  return(cicero_obj)
}


### modiefied function from cicero to plot gradient plots of TF motif enrichment or other features along pseudotime
plot_traj_grad = function (cds_subset, df_order = NULL, breaks = 50, expression_family = "binomialff", size = 10, option = 'viridis', title='') 
{
  assertthat::assert_that(is(cds_subset, "cell_data_set"))
  assertthat::assert_that(assertthat::is.count(breaks))
  assertthat::assert_that(breaks >= 2)
  assertthat::assert_that(nrow(fData(cds_subset)) <= 100, msg = paste("Too many sites to plot. Be sure you are", 
                                                                      "passing only a subset of your CDS.", collapse = " "))
  pData(cds_subset)$Pseudotime <- pseudotime(cds_subset)
  min_expr = 0
  fData(cds_subset)$site_name <- row.names(fData(cds_subset))
  df <- as.data.frame(as.matrix(exprs(cds_subset)))
  df$f_id <- row.names(df)
  cds_exprs <- tidyr::gather(df, "Cell", "expression", -f_id)
  cds_exprs$expression <- as.numeric(cds_exprs$expression > 0)
  cds_exprs <- merge(cds_exprs, as.data.frame(fData(cds_subset)), 
                     by.x = "f_id", by.y = "row.names")
  cds_exprs <- merge(cds_exprs, as.data.frame(pData(cds_subset)), 
                     by.x = "Cell", by.y = "row.names")
  trend_formula <- "expression~ VGAM::sm.ns(Pseudotime, df=3)"
  merged_df_with_vgam <- plyr::ddply(cds_exprs, plyr::.(f_id), 
                                     function(x) {
                                       fit_res <- tryCatch({
                                         vg <- suppressWarnings(VGAM::vgam(formula = as.formula(trend_formula), 
                                                                           family = expression_family, data = x, maxit = 30, 
                                                                           checkwz = FALSE))
                                         res <- predict(vg, type = "response")
                                         res[res < min_expr] <- min_expr
                                         res
                                       }, error = function(e) {
                                         print("Error! Curve fit failed!")
                                         print(e)
                                         res <- rep(NA, nrow(x))
                                         res
                                       })
                                       expectation <- as.numeric(fit_res)
                                       data.frame(Pseudotime = x$Pseudotime, expectation = expectation)
                                     })
  cds_exprs$br <- cut(cds_exprs$Pseudotime, breaks = breaks)
  df <- as.data.frame(with(cds_exprs, tapply(expression, list(br,f_id), mean)))
  df$Var1 <- row.names(df)
  mean.wt <- tidyr::gather(df, "Var2", "value", -Var1)
  names(mean.wt) <- c("Var1", "Var2", "value")
  mean.wt <- cbind(mean.wt, stringr::str_split_fixed(mean.wt$Var1, ",", 2))
  names(mean.wt) <- c("interval", "feature_label", "mean", 
                      "int_start", "int_end")
  mean.wt$int_start <- as.numeric(as.character(gsub("\\(", 
                                                    "", mean.wt$int_start)))
  merged_df_with_vgam$feature_label <- as.factor(merged_df_with_vgam$f_id)
  mean.wt$mean <- mean.wt$mean * 100
  merged_df_with_vgam$expectation <- merged_df_with_vgam$expectation * 100
  
  if(is.null(df_order)){
    rank_s <- subset(merged_df_with_vgam, Pseudotime == min(merged_df_with_vgam$Pseudotime)) %>% arrange(expectation)
    rank_s$rank <- as.factor(1:nrow(rank_s))
    rank_s <- rank_s[,c('feature_label', 'rank')]
  }else{
    rank_s <- df_order
  }
  
  merged_df_with_vgam <- left_join(merged_df_with_vgam, rank_s, by ='feature_label' )
  
  if(is.null(df_order)){
    df_order_out <<- data.frame(feature_label = unique(merged_df_with_vgam$f_id), rank = unique(merged_df_with_vgam$rank))
  }
  #merged_df_with_vgam$rank <- mean(merged_df_with_vgam$expectation)
  #merged_df_with_vgam <- merged_df_with_vgam %>% arrange(rank)
  #merged_df_with_vgam$feature_label <- factor(merged_df_with_vgam$feature_label ,  levels = merged_df_with_vgam$feature_label[order(merged_df_with_vgam$rank)])
  
  
  
  g_plot <- ggplot2::ggplot(merged_df_with_vgam, aes(x=Pseudotime, y=(rank)), position = 'identity') + 
    geom_line(aes(color=expectation), size = size) +
    labs(title = title, x='Pseudotime', y='', color='Expectation')+
    scale_color_viridis_c(option = option) + 
    theme_minimal() + theme(plot.title = element_text(size=12, face = "bold")) +
    theme(legend.position = 'bottom') 
  
  return(g_plot)
}


### Function to plot a 2x2 panel of UMAP plots indicating specific TF enrichment, pseudotime, 
###     as well as binned TF enrichment for single transition steps and the TF motif sequence logo.
plotTFenrichment <- function(cds, TF, pt_cut_off=0, cut_off, breaks = 10, legend_position='bottom', explore=0, save = F, path){
  library(grid)
  library(gridExtra)
  library(lemon)
  library(ggpubr)
  library(TFBSTools)
  library(arules)
  
  if(explore){
    TFsss <- x.sp@mmat[1:2,] %>% as.data.frame() %>% dplyr::select(contains(TF)) %>% colnames()
    # loop over TFsss
    
    for(i in TFsss){
      
      # get the data
      S_matrix <- reducedDims(cds)$UMAP
      data_df <- data.frame(S_matrix[, 1:2])
      colnames(data_df) <- c("data_dim_1", "data_dim_2")
      data_df$sample_name <- row.names(data_df) %>% as.numeric()
      data_df <- as.data.frame(cbind(data_df, colData(cds))) %>%
        #separate(col = agg_cell, into = c('agg', 'row_num'), sep = 3, convert = T) %>% 
        select(data_dim_1, data_dim_2, sample_name, Size_Factor) 
      
      # prepare helper TF-enrichment scores
      TFmat <- x.sp@mmat %>% as.data.frame() %>% dplyr::select(TF = i)
      TFmat$sample_name <- 1:length(TFmat$TF) %>% as.numeric()
      TFmat <- subset(TFmat, sample_name %in% data_df$sample_name)
      TFmat$cell_index <- row.names(TFmat)
      colnames(TFmat)[1] = c('motif_enrich')
      tCDS <- full_join(data_df, TFmat, by = "sample_name") 
      row.names(tCDS) <- row.names(data_df)
      
      # plot cell-trajectories
      p = ggplot(tCDS, aes(x = data_dim_1, 
                           y= data_dim_2, 
                           color = motif_enrich)) + 
        geom_point() + 
        scale_colour_viridis_c(option = 'plasma') +
        labs(x = 'UMAP1', y = 'UMAP2', color = i) +
        theme_bw(base_size = 15); 
      p = reposition_legend(p,'bottom left')
    } # close for-loop
  } # close else
  else
  {
    
    # retrieve list of TFs queried and present in the mmat
    TFsss <- x.sp@mmat[1:2,] %>% as.data.frame() %>% dplyr::select(contains(TF)) %>% colnames()
    
    # loop over TFsss
    for(i in TFsss){
      
      # get the data
      S_matrix <- reducedDims(cds)$UMAP
      data_df <- data.frame(S_matrix[, 1:2])
      colnames(data_df) <- c("data_dim_1", "data_dim_2")
      data_df$sample_name <- row.names(data_df) %>% as.numeric()
      data_df <- as.data.frame(cbind(data_df, colData(cds))) %>%
        #separate(col = agg_cell, into = c('agg', 'row_num'), sep = 3, convert = T) %>% 
        select(data_dim_1, data_dim_2, sample_name, Size_Factor)  
      
      # prepare helper TF-enrichment scores
      TFmat <- x.sp@mmat %>% as.data.frame() %>% dplyr::select(TF = i)
      TFmat$sample_name <- 1:length(TFmat$TF) %>% as.numeric()
      TFmat <- subset(TFmat, sample_name %in% data_df$sample_name)
      TFmat$cell_index <- row.names(TFmat)
      colnames(TFmat)[1] = c('motif_enrich')
      tCDS <- full_join(data_df, TFmat, by = "sample_name") 
      row.names(tCDS) <- row.names(data_df) %>% as.numeric()
      pt <- data.frame(pt = pseudotime(cds))
      pt$sample_name <- row.names(pt) %>% as.numeric()
      tCDS <- full_join(tCDS, pt, by = "sample_name") 
      
      
      
      # define trajectories
      if(pt_cut_off==0){
        tCDS$traj <- ifelse(tCDS$data_dim_2 > cut_off, 1,2) %>% factor(levels = c(2,1))
      }else if(pt_cut_off==1){
        tCDS$traj <- ifelse(tCDS$pt > cut_off, 1,2) %>% factor(levels = c(2,1))
      }else if(pt_cut_off==2){
        tCDS$traj <- discretize(tCDS$pt, method = 'interval', breaks = breaks, labels = c(1:breaks))
      }
      
      # plot cell-trajectories
      if(pt_cut_off==0){
        p = ggplot(tCDS, aes(x = data_dim_1, 
                             y= data_dim_2, 
                             color = motif_enrich)) + 
          geom_point() + 
          scale_colour_viridis_c(option = 'plasma') +
          geom_segment(aes(y=cut_off, yend = cut_off, 
                           x= min(data_dim_1)-0.25, xend = max(data_dim_1)+0.25), 
                       linetype = 2, color = 'darkgrey') +
          labs(x = 'UMAP 1', y = 'UMAP 2', color = i) +
          theme_bw(base_size = 13) +
          theme(legend.position = legend_position)
        invisible(p)
      }else{
        p = ggplot(tCDS, aes(x = data_dim_1, 
                             y= data_dim_2, 
                             color = motif_enrich)) + 
          geom_point() + 
          scale_colour_viridis_c(option = 'plasma') +
          labs(x = 'UMAP 1', y = 'UMAP 2', color = i) +
          theme_bw(base_size = 13) +
          theme(legend.position = legend_position)
        #p = reposition_legend(p,legend_position)
        invisible(p)
      }
      if(pt_cut_off==2){
        pt_plot <- ggplot(tCDS, aes(x = data_dim_1, 
                                    y= data_dim_2, 
                                    color = traj)) + 
          geom_point() + 
          scale_colour_viridis_d(option = 'cividis') +
          labs(x = 'UMAP 1', y = 'UMAP 2', color = 'Pseudotime') +
          theme_bw(base_size = 13) +
          theme(legend.position = legend_position)
        #pt_plot = reposition_legend(pt_plot,legend_position)
        invisible(pt_plot)
      }else{
        pt_plot <- ggplot(tCDS, aes(x = data_dim_1, 
                                    y= data_dim_2, 
                                    color = pt)) + 
          geom_point() + 
          scale_colour_viridis_c(option = 'magma') +
          labs(x = 'UMAP 1', y = 'UMAP 2', color = 'Pseudotime') +
          theme_bw(base_size = 13) +
          theme(legend.position = legend_position)
        invisible(pt_plot)
      }
      
      # boxplot-comparisons - wilcox with signif. statements 
      if(pt_cut_off==2){
        b <- ggboxplot(tCDS, 
                       x = "traj",
                       xlab = 'Pseudotime bin',
                       y = c("motif_enrich"),
                       ylab = paste0("TFME ", i), 
                       palette = viridis::viridis(option = 'cividis',n = breaks),
                       add = "jitter",
                       color = 'traj',
                       ggtheme = theme_bw(base_size = 13)) +
          stat_compare_means(hide.ns = T, method = 'wilcox',
                             label = "p.signif", ref.group = '1'
          ) 
        b = ggpar(b, legend = 'none') 
      }else{
        b <- ggboxplot(tCDS, 
                       x = "traj",
                       xlab = 'Pseudotime bin',
                       y = c("motif_enrich"),
                       ylab = paste0("TFME ", i), 
                       palette = "jco",
                       add = "jitter",
                       color = 'traj',
                       ggtheme = theme_bw(base_size = 13)) +
          stat_compare_means(hide.ns = T, 
                             label = "p.signif", 
                             label.x.npc = 'center') 
        b = ggpar(b, legend = 'none') 
      }
      
      
      # retrieve seq-logo
      suppressMessages(library(JASPAR2018))
      suppressMessages(library(ggseqlogo))
      x = strsplit(i, '_')[[1]][1]
      pfm <- getMatrixByID(JASPAR2018, ID=x) %>% as.matrix()
      
      m = ggplot() + 
        geom_logo(pfm) + 
        theme_logo(base_size = 13)
      
      # merge both into one double figure       
      figure1 <- ggarrange(pt_plot, p,
                           labels = c("a", "b"),
                           ncol = 2, 
                           align = 'h',
                           widths = c(1,1))
      
      figure2 <- ggarrange(b, m,
                           labels = c("c", "d"),
                           ncol = 2, 
                           align = 'h',
                           heights = c(1,1))
      figure <- ggarrange(figure1, figure2,
                          nrow = 2, 
                          align = 'hv',
                          heights = c(1,0.8)
      )
      
      #return(tCDS)
      if(save){
        tiff(paste0(path,i,'.tiff'), units = 'px',compression = 'lzw', res = 310, width = 3000, height = 3000)
        print(figure)
        dev.off()
      } else {
        print(figure)
      }
    } # close for-loop
  } # close else
} 


### extract pseudotime and trajectory from CellDataSet
extractPtRootWeigths <- function(cds){
    # Preparation steps adhering to tradeSeq's vignette for integrating monocle3-trajectories:
    # Get the closest vertices for every cell
    y_to_cells <-  principal_graph_aux(cds)$UMAP$pr_graph_cell_proj_closest_vertex %>%
      as.data.frame()
    y_to_cells$cells <- rownames(y_to_cells)
    y_to_cells$Y <- y_to_cells$V1
    
    # Get the root vertices
    # It is the same node as above
    root <- cds@principal_graph_aux$UMAP$root_pr_nodes
    
    # Get the other endpoints
    mst <- principal_graph(cds)$UMAP
    endpoints <- names(which(igraph::degree(mst) == 1))
    endpoints <- endpoints[!endpoints %in% root]
    
    # For each endpoint
    cellWeights <- lapply(endpoints, function(endpoint) {
      # We find the path between the endpoint and the root
      path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
      path <- as.character(path)
      # We find the cells that map along that path
      df <- y_to_cells[y_to_cells$Y %in% path, ]
      df <- data.frame(weights = as.numeric(colnames(cds) %in% df$cells))
      colnames(df) <- endpoint
      return(df)
    }) %>% do.call(what = 'cbind', args = .) %>%
      as.matrix()
    rownames(cellWeights) <- colnames(cds)
    pseudotime <- matrix(pseudotime(cds), ncol = ncol(cellWeights),
                         nrow = ncol(cds), byrow = FALSE)
    meta <- list(pseudotime = pseudotime, cellWeights = cellWeights)
    return(meta)
}


#### select and annotate TFs that are present in mmat
getTFannotations <- function(motifs, genes = F){
  
  library(ensembldb)
  edb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
  
  if(!genes){
    Tfs <- separate(as.data.frame(motifs), col = 1, sep = '_', into = c('dc', 'gene_name'))
    Tfs$motif <- motifs
    annot <- ensembldb::genes(edb, columns = c("gene_name", "entrezid", 'gene_biotype'),  
                              filter = GeneNameFilter(Tfs$gene_name), return.type = 'data.frame')
    annot <- left_join(annot,Tfs[,c('motif', 'gene_name')], by = 'gene_name')
    rownames(annot) <- names(annot$entrezid)
    
    }else{
      
    Tfs = data.frame(gene_name = motifs)
    annot <- ensembldb::genes(edb, columns = c("gene_name", "entrezid", 'gene_biotype'),  
                              filter = GeneNameFilter(Tfs$gene_name), return.type = 'data.frame')
    rownames(annot) <- names(annot$entrezid) 
    }
    
    return(annot)
}
