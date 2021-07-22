#' Visualization
#'
#' @param obj A snap object.
#' @param method Visulization method c("umap").
#' @param point.size Point size [1].
#' @param point.shape Point shape type [19].
#' @param point.alpha Point transparancy level [0.8].
#' @param point.color Color of point. 
#' @param text.add Whether to add cluster label text at the centroid of each cluster [TRUE].
#' @param text.color Cluster label text color ["black"].
#' @param legend.pos Position of the legend c("bottomleft", "bottom", "left", "topleft", "top", "topright", "right", "center").
#'
#' 
plotViz2 <- function(obj, 
                     main = NULL,
                    method = 'UMAP', 
                    point.size = 1.5, 
                    point.alpha = 0.8, 
                    point.color = NULL,
                    point.color.text = NULL,
                    text.add = NULL, 
                    text.size = 4, 
                    legend.add = F,
                    legend.pos = c("bottomleft", "bottom", "left", "topleft", "top", "topright", "right", "center"),
                    feature.mode = NULL,
                    scale_color = NULL,
                    stroke = 0.1,
                    threeD = F,
                    hex.bin = F,
                    bin.size = 100
                  ){	
if(!threeD){
  if(class(obj) == 'SingleCellExperiment'){
      paste('sce-plot')
      ncell = ncol(obj)
      data.use = as.data.frame(reducedDims(obj)$UMAP)
      
      if(is.null(feature.mode)){
        data.use$col = as.factor(obj@metadata[[point.color]])
      }else{
        data.use$col = feature.mode
      }
      
    }else if(class(obj) == "snap"){
      ncell = nrow(obj);
      data.use = as.data.frame(slot(obj, method));
      colnames(data.use) = c("UMAP1", "UMAP2")
      xlims = c(-max(abs(data.use[,1])) * 1.05, max(abs(data.use[,1])) * 1.2);
      ylims = c(-max(abs(data.use[,2])) * 1.05, max(abs(data.use[,2])) * 1.05);
      
      if(is.null(feature.mode)){
      data.use$col = factor(obj@metaData[,point.color]);
      }else{
        data.use$col = feature.mode
      }
    }
  
  legend.pos = match.arg(legend.pos);
  if(!hex.bin){
  pt_plot <- ggplot(data.use, 
                    aes(x = UMAP1, 
                        y= UMAP2, 
                        fill = col),color = 'grey42') + 
    geom_point(size = point.size, 
               alpha = point.alpha, 
               shape=21, stroke = stroke) + 
    labs(title = main, 
         x = 'UMAP1', 
         y = 'UMAP2 ', 
         fill = point.color.text)  +
    theme_void(base_family="Helvetica") + theme(plot.title = element_text(hjust = 0.5))
  
    if(!is.null(feature.mode)){
      pt_plot <- pt_plot + scale_color_viridis_c()
    }
    if(!is.null(scale_color)){
      pt_plot <- pt_plot + scale_color_manual(values = scale_color)
    }
  }else{
     pt_plot <- ggplot(data.use, 
                      aes(x = UMAP1, 
                          y= UMAP2, 
                          fill = as.factor(col)),color = 'grey42') + 
      geom_hex(bins = bin.size, color = 'grey22', size = stroke) + 
      labs(title = main, 
           x = 'UMAP1', 
           y = 'UMAP2 ', 
           fill = point.color.text)  +
      theme_void(base_family="Helvetica") + theme(plot.title = element_text(hjust = 0.5)) + 
      scale_fill_manual(values = scale_color)
  }
  
  
  if(legend.add){
    pt_plot <- pt_plot + theme(legend.position=legend.pos) 
  }else{
    pt_plot <- pt_plot + theme(legend.position = "none") 
  }

  if(text.add){
    centroid <- aggregate(cbind(data.use$UMAP1, data.use$UMAP2) ~ col, data=data.use, FUN=mean)
    pt_plot <-  pt_plot + geom_label(data = centroid, mapping = aes(x=V1, y=V2, label=factor(col)), color = 'black', size = text.size) 
  }
  return(pt_plot)
}else{
  if(class(obj) == 'SingleCellExperiment'){
    paste('sce-plot')
    ncell = ncol(obj)
    data.use = as.data.frame(reducedDims(obj)$UMAP)
    
    if(is.null(feature.mode)){
      data.use$col = as.factor(obj@metadata[[point.color]])
    }else{
      data.use$col = feature.mode
    }
    
  }else if(class(obj) == "snap"){
    ncell = nrow(obj);
    data.use = as.data.frame(slot(obj, method));
    colnames(data.use) = c("UMAP1", "UMAP2", "UMAP3")
    xlims = c(-max(abs(data.use[,1])) * 1.05, max(abs(data.use[,1])) * 1.2);
    ylims = c(-max(abs(data.use[,2])) * 1.05, max(abs(data.use[,2])) * 1.05);
    zlims = c(-max(abs(data.use[,2])) * 1.05, max(abs(data.use[,2])) * 1.05);
    if(is.null(feature.mode)){
      data.use$col = factor(obj@metaData[,point.color]);
    }else{
      data.use$col = feature.mode
    }
  }
  legend.pos = match.arg(legend.pos);
  
  pt_plot <- ggplot(data.use, 
                    aes(x = UMAP1, 
                        y= UMAP2, 
                        z = UMAP3,
                        fill = col),color = 'grey42') + 
    geom_point(size = point.size, 
               alpha = point.alpha, 
               shape=21, stroke = stroke) + 
    labs(title = main, 
         x = 'UMAP1', 
         y = 'UMAP2',
         z = 'UMAP2',
         fill = point.color.text)  +
    theme_void(base_family="Helvetica") + theme(plot.title = element_text(hjust = 0.5))
  
  if(!is.null(feature.mode)){
    pt_plot <- pt_plot + scale_color_viridis_c()
  }
  if(!is.null(scale_color)){
    pt_plot <- pt_plot + scale_color_manual(values = scale_color)
  }
  
  if(legend.add){
    pt_plot <- pt_plot + theme(legend.position=legend.pos) 
  }else{
    pt_plot <- pt_plot + theme(legend.position = "none") 
  }
  
  if(text.add){
    centroid <- aggregate(cbind(data.use$UMAP1, data.use$UMAP2, data.use$UMAP3) ~ col, data=data.use, FUN=mean)
    pt_plot <-  pt_plot + geom_label(data = centroid, mapping = aes(x=V1, y=V2, z=V3, label=factor(col)), color = 'black', size = text.size) 
  }
  return(pt_plot)
}
  
}


plotHeat <- function(mat = x.sp.ast@mmat, 
                     filter_obj = TF,
                     facet_x = x.sp.ast@cluster, 
                     xlab = 'TF',
                     ylab = 'TF enrichment',
                     table = F,
                     save = T,
                     save_name = 'TF_general',
                     width = NULL,
                     height = NULL,
                     res = NULL,
                     scale =T,
                     log2 = T,
                     return_mat =F){
 
  if(table){
    mat <- as.matrix(mat[,unique(filter_obj)]) %>% cbind(.,clus = as.factor(facet_x))
    mat <- aggregate(x = mat, by = list(facet_x), FUN = mean)
    mat <- reshape2::melt(select(mat,-'clus'))
    colnames(mat)[2] <- 'gene'
    mat <- left_join(mat, filter_obj, by = 'gene')
    mat$gene <- as.factor(mat$gene)
    mat$celltype <- as.factor(mat$celltype) 
    mat <- mat[order(mat$celltype),]
    
    
    g <- ggplot(mat, aes(x=Group.1, y=gene, size = value, fill = value, group = celltype)) + 
      scale_fill_viridis_c(option = 'inferno') +
      geom_tile(color="white", size=0.1) + 
      coord_equal() + labs(x = 'Cluster', y = 'Gene') +
      ggthemes::theme_tufte(base_family="Helvetica") +
      theme(legend.position="bottom") +
      scale_y_discrete(limit = unique(mat$gene))
  
  }else{
   motif_i <-mat[1:2,] %>% as.data.frame() %>% dplyr::select(contains(filter_obj)) %>% colnames()
    mat <- as.matrix(mat[,unique(motif_i)]) %>% base::cbind(.,clus = as.factor(facet_x))
    mat <- aggregate(x = mat, by = list(facet_x), FUN = mean) 
    temp <- mat[,c(1,ncol(mat))]
    mat_2 <- mat[,-c(1,ncol(mat))]
    if(scale){
      mat_2 <- wordspace::scaleMargins(as.matrix(mat_2), 
                                       cols = 1:ncol(mat_2)
                                       #rows = 1:nrow(mat_2)
                                       )
      
    }
    if(log2){
      mat_2 <- log(mat_2+1)
      
    }
    mat <- base::cbind(mat_2,temp)
    mat <- reshape2::melt(dplyr::select(mat,-'clus'))
    
    set.seed(42)
    hc <- hclust(dist(t(mat_2)))
    ord <- colnames(mat_2)[hc$order]
    
    mat$variable <- factor(mat$variable, levels = ord)
    mat <- mat[order(mat$variable),]
    
    g <- ggplot(mat, aes(x=Group.1, y=variable, size = value, fill = value)) + 
      scale_fill_viridis_c(option = 'inferno') +
      geom_tile(color="white", size=0.1) + 
      coord_equal() + labs(x = xlab, y = ylab) +
      ggthemes::theme_tufte(base_family="Helvetica") +
      theme(legend.position="bottom") +
      scale_y_discrete(limit = unique(mat$variable)) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  if(save){
    tiff(paste0('~/projs/scATACseq_MLL/output/',save_name,'.tiff'), units = 'px',compression = 'lzw', res = res, width = width, height = height)  
    plot(g)
    dev.off()
  }
  if(return_mat){
    heat_mat <<- mat
  }
  return(plot(g))
}


plotHeatGene <- function(mat = x.sp@gmat, 
                     filter_obj,
                     facet_x = x.sp@cluster, 
                     xlab = 'Cluster',
                     ylab = 'Gene',
                     aggregate_fun = mean,
                     hc = NULL,
                     save = T,
                     save_name = 'Gene_heatmap',
                     width = NULL,
                     height = NULL,
                     res = NULL){

  mat <- as.matrix(mat[,unique(filter_obj$gene)]) %>% cbind(.,clus = as.factor(facet_x))
  mat <- stats::aggregate(x = mat, by = list(facet_x), FUN = aggregate_fun)
  mat <- reshape2::melt(select(mat,-'clus'))
  colnames(mat)[2] <- 'gene'
  mat <- left_join(mat, filter_obj, by = 'gene')
  mat$gene <- as.factor(mat$gene)
  mat$celltype <- as.factor(mat$celltype)
  mat$ct_gene <- as.factor(paste0(mat$celltype,':',mat$gene))
  mat <- mat[order(mat$celltype),]
  
  mat$Group.1 <- factor(mat$Group.1, levels = mat$Group.1[hc$order])
  
  g <- ggplot(mat, aes(x=Group.1, y=ct_gene, size = value, fill = log2(value), group = celltype)) + 
    scale_fill_viridis_c(option = 'viridis') +
    geom_tile(color="white", size=0.1) + 
    coord_equal() + labs(x = xlab, y = ylab) +
    ggthemes::theme_tufte(base_family="Helvetica") +
    theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_discrete(limit = unique(mat$ct_gene)) 
    
  plot(g)
  return(g)

  if(save){
    tiff(paste0('~/projs/scATACseq_MLL/output/',save_name,'.tiff'), units = 'px',compression = 'lzw', res = res, width = width, height = height)  
    plot(g)
    dev.off()
  }
}


plotTriang <- function(mat = x.sp@mmat, 
                       mat_bg = x.sp.ctrl@mmat, 
                       filter_obj = unique(c(xplan, explanatory_c, explanatory_p)),
                       facet_x = 'subtype', 
                       xlab = '',
                       ylab = 'TF motif enrichment',
                       title = 'Triang-plot: exclusive TF motifs', 
                       aggregate_fun = median,
                       order = 'hc',
                       save = T,
                       write_table = F,
                       save_name = 'triang_plot',
                       width = NULL,
                       height = NULL,
                       res = NULL){

  levs <- levels(x.sp@metaData[,facet_x])
  temp_list <- list()
  temp_list <- mclapply(levs, mc.cores = 4, mc.set.seed = 42, function(dis){
    mat <- x.sp@mmat[which(x.sp@metaData[,facet_x] %in% dis),] 
    mat <- t(mat[,filter_obj])
    mat_bg <- t(mat_bg[,filter_obj])
    x <- matrixTests::row_wilcoxon_twosample(as.matrix(mat), as.matrix(mat_bg)) 
    x$FDR <- p.adjust(x$pvalue, method = 'BH', n = length(x$pvalue)) # Benjamini-Hochberg False discovery correction applied
    x$FDR[which(x$FDR>0.05)] = NA # mark those, which do not reach sign. level
    x$median.diff <- matrixStats::rowMedians(as.matrix(mat)) - matrixStats::rowMedians(as.matrix(mat_bg))
    x$dis <- dis
    x$tf <- rownames(x)
    x$updown <- ifelse(x$median.diff >= 0, 'up', 'down')
    x$comp <- paste0(dis, ' vs. Ctrl - adjust. method: Benjamini-Hochmerg')
    temp_list[[dis]] <- x
  } )

  names(temp_list) = levs
  TFs_table <- data.table::rbindlist(temp_list) %>% .[which(complete.cases(.)),]
    if(write_table){
    write.csv(TFs_table, paste0('output/',save_name,'.csv'))
    }
  
  mat <- TFs_table %>% dplyr::select('dis', 'tf', 'median.diff', 'FDR') %>% pivot_wider(names_from = 'dis', values_from = c('median.diff', 'FDR'))
  mat <- column_to_rownames(mat, var='tf')
  if(order == 'hc'){
  mat[is.na(mat)] <- 1000 # to get dist function working
  set.seed(42)
  hc <- hclust(dist(mat))
  ord <- rownames(mat)[hc$order]
  }else{
    ord = order
  }
  
  TFs_table$tf <- factor(TFs_table$tf, levels = ord)
  TFs_table <- TFs_table[order(TFs_table$tf),]

  g <- ggplot(TFs_table, aes(x=dis, y=tf, size = abs(median.diff), fill = -log(FDR), group = dis, shape = updown)) + 
  scale_fill_gradient2(low="blue",mid='grey',high = "red" ) +
  geom_point() + scale_shape_manual(values = c(25, 24)) + scale_size_continuous(limits = c(0,0.5)) + scale_y_discrete(limits = rev(levels(TFs_table$tf))) +
  coord_equal() + labs(title= title, 
                       x = '', y = '',shape='Direction \n of diff.', size = 'Abs. diff. of medians \n (Group of Interest-Ctrl)') +
  ggthemes::theme_tufte(base_family="Helvetica", base_size = 12) +
  theme(legend.position="right", axis.text.x = element_text(angle = 45, hjust = 1))


  if(save){
    tiff(paste0('output/',save_name,'.tiff'), units = 'px',compression = 'lzw', res = res, width = width, height = height)  
    plot(g)
    dev.off()
  }
  
  return(plot(g))
}


