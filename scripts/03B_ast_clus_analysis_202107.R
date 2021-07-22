###########################.
## TITLE: Astrocyte cluster analysis
## Author: Nils Briel, Feodor-Lynen-Strasse 23, 81377 Munich, Bavaria
## Date: "22/07/2021"
## Description: Here we investigate astrocyte heterogeneity and subject these findings to TFME analysis between all three disease groups.
##              We do also query DARs to rGREAT for GO analysis. 
###########################.

library(SnapATAC)
library(monocle3)
library(cicero)
library(GenomicRanges)
library(viridisLite)
library(tidyverse)
library(SummarizedExperiment)
library(rGREAT)
library(ggpubr)
library(cowplot)
library(grid)
library(parallel)
library(lemon)
library(ggrepel)
setwd("..")
source('scripts/plotViz2_v3.R')
source('scripts/04_cicero_analysis_functions.R')

# create new output sub-folder
ifelse(dir.exists(paste0('output/Ast_rev/')),
       stop(paste0("A storage directory for this project 'Ast_rev' already exists. Please edit first.")),
       dir.create('output/Ast_rev/'))
###########################.

  x.sp <- readRDS('psp_cbd_ctrl_peaks_rev.Rds')
  
  x.sp@metaData$disease <- factor(x.sp@metaData$disease, levels = c("Ctrl", "PSP", "CBD"))
  dis_cols <- c("#377EB8","#4DAF4A","#E41A1C")
  
  clus_cols <- ggpubr::get_palette(palette = "Dark2",length(levels(x.sp@cluster)))
  

#### 3B Astrocyte cluster analysis ####
  x.sp <- x.sp[which(x.sp@cluster==4),]

  x.sp@metaData$sample <- x.sp@sample
  set.seed(42)
 
  # Use CICERO for  prepare data set 
  # Preprocessing and Traj Plot
  # Run PCA then UMAP on the data
  set.seed(42)
  cicero_cds <- snapToCicero(x.sp, preprocces_cic = F, preprocess_traj = F)
  cicero_cds@metadata$line <- 'pspcbdctrl'
  cicero_cds <- preprocess_cds(cicero_cds, method = "PCA")
  cicero_cds <- align_cds(cicero_cds, alignment_group = c('case'), preprocess_method = 'PCA')
  cicero_cds <- reduce_dimension(cicero_cds, preprocess_method = "Aligned",
                                 umap.metric = "cosine",umap.min_dist = 0.33, umap.n_neighbors = 18,
                                 cores = 8) 
  monocle3::plot_cells(cicero_cds, color_cells_by = 'disease', cell_size = 1) + scale_color_manual(values = dis_cols)
  saveRDS(cicero_cds, 'output/Ast_rev/cicero_cds_ast_rev3.Rds')
  
  # transferring Ast data from cicero_cds to snap-object, and perform k-means clustering
  rown <- rownames(cicero_cds@int_colData@listData$reducedDims@listData$UMAP)
  cicero_cds@int_colData@listData$reducedDims@listData$UMAP -> x.sp@umap[,]
 
  # orientation subcluster
  x.sp@metaData$clusterkmeans <- as.factor(kmeans(x.sp@umap, centers = 6, nstart = 42, iter.max = 1000)$cluster)
  x.sp@cluster <- x.sp@metaData$clusterkmeans
  plot(cowplot::plot_grid(plotlist = list(
  plotViz2(obj=x.sp, method="umap", main="Entity",
           point.color='disease', point.size=1.8,text.add=F, text.size=4.5,legend.add=T, legend.pos = 'bottom',stroke = 0.3)+
    coord_equal() + scale_fill_manual(values = dis_cols) + coord_fixed(xlim = c(-10,10), ylim = c(-10,10)),
  plotViz2(obj=x.sp, method="umap", main="Case IDs",
           point.color='case', point.size=1.8, text.add=F, text.size=4.5,legend.add=T, legend.pos = 'bottom', stroke = 0.3)+
    coord_equal() + scale_fill_manual(values = colorspace::qualitative_hcl(n=13, palette = 'Dynamic')) +  coord_fixed(xlim = c(-10,10), ylim = c(-10,10)),
  plotViz2(obj=x.sp, method="umap", main="Clusters",
           point.color='clusterkmeans', point.size=1.8,text.add=F, text.size=4.5,legend.add=T,legend.pos = 'bottom', stroke = 0.3)+
    coord_equal() + scale_fill_manual(values = colorspace::sequential_hcl(n=length(levels(x.sp@metaData$clusterkmeans)), palette = 'RedOr'))+ coord_fixed(xlim = c(-10,10), ylim = c(-10,10))
  ), ncol = 3, nrow = 1, axis = 'tblr', align = 'hv'))
  x.sp@metaData$`Thal-Phase`<- as.factor(x.sp@metaData$`Thal-Phase`)
  
  # Plot Fig4a
  library(colorspace)
   p1 <- plotViz2(obj=x.sp[which(x.sp@metaData$disease %in% 'PSP'),], method="umap", main="PSP", 
             point.color='disease',point.size=1.2,text.add=FALSE,legend.add=TRUE, stroke = 0.09, point.alpha = 1
    ) + scale_fill_manual(values = '#4DAF4A')+coord_equal()+theme(title = element_text(size = 14), legend.position = 'none')+  coord_fixed(xlim = c(-10,10), ylim = c(-10,10))
    p2 <-  plotViz2(obj=x.sp[which(x.sp@metaData$disease %in% 'CBD'),],method="umap", main="CBD",
             point.color='disease',point.size=1.2,text.add=FALSE,legend.add=TRUE, stroke = 0.09, point.alpha = 1
    ) + scale_fill_manual(values = '#E41A1C')+coord_equal()+theme(title = element_text(size = 14), legend.position = 'none')+  coord_fixed(xlim = c(-10,10), ylim = c(-10,10))
    p3 <- plotViz2(obj=x.sp[which(x.sp@metaData$disease %in% 'Ctrl'),], method="umap", main="Ctrl",
             point.color='disease', point.size=1.2,text.add=FALSE,legend.add=TRUE, stroke = 0.09, point.alpha = 1
    ) + scale_fill_manual(values = '#377EB8')+coord_equal()+theme(title = element_text(size = 14), legend.position = 'none')+  coord_fixed(xlim = c(-10,10), ylim = c(-10,10))
    p4 <- plotViz2(obj=x.sp, method="umap", main="Clusters",
             point.color='clusterkmeans', point.size=1.2,point.alpha=1,
             text.add=TRUE, text.size=4.5,legend.add=FALSE, stroke = 0.09)+coord_equal() + scale_fill_manual(values = sequential_hcl(n=length(levels(x.sp@metaData$clusterkmeans)), palette = 'Reds')) +
      theme(title = element_text(size = 14), legend.text = element_text(size=10), legend.position = 'none')+  coord_fixed(xlim = c(-10,10), ylim = c(-10,10))
    
  tiff(paste0('output/Fig4_ast_metadata_cluster_viz_rev20210518.tiff'), units = 'px',compression = 'lzw', res = 500, width = 7000, height = 7000)
  plot(cowplot::plot_grid(plotlist = list(p1,p2,p3,p4), 
  ncol = 4, nrow = 1, axis = 'tblr', align = 'hv') )
  dev.off()
  cairo_pdf(paste0('output/Fig4_ast_metadata_cluster_viz_rev20210518.pdf'), width = 16, height = 4)
  plot(cowplot::plot_grid(plotlist = list(p1,p2,p3,p4), 
                          ncol = 4, nrow = 1, axis = 'tblr', align = 'hv') )
  dev.off()
  
  ################.
  DARs.list <- list()
  DARs.list = parallel::mclapply(levels(x.sp@cluster), mc.cores = 4, function(cluster_i){
    
    DARs = findDAR(obj=x.sp,input.mat="pmat",
      cluster.pos=cluster_i,cluster.neg.method="knn",
      bcv=0.25,test.method="exactTest", seed.use=10
      );
    DARs$FDR = p.adjust(DARs$PValue, method="BH");
    idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
    #DARs$clid <- new.cluster.ids[as.integer(cluster_i)]
    DARs.list[[cluster_i]] <<- DARs[idy,]
    DARs.list[[cluster_i]] <<- cbind(DARs.list[[cluster_i]], x.sp@peak[idy,])
  })
  
  names(DARs.list) = new.cluster.ids
  
  DARs.table <- data.table::rbindlist(DARs.list[1:5])
  
  saveRDS(DARs.list, 'output/Ast_rev/DARs.list.ast.Rds')
  
  ##########.
#### 3B.1 Astclust: chromVAR-motif Group-wise Comparisons ####
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg19)
library(gridExtra)

  library(matrixTests)
  set.seed(42)
  cpsp <- col_wilcoxon_twosample(x.sp@mmat[which(x.sp@metaData$disease == 'Ctrl'),], x.sp@mmat[which(x.sp@metaData$disease == 'PSP'),])
  ccbd <- col_wilcoxon_twosample(x.sp@mmat[which(x.sp@metaData$disease == 'Ctrl'),], x.sp@mmat[which(x.sp@metaData$disease == 'CBD'),])
  pspcbd <- col_wilcoxon_twosample(x.sp@mmat[which(x.sp@metaData$disease == 'CBD'),], x.sp@mmat[which(x.sp@metaData$disease == 'PSP'),])
  
  cpsp$bh <- p.adjust(cpsp$pvalue, 'BH', n = nrow(cpsp))
  ccbd$bh <- p.adjust(ccbd$pvalue, 'BH', n = nrow(ccbd))
  pspcbd$bh <- p.adjust(pspcbd$pvalue, 'BH', n = nrow(pspcbd))
  cpsp$feature <- rownames(cpsp)
  ccbd$feature <- rownames(ccbd)
  pspcbd$feature <- rownames(pspcbd)
  cpsp_s <- cpsp[which(cpsp$bh < 0.05),]
  ccbd_s <- ccbd[which(ccbd$bh < 0.05),]
  pspcbd_s <- pspcbd[which(pspcbd$bh < 0.05),] 
  cpsp_s$comparison <- 'PSP vs. Ctrl'
  ccbd_s$comparison <- 'CBD vs. Ctrl'
  pspcbd_s$comparison <- 'PSP vs. CBD'
  cpsp_s <- separate(cpsp_s, col = 'feature', into =c('prefix', 'TF'), sep = c("_"), remove = F)
  ccbd_s <- separate(ccbd_s, col = 'feature', into =c('prefix', 'TF'), sep = c("_"), remove = F)
  pspcbd_s <- separate(pspcbd_s, col = 'feature', into =c('prefix', 'TF'), sep = c("_"), remove = F)
  TF_mat_3comp <- data.table::rbindlist(list(cpsp_s, ccbd_s, pspcbd_s))
  TF_mat_3comp <- TF_mat_3comp[order(TF_mat_3comp$bh),]
  write.csv(TF_mat_3comp, 'output/scatac_tfs_psp_cbd_allcomparisons.csv')
  
  cpsp_s <- cpsp[which(cpsp$bh < 0.05 & ccbd$bh > 0.05 & pspcbd$bh > 0.05),] 
  ccbd_s <- ccbd[which(ccbd$bh < 0.05 & cpsp$bh > 0.05 & pspcbd$bh > 0.05),]
  pspcbd_s <- pspcbd[which(pspcbd$bh < 0.05 & cpsp$bh > 0.05 & ccbd$bh > 0.05),]
  cpsp_s$comparison <- 'PSP vs. Ctrl'
  ccbd_s$comparison <- 'CBD vs. Ctrl'
  pspcbd_s$comparison <- 'PSP vs. CBD'
  cpsp_s <- separate(cpsp_s, col = 'feature', into =c('prefix', 'TF'), sep = c("_"), remove = F)
  ccbd_s <- separate(ccbd_s, col = 'feature', into =c('prefix', 'TF'), sep = c("_"), remove = F)
  pspcbd_s <- separate(pspcbd_s, col = 'feature', into =c('prefix', 'TF'), sep = c("_"), remove = F)
  TF_mat_3exclusive<- data.table::rbindlist(list(cpsp_s, ccbd_s, pspcbd_s))
  write.csv(TF_mat_3exclusive, 'output/scatac_tfs_psp_cbd_exclcomparisons.csv')

  ## TF matrix homogeneity between clusters
  library(corrplot)
  col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                            "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                            "#4393C3", "#2166AC", "#053061"))
  x <- cbind(x.sp@mmat[,TF_mat_3exclusive$feature], x.sp@metaData$subtype)
  colnames(x)[ncol(x)] <- 'subtype'
  x <- aggregate(. ~ subtype, data = as.data.frame(x), median) %>% t() %>% .[2:nrow(.),]
  M <- cor(x, method = 'spearman')
  colnames(M) <- levels(x.sp@metaData$subtype)
  rownames(M) <- levels(x.sp@metaData$subtype)
  res1 <- cor.mtest(x, conf.level = 0.999,method= 'spearman', exact = T)
  corrplot(M, col = col2(200), method = 'shade', p.mat = res1$p, sig.level = 0.001, insig = 'blank', type = 'lower', tl.srt = 0, tl.col = "black")

  ##### GO analysis of significantly altered TFs (triangular comparison) ####
  library("pathfindR")
  library(colorspace)
  source('scripts/04_cicero_analysis_functions.R')
  annot <- getTFannotations(unique(colnames(x.sp@mmat)))
  
  levs <- levels(x.sp@metaData[,'disease'])[2:3]
  temp_list <- list()
  temp_list <- mclapply(levs, mc.cores = 4, mc.set.seed = 42, function(dis){
    mat <- x.sp@mmat[which(x.sp@metaData[,'disease'] %in% dis),] 
    mat <- t(mat[,annot$motif])
    mat_bg <- t(x.sp.ctrl@mmat[,annot$motif])
    x <- matrixTests::row_wilcoxon_twosample(as.matrix(mat), as.matrix(mat_bg)) # Wilcoxon rank-sum test
    x$FDR <- p.adjust(x$pvalue, method = 'BH', n = length(x$pvalue)) # Benjamini-Hochberg False discovery correction applied
    x$median.diff <- matrixStats::rowMedians(as.matrix(mat)) - matrixStats::rowMedians(as.matrix(mat_bg))
    x$comp <- paste0(dis, ' vs. Ctrl - adjust. method: Benjamini-Hochmerg')
    x$motif <- rownames(x)
    temp_list[[dis]] <- x
  } )
  
  # compare the extent of TFME deviations from Ctrl values between PSP and CBD:
  comp_tfme_scale <- rbind(temp_list[[1]][,c('comp', "median.diff")],
                           temp_list[[2]][,c('comp', "median.diff")])
  comp_tfme_scale$comp <- ifelse(comp_tfme_scale$comp == 'PSP vs. Ctrl - adjust. method: Benjamini-Hochmerg', 'PSP - TFs', 'CBD - TFs')
  cairo_pdf(file = "output/supplFig08_TFME_comp.pdf", width = 7, height = 7)
  plot(ggviolin(comp_tfme_scale, x = 'comp', y = 'median.diff', fill = 'comp', color = 'black', xlab = '', ylab = 'TFME(PSP|CBD) - TFME(Ctrl)',
                add = "jitter", add.params = list(color = 'black', alpha = 0.5)) +
         stat_compare_means(label = "p.format", method = "wilcox.test", hide.ns = T)+
         labs(fill = 'Entity')+ geom_hline(yintercept = 0, linetype = 2, color = 'grey40') + 
         theme_bw(base_family="Arial") +theme(legend.position="bottom")+ scale_fill_manual(values=c("#377EB8", "#E41A1C")))
  dev.off()
  
  pathPSP <- left_join(annot,temp_list[[1]], by = 'motif') %>% dplyr::select('gene_name', 'median.diff', 'FDR')
  pathCBD <- left_join(annot,temp_list[[2]], by = 'motif') %>% dplyr::select('gene_name', 'median.diff', 'FDR')
  output_psp <- run_pathfindR(pathPSP, gene_sets = 'GO-All', output_dir = 'PSP') 
  clustered_psp <- cluster_enriched_terms(output_psp, use_description = T)
  output_cbd <- run_pathfindR(pathCBD, gene_sets = 'GO-All', output_dir = 'CBD')
  clustered_cbd <- cluster_enriched_terms(output_cbd, use_description = T)
  
  enrichment_chart(result_df = clustered_cbd[order(clustered_cbd$Fold_Enrichment),][1:25,],plot_by_cluster = T, 
                   top_terms = 15)
  visualize_terms(result_df = output_psp, hsa_KEGG = FALSE, pin_name_path = "STRING")
  visualize_terms(result_df = output_cbd, hsa_KEGG = FALSE, pin_name_path = "STRING")
  set.seed(42)
  heatp <- term_gene_heatmap(output_psp[order(output_psp$Fold_Enrichment),][1:25, ], use_description = T, genes_df = pathPSP, 
                             low = diverge_hcl(n = 3, palette = 'Blue-Red 3')[1],
                             mid = diverge_hcl(n = 3, palette = 'Blue-Red 3')[2],
                             high = diverge_hcl(n = 3, palette = 'Blue-Red 3')[3])
  heatc <- term_gene_heatmap(output_cbd[order(output_cbd$Fold_Enrichment),][1:25, ], use_description = T, genes_df = pathCBD,
                             low = diverge_hcl(n = 3, palette = 'Blue-Red 3')[1],
                             mid = diverge_hcl(n = 3, palette = 'Blue-Red 3')[2],
                             high = diverge_hcl(n = 3, palette = 'Blue-Red 3')[3])
  
  upsp <- UpSet_plot(output_psp[order(output_psp$Fold_Enrichment),][1:10, ],use_description = T,
                     low = viridis(option = 'inferno',n = 3)[1],
                     mid = viridis(option = 'inferno',n = 3)[2],
                     high = viridis(option = 'inferno',n = 5)[3])
  ucbd <- UpSet_plot(output_cbd[order(output_cbd$Fold_Enrichment),][1:10, ],use_description = T,
                     low = viridis(option = 'inferno',n = 3)[1],
                     mid = viridis(option = 'inferno',n = 3)[2],
                     high = viridis(option = 'inferno',n = 5)[3])
  
  cairo_pdf(file = "output/supplFig09.pdf", width = 20, height = 12)
  plot(plot_grid(heatp, upsp, heatc, ucbd, ncol = 2, labels = c('PSP','','CBD','')))
  dev.off()
  
  saveRDS(output_psp, 'output/output_psp.pathfindR.Rds')
  saveRDS(output_cbd, 'output/output_cbd.pathfindR.Rds')
  

#### 3B.2 Astclust: rGREAT ####
saveRDS(x.sp, 'output/Ast_rev/psp_cbd_ast.Rds')  
saveRDS(list(TF_mat_3comp, TF_mat_3exclusive), 'output/Ast_rev/TF_mat_comparison_snap.Rds')  

great.ls <- list()
great.ls <- mclapply(as.numeric(levels(x.sp@cluster)), mc.cores = 4, function(i){
  DARs <- DARs.list[[i]]
  
  DARs$FDR = p.adjust(DARs$PValue, method="BH");
  print(paste0("Cluster ", i, " has its highest significance at: ", min(DARs$FDR)))
  
  idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
  
  if(dim(DARs.list[[i]])[1] > 0){
    
    job = rGREAT::submitGreatJob(
      gr                    = x.sp@peak[as.numeric(idy)],
      bg                    = NULL,
      species               = "hg19",
      includeCuratedRegDoms = TRUE,
      rule                  = "basalPlusExt",
      adv_upstream          = 5.0,
      adv_downstream        = 1.0,
      adv_span              = 1000.0,
      adv_twoDistance       = 1000.0,
      adv_oneDistance       = 1000.0,
      request_interval = 10,
      max_tries = 20,
      version = "default",
      base_url = "http://great.stanford.edu/public/cgi-bin"
    )
    
    job
    
    tb = getEnrichmentTables(job);
    
    # data cleaning
    df <- tb$'GO Cellular Component'[1:10,] %>% dplyr::select(Group = 'name', value = 'Hyper_Observed_Gene_Hits') %>% 
      mutate(Group = as.factor(Group),
             cumulative = cumsum(value),
             midpoint = cumulative - value / 2,
             label = paste0(round(value / sum(value) * 100, 1), "%"))
    
    GBP.1 <- tb[["GO Molecular Function"]] %>% .[order(.[,"Binom_Adjp_BH"]),] #%>% .[1:10,]
    GBP.2 <- tb[["GO Biological Process"]] %>% .[order(.[,"Binom_Adjp_BH"]),] #%>% .[1:10,]
    GBP.3 <- tb[["GO Cellular Component"]] %>% .[order(.[,"Binom_Adjp_BH"]),] #%>% .[1:10,]
    
    GBP.1$list <- "MF"
    GBP.2$list <- "BP"
    GBP.3$list <- "CC"
    
    great.ls[[i]] <- list(GBP.1, GBP.2, GBP.3)
    
  } else {print(paste0('Cluster ', i,' excluded from rGREAT analysis, because it lacks significant DAR-hits.'))}
})

names(great.ls) <- new.cluster.ids
names(great.res) <- new.cluster.ids
great.ls[[3]][[1]]

new.cluster.ids <- 1:6

# plot tile.plots 
plot_ls <- list()
for(k in 1:3){
  
  df <- data.frame(great.ls[[3]][k])
  df$id <- as.factor(new.cluster.ids[3])
  
   #temp <- data.frame(great.ls[[4]][k])
   #temp$id <- as.factor(new.cluster.ids[4])
      
    #df <- rbind(df, temp)
    
    #rm(temp)
  
  plot_ls[[k]] <- df
}
names(plot_ls) <- c("GO Molecular Function", "GO Biological Process", "GO Cellular Component")
saveRDS(plot_ls,'output/rgreat_results_ast_v3.Rds')
saveRDS(great.ls, 'output/Ast_rev/great.ls_ast_v3.Rds')

################## finished #####################.
print("Part 3B is done: now continue with '04_cicero_analysis.R'")
sessionInfo()