###########################.
## TITLE: SnapATAC analysis - Part 2
## Author: Nils Briel, Feodor-Lynen-Strasse 23, 81377 Munich, Bavaria
## Date: "22/07/2021"
## Description: Here we assess differentially acessible regions (DARs) of each cluster and subject these findings to a Gene Ontology database fetching package.
##              We do also apply chromVAR-motif to extract TF motif enrichments in the peaks, and finally calculate z-scores of protein degradation involvement. 
###########################.

library(SnapATAC)
library(GenomicRanges)
library(viridisLite)
library(tidyverse)
library(SummarizedExperiment)
library(ggpubr)
library(data.table)
setwd("..")

### 2.8.5 Add cell-by peak matrix to R.sp object
    x.sp <- readRDS("psp_cbd_ctrl_peaks_rev.Rds")
    
    # re-add the cell-by-bin matrix to the snap object;
    # CAVE: read the issue https://github.com/r3fang/SnapATAC/issues/154 and fkeramati's approach to correct the 'b'chr issue in the temporary fragments and bed files
    #df <- x.sp@file %>% as.data.frame() %>% separate('.', sep = 19, into = c('rm', 'keep')) %>% select(-1)
    #x.sp@file = paste0('/home/nes/',df$keep) ## paste here the specific path to the snap files
    
    x.sp <- addPmatToSnap(x.sp, do.par=F, num.cores=4);
    x.sp <- makeBinary(x.sp, mat="pmat");
    x.sp
    
    x.sp@metaData$disease <- factor(x.sp@metaData$disease, levels = c("Ctrl", "PSP", "CBD"))
    dis_cols <- c("#377EB8","#4DAF4A","#E41A1C")

##### 2.9 Identify differentially accessible regions ####
  ### 2.9.1 Define DARs in clusters by diff. analysis. 
    for(cluster_i in levels(x.sp@cluster)){
        print(cluster_i)
        DARs = findDAR(
        obj=x.sp,
        input.mat="pmat",
        cluster.pos=cluster_i,
        cluster.neg.method="knn",
        test.method="exactTest",
        bcv=0.25, 
        seed.use=10
      )
      
      DARs$FDR = p.adjust(DARs$PValue, method="BH");
      idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
      
      if(length(idy) > 0){
        
        tiff(paste0('output/2.9_DARs_cluster',cluster_i,'.tiff'), units = 'px',compression = 'lzw', res = 310, width = 3000, height = 1600) 
        par(mfrow = c(1, 2));
        plot(DARs$logCPM, DARs$logFC, 
             pch=19, cex=0.1, col="grey", 
             ylab="logFC", xlab="logCPM",
             main=paste0("Cluster ", cluster_i)
        )
        points(DARs$logCPM[idy], 
               DARs$logFC[idy], 
               pch=19, 
               cex=0.5, 
               col="red"
        )
        abline(h = 0, lwd=1, lty=2);
        
        covs = Matrix::rowSums(x.sp@pmat);
        vals = Matrix::rowSums(x.sp@pmat[,idy]) / covs;
        vals.zscore = (vals - mean(vals)) / sd(vals);
        
        plotFeatureSingle(
          obj=x.sp,
          feature.value=vals.zscore,
          method="umap", 
          main=paste0("Cluster ", cluster_i),
          point.size=0.1, 
          point.shape=19, 
          quantiles=c(0.01, 0.99),
        )
        dev.off()
      } else {print(paste0('Cluster ', cluster_i,' excluded from DAR visualization, because it lacks significant DAR-hits.'))}
    }


### 2.9.2 find DARs for each cluster
    idy.ls = lapply(levels(x.sp@cluster), function(cluster_i){
      
      DARs = findDAR(
        obj=x.sp,
        input.mat="pmat",
        cluster.pos=cluster_i,
        cluster.neg = NULL,
        cluster.neg.method="knn",
        bcv=0.25,
        test.method="exactTest",
        seed.use=10
      );
      DARs$FDR = p.adjust(DARs$PValue, method="BH");
      idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
    })
    
    DARs.list <- list()
    DARs.list = lapply(levels(x.sp@cluster), function(cluster_i){
      
      DARs = findDAR(
        obj=x.sp,
        input.mat="pmat",
        cluster.pos=cluster_i,
        cluster.neg = NULL,
        cluster.neg.method="knn",
        bcv=0.25,
        test.method="exactTest",
        seed.use=10
      );
      DARs$FDR = p.adjust(DARs$PValue, method="BH");
      idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
      DARs.list[[cluster_i]] <<- DARs[idy,]
      DARs.list[[cluster_i]] <<- cbind(DARs.list[[cluster_i]], x.sp@peak[idy,])
    })
    
    names(idy.ls) = levels(x.sp@cluster);
    names(DARs.list) = levels(x.sp@cluster);
    
    tiff(paste0('output/2.9_DARs.tiff'), units = 'px',compression = 'lzw', res = 400, width = 4500, height = 6000)  
    par(mfrow = c(4, 3));
    for(cluster_i in levels(x.sp@cluster)){
      
      print(cluster_i)
        idy =  unique(idy.ls[[cluster_i]]);
      
      if(length(idy) > 0){
      vals = Matrix::rowSums(x.sp@pmat[,idy]) / covs;
      vals.zscore = (vals - mean(vals)) / sd(vals);
      
      plotFeatureSingle(
        obj=x.sp,
        feature.value=vals.zscore,
        method="umap", 
        main=cluster_i,
        point.size=0.1, 
        point.shape=19, 
        quantiles=c(0.01, 0.99)
      )
      } else {print(paste0('Cluster ', cluster_i,' excluded from DAR visualization, because it lacks significant DAR-hits.'))}
    }
    dev.off()


saveRDS(DARs.list, 'DARs.list_rev.Rds')

# Reduce memory
if(exists('x.sp')){x.sp <- x.sp}else{x.sp <- readRDS("psp_cbd_ctrl_peaks_rev.Rds")}   
rm(list=setdiff(ls(), c("x.sp", 'DARs.list')))


##### 2.10 chromVAR-motif ####
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg19)

    register(MulticoreParam(6))
    x.sp = makeBinary(x.sp, "pmat");
    idy <- which(x.sp@peak@seqnames != 's37d5') #exclude unknown sequence
    x.sp <- x.sp[,idy,"pmat"]
    
    x.sp@mmat = suppressWarnings(SnapATAC::runChromVAR(obj=x.sp,  input.mat="pmat",  genome=BSgenome.Hsapiens.UCSC.hg19,  min.count=10,  species="Homo sapiens"))
    
   saveRDS(x.sp, 'psp_cbd_ctrl_motifs_rev.Rds')


# 2.10.1 chromVAR-motif visualization

# 2.10.2 cluster TF activity - pheatmap
  library(RColorBrewer)
  library(pheatmap)
  mat <- t(x.sp@mmat[,Matrix::colSums(x.sp@mmat) > 2.5]) %>% scale()
  colnames(mat) = paste0('cell_',1:length(x.sp@barcode))
  mat_col <- data.frame(Celltype = x.sp@metaData$celltype, Disease = x.sp@metaData$disease, Case = x.sp@metaData$case)
  rownames(mat_col) <- colnames(mat)
  mat_col <- mat_col[order(as.character(mat_col$Celltype)),]
  mat <- mat[,rownames(mat_col)]
  mat[is.na(mat)] <- 0
  idy <- which(!is.na(rownames(mat)))
  mat <- mat[idy,]
  
  # List with colors for each annotation.
  mat_colors <- list(Celltype = colorRampPalette(brewer.pal(8, "Set2"))(length(levels(mat_col$Celltype))),
                     Disease = c("#4DAF4A","#E41A1C","#377EB8"))
  names(mat_colors$Celltype) <- unique(mat_col$Celltype)
  names(mat_colors$Disease) <- unique(mat_col$Disease)

  tiff(paste0('output/2.10.4_cluster_TF_activity_heatmap_rev.tiff'), units = 'px',compression = 'lzw', res = 600, width = 8000, height = 5000)
  plot(pheatmap(
    mat               = mat,
    color             = inferno(25),
    border_color      = '',
    show_colnames     = FALSE,
    cluster_cols      = F,
    clustering_method = 'ward.D2',
    show_rownames     = T,
    scale = 'column',
    fontsize_row      = 7,
    annotation_col    = mat_col,
    annotation_colors = mat_colors,
    drop_levels       = TRUE,
    fontsize          = 10,
    main              = "TFME - all cell types"
  ))
  dev.off()

# 2.10.3 TF matrix homogeneity among clusters/celltypes 
  ## Next, we want to see the correlation of all or certain TF medians between cell types in a correlation matrix 
  library(corrplot)
  # all TFs
  col2 <- c("#BB3754FF","white", 'lightblue')
  # all included genes
  tf_mat_snap <- readRDS("TF_mat_comparison_snap.Rds") %>% .[[2]]
  x <- cbind(x.sp@mmat[,tf_mat_snap$feature], x.sp@metaData$celltype)
  colnames(x)[ncol(x)] <- 'celltype'
  x <- aggregate(. ~ celltype, data = as.data.frame(x), median) %>% t() %>% .[2:nrow(.),]
  M <- cor(x, method = 'spearman')
  colnames(M) <- levels(x.sp@metaData$celltype)
  rownames(M) <- levels(x.sp@metaData$celltype)
  res1 <- cor.mtest(x, conf.level = 0.95, method= 'spearman', exact = T)
  colnames(res1$p) <- levels(x.sp@metaData$celltype)
  rownames(res1$p) <- levels(x.sp@metaData$celltype)
  
  tiff(paste0('output/corr_tfactivities_sel.tiff'), units = 'px',compression = 'lzw', res = 300, width = 3000, height = 3000)
  plot(ggcorrplot::ggcorrplot(M, p.mat = res1$p, show.diag = F, lab = T, digits = 1,lab_size=4, colors = col2, insig = 'blank'))
  dev.off()
  
  x.sp.ctrl <- x.sp[which(x.sp@metaData$disease %in% 'Ctrl'),]
  
  temp_list <- list()
  temp_list <- mclapply(c('PSP', 'CBD'), mc.cores = 4, mc.set.seed = 42, function(dis){
    mat <- t(x.sp@mmat[which(x.sp@metaData$disease %in% dis),])
    #mat <- t(mat[,unique(c(xplan, explanatory_c, explanatory_p))])
    mat_bg <- t(x.sp.ctrl@mmat[,])
    #mat_bg <- t(mat_bg[,unique(c(xplan, explanatory_c, explanatory_p))])
    x <- row_wilcoxon_twosample(as.matrix(mat), as.matrix(mat_bg)) 
    x$FDR <- p.adjust(x$pvalue, method = 'BH', n = length(x$pvalue)) # Benjamini-Hochberg False discovery correction applied
    x$FDR[which(x$FDR>0.05)] = NA # mark those, which do not reach sign. level
    x$median.diff <- matrixStats::rowMedians(as.matrix(mat)) - matrixStats::rowMedians(as.matrix(mat_bg))
    x$dis <- dis
    x$tf <- rownames(x)
    x$updown <- ifelse(x$median.diff >= 0, 'up', 'down')
    x$comp <- paste0(dis, ' vs. Ctrl - adjust. method: Benjamini-Hochberg')
    temp_list[[dis]] <- x
  } )
  names(temp_list) = c('PSP', 'CBD')
  TFs_table <- data.table::rbindlist(temp_list) %>% .[which(complete.cases(.)),]
  # create new output sub-folder
  ifelse(dir.exists(paste0('output/Ast_rev/')),
         stop(paste0("A storage directory for this project 'Ast_rev' already exists. Please edit first.")),
         dir.create('output/Ast_rev/'))
  write.csv(TFs_table, 'output/Ast_rev/scatac_tfs_psp_cbd.csv')

##### 2.11 GREAT analysis ####
    library(ggpubr)
    library(cowplot)
    library(grid)
    library(rGREAT)
      rm(list=setdiff(ls(), c("x.sp", 'DARs.list')))
      x.sp <- SnapATAC::rmBmatFromSnap(x.sp)
      
    great_ls <- list()
### 2.11.1 rGREAT
    for(i in levels(x.sp@cluster)){
      DARs <- DARs.list[[i]]
      
      print(paste0("Cluster ", i, " has its highest significance at: ", min(DARs$FDR)))
     
       idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
      
      if(length(idy) > 0){
        
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
      
      great_ls[[paste0('cl',i)]] <- tb
     
      } else {print(paste0('Cluster ', i,' excluded from rGREAT analysis, because it lacks significant DAR-hits.'))}
    }
    
    current.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7, 9, 10)
    new.cluster.ids <- c("Oli #1 MOBP", #1
                         "Oli #2", #2
                         "Exc. ULN", #3
                         "Ast", #4
                         "Exc. DLN", #5
                         "Inh. N. SST|PVALB", #6 
                         "Mic", #7
                         "Inh. N. VIP", #9 
                         "OPC" #10
    )
    
    names(great_ls) <- plyr::mapvalues(x =names(great_ls), from = current.cluster.ids, to = new.cluster.ids)
      
    # plot tile.plots 
    plot_ls <- list()
    for(k in c("GO Molecular Function", "GO Biological Process", "GO Cellular Component")){
    
    df <- data.frame(great_ls[[1]][k])
    df$id <- names(great_ls)[1]
    
    for(i in names(great_ls)[2:length(great_ls)]){
      temp <- data.frame(great_ls[[i]][k])
      if(length(temp) > 1){
      temp$id <- as.factor(i)
      
      df <- rbind(df, temp)
      }
      rm(temp)
    }
    plot_ls[[k]] <- df
    }
    
    saveRDS(plot_ls,'output/rgreat_results_all.Rds')

  ##### 2.12 Assess protein degradation on system-level #####

    library(pheatmap)
    library(ggpubr)
    library(matrixTests)
    set.seed(42)
    x.sp <- readRDS("psp_cbd_ctrl_motifs_rev.Rds")
    
    cma_amiGO <- read.delim("refs/cma_amiGO.txt", header=FALSE, col.names = c('uniprotkb', 'gene')) %>% as.data.table() %>% dplyr::filter(gene %in% colnames(x.sp@gmat)) 
    upr_amiGO <- read.delim("refs/upr_amiGO.txt", header=FALSE, col.names = c('uniprotkb', 'gene')) %>% as.data.table() %>% dplyr::filter(gene %in% colnames(x.sp@gmat)) 
    ups_amiGO <- read.delim("refs/ups_genes_amiGO.txt", header=FALSE, col.names = c('uniprotkb', 'gene')) %>% as.data.table() %>% dplyr::filter(gene %in% colnames(x.sp@gmat)) 
    cma_amiGO$data <- c('CMA')
    upr_amiGO$data <- c('UPR')
    ups_amiGO$data <- c('UPS')
    int <- data.table::rbindlist(list(cma_amiGO,upr_amiGO, ups_amiGO)) %>% dplyr::filter(gene %in% colnames(x.sp@gmat))
    
    genes = read.table("refs/hsapiens.hg19.genes.bed");
    genes.gr = GRanges(genes[,1],IRanges(genes[,2], genes[,3]), name=genes[,7]);
    genes.sel.gr <- genes.gr[which(genes.gr$name %in% unique(int$gene))]
    set.seed(42)
    data_list <- list()
    for(ct in levels(x.sp@metaData$celltype)[levels(x.sp@metaData$celltype) != 'Oli #3']){
      set.seed(100)
      x.sp.sub <- x.sp[which(x.sp@metaData$celltype %in% ct),]
      x.sp.sub = suppressWarnings(createGmatFromMat(obj=x.sp.sub, input.mat="pmat",
                                                    genes=genes.sel.gr,do.par=TRUE,num.cores=4))
      # normalize the cell-by-gene matrix
      x.sp.sub = scaleCountMatrix(obj=x.sp.sub, mat="gmat",
                                  cov=rowSums(x.sp.sub, mat="pmat"),
                                  method = "RPM");
      # smooth the cell-by-gene matrix
      myRunMagic <- function (obj, input.mat, step.size) {
        A = obj@graph@mat;
        data.use = obj@gmat;
        # smooth
        A = A + t(A);
        A = A / DelayedMatrixStats::rowSums2(A);
        data.use.smooth = as.matrix(A) %*% as.matrix(data.use);
        if(step.size > 1){
          for(i in 1:step.size){
            data.use.smooth = A %*% data.use.smooth;
          }
        }
        slot(obj, input.mat) = Seurat::as.sparse(data.use.smooth);    
        return(obj)
      }
      x.sp.sub = myRunMagic(obj=x.sp.sub,input.mat="gmat",step.size=3);
      
      ## box plots ##
      plot_list <- list(ups_amiGO, cma_amiGO,upr_amiGO)
      for(i in 1:3){
        # assign degradation pathway-specific gene set to 'p'
        p <- plot_list[[i]] 
        # filter for significantly altered genes between each tauopathy and ctrls
        ref1 <- col_t_welch(x = as.matrix(x.sp.sub@gmat[which(x.sp.sub@metaData$disease == 'PSP'),p$gene]),
                            y = as.matrix(x.sp.sub@gmat[which(x.sp.sub@metaData$disease == 'Ctrl'),p$gene]))
        ref1$bonf <- p.adjust(ref1$pvalue, method = 'bonferroni', n = length(ref1$pvalue)) 
        ref2 <- col_t_welch(x = as.matrix(x.sp.sub@gmat[which(x.sp.sub@metaData$disease == 'CBD'),p$gene]),
                            y = as.matrix(x.sp.sub@gmat[which(x.sp.sub@metaData$disease == 'Ctrl'),p$gene]))
        ref2$bonf <- p.adjust(ref2$pvalue, method = 'bonferroni', n = length(ref2$pvalue)) 
        idy <- which(ref1$bonf < 5e-2 | ref2$bonf < 5e-2 )
        # calculate z-score by 'means'
        x <- x.sp.sub@gmat[,p$gene[idy]] %>% 
          .[rowSums(.)>0,] %>% 
          rowMeans() %>% 
          as.data.table(keep.rownames = T)
        x$dis <- x.sp.sub@metaData$disease
        x$path <- plot_list[[i]]$data[1]
        x$ct <- ct
        data_list[[paste0(ct,'_', plot_list[[i]]$data[1])]] <- x
        
      }
      
     }
    
    saveRDS(data_list, 'degra_data_list.Rds')
    
    Degra_table <- data.table::rbindlist(data_list) %>% .[which(complete.cases(.)),] 
    for(pathway in levels(as.factor(Degra_table$path))){
      degra_ref <- subset(Degra_table, path %in% pathway)
      df <- data.frame()
      for(ctype in levels(as.factor(degra_ref$ct))){
        p <- subset(degra_ref, ct %in% ctype) %>% group_by(dis) %>%
          summarise(mean = mean(.), sd = sd(.))
        p$ct <- ctype
        p$mean <- p$mean - p$mean[p$dis == 'Ctrl']
        v2 <- t.test(subset(degra_ref, ct %in% ctype & dis == 'PSP')$., subset(degra_ref, ct %in% ctype & dis == 'Ctrl')$.) 
        v3 <- t.test(subset(degra_ref, ct %in% ctype & dis == 'CBD')$., subset(degra_ref, ct %in% ctype & dis == 'Ctrl')$.) 
        p$pvalue = c(1, v2$p.value, v3$p.value)
        df <- rbind(df, p)
      }
      df$adj.pvalue <- p.adjust(df$pvalue, method = 'fdr', n = length(df$pvalue)) 
      df$updown <- ifelse(df$mean > 0, 'up', 'down')
      
      plot_list[[paste0(pathway,'_heat')]] <- plot(ggplot(df, aes(
        x = ct,
        y = dis,
        fill = mean)) + labs(fill = 'Z-score', y = '', x = '') + 
          colorspace::scale_fill_continuous_diverging("Blue-Red 3") + geom_tile(color = 'black',size = 0.5) +
          geom_text(aes(label=ifelse(pvalue <= 1e-1, signif((pvalue), 2),'')), nudge_y = 0.0, size = 3, fontface=2)+
          ggthemes::theme_tufte(base_family = 'arial', base_size = 12) + 
          theme(legend.position = 'bottom', axis.text.x =element_text(angle = 30, hjust = 1))
      )
    }

## --> see also plot_script_svg_202107.R (plotting)

################## finished #####################.
print("Part 3 is done: now continue with '03A_gchromVAR_202107.R'")
sessionInfo()

