###########################.
## TITLE: SNAPATAC ANALYSIS - Part 01
## Author: Nils Briel, Feodor-Lynen-Strasse 17, 81377 Munich, Bavaria
## Date: "22/07/2021"
## Description: This R script processes snaptools objects to R-compatible Snap-objects (1) and computes the actual analysis steps (2).
##              Here dimensionality is reduced, batch effect correction applied and cell types identified by marker gene projection.
##              GA-associations with tauopathy and other disease-gene-patterns are computed. Then peak calling is performed.
###########################.

Sys.setenv(RETICULATE_PYTHON = "/home/nes/miniconda3/bin/python") ## modify here to navigate to the specific directory
library(reticulate)
path_to_python <- "/home/nes/miniconda3/bin/python" ## modify here to navigate to the specific directory
use_python(path_to_python, required = T)
py_config()
library(SnapATAC)
library(GenomicRanges)
library(Biobase)
library(BiocParallel)
library(leiden)
library(umap)
library(viridisLite)
library(tidyverse)
library(SummarizedExperiment)
library(Matrix)
library(readxl)
library(ggpubr)
library(RColorBrewer)
library(pheatmap)
library(matrixTests)
library(parallel)
setwd("..")
source('scripts/plotViz2.R')

#### 1 Quality control & processing ####
  ##### 1.1 Merge multiple samples from snap to R .sp ####
    file.list = c("snaps/num102.snap", 
      "snaps/num104.snap", 
      "snaps/num105.snap", 
      "snaps/num108.snap", 
      "snaps/num112.snap", 
      "snaps/num113.snap", 
      "snaps/num114.snap", 
      "snaps/num115.snap",
      "snaps/numC1.snap", 
      "snaps/numC2.snap", 
      "snaps/numC3.snap", 
      "snaps/numC4.snap", 
      "snaps/numC5.snap" 
      )
    sample.list = c("num102", "num104", "num105", "num108", "num112", "num113", "num114", "num115", 
                    "numC1", "numC2", "numC3", "numC4", "numC5")
      
    x.sp = createSnap(file=file.list, sample=sample.list)
    
    x.sp.ls = lapply(seq(file.list), function(i){
      x.sp = createSnap(file=file.list[i], sample=sample.list[i]);
      x.sp
    })
    names(x.sp.ls) = sample.list;
    x.sp.ls
    
    dir.create('output/')

  ##### 1.2 Define barcode selection ####
      for(i in sample.list){
        m = strsplit(i, 'm')[[1]][2]
        barcodes = read.csv(
          #paste0("snaps/singlecell",m,".csv"), 
          paste0("",i,"/outs/singlecell.csv"),
          head=TRUE
        );
        
        barcodes = barcodes[2:nrow(barcodes),];
        promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
        
        UMI = log(barcodes$passed_filters+1, 10);
        data = data.frame(UMI=UMI, promoter_ratio=promoter_ratio);
        barcodes$promoter_ratio = promoter_ratio;
        
        p1 = ggplot(
          data, 
          aes(x= UMI, y= promoter_ratio)) + 
          geom_point(size=0.1, col="grey") +
          theme_classic() +
          ggtitle(paste0('QC:', i)) +
          ylim(0, 1) + xlim(0, 6) +
          labs(x = "log10(UMI)", y="promoter ratio") 
        
        tiff(paste0('output/1.2_barcode_sel_',i,'.tiff'), units = 'px',compression = 'lzw', res = 310, width = 1700, height = 1600)
        plot(p1)
        dev.off()
        
        barcodes.sel = barcodes[which(UMI >= 3 & UMI <= 6 & promoter_ratio >= 0.1 & promoter_ratio <= 0.7),]; 
        
        rownames(barcodes.sel) = barcodes.sel$barcode;
        write.table(barcodes.sel, paste0("snaps/",i,'barcodes.txt'))
      }
  
  
  for(i in sample.list){
  barcodes = read.csv(
    paste0("snaps/",i,"/singlecell.csv"), 
    #paste0("",i,"/outs/singlecell.csv"),
    head=TRUE
  );
  
  barcodes = barcodes[2:nrow(barcodes),];
  promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
  
  UMI = log(barcodes$passed_filters+1, 10);
  data = data.frame(UMI=UMI, promoter_ratio=promoter_ratio);
  barcodes$promoter_ratio = promoter_ratio;
  
  sel <- subset(x.after.sp@metaData, case %in% i)$barcode 
  barcodes.sel <- subset(barcodes, barcode %in% sel)
  
  rownames(barcodes.sel) = barcodes.sel$barcode;
  
  write.table(barcodes.sel, paste0("snaps/",i,'barcodes.txt'))
  }
  
  ##### 1.3. Select barcodes ####
  barcode.file.list = c(paste0("snaps/",sample.list[1],"barcodes.txt"),
                        paste0("snaps/",sample.list[2],"barcodes.txt"),
                        paste0("snaps/",sample.list[3],"barcodes.txt"),
                        paste0("snaps/",sample.list[4],"barcodes.txt"),
                        paste0("snaps/",sample.list[5],"barcodes.txt"),
                        paste0("snaps/",sample.list[6],"barcodes.txt"),
                        paste0("snaps/",sample.list[7],"barcodes.txt"),
                        paste0("snaps/",sample.list[8],"barcodes.txt"),
                        paste0("snaps/",sample.list[9],"barcodes.txt"),
                        paste0("snaps/",sample.list[10],"barcodes.txt"),
                        paste0("snaps/",sample.list[11],"barcodes.txt"),
                        paste0("snaps/",sample.list[12],"barcodes.txt"),
                        paste0("snaps/",sample.list[13],"barcodes.txt")
    );
    barcode.list = lapply(barcode.file.list, function(file){
      read.table(file)[,1];
    })
    x.sp.list = lapply(seq(x.sp.ls), function(i){
      x.sp = x.sp.ls[[i]];
      x.sp = x.sp[x.sp@barcode %in% barcode.list[[i]],];
      
    })
    names(x.sp.list) = sample.list;
    
    x.sp.list = lapply(seq(x.sp.list), function(i){
      x.sp = addBmatToSnap(x.sp.list[[i]], bin.size=5000);
      x.sp
    })

  
  ##### 1.4 Combining snap objects ####
  bin.shared = Reduce(intersect, lapply(x.sp.list, function(x.sp) x.sp@feature$name));
  x.sp.list <- lapply(x.sp.list, function(x.sp){
    idy = match(bin.shared, x.sp@feature$name);
    x.sp[,idy, mat="bmat"];
  })
  x.sp = Reduce(snapRbind, x.sp.list);
  rm(x.sp.list); # free memory
  gc();
  table(x.sp@sample);

  
  ##### 1.5 Matrix binarization & bin filtering ####
  x.sp = makeBinary(x.sp, mat="bmat");
  
  # First, we filter out any bins overlapping with the ENCODE blacklist to prevent from potential artifacts.
  black_list = read.table("refs/wgEncodeHg19ConsensusSignalArtifactRegions.bed");
  black_list.gr = GRanges(
    black_list[,1], 
    IRanges(black_list[,2], black_list[,3])
  );
  idy = queryHits(findOverlaps(x.sp@feature, black_list.gr));
  if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]}
  
  # Second, we remove unwanted chromosomes.
  chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))];
  idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature);
  if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};
  
  # Third, the bin coverage roughly obeys a log normal distribution. We remove the top 5% bins that overlap with invariant features such as promoters of the house keeping genes.
  bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
  bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
  idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
  x.sp = x.sp[, idy, mat="bmat"];
  x.sp
  

#### 2 SnapATAC analysis ####
  
  ##### 2.1 Dimensionality reduction - Nyström landmark (incompatible) ####
  x.sp@metaData$logUMI <- log(x.sp@metaData$UM)
  row.covs = log10(Matrix::rowSums(x.sp@bmat)+1);
  row.covs.dens = density(
    x = row.covs, 
    bw = 'nrd', adjust = 1
  );
  sampling_prob = 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = row.covs)$y + .Machine$double.eps); 
  set.seed(1);
  idx.landmark.ds = sort(sample(x = seq(nrow(x.sp)), size = 10000, prob = sampling_prob));
  x.landmark.sp = x.sp[idx.landmark.ds,];
  x.query.sp = x.sp[-idx.landmark.ds,];
  x.landmark.sp = runDiffusionMaps(
    obj= x.landmark.sp,
    input.mat="bmat", 
    num.eigs=50
  );
  x.query.sp = runDiffusionMapsExtension(
    obj1=x.landmark.sp, 
    obj2=x.query.sp,
    input.mat="bmat"
  );
  x.landmark.sp@metaData$landmark = 1;
  x.query.sp@metaData$landmark = 0;
  x.sp = snapRbind(x.landmark.sp, x.query.sp);
  ## combine landmarks and query cells;
  x.sp = x.sp[order(x.sp@sample),]; # IMPORTANT
  rm(x.landmark.sp, x.query.sp); # free memory
  

  ##### 2.2 Determine significant components, run graph-based clustering on uncorrected dataset & perform batch normalization ####
  set.seed(10)
  plotDimReductPW(obj=x.sp, eigs.dims=1:50, 
                  point.size=0.3,point.color="grey", point.shape=19,point.alpha=0.6,
                  down.sample=10000,pdf.file.name=NULL,pdf.height=7,pdf.width=7
                  );
  
  x.sp = runKNN(obj=x.sp,
                eigs.dims=1:25,k=15);
  
  x.sp=runCluster(obj=x.sp,tmp.folder=tempdir(),
                  louvain.lib="leiden", # sophisticated community detection algorithm
                  seed.use=10,resolution=1
                  );
  x.sp@metaData$cluster = x.sp@cluster
  
  x.sp = runViz(obj=x.sp,tmp.folder=tempdir(),
                dims=2, eigs.dims=1:25, method="umap",
                seed.use=42,num.cores	= 4, n_neighbors = 3
                );
  
  ##### 2.3 Batch effect normalization ####
  library(harmony);
  x.after.sp = runHarmony(obj=x.sp,eigs.dim=1:25, 
                          meta_data=x.sp@sample # sample index
                          );
  
  ##### 2.4 Graph-based clustering ####
  x.after.sp = runKNN(obj=x.after.sp,
                      eigs.dims=1:25,k=15
                      );
  
  x.after.sp=runCluster(obj=x.after.sp,tmp.folder=tempdir(),
                        louvain.lib="leiden", # sophisticated community detection algorithm
                        seed.use=10,resolution=1
                        );
  x.after.sp@metaData$cluster = x.after.sp@cluster
  
  # update metadata with qc metrics
  barcodes.df = data.frame()
  for (i in barcode.file.list){
    x <- read.table(i)
    df <- data.frame(x)
    barcodes.df <- rbind(barcodes.df,df)
  }
  
  barcodes.df <- barcodes.df[x.after.sp@barcode,]
  barcodes.df$case = x.after.sp@sample
  meta_df <- cbind(x.after.sp@metaData[,1:9], barcodes.df[,-1]) 
  
  df <- data.frame(case = c('num104',"num108", "num112", "num113", 'num102',"num105", "num114", "num115", 'numC1',"numC2", "numC3", "numC4", "numC5"), 
                   disease = c(rep('CBD', 4), rep('PSP',4), rep('Ctrl', 5)))

  tissue_list <- readxl::read_excel("refs/cryotissue_scatacseq_pspcbd_nils.xlsx") %>% select(.,-'NPdiagnosis') %>% as.data.frame() 
  tissue_list[is.na(tissue_list)] <- 0
  tissue_list$case <- paste0('num',tissue_list$Code)
  
  meta_df <- left_join(meta_df, tissue_list, by = 'case', keep=F) 
  meta_df<- left_join(meta_df, df, by = 'case') 
  meta_df[is.na(meta_df)] <- 0
  
  x.after.sp@metaData = subset(meta_df, barcode %in% x.after.sp@barcode)
   
   
  ##### 2.5 Metadata and cluster visualization ####
  x.after.sp = runViz(obj=x.after.sp, tmp.folder=tempdir(),
                      dims=2,eigs.dims=1:25,method="umap",
                      seed.use=42,num.cores	= 8,n_neighbors = 3
                      );
  
  # comparison before and after batch correction
  tiff(paste0('output/2.5_metadata_harmony_viz.tiff'), units = 'px',compression = 'lzw', res = 400, width = 8000, height = 3800)
  par(mfrow = c(1,2));
  plotViz(obj=x.sp, method="umap", main="Case Identifiers - before harmony",
          point.color=x.sp@cluster, point.size=0.25, text.add=FALSE, legend.add=TRUE
          );
  plotViz(obj=x.after.sp, method="umap", main="Case Identifiers - after harmony",
          point.color=x.after.sp@sample, text.add=FALSE, legend.add=TRUE
          );
  dev.off()
  
  # define colors
  dis_cols <- ggpubr::get_palette(palette = "Set1",length(levels(x.after.sp@metaData$disease)))
  clus_cols <- ggpubr::get_palette(palette = "Dark2",length(levels(x.after.sp@cluster)))

  # Metadata
  tiff(paste0('output/2.5_metadata_cluster_viz.tiff'), units = 'px',compression = 'lzw', res = 310, width = 4000, height = 3800)
  plot(cowplot::plot_grid(plotlist = list(
  plotViz2(obj=x.after.sp, method="umap", main="PSP/CBD/Ctrl Cluster",
          point.color='cluster', point.size=0.2,
          text.add=TRUE, text.size=3,legend.add=FALSE
          ),
  plotViz2(obj=x.after.sp, method="umap", main="Gender",text.add = F,
           point.size=0.2,point.color='gender',legend.add=TRUE),
  plotViz2(obj=x.after.sp, feature.mode = x.after.sp@metaData$Age, method="umap",  main="Age at death",
           point.size=0.2,text.add=FALSE, legend.add=TRUE
                    ),
  plotViz2(obj=x.after.sp,feature.mode=x.after.sp@metaData$pmi,method="umap", main="Post mortem interval",
           point.size=0.2,text.add=FALSE,legend.add=TRUE
                    ),
  plotViz2(obj=x.after.sp,method="umap", main="Braak & Braak stage",
          point.color='Braak&Braak (NFT)',point.size=0.2,text.add=FALSE,legend.add=TRUE
          ),
  plotViz2(obj=x.after.sp,method="umap", main="Thal-Phase",
          point.color='Thal-Phase', point.size=0.2,text.add=FALSE,legend.add=TRUE
          ),
  plotViz2(obj=x.after.sp[which(x.after.sp@metaData$disease %in% 'PSP'),], method="umap", main="PSP", 
          point.color='disease',point.size=0.2,text.add=FALSE,legend.add=TRUE
  ) + scale_color_manual(values = '#4DAF4A'),
  plotViz2(obj=x.after.sp[which(x.after.sp@metaData$disease %in% 'CBD'),],method="umap", main="CBD",
          point.color='disease',point.size=0.2,text.add=FALSE,legend.add=TRUE
  ) + scale_color_manual(values = '#E41A1C'),
  plotViz2(obj=x.after.sp[which(x.after.sp@metaData$disease %in% 'Ctrl'),], method="umap", main="Ctrl",
           point.color='disease', point.size=0.2,text.add=FALSE,legend.add=TRUE
  ) + scale_color_manual(values = '#377EB8')), ncol = 3, nrow = 3, axis = 'tblr', align = 'hv') )
  dev.off()
  
  # technical covariates
  tiff(paste0('output/2.5_technical_cluster_viz.tiff'), units = 'px',compression = 'lzw', res = 310, width = 4000, height = 2800)
  par(mfrow = c(2, 3));
  
  x.after.sp@metaData$readdepth <- log(x.after.sp@metaData[,"passed_filters"]+1,10)
  x.after.sp@metaData$frip <- x.after.sp@metaData$peak_region_fragments / x.after.sp@metaData$passed_filters %>% replace_na(0)
  x.after.sp@metaData$dupli <- x.after.sp@metaData$duplicate / x.after.sp@metaData$total
  
  plotFeatureSingle(obj=x.after.sp,feature.value=x.after.sp@metaData$readdepth,
                    method="umap", main="PSP/CBD/Ctrl Read Depth",
                    point.size=0.25, point.shape=19, quantiles=c(0.01, 0.99)
                    ); 
  plotFeatureSingle(obj=x.after.sp,feature.value=x.after.sp@metaData$frip,
                    method="umap", main="PSP/CBD/Ctrl FRiP",
                    point.size=0.25, point.shape=19, quantiles=c(0.01, 0.99)
                    ); 
  plotFeatureSingle(obj=x.after.sp,feature.value=x.after.sp@metaData$dupli,
                    method="umap", main="PSP/CBD/Ctrl Duplicate",
                    point.size=0.25, point.shape=19, quantiles=c(0.01, 0.99)
                    ); 
  plotFeatureSingle(obj=x.after.sp,feature.value=x.after.sp@metaData$promoter_ratio,
                    method="umap", main="PSP/CBD/Ctrl Promoter Ratio",
                    point.size=0.25, point.shape=19, quantiles=c(0.01, 0.99)
                    ); 
  plotFeatureSingle(obj=x.after.sp,feature.value=x.after.sp@metaData$enhancer_region_fragments,
                    method="umap", main="PSP/CBD/Ctrl Enhancer Ratio",
                    point.size=0.25, point.shape=19, quantiles=c(0.01, 0.99)
                    ); 
  dev.off()
  
  tiff(paste0('output/2.5_technical_pie.tiff'), units = 'px',compression = 'lzw', res = 310, width = 4000, height = 3800)
  dfr_prop <- x.after.sp@metaData %>% 
    dplyr::count(x.after.sp@metaData$disease) %>%     
    mutate(prop = prop.table(n))
  colnames(dfr_prop) <- c('DisEnt', 'value')
  plot(ggplot(dfr_prop, aes(x="", y=value, fill=DisEnt))+
    geom_bar(width = 1, stat = "identity") + theme_void() + labs(title = ' Proportions of cells in anaylsis - disease entity') +
    coord_polar("y", start=0) + scale_fill_manual(values = dis_cols))
  dev.off()

  
  ##### 2.6 Gene-based annotation ####
  # system("wget http://cf.10xgenomics.com/samples/cell-atac/1.1.0/gencode.v19.annotation.gene.bed"); # insert human reference
  genes = read.table("refs/hsapiens.hg19.genes.bed");
  genes.gr = GRanges(genes[,1], 
                     IRanges(genes[,2], genes[,3]), name=genes[,7]
  );
  
  # define TF-patterns here to visualize tau-populations
  brain_cell_types_mckenzie <- read_xlsx("refs/brain_cell_types_mckenzie.xlsx", 
                                         sheet = "top_human_enrich", skip = 1) %>% as.data.frame()
  ast.marker <- subset(brain_cell_types_mckenzie, Celltype =='ast')[1:5,]
  neu.marker <- subset(brain_cell_types_mckenzie, Celltype =='neu')[1:5,]
  oli.marker <- subset(brain_cell_types_mckenzie, Celltype =='oli')[1:5,]
  mic.marker <- subset(brain_cell_types_mckenzie, Celltype =='mic')[1:5,]
  end.marker <- subset(brain_cell_types_mckenzie, Celltype =='end')[1:5,]
  brain.markers <- c(ast.marker$gene, neu.marker$gene, oli.marker$gene, mic.marker$gene, end.marker$gene)
  
  marker.genes = data.frame(
    'Excitatory_neurons' = c("SLC17A7", "SATB2", "BCL11B", NA, "EOMES", NA, NA, NA, NA),
    'Inhibitory_neurons' = c("VIP", "PVALB", "SST", NA, NA, NA, NA, NA, NA), 
    'Neurons_undiff.' = c(NA,  NA , "CNR1", "SYNPR", NA, NA, NA, NA, NA),
    #'Endothelial_cells' = c(NA, NA "COBLL1", "APOLD1" , "SDPR",  "RGS5", "TM4SF1", "ABCG2", NA), 
    'Astrocytes' = c(NA, NA, NA, "GPR98", "AQP4", "GJA1", "ETNPPL", NA, "GFAP"), 
    'Oligodendrocytes' = c("MOBP", NA, "PLP1" ,  "TF", NA,NA, NA, NA, NA),  
    'OPCs' = c("PCDH15", "CSPG4", NA, NA, NA, NA, NA, NA, NA), #CSPG4 gene is NG2 protein
    'Microglia' = c(NA, "CXCR3", "TREM2", "CD68", "OLR1", "CCL4", NA, "GPR183", NA) 
  ) %>% pivot_longer(cols = 1:ncol(.), 
                 names_to = 'celltype', 
                 values_to = 'gene') %>% 
    drop_na() 
  marker.genes$gene <- as.vector(marker.genes$gene)
  marker.genes$celltype <- as.factor(marker.genes$celltype)
  marker.genes <- marker.genes[order(marker.genes$celltype),]
  
  # Diff Layer-specific neurons (transcriptome):
    ex_neu <- data.frame(
      'Ex1' = c('CBLN2', 'RASGRF2', NA, NA, NA, NA),
      'Ex2' = c('CUX2', 'NEFM', 'RORB', NA, NA, NA),
      'Ex3' = c('CUX2', 'NEFM', 'RORB', 'GLIS3', NA, NA),
      'Ex4' = c('RORB', 'ILRAPL2', 'TSHZ2', 'FOXP2', NA, NA),
      'Ex5' = c('RORB', 'ILRAPL2', 'FOXP2', 'PCP4', 'PDE1C', 'HS3ST5'),
      'Ex6' = c('FOXP2', 'PCP4', 'PDE1C', 'HTR2C', 'TLE4', NA)
      )%>% pivot_longer(cols = 1:ncol(.), 
                   names_to = 'celltype', 
                   values_to = 'gene') %>% 
      drop_na() 
    ex_neu$gene <- as.vector(ex_neu$gene)
    ex_neu$celltype <- as.factor(ex_neu$celltype)
    ex_neu <- ex_neu[order(ex_neu$celltype),]
  
  in_neu <- data.frame(
    'In1' = c('CCK', 'CNR1', 'RELN', 'VIP', 'TAC3', NA),
    'In2' = c('CCK', 'CNR1', 'RELN', 'CALB2', NA, NA),
    'In3' = c('TSHZ2', NA, NA, NA, NA, NA),
    'In4' = c('CCK', 'CNR1', 'RELN', 'COL5A', 'SV2C', 'EYA4'),
    'In6' = c('CA8', 'PVALB', 'RYR1', 'TAC1', NA, NA)
    ) %>%  pivot_longer(cols = 1:ncol(.), 
               names_to = 'celltype', 
               values_to = 'gene') %>% 
  drop_na() 
  in_neu$gene <- as.vector(in_neu$gene)
  in_neu$celltype <- as.factor(in_neu$celltype)
  in_neu <- in_neu[order(in_neu$celltype),]

  
  # McKenzie: Novel Marker Genes:
  atac_nov_marker <- data.frame(
  'Astrocyte' = c('ADGRV1', 'CLDN10', 'ETNPPL', 'PRSS35', 'RNF219-AS1', 'STON2', 'TPD52L1', NA, NA),
  'Microglia' = c('ARHGAP25', 'ATP8B4', 'FAM105A', 'HLA-DPA1', 'MS4A14', 'PIK3AP1', 'SH2B3', 'TRIB1', 'UBAC2'),
  'Neurons_undiff.' = c('BTBD11', 'DISP2', 'GALNTL6', 'SERTM1', 'SYT13', 'VSTM2A', 'ZMAT4', NA, NA),
  'Oligodendrocytes' = c('ANLN', 'CARNS1', 'CLCA4', 'CTNNA3', 'PAIP2B', 'QDPR', 'SLAIN1', 'SOX2-OT', 'TMEM144')
  )%>%  pivot_longer(cols = 1:ncol(.), 
                     names_to = 'celltype', 
                     values_to = 'gene') %>% 
    drop_na() 
  atac_nov_marker$gene <- as.vector(atac_nov_marker$gene)
  atac_nov_marker$celltype <- as.factor(atac_nov_marker$celltype)
  atac_nov_marker <- atac_nov_marker[order(atac_nov_marker$celltype),]

  ### 2.6.2 Marker gene enrichment visualization
  vec <- c(brain.markers, marker.genes$gene, ex_neu$gene, in_neu$gene, atac_nov_marker$gene)
  
  # plus disease associated genes from disgenet.org
  gda <- read_excel("refs/C0393570__C4551862__C0751072__C0494463__C0276496__C0949664__C0026769__C0017639_disease_gda_summary.xlsx") %>% as.data.frame()
  gda$Disease <- ifelse(gda$Disease == "Ophthalmoplegia, Progressive Supranuclear", "Progressive Supranuclear Palsy", gda$Disease)
  gda$Disease <- as.factor(gda$Disease)

  vec <- c(vec, gda$Gene) 
  genes.sel.gr <- genes.gr[which(genes.gr$name %in% vec)]
  
  marker.genes2 <- rbind(marker.genes, ex_neu, in_neu, atac_nov_marker) %>% .[order(.$celltype),]
  marker.genes2 = subset(marker.genes2, gene %in% colnames(x.after.sp@gmat))

  # create Gene Matrix
  x.after.sp = addBmatToSnap(x.after.sp, do.par=T, num.cores = 8);
  x.after.sp = suppressWarnings(createGmatFromMat(obj=x.after.sp, input.mat="bmat",
                                                  genes=genes.sel.gr,do.par=TRUE,num.cores=10
                                                  ))
  
  # normalize the cell-by-gene matrix
  x.after.sp = scaleCountMatrix(obj=x.after.sp, mat="gmat",
                                cov=rowSums(x.after.sp, mat="bmat"),
                                method = "RPM"
                                );
  
  # smooth the cell-by-gene matrix
  myRunMagic <- function (obj, input.mat, step.size) {
    A = obj@graph@mat;
    data.use = obj@gmat;
    
    # smooth
    A = A + t(A);
    A = A / Matrix::rowSums(A);
    data.use.smooth = A %*% data.use;
    if(step.size > 1){
      for(i in 1:step.size){
        data.use.smooth = A %*% data.use.smooth;
      }
    }
    
    slot(obj, input.mat) = data.use.smooth;    
    return(obj)
  }
  
  x.after.sp = myRunMagic(obj=x.after.sp,input.mat="gmat",step.size=3
                          );
  
  saveRDS(x.after.sp, file="psp_cbd_motifs_rev.Rds")
  
  # plot GA projected on cells as color in UMAP
  dir.create('output/marker_genes')
  
  atac_marker <- data.frame(
    'Astrocyte' = c('CLDN10', 'ETNPPL', 'STON2', 'TPD52L1', 'GJA1', 'GFAP'),
    'Microglia' = c('TREM2', 'CD68', 'OLR1', 'HLA-DPA1', 'CCL4', 'PIK3AP1'),
    'Exc. ULN' = c('CBLN2', 'CUX2', 'RASGRF2', NA, NA, NA),
    'Exc. DLN' = c('RORB', 'FOXP2', 'TSHZ2', NA, NA, NA),
    'Exc. Neu.' = c('NEFM', 'BCL11B', 'SATB2', NA, NA, NA),
    'Inh. VIP' = c('VIP', 'TAC3', 'CALB2', NA, NA, NA),
    'Inh. PVALB|SST' = c('PVALB', 'SST', 'TAC1', NA, NA, NA),
    'Oligodendro' = c('MOBP', 'CLCA4', 'PAIP2B', 'QDPR', 'SOX2-OT', 'TMEM144')
  )%>%  pivot_longer(cols = 1:ncol(.), 
                     names_to = 'celltype', 
                     values_to = 'gene') %>% 
    drop_na() 
  
  atac_marker$gene <- as.vector(atac_marker$gene)
  atac_marker$celltype <- as.factor(atac_marker$celltype)
  atac_marker <- atac_marker[order(atac_marker$celltype),]
  saveRDS(atac_marker, 'output/marker_genes/atac_marker.Rds')
  write_csv(atac_marker, 'output/marker_genes/atac_marker.csv')
  
  marker_genes_ls <- list()
  for(k in levels(atac_marker$celltype)){
    temp <- list()
    mg <- subset(atac_marker, celltype == k)$gene
    for(i in mg){
      celltype = subset(atac_marker, gene %in% i)$celltype
      tiff(paste0('output/marker_genes/2.6_marker_gene_',i,'.tiff'), units = 'px',compression = 'lzw', res = 310, width = 1700, height = 1600)  
      plotFeatureSingle(obj=x.after.sp,feature.value=x.after.sp@gmat[, i],
                        method="umap", main=paste0(k,': ', i),
                        point.size=0.1, point.shape=19, quantiles=c(0, 1)
                        )
      dev.off()
      p <- plotViz2(obj=x.after.sp,feature.mode = x.after.sp@gmat[, i],
                    method="umap", main=paste0(k,': ', i),
                    point.size=0.1, text.add = F, legend.add = T, legend.pos = ) + coord_equal()
      temp[[i]] <- p
      }
    marker_genes_ls[[k]] <- temp 
    rm(temp)
  }
  
  library(cowplot)
  for(k in names(marker_genes_ls)){
    temp <- marker_genes_ls[[k]]
    if(length(temp)==3){vec = c(3000, 900)}else{vec = c(3000, 1900)}
    tiff(paste0('output/marker_genes/2.6.0_marker_genes_grid_',k,'.tiff'), units = 'px',compression = 'lzw', 
         res = 350, width = vec[1], height = vec[2])
    plot(cowplot::plot_grid(plotlist = temp, axis = 'tblr', ncol = 3))
    dev.off()
  }
  
  ##### 2.7 Hierarchical clustering ####
  set.seed(42)
  rm(list=setdiff(ls(), c("x.after.sp", 'dis_cols', 'atac_marker', 'gda', 'syn')))
  source('scripts/plotViz2.R')
  # calculate the ensemble signals for each cluster
  ensemble.ls = mclapply(split(seq(length(x.after.sp@cluster)), x.after.sp@cluster), mc.cores = 1, function(x){
    SnapATAC::colMeans(x.after.sp[x,], mat="bmat");
  })
  saveRDS(ensemble.ls, 'ensemble.ls.Rds')
  # cluster using 1-cor as distance  
  set.seed(42)
  dist_mat <- 1 - cor(t(do.call(rbind, ensemble.ls)))
  hc <- hclust(as.dist(dist_mat,diag = TRUE, upper = TRUE), method="centroid")
  sil_cl <- cluster::silhouette(cutree(hc,k = 2) ,as.dist(dist_mat), title=title(main = 'Good'))
  rownames(sil_cl) <- rownames(dist_mat)
  rm(ensemble.ls)
  
  library(factoextra)
  set.seed(42)
  idcells <- as.integer(sample(rownames(x.after.sp@metaData),1000))
  idy <- which(Matrix::colMeans(x.after.sp@bmat) == quantile(Matrix::colMeans(x.after.sp@bmat))[3] | 
                 Matrix::colMeans(x.after.sp@bmat) == quantile(Matrix::colMeans(x.after.sp@bmat))[4] | 
                 Matrix::colMeans(x.after.sp@bmat) == quantile(Matrix::colMeans(x.after.sp@bmat))[5]  ) # quantile-based and column-sums-based selection of most important 
  # get bmat
  bmat <- as.matrix(x.after.sp@bmat[idcells,idy])
  res.hc <- eclust(as.matrix(bmat), 
                   FUNcluster="hclust",hc_method = 'centroid',k.max = 15,
                   seed = 42)
  ct_cols <- scales::hue_pal()(length(levels(x.after.sp@metaData$cluster)))
  dend <- fviz_dend(hc,k = 11, show_labels = T, k_colors = 'black') 
  gap <- fviz_gap_stat(res.hc$gap_stat) + geom_vline(xintercept = 11) + ggthemes::theme_base()

  stage_cells <- list(dend = dend, gap = gap)
  for(cellnum in c(500, 1000, 5000, 1e4)){
  idcells <- as.integer(sample(rownames(x.after.sp@metaData),cellnum))
  dat <- cbind(as.data.frame(x.after.sp@umap)[idcells,],
                       x.after.sp@metaData[idcells,c('celltype', 'disease')])
  p <- ggpubr::ggscatter(dat, "umap-1", "umap-2", fill = "celltype",color = 'grey42',title = paste0('n: ',cellnum,' cells'),
                         shape = 21, size = 1.5, repel = FALSE, 
                         ellipse = TRUE, ellipse.type = "convex", 
                         ellipse.level = 0.5, ellipse.alpha = 0.2, labelsize = 10, main = "Cluster plot", 
                         xlab = NULL, ylab = NULL, outlier.color = "black",facet.by = 'disease',
                         outlier.pointsize = 0.1, outlier.labelsize = 0.1, ellipse.border.remove = T,
                         ggtheme = ggthemes::theme_base()) + coord_equal()
  stage_cells[[paste('plot',cellnum)]] <- p
  }
  saveRDS(stage_cells, 'output/stage_cells_fig2.Rds')
  tiff(paste0('output/2.7_hc_clust.tiff'), units = 'px',compression = 'lzw', res = 250, width = 1700, height = 1600)  
  plot(hc, hang=-1, xlab="")
  dev.off()
  
 
  # 2.7.2 cell types
  current.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
  new.cluster.ids <- c("Oli #1 MOBP", #1
                       "Oli #2", #2
                       "Exc. ULN", #3
                       "Ast", #4
                       "Exc. DLN", #5
                       "Inh. N. SST|PVALB", #6 
                       "Mic", #7
                       "Exc. N. NEFM|BCL11B", #8
                       "Inh. N. VIP", #9 
                       "OPC", #10
                       "Oli #3" #1
                       )
  
  x.after.sp@metaData$celltype <- plyr::mapvalues(x = x.after.sp@cluster, from = current.cluster.ids, to = new.cluster.ids)
  x.after.sp@metaData$celltype <- x.after.sp@metaData$celltype %>% as.character() %>% as.factor()
  ct_cols <- scales::hue_pal()(length(levels(x.after.sp@metaData$celltype)))


  tiff(paste0('output/2.7.2_celltypes.tiff'), units = 'px',compression = 'lzw', res = 450, width = 3400, height = 3200)
  plotViz2(obj = x.after.sp,method = "umap",main = "PSP/CBD/Ctrl cell types",
           point.color = 'celltype',point.color.text = x.after.sp@metaData$celltype,legend.add = F,text.add = T)+
  coord_fixed()+ scale_color_manual(values = ct_cols)
  dev.off()
  

  ### 2.7.2.2 Celltype_frequencies - rel-boxplots
    set.seed(42)
    x.cbd <- subset(x.after.sp@metaData, disease == 'CBD')
    x.psp <- subset(x.after.sp@metaData, disease == 'PSP')
    x.ctrl <- subset(x.after.sp@metaData, disease == 'Ctrl')

    x.ctrl <- x.ctrl %>% dplyr::count(as.character(.$celltype)) %>%     
      mutate(prop = prop.table(n))
    
    x.psp_list <- list()
    x.psp_list <- mclapply(levels(as.factor(x.psp$case)), mc.cores = 4, function(i){
      x <- subset(x.psp, case == i)
      x <- x %>% dplyr::count(as.character(.$celltype)) %>%     
        mutate(prop = prop.table(n))
      if(nrow(x) == nrow(x.ctrl)){
        x$prop_to_ctrl <- x$prop/x.ctrl$prop
        x$prop_ctrl <- x.ctrl$prop
        
      }else{
        x <- left_join(select(x.ctrl, c(-prop, -n)), x, by = 'as.character(.$celltype)')
        x[is.na(x)] <- 0
        x$prop_to_ctrl <- x$prop/x.ctrl$prop
        x$prop_ctrl <- x.ctrl$prop
      }
      colnames(x)[1] <- 'celltype'
      x$case <- i
      x.psp_list[[i]] <- x
    }
    )
    x.psp <- data.table::rbindlist(x.psp_list)
    
    x.cbd_list <- list()
    x.cbd_list <- mclapply(levels(as.factor(x.cbd$case)), 
                           mc.cores = 4, 
                           function(i){
      x <- base::subset(x.cbd, case == i)
      x <- x %>% dplyr::count(as.character(.$celltype)) %>%     
        dplyr::mutate(prop = prop.table(n))
      if(nrow(x) == nrow(x.ctrl)){
        x$prop_to_ctrl <- x$prop/x.ctrl$prop
        x$prop_ctrl <- x.ctrl$prop
      
      }else{
        x <- left_join(dplyr::select(x.ctrl, c(-prop, -n)), x, by = 'as.character(.$celltype)')
        x[is.na(x)] <- 0
        x$prop_to_ctrl <- x$prop/x.ctrl$prop
        x$prop_ctrl <- x.ctrl$prop
      }
      colnames(x)[1] <- 'celltype'
      x$case <- i
      x.cbd_list[[i]] <- x
     }
    )
    
    x.cbd <- data.table::rbindlist(x.cbd_list)
    
    x <- data.frame(table(x.after.sp@metaData$celltype))
    colnames(x)[1] <- 'Celltype'
    saveRDS(list(x.psp = x.psp, x.cbd = x.cbd, tot_ct=x), 'output/celltype_freq.Rds')
      

  ## 2.7.4 cluster gene activity - pheatmap
    vec <- subset(atac_marker, gene %in% colnames(x.after.sp@gmat))$gene
    mat <- t(x.after.sp@gmat[,vec]) %>% scale()
    colnames(mat) = paste0('cell_',1:length(x.after.sp@barcode))
    
    mat_col <- data.frame(Diagnosis = x.after.sp@metaData$disease, 
                          Celltype = x.after.sp@metaData$celltype)
 
    rownames(mat_col) <- colnames(mat)
    mat_col <- mat_col[order(as.character(mat_col$Celltype)),]
    mat <- mat[,rownames(mat_col)]
    newnames <- lapply(
      rownames(mat),
      function(x) bquote(italic(.(x)))) #credits and thanks to Kevin Blighe https://www.biostars.org/p/400381/
    
    mat_colors <- list(Celltype = ct_cols,
                       Diagnosis = c("#4DAF4A","#E41A1C", "#377EB8"))
    names(mat_colors$Celltype) <- unique(mat_col$Celltype)
    names(mat_colors$Diagnosis) <- unique(mat_col$Diagnosis)
    
    tiff(paste0('output/2.7.4_cluster_gene_activity_heatmap_20210518_2.tiff'), 
          units = 'px',compression = 'lzw', res = 450, width = 3700, height = 2700)
      
    plot(
      pheatmap(
        mat               = as.data.frame((mat)), 
        color             = colorspace::diverge_hcl(256,palette = 'Blue-Red 3', power = 0.6),
        scale = 'row',
        border_color      = NA,
        show_colnames     = F,
        cluster_cols      = F,
        clustering_method = 'ward.D2',
        clustering_distance_rows = "correlation",
        show_rownames     = T,
        labels_row = as.expression(newnames),
        fontsize_row      = 10,
        annotation_col    = mat_col,
        annotation_colors = mat_colors,
        drop_levels       = TRUE,
        fontsize          = 10,
        main              = "Gene Accessibility - Marker Genes"
      )
    )
    dev.off()
      

  #### 3 Disease-associated GA patterns: ####
    
    # create new output sub-folder
    ifelse(dir.exists(paste0('output/disease_ass/')),
           stop(paste0("A storage directory for this project 'gda' already exists. Please edit first.")),
           dir.create('output/disease_ass/'))

    # generate disease-wise and -wise celltype X gene-acitvity matrices for heatmaps

    plot_list=list()
    table_list=list()
    
    x.sp.ctrl <- x.after.sp[which(x.after.sp@metaData$disease %in% 'Ctrl'),]
    
    for(i in c('PSP', 'CBD')){
      x.sp.dis <- x.after.sp[which(x.after.sp@metaData$disease %in% i),]
          gda_s <- subset(gda, Disease %in% 'Tauopathies' & Gene %in% colnames(x.after.sp@gmat))
          if(nrow(gda_s) > 50){
            gda_s <- gda_s[1:50,]
          }
          
          # compute wilcox between PSP or CBD against Ctrls among a defined set of GDA-genes  
          
          temp_list <- list()
          temp_list <- mclapply(levels(x.sp.ctrl@metaData$celltype), mc.cores = 4, mc.set.seed = 42, function(ct){
            
            mat <- x.sp.dis@gmat[which(x.sp.dis@metaData$celltype %in% ct),] 
            mat <- t(mat[,gda_s$Gene])
            
            mat_bg <- x.sp.ctrl@gmat[which(x.sp.ctrl@metaData$celltype %in% ct),] 
            mat_bg <- t(mat_bg[,gda_s$Gene])
            
            x <- row_wilcoxon_twosample(as.matrix(mat), as.matrix(mat_bg)) 
            x$FDR <- p.adjust(x$pvalue, method = 'bonferroni', n = length(x$pvalue)) # Bonferroni False discovery correction applied
            x$FDR[which(x$FDR>0.05)] = NA # mark those, which do not reach sign. level
            x$median.diff <- matrixStats::rowMedians(as.matrix(mat)) - matrixStats::rowMedians(as.matrix(mat_bg))
            x$ct <- ct
            x$gene <- rownames(x)
            x$updown <- ifelse(x$median.diff >= 0, 'up', 'down')
            x$comp <- paste0(ct, ' vs. Ctrl - adjust. method: Bonferroni')
            x$median.diff[which(abs(x$median.diff)<0.1)] = NA # mark those, that do not reach sign. level
            temp_list[['Tauopathies']] <- x
            } )
          names(temp_list) = levels(x.sp.ctrl@metaData$celltype)
          
          gda_table <- data.table::rbindlist(temp_list) #%>% .[which(complete.cases(.)),]
          gda_table$dis <- i
          gda_table$ct_gene <- paste0(gda_table$ct,':',gda_table$gene)
          table_list[[paste0(i,':Tauopathies')]] <- gda_table
          
          # plot the significances and mean differences 
          g <- ggplot(gda_table, aes(x=ct, y=gene, fill = ifelse(-log(FDR)>-log(0.05), median.diff, 0), group = ct)) + 
            colorspace::scale_fill_continuous_diverging(palette = 'Blue-Red', p1 = 0.5) +
            geom_tile(color = 'grey22', size = 0.25) + scale_shape_manual(values = c(25, 24)) + scale_size_continuous(limits=c(0,5))  +
            coord_equal() + labs(title= paste0(i,':Tauopathies'), x = '', y = '',fill = 'GA difference\n (Dis.Ent.-Ctrl)') +
            ggthemes::theme_tufte(base_family="Arial") + scale_y_discrete(limits=rev) +
            theme(legend.position="right", axis.text.x = element_text(angle = 45, hjust = 1),  axis.text.y = element_text(face = 'italic'))
        
          plot_list[[paste0(i,':Tauopathies')]] <- g
          tiff(paste0('output/disease_ass/2.7.4_dag_',paste0(i,':Tauopathies'),'_heat_v4.tiff'), units = 'px',compression = 'lzw', res = 450, width = 5000, height = 3000)
          plot(g)
          dev.off()
          }
    
    saveRDS(table_list, 'output/disease_ass/dis_ass_table_list_v4.Rds')
    
    tab <- data.table::rbindlist(table_list[c('CBD:Tauopathies', 'PSP:Tauopathies')])
    write.csv(tab, 'scatac_ga_psp_cbd_v4.csv')
    saveRDS(plot_list, 'output/disease_ass/dis_ass_plot_list_v4.Rds')
    
    
  ##### Microglial tauopathy genes? ####
    library("pathfindR")
    source('scripts/04_cicero_analysis_functions.R')
    gda_s <- subset(gda, Disease %in% 'Tauopathies' & Gene %in% colnames(gmat))#[1:50,]
    annot <- getTFannotations(colnames(x.after.sp@gmat)[colnames(x.after.sp@gmat) %in% gda_s$Gene], genes = T)
    levs <- levels(x.after.sp@metaData[,'disease'])[2:3]
    temp_list <- list()
    temp_list <- mclapply(levs, mc.cores = 4, mc.set.seed = 42, function(dis){
      mat <- x.after.sp@gmat[which(x.after.sp@metaData[,'disease'] %in% dis & x.after.sp@metaData[,'celltype'] %in% 'Mic'),] 
      mat <- t(mat[,gda_s$Gene])
      mat_bg <- t(x.after.sp@gmat[which(x.after.sp@metaData[,'disease'] %in% 'Ctrl'& x.after.sp@metaData[,'celltype'] %in% 'Mic'),gda_s$Gene])
      x <- matrixTests::row_wilcoxon_twosample(as.matrix(mat), as.matrix(mat_bg)) # Wilcoxon rank-sum test
      x$FDR <- p.adjust(x$pvalue, method = 'BH', n = length(x$pvalue)) # Benjamini-Hochberg False discovery correction applied
      x$median.diff <- matrixStats::rowMedians(as.matrix(mat)) - matrixStats::rowMedians(as.matrix(mat_bg))
      x$comp <- paste0(dis, ' vs. Ctrl - adjust. method: Benjamini-Hochmerg')
      x$gene_name <- rownames(x)
      temp_list[[dis]] <- x
    } )
    
    pathPSP <- left_join(annot,temp_list[[1]], by = 'gene_name') %>% dplyr::select('gene_name', 'median.diff', 'FDR')
    pathCBD <- left_join(annot,temp_list[[2]], by = 'gene_name') %>% dplyr::select('gene_name', 'median.diff', 'FDR')
    pathPSP$FDR[is.na(pathPSP$FDR)] <- 1
    pathCBD$FDR[is.na(pathCBD$FDR)] <- 1
    output_psp <- run_pathfindR(pathPSP, gene_sets = 'GO-All', output_dir = 'PSP_mic') 
    output_cbd <- run_pathfindR(pathCBD, gene_sets = 'GO-All', output_dir = 'CBD_mic')
    enrichment_chart(result_df = clustered_cbd[order(clustered_cbd$Fold_Enrichment),][1:25,],plot_by_cluster = T, 
                     top_terms = 15)
    visualize_terms(result_df = output_psp, hsa_KEGG = FALSE, pin_name_path = "STRING")
    visualize_terms(result_df = output_cbd, hsa_KEGG = FALSE, pin_name_path = "STRING")
    set.seed(42)
    
    ### --> no GO enrichment to report when taking all 263 tauopathy genes into account!
    genes = read.table("refs/hsapiens.hg19.genes.bed");
    genes.gr = GRanges(genes[,1], 
                       IRanges(genes[,2], genes[,3]), name=genes[,7])
    #http://amigo.geneontology.org/amigo/search/bioentity?q=*:*&fq=regulates_closure:%22GO:0001774%22&sfq=document_category:%22bioentity%22
    mic_act <- read.delim("refs/mg_positive_activation_amiGO.tsv", header=FALSE)
    genes.sel.gr <- genes.gr[which(genes.gr$name %in% mic_act$V2)]
    # create Gene Matrix
    x.after.sp = suppressWarnings(createGmatFromMat(obj=x.after.sp, input.mat="bmat",
                                                    genes=genes.sel.gr,do.par=TRUE,num.cores=8))
    # normalize the cell-by-gene matrix
    x.after.sp = scaleCountMatrix(obj=x.after.sp, mat="gmat",
                                  cov=rowSums(x.after.sp, mat="bmat"),
                                  method = "RPM")
    # smooth the cell-by-gene matrix
    myRunMagic <- function (obj, input.mat, step.size) {
      A = obj@graph@mat;
      data.use = obj@gmat;
      A = A + t(A);
      A = A / Matrix::rowSums(A);
      data.use.smooth = A %*% data.use;
      if(step.size > 1){
        for(i in 1:step.size){
          data.use.smooth = A %*% data.use.smooth;
        }
      }
      slot(obj, input.mat) = data.use.smooth;    
      return(obj)
    }
    x.after.sp = myRunMagic(obj=x.after.sp,input.mat="gmat",step.size=3)
    df <- x.after.sp@gmat[x.after.sp@metaData$celltype == 'Mic',]
    df <- base::cbind(as.matrix(df), x.after.sp@metaData[x.after.sp@metaData$celltype == 'Mic', 'disease'])
    pl_df <- cbind(as.data.table(df), dis = x.after.sp@metaData[x.after.sp@metaData$celltype == 'Mic', 'disease'])
    pl_df$mean_score <- (rowMeans2(df[,-ncol(df)]) - mean(rowMeans2(df[,-ncol(df)])))/sd(rowMeans2(df[,-ncol(df)])) # scoring-based assessment of positive regulation of MG activation 
    
    cairo_pdf('output/suppl_mic_activation_boxplots.pdf', width = 6, height = 6)
    ggboxplot(pl_df, x= 'dis', y='mean_score', fill = c("#377EB8","#4DAF4A","#E41A1C"),xlab = '', ylab='Z-score (microglial activation)', 
              add = 'jitter', add.params = list(alpha=0.5), alpha = 0.75, ggtheme = theme_bw()) + stat_compare_means(ref.group = 'Ctrl',method = 't.test') + theme(text = element_text(family = 'Arial'))
    dev.off()
    
   
  ##### 2.8 Identify peaks ####
  rm(list=ls())
  .rs.restartR()
  
  library(SnapATAC)
  x.after.sp <- readRDS("psp_cbd_ctrl_peaks_rev.Rds")
  library(reticulate)
  library(GenomicRanges)
  system("which snaptools")
  system("which macs2")
 
  ### 2.8.2 MACS2 Peak-calling - all clusters
  dir.create('peak_calling')
  set.seed(42)
  # call peaks for all cluster with more than 200 cells
  clusters.sel = as.integer(names(table(x.after.sp@cluster))[which(table(x.after.sp@cluster) > 200)])
  
  peaks.ls <- parallel::mclapply(seq(clusters.sel), function(i){
    print(clusters.sel[i]);
    peaks <<- runMACS(
      obj=x.after.sp[which(x.after.sp@cluster==clusters.sel[i]),], 
      output.prefix=paste0("psp_cbd_ctrl", gsub(" ", "_", clusters.sel)[i]),
      path.to.snaptools="/home/nes/miniconda3/bin/snaptools",  ## modify here to navigate to the specific directory
      path.to.macs="/home/nes/miniconda3/bin/macs2", ## modify here to navigate to the specific directory
      gsize="hs", # mm, hs, etc
      buffer.size=500, 
      num.cores=1, #don't change here, but...
      macs.options="--nomodel --shift -75 --ext 150 --qval 5e-2 -B --SPMR --call-summits",
      tmp.folder=tempdir()
    )
    }, mc.cores=4); # ... here, if required!
  
  # assuming all .narrowPeak files in the current folder are generated from the clusters
  peaks.names = system("ls | grep narrowPeak", intern=TRUE);
  peak.gr.ls = BiocGenerics::lapply(peaks.names, function(x){
    peak.df = read.table(x) 
    GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
  })
  peak.gr = GenomicRanges::reduce(Reduce(c, peak.gr.ls))
  peak.gr
  
  
  ### 2.8.3 Create a peak bed-file
  library(tidyverse)
  peaks.df = as.data.frame(peak.gr)[,1:3] %>% tidyr::separate(seqnames, into = c('remove','seqnames')) %>% .[,-1] ## added due to bug, see: https://github.com/r3fang/SnapATAC/issues/154
  write.table(peaks.df,file = "peaks.combined.bed",append=FALSE,
              quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
              fileEncoding = "")
  
  ### 2.8.4 metadata correlations
  x.after.sp@metaData$cluster <- as.numeric(x.after.sp@metaData$cluster)
  x <- select(x.after.sp@metaData, c(is.numeric, 'case')) %>% select(-c('frip', 'dupli', 'is__cell_barcode'))
  x$case <- as.factor(x$case)
  x <- x %>% dplyr::group_by(case) %>% summarise_all(mean) %>% .[,9:ncol(.)]
  c <- psych::corr.test(x, method = 'pearson')
  p <- ggcorrplot::cor_pmat(x, method = 'pearson')
  cairo_pdf(paste0('output/2.1_metadata_correlations.pdf'),width = 8, height = 8)
  plot(ggcorrplot::ggcorrplot(c$r, p.mat = c$p, type = 'upper', show.diag = T, lab = T, digits = 1,lab_size=2.5,insig = 'blank', legend.title = "Pearson's R") + 
    ggthemes::theme_tufte(base_family = 'helvetica') + theme(axis.text.x = element_text(angle = 30, vjust=0.99, hjust=1)) + labs(x='',y='')+
      scale_x_discrete(labels =c("total peaks", "duplicate" ,"chimeric" , 'lowmapq',"unmapped reads" , "mitochondrial reads", "passed filters", 
                                 "TSS fragments", "DNase sensitive region fragments", "enhancer region fragments",  "promoter region fragments", 
                                 "on_target_fragments", "blacklist region fragments" ,"peak region fragments",   "peak region cutsites", 
                                 "promoter ratio", "Age", "PMI" , "read depth"  )) +
      scale_y_discrete(labels =c("total peaks", "duplicate" ,"chimeric" ,"unmapped reads" , 'lowmapq', "mitochondrial reads", "passed filters", 
                                 "TSS fragments", "DNase sensitive region fragments", "enhancer region fragments",  "promoter region fragments", 
                                 "on_target_fragments", "blacklist region fragments" ,"peak region fragments",   "peak region cutsites", 
                                 "promoter ratio", "Age", "PMI" , "read depth"  ))
    ) 
  dev.off()
  
  
################## finished #####################
print("Part 2 is done: now execute 'bash scripts/02.8_bash_step.sh' in a terminal, please, and continue with '03_snapatac_analysis_p2_202107.R' afterwards.")
sessionInfo()
