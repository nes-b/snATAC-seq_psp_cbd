################.
## Author: Nils Briel, Feodor-Lynen-Strasse 23, 81377 Munich, Bavaria
## Date: "22/07/2021"
## Main & Supplemental plots (overview) script
################.

setwd("..")
library(SnapATAC)
library(monocle3)
library(cicero)
library(viridisLite)
library(viridis)
library(tidyverse)
library(export)
library(ggthemes)
library(ggpubr)
library(cowplot)
library(grid)
library(parallel)
library(EnhancedVolcano)
library(pheatmap)
library(colorspace)
source('scripts/plotViz2.R')

#####.

x.after.sp <- readRDS("psp_cbd_ctrl_peaks_rev.Rds")
set.seed(42)
idy <- sample(x.after.sp@barcode, 1e4)
x.after.sp <- x.after.sp[which(x.after.sp@barcode %in% idy),]

ct_cols <- c('#fa7878',  '#afdb85',  '#85db8b',  '#1a8239',  '#00b8b1',  '#3efae7',  '#fcc44c',  '#458cff',  '#7370ff',  '#9970ff',  '#e695dc')
clus_cols <- c("#458cff", "#7370ff", "#1a8239" ,"#fa7878", "#afdb85", "#00b8b1", "#fcc44c", "#85db8b", "#3efae7", "#e695dc", "#9970ff")


#_________________________________________________________________________

##########################################################################.
############################## MAIN FIGURES ##############################
##########################################################################.


##### FIGURE 1 #####
#created with vector graph editor inkscape ###


#_________________________________________________________________________


##### FIGURE 2 ##### 
# Fig2a UMAPs 
tiff(paste0('output/Fig2_metadata_cluster_viz_rev20210609.tiff'), units = 'px',compression = 'lzw', res = 550, width = 6000, height = 8000, type = 'cairo')
plot(cowplot::plot_grid(plotlist = list(
  plotViz2(obj=x.after.sp, method="umap", main="", text.add = F, 
                point.size=1,point.color='celltype',legend.add=F, stroke = 0.3)+ scale_fill_manual(values = ct_cols)+ theme_bw(base_family = 'Arial', base_size = 14)+coord_equal()+
         theme(title = element_text(size = 14), legend.position = 'none'),
  
  plotViz2(obj=x.after.sp,feature.mode=round(x.after.sp@metaData$pmi, digits = -1),method="umap", main="Post mortem interval",
           point.size=1,text.add=FALSE,legend.add=TRUE,  hex.bin = T,stroke=0.1,
           scale_color = sequential_hcl(length(levels(as.factor(round(x.after.sp@metaData$pmi, digits = -1)))),rev = T, palette = 'reds'),bin.size = 125)+ 
    coord_equal() + theme_bw()+ theme(title = element_text(size = 14), legend.text = element_text(size=12), legend.position = 'top'),
  
  plotViz2(obj=x.after.sp, method="umap", main="Sex",text.add = F,hex.bin = T,bin.size = 150,stroke=0.1,
           point.size=1,point.color='gender',legend.add=TRUE)+coord_equal()+scale_fill_manual(values = c('#CC1A36','#458cff'))+ theme_bw()+
    theme(title = element_text(size = 14), legend.text = element_text(size=12), legend.position = 'top'),
  
  plotViz2(obj=x.after.sp, feature.mode = round(x.after.sp@metaData$Age, digits = -1), method="umap",  main="Age at death",
           point.size=1,text.add=FALSE, legend.add=TRUE, hex.bin = T,stroke=0.1,
           scale_color = sequential_hcl(length(levels(as.factor(round(x.after.sp@metaData$Age, digits = -1)))),rev = T, palette = 'reds'),bin.size = 125)+ 
    coord_equal() + theme_bw()+ theme(title = element_text(size = 14), legend.text = element_text(size=12), legend.position = 'top'),
  
  plotViz2(obj=x.after.sp,method="umap", main="Braak & Braak stage",
           point.color='Braak&Braak (NFT)',point.size=0.5,text.add=FALSE,stroke=0.1,legend.add=TRUE, hex.bin = T, bin.size = 125)+coord_equal()+ theme_bw()+
    theme(title = element_text(size = 14), legend.text = element_text(size=12), legend.position = 'top')+ 
    scale_fill_manual(values = c('grey', '#458cff', '#CC1A36')),
  
  plotViz2(obj=x.after.sp,method="umap", main="Thal-Phase",hex.bin = T, bin.size = 125,stroke = 0.1,
           point.color='Thal-Phase', point.size=1,text.add=FALSE,legend.add=TRUE)+coord_equal()+ theme_bw()+
    theme(title = element_text(size = 14), legend.text = element_text(size=12), legend.position = 'top')+ 
    scale_fill_manual(values = c('grey', 'grey42','#458cff', '#CC1A36'))
  
), 
ncol = 2, nrow = 3, axis = 'tblr', align = 'hv') )
dev.off()

# Fig2b&c celltypes
x.after.sp <- readRDS("psp_cbd_ctrl_peaks_rev.Rds")
ct_list_psp = list()
meta_psp <- subset(x.after.sp@metaData, disease == 'PSP')

for(i in levels(as.factor(meta_psp$case))){
  ct_prop <- subset(meta_psp, case == i) %>% 
    dplyr::count(.$celltype) %>%     
    mutate(prop = prop.table(n))
  colnames(ct_prop)[1] = 'ct'
  ct_prop$ct <- as.character(ct_prop$ct)
  ct_prop <- ct_prop[order(ct_prop$ct),]
  
  bp <- ggplot(ct_prop, aes(x='', y=prop, fill=ct))+
    geom_bar(width = 2, color = 'black', stat = "identity") + labs(fill = '', x = i, y = '') + scale_fill_manual(values = ct_cols)+ 
    theme_minimal(base_family="Arial", base_size = 14) +  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 30))
  ct_list_psp[[i]] <- bp
}
ct_list_psp[['psp']] <- ggpubr::ggarrange(ct_list_psp[['num102']], ct_list_psp[['num105']], ct_list_psp[['num114']], ct_list_psp[['num115']],
                                          ncol = 4, common.legend = T, legend = 'none')

ct_list_cbd = list()
meta_cbd <- subset(x.after.sp@metaData, disease == 'CBD')

for(i in levels(as.factor(meta_cbd$case))){
  ct_prop <- subset(meta_cbd, case == i) %>% 
    dplyr::count(.$celltype) %>%     
    mutate(prop = prop.table(n))
  colnames(ct_prop)[1] = 'ct'
  ct_prop$ct <- as.character(ct_prop$ct)
  ct_prop <- ct_prop[order(ct_prop$ct),]
  
  bp <- ggplot(ct_prop, aes(x='', y=prop, fill=ct))+
    geom_bar(width = 2, color = 'black',stat = "identity") + labs(fill = '', x = i, y = '') + scale_fill_manual(values = ct_cols)+ 
    theme_minimal(base_family="Arial", base_size = 14) +  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 30))
  ct_list_cbd[[i]] <- bp
}

ct_list_cbd[['num112']] <- ct_list_cbd[['num112']] + scale_fill_manual(values= c('#fa7878',  '#85db8b',  '#1a8239',  '#fcc44c',  '#458cff',  '#7370ff'))
ct_list_cbd[['cbd']] <- ggpubr::ggarrange(ct_list_cbd[['num104']], ct_list_cbd[['num108']], ct_list_cbd[['num112']], ct_list_cbd[['num113']],
                                          ncol = 4, common.legend = T, legend = 'none')

ct_list_ctrl = list()
meta_ctrl <- subset(x.after.sp@metaData, disease == 'Ctrl')

for(i in levels(as.factor(meta_ctrl$case))){
  ct_prop <- subset(meta_ctrl, case == i) %>% 
    dplyr::count(.$celltype) %>%     
    mutate(prop = prop.table(n))
  colnames(ct_prop)[1] = 'ct'
  ct_prop$ct <- as.character(ct_prop$ct)
  ct_prop <- ct_prop[order(ct_prop$ct),]
  
  bp <- ggplot(ct_prop, aes(x='', y=prop, fill=ct))+
    geom_bar(width = 2, color = 'black',stat = "identity") + labs(fill = '', x = i, y = '') + scale_fill_manual(values = ct_cols)+ 
    theme_minimal(base_family="Arial", base_size = 14) +  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 30))
  ct_list_ctrl[[i]] <- bp
}
ct_list_ctrl[['ctrl']] <- ggpubr::ggarrange(ct_list_ctrl[['numC1']], ct_list_ctrl[['numC2']], ct_list_ctrl[['numC3']], ct_list_ctrl[['numC4']], ct_list_ctrl[['numC5']],
                                            ncol = 5, common.legend = T, legend = 'right')
###
ct_final <- ggpubr::ggarrange(ct_list_psp[['psp']], ct_list_cbd[['cbd']], ct_list_ctrl[['ctrl']], labels = c('PSP', 'CBD', 'Ctrl'),
                              ncol = 3, common.legend = T, legend = 'right', align = 'h', widths = c(0.8, 0.8, 1.65))


ct_frq <- readRDS('output/celltype_freq.Rds')
x.psp_p <- ggpubr::ggboxplot(ct_frq[['x.psp']], x = 'celltype', y = 'prop_to_ctrl', fill = 'celltype', 
                             xlab = '', ylab = 'Relative frequency') +
  geom_hline(yintercept = 1, linetype = 2, alpha = 0.5) + ylim(0,4.5) + scale_x_discrete(limits=rev) +
  coord_flip() + 
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = T)+
  ggthemes::theme_few(base_family="Arial", base_size = 14)  + scale_fill_manual(values = ct_cols) 

x.cbd_p <- ggpubr::ggboxplot(ct_frq[['x.cbd']], x = 'celltype', y = 'prop_to_ctrl', fill = 'celltype', 
                             xlab = '', ylab = 'Relative frequency') +
  geom_hline(yintercept = 1, linetype = 2, alpha = 0.5) + ylim(0,4.5) + scale_x_discrete(limits=rev) + 
  coord_flip() + stat_compare_means(label = "p.signif", method = "t.test",
                                    ref.group = ".all.", hide.ns = T) + scale_fill_manual(values = ct_cols) + 
  ggthemes::theme_few(base_family="Arial", base_size = 14) + theme(axis.text.y=element_blank())

tot_ct <- ggplot(ct_frq[["tot_ct"]], aes(x = Freq, y = Celltype, fill = Celltype)) + geom_col(color = 'black') + labs(x = 'Total frequency', y = '', fill = 'Cell type') +
  scale_fill_manual(values = ct_cols) + geom_text(aes(label=Freq),color = 'black',size = 4,hjust = -0.15) + scale_x_continuous(limits = c(0,1.1e4)) + scale_y_discrete(limits=rev) + 
  ggthemes::theme_tufte(base_family="Arial", base_size = 14) + theme(axis.text.y=element_blank())

ct_frqs <- ggpubr::ggarrange(plotlist = list(x.psp_p, x.cbd_p, tot_ct), common.legend = T, ncol = 3, 
                             legend = 'right', labels = c('PSP', 'CBD', ''),label.y = 1.03,label.x=c(0.37,0,0),widths = c(1.4,0.8, 0.6))

### Fig2d included in 02_snapatac_analysis_p1_202107.R - from l 702

cairo_pdf(file = "output/Fig2ct_final_20210518.pdf", width = 12)
plot(ct_final)
dev.off()
cairo_pdf(file = "output/Fig2ct_frqs_20210609.pdf", width = 12, height = 8)
plot(ct_frqs)
dev.off()


#_________________________________________________________________________


##### FIGURE 3 ##### 
## Fig3a gene disease associations
dis_ass <- readRDS('output/disease_ass/dis_ass_table_list_v4.Rds')
g1 <- ggplot(dis_ass$`PSP:Tauopathies`, aes(x=ct, y=gene, fill = ifelse(FDR<0.05, median.diff, 0), group = ct)) + 
  colorspace::scale_fill_continuous_diverging(palette = 'Blue-Red', p1 = 0.5, limits = c(-6,6), na.value="white") +
  geom_tile(color = 'grey22', size = 0.25) + 
  coord_equal() + labs(title= paste0('PSP:Tauopathies'), x = '', y = '',fill = 'GA difference\n (Dis.Ent.-Ctrl)') +
  ggthemes::theme_tufte(base_family="Arial") + scale_y_discrete(limits=rev) +
  theme(legend.position="right", axis.text.x = element_text(angle = 45, hjust = 1),  axis.text.y = element_text(face = 'italic'))

g2 <- ggplot(dis_ass$`CBD:Tauopathies`, aes(x=ct, y=gene, fill = ifelse(FDR<0.05, median.diff, 0), group = ct)) + 
  colorspace::scale_fill_continuous_diverging(palette = 'Blue-Red', p1 = 0.5, limits = c(-6,6), na.value="white") +
  geom_tile(color = 'grey22', size = 0.25) + 
  coord_equal() + labs(title= paste0('CBD:Tauopathies'), x = '', y = '',fill = 'GA difference\n (Dis.Ent.-Ctrl)') +
  ggthemes::theme_tufte(base_family="Arial") + scale_y_discrete(limits=rev) +
  theme(legend.position="right", axis.text.x = element_text(angle = 45, hjust = 1),  axis.text.y = element_text(face = 'italic'))

gda <- ggpubr::ggarrange(plotlist = list(g1,g2), ncol = 2, nrow = 1, 
                           vjust = c(0.5,-0.5), legend = 'bottom', common.legend = T)
cairo_pdf(paste0('output/Fig3_heatmap_gda_20210620.pdf'), width = 10, height = 10)
plot(gda)
dev.off()

# Fig3b included in 02_snapatac_analysis_p1_202107.R - from l 819

# Fig3c 
data_list <- readRDS('degra_data_list.Rds')
plot_list <- list()
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
      geom_text(aes(label=ifelse(pvalue <= 5e-2, signif((pvalue), 2),'')), nudge_y = 0.0, size = 3, fontface=2)+
      ggthemes::theme_tufte(base_family = 'arial', base_size = 12) + 
      theme(legend.position = 'bottom', axis.text.x =element_text(angle = 30, hjust = 1))
  )
}

cairo_pdf(paste0('output/Fig3_heatmap_degra_20210620.pdf'), width = 8, height = 6)
plot(ggarrange(plotlist = plot_list[1:3], ncol = 1, labels = c('CMA', 'UPR', 'UPS'), common.legend = T))
dev.off()

# Fig3d gchromVAR
zdf <- readRDS('output/gchromVAR_z_df_tauopathy.Rds')
gchrom <- ggplot(zdf, aes(
  x = ct,
  y = tr,
  fill = Zscore)) + labs(fill = 'Z-score', y = '', x = '') + 
  scale_fill_continuous_diverging(palette = 'Blue-Red') + coord_equal() + geom_tile() + 
    geom_text(aes(label=signif(adj.pvalue, 3)), nudge_y = 0.2, size = 3, fontface=2) + 
    geom_text(aes(label=signif(gchromVAR_pvalue, 3)), nudge_y = -0.2, size = 2.5, fontface=3)+
    ggthemes::theme_tufte(base_family = 'Helvetica', base_size = 12) + theme(legend.position = 'bottom', axis.text.x =element_text(angle = 30, hjust = 1))

cairo_pdf(paste0('output/Fig3_heatmap_gchrom.pdf'), width = 8, height = 6)
plot(gchrom)
dev.off()


#_________________________________________________________________________


##### FIGURE 4 ##### 
x.sp <- readRDS('output/Ast_rev/psp_cbd_ast.Rds')
dis_cols <- c("#377EB8","#4DAF4A","#E41A1C")

### Fig4a included in 03B_ast_clus_analysis_202107 - from l 80 

# Fig4b-d
cicero_cds <- readRDS('output/Ast_rev_cbd/cicero_cds_ast_cbd.Rds')
x.sp <- readRDS('output/Ast_rev_cbd/snap_cbd_ast_rev.Rds')
p1 <- plot_cells(cicero_cds, reduction_method = 'UMAP',color_cells_by = "Pseudotime",
                 group_cells_by = "cluster",graph_label_size = 5,label_branch_points = F,label_cell_groups = F,
                 cell_size = 1.1) + #geom_point(shape=21, color ='grey21',alpha =1, stroke = 0.01) +
  theme_void(base_family="Arial", base_size = 10) + theme(legend.position="bottom") + coord_equal() + 
  scale_color_viridis_c(option = 'cividis') + labs(color = 'Pseudotime')
cicero_cds@colData$EMX2 <- log2(x.sp@mmat[,colnames(dplyr::select(as.data.frame(x.sp@mmat),contains('EMX2')))]+1)
p2 <- plot_cells(cicero_cds,reduction_method = 'UMAP',color_cells_by = "EMX2",
                 group_cells_by = "cluster",group_label_size = 1.5,label_branch_points = F,label_cell_groups = F,
                 cell_size = 1.1) + theme_void(base_family="Arial", base_size = 10) + theme(legend.position="bottom") + coord_equal()+
  colorspace::scale_color_continuous_diverging(palette = 'Purp', p2 = .3,rev = T) + labs(color = "log2(TFME:EMX2)")

p3 <- plot_cells(cicero_cds, reduction_method = 'UMAP',color_cells_by = "disease",
                 group_cells_by = "cluster",graph_label_size = 5,label_branch_points = F, label_cell_groups = F,
                 cell_size = 1.1) + scale_color_manual(values =   c("#377EB8","#E41A1C")) + #geom_point(shape=21, color ='grey21',alpha =1, stroke = 0.01) +
  theme_void(base_family="Arial", base_size = 10) + theme(legend.position="bottom")+ coord_equal() + labs(fill = "Entity")
  
mid <- plot(cowplot::plot_grid(p2, p1, p3,ncol = 4, nrow = 1))

## tradeseq
library(tradeSeq)
library(colorspace)
## TFs 
# Fig4e-f
sce_tf <- readRDS('output/Ast_rev_cbd/fitGAM_tf.Rds')
ts_tf <- readRDS('output/Ast_rev_cbd/5.2.3_TF_int_tradeseq.Rds')
mat <- x.sp@mmat %>% scales::rescale(to = c(0,1), from = range(., na.rm = TRUE, finite = TRUE)) %>% t()
smoo1 <- plotSmoothers(sce_tf, gene = 'MA0099.2_FOS::JUN', counts = mat, alpha = 1, border = TRUE) + ggtitle(paste0('MA0099.2_FOS::JUN')) + labs(y = 'TFME')+
  theme_base(base_family="Arial", base_size = 12) 
smoo2 <- plotSmoothers(sce_tf, gene = 'MA0841.1_NFE2', counts = mat, alpha = 1, border = TRUE) + ggtitle(paste0('MA0841.1_NFE2')) + labs(y = 'TFME')+
  theme_base(base_family="Arial", base_size = 12)
# Fig4g
source('scripts/04_cicero_analysis_functions.R')
add <- x.sp@mmat %>% as.data.frame() %>% dplyr::select(contains(c('LHX9','SHOX','EMX','HESX','RFX4','IRF','TFEB','CREB'))) %>% colnames()
yhatSmooth <- tradeSeq::predictSmooth(sce_tf, gene =c(ts_tf[[1]],add), nPoints = 256, tidy = F)
heat1 <- plot(ggplotify::as.ggplot(pheatmap::pheatmap(t(scale(t(yhatSmooth))),
                        color = diverging_hcl(n = 256, palette = 'Blue-Red 3',power = 0.5),
                        fontsize_row = 7,
                        cluster_cols = FALSE,show_rownames = T,show_colnames = FALSE, scale = 'row', fontsize = 7)))
ts_1 <- plot_grid(smoo1, smoo2,nrow = 2, axis = 'tblr')
# suppl. figure --> part of Suppl.Fig.10
cairo_pdf(file = "output/supplFig9_hm_tfs.pdf", width = 18, height = 18)
yhatSmooth <- tradeSeq::predictSmooth(sce_tf, gene =ts_tf[[2]], nPoints = 256, tidy = F)
heat_tf_sup <- plot(ggplotify::as.ggplot(pheatmap::pheatmap(t(scale(t(yhatSmooth))),
                            color = diverging_hcl(n = 256, palette = 'Blue-Red 3',power = 0.5), fontsize_row = 7,
                            cluster_cols = FALSE,show_rownames = T,show_colnames = FALSE, scale = 'row', fontsize = 7)))
dev.off()


#_________________________________________________________________________


##### FIGURE 5 ##### 
## left half
x.sp <- readRDS('output/Ast_rev/psp_cbd_ast.Rds')
xgbTree_mod <- readRDS("output/modelling_rev/xgbTree_mod.Rds")
lime_explanation <- readRDS("output/modelling_rev/lime_explanation.Rds")
set.seed(42)
outcome <- as.numeric(as.factor(x.sp@metaData[,c('disease')])) %>% scales::rescale(c(0,1))
# define predictors
mat_t <- cbind(as.matrix(x.sp@mmat),outcome)
# train-test split
idy <- sample(nrow(mat_t), nrow(mat_t)*0.8)

# Suppl. Fig.12A
roc <- caTools::colAUC(as.numeric(as.factor(xgbTree_mod[['predy']])), as.factor(x.sp@metaData[-idy,c('disease')]), plotROC = TRUE)

# Fig5a-b
library(cvms)
data <- data.frame("target" = as.character(x.sp@metaData[-idy,c('disease')]),"prediction" = as.character(xgbTree_mod[['predy']]),
  stringsAsFactors = FALSE)
eval <- evaluate(data = data,target_col = "target",
  prediction_cols = "prediction",type = 'multinomial')

conf_mat <- plot_confusion_matrix(eval[["Confusion Matrix"]][[1]], add_sums = TRUE, add_row_percentages = F, add_col_percentages = F,
                      sums_settings = sum_tile_settings(label = "Total", tc_tile_border_color = "black"
                      )) + theme_tufte(base_family="Arial", base_size = 10)

# Fig5c
cairo_pdf(file = "output/Fig5_addon_suppl.pdf", width = 3, height = 4.5)
plot(caret::varImp(xgbTree_mod[['model']]), top = 25)
dev.off()

# Fig5d-f
best_pred_PSP <- subset(lime_explanation, label %in% 'PSP') %>% .[.$model_prediction == max(.$model_prediction),] %>% .$case %>% .[1]
best_pred_CBD <- subset(lime_explanation, label %in% 'CBD') %>% .[.$model_prediction == max(.$model_prediction),] %>% .$case %>% .[1]
best_pred_Ctrl <- subset(lime_explanation, label %in% 'Ctrl') %>% .[.$model_prediction == max(.$model_prediction),] %>% .$case %>% .[1]
fw_plot <- lime::plot_features(lime_explanation, ncol = 3, cases = c(best_pred_CBD,  best_pred_Ctrl, best_pred_PSP)) + 
  theme_base(base_family="Arial", base_size = 10) 
# Note: only barcodes derived from the targeted label (group entity) and no TFME ranges are depicted in the final figure

## FIGURE 05 -- right half 
# Fig5g (Diagram/flow chart) was created with vector graph editor inkscape
# Fig5g (Network-graph) is included in 06_RTN_202107.Rmd - from l 385 
# Fig4h is included in 06_RTN_202107.Rmd - from l 359  


#_________________________________________________________________________


##### FIGURE 6 ##### 
library(UpSetR)
library(ggplotify)
x <- readRDS('output/upset_list_all.Rds')
y <- readRDS('output/upset_list_psp_cbd.Rds')
x.sp <- readRDS("output/Ast_rev/psp_cbd_ast.Rds")  
int_1 <- Reduce(intersect, x)
int_2 <- Reduce(intersect, list(x[['XGB_model']], x[['trajectory']], x[['Pw.TFME']])) %>% .[!(. %in% int_1)]
int_3 <- Reduce(intersect, list(x[['XGB_model']], x[['RTN']], x[['Pw.TFME']])) %>% .[!(. %in% int_1 | . %in% int_2)]
int_4 <- Reduce(intersect, list(x[['trajectory']], x[['RTN']], x[['Pw.TFME']])) %>% .[!(. %in% int_1 | . %in% int_2 | . %in% int_3)]
xplan <-  as.data.frame(x.sp@mmat[1:2,]) %>% dplyr::select(contains(c(int_1, int_2, int_3, int_4))) %>% colnames() 

tf_mat_snap <- readRDS("TF_mat_comparison_snap.Rds") %>% .[[1]] 
tf_mat_snap <- subset(tf_mat_snap, comparison %in% "PSP vs. CBD") %>% .$TF

int_5 <- Reduce(intersect, y)
int_6 <- Reduce(intersect, list(y[['XGB_CBD']], y[['Pw.TFME_psp.vs.cbd']]))
int_7 <- Reduce(intersect, list(y[['XGB_PSP']], y[['Pw.TFME_psp.vs.cbd']]))
explanatory_c = as.data.frame(x.sp@mmat[1:2,]) %>% dplyr::select(contains(int_6)) %>% colnames()
explanatory_p = as.data.frame(x.sp@mmat[1:2,]) %>% dplyr::select(contains(int_7)) %>% colnames()

# Fig6a
cairo_pdf(file = "output/Fig6_up1.pdf", width = 10, height = 10)
upset(fromList(x), order.by = 'degree', sets.bar.color = 'lightblue')
dev.off()

# Fig6b
source('scripts/plotViz2.R')
x.sp.ctrl <- x.sp[which(x.sp@metaData$disease %in% 'Ctrl'),]
triang1 <- plotTriang(mat = x.sp@mmat, 
                       mat_bg = x.sp.ctrl@mmat, 
                       filter_obj = unique(xplan),
                       facet_x = 'disease', 
                       xlab = '',
                       ylab = 'TFME',
                       title = '', 
                       aggregate_fun = median,
                      order = unique(xplan),
                       save = F,
                       write_table = T)
# Fig6c
cairo_pdf(file = "output/Fig6_up2.pdf", width = 10, height = 10)
upset(fromList(y), order.by = 'degree', sets.bar.color = 'darkorange')
dev.off()

# Fig6d-e
triang2 <- plotTriang(mat = x.sp@mmat, 
                      mat_bg = x.sp.ctrl@mmat, 
                      filter_obj = unique(explanatory_c),
                      facet_x = 'disease', 
                      xlab = '',
                      ylab = 'TFME',
                      title = '', 
                      aggregate_fun = median,
                      order = unique(explanatory_c),
                      save = F,
                      write_table = F)
triang3 <- plotTriang(mat = x.sp@mmat, 
                      mat_bg = x.sp.ctrl@mmat, 
                      filter_obj = unique(explanatory_p),
                      facet_x = 'disease', 
                      xlab = '',
                      ylab = 'TFME',
                      title = '', 
                      aggregate_fun = median,
                      order = unique(explanatory_p),
                      save = F,
                      write_table = F)

tri <- ggarrange(triang3, triang2, common.legend = T, ncol = 1, legend = 'none')
cairo_pdf(file = "output/Fig6_triang.pdf", width = 8, height = 10)
plot(ggarrange(triang1, tri, common.legend = T, ncol = 2, legend = 'right'))
dev.off()

# Fig.6f-g is included in 07_branch_integration_202107.R - from l 168

# Fig.6h was created with vector graph editor inkscape


#_________________________________________________________________________


##########################################################################.
########################## SUPPLEMENTAL FIGURES ##########################
##########################################################################.

#####  Suppl.Fig. 1 included in 02_snapatac_analysis_p1_202107 - from l 562 ##### 
#####  Suppl.Fig. 2 included in 02_snapatac_analysis_p1_202107 - from l 291 and from 941 ##### 
#####  Suppl.Fig. 3 included in 02_snapatac_analysis_p1_202107 - from l 508 ##### 
#####  Suppl.Fig. 4 included in 02_snapatac_analysis_p1_202107 - from l 508 ##### 

#_________________________________________________________________________

#####  Suppl.Fig. 5 ##### 
### GREAT results
### all cell types
plot_ls <- readRDS('output/rgreat_results_all.Rds')

# MF
df <- plot_ls[[1]]
df <- df[order(df$GO.Molecular.Function.Binom_Adjp_BH),]
df <- subset(df, GO.Molecular.Function.Binom_Adjp_BH < 5e-2)
gr1 <- ggplot(df[1:100,], aes(x = GO.Molecular.Function.name, y = id, color = -log(GO.Molecular.Function.Binom_Adjp_BH), fill = GO.Molecular.Function.Binom_Fold_Enrichment)) + 
       scale_color_distiller(palette = 'Greens', direction = 1) + scale_fill_distiller(palette = 'Oranges', direction = 1) +
       geom_tile(size=1, width=0.7, height=0.7) + theme_base(base_size = 12,base_family = 'Arial') + coord_fixed() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
       labs(x='',y='',color= '-log(Adj. P)', fill = 'Binom. fold \nenrichment')

df <- plot_ls[[2]]
df <- df[order(df$GO.Biological.Process.Binom_Adjp_BH),]
df <- subset(df, GO.Biological.Process.Binom_Adjp_BH < 5e-2)
gr2 <- ggplot(df[1:100,], aes(y = id, x = GO.Biological.Process.name, color = -log(GO.Biological.Process.Binom_Adjp_BH), fill = GO.Biological.Process.Binom_Fold_Enrichment)) +
       scale_color_distiller(palette = 'Greens', direction = 1) + scale_fill_distiller(palette = 'Oranges', direction = 1) +
       geom_tile(size=1, width=0.7, height=0.7) + theme_base(base_size = 12,base_family = 'Arial') + coord_fixed() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
       labs(x='',y='',color= '-log(Adj. P)',fill = 'Binom. fold \nenrichment')

df <- plot_ls[[3]]
df <- df[order(df$GO.Cellular.Component.Binom_Adjp_BH),]
df <- subset(df, GO.Cellular.Component.Binom_Adjp_BH < 5e-2)
gr3 <- ggplot(df[1:100,], aes(y = id, x = GO.Cellular.Component.name, color = -log(GO.Cellular.Component.Binom_Adjp_BH), fill = GO.Cellular.Component.Binom_Fold_Enrichment)) + 
       scale_color_distiller(palette = 'Greens', direction = 1) + scale_fill_distiller(palette = 'Oranges', direction = 1) +
       geom_tile(size=1, width=0.7, height=0.7) + theme_base(base_size = 12,base_family = 'Arial') + coord_fixed() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
       labs(x='',y='',color= '-log(Adj. P)',fill = 'Binom. fold \nenrichment')

cairo_pdf(file = "output/supplFig3great.pdf", width = 23, height = 23)
plot(ggpubr::ggarrange(plotlist = list(gr1,gr2,gr3), common.legend = T, ncol = 1, nrow = 3, 
                           legend = 'bottom', align = 'hv'))
dev.off()


#_________________________________________________________________________

#####  Suppl.Fig. 6 included in 03_snapatac_analysis_p2_202107 - from l 173 (A) #####
#                             & 03_snapatac_analysis_p2_202107 - from l 212 (B) 

#_________________________________________________________________________

#####  Suppl.Fig. 7 ##### 
## Volcano plots of Gene-Disease and cell type associations
dis_ass <- readRDS('output/disease_ass/dis_ass_table_list_v4.Rds')
v1 <- EnhancedVolcano(dis_ass$`PSP:Tauopathies`,lab = dis_ass$`PSP:Tauopathies`$ct_gene,
                      x = 'median.diff',y = 'pvalue',#legendVisible = F,
                      title = '',
                      xlab = paste0('Difference of GA-medians: PSP - Ctrl'),ylab = '-log(FDR)', subtitle = '',hline = -log(0.05/nrow(dis_ass$`PSP:Tauopathies`)),
                      FCcutoff = 1, pointSize = 3.0,labSize = 3.4,xlim = c(-6,6),colAlpha = 0.9, 
                      pCutoff = 10e-4, boxedLabels = F) +
  theme_base(base_family="Helvetica", base_size = 12) 
v2 <- EnhancedVolcano(dis_ass$`CBD:Tauopathies`,lab = dis_ass$`CBD:Tauopathies`$ct_gene,
                      x = 'median.diff',y = 'pvalue',#legendVisible = F,
                      title = '',
                      xlab = paste0('Difference of GA-medians: CBD - Ctrl'),ylab = '-log(FDR)', subtitle = '',hline = -log(0.05/nrow(dis_ass$`CBD:Tauopathies`)),
                      FCcutoff = 1, pointSize = 3.0,labSize = 3.5,xlim = c(-6,6),colAlpha = 0.9,
                      pCutoff = 10e-4, boxedLabels = F) +
  theme_base(base_family="Helvetica", base_size = 12) 


### all asts (Ctrl, PSP, CBD)
plot_ls <- readRDS('output/rgreat_results_ast_v3.Rds')
# MF
df <- plot_ls[[1]]
df <- df[order(df$Hyper_Raw_PValue),]
df <- subset(df, Hyper_Raw_PValue < 1e-1)
gr1 <- ggplot(df[1:50,], aes(x = Binom_Fold_Enrichment, 
                             y = reorder(name,-Binom_Fold_Enrichment), 
                             fill = -log(Hyper_Raw_PValue))) + 
  scale_fill_distiller(palette = 'Blues', direction = 1) +
  geom_col(size=0.5, width=0.8, color = 'black') + theme_base(base_size = 12,base_family = 'Arial') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x='-log(Adj. P)',y='', fill = 'Binom. fold \nenrichment')

df <- plot_ls[[2]]
df <- df[order(df$Hyper_Raw_PValue),]
df <- subset(df, Hyper_Raw_PValue < 1e-1)
gr2 <- ggplot(df[1:50,], aes(x = Binom_Fold_Enrichment, 
                             y = reorder(name,-Binom_Fold_Enrichment), 
                             fill = -log(Hyper_Raw_PValue))) + 
  scale_fill_distiller(palette = 'Greens', direction = 1) +
  geom_col(size=0.5, width=0.8, color = 'black') + theme_base(base_size = 12,base_family = 'Arial') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x='-log(Adj. P)',y='', fill = 'Binom. fold \nenrichment')

df <- plot_ls[[3]]
df <- df[order(df$Hyper_Raw_PValue),]
df <- subset(df, Hyper_Raw_PValue < 1e-1)
gr3 <- ggplot(df[1:50,], aes(x = Binom_Fold_Enrichment, 
                             y = reorder(name,-Binom_Fold_Enrichment), 
                             fill = -log(Hyper_Raw_PValue))) + 
  scale_fill_distiller(palette = 'Reds', direction = 1) +
  geom_col(size=0.5, width=0.8, color = 'black') + theme_base(base_size = 12,base_family = 'Arial') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x='-log(Adj. P)',y='', fill = 'Binom. fold \nenrichment')

cairo_pdf(file = "output/supplFig4great_allasts.pdf", width = 23, height = 7.5)
plot(ggpubr::ggarrange(plotlist = list(gr1,gr2,gr3), common.legend = F, ncol = 3, nrow = 1, 
                          legend = 'bottom', align = 'hv', labels = c('Molecular Function','Biological Process','Cellular Component'),
                          label.y = 1))
dev.off()

#_________________________________________________________________________

#####  Suppl.Fig. 8 ##### 

# Part A is included in 03B_ast_clus_analysis_202107.R - from l 256
# Part B
# PSP ast overview-plots
cicero_cds <- readRDS("output/Ast_rev_psp/cicero_cds_ast_psp.Rds")
p1 <- plot_cells(cicero_cds, reduction_method = 'UMAP',color_cells_by = "Pseudotime",
                 group_cells_by = "cluster",graph_label_size = 3,label_branch_points = F,label_cell_groups = F,
                 cell_size = 0.75) + 
  theme_void() + theme(legend.position="bottom") + coord_equal() + scale_color_viridis_c(option = 'cividis') + labs(color = 'Pseudotime')
cicero_cds@colData$EMX2 <- x.sp@mmat[,colnames(dplyr::select(as.data.frame(x.sp@mmat),contains('EMX2')))]
p2 <- plot_cells(cicero_cds,reduction_method = 'UMAP',color_cells_by = "EMX2",
                 group_cells_by = "cluster",group_label_size = 3,label_branch_points = F,label_cell_groups = F,
                 cell_size = 0.75) + 
  theme_void() + theme(legend.position="bottom") + coord_equal() + 
  scale_colour_gradient2(low = "white",mid = 'grey90', high = "blue", limits = c(-0.5,0.7)) + labs(color = "TFME: EMX2")
p3 <- plot_cells(cicero_cds, reduction_method = 'UMAP',color_cells_by = "disease",
                 group_cells_by = "cluster",graph_label_size = 3,label_branch_points = F, label_cell_groups = F,
                 cell_size = 0.75) + 
  scale_color_manual(values = dis_cols) + theme_void() + theme(legend.position="bottom")+ coord_equal()
cairo_pdf('output/supplFig08_pt_sub_overview_part_psp.pdf',  width = 7, height = 18) 
plot(cowplot::plot_grid(p2, p1, p3, ncol = 1, nrow = 3, greedy = T, rel_widths = c(1, 1, 1)))
dev.off()

#_________________________________________________________________________

#####  Suppl.Fig. 9 is included in 03B_ast_clus_analysis_202107.R - from l 186 ####

#_________________________________________________________________________

#####  Suppl.Fig. 10 ##### 
# Part A is included in 04_cicero_analysis_202107.R - from l 332  

# Part B
x.sp <- readRDS('output/Ast_rev_cbd/snap_cbd_ast_rev.Rds')
sce <- readRDS('output/Ast_rev_cbd/fitGAM_ga.Rds')
ts_ga <- readRDS('output/Ast_rev_cbd/5.1.3_GA_int_tradeseq.Rds')
idy <- unique(colnames(x.sp@gmat))
mat <- as.matrix(t(x.sp@gmat[,idy]))
cairo_pdf(file = "output/supplFig10_hm_ga.pdf", width = 18, height = 18)
yhatSmooth <- tradeSeq::predictSmooth(sce, gene =ts_ga[[3]], nPoints = 256, tidy = F)
heat_ga_sup <- plot(ggplotify::as.ggplot(pheatmap::pheatmap(t(scale(t(yhatSmooth))),
                                                            color = diverging_hcl(n = 256, palette = 'Blue-Red',power = 0.75),
                                                            cluster_cols = FALSE,show_rownames = T,show_colnames = FALSE, scale = 'row', fontsize_row = 7)))
dev.off()

bott <- plot_grid(ts_1, heat1,
                  ncol = 2)

final_3 <- plot_grid(top,mid,bott, nrow = 3, rel_heights = c(1, 1, 1))

cairo_pdf(file = "output/Fig4_updated.pdf", width = 18, height = 18)
plot(final_3)
dev.off()

#_________________________________________________________________________

#####  Suppl.Fig. 11 is included in 04_cicero_analysis_202107.R - from l 103 ####

#_________________________________________________________________________

#####  Suppl.Fig. 12 ####
# Part A, see l 338

#_________________________________________________________________________

#####  Suppl.Fig. 13 is included in 06_RTN_202107.R - from l 398 ####

#_________________________________________________________________________

#####  Suppl.Fig. 14 is included in 08_branch_integration_202107.R - from l 142 ####


#################@SessionInfo##################
sessionInfo()
































