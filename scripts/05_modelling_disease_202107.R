###########################
## TITLE: 05 modelling disease-specific TF dynamics 
## Author: Nils Briel, Feodor-Lynen-Strasse 17, 81377 Munich, Bavaria
## This Rscript establishes a logistic regression model to predict disease state from sc-data and to understand the contributions to this prediction
##############################
library(SnapATAC)
library(caret)
library(GenomicRanges)
library(viridisLite)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(grid)
library(parallel)
library(fastglm)
library(ggrepel) 
library(lime)
set.seed(42)
setwd("..")
##############################
# make directory, if not existing
ifelse(dir.exists(paste0('output/modelling_rev')),
       stop(paste0("A storage directory for this project already exists. Please edit first.")),
       dir.create(paste0('output/modelling_rev')))
# load astrocytic dataset
x.sp <- readRDS("output/Ast_rev/psp_cbd_ast.Rds")  
# define label data
outcome <- as.numeric(as.factor(x.sp@metaData[,c('disease')])) %>% scales::rescale(c(0,1))
# define predictors
mat_t <- cbind(as.matrix(x.sp@mmat),outcome)
# train-test split
idy <- sample(nrow(mat_t), nrow(mat_t)*0.8)
mat_tr <- mat_t[idy,]
mat_te <- mat_t[-idy,]

# Train Xtreme Gradient Boosting Machine
  fitControl <- trainControl(method = "repeatedcv", number=10, repeats=3, classProbs = T)
  outcome <- as.factor(x.sp@metaData[,c('disease')])

  ## xgb Tree
  xgbTree_mod <- list()
  set.seed(825)
  mod <- train(y = outcome[idy], x = mat_tr[,-ncol(mat_tr)], 
               method = "xgbTree", 
               trControl = fitControl, 
               verbose = T, 
               metric = "ROC")
  
  xgbTree_mod[['summa']] <- summary(mod)
  xgbTree_mod[['predy']] <- stats::predict(mod, mat_te[,-ncol(mat_te)], type = 'raw')
  summary(xgbTree_mod[['predy']])
  xgbTree_mod[['confmat']] <- caret::confusionMatrix(xgbTree_mod[['predy']], as.factor(x.sp@metaData[-idy,c('disease')]))
  xgbTree_mod[['ROC_curves']] <- caTools::colAUC(as.numeric(as.factor(xgbTree_mod[['predy']])), as.factor(x.sp@metaData[-idy,c('disease')]), plotROC = TRUE)
  xgbTree_mod[['model']] <- mod


## Assess performance and confusion matrix
  tiff(paste0('output/modelling_rev/5_xgbTreeaccroc_dis_models.tiff'), units = 'px',compression = 'lzw', res = 300, width = 2300, height = 1600) 
  caTools::colAUC(as.numeric(as.factor(xgbTree_mod[['predy']])), as.factor(x.sp@metaData[-idy,c('disease')]), plotROC = TRUE)
  dev.off()

  library(cvms)
  data <- data.frame(
    "target" = as.character(x.sp@metaData[-idy,c('disease')]),
    "prediction" = as.character(xgbTree_mod[['predy']]),
    stringsAsFactors = FALSE
  )
  eval <- evaluate(
    data = data,
    target_col = "target",
    prediction_cols = "prediction",
    type = 'multinomial'
  )
  
  tiff(paste0('output/modelling_rev/5_xgb_confmat.tiff'), units = 'px',compression = 'lzw', res = 300, width = 1400, height = 1400) 
    plot_confusion_matrix(eval[["Confusion Matrix"]][[1]], add_sums = TRUE, add_row_percentages = F, add_col_percentages = F,
                        sums_settings = sum_tile_settings(
                          label = "Total",
                          tc_tile_border_color = "black"
                          ))
    dev.off()
    pdf("output/modelling_rev/5_xgb_confmat.pdf")       # Export PDF
    x <- as.data.frame(eval[1:8]) %>% t()
    x[,1] <- round(x[,1], 3)
    gridExtra::grid.table(x)
    dev.off()
    
    # Create ROC-curves
    library(ROCR)
    x <- stats::predict(xgbTree_mod[['model']], mat_te[,-ncol(mat_te)], type = 'prob')
    y_list <- list(PSP = ifelse(outcome[-idy]=='PSP',1,0),
                   CBD = ifelse(outcome[-idy]=='CBD',1,0),
                   Ctrl = ifelse(outcome[-idy]=='Ctrl',1,0))
    df_list <- list(PSP = data.frame(predictions = x$PSP, labels = y_list$PSP),
                    CBD = data.frame(predictions = x$CBD, labels = y_list$CBD),
                    Ctrl = data.frame(predictions = x$Ctrl, labels = y_list$Ctrl))
    for(i in c('PSP', 'CBD', 'Ctrl')){
    pred <- prediction(df_list[[i]]$predictions, df_list[[i]]$labels)
    perf <- performance(pred,"tpr","fpr")
    df_list[[i]] <- ggplotify::as.ggplot(~plot(perf,colorize=F, main=paste0('ROC: ',i))) + theme_classic()
    }
    pdf("output/modelling_rev/5_xgb_rocs.pdf")
    plot(plot_grid(plotlist = df_list, ncol = 3))
    dev.off()
    

# Use LIME to explain the model
  # Create explainer object
  components_lime <- lime(
    x = as.data.frame(mat_tr),
    model = xgbTree_mod[['model']], 
    n_bins = 10
  )
  summary(components_lime)
  
  psp_ob <- subset(mat_te, outcome = 1)[1,] %>% as.data.frame()
  ctrl_ob <- subset(mat_te, outcome = 0.5)[1,] %>% as.data.frame()
  cbd_ob <- subset(mat_te, outcome = 0)[1,] %>% as.data.frame()
  
  lime_explanation <- lime::explain(
    x = as.data.frame(mat_te), 
    explainer = components_lime, 
    n_permutations = 100,
    dist_fun = "gower",
    bin_continuous = T,
    kernel_width = 0.5,
    n_features = 10, 
    n_labels = 3,
    feature_select = "highest_weights"
  )

best_pred_PSP <- subset(lime_explanation, label %in% 'PSP') %>% .[.$model_prediction == max(.$model_prediction),] %>% .$case %>% .[1]
best_pred_CBD <- subset(lime_explanation, label %in% 'CBD') %>% .[.$model_prediction == max(.$model_prediction),] %>% .$case %>% .[1]
best_pred_Ctrl <- subset(lime_explanation, label %in% 'Ctrl') %>% .[.$model_prediction == max(.$model_prediction),] %>% .$case %>% .[1]

tiff(paste0('output/modelling_rev/5_lime_bestmodelxgb_featimportance.tiff'), units = 'px',compression = 'lzw', res = 450, width = 8000, height = 3600) 
plot_features(lime_explanation, ncol = 3, cases = c(best_pred_CBD, best_pred_Ctrl, best_pred_PSP))
dev.off()

top_TFs_ast <- list()
top_TFs_ast[['top_pred_PSP']] <- subset(lime_explanation, label %in% 'PSP') %>% .[(.$model_prediction > 0.6 & .$feature_weight > 0.1),] %>% .[order(desc(.$feature_weight)),] %>% .[1:5, 'feature']
top_TFs_ast[['top_pred_CBD']] <- subset(lime_explanation, label %in% 'CBD') %>% .[(.$model_prediction > 0.6 & .$feature_weight > 0.1),]  %>% .[order(desc(.$feature_weight)),] %>% .[1:5, 'feature']
top_TFs_ast[['top_pred_Ctrl']] <- subset(lime_explanation, label %in% 'Ctrl') %>% .[(.$model_prediction > 0.6 & .$feature_weight > 0.1),] %>% .[order(desc(.$feature_weight)),] %>% .[1:5, 'feature']

tiff(paste0('output/modelling_rev/5_featureplot.tiff'), units = 'px',compression = 'lzw', res = 450, width = 6000, height = 3600) 
plot(caret::featurePlot(mat_t[,as.vector(unlist(top_TFs_ast))], outcome, "box", jitter = F))
dev.off()

top_TFs_ast <- list()
top_TFs_ast[['top_pred_PSP']] <- subset(lime_explanation, label %in% 'PSP') %>% .[(.$model_prediction > 0.6 & .$feature_weight > 0.1),] %>% .[order(desc(.$feature_weight)),] %>% .[1:10, c('feature','feature_desc')]
top_TFs_ast[['top_pred_CBD']] <- subset(lime_explanation, label %in% 'CBD') %>% .[(.$model_prediction > 0.6 & .$feature_weight > 0.1),]  %>% .[order(desc(.$feature_weight)),] %>% .[1:10, c('feature','feature_desc')]
top_TFs_ast[['top_pred_Ctrl']] <- subset(lime_explanation, label %in% 'Ctrl') %>% .[(.$model_prediction > 0.6 & .$feature_weight > 0.1),] %>% .[order(desc(.$feature_weight)),] %>% .[1:10, c('feature','feature_desc')]

saveRDS(top_TFs_ast, 'output/modelling_rev/top_TF_ast.Rds')
saveRDS(components_lime, 'output/modelling_rev/components_lime.Rds')
saveRDS(lime_explanation, 'output/modelling_rev/lime_explanation.Rds')
saveRDS(xgbTree_mod, 'output/modelling_rev/xgbTree_mod.Rds')

################## finished #####################.    
print("Part 05 is done: The end of disease modelling. Now continue with '06_RTN_202107.Rmd'")
sessionInfo()
