rm(list=ls())

library(plyr)
library(car)

the_mri_tp <- 2
the_mri_glcm_binWidth <- 5
the_pet_bin_width <- 0.1
rootdir <- '/Users/shuang/Documents/Proj_Radiomics/Data/her2/her2_Analysis'
data_dir <-  sprintf('%s/PETMRI/PETbinwidth%s_MRItp%s_binwidth%s',rootdir,the_pet_bin_width, the_mri_tp, the_mri_glcm_binWidth)
fname <- sprintf('%s/data_all.csv',data_dir)

df_data <- read.csv(fname)
colnames(df_data)[which(colnames(df_data) == 'Sjoerd_Grade')] <- 'Tumor_Grade'
colnames(df_data)[which(colnames(df_data) == 'Marjan_Histology')] <- 'Tumor_Histology'
# df_data$T_stage_2 <- ifelse(df_data$T_stage == 0 | df_data$T_stage == 1, 1, df_data$T_stage)
# df_data$Overall_stage_2 <- recode(df_data$Overall_stage, 'c(0,1)=1;c(3,4)=3;else=2')
# df_data$Overall_stage_2 <- ifelse(df_data$Overall_stage == 3 | df_data$Overall_stage == 4, 3, df_data$Overall_stage)
# df_data$BCsubtype_2 <- ifelse(df_data$BC_subtype == 0 | df_data$BC_subtype == 1, 1, df_data$BC_subtype)

# tab <- table(df_data$cs_class,df_data$T_stage)
# print(assocstats(tab))

colname_df_data <- colnames(df_data)
# feature_var_names <- c(grep('ShapeSize',colname_df_data, value=TRUE),grep('FOstats',colname_df_data,value=TRUE), grep('texture',colname_df_data, value=TRUE))
feature_var_names <- c(grep('ShapeSize_',colname_df_data, value=TRUE),grep('FOstats_',colname_df_data,value=TRUE), 
                       grep('texture_',colname_df_data, value=TRUE), grep('MRI_',colname_df_data, value=TRUE),
                       grep('PET_',colname_df_data, value=TRUE))

# normalize all feature 
df_data[,feature_var_names] <- scale(df_data[,feature_var_names])

is_order_outcome <- 0 #0 if not ordered outcome, 1 if ordered outcome
if (is_order_outcome == 1){
  outcome_name_vec <- c('Tumor_Grade', 'T_stage', 'N_stage', 'Overall_stage','Overall_stage_2')
  fname_out <- sprintf('%s/assoc_spearmancorr_all_v2.csv',data_dir)
} else{
  # outcome_name_vec <- c('Recurrence_Type', 'Diseasefree_5yr', 'TripleNeg','Tumor_Histology')
  outcome_name_vec <- c('Recurrence','Recurrence_Type','TripleNeg','Tumor_Histology_updated','BoneMetsOrNot','DF_1yr',
                        'OS_1yr','DF_2yr','OS_2yr','DF_3yr','OS_3yr','DF_4yr','OS_4yr','DF_5yr','OS_5yr','BC_subtype')
  fname_out <- sprintf('%s/assoc_corr_all_v2.csv',data_dir)
}

# get correlation coeff between an image feature to an outcome (native values)
df_all <- data.frame()
for (ov in outcome_name_vec){
  print(ov)
  if (is_order_outcome == 0){
    # determine the correlation coefficient using lm
    lm_list <- lapply(feature_var_names, function(x) lm(as.formula(paste(x, '~ factor(', ov, ')')), data = df_data))
    cc_mat <- sapply(lm_list, function(x) summary(x)$r.squared)
  } else{
    # determine the spearman rank correlation (check to see if the relationship is monotonic)
    sprcc_mat <- sapply(feature_var_names, function(x) cor(df_data[,x],df_data[,ov],use='pairwise.complete.obs',method='spearman'))
  }
  
  if (length(unique(df_data[,ov])) == 2){
    print('wilcoxon test!')
    assoc_test_list <- lapply(feature_var_names, function(x) wilcox.test(as.formula(paste(x, '~ factor(', ov, ')')), data = df_data, conf.int = TRUE))
    tmp_mat <- t(sapply(assoc_test_list, function(x) c(x$p.value, x$conf.int[1], x$conf.int[2])))
    if (is_order_outcome == 0){
      tmp_df <- data.frame('feature' = feature_var_names, 'outcome' = ov, 'pval' = tmp_mat[,1], 'CI1' = tmp_mat[,2], 'CI2' = tmp_mat[,3], 'test_type' = c('wilcox'), 'corr_coeff' = cc_mat)
    } else{
      tmp_df <- data.frame('feature' = feature_var_names, 'outcome' = ov, 'pval' = tmp_mat[,1], 'CI1' = tmp_mat[,2], 'CI2' = tmp_mat[,3], 'test_type' = c('wilcox'), 'corr_coeff' = sprcc_mat)
    }
  }
  else{
    print('kruskal test!')
    assoc_test_list <- lapply(feature_var_names, function(x) kruskal.test(as.formula(paste(x, '~ factor(', ov, ')')), data = df_data))
    tmp_mat <- t(sapply(assoc_test_list, function(x) c(x$p.value, x$parameter[[1]], x$statistic[[1]])))
    
    if (is_order_outcome == 0){
      tmp_df <- data.frame('feature' = feature_var_names, 'outcome' = ov, 'pval' = tmp_mat[,1], 'df' = tmp_mat[,2], 'chi2_stats' = tmp_mat[,3], 'test_type' = c('kruskal'), 'corr_coeff'=cc_mat)
    } else{
      tmp_df <- data.frame('feature' = feature_var_names, 'outcome' = ov, 'pval' = tmp_mat[,1], 'df' = tmp_mat[,2], 'chi2_stats' = tmp_mat[,3], 'test_type' = c('kruskal'), 'corr_coeff'=sprcc_mat)
    }
  }
  df_all <- rbind.fill(df_all, tmp_df)
}

write.table(df_all, file = fname_out, sep = ',', row.names = FALSE)





