rm(list=ls())

library(plyr)
library(car)

the_mri_tp <- 2
the_mri_glcm_binWidth <- 5
the_pet_bin_width <- 0.1
rootdir <- '/Users/shuang/Documents/Proj_Radiomics/Data/her2/her2_Analysis'
data_dir <-  sprintf('%s/PETMRI/PETbinwidth%s_MRItp%s_binwidth%s',rootdir,the_pet_bin_width, the_mri_tp, the_mri_glcm_binWidth)
fname <- sprintf('%s/data_oncotype.csv',data_dir)

df_data <- read.csv(fname)

colname_df_data <- colnames(df_data)
print(colname_df_data)
# feature_var_names <- c(grep('ShapeSize',colname_df_data, value=TRUE),grep('FOstats',colname_df_data,value=TRUE), grep('texture',colname_df_data, value=TRUE))
feature_var_names <- c(grep('ShapeSize_',colname_df_data, value=TRUE),grep('FOstats_',colname_df_data,value=TRUE), 
                       grep('texture_',colname_df_data, value=TRUE), grep('MRI_',colname_df_data, value=TRUE),
                       grep('PET_',colname_df_data, value=TRUE))

# normalize all feature 
df_data[,feature_var_names] <- scale(df_data[,feature_var_names])

outcome_name_vec <- c('ODx_score')
fname_out <- sprintf('%s/assoc_spearmancorr_oncotype.csv',data_dir)


# get correlation coeff between an image feature to an outcome (native values)
df_all <- data.frame()
for (ov in outcome_name_vec){
  print(ov)
  # determine the spearman rank correlation (check to see if the relationship is monotonic)
  sprcc_mat <- sapply(feature_var_names, function(x) cor(df_data[,x],df_data[,ov],use='pairwise.complete.obs',method='spearman'))
  
  if (length(unique(df_data[,ov])) == 2){
    print('wilcoxon test!')
    assoc_test_list <- lapply(feature_var_names, function(x) wilcox.test(as.formula(paste(x, '~ factor(', ov, ')')), data = df_data, conf.int = TRUE))
    tmp_mat <- t(sapply(assoc_test_list, function(x) c(x$p.value, x$conf.int[1], x$conf.int[2])))
    tmp_df <- data.frame('feature' = feature_var_names, 'outcome' = ov, 'pval' = tmp_mat[,1], 'CI1' = tmp_mat[,2], 'CI2' = tmp_mat[,3], 'test_type' = c('wilcox'), 'corr_coeff' = sprcc_mat)
  }
  else{
    print('kruskal test!')
    assoc_test_list <- lapply(feature_var_names, function(x) kruskal.test(as.formula(paste(x, '~ factor(', ov, ')')), data = df_data))
    tmp_mat <- t(sapply(assoc_test_list, function(x) c(x$p.value, x$parameter[[1]], x$statistic[[1]])))
    tmp_df <- data.frame('feature' = feature_var_names, 'outcome' = ov, 'pval' = tmp_mat[,1], 'df' = tmp_mat[,2], 'chi2_stats' = tmp_mat[,3], 'test_type' = c('kruskal'), 'corr_coeff'=sprcc_mat)
  }
  df_all <- rbind.fill(df_all, tmp_df)
}

write.table(df_all, file = fname_out, sep = ',', row.names = FALSE)





