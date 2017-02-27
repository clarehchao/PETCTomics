# clear all variables in the workspace
rm(list=ls())

library(ConsensusClusterPlus)

# # PET image feature data
# rootdir <- '/Users/shuang/Documents/Proj_Radiomics/Data/her2/her2_Analysis/PET/IsoVoxel_IMGBIN128_GLCMBIN64'
# fname <- sprintf('%s/PETdataAll_glcmNbin64_normNbin128.csv',rootdir)
# data <- read.csv(fname,colClasses = c('Anon_Accession_Num'='character','MRN'='character'))
# sub_df <- subset(data, select = -c(breast_side,pet_series_fdir,X,voxel_size_mm3,Anon_Accession_Num,Laterality,PRIMARY_ID,MRN,CR_AccessionSeq,pt_id,TripleNeg,Marjan_Histology))

# MRI image feature data
theTP <-  c(1,2,3)
theglcmBin <-  c(64,128,256)
combo <- do.call(expand.grid,list(tp = theTP, glcmbin = theglcmBin))
for (t in 1:dim(combo)[1]) {
  rootdir <- sprintf('/Users/shuang/Documents/Proj_Radiomics/Data/her2/her2_Analysis/MRI/IsoVoxel_TP%d_GLCMBIN%d',combo$tp[t], combo$glcmbin[t])
  print(rootdir)
  fname <- sprintf('%s/MRIdataAll_tp%d_Nbin%d.csv',rootdir,combo$tp[t], combo$glcmbin[t])
  data <- read.csv(fname,colClasses = c('MRN'='character'))
  sub_df <- subset(data, select = -c(X,voxel_size_mm3,Laterality,PRIMARY_ID,MRN,CR_AccessionSeq,TripleNeg,Marjan_Histology))
  row.names(sub_df) <- sub_df$ptid_side
  
  the_data_mat <- data.matrix(sub_df)
  
  # compute z-score
  # http://stats.stackexchange.com/questions/89809/is-it-important-to-scale-data-before-clustering
  # should scale by z-score before running clustering since the feature values are vastly different and may affect it's clustered
  zscore <- scale(the_data_mat)
  
  # select 80% item resampling (pItem), 100% feature resampling (pFeature)
  # maximum evaluated k of 6 ==> cluster counts of 2,3,4,5,6 are evaluated (maxK)
  # 50 resamplings (reps), agglomerative heirarchical clustering algorithm upon 1- euclidean distance (as distance)
  the_input <-  t(zscore)
  
  # go through several cluster options
  # results <- ConsensusClusterPlus(the_input,maxK = 6, reps = 50, pItem = 0.8, pFeature = 1.0, title = title, innerLinkage='complete',finalLinkage='complete',clusterAlg = 'hc', distance = 'pearson', plot = 'pdf')
  
  cluster_linkage_list <- list('hc' = c('complete','average','mcquitty','median','centroid'),'km' = NULL,
                               'kmdist' = NULL,'pam' = NULL)
  distance_vector <- c('pearson','spearman','euclidean')
  theKmax <- 5
  theReps <- 50
  thePitem <- 0.8
  thePfeature <- 1.0
  
  # conclusion: clustering methods that works well:
  # complete option with HC using pearson or spearman for the method
  # km, kmdist, pam with pearson as method
  # will keep trying systematically!
  
  alg_list <- names(cluster_linkage_list)
  for (i in 1:length(alg_list)) {
    clusteralg_name <- alg_list[i]
    print(clusteralg_name)
    for (j in 1:length(distance_vector)) {
      print(distance_vector[j])
      if (is.null(cluster_linkage_list[[clusteralg_name]])) {
        cc_dir <- sprintf('%s/ConsensusCluster_%s_%s',rootdir,clusteralg_name,distance_vector[j])
        print(sprintf('cc_dir %s', cc_dir))
        results <- ConsensusClusterPlus(the_input,maxK = theKmax, reps = theReps, pItem = thePitem, pFeature = thePfeature,
                                        title = cc_dir,clusterAlg = clusteralg_name, distance = distance_vector[j], plot = 'pdf')
        icl <-  calcICL(results, title = cc_dir, plot = 'pdf')
        cc_fname <- sprintf('%s/ClusterConsensus.csv',cc_dir)
        ic_fname <- sprintf('%s/ItemConsensus.csv', cc_dir)
        write.csv(data.frame(icl[['clusterConsensus']]),cc_fname)
        write.csv(data.frame(icl[['itemConsensus']]),ic_fname)
        
        for (m in 2:theKmax) {
          cs_mat <- results[[m]][['consensusMatrix']]
          cs_class <- results[[m]][['consensusClass']]
          cs_mat_fname <- sprintf('%s/ConsensusMatrix_kmax%d.csv',cc_dir,m)
          cs_df <- data.frame(cs_mat,row.names = names(cs_class))
          colnames(cs_df) <- names(cs_class)
          write.csv(cs_df,cs_mat_fname)
          cs_class_fname <- sprintf('%s/ConsensusClass_kmax%d.csv',cc_dir,m)
          write.csv(data.frame(cs_class),cs_class_fname)
        }
      }
      else {
        for (k in 1:length(cluster_linkage_list[[clusteralg_name]])) {
          cc_dir <- sprintf('%s/ConsensusCluster_%s_%s_%s',rootdir,clusteralg_name,cluster_linkage_list[[clusteralg_name]][k],distance_vector[j])
          results <- ConsensusClusterPlus(the_input,maxK = theKmax, reps = theReps, pItem = thePitem, pFeature = thePfeature,
                                          title = cc_dir,clusterAlg = clusteralg_name, innerLinkage = cluster_linkage_list[[clusteralg_name]][k], 
                                          finalLinkage = cluster_linkage_list[[clusteralg_name]][k], distance = distance_vector[j], plot = 'pdf')
          icl <-  calcICL(results, title = cc_dir, plot = 'pdf')
          cc_fname <- sprintf('%s/ClusterConsensus.csv',cc_dir)
          ic_fname <- sprintf('%s/ItemConsensus.csv',cc_dir)
          write.csv(data.frame(icl[['clusterConsensus']]),cc_fname)
          write.csv(data.frame(icl[['itemConsensus']]),ic_fname)
          for (m in 2:theKmax) {
            cs_mat <- results[[m]][['consensusMatrix']]
            cs_class <- results[[m]][['consensusClass']]
            cs_mat_fname <- sprintf('%s/ConsensusMatrix_kmax%d.csv',cc_dir,m)
            cs_df <- data.frame(cs_mat,row.names = names(cs_class))
            colnames(cs_df) <- names(cs_class)
            write.csv(cs_df,cs_mat_fname)
            cs_class_fname <- sprintf('%s/ConsensusClass_kmax%d.csv',cc_dir,m)
            write.csv(data.frame(cs_class),cs_class_fname)
          }
        }
      }
    }
  }
}


# row.names(data) <- data$ptid_side
# print(data[names(cs_class),'TripleNeg'])
# print(cs_class)

# # test out cluster consensus equation in Monti et al., 2003
# library('Matrix')
# for (i in 1:kmax) {
#   ki <- (cs_class == i)
#   mask <- kronecker(ki,t(ki)) & lower.tri(matrix(T,length(ki),length(ki)))
#   cc_ki <- sum(cs_mat[mask]) / (sum(ki)*(sum(ki)-1)/2)
#   print(sprintf('Kmax = %d, cc_k%d: %f',kmax,i,cc_ki))
# }
# 
# print(icl[['clusterConsensus']])
# 
# # compute item consensus
# items_list <- names(cs_class)
# j_vector <- 1:length(items_list)
# for (i in 1:kmax) {
#   ki <- (cs_class == i)
#   for (j in j_vector){
#    tmp <- cs_mat[j,]
#    mask <- ki & (j_vector != j)
#    ic_i_j <- sum(tmp[mask])/(sum(ki)-ifelse(ki[j]==i,1,0))
#    print(sprintf('cluster k=%d, item %s, ic_i_j: %f',i,items_list[j],ic_i_j))
#   }
# }
# 
# print(icl[['itemConsensus']])