## PETCTomics
- A collection of radiomics-oriented projects using medical images such as PET, MRI, or CT images to build clinical relevant predictive models
- Explore the power of medical images to solve clinical, research or other problems using machine learning/deep learning

### Contributors
- Shih-ying (Clare) Huang, Ph.D.
- Roy Harnish, MS.

### Software Setup:
- Python 2.7.11 :: Anaconda custom (64-bit)
- ITK 4.10.1 (ImageIO and 3D texture feature calculation)
- On /data/francgrp1, one would need to run the following to set up the software packages needed for this repository:

```
source .francgrp_cshrc
```

Radiomics Pipeline
--------------------

##### Segment Images for Feature Extraction

- Segment the tumors from the PET images semi-manulally via MeVisLab (by R. Harnish)
- Segment the tumors from the MR images using N. Hylton's lab in-house segmentation software

##### Image Feature Extraction

- GetImageFeatures_MRI.py: extract radiomic features from MR images
- GetImageFeatures_PET.py: extract radiomic features from PET images
- GetImageFeatures_otherImagetype.py: extract radiomic features from other type of images, such as the MAMMI-PET images
- FeatureAnalysis.py: initial exploration of radioimic features and further processing code
- to-do: make it a more general package to read any type of images (MR, PET, or others)


##### Data Exploration

- OptimalCluster.py: cluster the tumors based on its features and determine the optimal setup for consensus clustering
- feature_cluster_analysis.py: univariate feature analysis
- her2_ImageFeatures_ConsensusClustering_PETMRI_v2.R: consensus clustering was implemented in R
- tumor_cluster_outcome_barplot.py: plot proportional table for freq wrt tumor cluster class and freq wrt clinical outcome
- notebooks/Radiomics_PET_MR_correlation.ipynb: evaluate the correlation between the PET and MR radiomics
- Radiomics_PrimVsRecurTumors.py: examine the relationship between the radiomics generated from the primary tumors vs. the recurred tumors via bokeh annular_wedge plot


##### Data-driven Learning


- classification_cv.py: perform classification tasks via nested cross validation where the inner loop optimizes the classifier's parameters and the outerloop optimize the classifier's performance
- GetLearningOutput.py: report classifer's performance and estimate the AUC confidence interval based on assumption (auc_mean - 1.96*auc_std/np.sqrt(# of innerloops))
- GetLearningOutput_CI.py: report classifer's performance and estimate AUC confidence interval via brute force approach (N of classification runs=1000 at least)
- Final_feature_importance.ipynb: determine feature importance of ElasticNet logistic regression classifer
- Learning_outcome_plot.ipynb: display the classification performanc across all classifier algorithms for manuscript
- RFS_analysis_her2.ipynb: perform survival analysis of recurrence free survival from the her2 data (may not have enough No of patients for effective survival analysis)
 

Toolkit Functions
--------------------


- ImageIO
- ImagePreProcess
- DataExplore
- her2
- DB
- aglio_y_olio





