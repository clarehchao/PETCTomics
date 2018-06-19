## PETCTomics
- A collection of radiomics-oriented projects using medical images such as PET, MRI, or CT images to build clinical relevant predictive models
- Explore the power of medical images to solve clinical, research or other problems using machine learning/deep learning

### Contributors
- Shih-ying (Clare) Huang, Ph.D.
- Roy Harnish, MS.
- Benjamin Franc, M.D., M.S.

### Software Setup:
- Python 2.7.11 :: Anaconda custom (64-bit)
- ITK 4.10.1 (ImageIO and 3D texture feature calculation)
- On Terra server of PRL, use the below anaconda environment to run radiomics-related code

```
source activate py27-pydicom
```

- For any classification related code using scikit-learn package, run the code with the latest python distribution version 3.5

```
source activate py35
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

- her2_ImageFeatures_ConsensusClustering_PETMRI_v2.R: perform consensus clustering (implemented in R)
- OptimalCluster.py
    - cluster the tumors based on its features and determine the optimal setup for consensus clustering
    - create cluster map plots (Radiomics manuscript Figure 1a-4a) 
- tumor_cluster_outcome_barplot.py: plot proportional table for freq wrt tumor cluster class and freq wrt clinical outcome (Radiomics manuscript Figure 1b - 4b)
- feature_cluster_analysis.py: univariate feature analysis (Radiomics manuscript Figure 5)
- notebooks/Radiomics_PET_MR_correlation.ipynb: evaluate the correlation between the PET and MR radiomics
- GetImageFeatures_PET_RecurredTumors.py: computed radiomics features of the recurrent tumors segmented from the PET images (found 5 - 10 years later than the primary DX)
- Radiomics_PrimVsRecurTumors.py
    - examine the relationship between the radiomics generated from the primary tumors vs. the recurred tumors
    - generated via bokeh annular_wedge plot (Radiomics manuscript revision figure R1)
    - generated a facet grid line plot via seaborn of the Pearson correlation coeff of PET radiomics between the primary and recurrent tumors (Radiomics manuscript revision figure R2)

##### Data-driven Learning
- classification_cv.py: perform classification tasks via nested cross validation where the inner loop optimizes the classifier's parameters and the outerloop optimize the classifier's performance
- GetLearningOutput.py: report classifer's performance and estimate the AUC confidence interval based on assumption (auc_mean - 1.96*auc_std/np.sqrt(# of innerloops))
- GetLearningOutput_CI.py
    - report classifer's performance and estimate AUC confidence interval via brute force approach (N of classification runs=1000 at least)
    - generated a heatmap of classifier CV performance (Radiomics manuscript Figure 6)
- Feature_Importance_CV.py
    - determine feature importance of logistic regression classifer (either ElasticNet or any LogReg)
    - print out top-10 ranked features at predicting an outcome (Radiomics manuscript Table 3)
- RFS_analysis_her2.ipynb: perform survival analysis of recurrence free survival from the her2 data (may not have enough No of patients for effective survival analysis)


Toolkit Functions
--------------------
- ImageIO
    - dicom_series.py: read dicom images and determine the number of image series
    - dmi.py: a class to read in .dmi file and create image geometry info
    - image_geometry.py: determine image geometry info from a dicom series or an instance of ITKImage
    - ITKImagHelper.py: a collection of ITK-based image operation helper functions

- ImagePreProcess
    - GLCMTextureFeature.py: compute 3D GLCM image features
    - ImageFeature.py: a class that computes all radiomic features (1st order stats, shape & size, and GLCM features) for a given input image and mask image
    - ImageProcess.py: a collection of image processing functions (e.g. segmentation, image display, etc.)
    - ITKImageFilters.py: a collection of ITK-based image filters
- DataExplore
    - DataHelper.py: a collection of helper function to categorize and clean-up radiomic-related data
    - LearningTool.py: machine learning related helper function including nested cross validation and feature importance
    - StatsTool.py: a collection of statistics-based helper function for data exploration
    - VizPlot.py: a collection of data visualization help functions
- her2
- DB
- aglio_y_olio
    - pet_util.py: determine the conversion factor for converting raw PET images to SUV values





