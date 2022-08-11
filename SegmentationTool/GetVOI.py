#! /usr/bin/env python
#  -*- coding: utf-8 -*-
"""
Created on 4/28/16 1:14 PM

@author: shuang Shih-ying Huang

"""
import matplotlib.pyplot as plt
import sys
import Annotate as ant
import dicom_series
import series_classifier
import ITKImageHelper
import itk
import AxesSequence as axs
import glob
import pandas as pd

# need to decompress the image data....
IM_ROOT_DIR = '/Users/shuang/Documents/Proj_Radiomics/IM/breast_radiomics/her2/PETCT/'
all_series = glob.glob('{}/*/*'.format(IM_ROOT_DIR))

VOIdata = pd.DataFrame(columns=['DataDir','xyz0','nxyz'])

# see this for the backend change so key pressing works http://stackoverflow.com/questions/32129630/key-press-event-not-registered-properly-by-matplotlib
# only 'Qt4Agg' backend works with MacbookPro keyboard, may be different when running on the linux machine
plt.switch_backend('Qt4Agg')

for ii in range(len(all_series)):
    # Set up the image series directory => can put this
    seriesDir = all_series[ii]

    # # get CT and PET series list from the directory
    seriesList = dicom_series.read_files(seriesDir)
    ctSeries, petSeries = series_classifier.get_ct_and_pet_series_from_series_list(seriesList)
    bestPET = series_classifier.determine_best_pet_series(petSeries)
    print bestPET.filenames[0]
    if bestPET is None:
        print 'No suitable PET series was found. Ending program'
        sys.exit()

    print '___________________________________________________________________\n'
    print 'selected PET series is: ', bestPET.info.SeriesNumber
    print '___________________________________________________________________\n'


    # import the PET images into ITK format and orient in Axial direction
    pet_itkImage, pet_itkImageType = ITKImageHelper.generate_oriented_itkImage(bestPET)
    axial_pet_itkImage = ITKImageHelper.itkImage_orient_to_axial(pet_itkImage, pet_itkImageType)

    pet_connector = itk.PyBuffer[pet_itkImageType]
    axialPETIM = pet_connector.GetArrayFromImage(axial_pet_itkImage)



    axes = axs.AxesSequence()
    for i,ax in zip(range(axialPETIM.shape[0]),axes):
        ax.imshow(axialPETIM[i,:,:],cmap='hot',interpolation='bicubic')
        ax.set_title('zlice: {}'.format(i))
    axes.show()
    zslices = axes.selected_slices
    del axes


    Nslice = zslices[1]-zslices[0]+1
    center_zslice = int(zslices[0] + 0.5*Nslice)
    theIM = axialPETIM[center_zslice,:,:]
    plt.imshow(theIM,cmap='hot',interpolation='bicubic')
    a = ant.Annotate()
    plt.show()

    # plot the ROI image
    rx0,ry0,rx1,ry1=[int(val) for val in [a.x0,a.y0,a.x1,a.y1]]
    theROI = theIM[ry0:ry1,rx0:rx1]
    plt.imshow(theROI,cmap='hot',interpolation='bicubic')
    plt.show()
    print rx0,ry0,rx1,ry1

    # double check if the VOI covers the entire tumor
    for n in range(zslices[0],zslices[1],1):
        theROI = axialPETIM[n,ry0:ry1, rx0:rx1]
        plt.imshow(theROI,cmap='hot',interpolation='bicubic')
        plt.show()

    # record the VOI index information
    xyz0 = [rx0,ry0,zslices[0]]
    nxyz = [rx1-rx0+1,ry1-ry0+1,Nslice]
    VOIdata.loc[ii] = [seriesDir,xyz0,nxyz]

# save the VOI data to .csv
VOIdata_fname = '{}/VOIdata.csv'.format(IM_ROOT_DIR)
VOIdata.to_csv(VOIdata_fname)













