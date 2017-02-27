# -*- coding: utf-8 -*-
"""
Created on Wed May  4 21:37:12 2016

@author: rharnish
"""

#%%

import dicom
import image_geometry
import numpy as np

#%%

class DMI:
    
    def __init__(self,*args):
        
        self._dmi_path       = None
        self._image_geometry = None
        self._DS             = None
        self._PIX            = None
        
        if args:
            
            self._dmi_path = args[0]
            self._load_dmi()
            
    def _load_dmi(self):
        
        self._DS = dicom.read_file(self._dmi_path)
        self._determine_geometry()
        self._load_pixel_data()
        
    def _determine_geometry(self):

        tag_ImagePositionPatient  = dicom.tag.Tag(0x0020,0x0032)
        tag_ImageDirectionCosines = dicom.tag.Tag(0x0119,0x1030)
        tag_PixelSpacing          = dicom.tag.Tag(0x0028,0x0030)
        tag_Rows                  = dicom.tag.Tag(0x0028,0x0010) 
        tag_Cols                  = dicom.tag.Tag(0x0028,0x0011) 
        tag_SliceThickness        = dicom.tag.Tag(0x0018,0x0050) 
        tag_SpacingBetweenSlices  = dicom.tag.Tag(0x0018,0x0088)
        
        # TODO: make this cover every case
        # to determine how many slices
        tag_Planes   = dicom.tag.Tag(0x0028,0x0012)
        Planes       = self._DS[tag_Planes]
        n_slices     = Planes.value
        
        ImagePositionPatient  = self._DS[tag_ImagePositionPatient]
        ImageDirectionCosines = self._DS[tag_ImageDirectionCosines]
        PixelSpacing          = self._DS[tag_PixelSpacing]
        Rows                  = self._DS[tag_Rows]
        Cols                  = self._DS[tag_Cols]        
        SliceThickness        = self._DS[tag_SliceThickness]        
        SpacingBetweenSlices  = self._DS[tag_SpacingBetweenSlices]
        
        # instantiate and get reference to ImageGeometry
        self._image_geometry = image_geometry.ImageGeometry()
        IG = self._image_geometry
        
        # set origin
        IG._origin = ImagePositionPatient.value

        # set direction cosines
        IG._dir_as_col_idx_increases   = [ImageDirectionCosines.value[i]   for i in np.arange(3)]
        IG._dir_as_row_idx_increases   = [ImageDirectionCosines.value[i+3] for i in np.arange(3)]
        IG._dir_as_slice_idx_increases = [ImageDirectionCosines.value[i+6] for i in np.arange(3)]
              
        # sampling and size      
        IG._samplingRCS = [PixelSpacing.value[0],PixelSpacing.value[1],SpacingBetweenSlices.value]
        IG._shapeRCS    = [Rows.value,Cols.value,n_slices]

        # permute RCS --> SRC
        IG._shapeSRC    = [IG._shapeRCS[2]] + IG._shapeRCS[0:2]
        IG._samplingSRC = [IG._samplingRCS[2]] + IG._samplingRCS[0:2]                
          
        # compute derived quantities
        IG._determine_orientation()
        IG._determine_span()
        IG._generate_direction_cosine_matrix()
        
        
    def _load_pixel_data(self):
        
        '''
        forcing NumberOfFrames to be the number of slices in 
        DS since pydicom uses NumberOfFrames when determining
        shape of pixel data
        '''
        self._DS.NumberOfFrames = self._image_geometry.shapeRCS[2]
        self._PIX = self._DS.pixel_array

# import matplotlib.pyplot as plt
#
# test_dmi_path = '/data/francgrp1/breast_radiomics/her2/MRI/001/20051010/E26508/dce/001_E26508S004_sagL_tp1.dmi'
# test_dmi = DMI(test_dmi_path)
# test_dmi._image_geometry.print_self()
#
# PIX   = test_dmi._PIX
# SLICE = PIX[int(PIX.shape[0]/2),:,:]
# plt.imshow(SLICE)
# plt.axis('off')
# plt.show()
#
#
# #%%
#
# def get_midslice(dmi_path):
#     test_dmi = DMI(dmi_path)
#     PIX   = test_dmi._PIX
#     SLICE = PIX[int(PIX.shape[0]/2),:,:]
#     return SLICE


# testing = False
# #testing = True
#
# if testing:
#
#     import os
#     from glob import glob as glob
#     mri_dir = '/data/francgrp1/breast_radiomics/her2/MRI'
#     dce_dirs = glob(mri_dir + '/*/*/*/dce')
#     dce_dirs = sorted(dce_dirs)
#
#     for dce_dir in dce_dirs[:]:
#         print '\n\n' + dce_dir
#         exam_dir    = os.path.dirname(dce_dir)
#         visit_dir   = os.path.dirname(exam_dir)
#         subject_dir = os.path.dirname(visit_dir)
#         subject_id  = os.path.basename(subject_dir)
#         pdf_name    = mri_dir + '/TEST/' + subject_id + '.sagR.PDF'
#         dmis = glob(dce_dir + '/*sagR*dmi')
#         dmis = sorted(dmis)
#         if len(dmis) > 0:
#             fig, ax = plt.subplots(nrows=1, ncols=len(dmis))
#             axidx = 0
#             for dmi in dmis:
#                 SLICE = get_midslice(dmi)
#                 if axidx == 0:
#                     wmin = np.min(SLICE)
#                     wmax = np.max(SLICE) * 1.5
#                 ax[axidx].imshow(SLICE,vmin=wmin,vmax=wmax)
#                 ax[axidx].axis('off')
#                 axidx = axidx + 1
#             plt.savefig(pdf_name,format='pdf')
#             plt.show()
#             del fig, ax
            