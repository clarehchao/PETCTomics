#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 2/14/18

@author: shuang
@goal:
"""

# NOTE: when addiing RadiomicsToolBox to pycharm project, it ended up importaning another ITKImageHelper
# which doesn't contain the function i need! so BEWARE when adding projects to another in pycharm
# as it will import packages unintentionally
# need to TURN ON plt.switch_backend('Qt4Agg') in ImP.display_overlay_volume(...) when running with MAC OSX
# need to install pyqt4 via conda to make it work ==> left and right key works! w00t!

import ImageProcess as ImP
import ITKImageHelper as iih


rootdir = '/Users/shuang/Desktop/MIBG_pt9/3dslicer'
nrrd_fname1 = '{}/Segmentation-label-heart.nrrd'.format(rootdir)
nrrd_fname2 = '{}/Right_lung_label.nrrd'.format(rootdir)
nrrd_fname3 = '{}/3 CT WB EFOV.nrrd'.format(rootdir)

itk_img = iih.itkImage_read(nrrd_fname3)
itk_img_orient = iih.itkImage_orient_to_axial(itk_img)

mask_itk = iih.itkImage_read(nrrd_fname2)
mask_itk_orient = iih.itkImage_orient_to_axial(mask_itk)

ImP.ITK_Image_OverlayPlot(itk_img_orient,mask_itk_orient,'CT image with organ mask')

