# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 21:30:36 2016

@author: rharnish
"""

#%%

import itk
import itkConfig
itkConfig.ProgressCallback = None
import ITKImageHelper

#%%

'''
Gaussian Smoothing
'''

def itkImage_DiscreteGaussianImageFilter(itkImageIn, variance):
    FilterType = itk.DiscreteGaussianImageFilter[type(itkImageIn),itk.Image[itk.F,3]]
    gaussianFilter = FilterType.New()
    gaussianFilter.SetInput( itkImageIn )
    gaussianFilter.SetVariance( variance )
    gaussianFilter.UpdateLargestPossibleRegion()
    outputItkImage = gaussianFilter.GetOutput()
    return outputItkImage
    
#%%
    
'''
Image Subtraction
'''   

def itkImage_Subtract(itkImageIn1,itkImageIn2):
    FilterType = itk.SubtractImageFilter[type(itkImageIn1),type(itkImageIn2),type(itkImageIn1)]
    subtractFilter = FilterType.New()
    subtractFilter.SetInput1( itkImageIn1 )
    subtractFilter.SetInput2( itkImageIn2 )
    subtractFilter.UpdateLargestPossibleRegion()
    outputItkImage = subtractFilter.GetOutput()
    return outputItkImage
    
#%%
    
'''
Difference of Gaussian
'''
    
def itkImage_DifferenceOfGaussian(itkImageIn, var1, var2):
    smooth1 = itkImage_DiscreteGaussianImageFilter(itkImageIn, var1)
    smooth2 = itkImage_DiscreteGaussianImageFilter(itkImageIn, var2)
    itkImageOut = itkImage_Subtract(smooth1,smooth2)
    return itkImageOut
    
#%%
    
def itkImage_GrayScaleFillHoles(itkImageIn):
    FilterType = itk.GrayscaleConnectedClosingImageFilter[type(itkImageIn),type(itkImageIn)]
    gsFilter   = FilterType.New()
    gsFilter.SetInput(itkImageIn)
    gsFilter.UpdateLargestPossibleRegion()
    return gsFilter.GetOutput()
    
#%%
    
def itkImage_Threshold(itkImageIn,low,high):
    FilterType = itk.BinaryThresholdImageFilter[type(itkImageIn),type(itkImageIn)]
    tfilt = FilterType.New()
    tfilt.SetInput(itkImageIn)
    tfilt.SetLowerThreshold(low)
    tfilt.SetUpperThreshold(high)
    tfilt.SetInsideValue(0)
    tfilt.SetOutsideValue(1)
    tfilt.UpdateLargestPossibleRegion()
    return tfilt.GetOutput()
    
def itkImage_MedianFilter(itkImageIn,radius):
    FilterType = itk.MedianImageFilter[type(itkImageIn),type(itkImageIn)]
    tfilt = FilterType.New()
    tfilt.SetInput(itkImageIn)
    tfilt.SetRadius(radius)
    tfilt.UpdateLargestPossibleRegion()
    return tfilt.GetOutput()
    
def itkImage_DistanceMap(itkImageIn):
    FilterType = itk.SignedDanielssonDistanceMapImageFilter[type(itkImageIn),itk.Image[itk.F,3]]
    dmap = FilterType.New()
    dmap.SetInput(itkImageIn)
    dmap.UpdateLargestPossibleRegion()
    return dmap.GetOutput()

def ITKImageResample(input_itk_img,output_img_spacing,is_mask):
    """
    :param input_itk: an ITK format image volume
           If read ITK image with ITKImageHelper.generate_oriented_itkImage()..
           need to make sure to run the output ITK with ITKImageHelper.itkImage_orient_to_[xxxxx] etc.
           or else, the output would be all 0's
    :param output_img_spacing: the voxel spacing in millimeter in row,col,slice direction
    :return: a resampled ITK format image volume
    """

    # get input image info
    input_img_type = type(input_itk_img)
    input_img_spacing = input_itk_img.GetSpacing()
    input_img_size = input_itk_img.GetLargestPossibleRegion().GetSize()
    input_img_dim = len(input_img_size)
    # ITKImageHelper.itkImage_print_image_geometry(input_itk_img)

    # set up a transform (required by ITK ResampleImageFilter
    transform = itk.IdentityTransform[itk.D, input_img_dim].New()
    transform.SetIdentity()

    # interpolation
    if is_mask:
        Interpolator = itk.NearestNeighborInterpolateImageFunction[input_img_type, itk.D].New()
    else:
        Interpolator = itk.LinearInterpolateImageFunction[input_img_type, itk.D].New()

    # bi-cubic interplator didn't work with the image type, UC3
    # Interpolator = itk.BSplineInterpolateImageFunction[input_img_type,itk.D].New()

    # setup output image size
    output_img_size = itk.Size[input_img_dim]()
    # output_img_size[0] = long((input_img_size[0]-1) * input_img_spacing[0] / output_img_spacing[0])
    # output_img_size[1] = long((input_img_size[1] -1) * input_img_spacing[1] / output_img_spacing[1])
    # output_img_size[2] = long((input_img_size[2]-1) * input_img_spacing[2] / output_img_spacing[2])
    output_img_size[0] = long(input_img_size[0]* input_img_spacing[0] / output_img_spacing[0])
    output_img_size[1] = long(input_img_size[1] * input_img_spacing[1] / output_img_spacing[1])
    output_img_size[2] = long(input_img_size[2] * input_img_spacing[2] / output_img_spacing[2])

    resampler = itk.ResampleImageFilter[input_img_type, input_img_type].New()
    resampler.SetInput(input_itk_img)
    resampler.SetTransform(transform)
    resampler.SetInterpolator(Interpolator)
    resampler.SetOutputSpacing(output_img_spacing)
    resampler.SetSize(output_img_size)
    resampler.SetOutputOrigin(input_itk_img.GetOrigin())
    resampler.SetOutputDirection(input_itk_img.GetDirection())
    resampler.UpdateOutputInformation()
    resampler.Update()

    output_itk_img = resampler.GetOutput()
    # ITKImageHelper.itkImage_print_image_geometry(output_itk_img)
    # print 'ImageProcess.ITKImageResample: Set the output ITK image size to {} with voxel spacing {} mm^3'.format(output_itk_img.GetLargestPossibleRegion().GetSize(), output_itk_img.GetSpacing())

    return output_itk_img

def ITKStandardizeImageIntensity(input_itkimg,Nbin):
    img_type = type(input_itkimg)

    minmaxcalcType = itk.MinimumMaximumImageCalculator[img_type]
    minmaxcalcFilter = minmaxcalcType.New()
    minmaxcalcFilter.SetImage(input_itkimg)
    minmaxcalcFilter.Compute()
    img_min = minmaxcalcFilter.GetMinimum()
    img_max = minmaxcalcFilter.GetMaximum()
    print 'ITKStandardizeImageIntensity:: (min,max) = ({},{})'.format(img_min,img_max)

    # rescale image intensity to 0-Nbin based on its instrinsic min and max
    intensitywindowfilterType = itk.IntensityWindowingImageFilter[img_type, img_type]
    intensitywindowFilter = intensitywindowfilterType.New()
    intensitywindowFilter.SetInput(input_itkimg)
    intensitywindowFilter.SetWindowMinimum(img_min)
    intensitywindowFilter.SetWindowMaximum(img_max)
    intensitywindowFilter.SetOutputMinimum(0)
    intensitywindowFilter.SetOutputMaximum(Nbin)
    intensitywindowFilter.Update()
    rescale_itk_img = intensitywindowFilter.GetOutput()

    return rescale_itk_img

def ITKImageMultiplyConstant(input_itkimg,val):
    img_type = type(input_itkimg)

    multiplyconstType = itk.MultiplyImageFilter[img_type,img_type,img_type]
    multiplyconstFilter = multiplyconstType.New()
    multiplyconstFilter.SetInput(input_itkimg)
    multiplyconstFilter.SetConstant(val)
    multiplyconstFilter.Update()
    multiply_const_itk_img = multiplyconstFilter.GetOutput()

    return multiply_const_itk_img