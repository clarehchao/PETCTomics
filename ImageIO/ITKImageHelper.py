# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 13:36:28 2015

@author: rharnish
"""

import itk
import itkConfig
itkConfig.ProgressCallback = None
  
#%%

from image_geometry import ImageGeometry
import dicom_series
import dmi
import numpy as np


#%%

import matplotlib.pyplot as plt

def itkImage_show_slice(itkImage,window=[0,1],title=''):
    IG = itkImage_get_image_geometry(itkImage)
    connector = itk.PyBuffer[type(itkImage)]
    PIX  = connector.GetArrayFromImage(itkImage)
    low  = np.min(PIX)
    high = np.max(PIX)
    win  = [float(low) + (float(high)-low)*window[0], 
            float(low) + (float(high)-low)*window[1]]
    plt.imshow(PIX[PIX.shape[0]/2,:,:],
               vmin=win[0],vmax=win[1],
               aspect=(IG.samplingSRC[2]/IG.samplingSRC[1]))
    plt.title(title)
    plt.show()
    
#%%

'''
I/O
'''

def imageIO_get_itk_ctype(imageIO):
    
    cType       = imageIO.GetComponentType()
    cTypeString = imageIO.GetComponentTypeAsString(cType)
    print 'imageIO dtype: ', cTypeString
    
    cType_dict = {}
    cType_dict['unsigned_char']  = itk.UC
    cType_dict['char']           = itk.SC
    cType_dict['short']          = itk.SS
    cType_dict['unsigned_short'] = itk.US
    cType_dict['float']          = itk.F
    cType_dict['double']         = itk.D

    # # TODO: not sure why the below datatype doesn't work with itk.Image[...type,dim]
    # cType_dict['unsigned_int'] = itk.UI
    # cType_dict['unsigned_long'] = itk.UL
    # cType_dict['unsigned_long_long'] = itk.ULL

    
    if cTypeString in cType_dict.keys():
        return cType_dict[cTypeString]
    else:
        print 'NEED TO ADD MORE CTYPES!!!!'


def itkImage_write(itkImage,output_path):
    writerType = itk.ImageFileWriter[ type(itkImage) ]
    writer = writerType.New()
    writer.SetInput( itkImage )
    writer.SetFileName( output_path )
    writer.Update()
    
    
def itkImage_read(input_path):    
    imageIO = itk.ImageIOFactory.CreateImageIO(input_path,itk.ImageIOFactory.ReadMode)
    imageIO.SetFileName(input_path)
    imageIO.ReadImageInformation()
    cType           = imageIO_get_itk_ctype(imageIO)
    nDims           = imageIO.GetNumberOfDimensions()
    itkImageType    = itk.Image[ cType, nDims]
    readerType      = itk.ImageFileReader[ itkImageType ]
    reader          = readerType.New()
    reader.SetFileName(input_path)
    reader.Update()
    itkImage        = reader.GetOutput()
    return itkImage

#%%

def itkImage_orient_to(itkImage, output_orientation = 'axial'):

    itkImageType = type(itkImage)
    input_orientation = itkImage_determine_orientation(itkImage)    
    
    orienterType = itk.OrientImageFilter[ itkImageType, itkImageType ]
    orienter     = orienterType.New()
    orienter.UseImageDirectionOn()
     
    if input_orientation == output_orientation:
        orientations = ['axial','coronal','sagittal']
        print 'same orient'
        idx     = orientations.index(input_orientation)
        new_idx = (idx + 1) % 3
        intermediate_orientation = orientations[new_idx]
        intermediate_orienter    = orienterType.New()
        print 'intermediate_orientation: {}'.format(intermediate_orientation)
        intermediate_orienter.UseImageDirectionOn()
        print 'after useImageDirectionOn'
        intermediate_orienter.SetInput( itkImage )
        print 'after SetupInput'
        
        if intermediate_orientation == 'axial':
            intermediate_orienter.SetDesiredCoordinateOrientationToAxial()
        if intermediate_orientation == 'coronal':
            intermediate_orienter.SetDesiredCoordinateOrientationToCoronal()
        if intermediate_orientation == 'sagittal':
            intermediate_orienter.SetDesiredCoordinateOrientationToSagittal()
        print 'intermediate_orienter.SetDesiredCoordinateOrientationTo...'
            
        intermediate_orienter.Update()
        print 'after IO update'
        orienter.SetInput( intermediate_orienter.GetOutput() )
    else:
        orienter.SetInput( itkImage )
        
    if output_orientation == 'axial':
        orienter.SetDesiredCoordinateOrientationToAxial()
    if output_orientation == 'coronal':
        orienter.SetDesiredCoordinateOrientationToCoronal()
    if output_orientation == 'sagittal':
        orienter.SetDesiredCoordinateOrientationToSagittal()
        
    orienter.Update()
    return orienter.GetOutput()

def itkImage_orient_to_axial(itkImage):
    return itkImage_orient_to(itkImage,output_orientation = 'axial')
    
def itkImage_orient_to_coronal(itkImage):
    return itkImage_orient_to(itkImage,output_orientation = 'coronal')
                              
def itkImage_orient_to_sagittal(itkImage):
    return itkImage_orient_to(itkImage,output_orientation = 'sagittal')
  
# TODO: probably many things... maybe explcit interpolator  
def itkImage_resample_to_template(sourceItkImage, templateItkImage):

    sourceItkImageType = type(sourceItkImage)    
    
    resampleFilterType = itk.ResampleImageFilter[sourceItkImageType, sourceItkImageType]
    resampleFilter = resampleFilterType.New()
    resampleFilter.UseReferenceImageOn()
    resampleFilter.SetReferenceImage(templateItkImage)
    resampleFilter.SetDefaultPixelValue(0)
 
    # determine gross orientation of src and template images
    src_orientation      = itkImage_determine_orientation(sourceItkImage)  
    template_orientation = itkImage_determine_orientation(templateItkImage) 
         
    if src_orientation == template_orientation:
        print 'source and template have same gross orientation'
        
        orienterType = itk.OrientImageFilter[ sourceItkImageType, sourceItkImageType ]
        orienter     = orienterType.New()
        orienter.UseImageDirectionOn()      
        
        orientations = ['axial','coronal','sagittal']
        idx     = orientations.index(src_orientation)
        new_idx = (idx + 1) % 3
        intermediate_orientation = orientations[new_idx]
        print 'intermediate_orientation: {}'.format(intermediate_orientation)
                
        intermediate_orienter    = orienterType.New()
        intermediate_orienter.UseImageDirectionOn()
        intermediate_orienter.SetInput( sourceItkImage )

        if intermediate_orientation == 'axial':
            intermediate_orienter.SetDesiredCoordinateOrientationToAxial()
        if intermediate_orientation == 'coronal':
            intermediate_orienter.SetDesiredCoordinateOrientationToCoronal()
        if intermediate_orientation == 'sagittal':
            intermediate_orienter.SetDesiredCoordinateOrientationToSagittal()
            
        intermediate_orienter.Update()
        resampleFilter.SetInput( intermediate_orienter.GetOutput() )
    else:
        resampleFilter.SetInput(sourceItkImage)
        
    resampleFilter.UpdateLargestPossibleRegion()
    outputImage = resampleFilter.GetOutput()
    return outputImage  
    

def itkImage_extract_region(itkImage, index, size):
    
    itkImageType = type(itkImage)

    #TODO add origin shift since this will keep origin the
    #same in output as original   

    largestPossibleRegion = itkImage.GetLargestPossibleRegion()    
    print 'input index: {}'.format(index)
    print 'input size:  {}'.format(size)
    print 'largest possible region: '
    print largestPossibleRegion
    
    itkRegion = itk.ImageRegion[3]
    desiredRegion = itkRegion(0)
    desiredRegion.SetIndex(index)
    desiredRegion.SetSize(size)
    
    extractType = itk.ExtractImageFilter[itkImageType,itkImageType]
    extract = extractType.New()
    extract.SetExtractionRegion(desiredRegion)
    extract.SetInput(itkImage)
    extract.SetDirectionCollapseToIdentity()
    extract.Update()
    itkImageOut = extract.GetOutput()
    
    return itkImageOut

#%%    
    
def itkImage_print_direction_cosine_matrix(itkImage):
    vnl_matrix = itkImage.GetDirection().GetVnlMatrix()
    for i in range(3):
        for j in range(3):
            print "{:>8.4f}".format(vnl_matrix.get(i,j)),
        print

def itkImage_set_direction_cosine_matrix(itkImage,direction_cosine_matrix):
    """ 
    direction_cosine_matrix is assumed to be
    of the form
    
    R0 C0 S0
    R1 C1 S1
    R2 C2 S2

    so we switch R and C to accomodate ITK's column
    major order.    
    """
    for i in range(3):
        for j in range(3):
            if j == 0:
                col_idx = 1
            if j == 1:
                col_idx = 0
            if j == 2:
                col_idx = 2              
            itkImage.GetDirection().GetVnlMatrix().set(i,j,direction_cosine_matrix[i,col_idx])

def itkImage_print_image_geometry(itkImage):
    origin  = itkImage.GetOrigin()
    spacing = itkImage.GetSpacing()
    region  = itkImage.GetLargestPossibleRegion()
    size    = region.GetSize()
    print 'itkImage dims    :',  size
    print 'itkImage origin  : ', origin
    print 'itkImage spacing : ', spacing
    print 'itkImage direction cosine matrix: '
    itkImage_print_direction_cosine_matrix(itkImage)
    
def itkImage_get_image_geometry(itkImage):

    IG = ImageGeometry()
    
    # use accessors to get info from itkImage
    origin     = itkImage.GetOrigin()
    spacing    = itkImage.GetSpacing()
    region     = itkImage.GetLargestPossibleRegion()
    size       = region.GetSize()
    vnl_matrix = itkImage.GetDirection().GetVnlMatrix()         
    
    # ITK image gives origin as XYZ
    IG._origin = [float(i) for i in origin]
         
    # TODO: really get this straightened out     
    # columns of ITK vnl_matrix represent directions of R, C, then S
    IG._dir_as_row_idx_increases   = [vnl_matrix.get(i,0) for i in np.arange(3)]        
    IG._dir_as_col_idx_increases   = [vnl_matrix.get(i,1) for i in np.arange(3)]
    IG._dir_as_slice_idx_increases = [vnl_matrix.get(i,2) for i in np.arange(3)]
          
    # sampling and size given in RCS order by ITK images 
#    spacing = [float(i) for i in spacing]
#    size    = [float(i) for i in size]
    IG._samplingRCS = [spacing[0], spacing[1], spacing[2]]
    IG._shapeRCS    = [size[1]   , size[0]   , size[2]   ]

    # permute RCS --> SRC
    IG._shapeSRC    = [IG._shapeRCS[2]]    + IG._shapeRCS[0:2]
    IG._samplingSRC = [IG._samplingRCS[2]] + IG._samplingRCS[0:2]                
            
    IG._determine_orientation()
    IG._determine_span()
    IG._generate_direction_cosine_matrix()
    
    return IG
    
def itkImage_determine_orientation(itkImage):
    IG = itkImage_get_image_geometry(itkImage)
    orientation = IG._determine_orientation()
    return orientation
    

def generate_oriented_itkImage(otherImage=None,pixarray=None,ig=None):
    if isinstance(otherImage, dicom_series.DicomSeries):
        PIX = otherImage.get_pixel_array()
        imageGeometry = ImageGeometry(otherImage)

    if isinstance(otherImage, dmi.DMI):
        PIX = otherImage._PIX
        imageGeometry = otherImage._image_geometry

    if pixarray != None and ig != None: # when pixarray and ig are both specified
        PIX = pixarray
        imageGeometry = ig

    pixelComponentType = ndarray_get_itk_ctype(PIX)
    print 'itk pixel component type: ', pixelComponentType
    numberOfDimensions = 3
    
    imageType = itk.Image[ pixelComponentType, numberOfDimensions ]
    print imageType
    connector = itk.PyBuffer[imageType]
    ITKIMG = connector.GetImageFromArray(PIX)
    
    # set origin
    ITKIMG.SetOrigin(imageGeometry.origin)
    
    # set spacing
    samplingRCS = imageGeometry.samplingRCS
    ITKIMG.SetSpacing( [samplingRCS[0], samplingRCS[1], samplingRCS[2]] )
        
    # use ImageGeometry.direction_cosine_matrix to set ITK image diection cosines    
    itkImage_set_direction_cosine_matrix(ITKIMG,imageGeometry.direction_cosine_matrix)
    
#    itkImage_print_image_geometry(ITKIMG)
    
    return ITKIMG
    
    
def ds_get_itk_ctype(series): 
    #TODO a field to DicomSeries that gives rescaled pixel type
    pix = series.get_pixel_array()
    dType = ndarray_get_itk_ctype(pix)
    return dType
    

#%%

'''
NumPy I/O
'''

def ndarray_get_itk_ctype(ndarray):
    dtype = ndarray.dtype
    print 'ndarry dtype: ', dtype
    
    # floats
    if dtype == 'float32':
        return itk.F

    if dtype == 'float64':
        return itk.D

    if dtype == 'int16':
        return itk.SS

    if dtype == 'uint16':
        return itk.US

    if dtype == 'int8':
        return itk.SC

    if dtype == 'uint8':
        return itk.UC

    if dtype == 'bool':
        return itk.B

    return None

def load_itkImage_to_numpy(imgPathIn):
    imageIO = itk.ImageIOFactory.CreateImageIO(imgPathIn,itk.ImageIOFactory.ReadMode)
    imageIO.SetFileName(imgPathIn)
    imageIO.ReadImageInformation()
    
    itkImageType = itk.Image[imageIO_get_itk_ctype(imageIO), imageIO.GetNumberOfDimensions()]
    readerType = itk.ImageFileReader[ itkImageType ]
    reader = readerType.New()
    reader.SetFileName( imgPathIn )
    reader.Update()
    
    itkImage = reader.GetOutput()
    
    connector = itk.PyBuffer[ itkImageType ]
    PIX = connector.GetArrayFromImage( itkImage )
    
    return PIX

def ndarray_to_itkImage(ndarray_in, IG=None):
    
    CType = ndarray_get_itk_ctype(ndarray_in)
    Dimensionality = len(ndarray_in.shape)
    connector = itk.PyBuffer[itk.Image[CType,Dimensionality]]
    itkImage = connector.GetImageFromArray(ndarray_in)
    
    if isinstance(IG,ImageGeometry):
        itkImage.SetOrigin(IG.origin)
        itkImage.SetSpacing(IG.samplingRCS)
        itkImage_set_direction_cosine_matrix(itkImage,IG.direction_cosine_matrix)
        
    
    return itkImage
    
def itkImage_to_ndarray(itkImageIn):
    connector = itk.PyBuffer[ type(itkImageIn) ]
    PIX = connector.GetArrayFromImage( itkImageIn )
    return PIX

#%%


#class ITKImageHelper:
    
    