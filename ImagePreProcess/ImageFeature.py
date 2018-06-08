# -*- coding: utf-8 -*-
"""
Created on 6/14/16 12:53 PM

@author: shuang Shih-ying Huang
@goal: extract image features including shape, volume, and textural features via itk

"""

import numpy as np
import itk
import GLCMTextureFeature as gtf
import image_geometry
import scipy.stats as ss
from skimage import measure as skm
import numpy.ma as ma
import ITKImageFilters as itkif
import ITKImageHelper
import matplotlib.pylab as plt
import pandas as pd


class ImageFeature:
    """ImageFeature"""

    def __init__(self, *args):

        self._inputImage = None
        self._maskImage = None
        self._cast_maskImage = None
        self._GLCM = None
        self._inputImageType = None
        self._inputImageDimension = None
        self._maskImageType = None
        self._GLCM_feature_list = None
        self._foutput_format = None
        self._dict_feature_output = None
        self._df_feature_output = None
        self._inputpix_min = None
        self._inputpix_max = None
        self._maskpix_min = None
        self._maskpix_max = None
        self._bin_width = 0.1 #default
        self._Nbin = 128 #Nbin
        # self._glcm_GSnorm_Nbin = 256 #default number of gray level to standardize the input image intensity
        self._IG = None
        self._maskImage_ndarray = None
        self._inputImage_ndarray = None

        if args:
            self._inputImage = args[0]

            #TODO: check if inputImage and maskImage is itkImage type
            if len(args) > 1:
                self._GLCM_feature_list = args[1]
                print 'feature list is defined!'
            if len(args) > 2:
                #TODO: need to consider multiple VOIs or mask (or multi-value mask)
                self._maskImage = args[2]
                print 'mask image is defined!'
                self._foutput_format = args[3]
                print 'feature output format: {}'.format(self._foutput_format)

            self.get_image_info()
        else:
            print '::Oh NO O_O:: No itk image is defined and default texture feature list will be used!'


    @property
    def GLCM(self):
        """ """
        return self._GLCM

    @property
    def inputImageType(self):
        return self._inputImageType

    @property
    def inputImageDimension(self):
        return self._inputImageDimension

    @property
    def maskImageType(self):
        return self._maskImageType

    @property
    def feature_output(self):
        if self._foutput_format == 'dict':
            return self._dict_feature_output
        elif self._foutput_format == 'df':
            return self._df_feature_output
        else:
            print 'invalid feature output format'

    @property
    def glcm_Nbin(self):
        return self._glcm_Nbin

    @property
    def input_pix_min(self):
        return self._inputpix_min

    @property
    def input_pix_max(self):
        return self._inputpix_max

    @property
    def IG(self):
        return self._IG

    @property
    def inputImage_ndarray(self):
        return self._inputImage_ndarray

    @property
    def maskImage_ndarray(self):
        return self._maskImage_ndarray

    def get_image_info(self):
        self._inputImageType = type(self._inputImage)
        self._inputImageDimension = self._inputImage.GetImageDimension()

        # get the numpy array from the image and make sure the images are in the right orientation
        connector = itk.PyBuffer[self._inputImageType]
        self._IG = image_geometry.ImageGeometry(self._inputImage)
        print 'ImageFeature class: image voxel size {}'.format(self._IG.samplingRCS)
        tmp = connector.GetArrayFromImage(self._inputImage)
        self._inputImage_ndarray = np.zeros(tmp.shape,dtype=tmp.dtype)
        self._inputImage_ndarray[tmp != 0] = tmp[tmp !=0]
        self._inputpix_min = np.amin(self._inputImage_ndarray)
        self._inputpix_max = np.amax(self._inputImage_ndarray)
        print '::ImageFeature:: input image min and max: {}, {}'.format(self._inputpix_min, self._inputpix_max)
        # del connector
        # del tmp

        if self._maskImage:
            self._maskImageType = type(self._maskImage)

            # cast the mask image to the same image type as the input image (image type needs to be the same when passing mask image to an itk glcm object)
            castImageFilter = itk.CastImageFilter[self._maskImageType,self._inputImageType].New()
            castImageFilter.SetInput(self._maskImage)
            castImageFilter.Update()
            self._cast_maskImage = castImageFilter.GetOutput()

            connector = itk.PyBuffer[self._inputImageType]
            tmp = connector.GetArrayFromImage(self._cast_maskImage)
            self._maskImage_ndarray = np.zeros(tmp.shape)
            self._maskImage_ndarray[tmp != 0] = tmp[tmp != 0]
            self._inside_pixel_val = int(np.unique(self._maskImage_ndarray[np.nonzero(self._maskImage_ndarray)])[0])
            self._maskImage_ndarray = self._maskImage_ndarray.astype('bool')
            print(self._maskImage_ndarray.shape)
            print(self._inputImage_ndarray.shape)

            # compute the image min and max within the mask
            self._maskpix_min = np.amin(self._inputImage_ndarray[self._maskImage_ndarray])
            self._maskpix_max = np.amax(self._inputImage_ndarray[self._maskImage_ndarray])
            print '::ImageFeature:: mask min and max: {}, {}'.format(self._maskpix_min, self._maskpix_max)
            # del connector
            # del tmp

        print '::ImageFeature:: complete get_image_info!'


    def _compute_texture_features(self, Rneighbor, binWidth=None):
        """Compute the grey-level co-occurrence matrix in 3D fashion
           Assume inputImage and inputMask are in ITK format
        """
        # # Rescale the ITK volume into the range of 0 - Nbin
        # if GSnorm_Nbin:
        #     self._glcm_GSnorm_Nbin = GSnorm_Nbin
        # inputImage_GSnorm = itkif.ITKStandardizeImageIntensity(self._inputImage,self._glcm_GSnorm_Nbin)

        # set up itk glcm filter
        # ImageToCoOccuranceType = itk.ScalarImageToCooccurrenceMatrixFilter[type(inputImage_GSnorm)]

        ImageToCoOccuranceType = itk.ScalarImageToCooccurrenceMatrixFilter[type(self._inputImage)]
        glcmGenerator = ImageToCoOccuranceType.New()
        # print glcmGenerator.GetMin(), glcmGenerator.GetMax(),glcmGenerator.GetNumberOfBinsPerAxis()
        # self._Nbin = glcmGenerator.GetNumberOfBinsPerAxis()

        #binWidth is more robust than setting the Nbin as sometimes certain Nbin doesn't produce enough quantification for a glcm
        if binWidth:
            self._bin_width = binWidth
        self._Nbin = int(np.round(np.ceil((self._inputpix_max - self._inputpix_min)/self._bin_width)))
        # self._Nbin = int(np.round((self._maskpix_max - self._maskpix_min) / self._bin_width))
        print '::IMageFeature:: set the GLCM Nbin to {}'.format(self._Nbin)

        glcmGenerator.SetNumberOfBinsPerAxis(self._Nbin)
        # glcmGenerator.SetPixelValueMinMax(int(self._maskpix_min), int(self._maskpix_max))
        glcmGenerator.SetPixelValueMinMax(int(self._inputpix_min),int(self._inputpix_max))
        # glcmGenerator.SetPixelValueMinMax(0,self._glcm_GSnorm_Nbin)
        # glcmGenerator.SetInput(inputImage_GSnorm)
        glcmGenerator.SetInput(self._inputImage)
        if self._cast_maskImage:
            glcmGenerator.SetMaskImage(self._cast_maskImage)
            if self._inside_pixel_val != 1:
                print '::ImageFeature:: set inside pixel value to {}'.format(self._inside_pixel_val)
                glcmGenerator.SetInsidePixelValue(self._inside_pixel_val)

        NeighborhoodType = itk.Neighborhood.F3
        neighborhood = NeighborhoodType()
        neighborhood.SetRadius(Rneighbor)
        centerIndx = neighborhood.GetCenterNeighborhoodIndex()

        self._GLCM = np.zeros((centerIndx,self._Nbin,self._Nbin))
        offset_list = []
        dict_foutput_tmp = []
        for d in range(centerIndx):
            offset = neighborhood.GetOffset(d)
            offset_tuple = tuple([int(s) for s in tuple(offset)])
            offset_list.append(offset_tuple)
            glcmGenerator.SetOffset(offset)
            glcmGenerator.Update()
            itkhistogram = glcmGenerator.GetOutput()

            for i in range(0, itkhistogram.GetSize()[0]):
                for j in range(0, itkhistogram.GetSize()[1]):
                    self._GLCM[d,i,j] = itkhistogram.GetFrequency((i, j))

            # print 'd = {}, glcm: {}'.format(d, self._GLCM[d] / np.sum(self._GLCM[d]))
            if d == 0:
                thegtf = gtf.GLCMTextureFeature(self._GLCM[d], self._GLCM_feature_list)
            else:
                thegtf.update_p_matrix(self._GLCM[d])

            thegtf.compute_features()
            thefeature = thegtf.feature_dict


            if d == 0:
                if self._foutput_format == 'df':
                    dict_foutput_tmp = thefeature
                elif self._foutput_format == 'dict':
                    self._dict_feature_output = thefeature
            else:
                # append features for each direction as a list
                if self._foutput_format == 'df':
                    for k,val in thefeature.items():
                        dict_foutput_tmp[k].append(val[0])
                elif self._foutput_format == 'dict':
                    for k, val in thefeature.items():
                        # print k, val
                        self._dict_feature_output[k].append(val[0])
            # del thegtf
            # del thefeature
        if self._foutput_format == 'df':
            # add the offset list to the dictionary
            dict_foutput_tmp['glcm_offset'] = offset_list

            # make sure the key string is prefixed with 'texture_' to distinguish from other features
            for k,val in dict_foutput_tmp.items():
                new_k = 'texture_' + k
                dict_foutput_tmp[new_k] = dict_foutput_tmp.pop(k)

            # convert to pandas dataframe
            self._df_feature_output = pd.DataFrame(dict_foutput_tmp)
        elif self._foutput_format == 'dict':
            self._dict_feature_output['glcm_offset'] = offset_list
            # make sure the key string is prefixed with 'texture_' to distinguish from other features
            for k, val in self._dict_feature_output.items():
                new_k = 'texture_' + k
                self._dict_feature_output[new_k] = self._dict_feature_output.pop(k)

        print '::ImageFeature:: complete compute_texture_features!'


    def _compute_shape_size_features(self):
        """ compute volume under a given mask"""

        # mask volume
        vox_size_SRC = 0.1*np.array(self._IG.samplingSRC)  # unit: cm
        vox_vol = np.prod(vox_size_SRC) # unit: cm^
        vol_cm3 = np.sum(self._maskImage_ndarray)*vox_vol# unit: ml or cm^3
        if self._foutput_format == 'df':
            self._df_feature_output['ShapeSize_vol_cm3'] = vol_cm3
        elif self._foutput_format == 'dict':
            self._dict_feature_output['ShapeSize_vol_cm3'] = vol_cm3


        # mask surface area
        the_tumor_vol = np.zeros(self._inputImage_ndarray.shape)
        the_tumor_vol[self._maskImage_ndarray] = self._inputImage_ndarray[self._maskImage_ndarray]
        verts, faces, _, _ = skm.marching_cubes_lewiner(the_tumor_vol, 0.0, tuple(vox_size_SRC))
        # verts, faces = skm.marching_cubes(the_tumor_vol, 0.0,tuple(vox_size_SRC))
        surf_area_cm2 = skm.mesh_surface_area(verts, faces) # unit: cm^2

        if self._foutput_format == 'df':
            self._df_feature_output['ShapeSize_surf_area_cm2'] = surf_area_cm2
            self._df_feature_output['ShapeSize_compactness1'] = vol_cm3/(np.sqrt(np.pi)*surf_area_cm2**(2./3.))
            self._df_feature_output['ShapeSize_compactness2'] = (36.*np.pi*vol_cm3**2)/(surf_area_cm2**3)
            self._df_feature_output['ShapeSize_sphericity'] = (np.pi**(1./3.))*(6*vol_cm3)**(2./3.)/surf_area_cm2
            self._df_feature_output['ShapeSize_surface2volratio'] = surf_area_cm2/vol_cm3
        elif self._foutput_format == 'dict':
            self._dict_feature_output['ShapeSize_surf_area_cm2'] = surf_area_cm2
            self._dict_feature_output['ShapeSize_compactness1'] = vol_cm3 / (np.sqrt(np.pi) * surf_area_cm2 ** (2. / 3.))
            self._dict_feature_output['ShapeSize_compactness2'] = (36. * np.pi * vol_cm3 ** 2) / (surf_area_cm2 ** 3)
            self._dict_feature_output['ShapeSize_sphericity'] = (np.pi ** (1. / 3.)) * (6 * vol_cm3) ** (2./3.)/surf_area_cm2
            self._dict_feature_output['ShapeSize_surface2volratio'] = surf_area_cm2 / vol_cm3

        # maximum 3D euclidean distance (or diameter?)
        print '::ImageFeature:: # of verts is {}'.format(len(verts))
        eucdis_tmp = []
        for n in range(len(verts) - 1):
            px1 = verts[n + 1:]
            px2 = np.tile(verts[n], (px1.shape[0], 1))
            eucdis = np.sqrt(np.sum((px1 - px2) ** 2, axis=1))
            eucdis_tmp.append(np.amax(eucdis))
        max_euc_dis = max(eucdis_tmp)
        # print '::ImageFeature:: max_euc_dis based on mesh verts: {}'.format(max_euc_dis)
        if self._foutput_format == 'df':
            self._df_feature_output['ShapeSize_max_euc_dis'] = max_euc_dis
        elif self._foutput_format == 'dict':
            self._dict_feature_output['ShapeSize_max_euc_dis'] = max_euc_dis


        # kk,ii,jj = np.where(self._maskImage_ndarray == True)
        # print '::ImageFeature:: len(kk) = {}'.format(len(kk))

        # # skip computing max euclidean distance if the mask is too big such as including skin and etc.
        # if len(kk) > 300000 and vol_cm3 > 70.:
        #     print '::ImageFeature:: compute_shape_size_feature, tumor mask is TOO big (vol: {} ml)! will not compute euc max distance :/'.format(vol_cm3)
        # else:
        #     print '::ImageFeature:: max euc distance # of voxels to go through: {}, vol = {} cm3'.format(len(kk),vol_cm3)
        #     eucdis_tmp = np.zeros((len(kk),1))
        #     for n in range(len(kk)-1):
        #         px1 = np.column_stack((kk[n+1:],ii[n+1:],jj[n+1:]))
        #         px2 = np.tile(np.array([kk[n],ii[n],jj[n]]),(px1.shape[0],1))
        #         vox_size_tile = np.tile(vox_size_SRC,(px1.shape[0],1))
        #         eucdis = np.sqrt(np.sum(((px1 - px2)*vox_size_tile)**2,axis=1))
        #         eucdis_tmp[n] = np.amax(eucdis)
        #     max_euc_dis = np.amax(eucdis_tmp)
        #     print '::ImageFeature:: mac euc dis (original voxels) = {}'.format(max_euc_dis)
        #     if self._foutput_format == 'df':
        #         self._df_feature_output['ShapeSize_max_euc_dis'] = max_euc_dis
        #     elif self._foutput_format == 'dict':
        #         self._dict_feature_output['ShapeSize_max_euc_dis'] = max_euc_dis

        # R is the radius of the sphere with the same volume as the tumor
        tumor_sphere_R = (3*vol_cm3/(4*np.pi))**(1./3)
        if self._foutput_format == 'df':
            self._df_feature_output['ShapeSize_spherical_disproportion'] = surf_area_cm2/(4*np.pi*tumor_sphere_R**2)
        elif self._foutput_format == 'dict':
            self._dict_feature_output['ShapeSize_spherical_disproportion'] = surf_area_cm2/(4*np.pi*tumor_sphere_R**2)

        print '::ImageFeature:: complete compute_shape_size_features!'

    def _compute_first_order_stats(self,binWidth=None):
        """
            compute first-order statistics of the data
            Nbin: the number of discrete intensity levels to compute the data histogrm, default Nbin = 128
        """
        if self._maskImage_ndarray is None:  # no mask is defined, so use all the entire image data
            data = self._inputImage_ndarray.flatten()
        else:
            data = self._inputImage_ndarray[self._maskImage_ndarray]

        #TODO: here is the bottleneck when data is all zero's
        # data_stats = dict(ss.describe(data)._asdict())
        # data_stats.pop('nobs')  # delete the dict entry of 'nobs'
        data_stats = {}
        data_stats['minmax'] = (np.min(data), np.max(data))
        data_stats['mean'] = np.mean(data)
        data_stats['variance'] = np.var(data)
        data_stats['kurtosis'] = ss.kurtosis(data, fisher=False)
        data_stats['skewness'] = ss.skew(data)

        # make sure each key is prefixed with 'FOstats_'
        for k, val in data_stats.items():
            new_k = 'FOstats_' + k
            data_stats[new_k] = data_stats.pop(k)

        if self._foutput_format == 'df':
            self._df_feature_output.update(data_stats)
        elif self._foutput_format == 'dict':
            self._dict_feature_output.update(data_stats)

        # compute histogram-related stats
        # density = True ==> the integral of p_data = 1.0, i.e. np.sum(p_data*np.diff(p_bin)) = 1.0
        # p_data,p_bin = np.histogram(data,bins=Nbin,density=True)

        if binWidth:
            self._bin_width = binWidth
        lowBound = np.min(data) - (np.min(data) % self._bin_width)
        # Add + binwidth to ensure the maximum value is included in the range generated by numpu.arange
        highBound = np.max(data) + self._bin_width

        bin_edge = np.arange(lowBound,highBound,self._bin_width)
        p_data = np.histogram(data, bin_edge)[0]
        eps = np.spacing(1)
        p_data = p_data / (np.sum(p_data) + eps)

        tmp = (-1)*np.sum(p_data*ma.log2(p_data))
        if tmp is ma.masked:
            print '::Oh NO O_O:: FOstats_entropy is a masked constant!!'
        else:
            if self._foutput_format == 'df':
                self._df_feature_output['FOstats_entropy'] = tmp
            elif self._foutput_format == 'dict':
                self._dict_feature_output['FOstats_entropy'] = tmp

        if self._foutput_format == 'df':
            self._df_feature_output['FOstats_energy'] = np.sum(data**2)
            self._df_feature_output['FOstats_uniformity'] = np.sum(p_data**2)
        elif self._foutput_format == 'dict':
            self._dict_feature_output['FOstats_energy'] = np.sum(data ** 2)
            self._dict_feature_output['FOstats_uniformity'] = np.sum(p_data**2)

        print '::ImageFeature:: complete compute_first_order_stats!'


















