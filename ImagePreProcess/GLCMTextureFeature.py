# -*- coding: utf-8 -*-
"""
Created on 6/22/16 1:23 PM

@author: shuang Shih-ying Huang
NOTE: the GLCM feature definition is based on the following papers:
- "Ultrasound GLCM texture analysis of radiation-induced parotid-gland injury in head-and-neck cancer radiotherapy: An
in viv o study of late toxicity" Yang et al., Med Phys 2012
- "Decoding tumour phenotype by noninvasive imaging using a quantitative radiomics approach" Aerts et al., Nature communications 2014 and its supplemantary document

"""

import numpy as np
import numpy.ma as ma


def is_mask_constant(maskval,outstr):
    if maskval is ma.masked:
        print '::OH NO O_O:: {} is a masked constant!'.format(outstr)
        return True
    else:
        return False

#TODO: the p matrix is not appropriately updated!!!!!!!
class GLCMTextureFeature:
    """GLCM Texture feature"""

    def __init__(self, *args):
        self._p = None # grey-level co-occurrence matrix
        self._is_p_zeros = False
        self._ux = None
        self._uy = None
        self._sigx = None
        self._sigy = None
        self._px = None
        self._py = None
        self._k1 = None
        self._k2 = None
        self._p_xplusy = None
        self._p_xminusy = None
        self._Hx = None
        self._Hy = None
        self._Hxy = None
        self._Hxy1 = None
        self._Hxy2 = None
        self._ii = None
        self._jj = None
        self._autocor = None
        self._cluster_prominence = None
        self._cluster_tendency = None
        self._cluster_shade = None
        self._contrast = None
        self._correlation = None
        self._diff_entropy = None
        self._dissimilarity = None
        self._energy = None
        self._entropy = None
        self._homogeneity1 = None
        self._homogeneity2 = None
        self._idmn = None
        self._idn = None
        self._inverse_var = None
        self._max_prob = None
        self._sum_avg = None
        self._sum_entropy = None
        self._sum_var = None
        self._diff_var = None
        self._diff_avg = None
        self._imc1 = None
        self._imc2 = None
        self._sum_squares = None
        self._avg_intensity = None
        self._feature_list = ['contrast']  # default feature list
        self._feature_dict = None
        self._compute_feature_map = None
        self._get_feature_map = None

        if args:
            self._p = args[0]
            self.compute_basic_stats()
            if len(args) > 1:
                self._feature_list = args[1]
        else:
            print 'RAR! no gray-level co-occurence matrix is defined!'


    # ---------------------------------------------------------------------
    # Properties
    # ---------------------------------------------------------------------
    @property
    def autocorrelation(self):
        return self._autocor

    @property
    def clusterprom(self):
        # measure image asymmetry; when value is high, the image is less symmetric; when the value is low, there is a peak in the GLCM matrix around the mean values
        return self._cluster_prominence

    @property
    def clustershade(self):
        # measure skewness of a matrix and gauge the concept of uniformity; when value is high, the image is asymetric
        return self._cluster_shade

    @property
    def clustertendency(self):
        return self._cluster_tendency

    @property
    def contrast(self):
        # measure the local variation in the image; the measure of contrast favors contributions from p(i,j) away from the diagonal
        # the value is high when there is a large amount of variation in an image
        return self._contrast

    @property
    def correlation(self):
        # measure the linear dependency of gray levels on the specific neighboring voxels
        # higher values can be obtained from similar gray-level regions
        return self._correlation

    @property
    def differenceentropy(self):
        return self._diff_entropy

    @property
    def dissimilarity(self):
        return self._dissimilarity

    @property
    def energy(self):
        # measure homogeneity of an image
        # higher values == textural uniformity
        return self._energy

    @property
    def entropy(self):
        # measure the randomness of the image texture (intensity distribution)
        # value is highest when all p(i,j)'s are equal and smaller when p(i,j) are not equal
        # homogenous images will have low entropy while heterogeneous images will have higher entropy value
        return self._entropy

    @property
    def homogeneity1(self):
        return self._homogeneity1

    @property
    def homogeneity2(self):
        # measure the local homogeneity of an image
        # low value for inhomogeneous images and high value for homogeneous images
        return self._homogeneity2

    @property
    def IDNM(self):
        return self._idmn

    @property
    def IDN(self):
        return self._idn

    @property
    def inverse_var(self):
        return self._inverse_var

    @property
    def maxprob(self):
        return self._max_prob

    @property
    def sum_avg(self):
        return self._sum_avg

    @property
    def sum_entropy(self):
        return self._sum_entropy

    @property
    def sum_var(self):
        return self._sum_var

    @property
    def imc1(self):
        return self._imc1

    @property
    def imc2(self):
        return self._imc2

    @property
    def diff_var(self):
        return self._diff_var

    @property
    def diff_avg(self):
        return self._diff_avg

    @property
    def avg_intensity(self):
        return self._avg_intensity

    @property
    def sum_squares(self):
        return self._sum_squares

    @property
    def feature_dict(self):
        return self._feature_dict

    def update_p_matrix(self,theglcm):
        self._p = theglcm
        self.compute_basic_stats()

    def compute_basic_stats(self):
        if np.count_nonzero(self._p) == 0:
            print '::O_O::GLCM matrix is all ZEROS! no need for any further computation!'
            self._is_p_zeros = True
        else:
            #TODO: make all math opearation working with 3D glcm matrix instead of 2D
            #Ensure the glcm matrix is normalized to probability matrix
            self._p = ma.true_divide(self._p,np.sum(self._p))
            # print 'p: {}'.format(self._p)
            self._px = np.sum(self._p,axis=1)
            self._py = np.sum(self._p,axis=0)

            # self._ii, self._jj = np.ogrid[0:self._p.shape[0], 0:self._p.shape[1]]
            self._ii, self._jj = np.ogrid[1:(self._p.shape[0]+1), 1:(1+self._p.shape[1])]
            self._ux = np.sum(self._ii * self._p)
            self._uy = np.sum(self._jj * self._p)
            # print 'after ux, uy'

            self._sigx = np.sum(((self._ii - self._ux) ** 2) * self._p) ** 0.5
            self._sigy = np.sum(((self._jj - self._uy) ** 2) * self._p) ** 0.5

            # print 'after sigx sigy'

            # self._k1 = np.arange(0, 2 * self._p.shape[0] - 1, 1)
            self._k1 = np.arange(2, 2 * self._p.shape[0] + 1)
            self._p_xplusy = np.zeros(self._k1.shape)
            for ik in range(len(self._k1)):
                tmp = (self._ii + self._jj == self._k1[ik])
                self._p_xplusy[ik] = ma.sum(self._p[tmp])
                del tmp
            # print 'after p_xplusy'

            # self._k2 = np.arange(0, self._p.shape[0], 1)
            self._k2 = np.arange(0, self._p.shape[0])
            self._p_xminusy = np.zeros(self._k2.shape)
            for ik in range(len(self._k2)):
                tmp = (np.abs(self._ii - self._jj) == self._k2[ik])
                self._p_xminusy[ik] = np.sum(self._p[tmp])
                del tmp
            # print 'after p xminusy'

            tmp =  -np.sum(self._px * ma.log2(self._px))
            if not is_mask_constant(tmp,'Hx'):
                self._Hx = tmp
            del tmp
            # print 'after Hx'

            tmp = -np.sum(self._py * ma.log2(self._py))
            if not is_mask_constant(tmp,'Hy'):
                self._Hy = tmp
            del tmp
            # print 'after Hy'

            tmp = -np.sum(self._p * ma.log2(self._p))
            if not is_mask_constant(tmp,'H'):
                self._Hxy = tmp
            del tmp

            pxx,pyy = np.meshgrid(self._px,self._py,indexing='ij')
            tmp = -np.sum(self._p * ma.log2(pxx * pyy))
            if not is_mask_constant(tmp,'Hxy1'):
                self._Hxy1 = tmp
            del tmp
            # print 'after Hxy1: {}'.format(self._Hxy1)

            tmp = -np.sum(pxx * pyy * ma.log2(pxx * pyy))
            if not is_mask_constant(tmp,'Hxy2'):
                self._Hxy2 = tmp
            del tmp
            # print 'after Hxy2'



    def compute_features(self):

        # set up the function call dictionary
        if not self._is_p_zeros: #only compute features when the matrix is not all 0's
            self._compute_feature_map = {'autocorrelation': self._compute_autocorrelation(),
                                         'cluster_prominence': self._compute_cluster_prominence(),
                                         'cluster_shade': self._compute_cluster_shade(),
                                         'cluster_tendency': self._compute_cluster_tendency(),
                                         'contrast': self._compute_contrast(), 'correlation': self._compute_correlation(),
                                         'diff_entropy': self._compute_difference_entropy(),
                                         'dissimilarity': self._compute_dissimilarity(), 'energy': self._compute_energy(),
                                         'entropy': self._compute_entropy(), 'homogeneity1': self._compute_homogeneity1(),
                                         'homogeneity2': self._compute_homogeneity2(), 'idmn': self._compute_IDMN(),
                                         'idn': self._compute_IDN(), 'inv_var': self._compute_inverse_var(),
                                         'maxprob': self._compute_max_prob(),
                                         'sum_avg': self._compute_sum_avg(), 'sum_entropy': self._compute_sum_entropy(),
                                         'sum_var': self._compute_sum_var(), 'imc1': self._compute_IMC1(),
                                         'imc2': self._compute_IMC2(), 'diff_avg': self._compute_diff_avg(),
                                         'diff_var': self._compute_diff_var(), 'avg_intensity': self._compute_avg_intensity(),
                                         'sum_squares': self._compute_sum_squares()}

            self._get_feature_map = {'autocorrelation': [self.autocorrelation],
                                     'cluster_prominence': [self.clusterprom],
                                     'cluster_shade': [self.clustershade],
                                     'cluster_tendency': [self.clustertendency],
                                     'contrast': [self.contrast], 'correlation': [self.correlation],
                                     'diff_entropy': [self.differenceentropy],
                                     'dissimilarity': [self.dissimilarity], 'energy': [self.energy],
                                     'entropy': [self.entropy], 'homogeneity1': [self.homogeneity1],
                                     'homogeneity2': [self.homogeneity2],
                                     'idmn': [self.IDNM], 'idn': [self.IDN],
                                     'inv_var': [self.inverse_var], 'maxprob': [self.maxprob],
                                     'sum_avg': [self.sum_avg],
                                     'sum_entropy': [self.sum_entropy], 'sum_var': [self.sum_var],
                                     'imc1': [self.imc1],'imc2': [self.imc2],
                                     'diff_avg': [self.diff_avg],'diff_var': [self.diff_var],
                                     'avg_intensity': [self.avg_intensity],'sum_squares': [self.sum_squares]}

        # compute all the features for a given feature list and popuulate the feature dictionary
        self._feature_dict = {}
        for ff in self._feature_list:
            if self._is_p_zeros:
                self._feature_dict[ff] = ['nan']
                print '{}: glcm matrix is all 0\'s'.format(ff)
            else:
                self._compute_feature_map[ff]
                self._feature_dict[ff] = self._get_feature_map[ff]
                # if ff is 'correlation':
                #     print self._sigx,self._sigy
                #     print self._p
                #     print 'glcm matrix is NOT all 0\'s, compute feature map: {},{}'.format(ff,self._feature_dict[ff])

    def _compute_autocorrelation(self):
        self._autocor = np.sum(self._p * (self._ii * self._jj))

    def _compute_cluster_prominence(self):
        self._cluster_prominence = np.sum(((self._ii + self._jj - self._ux - self._uy)**4) * self._p)

    def _compute_cluster_shade(self):
        self._cluster_shade = np.sum(((self._ii + self._jj - self._ux - self._uy)**3) * self._p)

    def _compute_cluster_tendency(self):
        self._cluster_tendency = np.sum(((self._ii + self._jj - self._ux - self._uy)**2) * self._p)

    def _compute_contrast(self):
        self._contrast = np.sum(((self._ii - self._jj)**2) * self._p)

    def _compute_correlation(self):
        tmp = ma.true_divide((np.sum(self._p * (self._ii * self._jj)) - (self._ux*self._uy)),self._sigx*self._sigy)
        if not is_mask_constant(tmp,'correlation'):
            self._correlation = tmp

    def _compute_difference_entropy(self):
        # tmp = np.sum(self._p_xminusy * ma.log2(self._p_xminusy))
        tmp = (-1)*np.sum(self._p_xminusy * ma.log2(self._p_xminusy))
        if not is_mask_constant(tmp,'diff_entropy'):
            self._diff_entropy = tmp

    def _compute_dissimilarity(self):
        self._dissimilarity = np.sum(np.abs(self._ii - self._jj)*self._p)

    def _compute_energy(self):
        self._energy = np.sum(self._p**2)

    def _compute_entropy(self):
        tmp = (-1)*np.sum(self._p * ma.log2(self._p))
        if not is_mask_constant(tmp,'entropy'):
            self._entropy = tmp

    def _compute_homogeneity1(self):
        tmp = np.sum(ma.true_divide(self._p,1 + np.abs(self._ii - self._jj)))
        if not is_mask_constant(tmp,'homogeneity1'):
            self._homogeneity1 = tmp

    def _compute_homogeneity2(self):
        tmp = np.sum(ma.true_divide(self._p,(1 + (self._ii - self._jj)**2)))
        if not is_mask_constant(tmp,'homogeneity2'):
            self._homogeneity2 = tmp


    def _compute_IMC1(self):
        """TODO: determine what HXY is in the equation in Aert et al."""
        div = np.max([self._Hx, self._Hy])
        tmp = ma.true_divide((self._Hxy - self._Hxy1), div)
        if not is_mask_constant(tmp, 'IMC1'):
            self._imc1 = tmp

    def _compute_IMC2(self):
        tmp = (1 - ma.exp(-2 * (self._Hxy2 - self._Hxy))) ** 0.5
        if not is_mask_constant(tmp, 'imc2'):
            self._imc2 = tmp

    def _compute_IDMN(self):
        denom = 1 + ma.true_divide((self._ii - self._jj)**2,self._p.shape[0]**2)
        tmp = np.sum(ma.true_divide(self._p,denom))
        if not is_mask_constant(tmp,'IDMN'):
            self._idmn = tmp

    def _compute_IDN(self):
        denom = 1 + ma.true_divide(np.abs(self._ii - self._jj), self._p.shape[0])
        tmp  = np.sum(ma.true_divide(self._p, denom))
        if not is_mask_constant(tmp,'IDN'):
            self._idn = tmp

    def _compute_inverse_var(self):
        tmp = np.sum(ma.true_divide(self._p,(self._ii - self._jj)**2))
        if not is_mask_constant(tmp,'inverse_var'):
            self._inverse_var = tmp

    def _compute_max_prob(self):
        self._max_prob = np.amax(self._p)

    def _compute_sum_avg(self):
        # print 'compute sum avg: {}, {}'.format(self._k1, self._p_xplusy)
        self._sum_avg = np.sum(self._k1 * self._p_xplusy)

    def _compute_sum_entropy(self):
        self._sum_entropy = -np.sum(self._p_xplusy * ma.log2(self._p_xplusy))

    def _compute_sum_var(self): # same as cluster tendency (definition has some variation as shown below)
        # self._sum_var = np.sum((self._k1 - self._sum_avg) ** 2 * self._p_xplusy)  # according to image biomarker standard
        self._sum_var = np.sum((self._k1 - self._sum_entropy)**2 * self._p_xplusy)  # according to http://earlglynn.github.io/RNotes/package/EBImage/Haralick-Textural-Features.html

    def _compute_diff_avg(self):
        self._diff_avg = np.sum(self._k2 * self._p_xminusy)

    def _compute_diff_var(self):
        diffavg = np.sum(self._k2 * self._p_xminusy)
        self._diff_var = np.sum(self._p_xminusy * (self._k2 - diffavg) ** 2)

    def _compute_avg_intensity(self):
        self._avg_intensity = self._ux

    def _compute_sum_squares(self):
        self._sum_squares = np.sum(self._p * ((self._ii - self._ux) ** 2))







