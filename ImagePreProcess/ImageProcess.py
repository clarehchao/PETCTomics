# -*- coding: utf-8 -*-
"""
Created on 6/22/16 9:50 AM

@author: shuang Shih-ying Huang
@Goal: a collection of functions for image segmentation, image display, and VOI extraction (built for breast MRI)

"""

import itk
import ITKImageHelper
import matplotlib.pyplot as plt
import numpy as np
import json
import skimage.morphology as skm
from skimage.feature import peak_local_max
from scipy import ndimage as ndi
from skimage import measure
import operator
import AxesSequence as axs
import matplotlib.patches as patches
from matplotlib.colors import colorConverter
import matplotlib as mpl
import re
import nrrd
import image_geometry


def display_volume(theVOI,fig_title,aspect=None,colormap=None,climarray=None):
    """ scroll through the slices in a VOI by pressing left or right arrow"""

    # plt.switch_backend('Qt4Agg')
    axes = axs.AxesSequence()
    for i,ax in zip(range(theVOI.shape[0]),axes):
        imgplot = ax.imshow(theVOI[i])
        if climarray:
            imgplot.set_clim(climarray)
        if colormap:
            imgplot.set_cmap(colormap)
        if aspect:
            ax.set_aspect(aspect)
        ax.set_title('{} slice # {}'.format(fig_title,i))
    axes.show()
    del axes

def display_overlay_volume(theVOI,theMask,fig_title,aspect=1.0):
    """ scroll through the slices in a VOI by pressing left or right arrow"""

    plt.switch_backend('Qt4Agg')

    # generate the colors for your colormap
    color1 = colorConverter.to_rgba('White')
    color2 = colorConverter.to_rgba('Green')
    cmap2 = mpl.colors.LinearSegmentedColormap.from_list('my_cmap2', [color1, color2], 256)
    cmap2._init()  # create the _lut array, with rgba values
    alphas = np.linspace(0, 0.5, cmap2.N + 3)
    cmap2._lut[:, -1] = alphas

    axes = axs.AxesSequence()
    for i,ax in zip(range(theVOI.shape[0]),axes):
        im1 = ax.imshow(np.squeeze(theVOI[i]),cmap=plt.cm.gray,aspect=aspect)
        ax.hold(True)
        im2 = ax.imshow(np.squeeze(theMask[i]).astype('uint8'), cmap=cmap2, aspect=aspect)
        # plt.savefig('{}/MPLCenterslice_overlay.pdf'.format(outputrootdir))
        ax.set_title('{} slice # {}'.format(fig_title,i))
    axes.show()
    del axes


def ComputeSER(data,timestamp_idx):
    """ Compute Signal Enhancement Ratio (SER) of the dynamic breast MRI images
        SER = (s1 - s0)/(s2-s0)

        Input:
        - data: a 4-D numpy array with shape of (# of time stamps, VOI.shape)
        - timestamp_idx: the time indx wrt to the dimension of data.shape[0] in a chronological order
    """
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        tmp = np.true_divide(data[timestamp_idx[1]] - data[timestamp_idx[0]],
                             data[timestamp_idx[2]] - data[timestamp_idx[0]])
        # #  nan_to_num replaces info with 1.0e+308 and NAN with 0.0
        # SER = np.nan_to_num(SER)

    tmp = np.squeeze(tmp)
    SER = np.zeros(tmp.shape)
    finite_mask = np.isfinite(tmp)
    SER[finite_mask] = tmp[finite_mask]

    return SER

def ComputePEMIP(data,timestamp_idx):
    """ Compute Percent Enhancement (PE) maximum intensity projection wrt time
        PE = 100*(s(t) - s0)/s0

        Input:
        - data: a 4-D numpy array with shape of (# of time stamps, VOI.shape)
        - timestamp_idx: the time indx wrt to the dimension of data.shape[0] in a chronological order

    """
    PE_all = np.zeros((data.shape[0] - 1, data.shape[1], data.shape[2], data.shape[3]))
    idx_select = timestamp_idx[1:]
    for ii in range(len(idx_select)):
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            tmp = 100.*np.true_divide(data[timestamp_idx[idx_select[ii]]] - data[timestamp_idx[0]],
                                        data[timestamp_idx[0]])
            tmp = np.squeeze(tmp)
            finite_mask = np.isfinite(tmp)
            # print 'ii = {}, min = {}, max = {}'.format(ii,np.amin(tmp[finite_mask]),np.amax(tmp[finite_mask]))
            PE_all[ii,finite_mask] = tmp[finite_mask]
            # PE_all[ii] = np.nan_to_num(np.squeeze(tmp))

    # perform maximum intensity projection of PE(t) over the 2 time points
    PE_MIP = np.amax(PE_all, axis=0)

    return PE_MIP

def TumorSegmentation(VOI_SER,VOI_PEMIP,thresh_ser,thresh_pe,aspect=1.0):
    """ Segment tumors based on SER and PE_MIP of a given dynamic breast MRI images """

    # Compute the central slice for testing
    center_slice = int(0.5 * VOI_SER.shape[0])

    # Create an initial mask based on a given SER and PEMIP threshold via global thresholding
    mask1 = np.logical_and(VOI_SER > thresh_ser[0], VOI_SER <= thresh_ser[1])
    mask2 = np.logical_and(VOI_PEMIP > thresh_pe[0], VOI_PEMIP <= thresh_pe[1])
    # mask3 = np.logical_and(mask1, mask2)
    # try to consider both mask result
    mask3 = np.logical_or(mask1, mask2)
    im = np.squeeze(mask3[center_slice])
    # plt.imshow(im, cmap='gray')
    # plt.show()

    # label the mask via connected component
    labels_all, nLabel = measure.label(mask3, connectivity=3, return_num=True)
    # print np.unique(labels_all)
    # display_volume(labels_all,'labels_all')

    # remove small objects, this may not be necessary since the histogram will be used to filter it out later...
    labels_filter = skm.remove_small_objects(labels_all, connectivity=3)
    # print np.unique(labels_filter)
    # display_volume(labels_filter, 'remove small objects')


    # compute a histogram of the labels to figure out which label has the largest # of voxels
    hist_label, bin_edge = np.histogram(labels_filter.ravel(), bins=np.arange(1, np.amax(labels_filter), 1), range=(1, np.amax(labels_filter)))
    # find the label with the largest # of voxels to pick out potential tumor
    theLabel_indx = np.argmax(hist_label)
    theLabel = bin_edge[theLabel_indx]
    # plt.bar(bin_edge[1:], hist_label)
    # plt.show()

    # Create a mask of the label with the largest # of voxels (hopefully the tumor!)
    mask4 = (labels_all == theLabel)
    # display_volume(mask4,'mask4')

    # Perform region growing via skimage watershed algorithm
    distance = ndi.distance_transform_edt(mask4)
    # plt.imshow(distance[center_slice], cmap='gray',aspect=aspect)
    # plt.show()

    local_maxi = peak_local_max(distance, indices=False)
    markers = measure.label(local_maxi)
    # plt.imshow(local_maxi[center_slice],cmap='gray')
    # plt.show()

    # A label of the tumor portion with enhanced signal (likely the vessel enhancement)
    labels_ws = skm.watershed(-distance, markers, mask=mask4)
    # display_volume(labels_ws, 'watershed seg')
    # plt.imshow(labels_ws[center_slice], cmap='gray')
    # plt.show()

    # Create a label with any potential holes filled in the mask to estimate tumor volume
    labels_rm = skm.binary_closing(labels_ws, np.ones((10, 10, 10)))
    # display_volume(labels_rm, 'binary closing')

    # TODO: random walker algorithm doesn't seem to work well or took a long computation time and didn't give any result...
    # labels_rw = segmentation.random_walker(mask4, markers)
    # plt.imshow(labels_rw[center_slice], cmap='gray')
    # plt.show()

    return labels_ws, labels_rm

def TumorSegmentation(theVol,theMask,theImGeo,theVOICenter,theVOIHalflen,aspect=1.0,is_plot=False):
    """ Segment tumors based on a given mask """

    # get all VOIs for faster computation
    theDatalist = [theVol,theMask]
    theVOIlist,VOI_indx,win_ixyz,win_size = GetVOI(theDatalist, theImGeo, theVOICenter, theVOIHalflen)
    theVol_VOI, theMask_VOI = theVOIlist

    center_slice = int(0.5 * theVol_VOI.shape[0])
    if is_plot:
        display_volume(theMask_VOI, 'the input mask VOI')
        # plt.imshow((theMask == 0)[center_slice])
        # plt.colorbar()
        # plt.show()
        # plt.imshow((theMask == 1)[center_slice])
        # plt.colorbar()
        # plt.show()

    # label the mask via connected component
    labels_all, nLabel = measure.label(theMask_VOI, connectivity=3, return_num=True)
    if is_plot:
        display_volume(labels_all,'the initial connected component label: labels_all')

    # remove small objects, this may not be necessary since the histogram will be used to filter it out later...
    print('before remove_small_object, # of labels: {}'.format(len(np.unique(labels_all))))
    labels_filter = skm.remove_small_objects(labels_all, connectivity=3,min_size=32)
    print('after remove_small_object, # of labels: {}'.format(len(np.unique(labels_filter))))

    if is_plot:
        # for ll in np.unique(labels_filter):
        #     plt.imshow((labels_filter==ll)[center_slice], cmap='gray', aspect=aspect)
        #     plt.title('label = {}'.format(ll))
        #     plt.show()
        display_volume(labels_filter, 'after remove small objects')

    if len(np.unique(labels_filter)) > 2:
        # compute a histogram of the labels to figure out which label has the largest # of voxels
        bin_max = np.amax(labels_filter)
        hist_label, bin_edge = np.histogram(labels_filter.ravel(), bins=np.linspace(1.,bin_max+1,bin_max+1), range=(1,bin_max+1))
        print hist_label, bin_edge
        # plt.bar(bin_edge[1:], hist_label)
        # plt.show()

        # # Create a mask of the label with the largest # of voxels (hopefully the tumor!)
        # theLabel_indx = np.argmax(hist_label)
        # theLabel = bin_edge[theLabel_indx]
        # print 'theLabel = {}'.format(theLabel)
        # initial_tumor_mask = (labels_filter == theLabel)

        # find the label with the largest # of voxels to pick out potential tumor regions
        the_max_count_label_indx = np.argmax(hist_label)

        # find the labels that are close to the maximum count
        the_max_count = hist_label[the_max_count_label_indx]
        count_thresh_ratio = 0.4
        the_largest_count_label_indx = np.where(hist_label >= count_thresh_ratio*the_max_count)
        the_largest_count_label = bin_edge[the_largest_count_label_indx]
        print('theLabel = {}'.format(the_largest_count_label))

        # Create a mask of the label with the largest # of voxels (hopefully the tumor!)
        initial_tumor_mask = np.in1d(labels_filter.ravel(),the_largest_count_label).reshape(labels_filter.shape)

    else:
        initial_tumor_mask = labels_filter

    if is_plot:
        display_volume(initial_tumor_mask, 'initial_tumor_mask')

    # Perform region growing via skimage watershed algorithm
    distance = ndi.distance_transform_edt(initial_tumor_mask)
    # bin_max = np.amax(distance)
    # hist, bin_edge = np.histogram(distance.ravel(), bins=np.linspace(0., bin_max, bin_max*3),range=(0, bin_max))
    # maxdist = bin_edge[np.argmax(hist)]
    # print np.diff(bin_edge)
    # print maxdist
    # plt.bar(bin_edge[1:], hist,width=np.diff(bin_edge)[0])
    # plt.show()

    if is_plot:
        plt.imshow(distance[center_slice], cmap='gray',aspect=aspect)
        plt.show()

    # set min_distance = 2 instead of 1
    local_maxi = peak_local_max(distance, min_distance=1.3,indices=False)
    initial_reg_grow_seed = measure.label(local_maxi)
    if is_plot:
        plt.imshow(initial_reg_grow_seed[center_slice],cmap='gray')
        plt.show()

    # A label of the tumor portion with enhanced signal (likely the vessel enhancement)
    tumor_labels_watershed = skm.watershed(-distance, initial_reg_grow_seed, mask=initial_tumor_mask)
    if is_plot:
        # display_volume(labels_watershed, 'watershed seg')
        plt.imshow(tumor_labels_watershed[center_slice], cmap='gray')
        plt.show()

    # Create a label with any potential holes filled in the mask to estimate tumor volume
    tumor_labels_watershed_binclose = skm.binary_closing(tumor_labels_watershed, np.ones((3, 3, 3)))
    if is_plot:
        display_volume(tumor_labels_watershed_binclose, 'binary closing')
        # # Visual checking of the segmented tumor via image overlay
        # display_overlay_volume(theVol_VOI,label_tumorall,'tumor label overlay',aspect=ar)


    # make an image mask with original image size
    theTumorMask = np.zeros(theVol.shape,dtype='uint8')
    theTumorMask[VOI_indx[0,0]:VOI_indx[0,1], VOI_indx[1,0]:VOI_indx[1,1],VOI_indx[2,0]:VOI_indx[2,1]] = tumor_labels_watershed_binclose.astype('uint8')
    return tumor_labels_watershed, tumor_labels_watershed_binclose,theTumorMask

def DisplayLabelOverlayImage(theVOI,theMask,voi_thresh,theImgGeo,outputrootdir,pt_id,is_plot=False):
    """ Overlay a label image with the original image and display for visual check"""

    if is_plot:
        center_slice = int(0.5 * theVOI.shape[0])
        ar = theImgGeo.samplingSRC[2]/theImgGeo.samplingSRC[1]

        # more details from http://stackoverflow.com/questions/10127284/overlay-imshow-plots-in-matplotlib
        # or http://matplotlib.org/examples/pylab_examples/layer_images.html

        # generate the colors for your colormap
        color1 = colorConverter.to_rgba('White')
        color2 = colorConverter.to_rgba('Green')

        # make the colormaps
        cmap2 = mpl.colors.LinearSegmentedColormap.from_list('my_cmap2', [color1, color2], 256)
        cmap2._init()  # create the _lut array, with rgba values
        alphas = np.linspace(0, 0.7, cmap2.N+3)
        cmap2._lut[:, -1] = alphas
        fig = plt.figure(frameon=False)
        im1 = plt.imshow(np.squeeze(theVOI[center_slice]),cmap=plt.cm.gray,clim=voi_thresh,aspect=ar)
        plt.hold(True)
        im2 = plt.imshow(np.squeeze(theMask[center_slice]).astype('uint8'),cmap=cmap2,aspect=ar)
        plt.savefig('{}/MPLCenterslice_overlay.pdf'.format(outputrootdir))
        plt.show()


    # Save the RGB overlay volume into a 24-bit unsigned char file (RGB 8-bit data) via ITK
    # cast the label image as unsigned char instead of bool for ITK I/O
    theMask_UC = theMask.astype('uint8')

    # cast the original VOI to float or unsigned char
    theVOI_UC = theVOI.astype('uint8')

    # overlay label with original images
    # NOTE: this only works with the ITKImageHelper on TumorSeg/shuang branch
    itkImage_voi, imageType_voi = ITKImageHelper.generate_oriented_itkImage(pixarray=theVOI_UC, ig=theImgGeo)
    ITKImageHelper.itkImage_print_image_geometry(itkImage_voi)
    itkImage_label, imageType_label = ITKImageHelper.generate_oriented_itkImage(pixarray=theMask_UC, ig=theImgGeo)
    ITKImageHelper.itkImage_print_image_geometry(itkImage_label)

    # adapted from c++ example code https://github.com/InsightSoftwareConsortium/ITKWikiExamples/blob/master/ImageProcessing/LabelOverlayImageFilter.cxx
    itkRGBtype = itk.Image.RGBUC3
    labeloverlayFilter = itk.LabelOverlayImageFilter[imageType_voi, imageType_label, itkRGBtype].New()
    labeloverlayFilter.SetInput(itkImage_voi)
    labeloverlayFilter.SetLabelImage(itkImage_label)
    labeloverlayFilter.SetOpacity(0.3)
    labeloverlayFilter.Update()
    labeloutput = labeloverlayFilter.GetOutput()

    # TODO: PyBuffer doesn't work itk.RGBUC3 type so have to go with another approach to display the overlay images (see below)
    # connector = itk.PyBuffer[itkRGBtype]
    # overlaypixarray = connector.GetArrayFromImage(labeloutput)
    # plt.imshow(overlaypixarray[center_slice])
    # plt.show()

    writerType = itk.ImageFileWriter[itkRGBtype]
    writer = writerType.New()
    writer.SetInput(labeloutput)
    itkoutftag = 'ITKTumorSegOverlay'
    writer.SetFileName('{}/{}.hdr'.format(outputrootdir, itkoutftag))
    writer.Update()

    if is_plot:
        # read the saved RGB files for display with matplotlib
        mpfile = '{}/{}.img'.format(outputrootdir, itkoutftag)
        img = np.fromfile(mpfile, dtype=np.uint8)

        # reshape the ITK saved file to display RGB properly (24-bit RGB, or each color is 8-bit)
        newshape = theVOI.shape + (3,)
        img = img.reshape(newshape)
        display_volume(img, 'pt_id: {}, tumor seg overlay'.format(pt_id))

def GetImInfo(jsonfile):
    """ Read the json file from brtool """
    with open(jsonfile) as f:
        data = json.load(f)
        VOI_LPS = data['VOI_LPS']
        jCENTER = VOI_LPS['CENTER']
        jHALFLENGTHS = VOI_LPS['HALFLENGTHS']
        jIMNAMES = data['IMAGE_NAMES']

    HL0 = [jHALFLENGTHS[0][0], jHALFLENGTHS[1][0], jHALFLENGTHS[2][0]]
    HL1 = [jHALFLENGTHS[0][1], jHALFLENGTHS[1][1], jHALFLENGTHS[2][1]]
    HL2 = [jHALFLENGTHS[0][2], jHALFLENGTHS[1][2], jHALFLENGTHS[2][2]]
    HALFLENGTHS = [HL0, HL1, HL2]
    CENTER = [jCENTER[0][0], jCENTER[1][0], jCENTER[2][0]]
    IMNAMES = str(jIMNAMES[0]).split('\\')

    # order time point
    time_list = {}
    for ii in range(len(IMNAMES)):
        check = re.search(r'.+_tp(\d).dmi',IMNAMES[ii])
        if check:
            tp = check.group(1)
            time_list[tp] = IMNAMES[ii]
        else:
            print 'cannot determine time point!'

    sorted_time_list = sorted(time_list.items(), key=operator.itemgetter(0))
    time_sorted_IMNAMES = np.array([sorted_time_list[ii][1] for ii in range(len(sorted_time_list))])

    return CENTER,HALFLENGTHS,time_sorted_IMNAMES

def GetVOIinfo(theCenter,theHalfLen,theImgGeo):
    for idx in range(3):
        side1 = np.add(theCenter, theHalfLen[idx])
        side2 = np.subtract(theCenter, theHalfLen[idx])
        tmp = np.vstack((side1, side2))
        if idx == 0:
            SIDE_POINTS = tmp
        else:
            SIDE_POINTS = np.vstack((SIDE_POINTS, tmp))

    SIDE_MIN = np.amin(SIDE_POINTS, axis=0)
    SIDE_MAX = np.amax(SIDE_POINTS, axis=0)

    minr, minc, mins = theImgGeo.coords_to_idx(SIDE_MIN[0], SIDE_MIN[1], SIDE_MIN[2])
    maxr, maxc, maxs = theImgGeo.coords_to_idx(SIDE_MAX[0], SIDE_MAX[1], SIDE_MAX[2])
    mr = min(minr, maxr)
    mc = min(minc, maxc)
    ms = min(mins, maxs)
    rlen = abs(minr - maxr)
    clen = abs(minc - maxc)
    slen = abs(mins - maxs)
    win_ixyz = [ms, mc, mr]
    win_size = [slen, clen, rlen]

    return win_ixyz,win_size


def GetVOI(theDatalist,theImgGeo,theCenter,VOI_ixyz,VOI_size,theThresh=None,is_display=False):
    """ Generate a VOI for a given patient dataset"""

    VOI_indx = np.array([[VOI_ixyz[0], VOI_ixyz[0] + VOI_size[0]], [VOI_ixyz[1],VOI_ixyz[1]+VOI_size[1]], [VOI_ixyz[2],VOI_ixyz[2]+VOI_size[2]]])

    # get the VOI
    theVOIlist = []
    for ii in range(len(theDatalist)):
        thePIX = theDatalist[ii]
        VOI = thePIX[VOI_indx[0,0]:VOI_indx[0,1], VOI_indx[1,0]:VOI_indx[1,1],VOI_indx[2,0]:VOI_indx[2,1]]

        # # determine initial threshold of VOI via skimage filters => this helps for TumorSegmenation(SER,PE_MIP...etc)
        # tmp = VOI.ravel()
        # data_thresh = tmp[np.logical_and(tmp < np.nanmax(VOI), tmp > 0.01)]
        #
        # # can try ostu, li, and yen method as well...
        # lo_thresh = skf.threshold_isodata(data_thresh)
        # hi_thresh = 0.2*np.nanmax(VOI)
        # theThresh = [lo_thresh, hi_thresh]

        # display for checking
        if is_display:
            cr, cc, cs = theImgGeo.coords_to_idx(theCenter[0], theCenter[1], theCenter[2])
            SLICE = thePIX[cs, :, :]
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
            color_map = plt.cm.gray
            imgplt = ax.imshow(SLICE, aspect=(theImgGeo.samplingSRC[2] / theImgGeo.samplingSRC[1]), cmap=color_map)
            imgplt.set_clim(theThresh)
            ax.add_patch(patches.Rectangle((VOI_ixyz[2], VOI_ixyz[1]), VOI_size[2], VOI_size[1], fill=False, color='red'))
            plt.show()
        theVOIlist.append(VOI)
        del thePIX
        del VOI

    return theVOIlist

def determine_dce_seriesid(pt_summary_file):

    dce_sid = None
    with open(pt_summary_file,'r') as f:
        lines = f.readlines()
        for ll in lines:
            ll_nospace = ''.join(ll.split()).lower()
            if ll_nospace.find('dyn') != -1:  # this should cover the case for 'dynamic' or 'dyn' in the text
                check = re.search(r'Series(\d+).+', ll_nospace, re.IGNORECASE)
                if check:
                    dce_sid = int(check.group(1))
                    break
    return dce_sid

def StandardizeImageIntensity(input_array,bin_info,method='fixed bin size'):
    """
    :param input_array: a N-D numpy array
    :param bin_info: bin_width if method = 'fixed bin size', number of bins if method = 'fixed bin number'
    :param method: the intensity standardization method
    :return: the standarized output numpy array
    """
    output_array = np.zeros(input_array.shape)
    # default standardization method: fixed bin size
    if method is 'fixed bin size':
        bin_width = bin_info
        output_array = np.round(np.true_divide(input_array,bin_width))
        output_array = output_array.astype('int16')
    elif method is 'fixed bin number':
        bin_number, max_percentile  = bin_info
        intensity_max = np.percentile(input_array,max_percentile)
        intensity_min = np.amin(input_array)
        intensity_res = (intensity_max - intensity_min)/bin_number
        output_array = np.round(np.true_divide(input_array - intensity_min,intensity_res))
        output_array[output_array > bin_number] = bin_number - 1
        output_array = output_array.astype('int16')

    return output_array



def GetITKVOI(input_itkimg,voi_size,voi_ixyz):
    """

    :param input_itkimg: input itk image
    :param voi_size: index defined in numpy fashion
    :param voi_ixyz: index order defined in numpy fashion
    :return: output roi itk image
    """
    img_dim = len(voi_size)
    img_type = type(input_itkimg)
    roiType = itk.RegionOfInterestImageFilter[img_type, img_type]
    roi = roiType.New()
    roi.SetInput(input_itkimg)

    window = itk.ImageRegion[img_dim]()
    size = itk.Size[img_dim]()
    size[0] = voi_size[2]
    size[1] = voi_size[1]
    size[2] = voi_size[0]
    window.SetSize(size)
    window.SetIndex(0, voi_ixyz[2])
    window.SetIndex(1, voi_ixyz[1])
    window.SetIndex(2, voi_ixyz[0])
    roi.SetRegionOfInterest(window)
    roi.Update()
    roi_itkimg = roi.GetOutput()

    return roi_itkimg

def NRRDImageCast(nrrd_fname,cast_dtype):
    """

    :param nrrd_fname: the file name of the .nrrd file
    :param cast_dtype: the datatype to cast the nrrd data to
    :return: a nrrd image file with the cast data type
    """
    # numpy to ITK type dictionary
    ctype2nptype_dict = {}
    ctype2nptype_dict['unsigned char'] = 'uint8' # 0 - 255
    ctype2nptype_dict['float'] = 'float32' # 3.4E +- 38
    ctype2nptype_dict['double'] = 'float64' #1.7E +- 308
    ctype2nptype_dict['short'] = 'int16' # -32768 - 32767
    ctype2nptype_dict['unsigned short'] = 'uint16' # 0 - 65535

    # read the nrrd file via python nrrd package
    data,options = nrrd.read(nrrd_fname)

    cast_dtype_nptype = ctype2nptype_dict[cast_dtype]
    if data.dtype != cast_dtype_nptype:
        # dtype_info = np.iinfo(cast_dtype_nptype)
        # norm_data = (float(dtype_info.max) * (data - np.amin(data)) / (np.amax(data) - np.amin(data))).astype(cast_dtype_nptype)
        norm_data = ((data - np.amin(data)) / (np.amax(data) - np.amin(data))).astype(cast_dtype_nptype)
        norm_options = options
        norm_options['type'] = cast_dtype
        fname_tag_check = re.search('(.+).nrrd',nrrd_fname)
        if fname_tag_check:
            fname_tag = fname_tag_check.group(1)
        else:
            print 'NRRDImageCast: nrrd_fname has not file extension .nrrd'

        cast_nrrd_fname = '{}_{}.nrrd'.format(fname_tag,cast_dtype.replace(' ','_'))
        nrrd.write(cast_nrrd_fname,norm_data,norm_options)
        print 'NRRDImageCast: write the cast file {}'.format(cast_nrrd_fname)
    else:
        print 'NRRDImageCast: no need to cast the image, same variable type!'


def ITK_Image_OverlayPlot(itk_img,itk_mask,fig_title=''):
    """

    :param itk_img: input itk image
    :param itk_mask: input itk mask
    :param fig_title: figure title
    :param aspect: figure aspect ratio
    NOTE: plot the overlay of itk image and the mask, plot only the images slices that the mask != 0
    """
    ig_img = image_geometry.ImageGeometry(itk_img)
    ar = ig_img.samplingRCS[0]/ig_img.samplingRCS[1]
    img_array = ITKImageHelper.itkImage_to_ndarray(itk_img)
    mask_array = ITKImageHelper.itkImage_to_ndarray(itk_mask)
    check = np.nonzero(mask_array)

    s1_idx = np.min(check[0])
    s2_idx = np.max(check[0])
    print(s1_idx, s2_idx)

    # # this works with the cellsite mouse data read in from dicom_series
    # s1_idx = np.min(check[2])
    # s2_idx = np.max(check[2])

    display_overlay_volume(img_array[s1_idx:s2_idx+1],mask_array[s1_idx:s2_idx+1],fig_title=fig_title,aspect=ar)