# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 15:14:12 2015

@author: rharnish
"""

from image_geometry import ImageGeometry


def get_ct_and_pet_series_from_series_list(seriesList):
    petSeries = []
    ctSeries  = []
    
    for series in seriesList:
        
        serNum = series.info.SeriesNumber
        
        # make sure series is 3D
        if len(series.shape) < 3:
            print 'not at least 3D: {}'.format(serNum)
            continue        
            
        # throw out series that are likely reports or scouts due to having 
        # number of minimal slices
        if series.shape[0] < 10:
            print 'not enough slices in series {} to use for analysis'.format(serNum)
            continue
        
        if series.uniform_orientation == False:
            print 'not all instances in series {} have same orientation, so not using'.format(serNum)
            continue
    
        # try to choose best PET and CT series to use     
        if series.info.Modality == 'PT':
            print 'adding series {} to PET list'.format(serNum)
            petSeries.append(series)
        if series.info.Modality == 'CT':
            print 'adding series {} to CT list'.format(serNum)
            ctSeries.append(series)
            
    return ctSeries, petSeries
    
def test_z_extent(IG,min_z_extent=700):
    # consider series to be candidate for use
    # as primary CT series if has "z" extent 
    # greater than some amount
    if not IG.span[2] > min_z_extent:
        print 'series does not have large enough z extent'
        return False
    else:
        return True
        
def test_axial(IG):
    # remove non-axial series from list
    if not IG.orientation == 'axial':
        print 'series is not axial'
        return False
    else:
        return True

def determine_best_ct_series(ctSeries):
    
    bestCT = None    
    highestRes = [999999,999999]
    
    for series in ctSeries:
        
        IG = ImageGeometry(series) 
        
        if not test_z_extent(IG):
            continue
        
        if not test_axial(IG):
            continue
        
        if series.sampling[1] * series.sampling[2] < highestRes[0] * highestRes[1]:
            highestRes = series.sampling[1:]
            bestCT = series
            
    return bestCT

def determine_best_pet_series(petSeries):
    
    bestPET = None
    highestRes = [999999,999999]
    
    for series in petSeries:
        
        IG = ImageGeometry(series)
        
        if not test_z_extent(IG):
            continue
        
        if not test_axial(IG):
            continue
    
        # TODO: this should be a function that more thoroughly determines
        # whether and how the series was ACed
        print
        try:
            acMethodField = series.info.AttenuationCorrectionMethod
        except:
            acMethodField = 'NONE'
        print 'AC: ', acMethodField
        if 'NONE' in acMethodField:
            print 'AC method description contains "NONE," removing...'
            continue
        if 'measured' in acMethodField:
            print '"measured" in AC method field description. Not getting rid of series...'
            
        if series.sampling[1] * series.sampling[2] < highestRes[0] * highestRes[1]:
            highestRes = series.sampling[1:]
            bestPET = series
            
    return bestPET

# put series info in db  
def getSeriesInfoForDB(series):
    seriesInfo = {}
    seriesInfo['SeriesNumber']   = series.info.SeriesNumber
    seriesInfo['Sampling']       = series.sampling
    seriesInfo['Shape']          = series.shape
    seriesInfo['SliceDirection'] = series.slice_direction
    return seriesInfo