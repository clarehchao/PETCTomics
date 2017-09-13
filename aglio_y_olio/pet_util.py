# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 18:51:09 2016

@author: rharnish, shuang
"""

#%%

#import dicom
from dicom.tag import Tag
import datetime
import math
import re

#%%


def get_SUV_multiplier(pet_series):
    '''
    pet_series should be DicomSeries object    
    
    SUV = C(T)/[decay-corrected injection dose (MBq)/patient's weight (kg)]
    
    TODO: check out
    http://www.turkupetcentre.net/petanalysis/model_suv.html
    
    '''
    info = pet_series.info
    
    # units
    tag_Units   = Tag(0x0054,0x1001)
    units_field = info[tag_Units]
    units       = units_field.value
    print 'image units:         {}'.format(units)
    
    # weight
    tag_Patient_Weight = Tag(0x0010,0x1030)
    pat_weight_field   = info[tag_Patient_Weight]
    pat_weight         = pat_weight_field.value
    print 'patient weight (kg): {}'.format(pat_weight)
    weight_units_factor = 1000.
    
    # injected dose
    tag_Radiopharmaceutical_Information_Sequence = Tag(0x0054,0x0016)
    tag_Radionuclide_Total_Dose                  = Tag(0x0018,0x1074)
    tag_Radiopharmaceutical_start_time           = Tag(0x0018,0x1072)
    tag_SeriesTime                               = Tag(0x0008,0x0031)
    tag_Radionuclide_Half_life                   = Tag(0x0018,0x1075)
    rad_pharm_seq        = info[tag_Radiopharmaceutical_Information_Sequence]

    try:
        rad_total_dose_field = rad_pharm_seq[0][tag_Radionuclide_Total_Dose]
    except:
        print '::O_O:: cannot get the total dose in the image DICOM header'
        return 0
    rad_total_dose       = rad_total_dose_field.value
    rad_half_life = rad_pharm_seq[0][tag_Radionuclide_Half_life].value  # unit: seconds
    print 'radionuclide total dose: {}'.format(rad_total_dose)

    # decay corrected the injected dose (or total dose)
    rad_pharm_start_time = rad_pharm_seq[0][tag_Radiopharmaceutical_start_time].value
    if re.search(r'\.',rad_pharm_start_time):
        dt_rad_pharm_start_time = datetime.datetime.strptime(rad_pharm_start_time,'%H%M%S.%f')
    else:
        dt_rad_pharm_start_time = datetime.datetime.strptime(rad_pharm_start_time, '%H%M%S')

    series_time = info[tag_SeriesTime].value
    if re.search(r'\.',series_time):
        dt_series_time = datetime.datetime.strptime(series_time,'%H%M%S.%f')
    else:
        dt_series_time = datetime.datetime.strptime(series_time,'%H%M%S')
    dt_time_elapse = dt_series_time - dt_rad_pharm_start_time
    time_elapse_sec = dt_time_elapse.seconds + dt_time_elapse.microseconds*1e-6  #unit: seconds
    rad_lambda = math.log(2)/rad_half_life  #unit: 1/sec
    decaycorr_rad_total_dose = rad_total_dose*math.exp(-time_elapse_sec*rad_lambda)

    print 'the decay-corrected total dose (MBq):    {}'.format(decaycorr_rad_total_dose)
    
    suv_multiplier = float(weight_units_factor*pat_weight) / float(decaycorr_rad_total_dose)
    
    return suv_multiplier 