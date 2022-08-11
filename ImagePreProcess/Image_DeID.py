# -*- coding: utf-8 -*-
"""
Created on 3/6/17 4:37 PM

@author: shuang Shih-ying Huang
@goal: de-identify images, e.g. MRN, patient name, and accession #

"""

import dicom
import glob
import os
import sys
import re
from dicom.errors import InvalidDicomError

# github source: https://github.com/darcymason/pydicom/blob/master/pydicom/examples/anonymize.py
def anonymize(filename, output_filename, new_person_name='anonymous',
              new_patient_id='id', new_acc_num = 'acc_num', remove_curves=True, remove_private_tags=True):
    """Replace data element values to partly anonymize a DICOM file.
    Note: completely anonymizing a DICOM file is very complicated; there
    are many things this example code does not address. USE AT YOUR OWN RISK.
    """

    # Define call-back functions for the dataset.walk() function
    def PN_callback(ds, data_element):
        """Called from the dataset "walk" recursive function for all data elements."""
        if data_element.VR == 'PN':
            # print 'found patient name!'
            data_element.value = new_person_name
        # else:
        #     print data_element.VR

    def curves_callback(ds, data_element):
        """Called from the dataset "walk" recursive function for all data elements."""
        if data_element.tag.group & 0xFF00 == 0x5000:
            del ds[data_element.tag]

    # Load the current dicom file to 'anonymize'
    dataset = dicom.read_file(filename)

    # Remove patient name and any other person names
    dataset.walk(PN_callback)

    # Change ID
    dataset.PatientID = new_patient_id

    # Change Accession Number
    dataset.AccessionNumber = new_acc_num

    # Remove data elements (should only do so if DICOM type 3 optional)
    # Use general loop so easy to add more later
    # Could also have done: del ds.OtherPatientIDs, etc.
    for name in ['OtherPatientIDs', 'OtherPatientIDsSequence']:
        if name in dataset:
            delattr(dataset, name)

    # Same as above but for blanking data elements that are type 2.
    for name in ['PatientBirthDate']:
        if name in dataset:
            dataset.data_element(name).value = ''

    # Remove private tags if function argument says to do so. Same for curves
    if remove_private_tags:
        dataset.remove_private_tags()
    if remove_curves:
        dataset.walk(curves_callback)

    # write the 'anonymized' DICOM out under the new filename
    dataset.save_as(output_filename)


if __name__ == "__main__":

    # Lassen
    rootdir = '/data/francgrp1'
    imdir = '{}/breast_radiomics/her2/MRI'.format(rootdir)
    all_dir = [d for d in glob.glob('{}/*'.format(imdir)) if os.path.isdir(d)]
    all_pt_id = [int(re.search(r'{}/(\d+)'.format(imdir), ss).group(1)) for ss in all_dir
                 if re.search(r'{}/(\d+)'.format(imdir), ss)]

    ids = [ii for ii in all_pt_id if ii != 88]  # don't have pt 88 json file for now...
    # ids ids = [1]
    # go through all .dmi in xxx/xxx/dce for each MRI_pt_id
    for pt_id in ids:
        series_dir = glob.glob('{}/{:0>3d}/*/*'.format(imdir, pt_id))[0]
        out_dir = '{}/dce_anonymized'.format(series_dir)
        if os.path.exists(out_dir):
            if not os.path.isdir(out_dir):
                raise IOError('O_O: Input is directory! output name exists but is not a directory')
        else:
            os.makedirs(out_dir)

        dce_dir = os.path.join(series_dir, 'dce')
        files = os.listdir(dce_dir)

        for ff in files:
            if ff.endswith('.dmi'):
                input_file = os.path.join(dce_dir, ff)
                output_file = os.path.join(out_dir, ff)
                print input_file
                try:
                    anonymize(input_file, output_file, remove_curves=False, remove_private_tags=False)
                except InvalidDicomError:
                    print("{} is NOT a valid dicom file, may need force=True on read_file\r".format(input_file))
                else:
                    print("done\r")

    # go through all mask dmi and anonymize
    mask_dir = '{}/MORE_INFO'.format(imdir)
    mask_out_dir = '{}/MASK_Anonymized'.format(imdir)
    if os.path.exists(mask_out_dir):
        if not os.path.isdir(mask_out_dir):
            raise IOError('O_O: Input is directory! output name exists but is not a directory')
    else:
        os.makedirs(mask_out_dir)

    files = os.listdir(mask_dir)
    for ff in files:
        if ff.endswith('.dmi'):
            input_file = os.path.join(mask_dir, ff)
            output_file = os.path.join(mask_out_dir, ff)
            print input_file
            try:
                anonymize(input_file, output_file, remove_curves=False, remove_private_tags=False)
            except InvalidDicomError:
                print("{} is NOT a valid dicom file, may need force=True on read_file\r".format(input_file))
            else:
                print("done\r")