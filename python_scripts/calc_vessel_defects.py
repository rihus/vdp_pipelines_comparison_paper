# -*- coding: utf-8 -*-
"""
Created on 18 Dec 2023

@author: Riaz Hussain, PhD

Code to generate coefficient of variation maps

"""
#%%Import libraries and define functions
import os
import csv
import glob
import numpy as np
import nibabel as nib

#%%Required functions
def process_nifti(in_array):
    """Load and process a nifti image"""
    nifti_read = nib.load(in_array)
    nifti_data = nifti_read.get_fdata()
    reorient_array = np.flip(np.rot90(nifti_data), 0)
    return reorient_array, nifti_read

def write_to_csv(filename, heading, data):
    """Output a csv file with heading provided and appends data each time"""
    # Check if the file exists and is empty
    file_exists = os.path.isfile(filename)
    file_is_empty = not file_exists or os.stat(filename).st_size == 0

    with open(filename, mode='a', newline='', encoding="utf-8") as file:
        writer_ = csv.writer(file)
        # Write heading only if the file is empty
        if file_is_empty:
            writer_.writerow(heading)
        writer_.writerow(data)

#%% Importing files, loading data, abd making a montage
### Path to directory containing subject folder
#SUBJ_DIR = r"C:\Users\HUSDQ4\Desktop\test\IRC740H-001" ##age_matched_spiral_healthy
PARENT_DIR = (r"C:\Users\HUSDQ4\OneDrive - cchmc\cincy_work\all_projects_data_work\vdp_analysis"
                r"\analysis_comparisons\age_matched_spiral_healthy")

SUBJECTS = ["ILD-HC-001", "ILD-HC-002","ILD-HC-006", "ILD-HC-009", "ILD-HC-012",
"ILD-HC-027", "ILD-HC-032", "ILD-HC-039", "ILD-HC-044", "ILD-HC-067",
"IRC740H-032c", "IRC740H-035c", "IRC740H-037c", "IRC740H-038c", "IRC740H-039c",
"IRC740H-040c", "IRC740H-041c", "IRC740H-042c", "IRC740H-043c", "IRC740H-044c",
"IRC740H-046c", "IRC740H-047c", "IRC740H-048c", "IRC740H-050c", "IRC740H-057c"]

# SUBJECTS = ["ILD-HC-001", "ILD-HC-002", "ILD-HC-009", "ILD-HC-012", "ILD-HC-027",
#             "ILD-HC-032", "ILD-HC-039", "ILD-HC-044", "IRC740H-043c", "IRC740H-047c"]

NUM=0
for SUBJECT_ID in SUBJECTS:
# SUBJECT_ID = "IRC740H-002"
    print(f"\nProcessing Subject: {SUBJECT_ID}\n")
    SUBJ_DIR = os.path.join(PARENT_DIR, SUBJECT_ID)
    ###Load Riaz's segmented defect mask
    vessel_msk, _ = process_nifti(os.path.join(SUBJ_DIR, "vessel_mask.nii.gz"))
    ##Load mask
    try:
        # Try to load the mask file file
        MSK_data, msk = process_nifti(glob.glob(os.path.join(SUBJ_DIR, "*4defects.nii.gz"))[0])
    except IndexError:
        MSK_data, msk = process_nifti(os.path.join(SUBJ_DIR, "img_ventilation_mask.nii.gz"))
    NUM = NUM + 1
    vessels_masked = vessel_msk * MSK_data
    ##Calculate reader scoring vdp
    vessel_vdp = (np.count_nonzero(vessels_masked)/np.count_nonzero(MSK_data))*100
    ##Save all vdp results
    vdp_csv_path = os.path.join(PARENT_DIR, "Vessels_vdp.csv")
    vdp_csv_heading = ["Subject_id", "vesselsVDP"]
    vdp_results = [SUBJECT_ID, vessel_vdp]
    # print(vdp_results)
    write_to_csv(vdp_csv_path, vdp_csv_heading, vdp_results)
print(f"\nTotal processed {NUM}")

#%%
