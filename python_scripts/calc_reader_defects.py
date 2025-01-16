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
                r"\analysis_comparisons\IRC740H_visit1_spiral_cf")

# SUBJECTS = ["ILD-HC-001", "ILD-HC-002","ILD-HC-006", "ILD-HC-009", "ILD-HC-012",
# "ILD-HC-027", "ILD-HC-032", "ILD-HC-039", "ILD-HC-044", "ILD-HC-067",
# "IRC740H-032c", "IRC740H-035c", "IRC740H-037c", "IRC740H-038c", "IRC740H-039c",
# "IRC740H-040c", "IRC740H-041c", "IRC740H-042c", "IRC740H-043c", "IRC740H-044c",
# "IRC740H-046c", "IRC740H-047c", "IRC740H-048c", "IRC740H-050c", "IRC740H-057c"]
SUBJECTS= ["IRC740H-001", "IRC740H-002", "IRC740H-003", "IRC740H-004", "IRC740H-005",
           "IRC740H-006", "IRC740H-007", "IRC740H-008", "IRC740H-009", "IRC740H-010",
           "IRC740H-012", "IRC740H-013", "IRC740H-014", "IRC740H-015", "IRC740H-016",
           "IRC740H-017", "IRC740H-018", "IRC740H-019", "IRC740H-020", "IRC740H-021",
           "IRC740H-022", "IRC740H-023", "IRC740H-024", "IRC740H-025", "IRC740H-027",
           "IRC740H-028", "IRC740H-030", "IRC740H-031", "IRC740H-033", "IRC740H-034",
           "IRC740H-045", "IRC740H-049", "IRC740H-052", "IRC740H-053", "IRC740H-054",
           "IRC740H-055", "IRC740H-056", "IRC740H-064"]
# SUBJECTS= ["IRC740H-028"]
# SUBJECTS = ["HYP-001-001", "HYP-001-002", "HYP-001-003", "HYP-001-004", "HYP-001-005",
#             "HYP-001-006", "HYP-001-007", "HYP-001-008", "HYP-001-009", "HYP-001-010",
#             "HYP-001-011", "HYP-002-001", "HYP-002-002", "HYP-002-003", "HYP-002-004",
#             "HYP-002-005", "HYP-002-006", "HYP-002-007", "HYP-002-008", "HYP-002-009",
#             "HYP-002-010", "HYP-003-002", "HYP-003-003", "HYP-003-004", "HYP-003-005",
#             "HYP-003-006", "HYP-003-007", "HYP-003-008", "HYP-003-009", "HYP-004-001",
#             "HYP-004-002", "HYP-004-003", "HYP-004-004", "HYP-004-005", "HYP-004-006",
#             "HYP-004-008", "HYP-004-009", "HYP-004-011"]

# SUBJECTS = ["ILD-HC-001", "ILD-HC-002", "ILD-HC-009", "ILD-HC-012", "ILD-HC-027",
#             "ILD-HC-032", "ILD-HC-039", "ILD-HC-044", "IRC740H-043c", "IRC740H-047c"]

# SUBJECTS = ["HYP-002-001", "HYP-002-002", "HYP-002-003", "HYP-002-004", "HYP-002-005",
#             "HYP-002-006", "HYP-002-007", "HYP-002-008", "HYP-002-009", "HYP-002-010"]
NUM=0
for SUBJECT_ID in SUBJECTS:
# SUBJECT_ID = "IRC740H-002"
    print(f"\nProcessing Subject: {SUBJECT_ID}\n")
    SUBJ_DIR = os.path.join(PARENT_DIR, SUBJECT_ID)
    ##Try to load Joey or Abood's segmented defect mask
    defects_masks = ["img_defect_mask_jp.nii.gz", "img_defect_mask_ab.nii.gz"]
    DEFECTS = None
    for mask in defects_masks:
        mask_path = os.path.join(SUBJ_DIR, mask)
        try:
            DEFECTS, _ = process_nifti(mask_path)
            break  # Exit loop if successful
        except OSError:
            print(f"Mask {mask} not found or could not be processed.")
    if DEFECTS is None:
        print("No valid defects mask (JP or AB) was found.")
    ###Load Riaz's segmented defect mask
    rh_defects, _ = process_nifti(os.path.join(SUBJ_DIR, "img_defect_mask_rh.nii.gz"))
    if DEFECTS is None:
        common_defects = rh_defects
    else:
        common_defects = rh_defects + DEFECTS == 2
    ##Load mask
    try:
        # Try to load the mask file file
        MSK_data, msk = process_nifti(glob.glob(os.path.join(SUBJ_DIR, "*4defects.nii.gz"))[0])
    except IndexError:
        MSK_data, msk = process_nifti(os.path.join(SUBJ_DIR, "img_ventilation_mask.nii.gz"))
    NUM = NUM + 1
    rdr_defect_masked = common_defects * MSK_data
    ##Calculate reader scoring vdp
    Reader_vdp = (np.count_nonzero(rdr_defect_masked)/np.count_nonzero(MSK_data))*100
    ##Save all vdp results
    vdp_csv_path = os.path.join(PARENT_DIR, "Reader_segmented_vdp.csv")
    vdp_csv_heading = ["Subject_id", "Reader"]
    vdp_results = [SUBJECT_ID, Reader_vdp]
    # print(vdp_results)
    write_to_csv(vdp_csv_path, vdp_csv_heading, vdp_results)
print(f"\nTotal processed {NUM}")

#%%
