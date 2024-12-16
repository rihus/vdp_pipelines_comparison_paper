# -*- coding: utf-8 -*-
"""
Created on 5 Mar 2024

@author: Riaz Hussain, PhD

"""
# %% Import required libraries and declare path
import os
import re
import csv
import numpy as np
import nibabel as nib
import pandas as pd

#%% Define necessary functions
def process_nifti(in_array):
    """Load and process a nifti image"""
    nifti_read = nib.load(in_array)
    nifti_data = nifti_read.get_fdata()
    reorient_array = np.flip(np.rot90(nifti_data), 0)
    return reorient_array, nifti_read

def write_to_csv(filename, heading, data):
    """Output a csv file with heading provided and appends data each time"""
    # #Check if the file exists and is empty
    file_exists = os.path.isfile(filename)
    file_is_empty = not file_exists or os.stat(filename).st_size == 0
    # #Save file
    with open(filename, mode='a', newline='', encoding="utf-8") as file:
        writer_ = csv.writer(file)
        # Write heading only if the file is empty
        if file_is_empty:
            writer_.writerow(heading)
        writer_.writerow(data)

def load_select_defects(subj_dir, f_name, sub_dir1=None, sub_dir2=None):
    """Select defects from input array"""
    if sub_dir1 is not None and sub_dir2 is not None:
        analysis_arr = np.load(os.path.join(subj_dir, sub_dir1, sub_dir2, f_name))
        defect_mask = analysis_arr == 1
        return defect_mask
    if sub_dir1 is not None:
        analysis_arr, _ = process_nifti(os.path.join(subj_dir, sub_dir1, f_name))
        defect_mask = analysis_arr == 3
        return defect_mask
    analysis_arr, _ = process_nifti(os.path.join(subj_dir, f_name))
    return analysis_arr

def calc_diagnostic_accuracy(rdr_defects, calc_defects, corr_anlys, msk,
                              subj_dir=None):
    """Calculate Sensitivity, Specificity, Positive/Negative Predictive Values,
        and Positive/Negative Likelihood Ratios of various analysis methods"""
    rdr = np.sum((rdr_defects==1) & (msk == 1))
    tp = np.sum((rdr_defects==1) & (calc_defects==1) & (msk == 1))
    fn = np.sum((rdr_defects==1) & (calc_defects==0) & (msk == 1))
    fp = np.sum((rdr_defects==0) & (calc_defects==1) & (msk == 1))
    tn = np.sum((rdr_defects==0) & (calc_defects==0) & (msk == 1))
    sensitivity = tp / (tp + fn)
    specificity = tn / (tn + fp)
    accuracy = (tp + tn) / (tp + fn + fp + tn)
    ppv = tp / (tp + fp)
    npv = tn / (tn + fn)
    print(f"Sensitivity of {corr_anlys[0]} {corr_anlys[1]}: {sensitivity}")
    print(f"Specificity of {corr_anlys[0]} {corr_anlys[1]}: {specificity}")
    print(f"Accuracy of {corr_anlys[0]} {corr_anlys[1]}: {accuracy}")
    print(f"Positive Predictive Value of {corr_anlys[0]} {corr_anlys[1]}: {ppv}")
    print(f"Negative Predictive Value of  {corr_anlys[0]} {corr_anlys[1]}: {npv}")
    print(f"Positive Likelihood Ratio of {corr_anlys[0]} {corr_anlys[1]}:"
          f"{sensitivity / (1 - specificity)}")
    print(f"Negative Likelihood Ratio of  {corr_anlys[0]} {corr_anlys[1]}:"
          f"{(1 - sensitivity) / specificity}")
    write_to_csv(os.path.join(os.path.dirname(subj_dir),
                f"rdr_{corr_anlys[0]}_{corr_anlys[1]}_sens_spec_acc_pvs_lrs.csv"),
                ["Subject_id", "TP", "TN","FP","FN", "Reader", "Sensitivity", "Specificity",
                 "Accuracy", "PPV", "NPV", "PLR", "NLR"], [os.path.basename(subj_dir),
                tp, tn, fp, fn, rdr, sensitivity, specificity, accuracy, ppv, npv,
                sensitivity / (1 - specificity), (1 - sensitivity) / specificity])

def process_csv_files(file_list, output_folder, correction):
    """Combine all csv files into a single"""
    combined_df = None

    for file in file_list:
        # Extract analysis name
        base_name = os.path.basename(file)
        method = base_name.split('_')[2] #.split('-')[0]
        # Read the CSV file into a DataFrame
        df = pd.read_csv(file)
        # Rename columns (except for 'Subject_id') to include the analysis name
        df.rename(columns={col: f"{col}_{method}" for col in df.columns
                           if col != 'Subject_id'}, inplace=True)
        # Merge with the combined DataFrame on 'Subject_id'
        if combined_df is None:
            combined_df = df
        else:
            combined_df = pd.merge(combined_df, df, on='Subject_id', how='inner')
            # Save the combined DataFrame to CSV in the specified output folder
    output_path = os.path.join(output_folder, f"{correction}_combined_sens_spec_acc_pvs_lrs.csv")
    combined_df.to_csv(output_path, index=False)

#%% Loading defect data, calculating dice coefficient
##Path to directory containing subject folders
PARENT_DIR = (r"C:\Users\HUSDQ4\OneDrive - cchmc\cincy_work\all_projects_data_work"
              r"\vdp_analysis\analysis_comparisons\IRC740H_visit1_spiral_cf")
###Specify options
ANALYSIS_MODE = "Batch" ##Single or Batch modes
SUBJECT_ID = "IRC740H-002"
CORR = "FA" # N4, FA
SUBJ_DIR_NAMES = [r'^IRC740H-\d{3}$', r'^IRC740H-\d{3}c$', r'^ILD-HC-\d{3}$'] #
## Single subject analysis mode
if ANALYSIS_MODE == "Single":
    SUBJECT_DIR = os.path.join(PARENT_DIR, SUBJECT_ID)
    ## Load defect arrays, select defects, and load reader defects
    md_defects = load_select_defects(SUBJECT_DIR, f"{CORR}_corr_glb-median_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "glb_median_analysis")
    mn_defects = load_select_defects(SUBJECT_DIR, f"{CORR}_corr_lb-mean_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "lb_mean_analysis")
    pct_defcts = load_select_defects(SUBJECT_DIR,f"{CORR}_corr_glb-percentile_defect_array.npy",
                                    f"{CORR}_corr_vdp_analysis", "glb_percentile_analysis")
    ak_defcts = load_select_defects(SUBJECT_DIR,f"{CORR}_corr_adaptive-kmeans_defect_array.npy",
                                    f"{CORR}_corr_vdp_analysis", "kmeans_adaptive_analysis")
    hk_defcts = load_select_defects(SUBJECT_DIR,
                                    f"{CORR}_corr_hierarchical-kmeans_defect_array.npy",
                                    f"{CORR}_corr_vdp_analysis", "kmeans_hierarchical_analysis")
    th_defcts = load_select_defects(SUBJECT_DIR, f"{CORR}_corr_thresholding_defect_array.npy",
                                    f"{CORR}_corr_vdp_analysis", "thresholding_analysis")
    rdr_rh = load_select_defects(SUBJECT_DIR,"img_defect_mask_rh.nii.gz")
    ##Try to load Joey or Abood's segmented defect mask
    defects_masks = ["img_defect_mask_jp.nii.gz", "img_defect_mask_ab.nii.gz"]
    DEFECTS = None
    for mask in defects_masks:
        mask_path = os.path.join(SUBJECT_DIR, mask)
        try:
            DEFECTS, _ = process_nifti(mask_path)
            break  # Exit loop if successful
        except OSError:
            print(f"Mask {mask} not found or could not be processed.")
    if DEFECTS is None:
        print("No valid defects mask (JP or AB) was found.")
    # rdr_jp = load_select_defects(SUBJECT_DIR,"img_defect_mask_jp.nii.gz")
    msk_arr = load_select_defects(SUBJECT_DIR, "img_ventilation_mask.nii.gz")
    # img_arr = load_select_defects(SUBJECT_DIR, "img_ventilation.nii.gz")
    # img_arr /= np.max(img_arr*msk_arr)
    rdrs_arr = (rdr_rh*msk_arr) + (DEFECTS*msk_arr) == 2
    # img_n_msk = [img_arr, msk_arr]
    calc_diagnostic_accuracy(rdrs_arr, md_defects, [CORR, "Median"], msk_arr, SUBJECT_DIR)
    calc_diagnostic_accuracy(rdrs_arr, mn_defects, [CORR, "Mean"], msk_arr, SUBJECT_DIR)
    calc_diagnostic_accuracy(rdrs_arr, pct_defcts, [CORR, "Percentile"], msk_arr, SUBJECT_DIR)
    calc_diagnostic_accuracy(rdrs_arr, ak_defcts, [CORR, "Adaptive"], msk_arr, SUBJECT_DIR)
    calc_diagnostic_accuracy(rdrs_arr,hk_defcts,[CORR, "Hierarchical"],msk_arr,SUBJECT_DIR)
    calc_diagnostic_accuracy(rdrs_arr, th_defcts, [CORR, "Thresholding"], msk_arr, SUBJECT_DIR)
## Batch analysis mode
elif ANALYSIS_MODE == "Batch":
    for dirpath, dirnames, filenames in os.walk(PARENT_DIR):
        for dirname in dirnames:
            # Check if any of the patterns match the dirname
            if any(re.match(pattern, dirname) for pattern in SUBJ_DIR_NAMES):
                print(f"\nSubject folder name: {dirname}\n")
                SUBJECT_DIR = os.path.join(PARENT_DIR, dirname)
                ## Load defect arrays, select defects, and load reader defects
                md_defects = load_select_defects(SUBJECT_DIR,
                                        f"{CORR}_corr_glb-median_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "glb_median_analysis")
                mn_defects = load_select_defects(SUBJECT_DIR,
                                        f"{CORR}_corr_lb-mean_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "lb_mean_analysis")
                pct_defcts = load_select_defects(SUBJECT_DIR,
                                        f"{CORR}_corr_glb-percentile_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "glb_percentile_analysis")
                ak_defcts = load_select_defects(SUBJECT_DIR,
                                        f"{CORR}_corr_adaptive-kmeans_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "kmeans_adaptive_analysis")
                hk_defcts = load_select_defects(SUBJECT_DIR,
                                        f"{CORR}_corr_hierarchical-kmeans_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "kmeans_hierarchical_analysis")
                th_defcts = load_select_defects(SUBJECT_DIR,
                                        f"{CORR}_corr_thresholding_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "thresholding_analysis")
                rdr_rh = load_select_defects(SUBJECT_DIR,"img_defect_mask_rh.nii.gz")
                ##Try to load Joey or Abood's segmented defect mask
                defects_masks = ["img_defect_mask_jp.nii.gz", "img_defect_mask_ab.nii.gz"]
                DEFECTS = None
                for mask in defects_masks:
                    mask_path = os.path.join(SUBJECT_DIR, mask)
                    try:
                        DEFECTS, _ = process_nifti(mask_path)
                        break  # Exit loop if successful
                    except OSError:
                        print(f"Mask {mask} not found or could not be processed.")
                if DEFECTS is None:
                    print("No valid defects mask (JP or AB) was found.")
                # rdr_jp = load_select_defects(SUBJECT_DIR,"img_defect_mask_jp.nii.gz")
                msk_arr = load_select_defects(SUBJECT_DIR, "img_ventilation_mask.nii.gz")
                # img_arr = load_select_defects(SUBJECT_DIR, "img_ventilation.nii.gz")
                # img_arr /= np.max(img_arr*msk_arr)
                rdrs_arr = (rdr_rh*msk_arr) +(DEFECTS*msk_arr) == 2
                # img_n_msk = [img_arr, msk_arr]
                calc_diagnostic_accuracy(rdrs_arr, md_defects, [CORR, "Median"],
                                         msk_arr, SUBJECT_DIR)
                calc_diagnostic_accuracy(rdrs_arr, mn_defects, [CORR, "Mean"],
                                         msk_arr, SUBJECT_DIR)
                calc_diagnostic_accuracy(rdrs_arr, pct_defcts, [CORR, "Percentile"],
                                         msk_arr, SUBJECT_DIR)
                calc_diagnostic_accuracy(rdrs_arr, ak_defcts, [CORR, "Adaptive"],
                                         msk_arr, SUBJECT_DIR)
                calc_diagnostic_accuracy(rdrs_arr, hk_defcts, [CORR, "Hierarchical"],
                                         msk_arr, SUBJECT_DIR)
                calc_diagnostic_accuracy(rdrs_arr, th_defcts, [CORR, "Thresholding"],
                                         msk_arr, SUBJECT_DIR)

#%%Collect all data into FA- and N4-corrected single files
# PARENT_DIR = (r"C:\Users\HUSDQ4\OneDrive - cchmc\cincy_work\all_projects_data_work"
#               r"\vdp_analysis\analysis_comparisons\IRC740H_visit1_spiral_cf")
analysis_fldr = os.path.join(PARENT_DIR,"vdp_analysis_results_December2024")
N4_CSVs = [os.path.join(analysis_fldr, "rdr_N4_Adaptive_sens_spec_acc_pvs_lrs.csv"),
            os.path.join(analysis_fldr,"rdr_N4_Hierarchical_sens_spec_acc_pvs_lrs.csv"),
            os.path.join(analysis_fldr,"rdr_N4_Thresholding_sens_spec_acc_pvs_lrs.csv"),
            os.path.join(analysis_fldr,"rdr_N4_Median_sens_spec_acc_pvs_lrs.csv"),
            os.path.join(analysis_fldr,"rdr_N4_Percentile_sens_spec_acc_pvs_lrs.csv"),
            os.path.join(analysis_fldr,"rdr_N4_Mean_sens_spec_acc_pvs_lrs.csv")]
process_csv_files(N4_CSVs, analysis_fldr, "N4")

FA_CSVs = [os.path.join(analysis_fldr,"rdr_FA_Median_sens_spec_acc_pvs_lrs.csv"),
           os.path.join(analysis_fldr,"rdr_FA_Percentile_sens_spec_acc_pvs_lrs.csv"),
           os.path.join(analysis_fldr,"rdr_FA_Adaptive_sens_spec_acc_pvs_lrs.csv"),
           os.path.join(analysis_fldr,"rdr_FA_Hierarchical_sens_spec_acc_pvs_lrs.csv"),
           os.path.join(analysis_fldr,"rdr_FA_Thresholding_sens_spec_acc_pvs_lrs.csv"),
           os.path.join(analysis_fldr,"rdr_FA_Mean_sens_spec_acc_pvs_lrs.csv")]
process_csv_files(FA_CSVs, analysis_fldr, "FA")

#%%
