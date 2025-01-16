#!/usr/bin/env python3
"""Document created by Riaz Hussain, PhD"""
#%% Libraries import
import os
import re
import numpy as np
import nibabel as nib
from skimage import io, img_as_ubyte
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

#%%Define Functions
def process_nifti(in_array):
    """Load and process a nifti image"""
    nifti_read = nib.load(in_array)
    nifti_data = nifti_read.get_fdata()
    reorient_array = np.flip(np.rot90(nifti_data), 0)
    return reorient_array, nifti_read

def create_montage_array(montage_in, slices=None, msk_file=None):
    """
    Creates a 2D N_pixel x (N_slice * N_pixel (x N_Channels)) array from a
    3D (4D) N_pixel x N_pixel x N_slice (x N_channel) array.
    Args:
        montage_in (ndarray): 3D/4D array, in form N_pixel x N_pixel x N_slice (x N_channel).
        msk_file (ndarray): 3D array with the same (3) dimensions as montage_in (optional).
        slices (list of ints): Slices to make into montage. Defaults to middle 7 slices.
        Other options: 'all' plots all the slices if msk_file=None. 
                        'all' plots non-zero(=masked) slices if mask provided.
    Returns:
        montage_output (ndarray): Numpy array for making montage.
    """
    if slices is None and msk_file is None:
        mid_im = montage_in.shape[2] // 2
        slices = [mid_im - 3, mid_im - 2, mid_im - 1, mid_im, mid_im + 1, mid_im + 2, mid_im + 3]
    elif slices == "all" and msk_file is None:
        slices = list(range(montage_in.shape[2]))
    else:
        slices = np.flatnonzero(np.sum(msk_file, axis=(0, 1)))
    if np.ndim(montage_in) == 3:
        montage_output = montage_in[:, :, slices[0]]
        for i in slices[1:]:
            montage_output = np.hstack((montage_output, montage_in[:, :, i]))
    elif np.ndim(montage_in) == 4:
        montage_output = montage_in[:, :, slices[0], :]
        for i in slices[1:]:
            montage_output = np.hstack((montage_output, montage_in[:, :, i, :]))
    return montage_output

def montage_plot_4d(montage_in):
    """To plot the 4D (3D image + color overlaid) montage image
        montage_in (ndarray): 4D array, in form N_pixel x N_pixel x N_slice (x N_channel),
        should be processed through create_montage_array() function already
    """
    fig, ax = plt.subplots(figsize=(16, 3), dpi=300)
    ax.imshow(montage_in, alpha=1)
    plt.axis('off')
    return fig

def save_image(fig, save_path, data_type=None, im_type=None):
    """To save image in the provided path
        save_path (str): path to directory for saving the image.
        data_type (str): type of data (None (default), Non_corr, N4_corr, FA_corr)
        im_type (str): optional image data type (e.g. raw, mask, defect, montage)
    """
    if data_type is not None:
        filename = os.path.join(save_path, f"{data_type}_img_{im_type}.png")
    else:
        ## benaam from Urdu means without/unknown name
        filename = os.path.join(save_path, 'benaam_img.png')
    i = 0
    while os.path.exists(filename):
        i += 1
        filename = f"{os.path.splitext(filename)[0]}_{i}.png"
    fig.savefig(filename, dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show(block=False)

def save_slicewise_imgs(image_in, output_dir, data_type=None):
    """Save slice-wise images of of a 3d/4d montage"""
    for an_img in range(image_in.shape[2]):
        out_file = (f"img_slice_{an_img}.png" if data_type is None
                    else f"{data_type}_img_slice_{an_img}.png")
        full_file_name = os.path.join(output_dir, out_file)
        io.imsave(full_file_name, img_as_ubyte(image_in[:,:,an_img,:]), check_contrast=False)

def overlay_images(norm_mr_image, defect_array):
    """
    Overlay 3D RGB image array on each 2D slice of 3D MRI image array and save as PNG image file.

    Args:
        mr_image (numpy.ndarray): normalized MRI image array of shape (height, width, depth).
        defect_array (np.ndarray): Array of shape (height, width, depth) with integer values
        representing different defect types (0 normal, 1 incomplete, 2 complete, and 4 hyperintense)
        mask_image (np.ndarray): Binary mask shape(height, width, depth) for normalizing the image.
        output_dir (str): Directory to save the output PNG image files.
        data_type (str): type of MR image data; Options: None (default), Non_corr, N4_corr, FA_corr
    Returns:
        4D array of Non_corr MR image overlayed by defects
    """
    defect_types = {0: 'Normal', 1: 'Defect', 4: 'Hyper'}
    #defect_dict = {defect_types[i]: (defect_array == i).astype(int) for i in defect_types}
    defect_dict = {defect_type: (defect_array == index).astype(int)
                    for index, defect_type in defect_types.items()}
    ##Define a color map
    cmap_list = [[[0, 0, 0], [1, 0, 0]], [[0, 0, 0], [0, 0, 1]], [[0, 0, 0], [1, 1, 1]]]
    cm_defect, cm_hyper, cm_normal = [ListedColormap(cm) for cm in cmap_list]
    # This part uses transpose to change the dimensions of each array,
    # applies color map function (cm_defect, cm_hyper, and cm_normal) to each.
    # Finally, the result is transposed again to obtain the desired shape of the output array.
    array_defect=cm_defect(defect_dict['Defect'].transpose(2,0,
                                                            1))[..., :3].transpose(1,2,0,3)
    array_hyper=cm_hyper(defect_dict['Hyper'].transpose(2,0,1))[..., :3].transpose(1,2,0,3)
    array_normal=cm_normal(defect_dict['Normal'].transpose(2,0,1))[..., :3].transpose(1,2,0,3)
    array4d_rgb = np.stack([norm_mr_image] * 3, axis=-1)
    defect_array_rgb = np.sum([array_defect, 4 * array_hyper], axis=0)
    defect_overlay_4d = ((defect_array_rgb + (array_normal * array4d_rgb)) * 255).astype(np.uint8)
    return defect_overlay_4d

def load_select_defects(subj_dir, f_name, sub_dir1=None, sub_dir2=None):
    """Select defects from input array"""
    if sub_dir1 is not None and sub_dir2 is not None:
        analysis_arr = np.load(os.path.join(subj_dir, sub_dir1, sub_dir2, f_name))
        defect_mask = analysis_arr == 1
        return defect_mask
    # if sub_dir1 is not None:
    #     analysis_arr, _ = process_nifti(os.path.join(subj_dir, sub_dir1, f_name))
    #     defect_mask = analysis_arr == 3
    #     return defect_mask
    analysis_arr, _ = process_nifti(os.path.join(subj_dir, f_name))
    return analysis_arr

#%%Load images, masks, and defect arry and overlay defect array onto the image
##Path to directory containing subject folders # age_matched_spiral_healthy
PARENT_DIR = (r"C:\Users\HUSDQ4\OneDrive - cchmc\cincy_work\all_projects_data_work"
              r"\vdp_analysis\analysis_comparisons\IRC740H_visit1_spiral_cf") 
###Specify options
ANALYSIS_MODE = "Batch" ##Single or Batch modes
SUBJECT_ID = "IRC740H-028"
CORR = "FA" # N4, FA
SUBJ_DIR_NAMES = [r'^IRC740H-\d{3}$', r'^IRC740H-\d{3}c$', r'^ILD-HC-\d{3}$'] #

## Single subject analysis mode
if ANALYSIS_MODE == "Single":
    SUBJECT_DIR = os.path.join(PARENT_DIR, SUBJECT_ID)
    ## Load defect arrays, select defects, and load reader defects
    md_ = load_select_defects(SUBJECT_DIR, f"{CORR}_corr_glb-median_defect_array.npy",
                                    f"{CORR}_corr_vdp_analysis", "glb_median_analysis")
    mn_ = load_select_defects(SUBJECT_DIR, f"{CORR}_corr_lb-mean_defect_array.npy",
                                    f"{CORR}_corr_vdp_analysis", "lb_mean_analysis")
    pct_ = load_select_defects(SUBJECT_DIR,f"{CORR}_corr_glb-percentile_defect_array.npy",
                                    f"{CORR}_corr_vdp_analysis", "glb_percentile_analysis")
    ak_ = load_select_defects(SUBJECT_DIR,f"{CORR}_corr_adaptive-kmeans_defect_array.npy",
                                    f"{CORR}_corr_vdp_analysis", "kmeans_adaptive_analysis")
    hk_ = load_select_defects(SUBJECT_DIR, f"{CORR}_corr_hierarchical-kmeans_defect_array.npy",
                                    f"{CORR}_corr_vdp_analysis", "kmeans_hierarchical_analysis")
    th_ = load_select_defects(SUBJECT_DIR, f"{CORR}_corr_thresholding_defect_array.npy",
                                    f"{CORR}_corr_vdp_analysis", "thresholding_analysis")
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
    ###Load Riaz's segmented defect mask
    rh_defects, _ = process_nifti(os.path.join(SUBJECT_DIR, "img_defect_mask_rh.nii.gz"))
    rdr_arr = rh_defects + DEFECTS == 2
    MSK_file = load_select_defects(SUBJECT_DIR, "img_ventilation_mask.nii.gz")
    IMG_file =  load_select_defects(SUBJECT_DIR, "img_ventilation.nii.gz")
    # IMG_masked = IMG_file * MSK_file
    IMG_file /= np.max(IMG_file)
    rdr_ = rdr_arr * MSK_file
    defects_dict = {'reader': rdr_, 'hierarchical': hk_, 'adaptive': ak_,
                    'percentile': pct_, 'mean': mn_, 'thresholding': th_,
                    'median': md_}
    for arr_name, arr in defects_dict.items():
        print(arr_name, np.count_nonzero(arr))
        defect_arr_4d = overlay_images(IMG_file, arr)
        defect_montage_arr = create_montage_array(defect_arr_4d, 'all', MSK_file)
        defect_montage_fig = montage_plot_4d(defect_montage_arr)
        ##Create directory for saving overlay images
        OUT_DIR = os.path.join(SUBJECT_DIR, "rdr_calcs_defect_overlays",
                            "defect_overlay_imgs")
        os.makedirs(OUT_DIR, exist_ok=True)
        if arr_name == 'reader':
            save_image(defect_montage_fig, OUT_DIR, data_type="Orig",
                        im_type=f"{arr_name}_defect_overlay")
            save_slicewise_imgs(defect_arr_4d, OUT_DIR,
                        data_type=f'{arr_name}_defect_overlay')
        else:
            save_image(defect_montage_fig, OUT_DIR, data_type=f"{CORR}",
                        im_type=f"{arr_name}_defect_overlay")
            save_slicewise_imgs(defect_arr_4d, OUT_DIR,
                        data_type=f'{CORR}_{arr_name}_defect_overlay')

## Batch analysis mode
elif ANALYSIS_MODE == "Batch":
    for dirpath, dirnames, filenames in os.walk(PARENT_DIR):
        for dirname in dirnames:
            # Check if any of the patterns match the dirname
            if any(re.match(pattern, dirname) for pattern in SUBJ_DIR_NAMES):
                print(f"\nSubject folder name: {dirname}\n")
                SUBJECT_DIR = os.path.join(PARENT_DIR, dirname)
                ## Load defect arrays, select defects, and load reader defects
                md_ = load_select_defects(SUBJECT_DIR,
                                        f"{CORR}_corr_glb-median_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "glb_median_analysis")
                mn_ = load_select_defects(SUBJECT_DIR,
                                        f"{CORR}_corr_lb-mean_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "lb_mean_analysis")
                pct_ = load_select_defects(SUBJECT_DIR,
                                        f"{CORR}_corr_glb-percentile_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "glb_percentile_analysis")
                ak_ = load_select_defects(SUBJECT_DIR,
                                        f"{CORR}_corr_adaptive-kmeans_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "kmeans_adaptive_analysis")
                hk_ = load_select_defects(SUBJECT_DIR,
                                        f"{CORR}_corr_hierarchical-kmeans_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "kmeans_hierarchical_analysis")
                th_ = load_select_defects(SUBJECT_DIR,
                                        f"{CORR}_corr_thresholding_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "thresholding_analysis")
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
                ###Load Riaz's segmented defect mask
                rh_defects, _ = process_nifti(os.path.join(SUBJECT_DIR,
                                                "img_defect_mask_rh.nii.gz"))
                rdr_arr = rh_defects + DEFECTS == 2
                MSK_file = load_select_defects(SUBJECT_DIR, "img_ventilation_mask.nii.gz")
                IMG_file =  load_select_defects(SUBJECT_DIR, "img_ventilation.nii.gz")
                # IMG_masked = IMG_file * MSK_file
                IMG_file /= np.max(IMG_file)
                rdr_ = rdr_arr * MSK_file
                print(f"reader: {np.count_nonzero(rdr_)}\n"
                      f"hierarchical: {np.count_nonzero(hk_)}\tadaptive: {np.count_nonzero(ak_)}\t"
                      f"percentile: {np.count_nonzero(pct_)}\nmean: {np.count_nonzero(mn_)}\t"
                      f"thresholding': {np.count_nonzero(th_)}\tmedian: {np.count_nonzero(md_)}\n")
                defects_dict = {'reader': rdr_, 'hierarchical': hk_, 'adaptive': ak_,
                                'percentile': pct_, 'mean': mn_, 'thresholding': th_,
                                'median': md_}
                for arr_name, arr in defects_dict.items():
                    print(arr_name, np.count_nonzero(arr))
                    defect_arr_4d = overlay_images(IMG_file, arr)
                    defect_montage_arr = create_montage_array(defect_arr_4d, 'all', MSK_file)
                    defect_montage_fig = montage_plot_4d(defect_montage_arr)
                    ##Create directory for saving overlay images
                    OUT_DIR = os.path.join(SUBJECT_DIR, "rdr_calcs_defect_overlays",
                                        "defect_overlay_imgs")
                    os.makedirs(OUT_DIR, exist_ok=True)
                    if arr_name == 'reader':
                        save_image(defect_montage_fig, OUT_DIR, data_type="Orig",
                                    im_type=f"{arr_name}_defect_overlay")
                        save_slicewise_imgs(defect_arr_4d, OUT_DIR,
                                    data_type=f'{arr_name}_defect_overlay')
                    else:
                        save_image(defect_montage_fig, OUT_DIR, data_type=f"{CORR}",
                                    im_type=f"{arr_name}_defect_overlay")
                        save_slicewise_imgs(defect_arr_4d, OUT_DIR,
                                    data_type=f'{CORR}_{arr_name}_defect_overlay')

#%%
