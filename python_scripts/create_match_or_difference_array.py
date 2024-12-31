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
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from skimage import io, img_as_ubyte

##Path to directory containing subject folders
PARENT_DIR = (r"C:\Users\HUSDQ4\OneDrive - cchmc\cincy_work\all_projects_data_work"
              r"\vdp_analysis\analysis_comparisons\IRC740H_visit1_spiral_cf")

#%% Define necessary functions
def process_nifti(in_array):
    """Load and process a nifti image"""
    nifti_read = nib.load(in_array)
    nifti_data = nifti_read.get_fdata()
    reorient_array = np.flip(np.rot90(nifti_data), 0)
    return reorient_array, nifti_read

def calc_dice(im1, im2, empty_score=1.0):
    """
    Computes the Dice coefficient, a measure of set similarity.
    Parameters
    ----------
    im1 : array-like, bool
        Any array of arbitrary size. If not boolean, will be converted.
    im2 : array-like, bool
        Any other array of identical size. If not boolean, will be converted.
    Returns
    -------
    dice : float
        Dice coefficient as a float on range [0,1].
        Maximum similarity = 1
        No similarity = 0
        Both are empty (sum eq to zero) = empty_score
        
    Notes
    -----
    The order of inputs for `dice` is irrelevant. The result will be
    identical if `im1` and `im2` are switched.
    """
    im1 = np.asarray(im1).astype(bool)
    im2 = np.asarray(im2).astype(bool)
    if im1.shape != im2.shape:
        raise ValueError("Shape mismatch: arrays must have the same shape.")
    im_sum = im1.sum() + im2.sum()
    if im_sum == 0:
        return empty_score
    # #Compute Dice coefficient
    intersection = np.logical_and(im1, im2)
    dice_coeff = 2. * intersection.sum() / im_sum
    return dice_coeff

def calc_jaccard(im1, im2, empty_score=1.0):
    """
    Computes the Jaccard index, a measure of set similarity.
    
    Parameters
    ----------
    im1 : array-like, bool
        Any array of arbitrary size. If not boolean, will be converted.
    im2 : array-like, bool
        Any other array of identical size. If not boolean, will be converted.
        
    Returns
    -------
    jaccard : float
        Jaccard index as a float on range [0,1].
        Maximum similarity = 1
        No similarity = 0
        Both are empty (sum equals zero) = empty_score
        
    Notes
    -----
    The order of inputs for `jaccard` is irrelevant. The result will be
    identical if `im1` and `im2` are switched.
    """
    im1 = np.asarray(im1).astype(bool)
    im2 = np.asarray(im2).astype(bool)
    if im1.shape != im2.shape:
        raise ValueError("Shape mismatch: arrays must have the same shape.")
    im_sum = im1.sum() + im2.sum()
    if im_sum == 0:
        return empty_score
    # Compute Jaccard index
    intersection = np.logical_and(im1, im2)
    union = np.logical_or(im1, im2)
    jaccard_idx = intersection.sum() / union.sum()
    return jaccard_idx

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
        montage_output (ndarray): Numpy array ready to make montage.
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
        for ii in slices[1:]:
            montage_output = np.hstack((montage_output, montage_in[:, :, ii]))
    elif np.ndim(montage_in) == 4:
        montage_output = montage_in[:, :, slices[0], :]
        for ii in slices[1:]:
            montage_output = np.hstack((montage_output, montage_in[:, :, ii, :]))
    return montage_output

def plot_binning_montage(binned_cmap_montage):
    """Montage for generalized linear binning maps"""
    # #Define color map with 0 values as black and other values assigned to colors
    my_cmap = ListedColormap(['#000000', '#FFEA00','#00FFFF','#FF00FF'])
    # #Create the montage with a fixed colorbar i.e. colors in cbar are fixed,
    # #regardless of the fact that they're present in plot or not.
    bin_fig, ax = plt.subplots(figsize=(16, 3), dpi=300)
    ax.imshow(binned_cmap_montage, cmap=my_cmap, alpha=1)
    plt.axis('off')
    sm = fixed_cbar()
    # # Create an axis for the colorbar using make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1.5%", pad="0.1%")
    cbar = plt.colorbar(sm, cax=cax)
    cbar.set_ticks([0, 0.5, 1.5, 2.5])
    cbar.set_ticklabels(['','FN', 'FP', 'TP'])
    cbar.ax.tick_params(size=0)
    plt.tight_layout()
    plt.show(block=False)
    return bin_fig

def fixed_cbar():
    """Generate a colorbar with fixed colors matching glb colors"""
    # # #Define the fixed colors and bounds for the color bar
    cbar_colors = ListedColormap(['#FFEA00','#00FFFF','#FF00FF'])
    cbar_bounds = [0, 1, 2, 3]
    # # #Create a ScalarMappable with the colormap
    s_map = ScalarMappable(cmap=cbar_colors, norm=BoundaryNorm(cbar_bounds, cbar_colors.N))
    return s_map

def save_image(figg, save_path, data_type=None, im_type=None):
    """To save image in the provided path
        save_path (str): path to directory for saving the image.
        data_type (str): Optional type of data being saved
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
    figg.savefig(filename, dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show(block=False)

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
        defect_array (np.ndarray): reader and any other defect overlay arr.
        integer valuesrepresenting different defect overlay types:
        (0 normal, 1 false negative, 2 false positive, and 3 true positive)
    Returns:
        4D array of Non_corr MR image overlayed by defects
    """
    defect_types = {0: 'Nrml', 1: 'FN', 2: 'FP', 3: 'TP'}
    defect_dict = {defect_type: (defect_array == index).astype(int)
                    for index, defect_type in defect_types.items()}
    ##Define a color map
    cmap_list = [[[0, 0, 0], '#FFEA00'], [[0, 0, 0], '#00FFFF'], [[0, 0, 0], '#FF00FF'],
                 [[0, 0, 0], [1, 1, 1]]]
    cm_fn, cm_fp, cm_tp, cm_nrml = [ListedColormap(cm) for cm in cmap_list]
    # This part uses transpose to change the dimensions of each array,
    # applies color map function (cm_fn, cm_fp, cm_tp, and cm_normal) to each.
    # Finally, the result is transposed again to obtain the desired shape of the output array.
    array_fn=cm_fn(defect_dict['FN'].transpose(2,0,1))[..., :3].transpose(1,2,0,3)
    array_fp=cm_fp(defect_dict['FP'].transpose(2,0,1))[..., :3].transpose(1,2,0,3)
    array_tp=cm_tp(defect_dict['TP'].transpose(2,0,1))[..., :3].transpose(1,2,0,3)
    array_normal=cm_nrml(defect_dict['Nrml'].transpose(2,0,1))[..., :3].transpose(1,2,0,3)
    defect_array_rgb = np.sum([array_fn, 2*array_fp, 3*array_tp], axis=0)
    defect_overlay_4d = ((defect_array_rgb + (array_normal * np.stack([norm_mr_image]*3,axis=-1))
                          )* 255).astype(np.uint8)
    return defect_overlay_4d

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

def calc_dice_match_diff_plot(rdr_defects, calc_defects, corr_anlys,
                              subj_dir=None, msk=None):
    """Calculate dice coefficient, find match and difference in reader defect mask
        and various analysis methods defect masks, save and print results"""
    dice_ = calc_dice(rdr_defects, calc_defects)
    jaccard_ = calc_jaccard(rdr_defects, calc_defects)
    rdr_and_calc = rdr_defects + (2 * calc_defects) ##Only reader=1, calc=2, common=3
    fn = (np.sum((rdr_defects==1) & (calc_defects==0) & (msk == 1))/np.sum(msk==1))*100
    fp = (np.sum((rdr_defects==0) & (calc_defects==1) & (msk == 1))/np.sum(msk==1))*100
    tp = (np.sum((rdr_defects==1) & (calc_defects==1) & (msk == 1))/np.sum(msk==1))*100
    tn = (np.sum((rdr_defects==0) & (calc_defects==0) & (msk == 1))/np.sum(msk==1))*100
    print(f"Rdr Defects missed by {corr_anlys[0]} {corr_anlys[1]}: {fn}%\n")
    print(f"Rdr Defects overestimated by {corr_anlys[0]} {corr_anlys[1]}: {fp}%\n")
    print(f"Rdr Defects common with {corr_anlys[0]} {corr_anlys[1]}: {tp}%\n")
    print(f"Rdr Vent common with {corr_anlys[0]} {corr_anlys[1]}: {tn}%\n")
    print(f"Dice: {dice_}\t Jaccard: {jaccard_}")
    data_file = os.path.join(os.path.dirname(subj_dir),
                             f"rdr_{corr_anlys[0]}_{corr_anlys[1]}_dice_overlay.csv")
    write_to_csv(data_file, ["Subject_id", "Dice", "Jaccard", "pFN", "pFP", "pTP", "pTN"],
                 [os.path.basename(subj_dir), dice_, jaccard_, fn, fp, tp, tn])
    if msk is not None:
        plot_binning_montage(create_montage_array(rdr_and_calc, msk_file=msk))
    else:
        plot_binning_montage(create_montage_array(rdr_and_calc, slices='all'))
    return rdr_and_calc

def master_calcs_plots(reader, auto_calc, corr_anlys, subj_dir, img_msk):
    """Calculate/plot dice coefficient and overlay b/w reader and anto-analysis methods"""
    defects_diff_arr = calc_dice_match_diff_plot(reader, auto_calc,
                                            corr_anlys, subj_dir, img_msk[1])
    # defects_overlay_4d = overlay_images(img_msk[0], defects_diff_arr)
    # defects_overlay_montage = create_montage_array(defects_overlay_4d, 'all', img_msk[1])
    # overlay_fig = plot_binning_montage(defects_overlay_montage)
    # outdir = os.path.join(subj_dir, "rdr_calcs_defect_overlays")
    # os.makedirs(outdir, exist_ok=True)
    # save_image(overlay_fig, outdir,
    #            f"rdr_{corr_anlys[0]}_{corr_anlys[1]}", "defects_ovrlay")
    # save_slicewise_imgs(defects_overlay_4d, outdir, f"rdr_{corr_anlys[0]}_{corr_anlys[1]}_ovrlay")

#%% Loading defect data, calculating dice coefficient
ANALYSIS_MODE = "Batch" ##Single or Batch modes
SUBJECT_ID = "IRC740H-003"
CORR = "N4" # N4, FA
SUBJ_DIR_NAMES = [r'^IRC740H-\d{3}$', r'^IRC740H-\d{3}c$', r'^ILD-HC-\d{3}$'] #
## Single subject analysis mode
if ANALYSIS_MODE == "Single":
    SUBJECT_DIR = os.path.join(PARENT_DIR, SUBJECT_ID)
    ## Load defect arrays, select defects, and load reader defects
    m_defects = load_select_defects(SUBJECT_DIR, f"{CORR}_corr_glb-median_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "glb_median_analysis")
    pct_defects = load_select_defects(SUBJECT_DIR,f"{CORR}_corr_glb-percentile_defect_array.npy",
                                    f"{CORR}_corr_vdp_analysis", "glb_percentile_analysis")
    ak_defects = load_select_defects(SUBJECT_DIR,f"{CORR}_corr_adaptive-kmeans_defect_array.npy",
                                    f"{CORR}_corr_vdp_analysis", "kmeans_adaptive_analysis")
    hk_defects = load_select_defects(SUBJECT_DIR,
                                    f"{CORR}_corr_hierarchical-kmeans_defect_array.npy",
                                    f"{CORR}_corr_vdp_analysis", "kmeans_hierarchical_analysis")
    th_defects = load_select_defects(SUBJECT_DIR, f"{CORR}_corr_thresholding_defect_array.npy",
                                    f"{CORR}_corr_vdp_analysis", "thresholding_analysis")
    rdr_rh = load_select_defects(SUBJECT_DIR,"img_defect_mask_rh.nii.gz")
    # rdr_jp = load_select_defects(SUBJECT_DIR,"img_defect_mask_jp.nii.gz")
    msk_arr = load_select_defects(SUBJECT_DIR, "img_ventilation_mask.nii.gz")
    img_arr = load_select_defects(SUBJECT_DIR, "img_ventilation.nii.gz")
    img_arr /= np.max(img_arr*msk_arr)
    rdr_arr = rdr_rh*msk_arr # +(rdr_jp*msk_arr) == 2
    img_n_msk = [img_arr, msk_arr]
    master_calcs_plots(rdr_arr, m_defects, [CORR, "glb-median"], SUBJECT_DIR, img_n_msk)
    master_calcs_plots(rdr_arr, pct_defects, [CORR, "glb-percentile"], SUBJECT_DIR, img_n_msk)
    master_calcs_plots(rdr_arr, ak_defects, [CORR, "adaptive-kmeans"], SUBJECT_DIR, img_n_msk)
    master_calcs_plots(rdr_arr, hk_defects, [CORR, "hierarchical-kmeans"], SUBJECT_DIR, img_n_msk)
    master_calcs_plots(rdr_arr, th_defects, [CORR, "thresholding"], SUBJECT_DIR, img_n_msk)
elif ANALYSIS_MODE == "Batch":
    for dirpath, dirnames, filenames in os.walk(PARENT_DIR):
        for dirname in dirnames:
            # Check if any of the patterns match the dirname
            if any(re.match(pattern, dirname) for pattern in SUBJ_DIR_NAMES):
                print(f"\nSubject folder name: {dirname}")
                SUBJECT_DIR = os.path.join(PARENT_DIR, dirname)
                ## Load defect arrays, select defects, and load reader defects
                m_defects = load_select_defects(SUBJECT_DIR,
                                        f"{CORR}_corr_glb-median_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "glb_median_analysis")
                pct_defects = load_select_defects(SUBJECT_DIR,
                                        f"{CORR}_corr_glb-percentile_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "glb_percentile_analysis")
                ak_defects = load_select_defects(SUBJECT_DIR,
                                        f"{CORR}_corr_adaptive-kmeans_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "kmeans_adaptive_analysis")
                hk_defects = load_select_defects(SUBJECT_DIR,
                                        f"{CORR}_corr_hierarchical-kmeans_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "kmeans_hierarchical_analysis")
                th_defects = load_select_defects(SUBJECT_DIR,
                                        f"{CORR}_corr_thresholding_defect_array.npy",
                                        f"{CORR}_corr_vdp_analysis", "thresholding_analysis")
                rdr_rh = load_select_defects(SUBJECT_DIR,"img_defect_mask_rh.nii.gz")
                # rdr_jp = load_select_defects(SUBJECT_DIR,"img_defect_mask_jp.nii.gz")
                msk_arr = load_select_defects(SUBJECT_DIR, "img_ventilation_mask.nii.gz")
                img_arr = load_select_defects(SUBJECT_DIR, "img_ventilation.nii.gz")
                img_arr /= np.max(img_arr*msk_arr)
                rdr_arr = rdr_rh*msk_arr # (rdr_rh*msk_arr)+(rdr_jp*msk_arr) == 2
                img_n_msk = [img_arr, msk_arr]
                master_calcs_plots(rdr_arr, m_defects, [CORR, "glb-median"],
                                   SUBJECT_DIR, img_n_msk)
                master_calcs_plots(rdr_arr, pct_defects, [CORR, "glb-percentile"],
                                   SUBJECT_DIR, img_n_msk)
                master_calcs_plots(rdr_arr, ak_defects, [CORR, "adaptive-kmeans"],
                                   SUBJECT_DIR, img_n_msk)
                master_calcs_plots(rdr_arr, hk_defects, [CORR, "hierarchical-kmeans"],
                                   SUBJECT_DIR, img_n_msk)
                master_calcs_plots(rdr_arr, th_defects, [CORR, "thresholding"],
                                   SUBJECT_DIR, img_n_msk)

# %%
