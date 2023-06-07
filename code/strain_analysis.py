import glob
from matplotlib import path as mpl_path
import matplotlib.patheffects as pe
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import image_analysis as ia
import numpy as np
from pathlib import Path
from typing import List, Union 


def box_to_mask(mask: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Given a mask (for dimensions) and a box. Will return a mask of the inside of the box."""
    p = mpl_path.Path([(box[0, 0], box[0, 1]), (box[1, 0], box[1, 1]), (box[2, 0], box[2, 1]), (box[3, 0], box[3, 1])])
    new_mask = np.zeros(mask.shape)
    for rr in range(0, mask.shape[0]):
        for cc in range(0, mask.shape[1]):
            new_mask[rr, cc] = p.contains_points([(rr, cc)])
    return new_mask


def shrink_box(box: np.ndarray, shrink_row: float = 0.1, shrink_col: float = 0.1) -> np.ndarray:
    """Given a box and shrink factors. Will make the box smaller by the given fraction."""
    r0, r1, c0, c1 = ia.box_to_bound(box)
    r0_new, r1_new = ia.shrink_pair(r0, r1, shrink_row)
    c0_new, c1_new = ia.shrink_pair(c0, c1, shrink_col)
    box = ia.bound_to_box(r0_new, r1_new, c0_new, c1_new)
    return box


def create_sub_domains(
    mask: np.ndarray,
    *,
    pillar_clip_fraction: float = 0.5,
    shrink_row: float = 0.1,
    shrink_col: float = 0.1,
    tile_dim_pix: int = 40,
    num_tile_row: int = 5,
    num_tile_col: int = 3,
    tile_style: int = 1,
    clip_columns: bool = True,
    clip_rows: bool = False,
    manual_sub: bool = False,
    sub_extents: List = None
) -> List:
    """Given a mask and sub-domain parameters. Will return a list of sub-domains define by 4 box coordinates.
    tile_style = 1 will fit as many square tiles of the given tile_dim_pix size in a grid
    tyle_style = 2 will create a num_tile_row x num_tile_col sized grid and adjust the tile_dim_pix as need
    """
    # clip pillars from the mask
    mask_removed_pillars = ia.remove_pillar_region(mask, pillar_clip_fraction,clip_columns, clip_rows)

    if manual_sub:
        r0_user, r1_user, c0_user, c1_user = sub_extents
        user_box = ia.bound_to_box(r0_user, r1_user, c0_user, c1_user)
        first_mask = box_to_mask(mask_removed_pillars, user_box)
        first_box_mask = ia.mask_to_box(first_mask)
        box_mask = shrink_box(first_box_mask, shrink_row, shrink_col)
    else:    
        # compute overall box
        first_box_mask = ia.mask_to_box(mask_removed_pillars)
        box_mask = shrink_box(first_box_mask, shrink_row, shrink_col)
    # tile sub-domains
    r0, r1, c0, c1 = ia.box_to_bound(box_mask)
    if tile_style == 1:
        num_tile_row = int(np.floor((r1 - r0) / tile_dim_pix))
        num_tile_col = int(np.floor((c1 - c0) / tile_dim_pix))
    elif tile_style == 2:
        tile_dim_pix = np.min([np.floor((r1 - r0) / num_tile_row), np.floor((c1 - c0) / num_tile_col)])
    # sub-divide into sub-domains
    shrink_row = 1.0 - num_tile_row * tile_dim_pix / (r1 - r0)
    shrink_col = 1.0 - num_tile_col * tile_dim_pix / (c1 - c0)
    box_mask_grid = shrink_box(box_mask, shrink_row, shrink_col)
    r0_box, _, c0_box, _ = ia.box_to_bound(box_mask_grid)
    tile_box_list = []
    for rr in range(0, num_tile_row):
        for cc in range(0, num_tile_col):
            tile_box = ia.bound_to_box(r0_box + rr * tile_dim_pix, r0_box + (rr + 1) * tile_dim_pix, c0_box + cc * tile_dim_pix, c0_box + (cc + 1) * tile_dim_pix)
            tile_box_list.append(tile_box)
    return tile_box_list, tile_dim_pix, num_tile_row, num_tile_col


def isolate_sub_domain_markers(tracker_row: np.ndarray, tracker_col: np.ndarray, sd_box: np.ndarray) -> List:
    """Given tracker row and column arrays and sub-domain box. Will return markers inside sub-domain."""
    
    all_num_sub_pts = [] 

    sd_tracker_row = []
    sd_tracker_col = []
    for jj in range(0, tracker_row.shape[0]):
        rr = tracker_row[jj, 0]
        cc = tracker_col[jj, 0]
        if ia.is_in_box(sd_box, rr, cc):
            sd_tracker_row.append(tracker_row[jj, :])
            sd_tracker_col.append(tracker_col[jj, :])
    sd_tracker_row = np.asarray(sd_tracker_row)
    sd_tracker_col = np.asarray(sd_tracker_col)
    num_sub_pts = len(sd_tracker_row)
    all_num_sub_pts.append(num_sub_pts)

    return sd_tracker_row, sd_tracker_col, all_num_sub_pts


def compute_F_from_Lambda_mat(Lambda_0: np.ndarray, Lambda_t: np.ndarray) -> np.ndarray:
    """Given Lambda (i.e., vectors connecting fiducial marker positions) matricies in positions 0 and t.
    Lambda matricies are units 2 x number of points.
    Will return the average deformation gradient F."""
    term_1 = np.dot(Lambda_t, np.transpose(Lambda_0))
    term_2 = np.linalg.inv(np.dot(Lambda_0, np.transpose(Lambda_0)))
    F = np.dot(term_1, term_2)
    # add in error handling if issues arise here --
    return F


def compute_Lambda_from_pts(row_pos: np.ndarray, col_pos: np.ndarray) -> np.ndarray:
    num_pts = row_pos.shape[0]
    pts_row_col = np.array([row_pos,col_pos])
    ii, jj = np.triu_indices(num_pts, k=1)
    Lambda_mat = pts_row_col[:,ii] - pts_row_col[:,jj]
    return Lambda_mat


def compute_sub_domain_strain(sd_tracker_row: np.ndarray, sd_tracker_col: np.ndarray) -> np.ndarray:
    """Given tracking point positions. Will return F at every frame with the first frame as the reference."""
    sd_F_list = []
    num_frames = sd_tracker_row.shape[1]
    Lambda_0 = compute_Lambda_from_pts(sd_tracker_row[:, 0], sd_tracker_col[:, 0])
    for kk in range(0, num_frames):
        Lambda_t = compute_Lambda_from_pts(sd_tracker_row[:, kk], sd_tracker_col[:, kk])
        F = compute_F_from_Lambda_mat(Lambda_0, Lambda_t)
        sd_F_list.append(F.reshape((-1, 1)))
    sd_F_arr = np.asarray(sd_F_list)[:, :, 0]
    return sd_F_arr


def get_box_center(box: np.ndarray) -> Union[float, int]:
    """Given a box. Will return the center."""
    center_row = np.mean(box[:, 0])
    center_col = np.mean(box[:, 1])
    return center_row, center_col


def compute_sub_domain_position(sd_tracker_row: List, sd_tracker_col: List, sd_box: np.ndarray) -> List:
    """Given sub-domain displacements and sub-domain box. Will return the box center displacement at every step."""
    sd_row = []
    sd_col = []
    center_row, center_col = get_box_center(sd_box)
    num_frames = sd_tracker_row.shape[1]
    for kk in range(0, num_frames):
        row_pos = center_row + np.mean(sd_tracker_row[:, kk] - sd_tracker_row[:, 0])
        col_pos = center_col + np.mean(sd_tracker_col[:, kk] - sd_tracker_col[:, 0])
        sd_row.append(row_pos)
        sd_col.append(col_pos)
    return np.asarray(sd_row), np.asarray(sd_col)


def compute_sub_domain_position_strain_all(tracker_row_all: List, tracker_col_all: List, sd_box_list: List) -> List:
    """Given all tracker rows, columns, and sub-domain box domains. Will return the markers within each sub-domain box."""
    sub_domain_F_all = []
    sub_domain_row_all = []
    sub_domain_col_all = []
    
    all_pts_sd = []
    for sd_box in sd_box_list:

        sd_tracker_row_all, sd_tracker_col_all, num_pts_sd = isolate_sub_domain_markers(tracker_row_all, tracker_col_all, sd_box)
        all_pts_sd.append(num_pts_sd)

        sd_tracker_row = sd_tracker_row_all
        sd_tracker_col = sd_tracker_col_all
        sd_F_array = compute_sub_domain_strain(sd_tracker_row, sd_tracker_col)
        sub_domain_F_all.append(sd_F_array)
        sd_row, sd_col = compute_sub_domain_position(sd_tracker_row, sd_tracker_col, sd_box)
        sub_domain_row_all.append(sd_row)
        sub_domain_col_all.append(sd_col)
    return sub_domain_F_all, sub_domain_row_all, sub_domain_col_all, all_pts_sd

def format_F_for_save(sub_domain_F_all: List) -> List:
    """Given sub_domain_F_all. Will return in a format where F_rr, F_rc, F_cr, and F_cc are formatted for saving."""
    """Subscript r refers to row and c to column: For example F_rr refers to F in the row-row direction."""

    num_subdomains = len(sub_domain_F_all)
    num_frames = sub_domain_F_all[0].shape[0]

    sub_domain_F_rr = np.zeros((num_subdomains, num_frames))
    sub_domain_F_rc = np.zeros((num_subdomains, num_frames))
    sub_domain_F_cr = np.zeros((num_subdomains, num_frames))
    sub_domain_F_cc = np.zeros((num_subdomains, num_frames))  
        
    for jj in range(0, num_subdomains):
        sub_domain_F_rr[jj,:] = sub_domain_F_all[jj][:, 0]
        sub_domain_F_rc[jj,:] = sub_domain_F_all[jj][:, 1]
        sub_domain_F_cr[jj,:] = sub_domain_F_all[jj][:, 2]
        sub_domain_F_cc[jj,:] = sub_domain_F_all[jj][:, 3]
        
    return sub_domain_F_rr, sub_domain_F_rc, sub_domain_F_cr, sub_domain_F_cc


def format_sd_row_col_for_save(sub_domain_row_all: List, sub_domain_col_all: List) -> List:
    """Given sub_domain_row_all and sub_domain_col_all. Will return in a format where sub_domain_row and sub_domain_col are formatted for saving."""

    num_subdomains = len(sub_domain_row_all)
    num_frames = sub_domain_row_all[0].shape[0]
    sub_domain_row = np.zeros((num_subdomains, num_frames))
    sub_domain_col = np.zeros((num_subdomains, num_frames))
    for jj in range(0, num_subdomains):
        sub_domain_row[jj, :] = sub_domain_row_all[jj]
        sub_domain_col[jj, :] = sub_domain_col_all[jj]

    return sub_domain_row, sub_domain_col


def save_sub_domain_strain(folder_path: Path, sub_domain_F_all: List, sub_domain_row_all: List, sub_domain_col_all: List, strain_sub_domain_info: np.ndarray, *, fname: str = "") -> List:
    """Given results of sub domain strain computation. Will save the strain results."""
    new_path = ia.create_folder(folder_path, "results")

    sub_domain_F_rr_all, sub_domain_F_rc_all, sub_domain_F_cr_all, sub_domain_F_cc_all = format_F_for_save(sub_domain_F_all)
    sub_domain_row_all, sub_domain_col_all = format_sd_row_col_for_save(sub_domain_row_all, sub_domain_col_all)
    saved_paths = []

    # save the sub domain positions
    file_path = new_path.joinpath("strain_" + fname + "_y.txt").resolve()
    saved_paths.append(file_path)
    np.savetxt(str(file_path), sub_domain_row_all)
    file_path = new_path.joinpath("strain_" + fname + "_x.txt").resolve()
    saved_paths.append(file_path)
    np.savetxt(str(file_path), sub_domain_col_all)
    # save the components of the sub domain deformation gradient
    file_path = new_path.joinpath("strain_" + fname + "_Fyy.txt").resolve()
    saved_paths.append(file_path)
    np.savetxt(str(file_path), sub_domain_F_rr_all)
    file_path = new_path.joinpath("strain_" + fname + "_Fyx.txt").resolve()
    saved_paths.append(file_path)
    np.savetxt(str(file_path), sub_domain_F_rc_all)
    file_path = new_path.joinpath("strain_" + fname + "_Fxy.txt").resolve()
    saved_paths.append(file_path)
    np.savetxt(str(file_path), sub_domain_F_cr_all)
    file_path = new_path.joinpath("strain_" + fname + "_Fxx.txt").resolve()
    saved_paths.append(file_path)
    np.savetxt(str(file_path), sub_domain_F_cc_all)
    # save the information about the sub domains (specifically: num_row_tile, num_col_tile, tile_pix_dim)
    # maybe also save information about rotation?
    file_path = new_path.joinpath("strain_" + fname + "_sub_domain_info.txt").resolve()
    np.savetxt(str(file_path), strain_sub_domain_info)
    saved_paths.append(file_path)
    return saved_paths


def load_sub_domain_strain(folder_path: Path, fname="") -> List:
    """Given folder path. Will load strain results. If there are none, will return an error."""
    res_folder_path = folder_path.joinpath("results").resolve()
    if res_folder_path.exists() is False:
        raise FileNotFoundError("tracking results are not present -- therefore strain results must not be present either")
    file_list = glob.glob(str(res_folder_path) + "/*strain*")
    if len(file_list) == 0:
        raise FileNotFoundError("strain results are not present")

    sub_domain_F_rr_all = np.loadtxt(str(res_folder_path) + "/strain_" + fname + "_Fyy.txt")
    sub_domain_F_rc_all = np.loadtxt(str(res_folder_path) + "/strain_" + fname + "_Fyx.txt")
    sub_domain_F_cr_all = np.loadtxt(str(res_folder_path) + "/strain_" + fname + "_Fxy.txt")
    sub_domain_F_cc_all = np.loadtxt(str(res_folder_path) + "/strain_" + fname + "_Fxx.txt")
    sub_domain_row_all = np.loadtxt(str(res_folder_path) + "/strain_" + fname + "_y.txt")
    sub_domain_col_all = np.loadtxt(str(res_folder_path) + "/strain_" + fname + "_x.txt")
    strain_info = np.loadtxt(str(res_folder_path) + "/strain_" + fname + "_sub_domain_info.txt")
    return sub_domain_F_rr_all, sub_domain_F_rc_all, sub_domain_F_cr_all, sub_domain_F_cc_all, sub_domain_row_all, sub_domain_col_all, strain_info


def F_to_E(F_rr: float, F_rc: float, F_cr: float, F_cc: float) -> float:
    """Given F, will compute E_cc (E_column_column) for visualization."""
    F_arr = np.asarray([[F_rr, F_rc], [F_cr, F_cc]])
    C = np.dot(F_arr.T, F_arr)
    E = 0.5 * (C - np.eye(2))
    return E[1, 1], E[1,0], E[0,0]


def F_to_E_all(sub_domain_F_rr_all: List, sub_domain_F_rc_all: List, sub_domain_F_cr_all: List, sub_domain_F_cc_all: List) -> List:

    num_sub_domains = sub_domain_F_rr_all.shape[0]
    num_frames = sub_domain_F_rr_all.shape[1]
   
    Ecc_arr = np.zeros((num_sub_domains, num_frames))
    Ecr_arr = np.zeros((num_sub_domains, num_frames))
    Err_arr = np.zeros((num_sub_domains, num_frames))
    for jj in range(0, num_sub_domains):
        for ii in range(0, num_frames):
            F_rr = sub_domain_F_rr_all[jj, ii]
            F_rc = sub_domain_F_rc_all[jj, ii]
            F_cr = sub_domain_F_cr_all[jj, ii]
            F_cc = sub_domain_F_cc_all[jj, ii]
            Ecc_arr[jj, ii], Ecr_arr[jj, ii], Err_arr[jj, ii]= F_to_E(F_rr, F_rc, F_cr, F_cc)
    return Ecc_arr, Ecr_arr, Err_arr


def get_text_str(row: int, col: int) -> str:
    """Given the row and the column of a location. Will return the row and column position string."""
    test_str = "%s%i" % (chr(row + 65), col + 1)
    return test_str


def png_sub_domains_numbered(
    folder_path: Path,
    example_tiff: np.ndarray,
    sub_domain_row: List,
    sub_domain_col: List,
    sub_domain_side: Union[float, int],
    num_sd_row: int,
    num_sd_col: int,
    fname: str = "strain_",
    col_map: object = plt.cm.rainbow
) -> List:
    """Given information to visualize the sub-domains. Will plot the subdomains and label them.
    Rows are labeled A, B, C, etc. -- columns are labeled 1, 2, 3, etc. """
    vis_folder_path = ia.create_folder(folder_path, "visualizations")
    pngs_folder_path = ia.create_folder(vis_folder_path, "strain_pngs")
    img_path = pngs_folder_path.joinpath(fname + "sub_domain_key.pdf").resolve()
    plt.figure()
    plt.imshow(example_tiff, cmap=plt.cm.gray)
    sds = sub_domain_side / 2.0
    for cc in range(0, num_sd_col):
        for rr in range(0, num_sd_row):
            idx = rr * num_sd_col + cc
            center_row = sub_domain_row[idx]
            center_col = sub_domain_col[idx]
            text_str = get_text_str(rr, cc)
            plt.text(center_col, center_row, text_str, color=col_map(idx / (num_sd_col * num_sd_row)), fontsize= int(sub_domain_side/12) ,horizontalalignment="center", verticalalignment="center",path_effects=[pe.Stroke(linewidth=1, foreground='k'), pe.Normal()])
            corners_rr = [center_row - sds, center_row - sds, center_row + sds, center_row + sds, center_row - sds]
            corners_cc = [center_col - sds, center_col + sds, center_col + sds, center_col - sds, center_col - sds]
            plt.plot(corners_cc, corners_rr, "k-", linewidth=1)
    plt.axis("off")
    plt.savefig(str(img_path),format='pdf')
    plt.close()
    return img_path


def png_sub_domain_strain_timeseries_all(
    folder_path: Path,
    sub_domain_strain_all: List,
    num_sd_row: int,
    num_sd_col: int,
    output: str,
    col_map: object = plt.cm.rainbow,
    *,
    fname: str = "strain_timeseries_Exx",
    xlabel: str = "frame",
    ylabel: str = "strain Exx"
) -> List:
    """Given strain timeseries. Will plot all timeseries on the same axis."""
    vis_folder_path = ia.create_folder(folder_path, "visualizations")
    main_pngs_folder_path = ia.create_folder(vis_folder_path, "strain_pngs")
    
    if output == "Ecc":
        pngs_folder_path = ia.create_folder(main_pngs_folder_path, "Exx")
    elif output == "Ecr":
        pngs_folder_path = ia.create_folder(main_pngs_folder_path, "Exy")  
    elif output == "Err":
        pngs_folder_path = ia.create_folder(main_pngs_folder_path, "Eyy")           
    path_list = []

    plt.figure(figsize=(15, 10))
    ax = plt.subplot(111)
    sub_domain_strain = sub_domain_strain_all
    for cc in range(0, num_sd_col):
        for rr in range(0, num_sd_row):
            lab = get_text_str(rr, cc)
            idx = rr * num_sd_col + cc
            lab = get_text_str(rr, cc)
            ax.plot(sub_domain_strain[idx, :], label=lab, color=col_map(idx/(num_sd_col * num_sd_row)))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # ax.set_title("beat %i" % (kk))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.25, box.width, box.height * 0.75])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=num_sd_col)
    img_path = pngs_folder_path.joinpath(fname + ".pdf").resolve()
    plt.savefig(str(img_path), format='pdf')
    plt.close()
    path_list.append(img_path)
    return path_list 

def pngs_sub_domain_strain(
    folder_path: Path,
    tiff_list: List,
    sub_domain_row_all: np.ndarray,
    sub_domain_col_all: np.ndarray,
    sub_domain_E_all: np.ndarray,
    sub_domain_side: Union[float, int],
    output: str,
    col_min: Union[float, int] = -0.025,
    col_max: Union[float, int] = 0.025,
    col_map: object = plt.cm.RdBu,
    fname: str = "strain",
    save_eps: bool = False
) -> List:
    """Given sub domain strain results. Will create pngs."""
    vis_folder_path = ia.create_folder(folder_path, "visualizations")
    main_pngs_folder_path = ia.create_folder(vis_folder_path, "strain_pngs")
    
    if output == "Ecc":
        pngs_folder_path = ia.create_folder(main_pngs_folder_path, "Exx")
        E_label = "Exx"

    elif output == "Ecr":
        pngs_folder_path = ia.create_folder(main_pngs_folder_path, "Exy")  
        E_label = "Exy"

    elif output == "Err":
        pngs_folder_path = ia.create_folder(main_pngs_folder_path, "Eyy") 
        E_label = "Eyy"
        
    path_list = []
    
    E = sub_domain_E_all
    num_frames = sub_domain_col_all.shape[1]
    for kk in range(num_frames):
        plt.figure()
        plt.imshow(tiff_list[kk], cmap=plt.cm.gray)
        plt.scatter(sub_domain_col_all[:,kk], sub_domain_row_all[:,kk], c=E[:,kk], s=50, cmap=col_map, vmin=col_min, vmax=col_max)
        # illustrate sub-domain borders
        sds = sub_domain_side / 2.0
        for ii in range(0, sub_domain_col_all.shape[0]):
            corners_rr = [sub_domain_row_all[ii, kk] - sds, sub_domain_row_all[ii, kk] - sds, sub_domain_row_all[ii, kk] + sds, sub_domain_row_all[ii, kk] + sds, sub_domain_row_all[ii, kk] - sds]
            corners_cc = [sub_domain_col_all[ii, kk] - sds, sub_domain_col_all[ii, kk] + sds, sub_domain_col_all[ii, kk] + sds, sub_domain_col_all[ii, kk] - sds, sub_domain_col_all[ii, kk] - sds]
            plt.plot(corners_cc, corners_rr, "k-", linewidth=0.1)
            
        plt.title("frame %i, %s" % (kk, E_label))
        cbar = plt.colorbar()
        cbar.ax.get_yaxis().labelpad = 15
        cbar.set_label("%s strain "%(E_label), rotation=270)
        plt.axis("off")
        path = pngs_folder_path.joinpath("%04d_" % (kk) + fname + ".png").resolve()
        plt.savefig(str(path))
        if save_eps:
            plt.savefig(str(path)[0:-4]+'.eps',format='eps')
        plt.close()
        path_list.append(path)
    return path_list


def create_gif(folder_path: Path, png_path_list: List, fname="sub_domain_strain") -> Path:
    """Given the pngs path list. Will create a gif."""
    img_list = []
    img = plt.imread(png_path_list[0])
    img_r, img_c,_ = img.shape
    fig, ax = plt.subplots(figsize=(img_c/100,img_r/100))
    plt.axis('off') 
    plt.tight_layout(pad=0.08, h_pad=None, w_pad=None, rect=None)
    for pa in png_path_list:
        img = ax.imshow(plt.imread(pa),animated=True)
        img_list.append([img])
    fn_gif = fname + ".gif"
    gif_path = folder_path.joinpath("visualizations").resolve().joinpath(fn_gif).resolve()
    ani = animation.ArtistAnimation(fig, img_list,interval=100)
    ani.save(gif_path,dpi=100)
    return gif_path


def run_sub_domain_strain_analysis(
    folder_path: Path,
    pillar_clip_fraction: float = 0.5,
    shrink_row: float = 0.1,
    shrink_col: float = 0.1,
    tile_dim_pix: int = 40,
    num_tile_row: int = 5,
    num_tile_col: int = 3,
    tile_style: int = 1,
    is_rotated: bool = True,
    clip_columns: bool = True,
    clip_rows: bool = False,
    manual_sub: bool = False,
    sub_extents: List = None,
    *,
    save_fname: str = ""
) -> List:
    """Given a folder path. Will perform strain analysis and save results as text files.
    Note that this function assumes that we have already run the tracking portion of the code."""
    # read images and mask file
    mask_file_path = folder_path.joinpath("masks").resolve().joinpath("tissue_mask.txt").resolve()
    mask = ia.read_txt_as_mask(mask_file_path)

    # load tracking results
    tracker_row_all, tracker_col_all = ia.load_tracking_results(folder_path=folder_path)
    
    # create sub-domains
    sub_domain_box_list, tile_dim_pix, num_tile_row, num_tile_col = create_sub_domains(mask, pillar_clip_fraction=pillar_clip_fraction, shrink_row=shrink_row, shrink_col=shrink_col, tile_dim_pix=tile_dim_pix, num_tile_row=num_tile_row, num_tile_col=num_tile_col, tile_style=tile_style,clip_columns=clip_columns,clip_rows=clip_rows,manual_sub= manual_sub,sub_extents=sub_extents)
    # compute strain in each sub-domain
        
    sub_domain_F_all, sub_domain_row_all, sub_domain_col_all, num_pts_sd = compute_sub_domain_position_strain_all(tracker_row_all, tracker_col_all, sub_domain_box_list)
    # save the sub-domain strains
    strain_sub_domain_info = np.asarray([[num_tile_row, num_tile_col], [tile_dim_pix, tile_dim_pix]])
    
    saved_paths = save_sub_domain_strain(folder_path, sub_domain_F_all, sub_domain_row_all, sub_domain_col_all, strain_sub_domain_info, fname=save_fname)
    return saved_paths


def visualize_sub_domain_strain(
    folder_path: Path,
    col_min: Union[float, int] = -0.025,
    col_max: Union[float, int] = 0.025,
    col_map: object = plt.cm.RdBu,
    fname: str = "strain",
    is_rotated: bool = True,
) -> List:
    """Given a folder path where strain analysis has already been run. Will save visualizations."""
    # read image files
    movie_folder_path = folder_path.joinpath("movie").resolve()
    name_list_path = ia.image_folder_to_path_list(movie_folder_path)
    tiff_list = ia.read_all_tiff(name_list_path)
    # load the strain results
    sub_domain_F_rr_all, sub_domain_F_rc_all, sub_domain_F_cr_all, sub_domain_F_cc_all, sub_domain_row_all, sub_domain_col_all, strain_info = load_sub_domain_strain(folder_path)
    
    num_sd_row = int(strain_info[0, 0])
    num_sd_col = int(strain_info[0, 1])
    sub_domain_side = strain_info[1, 0]

    # convert the strain results to Ecc for plotting
    sub_domain_Ecc_all, sub_domain_Ecr_all, sub_domain_Err_all = F_to_E_all(sub_domain_F_rr_all, sub_domain_F_rc_all, sub_domain_F_cr_all, sub_domain_F_cc_all)

    # create png visualizations
    png_path_list_Ecc = pngs_sub_domain_strain(folder_path, tiff_list, sub_domain_row_all, sub_domain_col_all, sub_domain_Ecc_all, sub_domain_side, "Ecc", col_min, col_max, col_map, fname, save_eps=False)
    png_path_list_Ecr = pngs_sub_domain_strain(folder_path, tiff_list, sub_domain_row_all, sub_domain_col_all, sub_domain_Ecr_all, sub_domain_side, "Ecr", col_min, col_max, col_map, fname, save_eps=False)
    png_path_list_Err = pngs_sub_domain_strain(folder_path, tiff_list, sub_domain_row_all, sub_domain_col_all, sub_domain_Err_all, sub_domain_side, "Err", col_min, col_max, col_map, fname, save_eps=False)    
    # create gif
    gif_path_Ecc = create_gif(folder_path, png_path_list_Ecc, fname="sub_domain_strain_Exx")
    gif_path_Ecr = create_gif(folder_path, png_path_list_Ecr, fname="sub_domain_strain_Exy")
    gif_path_Err = create_gif(folder_path, png_path_list_Err, fname="sub_domain_strain_Eyy")
    # create subdomain locations legend
    loc_legend_path = png_sub_domains_numbered(folder_path, tiff_list[0], sub_domain_row_all[:,0], sub_domain_col_all[:,0], sub_domain_side, num_sd_row, num_sd_col)
    # create subdomain strain timeseries plots
    timeseries_path_list_Ecc = png_sub_domain_strain_timeseries_all(folder_path, sub_domain_Ecc_all, num_sd_row, num_sd_col, output="Ecc", fname="strain_timeseries_Exx", ylabel="strain Exx")
    timeseries_path_list_Ecr = png_sub_domain_strain_timeseries_all(folder_path, sub_domain_Ecr_all, num_sd_row, num_sd_col, output="Ecr", fname="strain_timeseries_Exy", ylabel="strain Exy")
    timeseries_path_list_Err = png_sub_domain_strain_timeseries_all(folder_path, sub_domain_Err_all, num_sd_row, num_sd_col, output="Err", fname="strain_timeseries_Eyy", ylabel="strain Eyy")
    return png_path_list_Ecc, png_path_list_Ecr, png_path_list_Err, gif_path_Ecc, gif_path_Ecr, gif_path_Err, loc_legend_path, timeseries_path_list_Ecc, timeseries_path_list_Ecr, timeseries_path_list_Err

 


