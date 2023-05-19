import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
from skimage import io
import imageio
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("input_folder", help="the user input folder location")
args = parser.parse_args()
input_folder_str = args.input_folder

# specify file paths
self_path_file = Path(__file__)
self_path = self_path_file.resolve().parent
input_folder = self_path.joinpath(input_folder_str).resolve()
data_path = input_folder.joinpath("original_files").resolve()
new_data_path = input_folder.joinpath("prepared_files").resolve()
if new_data_path.exists() is False:
    os.mkdir(new_data_path)

pngs_data_path = input_folder.joinpath("visualize_files").resolve()
if pngs_data_path.exists() is False:
    os.mkdir(pngs_data_path)

# specify details of what to segment
# file_names = ["control1", "control2", "control3", "control4", "stim1", "stim2", "stim3", "stim4"]
file_names = ["control1"]

first_frame = [8]
last_frame = [98]
light_on_1 = [9]
light_on_2 = [39]
light_on_3 = [69]

# first_frame = [8, 146, 27, 197, 31, 356, 228, 267]
# last_frame = [98, 236, 117, 287, 121, 446, 318, 357]
# light_on_1 = [9, 147, 28, 198, 32, 357, 229, 268]
# light_on_2 = [39, 177, 58, 228, 62, 387, 258, 298]
# light_on_3 = [69, 207, 88, 258, 91, 417, 288, 327]



# adjust all given values for 0 based indexing
def adjust_0_based(li):
    li_new = []
    for kk in range(0, len(li)):
        li_new.append(li[kk] - 1)
    return li_new


first_frame = adjust_0_based(first_frame)
last_frame = adjust_0_based(last_frame)
light_on_1 = adjust_0_based(light_on_1)
light_on_2 = adjust_0_based(light_on_2)
light_on_3 = adjust_0_based(light_on_3)


# import and trim all movies and convert to grayscale (file sizes are quite large otherwise)
def rgb_to_grey(data):
    img_grey = 1.0 / 3.0 * data[:, :, :, 0] + 1.0 / 3.0 * data[:, :, :, 1] + 1.0 / 3.0 * data[:, :, :, 2]
    return img_grey


all_data = []
for kk in range(0, len(file_names)):
    data_rgb = io.imread(str(data_path) + "/" + file_names[kk] + ".tif")
    data_rgb = data_rgb[first_frame[kk]: last_frame[kk] + 1, :, :, :]
    data = rgb_to_grey(data_rgb)
    all_data.append(data)


# detect lit up frames
all_mean_intensity = []
median_intensity = []
for kk in range(0, len(all_data)):
    mean_intensity = []
    for jj in range(all_data[kk].shape[0]):
        mean_intensity.append(np.mean(all_data[kk][jj, :, :]))
    all_mean_intensity.append(mean_intensity)
    median_intensity.append(np.median(mean_intensity))

thresh = 2  # frames with brightness above this threshold
all_keep = []
for kk in range(0, len(all_data)):
    cutoff = median_intensity[kk] + thresh
    plt.plot(all_mean_intensity[kk], "k")
    keep_list = []
    for jj in range(0, len(all_mean_intensity[kk])):
        if all_mean_intensity[kk][jj] > cutoff:
            plt.plot(jj, all_mean_intensity[kk][jj], "ro")
            keep_list.append(0)
        else:
            keep_list.append(1)
    all_keep.append(keep_list)
plt.xlabel('Frames')
plt.ylabel('pixel value')
plt.savefig(str(pngs_data_path) + "/threshold.pdf", format='pdf')


# remove lit up frames
all_data_bright_removed = []
for kk in range(0, len(all_data)):
    bright_removed = all_data[kk][np.asarray(all_keep[kk]) > 0, :, :]
    all_data_bright_removed.append(bright_removed)


# save the final data in our modified fashion
for kk in range(0, len(file_names)):
    np.save(str(new_data_path) + "/" + file_names[kk] + ".npy", all_data_bright_removed[kk])

# save data in our modified fashion as individual TIF frames
for ff in file_names:
    frames_npy = np.load(str(new_data_path) + '/' + ff + '.npy')
    num_frames = len(frames_npy)
    data_res_folder = input_folder.joinpath(ff).resolve()
    if data_res_folder.exists() is False:
        os.mkdir(data_res_folder)
    imgs_folder = data_res_folder.joinpath('movie').resolve()
    if imgs_folder.exists() is False:
        os.mkdir(imgs_folder)
    for ii in range(num_frames):
        gray_scale_im = frames_npy[ii].astype('uint8')
        imageio.imwrite(str(imgs_folder) + '/%04d.TIF'%(ii), gray_scale_im)