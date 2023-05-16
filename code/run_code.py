import argparse
import image_analysis as ia
import strain_analysis as sa
from pathlib import Path
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument("input_folder", help="the user input folder location")
args = parser.parse_args()
input_folder_str = args.input_folder


self_path_file = Path(__file__)
self_path = self_path_file.resolve().parent
input_folder = self_path.joinpath(input_folder_str).resolve()

# run the tracking
ia.run_tracking(input_folder)

# run the tracking visualization
automatic_color_constraint = False # Put False if manual limits are to be specified
col_min_abs = 0
col_max_abs = 12
col_min_row = -4
col_max_row = 10
col_min_col = -4
col_max_col = 10
col_map = plt.cm.viridis
ia.run_visualization(input_folder, col_min_abs, col_max_abs, col_min_row, col_max_row, col_min_col, col_max_col, col_map = col_map)


# run the strain analysis (will automatically rotate based on the mask)
# pillar_clip_fraction = 0.2 #(used for control4)
pillar_clip_fraction = 0.1
clip_columns = True
clip_rows = True
shrink_row = 0.1
shrink_col = 0.1
tile_dim_pix = 100
# tile_dim_pix = 80 #(used for control4)
num_tile_row = 10
num_tile_col = 10
tile_style = 1 # or 2
manual_sub = False # or True
sub_extents = None # if manual_sub = True provide as [r0,r1,c0,c1]
sa.run_sub_domain_strain_analysis(input_folder, pillar_clip_fraction, shrink_row, shrink_col, tile_dim_pix, num_tile_row, num_tile_col, tile_style, is_rotated = False,clip_columns=clip_columns,clip_rows=clip_rows, manual_sub=manual_sub, sub_extents=sub_extents)

# visualize the strain analysis results
col_min = -0.02
col_max = 0.02
col_map = plt.cm.RdBu
sa.visualize_sub_domain_strain(input_folder, automatic_color_constraint, col_min, col_max, col_map, is_rotated = False)    


