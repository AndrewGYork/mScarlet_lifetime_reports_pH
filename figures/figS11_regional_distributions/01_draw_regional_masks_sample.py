from pathlib import Path
import numpy as np
import tifffile as tf
import pandas as pd
from skimage.measure import label, regionprops
from skimage.color import label2rgb
from skimage.morphology import binary_erosion, disk, remove_small_objects

# This script allows you to save the "region mask" to delineate cellular
# zones as performed in Grinstein JCB
# (https://doi.org/10.1083/jcb.201507112) to look at how pHlys varies by
# region. There is another script that generates and applies these masks
# in source_data; this is merely for generating a sample figure.

# I use "erode" to go a set distance in from the outside of the cell. By
# doing this sequentially, I identify shells/zones at the outside of the
# cell. I stayed pretty close to the grinstein analysis here. Note that
# this is imperfect where the cell is touching the edge, as I don't know
# if it is near another edge of the cell but just out of view.

# This script only processes a single sample image and outputs the
# region mask for illustrative purposes/overlaying in ImageJ.

# generate some paths
current_dir = Path.cwd()
manuscript_dir = current_dir.parents[1]
source_dir = manuscript_dir / 'source_data' / 'bafilomycin_endpoint'
cell_roi_path = source_dir / 'hand_cell_seg' / '2021-11-17_13-3_DMSO_mask.tif'

# load in the hand-drawn cell ROIs 
cell_mask = tf.imread(cell_roi_path).astype('bool')
# let's identify "shells" coming inward from the edge of the cell
eroded1 = binary_erosion(cell_mask, disk(14)) # 5 um from edge
eroded2 = binary_erosion(eroded1, disk(14))
eroded3 = binary_erosion(eroded2, disk(14))
# now let's mark out our regions
region_mask = np.zeros(cell_mask.shape, dtype='uint8')
region_mask[cell_mask ^ eroded1] = 1
region_mask[eroded1 ^ eroded2] = 2
region_mask[eroded2 ^ eroded3] = 3
region_mask[eroded3] = 4

tf.imwrite('2021-11-17_13-3_region_mask.tif', region_mask.astype('float32'),
           imagej=True)
