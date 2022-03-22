## Simple script to use skimage to get centroid position of nuclear masks
import imageio
from skimage.measure import label, regionprops, regionprops_table
from skimage import data, filters, measure, morphology
import sys
import pandas as pd

mask_file = sys.argv[1]
outfile = sys.argv[2]

mask_im = imageio.imread(mask_file)
labels = measure.label(mask_im)
props = regionprops_table(labels, properties=('centroid',
                                                 'orientation',
                                                 'axis_major_length',
                                                 'axis_minor_length'))
props_df = pd.DataFrame(props)
props_df.to_csv(outfile)