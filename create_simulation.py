import imageio
import glob
from PIL import Image
import numpy as np
import re

# Load images (sorted correctly)
# image_files = sorted(glob.glob('density_map_*.png'))  # Adjust pattern
# images = [imageio.imread(f) for f in image_files]

def extract_number(filename):
    match = re.search(r'density_map_([\d.]+) Myr\.png', filename)
    return float(match.group(1)) if match else float('inf')  # fallback for non-matching files

target_shape = (500, 500)  # Set your desired uniform size

# Resize all images
images = []
for file in sorted(glob.glob('density_map_*.png'), key=extract_number):
    img = Image.open(file).resize(target_shape)
    images.append(np.array(img))

# Save as GIF
imageio.mimsave('simulation.gif', images, fps=10)
