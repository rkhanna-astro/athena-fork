import imageio
import glob
from PIL import Image
import numpy as np

# Load images (sorted correctly)
image_files = sorted(glob.glob('density_map_*.png'))  # Adjust pattern
images = [imageio.imread(f) for f in image_files]

target_shape = (500, 500)  # Set your desired uniform size

# Resize all images
images = []
for file in sorted(glob.glob('density_map_*.png')):
    img = Image.open(file).resize(target_shape)
    images.append(np.array(img))

# Save as GIF
imageio.mimsave('simulation.gif', images, fps=10)