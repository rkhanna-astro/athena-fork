import yt
import numpy as np
import matplotlib.pyplot as plt
import imageio
from PIL import Image
import glob

# Load images (sorted correctly)
# image_files = sorted(glob.glob('density_map_*.png'))  # Adjust pattern
# images = [imageio.imread(f) for f in image_files]

target_shape = (500, 500)  # Set your desired uniform size

# Resize all images
images = []
for file in sorted(glob.glob('density_map_*.png')):
    img = Image.open(file).resize(target_shape)
    images.append(np.array(img))

# Save as GIF
imageio.mimsave('simulation.gif', images, fps=10)

# spacing = np.linspace(0, 500, 100, dtype=int)

# for space in spacing:
#     ds = yt.load(f'LinWave.out2.{space:05d}.athdf')
#     print(ds.field_list)

#     grid = 256
#     grid_x = np.linspace(0, 1, grid)
#     grid_y = np.linspace(0, 1, grid)

#     X,Y = np.meshgrid(grid_x, grid_y)

#     left_edge = ds.domain_left_edge
#     dims =  ds.domain_dimensions

#     print(ds.domain_left_edge)
#     print(ds.domain_dimensions)

#     grid = ds.covering_grid(level=0, left_edge=left_edge, dims=dims)

#     density = grid['rho'][:, :, 0]
#     print(density.shape)
#     print(f"X dtype: {X.dtype}, Y dtype: {Y.dtype}, density dtype: {density.dtype}")

#     plt.figure(figsize=(6, 5))
#     print(X.dtype, Y.dtype, density.dtype)
#     plt.pcolormesh(X, Y, density.T, cmap='hot', shading='auto')
#     plt.title(f'Density Map {space}')

#     plt.imshow(density.T, cmap='hot', 
#             extent=[X.min(), X.max(), Y.min(), Y.max()],
#             origin='lower', aspect='auto')
#     plt.colorbar(label="Density")
#     plt.savefig(f'density_map_{space:03d}.png', dpi=300, bbox_inches='tight')
#     plt.xlabel('X axis')
#     plt.ylabel('Y axis')


    # slc = yt.SlicePlot(ds, "z", ("parthenon", "prim_velocity_1")).save()
    # slc = yt.ProjectionPlot(ds, "x", ("parthenon", "prim_velocity_1")).save()
    # plot_1 = yt.LinePlot(
    #     ds, [("parthenon", "prim_velocity_1")], (0.0, 0.0, 0.0), (0.0, 3.0, 0.0), 256
    # )

    # plot_1.annotate_legend(("parthenon", "prim_velocity_1"))

    # # Save the line plot
    # plot_1.save()

    # plot_2 = yt.LinePlot(
    #     ds, [("parthenon", "prim_velocity_2")], (0.0, 0.0, 0.0), (0.0, 1.5, 0.0), 256
    # )

    # plot_2.annotate_legend(("parthenon", "prim_velocity_2"))

    # # Save the line plot
    # plot_2.save()

    # plot_4 = yt.LinePlot(
    #     ds, [("parthenon", "prim_velocity_3")], (0.0, 0.0, 0.0), (0.0, 1.5, 0.0), 256
    # )

    # plot_4.annotate_legend(("parthenon", "prim_velocity_3"))

    # # Save the line plot
    # plot_4.save()

    # plot_3 = yt.LinePlot(
    #     ds, [("parthenon", "prim_density")], (0.0, 0.0, 0.0), (1.0, 1.0, 0.0), 32
    # )

    # plot_3.annotate_legend(("parthenon", "prim_density"))

    # # Save the line plot
    # plot_3.save()

    # yt.SlicePlot(
    #     ds,
    #     "z",
    #     [("parthenon", "prim_density")]
    # ).save()

    # im, sc = yt.volume_render(ds, field=("parthenon", "prim_velocity_1"))