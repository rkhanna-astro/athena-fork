import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
import yt
from yt.units import mh

# def compute_mean_contrast(number_density_2d, cut_axis='x', smoothing_sigma=2.0):
contrasts = []

spacing = [150]

time_step = 0
grid = 512
grid_x = np.linspace(0, 1, grid)
grid_y = np.linspace(0, 1, grid)

X,Y = np.meshgrid(grid_x, grid_y)

for space in spacing:
    # For wavelength 2L - I need 2-pc/1.58-km/s to get desired properties - 15uG, 1.58 km/s and 200 cm-3    ?
    # For wavelength L - I need 1-pc/1.0-km/s to get desired properties - 15uG, 1.58 km/s and 200 cm-3      ?
    # For wavelength L - I need 1-pc/1.0-km/s to get desired properties - 15uG, 1.58 km/s and 200 cm-3      ?
    unit_base={"length_unit": (1.0,"pc"), "time_unit": (1.0,"1.0 * pc / (1.0 * km/s)"), "mass_unit": (2.225e34,"g")}

    ds_256 = yt.load(f'./time-varying-field/512/LinWave.out2.{space:05d}.athdf', units_override=unit_base)

    size = 512
    grid_x = np.linspace(0, 1, size)
    grid_y = np.linspace(0, 1, size)

    left_edge = ds_256.domain_left_edge
    right_edge = ds_256.domain_right_edge
    dims =  ds_256.domain_dimensions

    # print(left_edge, dims)

    data_256 = ds_256.covering_grid(level=0, left_edge=left_edge, dims=dims)

    accurate_number_density_256 = data_256['rho'].to("g/cm**3") / (2.34 * mh)
    number_density = accurate_number_density_256.to("cm**-3")[:, :, 0].T

    data = number_density

    extrema_coords_max = []
    extrema_coords_min = []

    for y in range(127, 128):
        raw_profile = data[y, :]
        
        smoothed_profile = gaussian_filter1d(raw_profile, sigma=3.0)

        # Find all extrema (both peaks and troughs)
        peaks, _ = find_peaks(smoothed_profile, distance=10)
        print(raw_profile[peaks])
        troughs, _ = find_peaks(-smoothed_profile, distance=10)
        extrema = np.sort(np.concatenate([peaks, troughs]))

        extrema_coords_max.extend([(grid_y[y], grid_x[i]) for i in peaks])
        extrema_coords_min.extend([(grid_y[y], grid_x[i]) for i in troughs])

        # Compute contrast between each successive extrema pair
        for j in range(len(extrema) - 1):
            a = extrema[j]
            b = extrema[j + 1]
            I1 = raw_profile[a]
            I2 = raw_profile[b]

            # Skip zero or negative intensity (can happen due to noise)
            if I1 + I2 == 0:
                continue

            contrast = 2 * abs(I1 - I2) / (I1 + I2)
            contrasts.append(contrast)

        plt.figure(figsize=(6, 5))
        plt.pcolormesh(X, Y, data, cmap='hot', shading='auto', vmin = 100, vmax = 250)
        ys, xs = zip(*extrema_coords_max)
        plt.scatter(xs, ys, color='cyan', marker='^', s=10, label='Crests')

        ys, xs = zip(*extrema_coords_min)
        plt.scatter(xs, ys, color='lime', marker='v', s=10, label='Troughs')
        plt.title(f'Number Density Map (time={1:.4f})')


        plt.colorbar(label=r'Number Density (cm$^{-3}$)')
        plt.xlabel('x (1 pc)')
        plt.ylabel('y (1 pc)')
        plt.savefig(f'overlap_map_{1:.4f}.png', dpi=300, bbox_inches='tight')

    mean_contrast = np.max(contrasts) if contrasts else 0.0
    print(mean_contrast)


# accurate_number_density_256 = data_256['rho'].to("g/cm**3") / (2.34 * mh)
#     number_density = accurate_number_density_256.to("cm**-3")[:, :, 0].T

#     velocity_x = data_256[('gas', 'velocity_x')].to("km/s")[:, : , 0].T

#     magnetic_field_y = data_256[('gas', 'magnetic_field_y')].to("uG")[:, : , 0].T

#     averaged_number_density = np.zeros((16, size))
#     velocity_plot = np.zeros((16, size))
#     magnetic_plot = np.zeros((16, size))

#     for x in range(16, size, 16):
#             averaged_value = (number_density[x-1] + number_density[x] + number_density[x+1]) / 3.0
#             averaged_number_density[int(x / 16)] = averaged_value

#             averaged_value = (velocity_x[x-1] + velocity_x[x] + velocity_x[x+1]) / 3.0
#             velocity_plot[int(x / 16)] = averaged_value

#             magnetic_plot[int(x / 16)] = (magnetic_field_y[x-1] + magnetic_field_y[x] + magnetic_field_y[x+1]) / 3.0