import yt
import numpy as np
import matplotlib.pyplot as plt
from yt.units import mh
from scipy.signal import windows
from scipy.ndimage import zoom
from scipy.signal.windows import hann




# spacing = np.linspace(0, 500, 100, dtype=int)

def downsample_to(target_shape, data):
    factor = np.array(target_shape) / np.array(data.shape)
    return zoom(data, factor, order=1)  # linear interpolation

def l2_error(ref, test):
    return np.sqrt(np.sum((ref - test)**2) / np.sum(ref**2))

def compute_power_spectrum_2d(data, apply_window=True):
    """Compute 2D power spectrum and return 1D radial average."""
    if hasattr(data, "value"):
        data = data.value  # strip unyt if needed

    # Remove mean
    data = data - np.mean(data)

    # Apply Hann window to suppress FFT edge ringing
    if apply_window:
        ny, nx = data.shape
        window = hann(ny)[:, None] * hann(nx)
        data *= window

    # 2D FFT and power spectrum
    fft2d = np.fft.fftshift(np.fft.fft2(data))
    power2d = np.abs(fft2d)**2

    # Get 1D radial profile
    y, x = np.indices(data.shape)
    center = np.array([(x.max() - x.min()) / 2.0, (y.max() - y.min()) / 2.0])
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(int)

    # Radial binning
    tbin = np.bincount(r.ravel(), power2d.ravel())
    nr = np.bincount(r.ravel())
    radial_profile = tbin / np.maximum(nr, 1)

    return radial_profile

spacing = [101]

time_step = 0

for space in spacing:
    # For wavelength 2L - I need 2-pc/1.58-km/s to get desired properties - 15uG, 1.58 km/s and 200 cm-3    ?
    # For wavelength L - I need 1-pc/1.0-km/s to get desired properties - 15uG, 1.58 km/s and 200 cm-3      ?
    # For wavelength L - I need 1-pc/1.0-km/s to get desired properties - 15uG, 1.58 km/s and 200 cm-3      ?
    unit_base={"length_unit": (1.0,"pc"), "time_unit": (1.0,"1.0 * pc / (1.0 * km/s)"), "mass_unit": (2.225e34,"g")}
    ds_512 = yt.load(f'./density-analysis/512/LinWave.out2.{space:05d}.athdf', units_override=unit_base)

    ds_128 = yt.load(f'./density-analysis/128/LinWave.out2.{space:05d}.athdf', units_override=unit_base)

    ds_256 = yt.load(f'./density-analysis/256/LinWave.out2.{space:05d}.athdf', units_override=unit_base)

    size = 512
    grid_x_512 = np.linspace(0, 1, 512)
    grid_x_256 = np.linspace(0, 1, 256)
    grid_x_128 = np.linspace(0, 1, 128)

    left_edge = ds_512.domain_left_edge
    dims =  ds_512.domain_dimensions

    data_512 = ds_512.covering_grid(level=0, left_edge=left_edge, dims=dims)

    left_edge = ds_256.domain_left_edge
    dims =  ds_256.domain_dimensions

    data_256 = ds_256.covering_grid(level=0, left_edge=left_edge, dims=dims)

    left_edge = ds_128.domain_left_edge
    dims =  ds_128.domain_dimensions

    data_128 = ds_128.covering_grid(level=0, left_edge=left_edge, dims=dims)

    accurate_number_density_512 = data_512['rho'].to("g/cm**3") / (2.34 * mh)
    number_density_512 = accurate_number_density_512.to("cm**-3")[:, :, 0].T

    accurate_number_density_256 = data_256['rho'].to("g/cm**3") / (2.34 * mh)
    number_density_256 = accurate_number_density_256.to("cm**-3")[:, :, 0].T

    accurate_number_density_128 = data_128['rho'].to("g/cm**3") / (2.34 * mh)
    number_density_128 = accurate_number_density_128.to("cm**-3")[:, :, 0].T

    number_density_512_to_256 = zoom(number_density_512, 128 / 512, order=3)
    number_density_256_to_128 = zoom(number_density_256, 128 / 256, order=3)

    # column_density = ds_128.integrate(("gas", "density"), axis=("index", "z")).to("g/cm**2")

    error_256 = l2_error(number_density_128, number_density_512_to_256 * number_density_256.units)

    error_128 = l2_error(number_density_128, number_density_256_to_128 * number_density_128.units)

    h128 = 1/128
    h256 = 1/256
    h512 = 1/512

    p = np.log(error_128 / error_256) / np.log(h256 / h512)
    print(f"Estimated convergence order: {p:.2f}")

    print(error_256, error_128)


    # ps128 = compute_power_spectrum_2d(number_density_128)
    # ps256 = compute_power_spectrum_2d(number_density_256)
    # ps512 = compute_power_spectrum_2d(number_density_512)

    # # k-axis: just index of radial bins, assuming uniform grid spacing
    # k128 = np.arange(len(ps128))
    # k256 = np.arange(len(ps256))
    # k512 = np.arange(len(ps512))

    # plt.figure(figsize=(8, 6))
    # plt.plot(k128, ps128, label="128²")
    # plt.plot(k256, ps256, label="256²")
    # plt.plot(k512, ps512, label="512²")
    # plt.xlabel("k (wavenumber)")
    # plt.ylabel("Power")
    # plt.title("1D Power Spectrum of Number Density")
    # plt.legend()
    # plt.grid(True, which="both", ls="--", alpha=0.5)
    # plt.tight_layout()
    # plt.show()

    mean_density_512 = np.mean(number_density_512)
    mean_density_256 = np.mean(number_density_256)
    mean_density_128 = np.mean(number_density_128)

    print(mean_density_512, mean_density_256, mean_density_128)


    

    # avg_density_512 = [0]*512
    # avg_density_256 = [0]*256
    # avg_density_128 = [0]*128

    # size = [128, 256, 512]

    # for s in size:
    #     for x in range(s):
    #         if s == 128:
    #             # avg_density_128[x] = (np.max(number_density_128[x]) - np.min(number_density_128[x])) * (100 / mean_density_128)
    #             avg_density_128[x] = np.mean(number_density_128[x])
    #         elif s == 256:
    #             # avg_density_256[x] = (np.max(number_density_256[x]) - np.min(number_density_256[x])) * (100 / mean_density_256)
    #             avg_density_256[x] = np.mean(number_density_256[x])
    #         else:
    #             # avg_density_512[x] = (np.max(number_density_512[x]) - np.min(number_density_512[x])) * (100 / mean_density_512)
    #             avg_density_512[x] = np.mean(number_density_512[x])
    
    # plt.plot(grid_x_128, avg_density_128)
    # plt.plot(grid_x_256, avg_density_256)
    # plt.plot(grid_x_512, avg_density_512)
    # plt.xlabel("X (pc)")
    # plt.ylabel("Number Density")
    # plt.title("Density analysis")
    # plt.grid()
    # plt.show()

    


    # averaged_data = number_density - np.mean(number_density)
    # window = windows.hann(averaged_data.shape[0])[:, None] * windows.hann(averaged_data.shape[1])
    # data_window = averaged_data * window


    # # print(number_density.shape)
    # # print("Accurate Number Density", accurate_number_density)

    # fft_image = np.fft.fft2(data_window)
    # fft_image = np.fft.fftshift(fft_image)
    # power_spectrum = np.abs(fft_image)**2


    # print("TOTAL Mass", grid[('gas', 'cell_mass')].to("g").sum())

#     print("Gamma", ds.gamma)
#     ds.gamma = 1.0

#     print(ds.derived_field_list)

#     pressure = grid[('gas', 'pressure')].to('dyne/cm**2')
#     print("Pressure", pressure)

#     iso_sound_speed = (pressure / density)**0.5
#     print("Sound speed", iso_sound_speed.to("km/s"))

#     time_step += 1

#     if time_step == 1:
#         break

#     if time_evolved > 4.8:
#         print("How Many Steps done", time_step)
#         break


    # plt.figure(figsize=(6, 5))
    # print(X.dtype, Y.dtype, number_density.dtype)
    # plt.pcolormesh(X, Y, number_density.T, cmap='hot', shading='auto')
    # plt.title(f'Density Map {space}')

    # plt.imshow(number_density.T, cmap='hot', 
    #         extent=[X.min(), X.max(), Y.min(), Y.max()],
    #         origin='lower', aspect='auto')
    # plt.colorbar(label="Density")
    # plt.savefig(f'density_map_{time_evolved:.4f}.png', dpi=300, bbox_inches='tight')
    # plt.xlabel('X axis')
    # plt.ylabel('Y axis')

    # plt.loglog(power_spectrum)
    

    # slc = yt.SlicePlot(ds, "z", ("p


#     [('athena_pp', 'Bcc1'), ('athena_pp', 'Bcc2'), ('athena_pp', 'Bcc3'), ('athena_pp', 'cell_volume'), ('athena_pp', 'dx'), ('athena_pp', 'dy'), ('athena_pp', 'dz'), ('athena_pp', 'path_element_x'), ('athena_pp', 'path_element_y'), ('athena_pp', 'path_element_z'), ('athena_pp', 'rho'), ('athena_pp', 'vel1'), ('athena_pp', 'vel2'), ('athena_pp', 'vel3'), ('athena_pp', 'volume'), ('athena_pp', 'x'), ('athena_pp', 'y'), ('athena_pp', 'z'), 
#      ('gas', 'alfven_speed'), ('gas', 'angular_momentum_magnitude'), ('gas', 'angular_momentum_x'), ('gas', 'angular_momentum_y'), ('gas', 'angular_momentum_z'), ('gas', 'averaged_density'), ('gas', 'cell_mass'), ('gas', 'cell_volume'), ('gas', 'density'), ('gas', 'density_gradient_magnitude'), ('gas', 'density_gradient_x'), ('gas', 'density_gradient_y'), ('gas', 'density_gradient_z'), ('gas', 'dx'), ('gas', 'dy'), 
#      ('gas', 'dynamical_time'), ('gas', 'dz'), ('gas', 'four_velocity_magnitude'), ('gas', 'four_velocity_t'), ('gas', 'four_velocity_x'), ('gas', 'four_velocity_y'), ('gas', 'four_velocity_z'), ('gas', 'kinetic_energy_density'), ('gas', 'lorentz_factor'), ('gas', 'mach_alfven'), ('gas', 'magnetic_energy_density'), ('gas', 'magnetic_field_los'), ('gas', 'magnetic_field_magnitude'), ('gas', 'magnetic_field_strength'), ('gas', 'magnetic_field_x'), ('gas', 'magnetic_field_y'), ('gas', 'magnetic_field_z'), ('gas', 'magnetic_pressure'), 
#      ('gas', 'mass'), ('gas', 'mean_molecular_weight'), ('gas', 'momentum_density_x'), ('gas', 'momentum_density_y'), ('gas', 'momentum_density_z'), ('gas', 'momentum_x'), ('gas', 'momentum_y'), ('gas', 'momentum_z'), ('gas', 'number_density'), ('gas', 'path_element_x'), ('gas', 'path_element_y'), ('gas', 'path_element_z'), ('gas', 'relative_magnetic_field_x'), ('gas', 'relative_magnetic_field_y'), ('gas', 'relative_magnetic_field_z'), ('gas', 'relative_velocity_x'), ('gas', 'relative_velocity_y'), ('gas', 'relative_velocity_z'), ('gas', 'shear'), ('gas', 'specific_angular_momentum_magnitude'), 
#      ('gas', 'specific_angular_momentum_x'), ('gas', 'specific_angular_momentum_y'), ('gas', 'specific_angular_momentum_z'), ('gas', 'velocity_los'), ('gas', 'velocity_magnitude'), ('gas', 'velocity_x'), ('gas', 'velocity_y'), ('gas', 'velocity_z'), ('gas', 'volume'), ('gas', 'vorticity_magnitude'), ('gas', 'vorticity_squared'), ('gas', 'vorticity_x'), ('gas', 'vorticity_y'), ('gas', 'vorticity_z'), ('gas', 'x'), ('gas', 'y'), ('gas', 'z'), 
#      ('index', 'cell_volume'), ('index', 'cylindrical_radius'), ('index', 'cylindrical_theta'), ('index', 'cylindrical_z'), ('index', 'dx'), ('index', 'dy'), ('index', 'dz'), ('index', 'grid_indices'), ('index', 'grid_level'), ('index', 'ones'), ('index', 'ones_over_dx'), ('index', 'path_element_x'), ('index', 'path_element_y'), ('index', 'path_element_z'), ('index', 'radius'), ('index', 'spherical_phi'), ('index', 'spherical_radius'), ('index', 'spherical_theta'), ('index', 'virial_radius_fraction'), ('index', 'volume'), ('index', 'x'), ('index', 'y'), ('index', 'z'), ('index', 'zeros')]