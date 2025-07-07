import yt
import numpy as np
import matplotlib.pyplot as plt
from yt.units import mh
from scipy.signal import windows
from scipy.ndimage import zoom
from scipy.signal.windows import hann
from scipy.ndimage import gaussian_filter1d
from scipy.fft import fft, fftshift

from astropy.timeseries import LombScargle
from scipy.signal import welch

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

    # # Apply Hann window to suppress FFT edge ringing
    # if apply_window:
    #     ny, nx = data.shape
    #     window = hann(ny)[:, None] * hann(nx)
    #     data *= window

    # 2D FFT and power spectrum
    fft2d = np.fft.fft2(data)
    power_2d = np.abs(fft2d)**2
    power_2d = np.fft.fftshift(power_2d)


    # Get 1D radial profile
    y, x = np.indices(data.shape)
    center = np.array([(x.max() - x.min()) / 2.0, (y.max() - y.min()) / 2.0])
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(int)

    # Radial binning
    tbin = np.bincount(r.ravel(), power_2d.ravel())
    nr = np.bincount(r.ravel())
    radial_profile = tbin / np.maximum(nr, 1)

    return radial_profile

spacing = [100]

time_step = 0

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

    velocity_x = data_256[('gas', 'velocity_x')].to("km/s")[:, : , 0].T

    magnetic_field_y = data_256[('gas', 'magnetic_field_y')].to("uG")[:, : , 0].T

    averaged_number_density = np.zeros((int(size / 16), size))
    velocity_plot = np.zeros((int(size / 16), size))
    magnetic_plot = np.zeros((int(size / 16), size))

    for x in range(16, size, 16):
            averaged_value = (number_density[x-1] + number_density[x] + number_density[x+1]) / 3.0
            averaged_number_density[int(x / 16)] = averaged_value

            averaged_value = (velocity_x[x-1] + velocity_x[x] + velocity_x[x+1]) / 3.0
            velocity_plot[int(x / 16)] = averaged_value

            magnetic_plot[int(x / 16)] = (magnetic_field_y[x-1] + magnetic_field_y[x] + magnetic_field_y[x+1]) / 3.0

    # print(averaged_number_density)



# min_freq = 0.1
# max_freq = 25

# freqs = np.linspace(min_freq, max_freq, 256)

# ls = LombScargle(freqs, averaged_number_density[1])
# density_power = ls.power(freqs)


for ind in range(16, size, 16):
    i = int(ind / 16)

    signal = averaged_number_density[i]
    vel_signal = velocity_plot[i]
    mag_signal = magnetic_plot[i]
    

    f1, P1 = welch(signal, fs=256, scaling='density')
    f2, P2 = welch(vel_signal, fs=256, scaling='density')
    f3, P3 = welch(mag_signal, fs=256, scaling='density')

    P_rho_norm = P1 / np.trapezoid(P1, f1)
    P_vel_norm = P2 / np.trapezoid(P2, f2)
    P_mag_norm = P3 / np.trapezoid(P3, f3)

    # first_10_indices = np.arange(1, 25)   # k = 1 to 10

    plt.figure(figsize=(8, 8))
    plt.plot(f1[:25], P_rho_norm[:25], label="density")
    plt.plot(f2[:25], P_vel_norm[:25], label="velocity")
    plt.plot(f3[:25], P_mag_norm[:25], label="magnetic")

    plt.legend()
    plt.xlabel("k (wavenumber)")
    plt.ylabel("Power")
    plt.title("1D Power Spectrum of Density, Velocity, Magnetic Field")
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()
    # plt.show()
    plt.savefig(f'power_spectra_comparison_{i}.png', dpi=300, bbox_inches='tight')


# print(freqs)

# velocity_signal = velocity_plot[1]
# fft_vals = fft(velocity_signal)
# fourier_shifted = fftshift(fft_vals)
# power_spectra = np.abs(fourier_shifted)**2 / 256

# vel_power = power_spectra[positive_freqs]
# print(velocity_signal)
# vel_fft_vals = np.fft.fft(velocity_signal)
# velocity_power_spectra = np.abs(vel_fft_vals[:256//2])**2

# magnetic_signal = magnetic_plot[1]
# mag_fft_vals = np.fft.fft(magnetic_signal)
# magnetic_power_spectra = np.abs(mag_fft_vals[:256//2])**2

# freqs = freqs[first_10_indices]

# den_power = power_spectra[first_10_indices]
# vel_power = velocity_power_spectra[first_10_indices]
# mag_power = magnetic_power_spectra[first_10_indices]

# den_power_norm = den_power / np.trapezoid(den_power, freqs)
# vel_power_norm = vel_power / np.trapezoid(vel_power, freqs)
# mag_power_norm = mag_power / np.trapezoid(mag_power, freqs)