import yt
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from yt.units import mh
# import imageio
# from PIL import Image
import glob
import matplotlib.patches as patches



spacing = np.linspace(0, 200, 200, dtype=int)
# spacing = [0,36]

time_step = 0

total_mass = []
mag_energy = []
kin_energy = []
time_in_myr = []

dcf = []
st = []
real = []

mag_x_low_y = []
mag_x_high_y = []
mean_backgroung_mag = []

plot_velocity = []
plot_dispersion = []
plot_alfven = []

grid = 512
grid_x = np.linspace(0, 1, grid)
grid_y = np.linspace(0, 1, grid)

X,Y = np.meshgrid(grid_x, grid_y)

for space in spacing:
    # For wavelength 2L - I need 2-pc/1.58-km/s to get desired properties - 15uG, 1.58 km/s and 200 cm-3    ?
    # For wavelength L - I need 1-pc/1.0-km/s to get desired properties - 15uG, 1.58 km/s and 200 cm-3      ?
    # For wavelength L - I need 1-pc/1.0-km/s to get desired properties - 15uG, 1.58 km/s and 200 cm-3      ?
    unit_base={"length_unit": (1.0,"pc"), "time_unit": (1.0,"1.0 * pc / (0.35 * km/s)"), "mass_unit": (3.532e34,"g")}
    ds = yt.load(f'LinWave.out2.{space:05d}.athdf', units_override=unit_base)

#     for unit in dir(ds.units):
#         print(unit)
#     print("Mass Unit", ds.mass_unit)
#     print("Density Unit", ds.density_unit)
#     print("Magnetic Unit", ds.magnetic_unit)
#     print("Time Unit", ds.time_unit)

    print(ds.field_list)

    left_edge = ds.domain_left_edge
    dims =  ds.domain_dimensions

    print(ds.domain_left_edge)
    print(ds.domain_dimensions)

    grid = ds.covering_grid(level=0, left_edge=left_edge, dims=dims)

    density = grid['rho'].to("g/cm**3")[:, :, 0].T
    print(f"X dtype: {X.dtype}, Y dtype: {Y.dtype}, density dtype: {density.dtype}")

    accurate_number_density = grid['rho'].to("g/cm**3") / (2.34 * mh)
    number_density = accurate_number_density.to("cm**-3")[:, :, 0].T
    print(number_density)
    print(number_density.shape)

    striations = number_density > 0

    physical_magnetic_field_y = grid[('gas', 'magnetic_field_y')].to("uG")[:, :, 0].T
    print("Magnetic Field Y", physical_magnetic_field_y[striations])

    if space == 0:
        mean_y = np.mean(physical_magnetic_field_y)

    physical_magnetic_field_x = grid[('gas', 'magnetic_field_x')].to("uG")[:, :, 0].T

    mag_x_low_y.append(np.max(physical_magnetic_field_x[5]))
    mag_x_high_y.append(np.max(physical_magnetic_field_x[250]))
    print("Magnetic Field X", np.max(physical_magnetic_field_x[254]))

    x_sel = X[striations]
    y_sel = Y[striations]

    velocity_x = grid[('gas', 'velocity_x')].to("km/s")[:, :, 0].T
    velocity_y = grid[('gas', 'velocity_y')].to("km/s")[:, :, 0].T

    print("Velocity_X", velocity_x[striations])
    print("Velocity_Y", velocity_y[striations])

    sigma_velocity = np.sqrt(np.sum(velocity_x[striations]**2 + velocity_y[striations]**2) / 256.0) 
    del_theta_dcf = (physical_magnetic_field_x[striations] / physical_magnetic_field_y[striations])
    sigma_angle_dcf = np.sqrt(np.sum(del_theta_dcf**2) / 256.0) 

    del_theta_st = ((physical_magnetic_field_y[striations] - mean_y) / mean_y)
    mean_backgroung_mag.append(np.mean(physical_magnetic_field_y[striations]))
    sigma_angle_st = np.sqrt(np.sum(del_theta_st**2) / 256.0)

    q = 0.5
    dcf_density = np.sqrt(4*3.14*np.mean(density[striations]))
    st_density = np.sqrt(2*3.14*np.mean(density[striations]))

    dcf_magnetic_field = ((q*dcf_density*sigma_velocity) / sigma_angle_dcf).to("uG")
    st_magnetic_field = 0
    if sigma_angle_st > 0:
        st_magnetic_field = ((st_density*sigma_velocity) / np.sqrt(sigma_angle_st)).to("uG")

    real.append(np.mean(physical_magnetic_field_y))
    dcf.append(dcf_magnetic_field)
    st.append(st_magnetic_field)

    # v_squared = (density * ((velocity_x)**2 + velocity_y**2)).sum() / density.sum()
    # v_rms = np.sqrt(v_squared).to("km/s")
    # plot_velocity.append(v_rms)
    # print("RMS velocity", v_rms)

    # # print("Velocity X", np.mean(velocity_xx))
    # # print("Velocity Y", np.mean(velocity_y))

    # sigma_x = np.std(velocity_x)
    # sigma_y = np.std(velocity_y)
    # sigma = np.sqrt(sigma_x**2 + sigma_y**2)
    # print("Velocity Dispersion", sigma)
    # plot_dispersion.append(sigma)

    # mean_velocity_x = np.sum(density * velocity_x) / np.sum(density)
    # mean_velocity_y = np.sum(density * velocity_y) / np.sum(density)
    
    # sigma_x = np.sqrt( (density * (velocity_x - mean_velocity_x)**2).sum()

#     mangetic_energy = grid[('gas', 'magnetic_energy_density')]
#     volume = grid[('gas', 'volume')]
    total_mag_energy = (physical_magnetic_field_x*(physical_magnetic_field_x/(8*3.14))).sum()
    mag_energy.append(total_mag_energy)

    kinetic_energy = grid[('gas', 'kinetic_energy_density')]
    total_kin_energy = (kinetic_energy).sum()
    kin_energy.append(total_kin_energy)
#     print("Magnetic Field X", total_mag_energy)

    time_evolved = ds.current_time.to("Myr")
    print("Time", ds.current_time.to("Myr"))


    alfven_speed = grid[('gas', 'alfven_speed')].to("km/s")[:,:, 0].T
    print("Alfven Speed", alfven_speed)
    print("Alfven Speed", np.mean(alfven_speed[0]))
    plot_alfven.append(np.mean(alfven_speed))


#     number_density = grid[('gas', 'number_density')].to("cm**-3")[:, :, 0]
#     print("Number Density", number_density)
    # print("Accurate Number Density", accurate_number_density)
    
    mass_of_mc = grid[('gas', 'cell_mass')].to("g").sum()
    total_mass.append(mass_of_mc)
    print("TOTAL Mass", mass_of_mc)

    time_in_myr.append(time_evolved)

#     print("Gamma", ds.gamma)
#     ds.gamma = 1.0

#     print(ds.derived_field_list)

#     pressure = grid[('gas', 'pressure')].to('dyne/cm**2')
#     print("Pressure", pressure)

#     iso_sound_speed = (pressure / density)**0.5
#     print("Sound speed", iso_sound_speed.to("km/s"))



    time_step += 1

#     if time_step < 20:
#         plt.figure(figsize=(6, 5))
#         plt.plot(grid_x, physical_magnetic_field_x[:, :, 0][128])
#         plt.xlabel('Time (MYr)')
#         plt.ylabel('Total Mass (g)') 
#         plt.show()
        # break
    plt.figure(figsize=(6, 5))
    print(X.dtype, Y.dtype, number_density.dtype)
    plt.pcolormesh(X, Y, number_density, cmap='gray', shading='auto', vmin = 100, vmax = 300)

    # plt.axhspan(0.3, 1.0, facecolor="red", alpha=0.15,
    #         edgecolor="black", linewidth=1.2)

    # Add label above the box (centered)
    # plt.text(
    #     0.5, 0.4, "Striations",    # x=0.5 (middle), y=0.75 (inside box)
    #     ha="center", va="center",
    #     fontsize=8,
    #     color="black",
    #     bbox=dict(facecolor="white", alpha=0.05, edgecolor="none")  # subtle background
    # )

    plt.title(f'Number Density Map (time={time_evolved:.4f})')
    plt.colorbar(label=r'Number Density (cm$^{-3}$)')
    plt.xlabel('x [pc]')
    plt.ylabel('y [pc]')
    plt.savefig(f'density_map_{time_evolved:.4f}.png', dpi=300, bbox_inches='tight')
    plt.close()

    # if time_step == 1:
    #     break

# Plot the total mass

# plt.figure(figsize=(6, 5))
# # plt.plot(time_in_myr, mag_x_low_y, label='y~0')
# # plt.plot(time_in_myr, mag_x_high_y, label='y~1')
# # plt.plot(time_in_myr, mean_backgroung_mag, label='mean Y')
# plt.plot(time_in_myr, real, label='real magnetic field')
# plt.plot(time_in_myr, dcf, label='DCF')
# plt.plot(time_in_myr, st, label='ST')
# # plt.yscale('log')
# plt.legend()
# plt.xlabel('Time (MYr)')
# plt.ylabel('Velocity RMS (km/s)')
# plt.show()

#     [('athena_pp', 'Bcc1'), ('athena_pp', 'Bcc2'), ('athena_pp', 'Bcc3'), ('athena_pp', 'cell_volume'), ('athena_pp', 'dx'), ('athena_pp', 'dy'), ('athena_pp', 'dz'), ('athena_pp', 'path_element_x'), ('athena_pp', 'path_element_y'), ('athena_pp', 'path_element_z'), ('athena_pp', 'rho'), ('athena_pp', 'vel1'), ('athena_pp', 'vel2'), ('athena_pp', 'vel3'), ('athena_pp', 'volume'), ('athena_pp', 'x'), ('athena_pp', 'y'), ('athena_pp', 'z'), 
#      ('gas', 'alfven_speed'), ('gas', 'angular_momentum_magnitude'), ('gas', 'angular_momentum_x'), ('gas', 'angular_momentum_y'), ('gas', 'angular_momentum_z'), ('gas', 'averaged_density'), ('gas', 'cell_mass'), ('gas', 'cell_volume'), ('gas', 'density'), ('gas', 'density_gradient_magnitude'), ('gas', 'density_gradient_x'), ('gas', 'density_gradient_y'), ('gas', 'density_gradient_z'), ('gas', 'dx'), ('gas', 'dy'), 
#      ('gas', 'dynamical_time'), ('gas', 'dz'), ('gas', 'four_velocity_magnitude'), ('gas', 'four_velocity_t'), ('gas', 'four_velocity_x'), ('gas', 'four_velocity_y'), ('gas', 'four_velocity_z'), ('gas', 'kinetic_energy_density'), ('gas', 'lorentz_factor'), ('gas', 'mach_alfven'), ('gas', 'magnetic_energy_density'), ('gas', 'magnetic_field_los'), ('gas', 'magnetic_field_magnitude'), ('gas', 'magnetic_field_strength'), ('gas', 'magnetic_field_x'), ('gas', 'magnetic_field_y'), ('gas', 'magnetic_field_z'), ('gas', 'magnetic_pressure'), 
#      ('gas', 'mass'), ('gas', 'mean_molecular_weight'), ('gas', 'momentum_density_x'), ('gas', 'momentum_density_y'), ('gas', 'momentum_density_z'), ('gas', 'momentum_x'), ('gas', 'momentum_y'), ('gas', 'momentum_z'), ('gas', 'number_density'), ('gas', 'path_element_x'), ('gas', 'path_element_y'), ('gas', 'path_element_z'), ('gas', 'relative_magnetic_field_x'), ('gas', 'relative_magnetic_field_y'), ('gas', 'relative_magnetic_field_z'), ('gas', 'relative_velocity_x'), ('gas', 'relative_velocity_y'), ('gas', 'relative_velocity_z'), ('gas', 'shear'), ('gas', 'specific_angular_momentum_magnitude'), 
#      ('gas', 'specific_angular_momentum_x'), ('gas', 'specific_angular_momentum_y'), ('gas', 'specific_angular_momentum_z'), ('gas', 'velocity_los'), ('gas', 'velocity_magnitude'), ('gas', 'velocity_x'), ('gas', 'velocity_y'), ('gas', 'velocity_z'), ('gas', 'volume'), ('gas', 'vorticity_magnitude'), ('gas', 'vorticity_squared'), ('gas', 'vorticity_x'), ('gas', 'vorticity_y'), ('gas', 'vorticity_z'), ('gas', 'x'), ('gas', 'y'), ('gas', 'z'), 
#      ('index', 'cell_volume'), ('index', 'cylindrical_radius'), ('index', 'cylindrical_theta'), ('index', 'cylindrical_z'), ('index', 'dx'), ('index', 'dy'), ('index', 'dz'), ('index', 'grid_indices'), ('index', 'grid_level'), ('index', 'ones'), ('index', 'ones_over_dx'), ('index', 'path_element_x'), ('index', 'path_element_y'), ('index', 'path_element_z'), ('index', 'radius'), ('index', 'spherical_phi'), ('index', 'spherical_radius'), ('index', 'spherical_theta'), ('index', 'virial_radius_fraction'), ('index', 'volume'), ('index', 'x'), ('index', 'y'), ('index', 'z'), ('index', 'zeros')]

#     plt.imshow(number_density.T, cmap='hot', 
#             extent=[X.min(), X.max(), Y.min(), Y.max()],
#             origin='lower', aspect='auto')
