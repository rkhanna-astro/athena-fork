import yt
import numpy as np
import matplotlib.pyplot as plt
from yt.units import mh
# import imageio
# from PIL import Image
import glob



spacing = np.linspace(0, 500, 100, dtype=int)
# spacing = [50]

time_step = 0

total_mass = []
mag_energy = []
kin_energy = []
time_in_myr = []

velocity_x = []

grid = 512
grid_x = np.linspace(0, 1, grid)
grid_y = np.linspace(0, 1, grid)

X,Y = np.meshgrid(grid_x, grid_y)

for space in spacing:
    # For wavelength 2L - I need 2-pc/1.58-km/s to get desired properties - 15uG, 1.58 km/s and 200 cm-3    ?
    # For wavelength L - I need 1-pc/1.0-km/s to get desired properties - 15uG, 1.58 km/s and 200 cm-3      ?
    # For wavelength L - I need 1-pc/1.0-km/s to get desired properties - 15uG, 1.58 km/s and 200 cm-3      ?
    unit_base={"length_unit": (1.0,"pc"), "time_unit": (1.0,"1.0 * pc / (1.0 * km/s)"), "mass_unit": (2.2254e34,"g")}
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

    density = grid['rho'].to("g/cm**3")[:, :, 0]
#     print(density)
    print(f"X dtype: {X.dtype}, Y dtype: {Y.dtype}, density dtype: {density.dtype}")

    physical_density = grid['rho'].to("g/cm**3")
#     print("Density", physical_density)

    physical_magnetic_field = grid[('gas', 'magnetic_field_y')].to("uG")
    print("Magnetic Field", physical_magnetic_field)

    physical_magnetic_field_x = grid[('gas', 'magnetic_field_x')].to("uG")
    print("Magnetic Field X", physical_magnetic_field_x)

#     mangetic_energy = grid[('gas', 'magnetic_energy_density')]
#     volume = grid[('gas', 'volume')]
    total_mag_energy = (physical_magnetic_field_x*(physical_magnetic_field_x/(8*3.14))).sum()
    mag_energy.append(total_mag_energy)

    kinetic_energy = grid[('gas', 'kinetic_energy_density')]
    total_kin_energy = (kinetic_energy).sum()
    kin_energy.append(total_kin_energy)
#     print("Magnetic Field X", total_mag_energy)

    vel_x = grid[('gas', 'velocity_y')].to('km/s')[:, :, 0]
    velocity_x.append(vel_x[128, 128])


    time_evolved = ds.current_time.to("Myr")
    print("Time", ds.current_time.to("Myr"))


#     alfven_speed = grid[('gas', 'alfven_speed')].to("km/s")
#     print("Alfven Speed", alfven_speed)


#     number_density = grid[('gas', 'number_density')].to("cm**-3")[:, :, 0]
#     print("Number Density", number_density)

    accurate_number_density = grid['rho'].to("g/cm**3") / (2.34 * mh)
    number_density = accurate_number_density.to("cm**-3")[:, :, 0]
    print(number_density.shape)
    print("Accurate Number Density", accurate_number_density)
    
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
    plt.pcolormesh(X, Y, number_density.T, cmap='hot', shading='auto', vmin = 100, vmax = 250)
    plt.title(f'Number Density Map (time={time_evolved:.4f})')


    plt.colorbar(label=r'Number Density (cm$^{-3}$)')
    plt.xlabel('x (1 pc)')
    plt.ylabel('y (1 pc)')
    plt.savefig(f'density_map_{time_evolved:.4f}.png', dpi=300, bbox_inches='tight')

    

# Plot the total mass
# plt.figure(figsize=(6, 5))
# # plt.plot(time_in_myr, mag_energy, label='magnetic energy density')
# plt.plot(time_in_myr, mag_energy, label='magnetic_energy_density')
# plt.yscale('log')
# plt.legend()
# plt.xlabel('Time (MYr)')
# plt.ylabel('Velocity_x (km/s)')
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
