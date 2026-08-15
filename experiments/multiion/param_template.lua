-- La II forest, 3850-3950 A. Placeholders: @BB@, @EXP@
sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = "../atom_multiion.hdf5"

grid_type    = "grid_1D_sphere"
model_file   = "../multiion.mod"
hydro_module = "homologous"

-- window 3790-3980 A: forest plus its 3000 km/s blue extension and margins.
-- 8 transport bins per 100 km/s Doppler width.
transport_nu_grid  = {7.50e14, 7.95e14, 4.17e-5, 1}
spectrum_nu_grid   = {7.50e14, 7.95e14, 1.0e-4, 1}
transport_radiative_equilibrium = 0
transport_steady_iterate        = 1

texp             = 86400.0
tstep_time_start = texp

core_n_emit      = 2e6
core_radius      = 8.64e12
core_temperature = 6000.0
core_luminosity  = 6.8937e37

opacity_grey_opacity        = 0
opacity_electron_scattering = 0
opacity_bound_bound         = @BB@
opacity_line_expansion      = @EXP@
opacity_epsilon             = 1
line_velocity_width         = 1.0e7

output_write_radiation = 0
