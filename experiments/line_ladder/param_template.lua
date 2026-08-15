-- Line-ladder experiment template. Placeholders: @ATOMFILE@, @BB@, @EXP@
sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = "../@ATOMFILE@"

grid_type    = "grid_1D_sphere"
model_file   = "../ladder.mod"
hydro_module = "homologous"

-- window: 10.2 eV line at 2.4662e15 Hz, forest extends 5% blueward,
-- troughs 1% further blue, margins beyond that.
transport_nu_grid  = {2.34e15, 2.70e15, 2.0e-4, 1}
spectrum_nu_grid   = {2.34e15, 2.70e15, 5.0e-4, 1}
transport_radiative_equilibrium = 0
transport_steady_iterate        = 1

texp             = 20*3600.0*24.0
tstep_time_start = texp

core_n_emit      = 2e6
core_radius      = 1.728e14
core_temperature = 2.0e4
core_luminosity  = 3.40e42

opacity_grey_opacity        = 0
opacity_electron_scattering = 0
opacity_bound_bound         = @BB@
opacity_line_expansion      = @EXP@
opacity_epsilon             = 1
line_velocity_width         = 1.0e7

output_write_radiation = 0
