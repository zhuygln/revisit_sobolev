-- Minimal 1-line experiment: RESOLVED bound-bound treatment.
-- Fixed T (no radiative equilibrium), pure absorption, blackbody core.
sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = sedona_home.."/data/2level_atomdata.hdf5"

grid_type    = "grid_1D_sphere"
model_file   = "../minimal_1line.mod"
hydro_module = "homologous"

-- narrow window around 1215.5 A (nu0 = 2.4662e15 Hz), +-6% covers v_max/c = 1%
transport_nu_grid  = {2.30e15, 2.62e15, 2.0e-4, 1}
spectrum_nu_grid   = {2.30e15, 2.62e15, 5.0e-4, 1}
transport_radiative_equilibrium = 0
transport_steady_iterate        = 1

texp             = 20*3600.0*24.0
tstep_time_start = texp

-- blackbody core just inside the shell: L = 4 pi r^2 sigma T^4 at 2e4 K
core_n_emit      = 2e6
core_radius      = 1.728e14
core_temperature = 2.0e4
core_luminosity  = 3.40e42

opacity_grey_opacity        = 0
opacity_electron_scattering = 0
opacity_bound_bound         = 1
opacity_line_expansion      = 0
opacity_epsilon             = 1
line_velocity_width         = 1.0e7   -- 100 km/s Gaussian-ish width

output_write_radiation = 0
