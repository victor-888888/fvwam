FOBJ=linked_list_mod.o        \
     array_mod.o              \
     hash_table_mod.o         \
     container.o              \
     timedelta_mod.o          \
     datetime_mod.o           \
     face_mod.o               \
     log_mod.o                \
     flogger.o                \
     fiona_mod.o              \
     string_mod.o             \
     string_actions_mod.o     \
     string_numerics_mod.o    \
     string.o                 \
     params_mod.o             \
     mesh_mod.o               \
     mask_mod.o               \
     forcing_mod.o            \
     state_mod.o              \
     time_scheme_mod.o        \
     wave_mod.o               \
     initial_mod.o            \
     read_forcing_mod.o       \
     const_mod.o              \
     time_mod.o               \
     history_mod.o            \
     fre_dir_mod.o            \
     wam_source_module.o      \
     wave_propag_mod.o        \
     restart_mod.o            \
     mod_mpi_interfaces.o     \
     mod_mpi_variables.o      \
     mod_mpi_reallocate.o     \
     mod_mpi_test_variables.o \
     mod_mpi_io_netcdf.o      \
     mod_mpi_test.o           \
     wave_main.o



wave.exe: $(FOBJ)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(FOBJ) $(LIBS)

linked_list_mod.o : linked_list_mod.F90
array_mod.o : array_mod.F90
hash_table_mod.o : hash_table_mod.F90 linked_list_mod.o
container.o : container.F90 linked_list_mod.o array_mod.o hash_table_mod.o 
timedelta_mod.o : timedelta_mod.F90
datetime_mod.o : datetime_mod.F90 timedelta_mod.o
face_mod.o : face_mod.F90
log_mod.o : log_mod.F90 hash_table_mod.o string.o face_mod.o
flogger.o : flogger.F90 log_mod.o
fiona_mod.o : fiona_mod.F90 flogger.o string.o container.o
string_mod.o : string_mod.F90 
string_actions_mod.o : string_actions_mod.F90 string_mod.o
string_numerics_mod.o : string_numerics_mod.F90 
string.o : string.F90 string_mod.o string_actions_mod.o string_numerics_mod.o
params_mod.o : params_mod.F90 const_mod.o mod_mpi_variables.o
mesh_mod.o : mesh_mod.F90 const_mod.o params_mod.o log_mod.o fiona_mod.o mod_mpi_interfaces.o
mask_mod.o : mask_mod.F90 params_mod.o mesh_mod.o fiona_mod.o log_mod.o
forcing_mod.o : forcing_mod.F90 const_mod.o mesh_mod.o log_mod.o read_forcing_mod.o time_mod.o params_mod.o string.o
state_mod.o : state_mod.F90 const_mod.o mesh_mod.o fre_dir_mod.o
time_scheme_mod.o : time_scheme_mod.F90 params_mod.o time_mod.o log_mod.o state_mod.o wave_propag_mod.o wam_source_module.o forcing_mod.o
wave_mod.o : wave_mod.F90 params_mod.o time_mod.o log_mod.o state_mod.o mesh_mod.o mask_mod.o history_mod.o time_scheme_mod.o forcing_mod.o initial_mod.o fre_dir_mod.o wam_source_module.o restart_mod.o 
initial_mod.o : initial_mod.F90  params_mod.o const_mod.o mesh_mod.o state_mod.o mask_mod.o forcing_mod.o fre_dir_mod.o restart_mod.o 
read_forcing_mod.o : read_forcing_mod.F90 const_mod.o params_mod.o
const_mod.o : const_mod.F90
time_mod.o : time_mod.F90  timedelta_mod.o datetime_mod.o hash_table_mod.o params_mod.o
history_mod.o : history_mod.F90 params_mod.o mesh_mod.o state_mod.o mask_mod.o forcing_mod.o fiona_mod.o log_mod.o string.o time_mod.o wave_propag_mod.o initial_mod.o wam_source_module.o const_mod.o fre_dir_mod.o mod_mpi_variables.o
fre_dir_mod.o : fre_dir_mod.F90 const_mod.o params_mod.o
wam_source_module.o : wam_source_module.F90 fre_dir_mod.o mesh_mod.o mask_mod.o const_mod.o wave_propag_mod.o initial_mod.o 
wave_propag_mod.o : wave_propag_mod.F90 const_mod.o fre_dir_mod.o params_mod.o mesh_mod.o  mask_mod.o
restart_mod.o : restart_mod.F90  params_mod.o mesh_mod.o time_mod.o fiona_mod.o log_mod.o string.o state_mod.o mod_mpi_interfaces.o 
mod_mpi_test_variables.o : mod_mpi_test_variables.F90 const_mod.o
mod_mpi_interfaces.o : mod_mpi_interfaces.F90 mod_mpi_test_variables.o mod_mpi_variables.o params_mod.o mod_mpi_reallocate.o fiona_mod.o
mod_mpi_variables.o : mod_mpi_variables.F90 const_mod.o
mod_mpi_reallocate.o : mod_mpi_reallocate.F90 mod_mpi_variables.o const_mod.o
mod_mpi_io_netcdf.o : mod_mpi_io_netcdf.F90 const_mod.o mod_mpi_test_variables.o
mod_mpi_test.o : mod_mpi_test.F90 mod_mpi_test_variables.o mod_mpi_variables.o mod_mpi_io_netcdf.o
ice_source_mod.o : ice_source_mod.F90 const_mod.o params_mod.o mesh_mod.o wave_propag_mod.o mod_mpi_interfaces.o mod_mpi_variables.o read_forcing_mod.o string.o time_mod.o log_mod.o
wave_main.o : wave_main.F90 params_mod.o log_mod.o wave_mod.o mod_mpi_test.o

