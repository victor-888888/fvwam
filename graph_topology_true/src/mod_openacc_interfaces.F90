MODULE mod_openacc_interfaces 
  use mesh_mod
  use fre_dir_mod
  use mask_mod
  use WAM_SOURCE_MODULE
  use forcing_mod
  use state_mod 
  use wave_propag_mod  
  use params_mod
  use history_mod
  implicit none
  contains
SUBROUTINE openacc_init
  !------- copy variables in mesh_mod to device --------------------!

  !------- Copy variables in fre_dir to device --------------------!

  !--------- Copy variables in mask to device ------------------------!
  
  !--------- Copy variables in WAM_SOURCE to device ------------------------!

  !----------------Copy variables in state_mod to device--------------------!

  !----------------Copy variables in wave_propag to device--------------------!
  
  !----------------Copy variables in history_mod to device--------------------!
  if (spectrum_output_option) then
  end if

  !----------------Copy variables in wam_source to device--------------------!

  !----------------update sclar variable--------------------!

  !----------------Copy variables which can not be declared with create--------------------!


END SUBROUTINE

SUBROUTINE openacc_finalize
END SUBROUTINE
END MODULE
