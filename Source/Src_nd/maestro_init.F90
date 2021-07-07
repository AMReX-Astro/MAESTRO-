module maestro_init_module

  use amrex_error_module
  use network, only: network_init, nspec, short_spec_names, short_aux_names
  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use amrex_fort_module, only: amrex_spacedim
  use model_parser_module
  use meth_params_module, only: rho_comp, rhoh_comp, spec_comp, temp_comp, pi_comp, &
       nscal, small_dens, small_temp
  use runtime_init_module
  implicit none

  private

contains

  subroutine maestro_network_init() bind(C, name="maestro_network_init")
    ! Binds to C function ``maestro_network_init``

    use actual_rhs_module, only: actual_rhs_init

    call network_init()
    call actual_rhs_init()

  end subroutine maestro_network_init

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine maestro_eos_init() bind(C, name="maestro_eos_init")
    ! Binds to C function ``maestro_eos_init``
    
    use eos_module, only: eos_init

    integer :: ioproc

    !---------------------------------------------------------------------
    ! other initializations
    !---------------------------------------------------------------------

    ! This is a routine which links to the C++ ParallelDescriptor class

    call bl_pd_is_ioproc(ioproc)

    !---------------------------------------------------------------------
    ! safety checks
    !---------------------------------------------------------------------

    if (small_dens <= 0.d0) then
        if (ioproc == 1) then
           call amrex_warning("Warning:: small_dens has not been set, defaulting to 1.d-200.")
        endif
        small_dens = 1.d-200
     endif
 
     if (small_temp <= 0.d0) then
        if (ioproc == 1) then
           call amrex_warning("Warning:: small_temp has not been set, defaulting to 1.d-200.")
        endif
        small_temp = 1.d-200
     endif
 
     ! Note that the EOS may modify our choices because of its
     ! internal limitations, so the small_dens and small_temp
     ! may be modified coming back out of this routine.
 
     call eos_init(small_dens=small_dens, small_temp=small_temp)

  end subroutine maestro_eos_init

end module maestro_init_module
