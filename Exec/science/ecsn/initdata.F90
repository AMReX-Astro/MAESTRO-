
module initdata_module

  use network, only: nspec, network_species_index
  use amrex_error_module
  use amrex_mempool_module, only: bl_allocate, bl_deallocate
  use amrex_constants_module
  use amrex_fort_module, only : amrex_spacedim, amrex_random
  use base_state_geometry_module, only: nr_fine, max_radial_level, center
  use meth_params_module, only: nscal, rho_comp, rhoh_comp, temp_comp, &
       spec_comp, pi_comp, prob_lo, prob_hi, sponge_start_factor, &
       sponge_center_density
  use eos_module
  use eos_type_module
  use probin_module, only: velpert_amplitude, velpert_steep, velpert_scale
  use fill_3d_data_module, only: put_1d_array_on_cart_sphr

  implicit none

  private

contains

  subroutine initdata(lev, time, lo, hi, &
       scal, scal_lo, scal_hi, nc_s, &
       vel, vel_lo, vel_hi, nc_v, &
       s0_init, p0_init, dx) bind(C, name="initdata")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: scal_lo(3), scal_hi(3), nc_s
    integer         , intent(in   ) :: vel_lo(3), vel_hi(3), nc_v
    double precision, intent(in   ) :: time
    double precision, intent(inout) :: scal(scal_lo(1):scal_hi(1), &
         scal_lo(2):scal_hi(2), &
         scal_lo(3):scal_hi(3), 1:nc_s)
    double precision, intent(inout) :: vel(vel_lo(1):vel_hi(1), &
         vel_lo(2):vel_hi(2), &
         vel_lo(3):vel_hi(3), 1:nc_v)
    double precision, intent(in   ) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(in   ) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dx(3)

    integer          :: i,j,k,r,n
    double precision :: rand, velpert_r_inner, velpert_r_outer

    ! random numbers between -1 and 1
    double precision :: alpha(3,3,3), beta(3,3,3), gamma(3,3,3)

    ! random numbers between 0 and 2*pi
    double precision :: phix(3,3,3), phiy(3,3,3), phiz(3,3,3)

    ! L2 norm of k
    double precision :: normk(3,3,3)

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             if (amrex_spacedim .eq. 1) then
                r = i
             else if (amrex_spacedim .eq. 2) then
                r = j
             else if (amrex_spacedim .eq. 3) then
                r = k
             end if

             ! set scalars using s0
             scal(i,j,k,rho_comp)  = s0_init(lev,r,rho_comp)
             scal(i,j,k,rhoh_comp) = s0_init(lev,r,rhoh_comp)
             scal(i,j,k,temp_comp) = s0_init(lev,r,temp_comp)
             scal(i,j,k,spec_comp:spec_comp+nspec-1) = &
                  s0_init(lev,r,spec_comp:spec_comp+nspec-1)

             ! initialize pi to zero for now
             scal(i,j,k,pi_comp) = 0.d0

          end do
       end do
    end do

    ! set velocity to zero
    vel(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nc_v) = 0.d0

  end subroutine initdata

  subroutine initdata_sphr(time, lo, hi, &
       scal, scal_lo, scal_hi, nc_s, &
       vel, vel_lo, vel_hi, nc_v, &
       s0_init, p0_init, &
       dx, r_cc_loc, r_edge_loc, &
       cc_to_r, ccr_lo, ccr_hi) bind(C, name="initdata_sphr")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: scal_lo(3), scal_hi(3), nc_s
    integer         , intent(in   ) :: vel_lo(3), vel_hi(3), nc_v
    integer         , intent(in   ) :: ccr_lo(3), ccr_hi(3)
    double precision, intent(in   ) :: time
    double precision, intent(inout) :: scal(scal_lo(1):scal_hi(1), &
                                            scal_lo(2):scal_hi(2), &
                                            scal_lo(3):scal_hi(3), nc_s)
    double precision, intent(inout) :: vel(vel_lo(1):vel_hi(1), &
                                           vel_lo(2):vel_hi(2), &
                                           vel_lo(3):vel_hi(3), nc_v)
    double precision, intent(in   ) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(in   ) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dx(3)
    double precision, intent(in   ) :: r_cc_loc  (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: cc_to_r   (ccr_lo(1):ccr_hi(1), &
         ccr_lo(2):ccr_hi(2),ccr_lo(3):ccr_hi(3))

    integer          :: i,j,k,r,n,comp
    double precision, pointer :: p0_cart(:,:,:,:)

    type (eos_t) :: eos_state

    double precision :: rand, velpert_r_inner, velpert_r_outer

    ! random numbers between -1 and 1
    double precision :: alpha(3,3,3), beta(3,3,3), gamma(3,3,3)

    ! random numbers between 0 and 2*pi
    double precision :: phix(3,3,3), phiy(3,3,3), phiz(3,3,3)

    ! L2 norm of k
    double precision :: normk(3,3,3)

    integer :: iloc, jloc, kloc

    ! cos and sin of (2*pi*kx/L + phix), etc
    double precision :: cx(3,3,3), cy(3,3,3), cz(3,3,3)
    double precision :: sx(3,3,3), sy(3,3,3), sz(3,3,3)

    ! location of center of star
    double precision :: xc(3)

    ! radius, or distance, to center of star
    double precision :: rloc

    ! the point we're at
    double precision :: xloc(3)

    ! perturbational velocity to add
    double precision :: vpert(3)

    integer :: ihe4

    ihe4 = network_species_index("helium-4")

    ! initialize the domain with the base state
    scal(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nc_s) = 0.d0

    ! if we are spherical, we want to make sure that p0 is good, since that is
    ! what is needed for HSE.  Therefore, we will put p0 onto a cart array and
    ! then initialize h from rho, X, and p0.
    call bl_allocate(p0_cart,lo,hi,1)

    ! initialize temp
    call put_1d_array_on_cart_sphr(lo,hi,scal(:,:,:,temp_comp),scal_lo,scal_hi,1, &
                                   s0_init(0,:,temp_comp),dx,0,0,r_cc_loc,r_edge_loc, &
                                   cc_to_r,ccr_lo,ccr_hi)

    ! initialize p0_cart
    call put_1d_array_on_cart_sphr(lo,hi,p0_cart,lo,hi,1,p0_init,dx,0,0,r_cc_loc,r_edge_loc, &
                                   cc_to_r,ccr_lo,ccr_hi)

    ! initialize species
    do comp = spec_comp, spec_comp+nspec-1
       call put_1d_array_on_cart_sphr(lo,hi,scal(:,:,:,comp),scal_lo,scal_hi,1, &
                                      s0_init(0,:,comp),dx,0,0,r_cc_loc,r_edge_loc, &
                                      cc_to_r,ccr_lo,ccr_hi)
    end do

    ! initialize rho as sum of partial densities rho*X_i
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             scal(i,j,k,rho_comp) = 0.d0
             do comp = spec_comp, spec_comp+nspec-1
                scal(i,j,k,rho_comp) = scal(i,j,k,rho_comp) + scal(i,j,k,comp)
             enddo
          enddo
       enddo
    enddo

    ! initialize (rho h) and T using the EOS
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state%T     = scal(i,j,k,temp_comp)
             eos_state%p     = p0_cart(i,j,k,1)
             eos_state%rho   = scal(i,j,k,rho_comp)
             eos_state%xn(:) = scal(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             call eos(eos_input_rp, eos_state)

             scal(i,j,k,rhoh_comp) = eos_state%rho*eos_state%h
             scal(i,j,k,temp_comp) = eos_state%T

          enddo
       enddo
    enddo

    call bl_deallocate(p0_cart)

    ! set velocity to zero
    vel(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nc_v) = 0.d0

    ! compute the radial bounds of the perturbation
    do i = 0, nr_fine-1
       if (s0_init(0,i,rho_comp) < sponge_start_factor*sponge_center_density) then
          velpert_r_outer = r_cc_loc(0,i)
          exit
       endif
    enddo

    do i = 0, nr_fine-1
       if (s0_init(0,i,spec_comp-1+ihe4)/s0_init(0,i,rho_comp) > 0.9d0) then
          velpert_r_inner = r_cc_loc(0,i)
          exit
       endif
    enddo

    ! adjust velpert_r_outer to be halfway between the base of the He layer and
    ! where the sponge turns on
    velpert_r_outer = velpert_r_inner + HALF*(velpert_r_outer - velpert_r_inner)

    velpert_r_outer = 7.5d6
    velpert_r_inner = 7.0d6
    if (amrex_spacedim .eq. 3) then

       ! generate random numbers
       ! random numbers are not currently supported
       ! use functions that result in numbers between (-1,1)
       do i=1,3
          do j=1,3
             do k=1,3
                rand = amrex_random()
                rand = 2.0d0*rand - 1.0d0
                ! rand = (.5)**i * (.7)**j * (.3)**k * (-1.)**i
                alpha(i,j,k) = rand
                rand = amrex_random()
                rand = 2.0d0*rand - 1.0d0
                ! rand = (.5)**i * (.3)**j * (.7)**k * (-1.)**j
                beta(i,j,k) = rand
                rand = amrex_random()
                rand = 2.0d0*rand - 1.0d0
                ! rand = (.3)**i * (.5)**j * (.7)**k * (-1.)**k
                gamma(i,j,k) = rand
                rand = amrex_random()
                ! rand = (.3)**i * (.7)**j * (.5)**k
                rand = 2.0d0*M_PI*rand
                phix(i,j,k) = rand
                rand = amrex_random()
                ! rand = (.7)**i * (.3)**j * (.5)**k
                rand = 2.0d0*M_PI*rand
                phiy(i,j,k) = rand
                rand = amrex_random()
                ! rand = (.7)**i * (.5)**j * (.3)**k
                rand = 2.0d0*M_PI*rand
                phiz(i,j,k) = rand
             enddo
          enddo
       enddo

       ! compute the norm of k
       do i=1,3
          do j=1,3
             do k=1,3
                normk(i,j,k) = sqrt(dble(i)**2+dble(j)**2+dble(k)**2)
             enddo
          enddo
       enddo

       ! now do the big loop over all points in the domain
       do iloc = lo(1),hi(1)
          do jloc = lo(2),hi(2)
             do kloc = lo(3),hi(3)

                ! set perturbational velocity to zero
                vpert = ZERO

                ! compute where we physically are
                xloc(1) = prob_lo(1) + (dble(iloc)+0.5d0)*dx(1) - center(1)
                xloc(2) = prob_lo(2) + (dble(jloc)+0.5d0)*dx(2) - center(2)
                xloc(3) = prob_lo(3) + (dble(kloc)+0.5d0)*dx(3) - center(3)

                ! compute distance to the center of the star
                rloc = ZERO
                do i=1,3
                   rloc = rloc + xloc(i)**2
                enddo
                rloc = sqrt(rloc)

                ! loop over the 27 combinations of fourier components
                do i=1,3
                   do j=1,3
                      do k=1,3
                         ! compute cosines and sines
                         cx(i,j,k) = cos(2.0d0*M_PI*dble(i)*xloc(1)/velpert_scale + phix(i,j,k))
                         cy(i,j,k) = cos(2.0d0*M_PI*dble(j)*xloc(2)/velpert_scale + phiy(i,j,k))
                         cz(i,j,k) = cos(2.0d0*M_PI*dble(k)*xloc(3)/velpert_scale + phiz(i,j,k))
                         sx(i,j,k) = sin(2.0d0*M_PI*dble(i)*xloc(1)/velpert_scale + phix(i,j,k))
                         sy(i,j,k) = sin(2.0d0*M_PI*dble(j)*xloc(2)/velpert_scale + phiy(i,j,k))
                         sz(i,j,k) = sin(2.0d0*M_PI*dble(k)*xloc(3)/velpert_scale + phiz(i,j,k))
                      enddo
                   enddo
                enddo

                ! loop over the 27 combinations of fourier components
                do i=1,3
                   do j=1,3
                      do k=1,3
                         ! compute contribution from perturbation velocity from each mode
                         vpert(1) = vpert(1) + &
                              (-gamma(i,j,k)*dble(j)*cx(i,j,k)*cz(i,j,k)*sy(i,j,k) &
                              +beta(i,j,k)*dble(k)*cx(i,j,k)*cy(i,j,k)*sz(i,j,k)) &
                              / normk(i,j,k)

                         vpert(2) = vpert(2) + &
                              (gamma(i,j,k)*dble(i)*cy(i,j,k)*cz(i,j,k)*sx(i,j,k) &
                              -alpha(i,j,k)*dble(k)*cx(i,j,k)*cy(i,j,k)*sz(i,j,k)) &
                              / normk(i,j,k)

                         vpert(3) = vpert(3) + &
                              ( -beta(i,j,k)*dble(i)*cy(i,j,k)*cz(i,j,k)*sx(i,j,k) &
                              +alpha(i,j,k)*dble(j)*cx(i,j,k)*cz(i,j,k)*sy(i,j,k)) &
                              / normk(i,j,k)
                      enddo
                   enddo
                enddo

                ! apply the cutoff function to the perturbational velocity
                do i=1,3
                   vpert(i) = velpert_amplitude*vpert(i) * &
                     HALF*(ONE + tanh((velpert_r_outer - velpert_steep - rloc)/velpert_steep))* &
                     HALF*(ONE + tanh((rloc - velpert_r_inner - velpert_steep)/velpert_steep))
                enddo

                ! add perturbational velocity to background velocity
                do i=1,3
                   vel(iloc,jloc,kloc,i) = vel(iloc,jloc,kloc,i) + vpert(i)
                enddo

             enddo
          enddo
       enddo

    endif

  end subroutine initdata_sphr


end module initdata_module
