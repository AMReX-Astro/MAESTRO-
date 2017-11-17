subroutine average(lev,lo,hi,phi,p_lo,p_hi,phisum) bind (C,name="average")

  use base_state_geometry_module, only: max_radial_level, nr_fine
  use amrex_fort_module, only: amrex_spacedim

  integer         , intent (in   ) :: lev, lo(3), hi(3)
  integer         , intent (in   ) :: p_lo(3), p_hi(3)
  double precision, intent (in   ) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
  double precision, intent (inout) :: phisum(0:max_radial_level,0:nr_fine-1)

  integer :: i,j,k

  if (amrex_spacedim .eq. 2) then
     do j=lo(2),hi(2)
        phisum(lev,j) = sum(phi(lo(1):hi(1),j,0))
     end do
  else if (amrex_spacedim .eq. 3) then
     do k=lo(3),hi(3)
        phisum(lev,k) = sum(phi(lo(1):hi(1),lo(2):hi(2),k))
     end do
  end if

end subroutine average


subroutine divide_phisum_by_ncell(phisum,ncell) bind (C,name="divide_phisum_by_ncell")

  use base_state_geometry_module, only: max_radial_level, nr_fine

  double precision, intent(inout) :: phisum(0:max_radial_level,0:nr_fine-1)
  integer         , intent(in   ) ::  ncell(0:max_radial_level)

  integer :: n

  do n=0,max_radial_level
     phisum(:,n) = phisum(:,n) / ncell(n)
  end do

end subroutine divide_phisum_by_ncell
