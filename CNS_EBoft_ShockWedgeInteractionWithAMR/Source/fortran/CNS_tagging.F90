module cns_tagging_module
  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: cns_tag_denerror

contains

  subroutine cns_tag_denerror (lo, hi, tag, tlo, thi, rho, rlo, rhi, flag, flo, fhi, &
       dengrad, tagval, clearval, dx) bind(c,name='cns_tag_denerror')
    use iso_c_binding, only : c_char
    use amrex_ebcellflag_module, only : is_regular_cell, is_single_valued_cell
    integer, dimension(3), intent(in) :: lo, hi, tlo, thi, rlo, rhi, flo, fhi
    character(kind=c_char), intent(inout) :: tag(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3))
    real(kind=amrex_real),  intent(in)    :: rho(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    integer,                intent(in)   :: flag(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    real(kind=amrex_real), intent(in) :: dengrad
    real(kind=amrex_real)	      :: wx, wy, mean, wtilde
    character(kind=c_char), intent(in) :: tagval, clearval
    real(kind=amrex_real) , intent(in) :: dx(3)
    real(kind=amrex_real) :: schlieren,  drhodx,  drhody, drhodz

    integer :: i,j,k
    real(amrex_real) :: ax, ay, az

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (is_regular_cell(flag(i,j,k))) then
		wx=abs(rho(i+1,j,k)-2.0d0*rho(i,j,k)+rho(i-1,j,k))
                wy=abs(rho(i,j+1,k)-2.0d0*rho(i,j,k)+rho(i,j-1,k))
                mean=dsqrt((rho(i+1,j,k)+2.0d0*rho(i,j,k)+rho(i-1,j,k)**2+(rho(i,j+1,k)+2.0d0*rho(i,j,k)+rho(i,j-1,k)**2)))
                wtilde=dsqrt(wx**2+wy**2)/(mean+1e-10)
                if ( wtilde.gt.0.01d0) then
                !tag(i,j,k) = tagval
                endif


		drhodx = (rho(i+1,j,k)-rho(i-1,j,k))/(2.0d0*dx(1))
		drhody = (rho(i,j+1,k)-rho(i,j-1,k))/(2.0d0*dx(2))
		drhodz = (rho(i,j,k+1)-rho(i,j,k-1))/(2.0d0*dx(3))

		schlieren = dsqrt(drhodx**2+drhody**2+drhodz**2)

		if(schlieren .ge. 300.0d0)then
			tag(i,j,k) = tagval
		endif



	     endif
          end do
       end do
    end do

  end subroutine cns_tag_denerror

end module cns_tagging_module
