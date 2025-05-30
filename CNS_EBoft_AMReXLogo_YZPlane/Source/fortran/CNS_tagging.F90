module cns_tagging_module
  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: cns_tag_denerror

contains

  subroutine cns_tag_denerror (lo, hi, tag, tlo, thi,& 
                               yvel, yvel_lo, yvel_hi,& 
                               zvel, zvel_lo, zvel_hi,& 
                               flag, flo, fhi, &
                                  dengrad, tagval, clearval, dx) bind(c,name='cns_tag_denerror')
    use iso_c_binding, only : c_char
    use amrex_ebcellflag_module, only : is_regular_cell, is_single_valued_cell
    integer, dimension(3), intent(in) :: lo, hi, tlo, thi, yvel_lo, yvel_hi, zvel_lo, zvel_hi, flo, fhi
    character(kind=c_char), intent(inout) :: tag(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3))
    real(kind=amrex_real),  intent(in)    :: yvel(yvel_lo(1):yvel_hi(1),yvel_lo(2):yvel_hi(2),yvel_lo(3):yvel_hi(3))
    real(kind=amrex_real),  intent(in)    :: zvel(zvel_lo(1):zvel_hi(1),zvel_lo(2):zvel_hi(2),zvel_lo(3):zvel_hi(3))
    integer,                intent(in)   :: flag(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    real(kind=amrex_real), intent(in) :: dengrad
    real(kind=amrex_real)          :: vortmag
    character(kind=c_char), intent(in) :: tagval, clearval

    integer :: i,j,k
    real(amrex_real) :: x, y, z, dvdz, dwdy
    real(amrex_real), intent(in) :: dx(3)

    do k = lo(3), hi(3)
       z = -5.0 + (k+0.5d0)*dx(3)    
       do j = lo(2), hi(2)
          y = -3.0d0 + (j+0.5d0)*dx(2)    
          do i = lo(1), hi(1)
          !x = 0.0d0 + (i+0.5d0)*dx(1)    
             if (is_regular_cell(flag(i,j,k))) then
                dvdz=(yvel(i,j,k+1)-yvel(i,j,k-1))/(2.0*dx(3))
                dwdy=(zvel(i,j+1,k)-zvel(i,j-1,k))/(2.0*dx(2))
                vortmag = abs(dvdz - dwdy)
                if ( z .ge. -3.0d0 .and. z .le. 12.0d0 .and. y.ge.1.0d0 .and. y .le. 3.0d0)then
                    !tag(i,j,k) = tagval
                endif
                !if(vortmag .ge. 300.0d0 .and. y .gt. -1.5d0 .and. y .lt. 5.5d0) then
                if(vortmag .ge. 300.0d0) then
                    tag(i,j,k) = tagval
                endif
             endif
          end do
       end do
    end do

  end subroutine cns_tag_denerror


end module cns_tagging_module
