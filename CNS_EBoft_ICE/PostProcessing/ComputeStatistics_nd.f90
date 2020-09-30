module transform_module

  use amrex_fort_module, only : amrex_real
  implicit none

  public

contains

  subroutine transform(lo, hi, dataIn, dataInlo, dataInhi, ncIn, dx, problo, volume_and_mass) bind(C, name="transform")

    integer, intent(in) :: lo(3),hi(3)
    integer, intent(in) :: dataInlo(3),dataInhi(3)
    integer, intent(in) :: ncIn

    real (kind=amrex_real),intent(in   ) :: dataIn(dataInlo(1):dataInhi(1),dataInlo(2):dataInhi(2),dataInlo(3):dataInhi(3),ncIn)
    real (kind=amrex_real), intent(in)   :: dx(3), problo(3)
    real (kind=amrex_real), intent(out) :: volume_and_mass(ncIn)
    integer :: i, j, k
    real(kind=amrex_real) :: x, y, z
    character(len=20) :: outfile

    volume_and_mass = 0.0d0
    do k=lo(3),hi(3)
	z = problo(3) + (float(k) + 0.5d0)*dx(3)
       do j=lo(2),hi(2)
	  y = problo(2) + (float(j) + 0.5d0)*dx(2)
          do i=lo(1),hi(1)
	     x = problo(1) + (float(i) + 0.5d0)*dx(1)
		volume_and_mass(1)=volume_and_mass(1)+dx(1)*dx(2)*dx(3)*dataIn(i,j,k,2)
		volume_and_mass(2)=volume_and_mass(2)+dataIn(i,j,k,1)*dx(1)*dx(2)*dx(3)*dataIn(i,j,k,2)
		volume_and_mass(3)=volume_and_mass(3)+dataIn(i,j,k,3)*dx(1)*dx(2)*dx(3)*dataIn(i,j,k,2)
          enddo
       enddo
    enddo

  end subroutine transform
end module transform_module
