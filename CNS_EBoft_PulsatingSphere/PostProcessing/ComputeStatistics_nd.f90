module transform_module

  use amrex_fort_module, only : amrex_real
  implicit none

  public

contains

  subroutine transform(lo, hi, dataIn, dataInlo, dataInhi, ncIn, sOut, sOutlo, sOuthi, ncOut, dx, problo, sigma_xvel, filenum) bind(C, name="transform")

    integer, intent(in) :: lo(3),hi(3)
    integer, intent(in) :: dataInlo(3),dataInhi(3)
    integer, intent(in) :: sOutlo(3),sOuthi(3)
    integer, intent(in) :: ncIn, ncOut, filenum

    real (kind=amrex_real),intent(in   ) :: dataIn(dataInlo(1):dataInhi(1),dataInlo(2):dataInhi(2),dataInlo(3):dataInhi(3),ncIn)
    real (kind=amrex_real),intent(inout) :: sOut(sOutlo(1):sOuthi(1),sOutlo(2):sOuthi(2),sOutlo(3):sOuthi(3),ncOut)
    real (kind=amrex_real), intent(in)   :: dx(3), problo(3)
    real (kind=amrex_real), intent(out) :: sigma_xvel(0:63)

    integer :: i, j, k
    real(kind=amrex_real) :: x, y, z
    character(len=20) :: outfile

    write(outfile, "(A8,I0.5,A4)") "pressure", filenum,".txt"
    open(unit=10,file=outfile,position='append')

    ! This is an example pointwise transformation
    ! Here, dataIn(i,j,k,1LncIn) contains data from the plotfile, stacked up in the order that the 
    ! user called the function with.  The output data has (hardwired) ncOut components
  
    sigma_xvel = 0.0d0
    do k=lo(3),hi(3)
	z = problo(3) + (float(k) + 0.5d0)*dx(3)
       do j=lo(2),hi(2)
	  y = problo(2) + (float(j) + 0.5d0)*dx(2)
          do i=lo(1),hi(1)
	     x = problo(1) + (float(i) + 0.5d0)*dx(1)
		if(x.gt.0.0d0.and.abs(x).lt.dx(1).and.y.gt.0.0d0.and.abs(y).lt.dx(2).and.z.ge.0.01d0)then
		write(10, *)z,dataIn(i,j,k,ncIn)
		endif
          enddo
       enddo
    enddo

	close(10)

  end subroutine transform
end module transform_module
