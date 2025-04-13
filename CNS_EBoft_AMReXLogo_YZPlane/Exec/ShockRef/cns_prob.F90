module probdata_module
  use amrex_fort_module, only : rt => amrex_real
  use cns_module, only : center, nvar, urho, umx, umy, umz, ueden, ueint, utemp
  implicit none
  real(rt), save :: p0   = 101325.0d0
  real(rt), save :: p1   = 101325.0d0
  real(rt), save :: rho0 = 1.226d0
  real(rt), save :: rho1 = 1.226d0
  real(rt), save :: v0   = 0.d0
  real(rt), save :: v1   = 0.d0
  real(rt), save :: x1   = 4.2
 
  real(rt), save :: inflow_state(nvar)
  real(rt), save :: interior_state(nvar)
 
end module probdata_module


subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)
  use amrex_fort_module, only : rt => amrex_real
  use amrex_parmparse_module
  use probdata_module
  implicit none
  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(*), probhi(*)
  type(amrex_parmparse) :: pp
  real(rt) :: vel = 50.0;
  call amrex_parmparse_build(pp,"prob")
  call pp%query("p0", p0)
  call pp%query("p1", p1)
  call pp%query("rho0", rho0)
  call pp%query("rho1", rho1)
  call pp%query("v0", v0)
  call pp%query("v1", v1)
  call pp%query("x1", x1)
  call amrex_parmparse_destroy(pp)

  inflow_state(urho) = 1.226d0
  inflow_state(umx) = 0.0d0
  inflow_state(umy) = 0.0d0
  inflow_state(umz) = 1.226d0*vel
  inflow_state(ueden) = 101325.0d0/0.4d0 + 0.5*1.226*vel*vel!+ 0.5d0*1.226d0*5.0d0**2
  inflow_state(ueint) = 101325.0d0/0.4d0
  inflow_state(utemp) = 300.0d0

  interior_state(urho) = 1.226d0
  interior_state(umx) = 0.0d0
  interior_state(umy) = 0.0d0
  interior_state(umz) = 0.0d0
  interior_state(ueden) = 101325.0d0/0.4d0
  interior_state(ueint) = 101325.0d0/0.4d0
  interior_state(utemp) = 300.0d0

end subroutine amrex_probinit


subroutine cns_initdata(level, time, lo, hi, q, qlo, qhi, dx, prob_lo) bind(C, name="cns_initdata")
  use amrex_fort_module, only : rt => amrex_real
  use cns_physics_module, only : gamma, cv
  use cns_module, only : center, nvar, urho, umx, umy, umz, ueden, ueint, utemp, volfraction
  use probdata_module, only : p0,p1,rho0,rho1,v0,v1,x1, interior_state, inflow_state
  implicit none
  integer, intent(in) :: level, lo(3), hi(3), qlo(3), qhi(3)
  real(rt), intent(in) :: time
  real(rt), intent(inout) :: q(qlo(1):qhi(1), qlo(2):qhi(2), qlo(3):qhi(3),nvar)
  real(rt), intent(in) :: dx(3), prob_lo(3)
  
  integer :: i,j,k
  real(rt) :: x, y, z, p_l, p_r, rho_l, rho_r, u_l, u_r, rho, w, p 
   
	rho = 1.226;
	p = 101325.0;
	
	do k = lo(3)-2, hi(3)+2
		z = prob_lo(2) + (k+0.5)*dx(2)/2**level
     	do j = lo(2)-2, hi(2)+2
 	       do i = lo(1)-2, hi(1)+2
	 
 			!if(z .lt. -1.0d0) then
			!	q(i,j,k,:) = inflow_state
			!else
				q(i,j,k,:) = interior_state	
			!endif

        end do
     end do
  end do

        
end subroutine cns_initdata
