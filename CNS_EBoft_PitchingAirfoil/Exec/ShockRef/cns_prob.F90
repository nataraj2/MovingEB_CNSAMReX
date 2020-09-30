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
  inflow_state(umx) = 1.226d0*256.9d0
  inflow_state(umy) = 0.0d0
  inflow_state(umz) = 0.0d0
  inflow_state(ueden) = 101325.0d0/0.4d0+ 0.5d0*1.226d0*256.9d0**2
  inflow_state(ueint) = 101325.0d0/0.4d0
  inflow_state(utemp) = 300.0d0

  interior_state(urho) = 1.226d0
  interior_state(umx) = 1.226d0*0.0d0
  interior_state(umy) = 0.d0
  interior_state(umz) = 0.d0
  interior_state(ueden) = 101325.0d0/0.4d0
  interior_state(ueint) = 101325.0d0/0.4d0
  interior_state(utemp) = 300.0d0

end subroutine amrex_probinit


subroutine cns_initdata(level, time, lo, hi, q, qlo, qhi, dx, prob_lo) bind(C, name="cns_initdata")
  use amrex_fort_module, only : rt => amrex_real
  use cns_physics_module, only : gamma, cv
  use cns_module, only : center, nvar, urho, umx, umy, umz, ueden, ueint, utemp, volfraction
  use probdata_module, only : p0,p1,rho0,rho1,v0,v1,x1, interior_state
  implicit none
  integer, intent(in) :: level, lo(3), hi(3), qlo(3), qhi(3)
  real(rt), intent(in) :: time
  real(rt), intent(inout) :: q(qlo(1):qhi(1), qlo(2):qhi(2), qlo(3):qhi(3),nvar)
  real(rt), intent(in) :: dx(3), prob_lo(3)
  
  integer :: i,j,k
  real(rt) :: x, y, z, p_l, p_r, rho_l, rho_r, u_l, u_r, rho, u, p 
   
  rho_l = 1.226d0
  rho_r = 0.944d0
  u_l = 459.54d0
  u_r = 289.86d0
  p_l = 101325.0d0
  p_r = 0.096d6

   do k = lo(3)-2, hi(3)+2
	z = prob_lo(3) + (k+0.5d0)*dx(3)
     do j = lo(2)-2, hi(2)+2
	y = prob_lo(2) + (j+0.5d0)*dx(2)
        do i = lo(1)-2, hi(1)+2
	   x = prob_lo(1) + (i+0.5d0)*dx(1)

		rho = 1.226d0
		u = 256.9d0
		p = 101325.0d0

             q(i,j,k,1) = rho
	     q(i,j,k,2) = rho*u
	     q(i,j,k,3) = 0.0d0
	     q(i,j,k,4) = 0.0d0
	     q(i,j,k,5) = p/(0.4d0) + 0.5d0*rho*u**2
	     q(i,j,k,6) = p/0.4d0
	     q(i,j,k,7) = p/(287.0d0*rho)
        end do
     end do
  end do

        
end subroutine cns_initdata
