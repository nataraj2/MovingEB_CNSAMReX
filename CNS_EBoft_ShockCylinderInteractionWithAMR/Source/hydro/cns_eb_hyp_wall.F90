module cns_eb_hyp_wall_module
  use amrex_fort_module, only : rt=>amrex_real
  use amrex_error_module, only : amrex_abort
  implicit none
  private
  public :: compute_hyp_wallflux

contains

  subroutine compute_hyp_wallflux (divw, i,j,k, rho, u, v, w, p, &
       axm, axp, aym, ayp, azm, azp,time)
    use cns_physics_module, only : gamma
    use cns_module, only : smallp, smallr, umx, umy, umz
    use riemann_module, only : analriem
    !use riemann_module_EB, only : analriem_EB 
    integer, intent(in) :: i,j,k
    real(rt), intent(in) :: rho, u, v, w, p, axm, axp, aym, ayp, azm, azp
    real(rt), intent(out) :: divw(5)
    real(rt), intent(in) :: time
    real(rt) :: x, y, z

    
    real(rt) :: apnorm, apnorminv, anrmx, anrmy, anrmz, un
    real(rt) :: flux(1,1,1,5)
    real(rt) :: velocity, uvel, vvel, wvel, ublade, vblade, wblade, velmag, rad, uwalldotn, twouwalldotn, theta

    apnorm = sqrt((axm-axp)**2 + (aym-ayp)**2 + (azm-azp)**2)

    if (apnorm .eq. 0.d0) then
       print *, "compute_hyp_wallflux: ", i,j,k, axm, axp, aym, ayp, azm, azp
       flush(6)
       call amrex_abort("compute_hyp_wallflux: we are in trouble.")
    end if

    apnorminv = 1.d0 / apnorm
    anrmx = (axm-axp) * apnorminv  ! pointing to the wall
    anrmy = (aym-ayp) * apnorminv
    anrmz = (azm-azp) * apnorminv

    un = u*anrmx + v*anrmy + w*anrmz

    ublade = 967.25d0
    vblade = 0.0d0
    wblade = 0.0d0

    twouwalldotn = 2.0d0*(ublade*anrmx+vblade*anrmy+wblade*anrmz) 
    uwalldotn = (ublade*anrmx+vblade*anrmy+wblade*anrmz) 


    call analriem(gamma, smallp, smallr, 1, 1, 1, 1, &
         [rho], [ un], [p], [0.d0], [0.d0], &  ! fluid
         [rho], [twouwalldotn-un], [p], [0.d0], [0.d0], &  ! body
        flux, [1,1,1], [1,1,1], 2,3,4)
	
	uvel = u - un*anrmx + uwalldotn*anrmx
	vvel = v - un*anrmy + uwalldotn*anrmy
	wvel = w - un*anrmz + uwalldotn*anrmz

    divw = 0.d0
    !divw(1)   = rho*(uvel*(axm-axp)+vvel*(aym-ayp)+wvel*(azm-azp))
    !divw(umx) = (axm-axp)*(p+rho*uvel**2) + (aym-ayp)*rho*uvel*vvel   + (azm-azp)*rho*uvel*wvel
    !divw(umy) = (axm-axp)*rho*vvel*uvel    + (aym-ayp)*(p+rho*vvel**2) + (azm-azp)*rho*vvel*wvel
    !divw(umz) = (axm-axp)*rho*wvel*uvel   + (aym-ayp)*rho*wvel*vvel   + (azm-azp)*(p+rho*wvel**2)
    !divw(5)   = (p+p/0.4d0+0.5d0*rho*(uvel**2+vvel**2+wvel**2))*(uvel*(axm-axp)+vvel*(aym-ayp)+wvel*(azm-azp))	

    divw(1)   = rho*(uvel*(axm-axp)+vvel*(aym-ayp)+wvel*(azm-azp))
    !divw(umx) = (axm-axp)*p + rho*uvel*ucyl1*(axm-axp)
    divw(umx) = (axm-axp)*p + uvel*flux(1,1,1,1)/apnorminv
    divw(umy) = (aym-ayp)*p + vvel*flux(1,1,1,1)/apnorminv
    divw(umz) = (azm-azp)*p + wvel*flux(1,1,1,1)/apnorminv
    divw(5)   = (p+p/0.4d0+0.5d0*rho*(uvel**2+vvel**2+wvel**2))*(uvel*(axm-axp)+vvel*(aym-ayp)+wvel*(azm-azp))

  end subroutine compute_hyp_wallflux

end module cns_eb_hyp_wall_module
