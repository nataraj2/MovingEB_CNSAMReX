module cns_eb_hyp_wall_module
  use amrex_fort_module, only : rt=>amrex_real
  use amrex_error_module, only : amrex_abort
  implicit none
  private
  public :: compute_hyp_wallflux

contains

  subroutine compute_hyp_wallflux (divw, i,j,k, rho, u, v, w, p, &
       axm, axp, aym, ayp, azm, azp,time, dx, level, RKstep)
    use cns_physics_module, only : gamma
    use cns_module, only : smallp, smallr, umx, umy, umz
    use riemann_module, only : analriem
    use riemann_module_EB, only : analriem_EB 
    integer, intent(in) :: i,j,k
    real(rt), intent(in) :: rho, u, v, w, p, axm, axp, aym, ayp, azm, azp
    real(rt), intent(out) :: divw(5)
    real(rt), intent(in) :: time
    real(rt) :: x, y, z
    real(rt), intent(in):: dx(3)
    integer, intent(in) :: level, RKstep

    
    real(rt) :: apnorm, apnorminv, anrmx, anrmy, anrmz, un
    real(rt) :: flux(1,1,1,5)
    real(rt) :: velocity, uvel, vvel, wvel, ublade, vblade, wblade, velmag, rad, uwalldotn, twouwalldotn, theta
    real(rt) :: zval, zvec(6), timevec(6) 

    zvec(1) = 0.758945
    zvec(2) = 2.45757
    zvec(3) = 4.65736
    zvec(4) = 6.29964
    zvec(5) = 7.6947
    zvec(6) = 9.31447

    timevec(1) = 0.00d0
    timevec(2) = timevec(1) + 813*5e-5
    timevec(3) = timevec(2) + 818*5e-5
    timevec(4) = timevec(3) + 794*5e-5
    timevec(5) = timevec(4) + 794*5e-5
    timevec(6) = timevec(5) + 794*5e-5

    apnorm = sqrt((axm-axp)**2 + (aym-ayp)**2 + (azm-azp)**2)

    if (apnorm .eq. 0.d0) then
       print *, "compute_hyp_wallflux: ", i,j,k, axm, axp, aym, ayp, azm, azp
       flush(6)
       !call amrex_abort("compute_hyp_wallflux: we are in trouble.")
    end if

    apnorminv = 1.d0 / apnorm
    anrmx = (axm-axp) * apnorminv  ! pointing to the wall
    anrmy = (aym-ayp) * apnorminv
    anrmz = (azm-azp) * apnorminv

    un = u*anrmx + v*anrmy + w*anrmz

    ublade = 0.0d0
    vblade = 0.0d0
    zval = -5.0 + (k+0.5)*dx(3)
    if(zval .lt. 0.0d0) then
        vblade =  -0.5d0*0.2d0*2.0d0*3.14159265359d0*24.375d0*sin(2.0d0*3.14159265359d0*24.375d0*time)
    else if(zval .gt. zvec(1) .and. zval .lt. zvec(2) .and. & 
            time .gt. timevec(1) .and. time .lt. timevec(2)) then 
        vblade = 83.1d0
    else if(zval .gt. zvec(2) .and. zval .lt. zvec(3) .and. &
            time .gt. timevec(2) .and. time .lt. timevec(3)) then 
        vblade = -83.1d0
    else if(zval .gt. zvec(3) .and. zval .lt. zvec(4) .and. &
            time .gt. timevec(3) .and. time .lt. timevec(4)) then 
        vblade = 83.1d0
    else if(zval .gt. zvec(4) .and. zval .lt. zvec(5) .and. &
            time .gt. timevec(4) .and. time .lt. timevec(5)) then 
        vblade = -83.1d0
    else if(zval .gt. zvec(5) .and. zval .lt. zvec(6) .and. &
            time .gt. timevec(5) .and. time .lt. timevec(6)) then 
        vblade = 83.1d0
    endif
    vblade =  -0.5d0*2.0*3.14159265359d0*30.00d0*sin(2.0d0*3.14159265359d0*30.00d0*time)

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

    if(i==-2)then
      !divw = 0.0d0
    endif


  end subroutine compute_hyp_wallflux

end module cns_eb_hyp_wall_module
