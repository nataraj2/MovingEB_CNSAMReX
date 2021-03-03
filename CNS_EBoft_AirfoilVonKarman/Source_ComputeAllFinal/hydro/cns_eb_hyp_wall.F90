module cns_eb_hyp_wall_module
  use amrex_fort_module, only : rt=>amrex_real
  use amrex_error_module, only : amrex_abort
  use amrex_fi_mpi
  implicit none
  private
  public :: compute_hyp_wallflux

contains

  subroutine compute_hyp_wallflux (divw, i,j,k, rho, u, v, w, p, &
       axm, axp, aym, ayp, azm, azp,time,dx,level)
    use cns_physics_module, only : gamma
    use cns_module, only : smallp, smallr, umx, umy, umz
    use riemann_module, only : analriem
    use riemann_module_EB, only : analriem_EB
    integer, intent(in) :: i,j,k,level
    real(rt), intent(in) :: rho, u, v, w, p, axm, axp, aym, ayp, azm, azp
    real(rt), intent(out) :: divw(5)
    real(rt), intent(in) :: dx(3), time
    real(rt) :: x, y, z
    

    
    real(rt) :: apnorm, apnorminv, anrmx, anrmy, anrmz, un
    real(rt) :: flux(1,1,1,5)
    real(rt) :: velocity, uvel, vvel, wvel, ublade, vblade, wblade, velmag, rad, uwalldotn, twouwalldotn, theta, omega
    integer  :: myrank, ierr
    character(len=30) :: filename1, filename2, filename3 

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

    x = -7.0d0 + (i+0.5d0)*dx(1)
    y = -5.0d0 + (j+0.5d0)*dx(2)
    z = -0.3125d0 + (k+0.5d0)*dx(3)
  
    rad = dsqrt(x**2+y**2)
    omega = 41.45d0
    velmag = -2.51d0*3.1415d0/180.0d0*omega*cos(omega*time)*rad
    !velmag = 0.0d0
    theta = atan2(y,x)

	
    
    ublade = -velmag*sin(theta)
    vblade = velmag*cos(theta)
    wblade = 0.0d0

    twouwalldotn = 2.0d0*(ublade*anrmx+vblade*anrmy+wblade*anrmz) 
    uwalldotn = (ublade*anrmx+vblade*anrmy+wblade*anrmz) 


	call analriem_EB(gamma, smallp, smallr, 1, 1, 1, 1, &
         [rho], [ un], [p], [0.d0], [0.d0], &  ! fluid
         [rho], [twouwalldotn-un], [p], [0.d0], [0.d0], &  ! body
        flux, [1,1,1], [1,1,1], 2,3,4)

	call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)

	WRITE(filename1,'(a,i4.4,a)') "pressure_top.",myrank,".txt"
        WRITE(filename2,'(a,i4.4,a)') "pressure_bottom.",myrank,".txt"
        open(unit=10,file=filename1,position='append')
        open(unit=20,file=filename2,position='append')
        if(abs(z-dx(3)/2.0d0).le.1e-5.and.y+0.04086342658d0*x.ge.0.0d0 .and.level .eq. 3)then
                write(10,*)x,p,flux(1,1,1,2)
        endif
        if(abs(z-dx(3)/2.0d0).le.1e-5.and.y+0.04086342658d0*x.le.0.0d0 .and.level.eq. 3)then
                write(20,*)x,p,flux(1,1,1,2)
        endif
        close(10)
        close(20)

	uvel = u - un*anrmx + uwalldotn*anrmx
	vvel = v - un*anrmy + uwalldotn*anrmy
	wvel = w - un*anrmz + uwalldotn*anrmz

    divw = 0.d0
    
    divw(1)   = rho*(uvel*(axm-axp)+vvel*(aym-ayp)+wvel*(azm-azp))
    divw(umx) = (axm-axp)*flux(1,1,1,2) + uvel*flux(1,1,1,1)/apnorminv
    divw(umy) = (aym-ayp)*flux(1,1,1,2) + vvel*flux(1,1,1,1)/apnorminv
    divw(umz) = (azm-azp)*flux(1,1,1,2) + wvel*flux(1,1,1,1)/apnorminv
    divw(5)   = (flux(1,1,1,2)+flux(1,1,1,2)/0.4d0+0.5d0*rho*(uvel**2+vvel**2+wvel**2))*(uvel*(axm-axp)+vvel*(aym-ayp)+wvel*(azm-azp))

  end subroutine compute_hyp_wallflux

end module cns_eb_hyp_wall_module
