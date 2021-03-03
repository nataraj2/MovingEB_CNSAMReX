module cns_eb_hyp_wall_module
  use amrex_fort_module, only : rt=>amrex_real
  use amrex_error_module, only : amrex_abort
  use amrex_fi_mpi
  implicit none
  private
  public :: compute_hyp_wallflux

contains

  subroutine compute_hyp_wallflux (divw, i,j,k, rho, u, v, w, p, &
       axm, axp, aym, ayp, azm, azp,time, dx, level, finest_level, RKstep)
    use cns_physics_module, only : gamma
    use cns_module, only : smallp, smallr, umx, umy, umz
    use riemann_module, only : analriem
    use riemann_module_EB, only : analriem_EB 
    integer, intent(in) :: i,j,k, level, finest_level
    real(rt), intent(in) :: rho, u, v, w, p, axm, axp, aym, ayp, azm, azp
    real(rt), intent(out) :: divw(5)
    real(rt), intent(in) :: dx(3), time
    real(rt) :: x, y, z

    
    real(rt) :: apnorm, apnorminv, anrmx, anrmy, anrmz, un
    real(rt) :: flux(1,1,1,5)
    real(rt) :: velocity, uvel, vvel, wvel, ublade, vblade, wblade, velmag, rad, uwalldotn, twouwalldotn, theta
    integer :: myrank, ierr
    character(len=100) ::filename1
    integer, intent(in) :: RKstep
    character (len=100) :: dir_name, file_name
    integer :: step
    logical :: dirExists

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

    x =  -0.15625d0 + (i+0.5d0)*dx(1)
    y = -5.0d0 + (j+0.5d0)*dx(2)
    z =  0.0d0 + (k+0.5d0)*dx(3)

    ublade = 0.0d0
    vblade = -0.2d0*0.4d0*2.0d0*3.14159265359d0*24.375d0*1.2d0*sin(2.0d0*3.14159265359d0*24.375d0*1.2d0*time)
    wblade = 0.0d0

    twouwalldotn = 2.0d0*(ublade*anrmx+vblade*anrmy+wblade*anrmz) 
    uwalldotn = (ublade*anrmx+vblade*anrmy+wblade*anrmz) 

    call analriem_EB(gamma, smallp, smallr, 1, 1, 1, 1, &
         [rho], [ un], [p], [0.d0], [0.d0], &  ! fluid
         [rho], [twouwalldotn-un], [p], [0.d0], [0.d0], &  ! body
        flux, [1,1,1], [1,1,1], 2,3,4)

	call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)

        WRITE(filename1,'(a,i4.4,a)') "pressure_top.",myrank,".txt"
   	if(abs(time-0.34188d0).le.1e-5)then
        open(unit=10,file=filename1,position='append')
        if(abs(x-dx(1)/2.0d0).le.1e-5.and.level .eq. finest_level .and. RKstep .eq. 2)then
                write(10,*)y,z,p,flux(1,1,1,2)
        endif
        close(10)
        endif



	step = int(time/5e-5)
	if(mod(step,1000000).eq.0.0d0)then
        step = int(time/5e-5)
        write (dir_name,"('dir',i4.4)") step
        !print*,"dir name is ",trim(dir_name)

        WRITE(file_name,'(a,a,i4.4,a)')trim(dir_name),"/pressure_surface.",myrank,".txt"
        inquire(file=file_name, exist=dirExists )
        if(dirExists.eqv..false.)call system('mkdir ' // dir_name)

        open(unit=10,file=file_name,position='append')
        if(abs(x-dx(1)/2.0d0).le.1e-5.and.level .eq. finest_level .and. RKstep .eq. 2)then
                write(10,*)y,z,p,flux(1,1,1,2),anrmy, anrmz, apnorm
        endif
        close(10)
        endif

	
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
    divw(umx) = (axm-axp)*flux(1,1,1,2) + uvel*flux(1,1,1,1)/apnorminv
    divw(umy) = (aym-ayp)*flux(1,1,1,2) + vvel*flux(1,1,1,1)/apnorminv
    divw(umz) = (azm-azp)*flux(1,1,1,2) + wvel*flux(1,1,1,1)/apnorminv
    divw(5)   = (flux(1,1,1,2)+flux(1,1,1,2)/0.4d0+0.5d0*rho*(uvel**2+vvel**2+wvel**2))*(uvel*(axm-axp)+vvel*(aym-ayp)+wvel*(azm-azp))

  end subroutine compute_hyp_wallflux

end module cns_eb_hyp_wall_module
