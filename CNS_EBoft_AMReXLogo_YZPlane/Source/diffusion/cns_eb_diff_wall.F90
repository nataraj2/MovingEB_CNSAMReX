
module cns_eb_diff_wall_module
  use amrex_fort_module, only : rt=>amrex_real
  use amrex_error_module, only : amrex_abort
  use cns_module, only : qvar, qu, qv, qw, qtemp
  use amrex_fi_mpi
  implicit none
  private
  public :: compute_diff_wallflux

contains

  !
  ! The wall is assumed to be adiabatic for now. Thus dTdn=0.
  ! Later we can support const T wall
  !
  ! We use no-slip boundary for velocities.
  !
  subroutine compute_diff_wallflux (divw, dxinv, i,j,k, &
       q, qlo, qhi, &
       lam, mu, xi, clo, chi, &
       bcent, blo, bhi, &
       centroid, centlo, centhi, &
       vfrac, vlo, vhi, &
       apx, axlo, axhi, &
       apy, aylo, ayhi, &
       apz, azlo, azhi, level, time, RKstep)
    integer, intent(in) :: i,j,k,qlo(3),qhi(3),clo(3),chi(3),axlo(3),axhi(3), &
         aylo(3),ayhi(3),azlo(3),azhi(3),blo(3),bhi(3), centlo(3), centhi(3), vlo(3), vhi(3)
    real(rt), intent(in) :: dxinv(3)
    real(rt), intent(out) :: divw(5)
    real(rt), intent(in) :: q  (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),qvar)
    real(rt), intent(in) :: lam(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(rt), intent(in) :: mu (clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(rt), intent(in) :: xi (clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(rt), intent(in) :: bcent( blo(1): bhi(1), blo(2): bhi(2), blo(3): bhi(3), 3)
    real(rt), intent(in) :: centroid(centlo(1):centhi(1),centlo(2):centhi(2),centlo(3):centhi(3),3)
    real(rt), intent(in) :: vfrac(vlo(1):vhi(1),vlo(2):vhi(2), vlo(3):vhi(3))
    real(rt), intent(in) :: apx  (axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(rt), intent(in) :: apy  (aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(rt), intent(in) :: apz  (azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))

    real(rt) :: dapx, dapy, dapz
    real(rt) :: apnorm, apnorminv, anrmx, anrmy, anrmz
    real(rt) :: xit, yit, zit, s
    integer :: ixit, iyit, izit, is
    real(rt) :: bct(3), d1, d2, ddinv, cxm, cx0, cxp, cym, cy0, cyp, czm, cz0, czp
    real(rt) :: u1, v1, w1, u2, v2, w2, dudn, dvdn, dwdn
    real(rt) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, divu
    real(rt) :: tauxx, tauyy, tauzz, tauxy, tauxz, tauyz, tautmp
    real(rt), parameter :: twoThirds = 2.d0/3.d0
    real(rt) :: u, v, w, ublade, vblade, wblade, uvel, vvel, wvel, uwalldotn, un 
    real(rt), intent(in) :: time
    real(rt) :: x, y, z
    character(len=100) ::filename1
    integer, intent(in) :: level
    integer, intent(in) :: RKstep
    integer :: myrank, ierr


    divw = 0.d0

    dapx = apx(i+1,j,k)-apx(i,j,k)
    dapy = apy(i,j+1,k)-apy(i,j,k)
    dapz = apz(i,j,k+1)-apz(i,j,k)

    apnorm = sqrt(dapx**2+dapy**2+dapz**2)

    if (apnorm .eq. 0.d0) then
       call amrex_abort("compute_diff_wallflux: we are in trouble.")
    end if

    apnorminv = 1.d0/apnorm
    anrmx = -dapx * apnorminv  ! unit vector pointing toward the wall
    anrmy = -dapy * apnorminv
    anrmz = -dapz * apnorminv

    ! The center of the wall
    bct = bcent(i,j,k,:)

    if (abs(anrmx).ge.abs(anrmy) .and. abs(anrmx).ge.abs(anrmz)) then
       ! y-z plane: x = const
       ! the equation for the line:  x = bct(1) - d*anrmx
       !                             y = bct(2) - d*anrmy
       !                             z = bct(3) - d*anrmz
       s = sign(1.0d0,-anrmx)
       is = nint(s)

       !
       ! the line intersects the y-z plane (x = s) at ...
       !
       d1 = (bct(1) - s) * (1.0d0/anrmx)  ! this is also the distance from wall to intersection
       yit = bct(2) - d1*anrmy
       zit = bct(3) - d1*anrmz
       iyit = j + nint(yit)
       izit = k + nint(zit)
       yit = yit - nint(yit)  ! shift so that the center of the nine cells are (0.,0.)
       zit = zit - nint(zit)
       !
       ! coefficents for quadratic interpolation
       cym = 0.5d0*yit*(yit-1.d0)
       cy0 = 1.d0-yit*yit
       cyp = 0.5d0*yit*(yit+1.d0)
       czm = 0.5d0*zit*(zit-1.d0)
       cz0 = 1.d0-zit*zit
       czp = 0.5d0*zit*(zit+1.d0)
       !
       ! interploation
       u1 = interp2d(cym,cy0,cyp,czm,cz0,czp, q(i+is,iyit-1:iyit+1,izit-1:izit+1,qu))
       v1 = interp2d(cym,cy0,cyp,czm,cz0,czp, q(i+is,iyit-1:iyit+1,izit-1:izit+1,qv))
       w1 = interp2d(cym,cy0,cyp,czm,cz0,czp, q(i+is,iyit-1:iyit+1,izit-1:izit+1,qw))

       !
       ! the line intersects the y-z plane (x = 2*s) at ...
       !
       d2 = (bct(1) - 2.d0*s) * (1.0d0/anrmx)
       yit = bct(2) - d2*anrmy
       zit = bct(3) - d2*anrmz
       iyit = j + nint(yit)
       izit = k + nint(zit)
       yit = yit - nint(yit)  ! shift so that the center of the nine cells are (0.,0.)
       zit = zit - nint(zit)
       !
       ! coefficents for quadratic interpolation
       cym = 0.5d0*yit*(yit-1.d0)
       cy0 = 1.d0-yit*yit
       cyp = 0.5d0*yit*(yit+1.d0)
       czm = 0.5d0*zit*(zit-1.d0)
       cz0 = 1.d0-zit*zit
       czp = 0.5d0*zit*(zit+1.d0)
       !
       ! interploation
       u2 = interp2d(cym,cy0,cyp,czm,cz0,czp, q(i+2*is,iyit-1:iyit+1,izit-1:izit+1,qu))
       v2 = interp2d(cym,cy0,cyp,czm,cz0,czp, q(i+2*is,iyit-1:iyit+1,izit-1:izit+1,qv))
       w2 = interp2d(cym,cy0,cyp,czm,cz0,czp, q(i+2*is,iyit-1:iyit+1,izit-1:izit+1,qw))

    else if (abs(anrmy).ge.abs(anrmx) .and. abs(anrmy).ge.abs(anrmz)) then
       ! z-x plane
       s = sign(1.0d0,-anrmy)
       is = nint(s)

       d1 = (bct(2) - s) * (1.0d0/anrmy)
       xit = bct(1) - d1*anrmx
       zit = bct(3) - d1*anrmz
       ixit = i + nint(xit)
       izit = k + nint(zit)
       xit = xit - nint(xit)
       zit = zit - nint(zit)

       cxm = 0.5d0*xit*(xit-1.d0)
       cx0 = 1.d0-xit*xit
       cxp = 0.5d0*xit*(xit+1.d0)
       czm = 0.5d0*zit*(zit-1.d0)
       cz0 = 1.d0-zit*zit
       czp = 0.5d0*zit*(zit+1.d0)

       u1 = interp2d(cxm,cx0,cxp,czm,cz0,czp, q(ixit-1:ixit+1,j+is,izit-1:izit+1,qu))
       v1 = interp2d(cxm,cx0,cxp,czm,cz0,czp, q(ixit-1:ixit+1,j+is,izit-1:izit+1,qv))
       w1 = interp2d(cxm,cx0,cxp,czm,cz0,czp, q(ixit-1:ixit+1,j+is,izit-1:izit+1,qw))

       d2 = (bct(2) - 2.d0*s) * (1.0d0/anrmy)
       xit = bct(1) - d2*anrmx
       zit = bct(3) - d2*anrmz
       ixit = i + nint(xit)
       izit = k + nint(zit)
       xit = xit - nint(xit)
       zit = zit - nint(zit)

       cxm = 0.5d0*xit*(xit-1.d0)
       cx0 = 1.d0-xit*xit
       cxp = 0.5d0*xit*(xit+1.d0)
       czm = 0.5d0*zit*(zit-1.d0)
       cz0 = 1.d0-zit*zit
       czp = 0.5d0*zit*(zit+1.d0)

       u2 = interp2d(cxm,cx0,cxp,czm,cz0,czp, q(ixit-1:ixit+1,j+2*is,izit-1:izit+1,qu))
       v2 = interp2d(cxm,cx0,cxp,czm,cz0,czp, q(ixit-1:ixit+1,j+2*is,izit-1:izit+1,qv))
       w2 = interp2d(cxm,cx0,cxp,czm,cz0,czp, q(ixit-1:ixit+1,j+2*is,izit-1:izit+1,qw))

    else
       ! x-y plane
       s = sign(1.d0,-anrmz)
       is = nint(s)

       d1 = (bct(3) - s) * (1.0d0/anrmz)
       xit = bct(1) - d1*anrmx
       yit = bct(2) - d1*anrmy
       ixit = i + nint(xit)
       iyit = j + nint(yit)
       xit = xit - nint(xit)
       yit = yit - nint(yit)

       cxm = 0.5d0*xit*(xit-1.d0)
       cx0 = 1.d0-xit*xit
       cxp = 0.5d0*xit*(xit+1.d0)
       cym = 0.5d0*yit*(yit-1.d0)
       cy0 = 1.d0-yit*yit
       cyp = 0.5d0*yit*(yit+1.d0)

       u1 = interp2d(cxm,cx0,cxp,cym,cy0,cyp, q(ixit-1:ixit+1,iyit-1:iyit+1,k+is,qu))
       v1 = interp2d(cxm,cx0,cxp,cym,cy0,cyp, q(ixit-1:ixit+1,iyit-1:iyit+1,k+is,qv))
       w1 = interp2d(cxm,cx0,cxp,cym,cy0,cyp, q(ixit-1:ixit+1,iyit-1:iyit+1,k+is,qw))
       
       d2 = (bct(3) - 2.d0*s) * (1.0d0/anrmz)
       xit = bct(1) - d2*anrmx
       yit = bct(2) - d2*anrmy
       ixit = i + nint(xit)
       iyit = j + nint(yit)
       xit = xit - nint(xit)
       yit = yit - nint(yit)

       cxm = 0.5d0*xit*(xit-1.d0)
       cx0 = 1.d0-xit*xit
       cxp = 0.5d0*xit*(xit+1.d0)
       cym = 0.5d0*yit*(yit-1.d0)
       cy0 = 1.d0-yit*yit
       cyp = 0.5d0*yit*(yit+1.d0)

       u2 = interp2d(cxm,cx0,cxp,cym,cy0,cyp, q(ixit-1:ixit+1,iyit-1:iyit+1,k+2*is,qu))
       v2 = interp2d(cxm,cx0,cxp,cym,cy0,cyp, q(ixit-1:ixit+1,iyit-1:iyit+1,k+2*is,qv))
       w2 = interp2d(cxm,cx0,cxp,cym,cy0,cyp, q(ixit-1:ixit+1,iyit-1:iyit+1,k+2*is,qw))
       
    end if

    ublade = 0.0d0
    !vblade = -83.1d0
	vblade =  -0.5d0*0.2d0*2.0d0*3.14159265359d0*24.375d0*sin(2.0d0*3.14159265359d0*24.375d0*time)
    wblade = 0.0d0
    
    uvel = ublade!u - un*anrmx + uwalldotn*anrmx
    vvel = vblade!v - un*anrmy + uwalldotn*anrmy
    wvel = wblade!w - un*anrmz + uwalldotn*anrmz


    !
    ! compute derivatives on the wall given that velocity is zero on the wall.
    !
    ddinv = 1.d0/(d1*d2*(d2-d1))
    dudn = -ddinv*(d2*d2*(u1-uvel)-d1*d1*(u2-uvel))  ! note that the normal vector points toward the wall
    dvdn = -ddinv*(d2*d2*(v1-vvel)-d1*d1*(v2-vvel))
    dwdn = -ddinv*(d2*d2*(w1-wvel)-d1*d1*(w2-wvel))

    !dudn = (uvel-u1)/d1
    !dvdn = (vvel-v1)/d1
    !dwdn = (wvel-w1)/d1
    !
    ! transform them to d/dx, d/dy and d/dz given transverse derivatives are zero
    dudx = dudn * anrmx
    dudy = dudn * anrmy
    dudz = dudn * anrmz
    !
    dvdx = dvdn * anrmx
    dvdy = dvdn * anrmy
    dvdz = dvdn * anrmz
    !
    dwdx = dwdn * anrmx
    dwdy = dwdn * anrmy
    dwdz = dwdn * anrmz

    divu = dudx+dvdy+dwdz
    tautmp = (xi(i,j,k)-twoThirds*mu(i,j,k))*divu
    tauxx = mu(i,j,k)*2.d0*dudx + tautmp
    tauyy = mu(i,j,k)*2.d0*dvdy + tautmp
    tauzz = mu(i,j,k)*2.d0*dwdz + tautmp
    tauxy = mu(i,j,k)*(dudy+dvdx)
    tauxz = mu(i,j,k)*(dudz+dwdx)
    tauyz = mu(i,j,k)*(dwdy+dvdz)

    x =  -0.15625d0 + (i+0.5d0)/dxinv(1)
    y = -5.0d0 + (j+0.5d0)/dxinv(1)
    z =  0.0d0 + (k+0.5d0)/dxinv(1)

    call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)

    WRITE(filename1,'(a,i4.4,a)') "skin_friction.",myrank,".txt"
    if(abs(time-0.205128d0).le.1e-5)then
    open(unit=10,file=filename1,position='append')
    if(abs(x-1.0d0/(dxinv(1)*2.0d0)).le.1e-5.and.level .eq. 3 .and. RKstep .eq. 2)then
        write(10,*)y,z,tauyy*dxinv(1),tauyz*dxinv(1),tauzz*dxinv(1), anrmy, anrmz
    endif
    close(10)
    endif


    ! dx == dy == dz
    divw(2) = dxinv(1) * (dapx*tauxx + dapy*tauxy + dapz*tauxz)
    divw(3) = dxinv(1) * (dapx*tauxy + dapy*tauyy + dapz*tauyz)
    divw(4) = dxinv(1) * (dapx*tauxz + dapy*tauyz + dapz*tauzz)
 
  end subroutine compute_diff_wallflux

  
  real(rt) function interp2d(cym,cy0,cyp,czm,cz0,czp,v)
    real(rt), intent(in) :: cym,cy0,cyp,czm,cz0,czp,v(3,3)
    interp2d = czm*(cym*v(1,1) + cy0*v(2,1) + cyp*v(3,1)) &
         +     cz0*(cym*v(1,2) + cy0*v(2,2) + cyp*v(3,2)) &
         +     czp*(cym*v(1,3) + cy0*v(2,3) + cyp*v(3,3))
  end function interp2d

end module cns_eb_diff_wall_module
