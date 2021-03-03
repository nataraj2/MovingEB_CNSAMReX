
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
       apz, azlo, azhi,level, finest_level, time, RKstep)
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
    real(rt) :: bct(3), cct(3), d1, d2, ddinv, cxm, cx0, cxp, cym, cy0, cyp, czm, cz0, czp
    real(rt) :: u1, v1, w1, u2, v2, w2, dudn, dvdn, dwdn
    real(rt) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, divu
    real(rt) :: tauxx, tauyy, tauzz, tauxy, tauxz, tauyz, tautmp
    real(rt), parameter :: twoThirds = 2.d0/3.d0
    real(rt) :: u, v, w, ublade, vblade, wblade, uvel, vvel, wvel, uwalldotn, un 
    integer :: myrank, ierr
    real(rt) :: x, y, z
    character(len=100) ::filename1
    integer, intent(in) :: level, finest_level
    real(rt), intent(in) :: time
    integer, intent(in) :: RKstep
    real(rt) :: dx(3)
    integer :: ii, jj, kk
    real(rt) :: a11, a12, a13, a21, a22, a23, a31, a32, a33, x0, y0, z0, xval, yval, zval
    real(rt) :: u_rhs1, u_rhs2, u_rhs3
    real(rt) :: v_rhs1, v_rhs2, v_rhs3
    real(rt) :: w_rhs1, w_rhs2, w_rhs3
    integer, parameter :: cluster_size = 3, matsize = 9, nvar = 3
    real(rt) :: weight
    real(rt) ::  ATA(matsize,matsize), ATb(matsize,nvar), dxi, dyi, dzi, vec(matsize), ebvel(nvar), xface, yface, zface
    integer  pivot(matsize), ok, iind, jind, var
    character (len=100) :: dir_name, file_name
    integer :: step
    logical :: dirExists

    dx = 1.0/dxinv
     
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

    ublade = 0.0d0
    vblade = 0.0d0 
    wblade = 50.0d0

    uvel = ublade
    vvel = vblade
    wvel = wblade

    dudx=0.0d0; dudy=0.0d0; dudz=0.0d0
    dvdx=0.0d0; dvdy=0.0d0; dvdz=0.0d0
    dwdx=0.0d0; dwdy=0.0d0; dwdz=0.0d0

    x0 =  (i+0.5d0)*dx(1) + bct(1)*dx(1)
    y0 =  (j+0.5d0)*dx(2) + bct(2)*dx(2)
    z0 =  (k+0.5d0)*dx(3) + bct(3)*dx(3)

	
    ATA = 0.0d0
    Atb = 0.0d0

    ebvel(1) = uvel
    ebvel(2) = vvel
    ebvel(3) = wvel
   

    do kk = k-cluster_size, k+cluster_size 
       do jj = j-cluster_size, j+cluster_size
          do ii = i-cluster_size, i+cluster_size
		if(ii.eq.i.and.jj.eq.j.and.kk.eq.k)then
		else
			if(vfrac(ii,jj,kk).ne.0.0d0)then
				if(vfrac(ii,jj,kk).lt.1-1e-10)then
    					cct =   centroid(ii,jj,kk,:)
					xval =  (ii+0.5d0)*dx(1) + cct(1)*dx(1) 
					yval =  (jj+0.5d0)*dx(2) + cct(2)*dx(2)
					zval =  (kk+0.5d0)*dx(3) + cct(3)*dx(3)

					bct = bcent(ii,jj,kk,:)
					xface = (ii+0.5d0)*dx(1) + bct(1)*dx(1)
					yface = (jj+0.5d0)*dx(2) + bct(2)*dx(2)
					zface = (kk+0.5d0)*dx(3) + bct(3)*dx(3)
				else
					xval =  (ii+0.5d0)*dx(1)
                                	yval =  (jj+0.5d0)*dx(2)
                                	zval =  (kk+0.5d0)*dx(3)
				endif

				weight = 1.0d0!/dsqrt((xval-x0)**2+(yval-y0)**2+(zval-z0)**2)

				dxi = (xval-x0)
				dyi = (yval-y0)
				dzi = (zval-z0)
		
				vec(1) = dxi;vec(2) = dyi;vec(3) = dzi			
				vec(4) = dxi**2;vec(5) = dyi**2;vec(6) = dzi**2
				vec(7) = dxi*dyi;vec(8) = dyi*dzi;vec(9) = dzi*dxi

		
				do iind = 1, matsize
			   		do jind = 1, matsize
						ATA(iind,jind) = ATA(iind,jind) + vec(iind)*vec(jind)	
			   		end do
				end do
			
				do var = 1, nvar
					do iind = 1, matsize
						ATb(iind,var) = ATb(iind,var) + (q(ii,jj,kk,var+1)-ebvel(var))*vec(iind)
					end do
				end do 

			endif
		endif
	   enddo
        enddo
     enddo

     call DGESV(matsize, nvar, ATA, matsize, pivot, ATb, matsize, ok)

     dudx = 0.0d0 
     dudy = 0.0d0
     dudz = 0.0d0

     dvdx = 0.0d0 
     dvdy = ATb(2,2) 
     dvdz = ATb(3,2)

     dwdx = 0.0d0 
     dwdy = ATb(2,3)
     dwdz = ATb(3,3)


    dudx = dudx*dx(1)
    dudy = dudy*dx(2)
    dudz = dudz*dx(3)

    dvdx = dvdx*dx(1)
    dvdy = dvdy*dx(2)
    dvdz = dvdz*dx(3)

    dwdx = dwdx*dx(1)
    dwdy = dwdy*dx(2)
    dwdz = dwdz*dx(3)


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
    if(abs(time-0.25d0).le.1e-5)then
    open(unit=10,file=filename1,position='append')		
    if(abs(x-1.0d0/(dxinv(1)*2.0d0)).le.1e-5.and.level .eq. finest_level .and. RKstep .eq. 2)then
        write(10,*)y,z,tauyy*dxinv(1),tauyz*dxinv(1),tauzz*dxinv(1), anrmy, anrmz
    endif
    close(10)
    endif

	
	step = int(time/2.5e-5)
    if(mod(step,100000).eq.0.0d0)then
	step = int(time/5e-5)
	write (dir_name,"('dir',i4.4)") step
        !print*,"dir name is ",trim(dir_name)


        WRITE(file_name,'(a,a,i4.4,a)')trim(dir_name),"/skin_friction.",myrank,".txt"
        inquire(file=file_name, exist=dirExists )
        if(dirExists.eqv..false.)call system('mkdir ' // dir_name)

    	open(unit=10,file=file_name,position='append')		
    	if(abs(x-1.0d0/(dxinv(1)*2.0d0)).le.1e-5.and.level .eq. finest_level .and. RKstep .eq. 2)then
        	write(10,*)y,z,tauyy*dxinv(1),tauyz*dxinv(1),tauzz*dxinv(1), anrmy, anrmz, apnorm
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
