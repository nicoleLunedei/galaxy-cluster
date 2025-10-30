!******************************************************************
!  Programma ottimizzato con BCG + gradT + turbolenza con una
!  sola velocità (250 km/s) e calcolo di rho0
!*****************************************************************
module constants     
    implicit none
    real*8, parameter :: msol = 1.989d33
    real*8, parameter :: cmkpc = 3.084e21
    real*8, parameter :: mu = 0.62d0
    real*8, parameter :: boltz = 1.38066e-16
    real*8, parameter :: guniv = 6.6720e-8
    real*8, parameter :: mp = 1.67265e-24

    real*8, parameter :: rho0nfw = 7.35d-26
    real*8, parameter :: rs = 435.7d0 * cmkpc
    real*8, parameter :: ticm = 8.9d7
    real*8, parameter :: rvir = 2797.d0 * cmkpc
    real*8, parameter :: fc = 1.138799d0
    real*8, parameter :: mvir = 1.3d15 * msol
    real*8, parameter :: ahern = 10.d0 * cmkpc / (1.d0 + sqrt(2.d0))
    real*8, parameter :: r500 = rvir / 2.d0
end module constants

program heq_prova
    use constants
    parameter(jmax=5000)
    implicit real*8 (a-h,o-z)
    integer :: j
    real*8, dimension(jmax) :: r, rr, vol, mnfw, mdark, mhern, grv, gradT, T
    real*8, dimension(jmax) :: rho_nobcg, rho_bcg, rho_thermal, rho_vt
    real*8, dimension(jmax) :: mg_nobcg, mg_bcg, mg_thermal, mg_vt
    real*8 :: fb_nobcg, fb_bcg, fb_thermal, fb_target
    real*8 :: vt_kms, mbgc_nonzero

    vt_kms = 250.0d0
    mbgc_nonzero = 1.d12 * msol

    !************************************************************
    ! rhobest per fb = 0.16??
    !************************************************************
        
    call rhotest(0.0d0, 0.0d0 , 0.0d0 , rho0_nobcg, fb_nobcg)
    call rhotest(mbgc_nonzero, 0.0d0 , 0.0d0 , rho0_bcg, fb_bcg)
    call rhotest(mbgc_nonzero, 1.0d0 , 0.0d0 , rho0_thermal, fb_thermal)
    call rhotest(mbgc_nonzero, 1.0d0 , vt_kms , rho0_vt, fb_vt)

    !****************************************************************************************************
    !end turb rho0test*************
    !****************************************************************************************************

    print *, 'rho0:', rho0_nobcg, rho0_bcg, rho0_thermal, rho0_vt

    !*******************************************************************************
    !DENSITIES
    !*******************************************************************************

    CALL rho_gas(0.0d0, 0.0d0, rho0_nobcg, 0d0, rho_nobcg, mg_nobcg)            !NO BCG
    CALL rho_gas(mbgc_nonzero, 0.0d0, rho0_bcg, 0d0, rho_bcg, mg_bcg)   !BCG+ISOTHERMAL
    CALL rho_gas(mbgc_nonzero, 1.0d0, rho0_thermal, 0d0, rho_thermal, mg_thermal)   !BCG+GRADT
    CALL rho_gas(mbgc_nonzero, 1.0d0, rho0_vt, vt_kms, rho_vt, mg_vt)   !BCG+GRADT+TURBOL

    open(20,file='density_prova.dat',status='unknown')
    do j=1, jmax
        write(20,1000)rr(j)/cmkpc,rho_nobcg(j),rho_bcg(j),rho_thermal(j),rho_vt(j)
    enddo
    close(20)
    1000 format(5(1pe12.4))

    open(20,file='gasmass_prova.dat')
    do j=1,jmax
        write(20,1006)r(j)/cmkpc,mg_nobcg(j)/msol,mg_bcg(j)/msol,mg_thermal(j)/msol, mg_vt(j)/msol
    enddo
    1006 format(5(1pe12.4))
    close(20)

    print*, 'baryon fractions:', fb_nobcg, fb_bcg, fb_thermal, fb_vt

stop
end program heq_prova

!********************************************************************************
! SUBROUTINE MAIN
! Calcola grid, profili di massa, temperatura e gravità
!********************************************************************************

SUBROUTINE main(mbcg, Tprof, r, rr, vol, mnfw, mdark, mhern, T, gradT, grv)
    use constants
    implicit real*8 (a-h,o-z)
    integer :: j
    parameter(jmax=5000)
    real*8, INTENT(IN):: mbcg, Tprof 
    real*8, INTENT(OUT) :: r(jmax), rr(jmax), vol(jmax), mnfw(jmax), mdark(jmax), &
                            mhern(jmax), T(jmax), gradT(jmax), grv(jmax)
    real*8 :: x, y, rhonfw(jmax)

    !------------------
    !!set the grid
    !------------------
    rmin = 0.*cmkpc
    rmax = 3000.*cmkpc
    do j=1,jmax 
        r(j)=rmin+(j-1)*rmax/(jmax-1)
    enddo
    do j=1,jmax-1
        rr(j)=r(j)+0.5*(r(j+1)-r(j))
    enddo
    rr(jmax)=rr(jmax-1)+(rr(jmax-1)-rr(jmax-2))
    open(10,file='grid.dat',status='unknown')
    do j=1,jmax
        write(10,*)real(r(j)/cmkpc),real(rr(j)/cmkpc)
    enddo
    close(10)

    vol(1)=4.1888*r(1)**3
    do j=2,jmax
        vol(j)=4.1888*(r(j)**3-r(j-1)**3)    !! centrato a rr(j-1) !!
    enddo

    !_______________________
    ! NFW and BCG profiles
    !_______________________
    do j=1,jmax
        x=rr(j)/rs
        rhonfw(j)=rho0nfw/(x*(1.+x)**2)
    enddo

    open(20,file='masse_prova.dat')
    mnfw(1)=0.
    do j=2,jmax
        x=r(j)/rs
        mnfw(j)=mnfw(j-1)+rhonfw(j-1)*vol(j)
        mdark(j)=mvir*(log(1.+x)-x/(1.+x))/fc
        mhern(j)=mbgc*r(j)**2/(r(j)+ahern)**2
        write(20,1001)r(j)/cmkpc,mnfw(j)/msol,mdark(j)/msol,mhern(j)/msol
    enddo
    1001 format(4(1pe12.4))
    close(20)

    !.........................
    ! Temperatures
    !.........................
    open(20, file='temperature_prova.dat')
    if (Tprof /= 0.0) then
        do j=1, jmax
            y=rr(j)/r500
            T(j)=ticm*1.35d0*((y/0.045d0)**1.9+0.45d0)/((y/0.045d0)**1.9+1)*(1+(y/0.6d0)**2)**(-0.45)
        end do
    else
        do j=1, jmax
            T(j)=ticm
        end do
    end if

    gradT(1)=0.0d0
    do j=2, jmax
        gradT(j) = log(T(j)) - log(T(j - 1))
    end do

    !::::::::::::
    !gravity
    !::::::::::::
    grv(1)=0.        
    do j=2,jmax
        grv(j)=guniv*(mnfw(j)+mhern(j))/r(j)**2
    enddo

return
end subroutine main


!********************************************************************************
! SUBROUTINE RHOTEST
! Calcola rho0_best per un dato Tprof, vt_kms, e BCG
!********************************************************************************

subroutine rhotest(mbcg, Tprof, vt_kms, rho0_best, fb)
    use constants
    implicit real*8 (a-h,o-z)
    parameter(jmax=5000)
    real*8,  INTENT(IN) :: mbcg, Tprof, vt_kms
    real*8, INTENT(OUT) :: rho0_best

    
    real*8 :: r(jmax), rr(jmax), vol(jmax), mnfw(jmax), mdark(jmax), mhern(jmax)
    real*8 :: T(jmax), gradT(jmax), grv(jmax)
    real*8 :: rho(jmax), lnd(jmax), mg(jmax)
    real*8 :: rho0min, rho0max, rho0test, fb_target, fb, tol
    real*8 :: mtotvir, mgasvir, vturbl
    integer :: i, j, jvir, maxiter, j_min
    real*8 :: dist_min

    fb_target = 0.16d0
    rho0min = 1.d-28
    rho0max = 1.d-24
    tol = 1.d-3
    maxiter = 10000

    vturbl = vt_kms * 1.0d5 !conversione in cm/s

    CALL main(mbcg, Tprof, r, rr, vol, mnfw, mdark, mhern, T, gradT, grv)
    
    dist_min = abs(r(1) - rvir)
    j_min = 1
    do j = 2, jmax
    if (abs(r(j) - rvir) < dist_min) then
        dist_min = abs(r(j) - rvir)
        j_min = j
    endif
    enddo
    jvir = j_min

    do i = 1, maxiter
        rho0test = 0.5d0 * (rho0min + rho0max)

        !density profile test
        grv(1)=0. 
        lnd(1) = log(rho0test)
        do j = 2, jvir
            grv(j)=guniv*(mnfw(j)+mhern(j))/r(j)**2
            lnd(j)=lnd(j-1)-(1.0d0+(vturbl**2/((1.5d4**2)*T(j))))**(-1)*(grv(j)*(mu*mp)*(rr(j)-rr(j-1))/(boltz*T(j))+ gradT(j))
        enddo

        do j = 1, jvir
            rho(j) = exp(lnd(j))
        enddo

        mg(1) = 0.d0
        do j = 2, jvir
            mg(j) = mg(j-1) + rho(j-1)*vol(j)
        enddo

        mtotvir = mnfw(jvir) + mhern(jvir) + mg(jvir)
        mgasvir = mg(jvir)
        fb = mgasvir / mtotvir

        if (abs(fb - fb_target) < tol) then
            print *, 'rho0best = ', rho0test, ' -> f_b =', fb, 'jvir =', jvir
            exit
        endif

        ! new range
        if ((fb - fb_target) > 0.d0) then
            rho0max = rho0test
        else
            rho0min = rho0test
        endif

        if (i == maxiter) then
            print *, 'max iterations'
        endif
    enddo

    rho0_best = rho0test

return
end subroutine rhotest


!********************************************************************************
! SUBROUTINE RHO_GAS
! Calcola densità e massa di gas per un set di parametri
!********************************************************************************
subroutine rho_gas(mbcg, Tprof, rho0, vt_kms, rho, mg)
    use constants
    parameter(jmax=5000)
    real*8,  INTENT(IN) :: mbcg, Tprof, rho0, vt_kms
    real*8, INTENT(OUT) :: rho(jmax), mg(jmax)
    real*8 :: r(jmax), rr(jmax), vol(jmax), mnfw(jmax), mdark(jmax), mhern(jmax)
    real*8 :: T(jmax), gradT(jmax), grv(jmax), lnd(jmax)
    vturbl = vt_kms * 1.0d5 !conversione in cm/s

    CALL main(mbcg, Tprof, r, rr, vol, mnfw, mdark, mhern, T, gradT, grv)
    !CALL rhotest(mbcg, Tprof, vt_kms, rho0_best, fb)
   
    lnd(1)=log(rho0_best)          !! mette il gas in eq. con il potenziale
    do j=2,jmax
        grv(1)=0.
        grv(j)=guniv*(mnfw(j)+mhern(j))/r(j)**2
        lnd(j)=lnd(j-1)-(1.0d0+(vturbl**2/((1.5d4**2)*T(j))))**(-1)*(grv(j)*(mu*mp)*(rr(j)-rr(j-1))/(boltz*T(j))+ gradT(j)) !cosidera il grad di temp mentre calcola la densità del gas
    end do

    !calcolo rho e mgas per diversi valori di vturbl
    mg(1)=0.d0
    do j = 1, jmax
        rho(j) = exp(lnd(j))
    end do
    do j = 2, jmax
        mg(j) = mg(j-1) + rho(j-1)*vol(j)
    end do

return
end subroutine rho_gas
   

