!*****************************************************************************
!  Programma con calcolo di rhobest e casi diversi di velocità di turbolenza
!*****************************************************************_***********

parameter(jmax=5000)
implicit real*8 (a-h,o-z)
real*8, dimension(jmax) :: r(jmax),rr(jmax),vol(jmax),mnfw(jmax),&
        rhost(jmax), rho(jmax), rhonobcg(jmax), mhern(jmax), rhonfw(jmax), mdark(jmax),&
        grvnfw(jmax), grvnobcg(jmax), lnd(jmax), mg(jmax), mgnobcg(jmax), mgthermal(jmax), rhothermal(jmax), T(jmax), lnT(jmax)
real*8 :: msol,mu,mp,rmin,rmax,mvir,rvir,mbgc,ahern
!test parameters
real*8 :: fb_target, fb, rho0min, rho0max, rho0test, mtotvir, mgasvir, tol, dist_min
integer :: i, maxiter, jvir, j_min
!turb parameters
integer, parameter :: nvt = 5
real*8, dimension(nvt) :: vt_values = (/ 200.0d0, 250.0d0, 300.0d0, 400.0d0, 500.0d0 /)  
real*8, dimension(jmax, nvt) :: Mgas_vt, rho_vt
real*8, dimension(nvt) :: fb_vt
real*8, dimension(jmax, nvt) :: mg_vt

!Constants
msol = 1.989d33
cmkpc = 3.084e21
mu=0.62
boltz=1.38066e-16
guniv=6.6720e-8
mp=1.67265e-24

!Test parameters
fb_target = 0.16d0
rho0min = 1.d-28
rho0max = 1.d-24
tol = 1.d-3
maxiter = 10000

!Problem's parameters
rho0nfw=7.35d-26
rs=435.7*cmkpc
ticm=8.9e7
rvir=2797.*cmkpc
r500=rvir/2.0d0
fc=1.138799
mvir=1.3e15*msol
mbgc=1.d12*msol
ahern=12.*cmkpc/(1.+2.**(0.5))


!Set the grid
rmin = 0.*cmkpc
rmax = 3000.*cmkpc
do j=1,jmax
   r(j)=rmin+(j-1)*rmax/(jmax-1)
end do
do j=1,jmax-1
   rr(j)=r(j)+0.5*(r(j+1)-r(j))
end do
rr(jmax)=rr(jmax-1)+(rr(jmax-1)-rr(jmax-2))

!Volume
vol(1)=4.1888*r(1)**3
do j=2,jmax
   vol(j)=4.1888*(r(j)**3-r(j-1)**3)    !! centrato a rr(j-1) !!
end do

do j=1,jmax
   x=rr(j)/rs
   rhonfw(j)=rho0nfw/(x*(1.+x)**2)
end do

mnfw(1)=0.
do j=2,jmax
   x=r(j)/rs
   mnfw(j)=mnfw(j-1)+rhonfw(j-1)*vol(j)
   mdark(j)=mvir*(log(1.+x)-x/(1.+x))/fc
   mhern(j)=mbgc*r(j)**2/(r(j)+ahern)**2
end do

!!Temperature
do j=1, jmax
    y=rr(j)/r500
    T(j)=ticm*1.35d0*((y/0.045d0)**1.9+0.45d0)/((y/0.045d0)**1.9+1)*(1+(y/0.6d0)**2)**(-0.45)
    lnT(j)=log(T(j))
end do

!Gravitational force
grvnfw(1)=0.         
do j=2,jmax
   grvnobcg(j)=guniv*(mnfw(j))/r(j)**2
   grvnfw(j)=guniv*(mnfw(j)+mhern(j))/r(j)**2
end do

!****************************************************************************************************
! TEST FOR RHOBEST so that fb = 0.16
!****************************************************************************************************

dist_min = abs(r(1) - rvir)
j_min = 1
do j = 2, jmax
   if (abs(r(j) - rvir) < dist_min) then
      dist_min = abs(r(j) - rvir)
      j_min = j
   endif
enddo
jvir = j_min

!!!!!!!!!!!!!!!!!!!!!BCG **************************************************
do i = 1, maxiter
   rho0test = 0.5d0 * (rho0min + rho0max)

   !density profile test
   lnd(1) = log(rho0test)
   do j = 2, jvir
      lnd(j) = lnd(j-1) - grvnfw(j)*mu*mp*(rr(j)-rr(j-1))/(boltz*ticm)
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
      print *, 'BCG only ', 'rho0best = ', rho0test, ' -> f_b =', fb
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

rho0 = rho0test

!!!!!!!!!!!!!!!!!!!!!!!!!! NO BCG ***********************************************
rho0min = 1.d-28
rho0max = 1.d-24
do i = 1, maxiter
   rho0test = 0.5d0 * (rho0min + rho0max)

   !density profile test
   lnd(1) = log(rho0test)
   do j = 2, jvir
      lnd(j) = lnd(j-1) - grvnobcg(j)*mu*mp*(rr(j)-rr(j-1))/(boltz*ticm)
   enddo

   do j = 1, jvir
      rhonobcg(j) = exp(lnd(j))
   enddo
   mg(1)=0.d0
   do j = 2, jvir
      mg(j) = mg(j-1) + rhonobcg(j-1)*vol(j)
   enddo

   mtotvir = mnfw(jvir) + mg(jvir)
   mgasvir = mg(jvir)
   fb = mgasvir / mtotvir

if (abs(fb - fb_target) < tol) then
      print *, 'No BCG ', 'rho0best = ', rho0test, ' -> f_b =', fb
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

rho0nobcg = rho0test

!!!!!!!!!!!!!!!!!!!!!!!!!!! BCG + gradT *******************************************************+
rho0min = 1.d-28
rho0max = 1.d-24
do i = 1, maxiter
   rho0test = 0.5d0 * (rho0min + rho0max)

   !density profile test
   lnd(1) = log(rho0test)
   do j = 2, jmax
      lnd(j)=lnd(j-1)-(grvnfw(j)*(mu*mp)*(rr(j)-rr(j-1))/(boltz*T(j))+ lnT(j)-lnT(j-1))
   enddo

   do j = 1, jmax
      rhothermal(j) = exp(lnd(j))
   enddo
   mg(1)=0.d0
   do j = 2, jmax
      mg(j) = mg(j-1) + rhothermal(j-1)*vol(j)
   enddo

   mtotvir = mnfw(jvir) + mhern(jvir) + mg(jvir)
   mgasvir = mg(jvir)
   fb = mgasvir / mtotvir

if (abs(fb - fb_target) < tol) then
      print *, 'BCG+gradT ', 'rho0best = ', rho0test, ' -> f_b =', fb
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

rho0thermal = rho0test

!!!!!!!!!!!!!!!!!!!!!! BCG + gradT + TURBOLENCE *********************************************
rho0min = 1.d-28
rho0max = 1.d-24
do i = 1, maxiter
   rho0test = 0.5d0 * (rho0min + rho0max)

   !density profile test
   do ivt = 1, nvt 
      vt_kms = vt_values(ivt)
      vturbl = vt_kms * 1.0d5 
   
      lnd(1)=log(rho0test)     
      do j=2,jmax
         gg=grvnfw(j)
         lnd(j)=lnd(j-1)-(1.0d0+(vturbl**2/((1.5d4**2)*T(j))))**(-1)*(gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*T(j))+ lnT(j)-lnT(j-1)) 
      end do

     !calcolo rho e mgas per diversi valori di vturbl
      mg_vt(1, ivt)=0.d0
      do j = 1, jmax
         rho_vt(j, ivt) = exp(lnd(j))
      end do
      do j = 2, jmax
         mg_vt(j, ivt) = mg_vt(j-1, ivt) + rho_vt(j-1, ivt)*vol(j)
      end do

      mtotvir = mnfw(jvir) + mhern(jvir) + mg_vt(jvir, ivt)
      mgasvir = mg_vt(jvir, ivt)
      fb = mgasvir/ mtotvir
      
      if (abs(fb - fb_target) < tol) then
         print *, 'Turbolence ', vt_values(ivt), 'rho0best = ', rho0test, ' -> f_b =', fb
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
   end do
end do

rho0_vt = rho0test


!*****************************
!end turb rho0test*************
!*****************************

!print *, 'densities0', rho0, rho0nobcg, rho0thermal, rho0_vt

stop
end


