!******************************************************************
!  this program solves the hydrostatic equilibrium equation
!  for an isothermal gas in a NFW halo
!*****************************************************************

parameter(jmax=5000)
implicit real*8 (a-h,o-z)
real*8, dimension(jmax) :: r(jmax),rr(jmax),vol(jmax),mnfw(jmax),&
        rhost(jmax), rho(jmax), rhonobcg(jmax), mhern(jmax), rhonfw(jmax), mdark(jmax),&
        grvnfw(jmax), grvnobcg(jmax), lnd(jmax), mg(jmax), mgnobcg(jmax), mgthermal(jmax), rhothermal(jmax),&
        T(jmax), lnT(jmax)
real*8 :: msol,mu,mp,rmin,rmax,mvir,rvir,mbgc,ahern
real*8, dimension(jmax, 4) :: mgturb, rhoturb, lndt
real(8), dimension(4) :: rho0turb, vturb


real*8 :: fb_target, fb, rho0min, rho0max, rho0test, mtotvir, mgasvir
real*8 :: tol
integer :: i, maxiter, jvir
real*8 :: dist_min
integer :: j_min


!  constants

msol = 1.989d33
cmkpc = 3.084e21
mu=0.62
boltz=1.38066e-16
guniv=6.6720e-8
mp=1.67265e-24
vturb = (/ 2.0d7, 3.0d7, 4.0d7, 5.0d7 /) !vt appartiene a intervallo 200-400 km/s già convertito in cm qua 

fb_target = 0.16d0
rho0min = 1.d-28
rho0max = 1.d-24
tol = 1.d-3
maxiter = 1000


!    set the grid

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

!  parametri del problema

rho0nfw=7.35d-26
rs=435.7*cmkpc
!rho0=2.882d-26
ticm=8.9e7

rvir=2797.*cmkpc
fc=1.138799
mvir=1.3e15*msol
mbgc=1.d12*msol
ahern=10.*cmkpc/(1.+2.**(0.5))
r500=rvir/2.0d0

do j=1,jmax
   x=rr(j)/rs
   rhonfw(j)=rho0nfw/(x*(1.+x)**2)
enddo

open(20,file='masse.dat')
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


!! temperature
open(20, file='temperature.dat')
do j=1, jmax
    y=rr(j)/r500
    T(j)=ticm*1.35d0*((y/0.045d0)**1.9+0.45d0)/((y/0.045d0)**1.9+1)*(1+(y/0.6d0)**2)**(-0.45)
    lnT(j)=log(T(j))
    write(20,1003)rr(j)/cmkpc, T(j), ticm
enddo
1003 format(3(1pe12.4))
close(20)

!gravity
open(20,file='grv.dat')
grvnfw(1)=0.          !! ok per alone NFW, isotermo o beta-model
do j=2,jmax
   grvnobcg(j)=guniv*(mnfw(j))/r(j)**2
   grvnfw(j)=guniv*(mnfw(j)+mhern(j))/r(j)**2
   write(20,1002)r(j)/cmkpc,grvnfw(j)/msol,grvnobcg(j)/msol
enddo
1002 format(3(1pe12.4))
close(20)


!************************************************************
! rhobest per fb = 0.16??
!************************************************************

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

rho0 = rho0test

!same for nobcg

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

rho0nobcg = rho0test

!same for thermal

rho0min = 1.d-28
rho0max = 1.d-24
do i = 1, maxiter
   rho0test = 0.5d0 * (rho0min + rho0max)

   !density profile test
   lnd(1) = log(rho0test)
   do j = 2, jvir
      lnd(j)=lnd(j-1)-(grvnfw(j)*(mu*mp)*(rr(j)-rr(j-1))/(boltz*T(j))+ lnT(j)-lnT(j-1))
   enddo

   do j = 1, jvir
      rhothermal(j) = exp(lnd(j))
   enddo
   mg(1)=0.d0
   do j = 2, jvir
      mg(j) = mg(j-1) + rhothermal(j-1)*vol(j)
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

rho0thermal = rho0test

!same for turb: 

do iturb = 1, 4  
   rho0min = 1.d-28 !qua ridefinisco dentro il ciclo perchè quando esce dal do precedente i valori sono ancora quelli dell'iturb precedente
   rho0max = 1.d-24
  

   do i = 1, maxiter
    rho0test = 0.5d0 * (rho0min + rho0max)  !qua riinializziamo per ogni iterazione fissato i turb
      lndt(1,iturb) = log(rho0test)

      do j = 2, jvir
         gg = grvnfw(j)
         lndt(j, iturb) = lndt(j-1, iturb) - (1.0d0 + (vturb(iturb)**2 / ((1.5d4**2) * T(j))))**(-1) * & !espressione di anna con tutti i contributi
         (gg * (mu * mp) * (rr(j) - rr(j-1)) / (boltz * T(j)) + lnT(j) - lnT(j-1))
      end do

      do j = 1, jvir
         rhoturb(j,iturb) = exp(lndt(j,iturb))
      end do

      mgturb(1,iturb) = 0.d0
      do j = 2, jvir
         mgturb(j,iturb) = mgturb(j-1,iturb) + rhoturb(j-1,iturb) * vol(j)
      end do

      mtotvir = mnfw(jvir) + mhern(jvir) + mgturb(jvir,iturb)  
      mgasvir = mgturb(jvir,iturb)  
      fb = mgasvir / mtotvir

      if (abs(fb - fb_target) < tol) then
         print *, 'rho0best for iturb ', iturb, ' = ', rho0test, ' -> f_b =', fb, ' jvir =', jvir
         exit
      endif

      if ((fb - fb_target) > 0.d0) then
         rho0max = rho0test
      else
         rho0min = rho0test
      endif

      if (i == maxiter) then
         print *, 'max iterations for iturb ', iturb
      endif
   enddo

   rho0turb(iturb) = rho0test !4dimensioni, una per ogni iturb (in ordine 200, 300, 400, 500 km/h)
enddo



!end

print *, 'densities0', rho0, rho0nobcg, rho0thermal

!     calculate the gas density, assuming ticm

lnd(1)=log(rho0)          !! mette il gas in eq. con il potenziale
do j=2,jmax
   gg=grvnfw(j)
   lnd(j)=lnd(j-1)-gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*ticm)
enddo

do j=1,jmax
   rho(j)=exp(lnd(j))
enddo

!same for nobcg
lnd(1)=log(rho0nobcg)          !! mette il gas in eq. con il potenziale
do j=2,jmax
   gg=grvnobcg(j)
   lnd(j)=lnd(j-1)-gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*ticm)
enddo

do j=1,jmax
   rhonobcg(j)=exp(lnd(j))
enddo

!!!end

!thermal g density

lnd(1)=log(rho0thermal)          !! mette il gas in eq. con il potenziale
do j=2,jmax
   gg=grvnfw(j)
   lnd(j)=lnd(j-1)-(gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*T(j))+ lnT(j)-lnT(j-1))
enddo

do j=1,jmax
   rhothermal(j)=exp(lnd(j))  
enddo

!end

!turbulent gas density

 do iturb = 1, 4
   lndt(1,iturb) = log(rho0turb(iturb))
 
 do j = 2, jmax
      gg = grvnfw(j)
      lndt(j, iturb) = lndt(j-1, iturb) - (1.0d0 + (vturb(iturb)**2 / ((1.5d4**2) * T(j))))**(-1) * & !espressione di anna con tutti i contributi
      (gg * (mu * mp) * (rr(j) - rr(j-1)) / (boltz * T(j)) + lnT(j) - lnT(j-1))
   end do

   do j = 1, jmax
      rhoturb(j,iturb) = exp(lndt(j,iturb))
   end do
end do
!end


open(20,file='density.dat',status='unknown')
do j=1,jmax
   write(20,1000)rr(j)/cmkpc,rho(j),rhonfw(j),rhonobcg(j),rhothermal(j), &
                   (rhoturb(j,iturb), iturb=1,4)
enddo
close(20)
1000 format(9(1pe12.4))

open(20,file='gasmass.dat')
mg(1)=0
mgthermal(1)=0
mgnobcg(1)=0
do iturb = 1, 4
   mgturb(1, iturb)=0
end do
do j=2,jmax
   mg(j)=mg(j-1)+rho(j-1)*vol(j)
   mgnobcg(j)=mgnobcg(j-1)+rhonobcg(j-1)*vol(j)
   mgthermal(j)=mgthermal(j-1)+rhothermal(j-1)*vol(j)
    do iturb = 1, 4
       mgturb(j,iturb) = mgturb(j-1,iturb) + rhoturb(j-1,iturb)*vol(j)
   end do

   write(20,1006) r(j)/cmkpc, mg(j)/msol, mgnobcg(j)/msol, mgthermal(j)/msol, &
                   (mgturb(j,iturb)/msol, iturb=1,4)
end do

1006 format(8(1pe12.4))
close(20)



stop
end

