!*******************************************************************************
! Codice finale con tutti i casi: NO BCG, BCG, BCG + gradT, BCG + gradT + turb
!*******************************************************************************

parameter(jmax=5000)
implicit real*8 (a-h,o-z)
real*8 :: r(jmax),rr(jmax),vol(jmax),mnfw(jmax),&
        rhost(jmax),rho(jmax),mhern(jmax),rhonfw(jmax),mdark(jmax),&
        grvnfw(jmax),lnd(jmax), mdm_analytic(jmax), T(jmax), T_mean(jmax), lnT(jmax),&
        Mgas(jmax)
real*8 :: msol,mu,mp, gamma,rmin,rmax,mvir,rvir,mbgc,ahern, BCG, gradT, turbol, vturbl

!constants
msol = 1.989d33
cmkpc = 3.084e21
mu=0.62
boltz=1.38066e-16
guniv=6.6720e-8
mp=1.67265e-24
gamma= 1.666667

!  parametri del problema
rho0nfw=7.35d-26
rs=435.7*cmkpc
ticm=8.9d7		      !!temperature of the ICM
rvir=2797.*cmkpc
r500=rvir/2.0d0
fc=1.138799
mvir=1.3e15*msol
mbgc=1.d12*msol
ahern=12.*cmkpc/(1.+2.**(0.5))

!!Deciding the case
print*, 'BCG? No=0, Yes=1'
read(*,*) BCG
print*, 'gradT? No=0, Yes=1'
read(*,*) gradT
print*, 'Turbolence? No=1, Yes=1'
read(*,*) turbol

!set the grid
rmin = 0.*cmkpc
rmax = 3000.*cmkpc
do j=1,jmax
   r(j)=rmin+(j-1)*rmax/(jmax-1)    !!r_j
enddo
do j=1,jmax-1
   rr(j)=r(j)+0.5*(r(j+1)-r(j))     !!r_j+1/2
enddo
rr(jmax)=rr(jmax-1)+(rr(jmax-1)-rr(jmax-2))

!!volume of shells
vol(1)=4.1888*r(1)**3
do j=2,jmax
   vol(j)=4.1888*(r(j)**3-r(j-1)**3)    !! centrato a rr(j-1) !!
enddo

!!DM density profile
do j=1,jmax
   x=rr(j)/rs
   rhonfw(j)=rho0nfw/(x*(1.+x)**2)
enddo

open(20,file='masse.dat')
mnfw(1)=0.		!!initial dark matter mass
mdark(1)=0.
mhern(1)=0.
mdm_analytic(1)=0.
mdm_analytic(1)=0.
do j=2,jmax
   x=r(j)/rs
   mnfw(j)=mnfw(j-1)+rhonfw(j-1)*vol(j)					                     !!DM mass discrete NFW formula in pdf centered in j-1/2
   mdm_analytic(j)=4*3.14*rho0nfw*(rs**3)*(log(1.+x)-(r(j)/(r(j)+rs)))     !!DM analytical
   mdark(j)=mvir*(log(1.+x)-x/(1.+x))/fc				                        !!DM analytic NFW formula
   !!analytic stellar mass
   if (BCG == 0) then      
      mhern(j)=0.                             !!no BCG case
   else if (BCG == 1) then 
      mhern(j)=mbgc*r(j)**2/(r(j)+ahern)**2  !!with BCG case
   end if				                        
enddo
do j=1, jmax
   write(20,1001)r(j)/cmkpc,mnfw(j)/msol,mdark(j)/msol,mhern(j)/msol, mdm_analytic(j)/msol
enddo
1001 format(5(1pe12.4))
close(20)

!! temperature
open(20, file='temperature_graT.dat')
do j=1, jmax
   y=rr(j)/r500
   if (gradT == 0) then
      T(j)=ticm
      lnT(j)=0.
   else if (gradT == 1) then
      T(j)=ticm*1.35d0*((y/0.045d0)**1.9+0.45d0)/((y/0.045d0)**1.9+1)*(1+(y/0.6d0)**2)**(-0.45)
      lnT(j)=log(T(j))
      write(20,1003)rr(j)/cmkpc, T(j), ticm
   end if
enddo
1003 format(3(1pe12.4))
close(20)

!!gravitational force
grvnfw(1)=0.          !! ok per alone NFW, isotermo o beta-model
do j=2,jmax
   grvnfw(j)=guniv*(mnfw(j)+mhern(j))/r(j)**2		!!already takes into account the mass of the central galaxy BCG
enddo

if (BCG == 0) then                     
   rho0=4.453d-26          !!no BCG + isothermal
   vturbl=0.
else if (BCG == 1 .and. gradT == 0) then
   rho0=1.007d-25          !!BCG + isothermal
   vturbl=0.
else if (BCG == 1 .and. gradT == 1 .and. turbol == 0) then
   rho0=2.032d-25          !!BCG + thermal grad
   vturbl=0.
else if (BCG == 1 .and. gradT == 1 .and. turbol == 1) then
   vturbl= 250.0d5         !!BCG + gradT + turbolence
   rho0=1.602d-25
end if

!calculate the gas density
lnd(1)=log(rho0)     
do j=2,jmax
   gg=grvnfw(j)
   T_j = 0.5 * (T(j) + T(j-1))
   lnd(j)=lnd(j-1)-(1.0d0 + 0.5555556*(vturbl**2)/((1.5d4**2)*T_j))**(-1)&
           *(gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*T_j) + lnT(j) - lnT(j-1))
enddo
do j=1,jmax
   rho(j)=exp(lnd(j))  
enddo

if (BCG == 0) then                     
   open(10,file='density_noBCG.dat',status='unknown')          !!no BCG + isothermal
else if (BCG == 1 .and. gradT == 0) then
   open(10,file='density_BCG.dat',status='unknown')          !!BCG + isothermal
else if (BCG == 1 .and. gradT == 1 .and. turbol == 0) then
   open(10,file='density_BCG_gradT.dat',status='unknown')          !!BCG + thermal grad
else if (BCG == 1 .and. gradT == 1 .and. turbol == 1) then
   open(10,file='density_BCG_gradT_turb.dat',status='unknown')          !!BCG + thermal grad + turbolence
end if
do j=1,jmax
   write(10,1000)rr(j)/cmkpc,rhonfw(j),rho(j)     !!restituisce i punti, la densità  numerica e la densità analitica
end do
close(10)
1000 format(3(1pe12.4))


stop
end
