!******************************************************************
! Codice per il caso BCG + T gradient + vari valori di turbolenza
!*****************************************************************

parameter(jmax=5000)
implicit real*8 (a-h,o-z)
real*8:: r(jmax),rr(jmax),vol(jmax),mnfw(jmax),&
        rhost(jmax),mhern(jmax),rhonfw(jmax),mdark(jmax),&
        grvnfw(jmax),lnd(jmax), mdm_analytic(jmax), T(jmax), lnT(jmax),&
        Mgas(jmax)
real*8, dimension(5, jmax) :: rho
real*8 :: msol,mu,mp,rmin,rmax,mvir,rvir,mbgc,ahern, T_j
real*8, dimension(5) :: vt_kms = (/ 200.0d0, 250.0d0, 300.0d0, 400.0d0, 500.0d0 /)
real*8, dimension(5) ::rho0= (/1.76d-25, 1.602d-25, 1.47d-25, 1.15d-25, 8.80d-26/)

!constants
msol = 1.989d33
cmkpc = 3.084e21
mu=0.62
boltz=1.38066e-16
guniv=6.6720e-8
mp=1.67265e-24

! parametri del problema
rho0nfw=7.35d-26
rs=435.7*cmkpc 		      !!central density/initial condition to be modified
ticm=8.9d7		      !!temperature of the ICM
rvir=2797.*cmkpc
r500=rvir/2.0d0
fc=1.138799
mvir=1.3e15*msol
mbgc=1.d12*msol
ahern=12.*cmkpc/(1.+2.**(0.5))

do i=1, 5
   vturbl = vt_kms(i) * 1.0d5

   !set the grid
   rmin = 0.*cmkpc
   rmax = 3000.*cmkpc
   do j=1,jmax
      r(j)=rmin+(j-1)*rmax/(jmax-1)    !!r_j
   end do
   do j=1,jmax-1
      rr(j)=r(j)+0.5*(r(j+1)-r(j))     !!r_j+1/2
   end do
   rr(jmax)=rr(jmax-1)+(rr(jmax-1)-rr(jmax-2))

   !!volume of shells
   vol(1)=4.1888*r(1)**3
   do j=2,jmax
      vol(j)=4.1888*(r(j)**3-r(j-1)**3)    !! centrato a rr(j-1) !!
   end do

   !!DM density profile
   do j=1,jmax
      x=rr(j)/rs
      rhonfw(j)=rho0nfw/(x*(1.+x)**2)
   end do

   !!Masses
   !open(20,file='masse.dat')
   mnfw(1)=0.		!!initial dark matter mass
   mdark(1)=0.
   mhern(1)=0.
   mdm_analytic(1)=0.
   mdm_analytic(1)=0.
   do j=2,jmax 		!!NB: starts at j=2 not j=1
      x=r(j)/rs
      mnfw(j)=mnfw(j-1)+rhonfw(j-1)*vol(j)					                     !!DM mass discrete NFW formula in pdf centered in j-1/2
      mdm_analytic(j)=4*3.14*rho0nfw*(rs**3)*(log(1.+x)-(r(j)/(r(j)+rs)))     !!DM analytical
      mdark(j)=mvir*(log(1.+x)-x/(1.+x))/fc				                        !!DM analytic NFW formula
      mhern(j)=mbgc*r(j)**2/(r(j)+ahern)**2				                        !!analytic stellar mass
   end do
   !do j=1, jmax
      !write(20,1001)r(j)/cmkpc,mnfw(j)/msol,mdark(j)/msol,mhern(j)/msol, mdm_analytic(j)/msol
   !enddo
   !1001 format(5(1pe12.4))
   !close(20)

   do j=1, jmax
      y=rr(j)/r500
      T(j)=ticm*1.35d0*((y/0.045d0)**1.9+0.45d0)/((y/0.045d0)**1.9+1)*(1+(y/0.6d0)**2)**(-0.45)
      lnT(j)=log(T(j))
   end do

   !!Gravitational force
   grvnfw(1)=0.          
   do j=2,jmax
      grvnfw(j)=guniv*(mnfw(j)+mhern(j))/r(j)**2	

   !Calculate the gas density
   lnd(1)=log(rho0(i))          !! mette il gas in eq. con il potenziale
   do j=2,jmax
      gg=grvnfw(j)
      T_j = 0.5 * (T(j) + T(j-1))
      lnd(j)=lnd(j-1)-(1.0d0 + 0.5555556*(vturbl**2)/((1.5d4**2)*T_j))**(-1)&
            *(gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*T_j) + lnT(j) - lnT(j-1))
   end do 

   do j=1,jmax
      rho(i, j)=exp(lnd(j))  
   end do
end do

open(20,file='density_turb_var.dat',status='unknown')
do j=1,jmax
   write(20,1000)rr(j)/cmkpc, rho(1,j),rho(2,j),rho(3,j),rho(4,j),rho(5,j), rhonfw(j)     
close(20)
1000 format(7(1pe12.4))

stop
end

