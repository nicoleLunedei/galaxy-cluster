parameter(jmax=5000)
implicit real*8 (a-h,o-z)
real*8, dimension(jmax) :: r(jmax),rr(jmax),vol(jmax),mnfw(jmax),&
        rhost(jmax),rho(jmax),mhern(jmax),rhonfw(jmax),mdark(jmax),&
        grvnfw(jmax),lnd(jmax), mdm_analytic(jmax), T(jmax), lnT(jmax),&
        Mgas(jmax)
real*8 :: msol,mu,mp,rmin,rmax,mvir,rvir,mbgc,ahern, fb
integer, parameter :: nvt = 5
real*8, dimension(nvt) :: vt_values
real*8, dimension(jmax, nvt) :: Mgas_vt, rho_vt
real*8, dimension(nvt) :: fb_vt




!constants

msol = 1.989d33
cmkpc = 3.084e21
mu=0.62
boltz=1.38066e-16
guniv=6.6720e-8
mp=1.67265e-24

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
open(10,file='grid_gradT_turbl.dat',status='unknown')
do j=1,jmax
   write(10,*)real(r(j)/cmkpc),real(rr(j)/cmkpc)
enddo
close(10)

!!volume of shells

vol(1)=4.1888*r(1)**3
do j=2,jmax
   vol(j)=4.1888*(r(j)**3-r(j-1)**3)    !! centrato a rr(j-1) !!
enddo

!  parametri del problema

rho0nfw=7.35d-26
rs=435.7*cmkpc
rho0=8d-26 		      !!central density/initial condition to be modified
ticm=8.9d7		      !!temperature of the ICM

rvir=2797.*cmkpc
r500=rvir/2.0d0
fc=1.138799
mvir=1.3e15*msol
mbgc=1.d12*msol
ahern=10.*cmkpc/(1.+2.**(0.5))

!!DM density profile

do j=1,jmax
   x=rr(j)/rs
   rhonfw(j)=rho0nfw/(x*(1.+x)**2)
enddo

open(20,file='masse_gradT_turbl.dat')
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
enddo
do j=1, jmax
   write(20,1001)r(j)/cmkpc,mnfw(j)/msol,mdark(j)/msol,mhern(j)/msol, mdm_analytic(j)/msol
enddo
1001 format(5(1pe12.4))
close(20)

!! temperature
open(20, file='temperature_gradT_turbl.dat')
do j=1, jmax
    y=rr(j)/r500
    T(j)=ticm*1.35d0*((y/0.045d0)**1.9+0.45d0)/((y/0.045d0)**1.9+1)*(1+(y/0.6d0)**2)**(-0.45) !no more temp const
    lnT(j)=log(T(j))
    write(20,1003)rr(j)/cmkpc, T(j), ticm !ticm is temp const
enddo
1003 format(3(1pe12.4))
close(20)

!!gravitational energy

open(20,file='grv_gradT_turbl.dat')
grvnfw(1)=0.          !! ok per alone NFW, isotermo o beta-model
do j=2,jmax
   grvnfw(j)=guniv*(mnfw(j)+mhern(j))/r(j)**2		!!already takes into account the mass of the central galaxy BCG
   write(20,1002)r(j)/cmkpc,grvnfw(j)/msol
enddo
1002 format(2(1pe12.4))
close(20)

!calculate the gas density
vt_values = (/ 200.0d0, 250.0d0, 300.0d0, 350.0d0, 400.0d0 /) !vt appartien a intervallo 200-400 km/s 
do ivt = 1, nvt !testare i valori di vt 
  vt_kms = vt_values(ivt) !conversione in m/s
  vturbl = vt_kms * 1.0d5 !conversione in m/s
   
  lnd(1)=log(rho0)          !! mette il gas in eq. con il potenziale
  do j=2,jmax
     gg=grvnfw(j)
     lnd(j)=lnd(j-1)-(1.0d0+(vturbl**2/((1.5d4**2)*T(j))))**(-1)*(gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*T(j))+ lnT(j)-lnT(j-1)) !cosidera il grad di temp mentre calcola la densità del gas
  end do

!calcolo rho e mgas per diversi valori di vturbl
  rho_vt(1, ivt) = exp(lnd(1))
  Mgas_vt(1, ivt) = 0.0d0
  do j=2,jmax
     rho_vt(j, ivt)=exp(lnd(j)) 
     Mgas_vt(j, ivt) = Mgas_vt(j-1, ivt) + rho_vt(j-1, ivt) * vol(j)
  end do

  fb_vt(ivt) = (mhern(jmax) + Mgas_vt(jmax, ivt))/(mnfw(jmax) + mhern(jmax) + Mgas_vt(jmax, ivt))!!baryon fraction
end do

open(20,file='density_gradT_turbl.dat',status='unknown')
do j=1,jmax
   write(20,1000)rr(j)/cmkpc,(rho_vt(j, ivt), ivt = 1, nvt),rhonfw(j)   !!restituisce i punti, la densità  numerica e la densità analitica, massa gas, raggio viriale e baryon fraction
end do
close(20)
1000 format(7(1pe12.4))

open(21, file='fb_gradT_turbl.dat', status='unknown')
do ivt = 1, nvt
   write(21,1004) vt_values(ivt), fb_vt(ivt)
enddo
close(21)
1004 format(2(1pe12.4))


stop
end
