!*****************************************************************2
!  Program for project 1- Calcolo 20121. 
!  Solves the hydrostatic equilibrium equation in a NFW+BCG potential
!  Then solves the Fe diffusion eq. (for Perseus)
!*****************************************************************

parameter(jmax=5000)
implicit real*8 (a-h,o-z)
real*8, dimension(jmax) :: r(jmax),rr(jmax),vol(jmax),mnfw(jmax),&
        rhost(jmax),rho(jmax),mhern(jmax),rhonfw(jmax),mdark(jmax),&
        grvnfw(jmax),lnd(jmax),tem(jmax),tem2(jmax),mgas(jmax),&
        fbarr(jmax),fgasr(jmax)
real*8 :: msol,mu,mp,rmin,rmax,mvir,rvir,mbcg,ahern,lsol,h,me,&
          ne0,alphast,alphasn,zfesn

real*8, dimension(jmax) :: u(jmax),flux(jmax),ne(jmax),zfe(jmax),& 
        kappa(jmax),lturb,rhofedot(jmax),rhofe(jmax),zfest(jmax),&
        amfeiniz(jmax),amfe(jmax),gradzfe(jmax),zfeobs(jmax),&
        amfeobs(jmax),rhofeobs(jmax)

!  constants

pi=3.14159265359
msol = 1.989d33
cmkpc = 3.084e21
years=3.156d7
mu=0.61
boltz=1.38066e-16
guniv=6.6720e-8
mp=1.67265e-24
zfesol=1.8e-3
zfesn=0.744/1.4
snu=0.5

tnow=13.7*1.e9*years
time0=tnow-5.*1.e9*years
time=time0

!    set the grid

rmin = 0*cmkpc
rmax = 2800.*cmkpc
do j=1,jmax
   r(j)=rmin+(j-1)*rmax/(jmax-1)
enddo
do j=1,jmax-1
   rr(j)=r(j)+0.5*(r(j+1)-r(j))
enddo
rr(jmax)=rr(jmax-1)+(rr(jmax-1)-rr(jmax-2))
open(10,file='grid.dat',status='unknown')
do j=1,jmax
   write(10,*)j,real(r(j)/cmkpc),real(rr(j)/cmkpc)
enddo
close(10)

vol(1)=4.1888*r(1)**3
do j=2,jmax
   vol(j)=4.1888*(r(j)**3-r(j-1)**3)    !! centered in rr(j-1) !!
enddo

!  problem parameters

rho0nfw=7.35d-26
rs=435.7*cmkpc
!!rho0=4.d-26   !! 5000 points, with BCG, isothermal gas !!old: 2.882d-26
!!rho0=9.d-26   !! 5000 points, with BCG !!
rho0=1.6d-25   !! 5000 points, with BCG and dT/dr !!
ticm=8.9e7

rvir=2797.*cmkpc
r500=rvir/2.
fc=1.138799
mvir=1.3e15*msol
mbcg=1.d12*msol
ahern=12.*cmkpc/(1.+2.**(0.5))
aml=7.5    !! this is the mass to light ratio 

do j=1,jmax
   x=rr(j)/rs
   rhonfw(j)=rho0nfw/(x*(1.+x)**2)
   rhost(j)=mbcg/(2.*pi)*(ahern/rr(j))/(ahern+rr(j))**3 !densta analitica
enddo

open(20,file='masse.dat')
mnfw(1)=0.
mhern(1)=0.
do j=2,jmax
   x=r(j)/rs
   mnfw(j)=mnfw(j-1)+rhonfw(j-1)*vol(j)
   mdark(j)=mvir*(log(1.+x)-x/(1.+x))/fc
   mhern(j)=mbcg*r(j)**2/(r(j)+ahern)**2
   write(20,1101)r(j)/cmkpc,mnfw(j)/msol,mdark(j)/msol,mhern(j)/msol
enddo
1101 format(4(1pe12.4))
close(20)

open(20,file='grv.dat')
grvnfw(1)=0.
do j=2,jmax
   grvnfw(j)=guniv*(mnfw(j)+ mhern(j))/r(j)**2   !! with BCG !!
   write(20,1002)r(j)/cmkpc,grvnfw(j)/msol
enddo
1002 format(2(1pe12.4))
close(20)
!
!inizia parte 2 del projet1
! Temperature profile 
!

!identifica costanti
  temp0=8.12e7 
  rtemp1=71.
  rtemp2=71.
  r0bill=40.
  t0bill=0.9 
  tmaxbill=4.8
  tcoeffbill=tmaxbill-t0bill
  pbill=1.8
  qbill=0.2  !!0.15
  sbill=1.6
  qplusp = qbill + pbill
  sinv=1./sbill
!
!profilo temperatura calcolato ma in che modo??
open(20,file='temperature.dat',status='unknown')
 do j=1,jmax
    rkpc=rr(j)/cmkpc
    roverr0 = rkpc/r0bill
    temp1 = t0bill + tcoeffbill * roverr0**pbill
    cut=0.0 !cutoff?
    temp2 = tmaxbill/(cut*(rkpc/60.)**2 + roverr0**qbill) !profilo di temperatura 2
    ttt = 1./( (1./temp1)**sbill + (1./temp2)**sbill )**sinv !ttt intervallo?
    tcorr = -0.12*exp(-rkpc/1.5) !t corente
    tem(j) = (ttt + tcorr)*1.e7    !! this is for NGC 5044 !!
    x=rr(j)/r500
    xx=x/0.045
    tem2(j)=ticm*1.35*(xx**1.9+0.45)/(xx**1.9+1.)* &   !! this is for Perseus !!
        1./(1.+(x/0.6)**2)**0.45
       tem(j)=tem2(j)  !!ticm --> for Perseus
    write(20,1003)rr(j)/cmkpc,tem(j),tem2(j),ticm
 enddo
close(20)
1003 format(4(1pe12.4))

!     calculate the gas density, assuming ticm

lnd(1)=log(rho0)          !! mette il gas in eq. con il potenziale
do j=2,jmax
   gg=grvnfw(j)
   temmed=0.5*(tem(j)+tem(j-1)) !temperatura media
!! isoth !!   lnd(j)=lnd(j-1)-gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*ticm)
   lnd(j)=lnd(j-1)-gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*temmed) &
          - (log(tem(j)) - log(tem(j-1))) 
enddo

do j=1,jmax
   rho(j)=exp(lnd(j))
enddo

open(20,file='density.dat',status='unknown')
do j=1,jmax
   write(20,1101)rr(j)/cmkpc,rho(j),rhonfw(j),rhost(j)
enddo
close(20)

open(20,file='mgas.dat',status='unknown')
mgas(1)=rho(j)*4.188*r(1)**3
do j=2,jmax
   mgas(j)=mgas(j-1)+rho(j-1)*vol(j)
   write(20,1100)r(j)/cmkpc,mgas(j)/msol,mnfw(j)/msol
enddo
close(20)

!!fbar=mgas(jmax-1)/(mnfw(jmax-1)+mgas(jmax-1))
fbar=(mhern(jmax-1)+mgas(jmax-1))/(mnfw(jmax-1)+mgas(jmax-1)+mhern(jmax-1))
fgas=(mgas(jmax-1))/(mnfw(jmax-1)+mgas(jmax-1)+mhern(jmax-1))
print*,'fbari, fgas = ',real(fbar), real(fgas)

open(20,file='barfrac.dat')
do j=2,jmax-1
   fbarr(j)=(mhern(j)+mgas(j))/(mnfw(j)+mgas(j)+mhern(j))
   fgasr(j)=(mgas(j))/(mnfw(j)+mgas(j)+mhern(j))
   write(20,1100)r(j)/cmkpc,fbarr(j),fgasr(j)
enddo
close(20)
1100 format(4(1pe12.4))

!***********************************************************************
!! At this point we have the gas density profile and we can proceed
!! with the integration of the diffusion equation for rhofe
!***********************************************************************

!! Set the initial abundance profile

zfeout=0.4*zfesol   !! this is the background abundance -->viene considerato solo ferro

do j=1,jmax
   x=rr(j)/(80.*cmkpc)
   zfeobs(j)=zfesol*0.3*1.4*1.15*(2.2+x**3)/(1+x**3)/1.15  !Perseus!
   zfeobs(j)=zfeobs(j) - zfeout   !! subtract z_Fe,out-->tolto il background
   zfeobs(j)=max(zfeobs(j),0.) !calcola il max-->per successivo controllo
   zfe(j)=0. !!zfeout !!zfeobs(j)  !! which initial zfe? !!
   rhofe(j)=rho(j)*zfe(j)/1.4 
   rhofeobs(j)=rho(j)*zfeobs(j)/1.4
enddo

 do j=1,jmax
    zfest(j)=1.*zfesol    !! ferro nella stella
 enddo

!! Calculate the initial excess of iron mass

amfeiniz(1)=rhofe(1)*vol(1) !massa del ferro iniziale teorica
amfeiniz(1)=rhofeobs(1)*vol(1) !massa del ferro iniziale numerica
!calcolo massa in modo discreto
do j=2,jmax
   amfeiniz(j)=amfeiniz(j-1)+rhofe(j-1)*vol(j)
   amfeobs(j)=amfeobs(j-1)+rhofeobs(j-1)*vol(j)
enddo

open(20,file='zfe_initial.dat')
do j=1,jmax
   write(20,1500)rr(j)/cmkpc,zfe(j)/zfesol,zfeobs(j)/zfesol, &
                 r(j)/cmkpc,amfeiniz(j)/msol,amfeobs(j)/msol
enddo
close(20)
1500 format(6(1pe12.4))

open(20,file='initial.dat',status='unknown')
do j=1,jmax
   write(20,3001)rr(j)/cmkpc+0.001,zfe(j)/zfesol,zfeobs(j)/zfesol,ne(j)
enddo
close(20)
3001  format(4(1pe12.4))

!! boundary conditions (outflows)

zfe(1)=zfe(2) !per non avere discontinutà agli estremi-->settare i punti corretti per successivi calcoli
zfe(jmax)=zfe(jmax-1)
rhofe(1)=rho(1)*zfe(1)/1.4
rhofe(jmax)=rho(jmax)*zfe(jmax)/1.4

!! Here start the time integration (use FTCS method) !ftcs unico metodo stabile per il problema

 print*,'dai ncycle'
 read(*,*)ncycle
 tend=tnow !per integrare fino al tempo corrente (oggi)

!!  set the diffusion coefficient kappa = C*v*l (for now constant)

 open(20,file='kappa.dat',status='unknown')! kappa è d cioè coeff. di diffusione
 vturb=260.e5   !! come Perseus !!
 lturb=15.*cmkpc  !! this is quite uncertain !!
 rscala=30.*cmkpc !raggio scalar

 kappa=0.11*vturb*lturb !coeff. diffusione

!! do j=1,jmax    !! for variable kappa
!!!    kappa(j)=0.11*vturb*lturb   !! constant !!
!!    kappa(j)=rhost(j)  !!0.333*vturb*lturb   !! constant !!
!!    kappa(j)=kappa(j)-0.6*kappa(j)*exp(-(r(j)/rscala)**2)
!!    write(20,*)real(r(j)/cmkpc),kappa(j)
!! enddo
 close(20)

      n=0
1000  continue      !! here start the main time cycle
      n=n+1

!! calculate the timestep (to be modified if the grid is non-uniform)

      dt=0.4*(r(5)-r(4))**2/(2.*kappa(5))  !! ok for Delta_r costant !!
      time=time+dt
!!    print*,'n,dt (yr),time (Gyr) = ',n,real(dt/years), &
!!            real(time/1.e9/years)

!! write the source terms (SNIa and stellar winds)

 slope=1.1
 alphast=4.7e-20*(time/tnow)**(-1.26)
 alphasn=4.436e-20*(snu/aml)*(time/tnow)**(-slope)

!!  print*,'alphast,sn = ',alphast,alphasn

 do j=2,jmax-2
    rhofedot(j)=(alphast*zfest(j)/1.4+alphasn*zfesn)*rhost(j) !rhost= rho stellar
 enddo

!! the equation to be solved is d(n*zfe)/dt = div(kappa*n*grad(zfe)) + S 
!! (according to Rebusco et al. 2006)
!! Use the FTCS scheme.

!!      goto776
!! source step

 do j=2,jmax-1
!!!    write(70,*)rhofe(j),dt*rhofedot(j)
!!!    if(j.eq.5)print*,'azz ',dt,rhofe(j),dt*rhofedot(j),rhofedot(j)
    rhofe(j)=rhofe(j) + dt*rhofedot(j) !cambiamento di rhofe nel tempo (dt)
    zfe(j)=rhofe(j)/rho(j) * 1.4
 enddo

!! set the boundary conditions (outflows) -->sono cambiati i valori rispetto a prima

      zfe(1)=zfe(2) 
      zfe(jmax)=zfe(jmax-1)
      rhofe(1)=rhofe(2)
      rhofe(jmax)=rhofe(jmax-1)
776   continue

!!  goto777
!  diffusive step   !  check the Fe conservation !

 do j=2,jmax-1
    gradzfe(j)=(zfe(j)-zfe(j-1))/(rr(j)-rr(j-1))  !! dZ/dr centered at "j" !!
 enddo
 gradzfe(1)=0.        !! B.C. !!
 gradzfe(jmax)=0.

 do j=2,jmax-1
    rhojp1=0.5*(rho(j+1)+rho(j))  !! vogliamo ricavre rho ricavata con la prima griglia. per farlo rho si calcola come punto medio di intervallo griglia 1--> rho centered at 'j+1' !!
    rhoj=0.5*(rho(j-1)+rho(j))    !! stesso ragionamento-->rho centered at j-1/2 !!
    !non è j del loop è il j della griglia 1
    rhofe(j)=rhofe(j) &
            + (dt/1.4)*(r(j+1)**2*kappa(j+1)*rhojp1*gradzfe(j+1) &
            -r(j)**2*kappa(j)*rhoj*gradzfe(j))   &
             / (0.33333333*(r(j+1)**3-r(j)**3))
         zfe(j)=1.4*rhofe(j)/rho(j)  !! update Z_Fe with the new rho_Fe !! !! fomula a pag 28 di project1 pdf
      enddo
2000  format(3(1pe12.4))

!! set the boundary conditions (outflows)

      zfe(1)=zfe(2)
      zfe(jmax)=zfe(jmax-1)
      rhofe(1)=rhofe(2)
      rhofe(jmax)=rhofe(jmax-1)
777   continue

      if (time.ge.tend) goto1001 !ge= Greater than or Equal
      if (n.ge.ncycle) goto1001 

      goto 1000

1001  continue

      do j=2,jmax
         write(99,*)real(log10(r(j)/cmkpc)),real(log10(rhofedot(j)))
      enddo

!! calcola la massa di Fe al tempo finale

      amfe(1)=rhofe(1)*vol(1) !massa ferro
      do j=2,jmax
         amfe(j)=amfe(j-1)+rhofe(j-1)*vol(j)
      enddo

      write(6,3002)amfe(jmax)/msol,amfeiniz(jmax)/msol,amfeobs(jmax)/msol
      write(6,3003)amfe(180)/msol,amfeiniz(180)/msol,amfeobs(180)/msol
3002  format('M_Fe(tot), M_Fein(tot) (Msol) = ',3(1pe12.4))
3003  format('M_Fe(100kpc), M_Fein(100kpc) (Msol) = ',3(1pe12.4))

!!      print*,'TIME (Gyr) = ',time/3.156e16

      open(21,file='diff.dat',status='unknown')
      do j=2,jmax
         write(21,3000)rr(j)/cmkpc,zfe(j)/zfesol
      enddo
      close(21)
3000  format(2(1pe12.4))

stop
end
   

