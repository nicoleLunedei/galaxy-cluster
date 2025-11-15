!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Nuova versione corretta: ciclo su n invece che su time, stampa 
!! su file ed evoluzione funzionante
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

parameter(jmax=5000)
implicit real*8 (a-h,o-z)
!!!Declaration of real arrays of dimension of jmax!!!
real*8, dimension(jmax) :: r,rr,vol,mnfw,&
        rhost,rho,mhern,rhonfw,mdark,&
        grvnfw,lnd, mdm_analytic, T, lnT,Mgas,&
        ne_reb, rho_reb, T_reb, zfe_obs,& !!!!Variables for diffusion program
        zfe_ex, zfe_s, rho_fe_ex, rhofe_s, z_fe, rho_fe, gradz_fe,rho_fe_obs


!!!!Declaration of real numerical values!!!!
real*8 :: msol,mu,mp,rmin,rmax,mvir,rvir,mbgc,ahern,fb,rkpc,zfesol,zfe_ground,&
          rho_jp1,rho_j, dt, D, t_final, time, n_cycle, v_t, l_t
!constants
msol = 1.989d33
cmkpc = 3.084e21
mu=0.62
boltz=1.38066e-16
guniv=6.6720e-8
mp=1.67265e-24
years=3.156d7

!set the grid
rmin = 0.*cmkpc
rmax = 3000.*cmkpc
!Original domain!
do j=1,jmax
   r(j)=rmin+(j-1)*rmax/(jmax-1)    !!r_j
enddo
!Shifted domain!
do j=1,jmax-1
   rr(j)=r(j)+0.5*(r(j+1)-r(j))     !!r_j+1/2
end do
rr(jmax)=rr(jmax-1)+(rr(jmax-1)-rr(jmax-2))

open(10,file='grid.dat',status='unknown')
do j=1,jmax
   write(10,*)real(r(j)/cmkpc),real(rr(j)/cmkpc)
end do
close(10)

!!volumes of centered shells on the shifted domain!!
vol(1)=4.1888*r(1)**3
do j=2,jmax
   vol(j)=4.1888*(r(j)**3-r(j-1)**3)   
end do

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
!!Cycle on the Shifted domain!!
do j=1,jmax
   x=rr(j)/rs
   rhonfw(j)=rho0nfw/(x*(1.+x)**2)
end do

open(20,file='masse_diff.dat')
mnfw(1)=0.		!!initial dark matter mass
mdark(1)=0.
mhern(1)=0.
mdm_analytic(1)=0.
mdm_analytic(1)=0.
!!Cycle on the Original domain!!
do j=2,jmax 		!!NB: starts at j=2 not j=1
   x=r(j)/rs
   mnfw(j)=mnfw(j-1)+rhonfw(j-1)*vol(j)					                     !!DM mass discrete NFW formula in pdf centered in j-1/2
   mdm_analytic(j)=4*3.14*rho0nfw*(rs**3)*(log(1.+x)-(r(j)/(r(j)+rs)))     !!DM analytical
   mdark(j)=mvir*(log(1.+x)-x/(1.+x))/fc				                        !!DM analytic NFW formula
   mhern(j)=mbgc*r(j)**2/(r(j)+ahern)**2				                        !!analytic stellar mass
end do
!!Cycle on the Original domain!!
do j=1, jmax
   write(20,1001)r(j)/cmkpc,mnfw(j)/msol,mdark(j)/msol,mhern(j)/msol, mdm_analytic(j)/msol
end do
1001 format(5(1pe12.4))
close(20)

!! Temperature profile from a given analytical function !!
open(20, file='temperature_diff.dat')
do j=1, jmax
    y=rr(j)/r500
    T(j)=ticm*1.35d0*((y/0.045d0)**1.9+0.45d0)/((y/0.045d0)**1.9+1)*(1+(y/0.6d0)**2)**(-0.45)
    lnT(j)=log(T(j))
    write(20,1003)rr(j)/cmkpc, T(j), ticm
end do
1003 format(3(1pe12.4))
close(20)

!!gravitational energy
open(20,file='grv_diff.dat')
grvnfw(1)=0.          !! ok per alone NFW, isotermo o beta-model
!!!!Cycle on the Original domain!!
do j=2,jmax
   grvnfw(j)=guniv*(mnfw(j)+mhern(j))/r(j)**2		!!already takes into account the mass of the central galaxy BCG
   write(20,1002)r(j)/cmkpc,grvnfw(j)/msol
end do
1002 format(2(1pe12.4))
close(20)

!calculate the gas density
lnd(1)=log(rho0)          !! mette il gas in eq. con il potenziale
!!Cycle on the Shifted domain for the density and the original domain for the temperature!!
do j=2,jmax
   gg=grvnfw(j)
   lnd(j)=lnd(j-1)-(gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*T(j))+ lnT(j)-lnT(j-1))
end do
!!Just passing from lnd(j) to rho(j)!!
do j=1,jmax
   rho(j)=exp(lnd(j))  
end do

!!Cycle on the Shifted domain!!
Mgas(1)=0.
do j=2,jmax
   Mgas(j)=Mgas(j-1)+rho(j-1)*vol(j)
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NOW THE REAL FE DIFFUSION PART STARTS!!!!!!!!!!!!!!!!!!!!!!!!!!!! ************************************************************************************************+
!Comparison between our profiles, of density and temperature, and the Rebusco profiles!
 open(10,file='comparison.dat')

 !!Cycle on the Original domain!!
  !!In this cycle we are calculating the gas density and the temperature profile used by Rebusco, from XMM Newton data!!
  zfesol= 1.8d-3
  zfe_ground = 0.4 * zfesol
  v_t = 200d6
  l_t = 20 * cmkpc
  KeV = 1.16d7
  !D = 0.333333333 * v_t * l_t   !Our value
  D=0.11*260d5*15*cmkpc          !Rebusco value

 do j = 1, jmax
  rkpc = rr(j)/cmkpc
  ne_reb(j) = ((4.6d-2)/(1 + (rkpc/57)**2)**1.8) + ((4.8d-3)/(1 + (rkpc/200)**2)**0.87)
  rho_reb(j) = 1.937d-27 * ne_reb(j)
  T_reb(j) = 7 * ((1+(rkpc/71)**3)/(2.3+(rkpc/71)**3)) * KeV
  write(10,1000)rkpc,rho(j), rho_reb(j), T_reb(j)
 end do 
 1000 format(4(1pe12.4))

!!!!!!!!!!!Saving the initial values!!!!!!!!!!!!
open(20, file='initial_zfe.dat' )

!!!!!!!!!!!Saving the final values!!!!!!!!
open(40, file= 'final_zfe.dat')

***************!Creating the initial condition data file!****************
do j = 1, jmax
 !Rescaled radius!
    x=rr(j)/(80.*cmkpc)
    
 !DIFFUSION CASE!**********
 !!Computing the Fe abundance profile in Perseus Cluster ! 
 !!Perseus profile: here I would remove the factor 1.15/1.15, because for definition is equal to 1!!
   zfe_obs(j) = zfesol * 0.3 * 1.4 * 1.15 * (2.2 + x**3)/(1 + x**3)/1.15  

 !!!Cleaning operation from the abundance_fe di background!
   !!!!zfe_ex is reffering to data of the Fe abundance without the background contribution!!  
   zfe_ex(j) = zfe_obs(j) - zfe_ground 
   zfe_ex(j)= max(zfe_ex(j),0.) !calcola il max-->per successivo controllo! Forse lo fa perchè la differenza potrebbe andare in negativo 
                                !quindi sostituisce l'eccesso negativo con lo 0 per scanso di equivoci
   !!!!zfe_ex will be the initial condition for the diffusion case
 !!Computing the Fe density: t = 0!!
   rho_fe_obs(j) = rho(j) * zfe_obs(j)/1.4    !!!From raw data!!!
   rho_fe_ex(j) = rho(j) * zfe_ex(j)/1.4      !!!From cleaned data!!! 

 !SOURCE CASE: responsible for the Fe-peak!******************
   !!!Initial condition: zfe(r,0) = 0 vs zfe(r,0) = zfe_ground!!!
   zfe_s(j) = 0.  
   rhofe_s(j)= 0.
 write(20,1004) rr(j)/cmkpc, zfe_obs(j), zfe_ex(j), zfe_s(j)
enddo
1004 format(4(1pe12.4))

!Print*, zfe_obs(1)
!Print*, zfe_obs(100)

!!!!!!!!!!!!!!!!!!Computing the initial mass of the iron profile""""""""""""""""""
M0_fe(1) = 0.
   do j=2,jmax
     M0_fe(j)=M0_fe(j-1)+rho_fe_ex(j-1)*vol(j) 
   enddo
        Print*, "Initial Fe Mass (< 100kpc)", M0_fe(166)/ msol, "masse solari"


!!Initializing rho_fe and z_fe: we assign their initial value at t=0 equal to the (Rebusco- background) profiles
do j=1, jmax
rho_fe(j)=rho_fe_ex(j)
z_fe(j)=zfe_ex(j)
end do

************************!!!!!!!!!TIME CYCLE!!!!!!!!!!!!!!!!!!!!******************************************
t_final = 1.e9 * years                  !!!!Change to see how it evolves in different timelines
dt = 0.4 * (r(5) - r(4))**2 / (2. * D)
n_cycle=t_final/dt
Print*, "Intervallo temporale ", dt
Print*, "t_final:", t_final
Print*, "Numero cicli:", int(n_cycle)
Print*, "Coefficiente di diffusione ", D
time = 0.

do n=1, int(n_cycle)
!!!!!!!!!!!!!!Computing the grad Zfe!!!!!!!!!!!!!!
 !! Cycle on the Shifted domain!!
  do  j=2,jmax-1
    gradz_fe(j)=(z_fe(j)-z_fe(j-1))/(rr(j)-rr(j-1))  !! dZ/dr centered at "j" !!
  enddo
  gradz_fe(1)= 0. 
  gradz_fe(jmax)= 0.

 !!!Medium values of the gas density in each shell!!!
  do j=2,jmax-1
  !! Remember we are computing tha average value of the density which has been defined on the grid
      rho_jp1 = 0.5 * (rho(j+1) + rho(j)) 
      rho_j = 0.5 * (rho(j-1) + rho(j))  
  !! So now we have the gas density on the grid r(j), and note that they are numeric values not vectors
      rho_fe(j) = rho_fe(j) &
            + (dt/1.4)*(r(j+1)**2*D*rho_jp1*gradz_fe(j+1) &
            -r(j)**2*D*rho_j*gradz_fe(j))   &
             / (0.33333333*(r(j+1)**3-r(j)**3))            
   !!!!!!!!!!!!Getting the abundance from the rho_fe!!!!!!!!!!!!!!!!!!!!!!
      z_fe(j)=1.4*rho_fe(j)/rho(j)  !! update Z_Fe with the new rho_Fe !! !! fomula a pag 28 di project1 pdf
  !!!!!!!!!!!So now we are happy because we have the Fe density !!!!!!!!!!!!!!!!!

  !!!!!!!!!!!Now we are writing the data of different times in different files!!!!!!
 if (n .EQ. int(n_cycle)) then  
   do j=1, jmax
      write (40, 1005) r(j)/cmkpc, rho_fe_obs(j), rho_fe(j), zfe_obs(j)/zfesol, zfe_ex(j)/zfesol, z_fe(j)/zfesol
   end do
   print*, "n=n_cycle"
     !!!!!!!!!!!!!!!!!!!!!!!computing the Fe mass!!!!!!!!!!!!!!!!!!!
  M_fe(1) = 0.
   do j=2,jmax
     M_fe(j)=M_fe(j-1)+rho_fe(j-1)*vol(j)
     write(50,1006) r(j)/cmkpc, M0_fe(j)/msol, M_fe(j)/msol 
   enddo
    Print*, "Final Fe Mass (< 100kpc)", M_fe(166)/ msol, "masse solari"
   print*, "n=n_cycle"
 end if
!!!!!!!!!!Radial cyle finished!!!!!!!!

!!!!! Setting the boundary condition, in order to guaranteeing continuity during the computing in the cycle with a correct managing the j index!
!!(1) = (2), (jmax) = (jmax-1)!!
z_fe(1)=z_fe(2)
z_fe(jmax)=z_fe(jmax-1)
rho_fe(1)=rho_fe(2)
rho_fe(jmax)=rho_fe(jmax-1)

!!Evolving time
time = time + dt
end do
!!!!!!!!!!Time cycle finished!!!!!!!!!
close(40)
1005 format(6 (1pe12.4))

!!!!!!! end abundance cycle!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!THE END :)
end 
