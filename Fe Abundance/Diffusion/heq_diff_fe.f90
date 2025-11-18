!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!Codice per il caso solo diffusione
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

parameter(jmax=5000)
implicit real*8 (a-h,o-z)
!!!Declaration of real arrays of dimension of jmax!!!
real*8, dimension(jmax) :: r,rr,vol,mnfw,&
        rho,mhern,rhonfw,mdark,&
        grvnfw,lnd, mdm_analytic, T, lnT,Mgas,&
        ne_reb, rho_reb, T_reb, zfe_obs,& !!!!Variables for diffusion program
        zfe_ex, rho_fe_ex, z_fe, rho_fe, gradz_fe, rho_fe_obs, M0_fe, M_fe
   

!!!!Declaration of real numerical values!!!!
real*8 :: msol,mu,mp,rmin,rmax,mvir,rvir,mbgc,ahern,&
          zfesol,zfe_ground,rho_jp1,rho_j, D, v_t, l_t,&
          t_final, time, n_cycle, dt

!constants and units of measure
msol = 1.989d33
cmkpc = 3.084e21
mu=0.62
boltz=1.38066e-16
guniv=6.6720e-8
mp=1.67265e-24
years=3.156d7
KeV = 1.16d7

!parameters
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
zfesol= 1.8d-3
zfe_ground = 0.4 * zfesol  !!background Z
v_t = 200d6
l_t = 20 * cmkpc
!!Diffusion coefficient
!D = 0.333333333 * v_t * l_t
D=0.11*260d5*15*cmkpc   !!Rebusco

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
enddo
rr(jmax)=rr(jmax-1)+(rr(jmax-1)-rr(jmax-2))

!!volumes of centered shells on the shifted domain!!
vol(1)=4.1888*r(1)**3
do j=2,jmax
   vol(j)=4.1888*(r(j)**3-r(j-1)**3)   
enddo

!!DM density profile
!!Cycle on the Shifted domain!!
do j=1,jmax
   x=rr(j)/rs
   rhonfw(j)=rho0nfw/(x*(1.+x)**2)
enddo

mnfw(1)=0.
mdark(1)=0.
mhern(1)=0.
mdm_analytic(1)=0.
mdm_analytic(1)=0.
!!Cycle on the Original domain!!
do j=2,jmax 	
   x=r(j)/rs
   mnfw(j)=mnfw(j-1)+rhonfw(j-1)*vol(j)					                     !!DM mass discrete NFW formula in pdf centered in j-1/2
   mdm_analytic(j)=4*3.14*rho0nfw*(rs**3)*(log(1.+x)-(r(j)/(r(j)+rs)))     !!DM analytical
   mdark(j)=mvir*(log(1.+x)-x/(1.+x))/fc				                        !!DM analytic NFW formula
   mhern(j)=mbgc*r(j)**2/(r(j)+ahern)**2				                        !!analytic stellar mass
enddo

!!Temperature profile from a given analytical function !!
do j=1, jmax
    y=rr(j)/r500
    T(j)=ticm*1.35d0*((y/0.045d0)**1.9+0.45d0)/((y/0.045d0)**1.9+1)*(1+(y/0.6d0)**2)**(-0.45)
    lnT(j)=log(T(j))
enddo

!!Gravitational energy
grvnfw(1)=0.
!!!!Cycle on the Original domain!!
do j=2,jmax
   grvnfw(j)=guniv*(mnfw(j)+mhern(j))/r(j)**2		!!already takes into account the mass of the central galaxy BCG
enddo

!Calculate the gas density
lnd(1)=log(rho0)       
!!Cycle on the Shifted domain for the density and the original domain for the temperature!!
do j=2,jmax
   gg=grvnfw(j)
   lnd(j)=lnd(j-1)-(gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*T(j))+ lnT(j)-lnT(j-1))
enddo
do j=1,jmax
   rho(j)=exp(lnd(j))  
enddo

!!Cycle on the Shifted domain!!
Mgas(1)=0.
do j=2,jmax
   Mgas(j)=Mgas(j-1)+rho(j-1)*vol(j)
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NOW THE REAL FE DIFFUSION PART STARTS!!!!!!!!!!!!!!!!!!!!!!!!!!!! ************************************************************************************************
!!!!!!!!!!!Files for the final values!!!!!!!!
open(20, file='final_zfe_diff1.dat')
open(30, file='final_zfe_diff2.dat')
open(40, file= 'final_zfe_diff5.dat')
open(50, file= 'Mass_fe_diff.dat')

do j = 1, jmax
 !Rescaled radius!
    x=rr(j)/(80.*cmkpc)

 !DIFFUSION CASE!
 !!Computing the Fe abundance profile in Perseus Cluster ! 
  !!Perseus: here I would remove the factor 1.15/1.15, because for definition is equal to 1!!
   zfe_obs(j) = zfesol * 0.3 * 1.4 * 1.15 * (2.2 + x**3)/(1 + x**3)/1.15  

   !!!Cleaning operation from the abundance_fe di background!
    !!!!zfe_ex is reffering to data of the Fe abundance without the background contribution!!  
   zfe_ex(j) = zfe_obs(j) - zfe_ground 
   zfe_ex(j)= max(zfe_ex(j),0.) !calcola il max-->per successivo controllo! Forse lo fa perchè la differenza potrebbe andare in negativo 
                                !quindi sostituisce l'eccesso negativo con lo 0 per scanso di equivoci
    !!!!zfe_ex will be the initial condition for the diffusion case
 !!Computing the Fe density: t = 0!!
   
  !!!From raw data!!!
   rho_fe_obs(j) = rho(j) * zfe_obs(j)/1.4

  !!!From cleaned data!!!
   rho_fe_ex(j) = rho(j) * zfe_ex(j)/1.4
end do
!!!!!!!!!!!!!!!!!!Computing the initial mass of the iron profile""""""""""""""""""
M0_fe(1) = 0.
do j=2,jmax
     M0_fe(j)=M0_fe(j-1)+rho_fe_ex(j-1)*vol(j) 
end do

!!!!! Setting the initial condition
do j=1, jmax
   rho_fe(j)=rho_fe_ex(j)
   z_fe(j)=zfe_ex(j)
enddo

!!!!!!!!!TIME CYCLE!!!!!!!!!!!!!!!!!!!!***************************************************************************
print *, 'Final time: [yr] (1-2-5d9)'
read(*, *) t_final
t_final = t_final* years
dt = 0.4 * (r(5) - r(4))**2 / (2. * D)
n_cycle=t_final/dt
time = 0.

do n=1, int(n_cycle)
!!!!!!!!!!!!!!Computing the grad Zfe!!!!!!!!!!!!!!
 !! Cycle on the Shifted domain!!
 !!!Gradient of Zfe!!!
  do  j=2,jmax-1
    gradz_fe(j)=(z_fe(j)-z_fe(j-1))/(rr(j)-rr(j-1))  !! dZ/dr centered at "j" !!
  end do
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

   !!!!!!!!!!!!!!!!!!!!computing the Fe mass!!!!!!!!!!!!!!!!!!!
      M_fe(1) = 0.
      M_fe(j)=M_fe(j-1)+rho_fe(j-1)*vol(j)

 end do
!!!!!!!!!!Radial cyle finished!!!!!!!!**************************

!!Boundary condition (1) = (2), (jmax) = (jmax-1)!!
!ABUNDANCE!
z_fe(1) = z_fe(2) 
z_fe(jmax) = z_fe(jmax-1)

!FE DENSITY!
rho_fe(1) = rho_fe(2)
rho_fe(jmax) = rho_fe(jmax-1)

time = time + dt
end do
!!!!!!!!!!Time cycle finished!!!!!!!!!**************************

!!to check if at the end time=t_final (note time is always <t_final that's way later we aren't imposing the equality)
print*, 'Final time', time/years    
!!To check MASS CONSERVATION
Print*, "Initial Fe Mass (< 100kpc)=", M0_fe(166)/ msol, "masse solari"     !!element 166 corresponds to radius 100 kpc

!Writing the final values rho_fe, z_fe, M_fe(<100 kpc)
!!Case diffusion after 1Gyr
if ((time/years) .le. 1d9) then     
   do j=1, jmax
      write (20, 1005) r(j)/cmkpc, rho_fe_obs(j), rho_fe(j), zfe_obs(j)/zfesol, zfe_ex(j)/zfesol, z_fe(j)/zfesol
   end do
   print*, 'Final Fe Mass(<100 kpc)=', M_fe(166)/msol, 'masse solari per time=', time/years 
end if
!!Case after 2 Gyr
if ((time/years) .ge. 1d9 .and. (time/years) .le. 2d9) then
   do j=1, jmax
      write (30, 1005) r(j)/cmkpc, rho_fe_obs(j), rho_fe(j), zfe_obs(j)/zfesol, zfe_ex(j)/zfesol, z_fe(j)/zfesol
   end do
   print*, 'Final Fe Mass(<100 kpc)=', M_fe(166)/msol, 'masse solari per time=', time/years 
end if
!!Case after 5 Gyr
if ((time/years) .ge. 2d9 .and. (time/years) .le. 5d9) then
   do j=1, jmax
      write (40, 1005) r(j)/cmkpc, rho_fe_obs(j), rho_fe(j), zfe_obs(j)/zfesol, zfe_ex(j)/zfesol, z_fe(j)/zfesol
   end do
   print*, 'Final Fe Mass(<100 kpc)=', M_fe(166)/msol, 'masse solari per time=', time/years 
end if
close(20)
close(30)
close(40)
1005 format(6 (1pe12.4))

end 
