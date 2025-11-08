parameter(jmax=5000)

real*8, dimension(jmax) :: rr, rho_rebusco, temp_rebusco, Te_keV, r
real*8 :: rmin, rmax, cmkpc, ne_rebusco

cmkpc = 3.084d21
rmin = 0.*cmkpc
rmax = 3000.d0*cmkpc
 
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

zfesol=1.8d-3

!griglia
do j=1,jmax
   r(j)=rmin+(j-1)*rmax/(jmax-1)
enddo
do j=1,jmax-1
   rr(j)=r(j)+0.5*(r(j+1)-r(j))
enddo
rr(jmax)=rr(jmax-1)+(rr(jmax-1)-rr(jmax-2))
!salva griglia
open(10,file='grid_ferrodiff.dat',status='unknown')
do j=1,jmax
   write(10,*)real(r(j)/cmkpc),real(rr(j)/cmkpc)
enddo
close(10)

!program diffusion_fe
!Comparison between our profiles of density and temperature and the Rebusco profiles!



! profilo Rebusco
open(20,file='density_reb.dat',status='unknown')
do j=1,jmax
   rkpc = rr(j)/cmkpc
   ! Densità elettronica ne_rebusco
   ne_rebusco = 4.6d-2 / (1.d0 + (rkpc/57.d0)**2.d0)**1.8d0 + 4.8d-3 / (1.d0 + (rkpc/200.d0)**2.d0)**0.87d0
   rho_rebusco(j) = 1.937d-24 * ne_rebusco  ! in g/cm^3
   ! Temperatura in keV 
   Te_keV(j) = 7.d0 * (1.d0 + (rkpc/71.d0)**3) / (2.3d0 + (rkpc/71.d0)**3)
   !Temp in Kelvin
   temp_rebusco(j) = Te_keV(j) * 1.16d7 
   write(20,1000) rr(j)/cmkpc, rho_rebusco(j), Te_keV(j), temp_rebusco(j)
enddo
close(20)

! perseus' profile
! Temperature profile

 do j=1,jmax
    rkpc=rr(j)/cmkpc
    roverr0 = rkpc/r0bill
    temp1 = t0bill + tcoeffbill * roverr0**pbill
    cut=0.0
    temp2 = tmaxbill/(cut*(rkpc/60.)**2 + roverr0**qbill)
    ttt = 1./( (1./temp1)**sbill + (1./temp2)**sbill )**sinv
    tcorr = -0.12*exp(-rkpc/1.5)
    tem(j) = (ttt + tcorr)*1.e7    !! this is for NGC 5044 !!
    x=rr(j)/r500
    xx=x/0.045
    tem2(j)=ticm*1.35*(xx**1.9+0.45)/(xx**1.9+1.)* &   !! this is for Perseus !!
        1./(1.+(x/0.6)**2)**0.45
    tem(j)=tem2(j)  !!ticm --> for Perseus
    write(20,1001)rr(j)/cmkpc, tem(j), tem2(j), ticm
 enddo
close(20)
1001 format(4(1pe12.4))

!     calculate the gas density, assuming ticm

lnd(1)=log(rho0)          !! mette il gas in eq. con il potenziale
do j=2,jmax
   gg=grvnfw(j)
   temmed=0.5*(tem(j)+tem(j-1))
!! isoth !!   lnd(j)=lnd(j-1)-gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*ticm)
   lnd(j)=lnd(j-1)-gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*temmed) &
          - (log(tem(j)) - log(tem(j-1))) 
enddo

do j=1,jmax
   rho(j)=exp(lnd(j))
enddo

open(20,file='density_irondiff.dat',status='unknown')
do j=1,jmax
   write(20,1101)rr(j)/cmkpc,rho(j),rhonfw(j),rhost(j)
enddo
close(20)

open(10, file='iron_diffusion.dat' )

!We are creating the data file for the diffusion!
zfe_ground=0.4*zfesol   !! this is the background abundance !!
do j = 1, jmax
!Rescaled radius!
 x=rr(j)/(80.*cmkpc)
!Computing the Fe abundance profile in Perseus Cluster ! 
  !!Perseus: here I would remove the factor 1.15/1.15, because for definition is equal to 1!!
 zfe_obs(j)=zfesol*0.3*1.4*1.15*(2.2+x**3)/(1+x**3)/1.15  !!Perseus' profile-->observed value of iron metallicity

!Cleaning operation from the abundance_fe di background!
  !!zfe_ex is reffering to data of the Fe abundance without the background contribution!!  
   zfe_ex(j)=zfe_obs(j) - zfe_ground !zfe_ex is the excess of iron adundance 
   zfe_ex(j)=max(zfe_obs_ex(j),0.) !calcola il max-->per successivo controllo!

!Here we are considering the different initial values for the Zfe!
  !!Investigating just diffusion:zfe_obs(j), meaning the abundance at radius corrisponding to the index j !!
   rho_fe_ex(j)=rho(j)*zfe_obs(j)/1.4
  !!Investigating only the action of SNIa, responsible for the Fe-peak!!
   !!!Initial condition: zfe(r,0) = 0 vs zfe(r,0) = zfe_ground!!!
   !zfe_s(j)= qui ci andrà l'implementazione della evoluzione temporale dell'abbondanza nel caso di sola sorgente!  
   rhofe_s(j)= 0.
   
enddo

!!!!! Setting the initial condition, in order to guaranteeing continuity during the computation in the cycle with a correct managing of the j index!
 !boundary conditions (outflows)
z_fe(1)=z_fe(2) 
z_fe(jmax)=z_fe(jmax-1)
!!!!!!!! correct values of density at the extremes!!!!
rho_fe(1)=rho(1)*z_fe(1)/1.4
rho_fe(jmax)=rho(jmax)*z_fe(jmax)/1.4

!!!!!!!!!!!!!!Computing the grad Zfe!!!!!!!!!!!!!!
!!!!!!!!!!!!Remember i = j +1/2!!!!!!!!!!!!!!!
 do  i=2,jmax-1
    gradz_fe(i)=(z_fe(i)-z_fe(i-1))/(rr(i)-rr(i-1))  !! dZ/dr centered at "j" !!
 enddo
 
 gradz_fe(1)= 0. 
 gradz_fe(jmax)= 0.

  do j=2,jmax-1
  !! Remember we are computing the average value of the density which has been defined on the grid in witch rr(i) where i = j + 1/2
    rho_jp1=0.5*(rho(j+1)+rho(j)) !! media tra j+1/2 e j+3/2 → centrata in j+1
    rho_j=0.5*(rho(j-1)+rho(j))  !! media tra j-1/2 e j+1/2 → centrata in j
    
  !! So now we have the gas density on the grid r(j), and note that they are numeric values not vectors
    rho_fe(j)=rho_fe(j) &
            + (dt/1.4)*(r(j+1)**2*D*rho_jp1*gradz_fe(j+1) &
            -r(j)**2*D*rho_j*gradz_fe(j))   &
             / (0.33333333*(r(j+1)**3-r(j)**3))
         z_fe(j)=1.4*rho_fe(j)/rho(j)  !! update Z_Fe with the new rho_Fe !! !! fomula a pag 28 di project1 pdf
  !!!!!!!!!!!So now we are happy because we have the Fe density !!!!!!!!!!!!!!!!!
   write (10, 1000 ) !!!!!!!!!!!!!Remember to change the format number
  enddo


!!!!!!! end abundance cycle!!!!!!!!!!!!!!!!!!
