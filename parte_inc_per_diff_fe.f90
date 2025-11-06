program diffusion_fe
!Comparison between our profiles, of density and temperature, and the Rebusco profiles!
 open(20,file='comparison')
 do j = 1, jmax
  rkpc = rkpc/cmkpc
  ne_re(j) = ((4.6d-2)/(1 + (rkpc/57)**2)**1.8) + ((4.8d-3)/(1 + (rkpc/200)**2)**0.87)
  rho_re(j) = 1.937d-27 * ne_re(j)
  T_re(j) = 
  write(20,1001)rkpc,rho(j), rho_re(j)
 enddo 

open(10, file='diffusion_fe.dat' )
!We are creating the data file for the diffusion!
do j = 1, jmax
!Rescaled radius!
    x=rr(j)/(80.*cmkpc)
!Computing the Fe abundance profile in Perseus Cluster ! 
  !!Perseus: here I would remove the factor 1.15/1.15, because for definition is equal to 1!!
   zfe_obs(j)=zfesol*0.3*1.4*1.15*(2.2+x**3)/(1+x**3)/1.15  

!Cleaning operation from the abundance_fe di background!
  !!zfe_ex is reffering to data of the Fe abundance without the background contribution!!  
   zfe_ex(j)=zfe_obs(j) - zfe_ground 
   !zfe_ex(j)=max(zfe_obs_ex(j),0.) !calcola il max-->per successivo controllo!

!Here we are considering the different initial values for the Zfe!
  !!Investigating just diffusion:zfe_obs(j), meaning the abundance at radius corrisponding to the index j !!
   rho_fe_ex(j)=rho(j)*zfe_obs(j)/1.4
  !!Investigating only the action of SNIa, responsible for the Fe-peak!!
   !!!Initial condition: zfe(r,0) = 0 vs zfe(r,0) = zfe_ground!!!
   !zfe_s(j)= qui ci andrà l'implementazione della evoluzione temporale dell'abbondanza nel caso di sola sorgente!  
   rhofe_s(j)= 0.
   
enddo

!!!!! Setting the initial condition, in order to guaranteeing continuity during the computing in the cycle with a correct managing the j index!

z_fe(1)=z_fe(2) 
z_fe(jmax)=z_fe(jmax-1)
!!!!!!!! correct values of density at the extremes!!!!
rho_fe(1)=rho(1)*z_fe(1)/1.4
rho_fe(jmax)=rho(jmax)*z_fe(jmax)/1.4
!!!!!!!!!!!!!!Computing the grad Zfe!!!!!!!!!!!!!!
!!!!!!!!!!!!Remember i = j +1/2!!!!!!!!!!!!!!!
 do  =2,jmax-1
    gradz_fe(i)=(z_fe(i)-z_fe(i-1))/(rr(i)-rr(i-1))  !! dZ/dr centered at "j" !!
 enddo
 gradz_fe(1)= 0. 
 gradz_fe(jmax)= 0.

  do j=2,jmax-1
  !! Remember we are computing tha average value of the density which has been defined on the grid
  !! rr(i) where i = j + 1/2
    rho_jp1=0.5*(rho(j+1)+rho(j)) 
    rho_j=0.5*(rho(j-1)+rho(j))  
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
