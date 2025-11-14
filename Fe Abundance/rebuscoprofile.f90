parameter(jmax=5000)

real*8, dimension(jmax) :: rr, rho_rebusco, temp_rebusco, Te_keV, r
real*8 :: rmin, rmax, cmkpc, ne_rebusco


integer :: j


cmkpc = 3.084d21
rmin = 0.*cmkpc
rmax = 3000.d0*cmkpc


!griglia
do j=1,jmax
   r(j)=rmin+(j-1)*rmax/(jmax-1)
enddo
do j=1,jmax-1
   rr(j)=r(j)+0.5*(r(j+1)-r(j))
enddo
rr(jmax)=rr(jmax-1)+(rr(jmax-1)-rr(jmax-2))
!salva griglia
open(10,file='grid_reb.dat',status='unknown')
do j=1,jmax
   write(10,*)real(r(j)/cmkpc),real(rr(j)/cmkpc)
enddo
close(10)

! profilo Rebusco
open(20,file='density_reb.dat',status='unknown')
do j=1,jmax
   rkpc = rr(j)/cmkpc
   ! Densità elettronica ne
   ne_rebusco = 4.6d-2 / (1.d0 + (rkpc/57.d0)**2.d0)**1.8d0 + 4.8d-3 / (1.d0 + (rkpc/200.d0)**2.d0)**0.87d0
   rho_rebusco(j) = 1.937d-24 * ne_rebusco  ! in g/cm^3
   ! Temperatura in keV 
   Te_keV(j) = 7.d0 * (1.d0 + (rkpc/71.d0)**3) / (2.3d0 + (rkpc/71.d0)**3) !errore commesso prima era nello scrivere 2.3 su fortran si scrive 2.3d0 
   !Temp in Kelvin
   temp_rebusco(j) = Te_keV(j) * 1.16d7 
  
   write(20,1000) rr(j)/cmkpc, rho_rebusco(j), Te_keV(j), temp_rebusco(j)
enddo
close(20)

1000 format(4((1pe12.4)))
end


