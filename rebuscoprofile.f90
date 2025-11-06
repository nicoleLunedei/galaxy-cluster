parameter(jmax=5000)
real*8, dimension(jmax) :: rr(jmax), rho_rebusco, temp_rebusco
real*8 ::rmin,rmax
real*8, dimension(jmax) :: r
real*8 :: cmkpc, ne_rebusco

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
open(10,file='grid_rebusco.dat',status='unknown')
do j=1,jmax
   write(10,*)real(r(j)/cmkpc),real(rr(j)/cmkpc)
enddo
close(10)

! profilo Rebusco
open(20,file='density_rebusco.dat',status='unknown')
do j=1,jmax
   rkpc = rr(j)/cmkpc
   ! Densità elettronica ne
   ne_rebusco = 4.6d-2 / (1.d0 + (rkpc/57.d0)**2.d0)**1.8d0 + 4.8d-3 / (1.d0 + (rkpc/200.d0)**2.d0)**0.87d0
   rho_rebusco(j) = 1.937d-24 * ne_rebusco  ! in g/cm^3
   ! Temperatura Te(r) from keV to Kelvin
   Te_keV = 7.d0 * (1.d0 + (rkpc/71.d0)**3) / (1.d0 + (rkpc/300.d0)**3)
temp_rebusco(j) = Te_keV * 1.16d7
!temp_rebusco(j) = 7.d0 * (1.d0 + (rkpc/71.d0)**3) / (2.d3 + (rkpc/71.d0)**3) * 1.16d7  ! in Kelvin perchè c'è * 1.16d7
  
   write(20,1000) rr(j)/cmkpc, rho_rebusco(j), temp_rebusco(j)
enddo
close(20)

1000 format(3(1pe13.5))
end


