program xelbd
!   test driver for "elbd"
real*8 PI,PIHALF
parameter (PI=3.1415926535897932384626433d0)
parameter (PIHALF=1.5707963267948966192313216916398d0)
real*8 dmc,dphi,mc,mm,phi,phic,b,d
real*4 rmc,rphi,rphic,rb,rd
integer jend,iend,j,i
!
jend=6
iend=5
dmc=1.d0/dble(jend-1)
dphi=PIHALF/dble(iend)
write(*,'(1x,2a10,2a25,2a15)') &
    'm','phi/PI','elb','eld','relb','reld'
do j=1,jend
    write(*,'(1x)')
    mc=dble(j-1)*dmc
    if(mc.le.0.d0) mc=1.21d-32
    rmc=mc
    mm=1.d0-mc
    do i=0,iend
        phi=dphi*dble(i)
        phic=dphi*dble(iend-i)
        rphi=phi
           rphic=phic
                call elbd(phi,phic,mc,b,d)
                call relbd(rphi,rphic,rmc,rb,rd)
                write(*,'(1x,0p2f10.5,0p2f25.15,0p2f15.7)') &
                mm,phi/PI,b,d,rb,rd
    enddo
enddo
end program xelbd
