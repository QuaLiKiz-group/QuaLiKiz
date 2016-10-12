program xcei
!
!   Test driver of complete elliptic integrals of first and second kind:
!   K(m), E(m), B(m), D(m), S(m)
!
!     References:
!
!       T. Fukushima, (2015), J. Comp. Appl. Math., 282, 71-76
!
!      "Precise and Fast Computation of General Complete Elliptic
!       Integral of Second Kind"
!
!       T. Fukushima, (2016), Astron. J., re-revised
!
!      "Zonal Toroidal Harmonic Expansions of External Gravitational
!       Fields for Ring-Like Objects"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
!
nend=20
dmc=1.d0/dble(nend)
!
write (*,"(a10,5a25)") "m","K(m)","E(m)","B(m)","D(m)","S(m)"
do n=1,nend
    fmc=dmc*dble(n)
    cek=ceik(fmc)
    cee=ceie(fmc)
    ceb=ceib(fmc)
    ced=ceid(fmc)
    ces=ceis(fmc)
    fm=1.d0-fmc
    write (*,"(0pf10.5,1p5e25.15)") fm,cek,cee,ceb,ced,ces
enddo
!
stop
end program xcei
