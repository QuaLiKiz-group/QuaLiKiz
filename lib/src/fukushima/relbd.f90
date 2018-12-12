subroutine relbd(phi,phic,mc,b,d)
!
!   Single precision general incomplete elliptic integrals of the second kind
!
!     Reference: T. Fukushima, (2011) J. Comp. Appl. Math., 235, 4140-4148
!        "Precise and Fast Computation of General Incomplete Elliptic Integral
!         of Second Kind by Half and Double Argument Transformations"       
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
!     Used subprograms: relsbd, relcbd, rcelbd, rserbd
!
!     Inputs: phi  = argument                0 <= phi  <= PI/2
!             phic = complementar argument   0 <= phic <= PI/2
!             mc   = complementary parameter 0 <= mc   <= 1
!
!     Outputs: b, d
!
!     CAUTION: phi and phic must satisfy condition, phi + phic = PI/2
!
real*4 phi,phic,mc,b,d
real*4 m,c,x,d2,z,k,bc,ec,dc,sz,v,t2

if(phi.lt.1.250) then
! write(*,*) "relsbd"
    call relsbd(sin(phi),mc,b,d)
! write(*,*) "relsbd:b,d=",b,d
else
    m=1.0-mc
! write(*,'(a10,1pe25.17)') "m=",m
    c=sin(phic)
    x=c*c
    d2=mc+m*x
    if(x.lt.0.9*d2) then
        z=c/sqrt(d2)
! write(*,*) "relsbd"
        call relsbd(z,mc,b,d)
! write(*,*) "relsbd:b,d=",b,d
        call rcelbd(mc,bc,dc)
        sz=z*sqrt(1.0-x)
        b=bc-(b-sz)
        d=dc-(d+sz)
    else
        v=mc*(1.0-x)
        if(v.lt.x*d2) then
! write(*,*) "relcbd"
            call relcbd(c,mc,b,d)
! write(*,*) "relcbd:b,d=",b,d
        else
! write(*,*) "relcbd"
            t2=(1.0-x)/d2
! write(*,*) "relcbd"
            call relcbd(sqrt(mc*t2),mc,b,d)
! write(*,*) "relcbd:b,d=",b,d
            call rcelbd(mc,bc,dc)
            sz=c*sqrt(t2)
            b=bc-(b-sz)
            d=dc-(d+sz)
        endif
    endif
endif
return
end
