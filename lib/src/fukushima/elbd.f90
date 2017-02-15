subroutine elbd(phi,phic,mc,b,d)
!
!   Double precision general incomplete elliptic integrals of the second kind
!
!     Reference: T. Fukushima, (2011) J. Comp. Appl. Math., 235, 4140-4148
!        "Precise and Fast Computation of General Incomplete Elliptic Integral
!         of Second Kind by Half and Double Argument Transformations"       
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
!     Used subprograms: elsbd, elcbd, celbd, serbd
!
!     Inputs: phi  = argument                0 <= phi  <= PI/2
!             phic = complementar argument   0 <= phic <= PI/2
!             mc   = complementary parameter 0 <= mc   <= 1
!
!     Outputs: b, d
!
!     CAUTION: phi and phic must satisfy condition, phi + phic = PI/2
!
real*8 phi,phic,mc,b,d
real*8 m,c,x,d2,z,k,bc,ec,dc,sz,v,t2

if(phi.lt.1.25d0) then
    call elsbd(sin(phi),mc,b,d)
else
    m=1.d0-mc
    c=sin(phic)
    x=c*c
    d2=mc+m*x
    if(x.lt.0.9d0*d2) then
        z=c/sqrt(d2)
        call elsbd(z,mc,b,d)
        call celbd(mc,bc,dc)
            sz=z*sqrt(1.d0-x)
            b=bc-(b-sz)
            d=dc-(d+sz)
      else
            v=mc*(1.d0-x)
            if(v.lt.x*d2) then
            call elcbd(c,mc,b,d)
        else
            t2=(1.d0-x)/d2
            call elcbd(sqrt(mc*t2),mc,b,d)
            call celbd(mc,bc,dc)
                  sz=c*sqrt(t2)
                  b=bc-(b-sz)
                  d=dc-(d+sz)
            endif
      endif
endif
return
end
