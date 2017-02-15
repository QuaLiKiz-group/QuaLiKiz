!Contains routines used to build the quasilinear fluxes based on the dispersion relation solution
!The integrals are done 
MODULE QLflux
  USE kind
  USE datmat
  USE datcal
  USE callpassQLints
  USE calltrapQLints

  IMPLICIT NONE

CONTAINS

  SUBROUTINE make_QLflux(p, nu, omega)
    !Extracts the QL linear response based on solution

    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX, INTENT(IN)  :: omega

    COMPLEX(KIND=DBL) :: fonctce, fonctcgte, fonctcgne, fonctcgue, fonctcce, fonctece, fonctvce, fonctecgte, fonctecgne, fonctecgue, fonctecce
    COMPLEX(KIND=DBL) :: fonctpe, fonctpgte, fonctpgne, fonctpgue, fonctpce, fonctepe, fonctvpe, fonctepgte, fonctepgne, fonctepgue, fonctepce
    COMPLEX(KIND=DBL), DIMENSION(nions) :: fonctci, fonctcgti, fonctcgni, fonctcgui, fonctcci, foncteci, fonctvci, fonctecgti, fonctecgni, fonctecgui, fonctecci
    COMPLEX(KIND=DBL), DIMENSION(nions) :: fonctpi, fonctpgti, fonctpgni, fonctpgui, fonctpci, fonctepi, fonctvpi, fonctepgti, fonctepgni, fonctepgui, fonctepci

    IF (fc(p)==0. .OR. REAL(mwidth) < d/3. .OR. (calccirc .EQV. .FALSE.) ) THEN
       fonctce = 0.
       fonctci(:) = 0.
       fonctcgte = 0.
       fonctcgti(:) = 0.
       fonctcgne = 0.
       fonctcgni(:) = 0.
       fonctcgue = 0.
       fonctcgui(:) = 0.
       fonctcce = 0.
       fonctcci(:) = 0.

       fonctece = 0.
       foncteci(:) = 0.

       fonctvce = 0.
       fonctvci(:) = 0.

       fonctecgte = 0.
       fonctecgti(:) = 0.
       fonctecgne = 0.
       fonctecgni(:) = 0.
       fonctecgue = 0.
       fonctecgui(:) = 0.
       fonctecce = 0.
       fonctecci(:) = 0.

    ELSE

       IF ( (rotflagarray(p) == 1) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) THEN
          CALL passQLintsrot( p, nu, omega, fonctce, fonctci, fonctcgte, fonctcgti, fonctcgne, fonctcgni, fonctcgue, fonctcgui, &
               & fonctcce, fonctcci, fonctece, foncteci, fonctecgte, fonctecgti, fonctecgne, fonctecgni, fonctecgue, fonctecgui, fonctecce, fonctecci, fonctvce, fonctvci)        
       ELSE
          CALL passQLints( p, nu, omega, fonctce, fonctci, fonctcgte, fonctcgti, fonctcgne, fonctcgni, &
               & fonctcce, fonctcci, fonctece, foncteci, fonctecgte, fonctecgti, fonctecgne, fonctecgni, fonctecce, fonctecci) 
          fonctcgue = 0.
          fonctcgui(:) = 0.
          fonctvce = 0.
          fonctvci(:) = 0.
          fonctecgue = 0.
          fonctecgui(:) = 0.
       END IF

    END IF

    IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
       fonctpe = 0.
       fonctpi(:) = 0. 
       fonctpgte = 0.
       fonctpgti(:) = 0. 
       fonctpgne = 0.
       fonctpgni(:) = 0. 
       fonctpgue = 0.
       fonctpgui(:) = 0.
       fonctpce = 0.
       fonctpci(:) = 0. 

       fonctepe = 0.
       fonctepi(:) = 0. 

       fonctvpe = 0.
       fonctvpi(:) = 0.

       fonctepgte = 0.
       fonctepgti(:) = 0. 
       fonctepgne = 0.
       fonctepgni(:) = 0. 
       fonctepgue = 0.
       fonctepgui(:) = 0. 
       fonctepce = 0.
       fonctepci(:) = 0. 

    ELSE      

       IF ( (rotflagarray(p) == 1) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) THEN
          CALL trapQLintsrot( p, nu, omega, fonctpe, fonctpi, fonctpgte, fonctpgti, fonctpgne, fonctpgni, fonctpgue, fonctpgui, &
               & fonctpce, fonctpci, fonctepe, fonctepi, fonctepgte, fonctepgti, fonctepgne, fonctepgni, fonctepgue, fonctepgui, fonctepce, fonctepci, fonctvpe, fonctvpi)
       ELSE
          CALL trapQLints( p, nu, omega, fonctpe, fonctpi, fonctpgte, fonctpgti, fonctpgne, fonctpgni, &
               & fonctpce, fonctpci, fonctepe, fonctepi, fonctepgte, fonctepgti, fonctepgne, fonctepgni,fonctepce, fonctepci)
          fonctpgue = 0.
          fonctpgui(:) = 0.
          fonctvpe = 0.
          fonctvpi(:) = 0.
          fonctepgue = 0.
          fonctepgui(:) = 0.
       END IF
    END IF

    IF (el_type == 2) THEN
       fonctce = fc(p) * Nex(p)
       fonctpe = ft(p) * Nex(p)
    ENDIF

    !If rot_flag == 2, then still set momentum transport according to case with symmetry breaking
    IF ((rho(p)<1.0) .AND. (rotflagarray(dimx) == 1) .AND. (ETG_flag(nu) .EQV. .FALSE.)) THEN
       IF (rot_flag == 2) THEN
          Machi=Machiorig; Aui=Auiorig; gammaE=gammaEorig;
          mwidth=mwidth_rot
          mshift=mshift_rot
          widthhat = ABS(mwidth)**2 / SQRT(REAL(mwidth**2))
          Athe=widthhat*cthe(p)/qRd 
          Athi(:)=widthhat*cthi(p,:)/qRd
       ENDIF
       CALL momtrapQLintsrot( p, nu, omega, fonctvpi)
       CALL mompassQLintsrot( p, nu, omega, fonctvci)
    ENDIF

    fonxad = Ac(p)

    fonxcirce = fonctce
    fonxpiege = fonctpe
    fonxcirci(:) = fonctci(:)
    fonxpiegi(:) = fonctpi(:)

    !To save memory space, this condition is in place
    IF (phys_meth .NE. 0.0) THEN
       fonxcircgte = fonctcgte
       fonxpieggte = fonctpgte
       fonxcircgti(:) = fonctcgti(:)
       fonxpieggti(:) = fonctpgti(:)
       fonxcircgne = fonctcgne
       fonxpieggne = fonctpgne
       fonxcircgni(:) = fonctcgni(:)
       fonxpieggni(:) = fonctpgni(:)
       fonxcircgue = fonctcgue
       fonxpieggue = fonctpgue
       fonxcircgui(:) = fonctcgui(:)
       fonxpieggui(:) = fonctpgui(:)
       fonxcircce = fonctcce
       fonxpiegce = fonctpce
       fonxcircci(:) = fonctcci(:)
       fonxpiegci(:) = fonctpci(:)
       !!
       IF (phys_meth == 2) THEN
          fonxecircgte = fonctecgte
          fonxepieggte = fonctepgte
          fonxecircgti(:) = fonctecgti(:)
          fonxepieggti(:) = fonctepgti(:)
          fonxecircgne = fonctecgne
          fonxepieggne = fonctepgne
          fonxecircgni(:) = fonctecgni(:)
          fonxepieggni(:) = fonctepgni(:)
          fonxecircgue = fonctecgue
          fonxepieggue = fonctepgue
          fonxecircgui(:) = fonctecgui(:)
          fonxepieggui(:) = fonctepgui(:)
          fonxecircce = fonctecce
          fonxepiegce = fonctepce
          fonxecircci(:) = fonctecci(:)
          fonxepiegci(:) = fonctepci(:)
       ENDIF
    ENDIF

    fonxecirce = fonctece
    fonxepiege = fonctepe
    fonxecirci(:) = foncteci(:)
    fonxepiegi(:) = fonctepi(:)     

    fonxvcirce = fonctvce
    fonxvpiege = fonctvpe
    fonxvcirci(:) = fonctvci(:)
    fonxvpiegi(:) = fonctvpi(:)     

  END SUBROUTINE make_QLflux

END MODULE QLflux
 
