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

    COMPLEX(KIND=DBL) :: fonctce, fonctcgte, fonctcgne, fonctcce, fonctece, fonctecgte, fonctecgne, fonctecce
    COMPLEX(KIND=DBL) :: fonctpe, fonctpgte, fonctpgne, fonctpce, fonctepe, fonctepgte, fonctepgne, fonctepce
    COMPLEX(KIND=DBL), DIMENSION(nions) :: fonctci, fonctcgti, fonctcgni, fonctcgui, fonctcci, foncteci, fonctvci, fonctecgti, fonctecgni, fonctecgui, fonctecci
    COMPLEX(KIND=DBL), DIMENSION(nions) :: fonctpi, fonctpgti, fonctpgni, fonctpgui, fonctpci, fonctepi, fonctvpi, fonctepgti, fonctepgni, fonctepgui, fonctepci

    IF (fc(p)==0. .OR. REAL(mwidth) < d/3. .OR. (calccirc .EQV. .FALSE.) ) THEN
       fonctce = 0.
       fonctci(:) = 0.
       fonctcgte = 0.
       fonctcgti(:) = 0.
       fonctcgne = 0.
       fonctcgni(:) = 0.
       fonctcgui(:) = 0.
       fonctcce = 0.
       fonctcci(:) = 0.

       fonctece = 0.
       foncteci(:) = 0.

       fonctvci(:) = 0.

       fonctecgte = 0.
       fonctecgti(:) = 0.
       fonctecgne = 0.
       fonctecgni(:) = 0.
       fonctecgui(:) = 0.
       fonctecce = 0.
       fonctecci(:) = 0.

    ELSE

       IF ( (rotflagarray(p) == 1) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) THEN
          CALL passQLintsrot( p, nu, omega, fonctce, fonctci, fonctcgte, fonctcgti, fonctcgne, fonctcgni, fonctcgui, &
               & fonctcce, fonctcci, fonctece, foncteci, fonctecgte, fonctecgti, fonctecgne, fonctecgni, fonctecgui, fonctecce, fonctecci, fonctvci)        
       ELSE
          CALL passQLints( p, nu, omega, fonctce, fonctci, fonctcgte, fonctcgti, fonctcgne, fonctcgni, &
               & fonctcce, fonctcci, fonctece, foncteci, fonctecgte, fonctecgti, fonctecgne, fonctecgni, fonctecce, fonctecci) 
          fonctcgui(:) = 0.
          fonctvci(:) = 0.
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
       fonctpgui(:) = 0.
       fonctpce = 0.
       fonctpci(:) = 0. 

       fonctepe = 0.
       fonctepi(:) = 0. 

       fonctvpi(:) = 0.

       fonctepgte = 0.
       fonctepgti(:) = 0. 
       fonctepgne = 0.
       fonctepgni(:) = 0. 
       fonctepgui(:) = 0. 
       fonctepce = 0.
       fonctepci(:) = 0. 

    ELSE      

       IF ( (rotflagarray(p) == 1) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) THEN
          CALL trapQLintsrot( p, nu, omega, fonctpe, fonctpi, fonctpgte, fonctpgti, fonctpgne, fonctpgni, fonctpgui, &
               & fonctpce, fonctpci, fonctepe, fonctepi, fonctepgte, fonctepgti, fonctepgne, fonctepgni, fonctepgui, fonctepce, fonctepci, fonctvpi)
       ELSE
          CALL trapQLints( p, nu, omega, fonctpe, fonctpi, fonctpgte, fonctpgti, fonctpgne, fonctpgni, &
               & fonctpce, fonctpci, fonctepe, fonctepi, fonctepgte, fonctepgti, fonctepgne, fonctepgni,fonctepce, fonctepci)
          fonctpgui(:) = 0.
          fonctvpi(:) = 0.
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
    fonxcirci(:) = fonctci(:)*ninorm(p,:) !renormalise coefi
    fonxpiegi(:) = fonctpi(:)*ninorm(p,:) !renormalise coefi

    !To save memory space, this condition is in place
    IF (phys_meth .NE. 0.0) THEN
       fonxcircgte = fonctcgte
       fonxpieggte = fonctpgte
       fonxcircgti(:) = fonctcgti(:)*ninorm(p,:) !renormalise coefi
       fonxpieggti(:) = fonctpgti(:)*ninorm(p,:) !renormalise coefi
       fonxcircgne = fonctcgne
       fonxpieggne = fonctpgne
       fonxcircgni(:) = fonctcgni(:)*ninorm(p,:) !renormalise coefi
       fonxpieggni(:) = fonctpgni(:)*ninorm(p,:) !renormalise coefi
       fonxcircgui(:) = fonctcgui(:)*ninorm(p,:) !renormalise coefi
       fonxpieggui(:) = fonctpgui(:)*ninorm(p,:) !renormalise coefi
       fonxcircce = fonctcce
       fonxpiegce = fonctpce
       fonxcircci(:) = fonctcci(:)*ninorm(p,:) !renormalise coefi
       fonxpiegci(:) = fonctpci(:)*ninorm(p,:) !renormalise coefi
       !!
       IF (phys_meth == 2) THEN
          fonxecircgte = fonctecgte
          fonxepieggte = fonctepgte
          fonxecircgti(:) = fonctecgti(:)*ninorm(p,:) !renormalise coefi
          fonxepieggti(:) = fonctepgti(:)*ninorm(p,:) !renormalise coefi
          fonxecircgne = fonctecgne
          fonxepieggne = fonctepgne
          fonxecircgni(:) = fonctecgni(:)*ninorm(p,:) !renormalise coefi
          fonxepieggni(:) = fonctepgni(:)*ninorm(p,:) !renormalise coefi
          fonxecircgui(:) = fonctecgui(:)*ninorm(p,:) !renormalise coefi
          fonxepieggui(:) = fonctepgui(:)*ninorm(p,:) !renormalise coefi
          fonxecircce = fonctecce
          fonxepiegce = fonctepce
          fonxecircci(:) = fonctecci(:)*ninorm(p,:) !renormalise coefi
          fonxepiegci(:) = fonctepci(:)*ninorm(p,:) !renormalise coefi
       ENDIF
    ENDIF

    fonxecirce = fonctece
    fonxepiege = fonctepe
    fonxecirci(:) = foncteci(:)*ninorm(p,:) !renormalise coefi
    fonxepiegi(:) = fonctepi(:)*ninorm(p,:) !renormalise coefi

    fonxvcirci(:) = fonctvci(:)*ninorm(p,:) !renormalise coefi
    fonxvpiegi(:) = fonctvpi(:)*ninorm(p,:) !renormalise coefi

  END SUBROUTINE make_QLflux

END MODULE QLflux
 
