!===============================================================================
real*8 function ceis(mc)
!
! Double precision minimax rational approximation of S(m):
! the special complete elliptic integral
!
!     Reference: T. Fukushima, (2016), Astron. J., re-revised
!
!      "Zonal Toroidal Harmonic Expansions of External Gravitational
!       Fields for Ring-Like Objects"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
real*8 mc
real*8 t,t2
if(mc.gt.1.d0) then
    write(*,"(a20,1pe15.7)") "(ceis) negative m: mc=",mc
elseif(mc.gt.0.601505d0) then
    t=2.5094417746772230517d0*mc-1.5094417746772230517d0
    t2=t*t
    ceis=((1026.17763535952307d0 &
    +t2*(485.677469611677641d0 &
    -t2*1.03090351230372971d0)) &
    +t*(1274.82870240508851d0 &
    +t2*(52.3124844458365222d0 &
    +t2*0.00739324234515133022d0)))/ &
    ((3634.91161208781745d0 &
    +t2*(3687.92980597191086d0 &
    +t2*72.4741783934871057d0)) &
    +t*(6150.43453804219490d0 &
    +t2*(908.926554921646573d0 &
    -t2)))
elseif(mc.gt.0.399058d0) then
    t=4.9395644292086323828d0*mc-1.9711727019911384214d0
    t2=t*t
    ceis=((86.7180733051785905d0 &
    +t2*(-36.6337761424801254d0 &
    -t2*1.38891143351659923d0)) &
    +t*(15.7273307842061126d0 &
    +t2*(-14.8291225844217536d0 &
    +t2*0.00361918439375104252d0)))/ &
    ((235.196742248774012d0 &
    +t2*(-88.9257271644547379d0 &
    -t2*15.3339604566201331d0)) &
    +t*(116.962483360689959d0 &
    +t2*(-71.2170053927684721d0 &
    -t2)))
elseif(mc.gt.0.219501d0) then
    t=5.569262128460600255d0*mc-1.2224586064592302166d0
    t2=t*t
    ceis=((74.7877514681889701d0 &
    +t2*(82.9949360943443884d0 &
    +t2*1.56773125483764527d0)) &
    +t*(134.483076192401992d0 &
    +t2*(20.2462516310088192d0 &
    -t2*0.00367147603407009175d0)))/ &
    ((143.112254337953479d0 &
    +t2*(267.680587804708068d0 &
    +t2*17.1824818506270368d0)) &
    +t*(321.171975806674627d0 &
    +t2*(101.689772698952843d0 &
    +t2)))
elseif(mc.gt.0.128079d0) then
    t=11.059255490920351242d0*mc-1.4275176395125080180d0
    t2=t*t
    ceis=((368.838713207397348d0 &
    +t2*(273.591392228690348d0 &
    +t2*3.11283723054624760d0)) &
    +t*(546.423411840019851d0 &
    +t2*(52.8283480246527158d0 &
    -t2*0.00267666613712733130d0)))/ &
    ((537.483931788056883d0 &
    +t2*(649.569549167506858d0 &
    +t2*24.2918898835263356d0)) &
    +t*(977.761239499760638d0 &
    +t2*(191.900997942387173d0 &
    +t2)))
elseif(mc.gt.0.071412d0) then
    t=19.090527280363483639d0*mc-1.4641861708220381047d0
    t2=t*t
    ceis=((737.51190337498540d0 &
    +t2*(509.23908689777924d0 &
    +t2*5.2334058551069839d0)) &
    +t*(1056.48245489063616d0 &
    +t2*(93.951657090929254d0 &
    +t2*0.00127302402727680993d0)))/ &
    ((848.68441046566805d0 &
    +t2*(911.51082544227172d0 &
    +t2*28.2934744382382281d0)) &
    +t*(1462.80391382404865d0 &
    +t2*(248.592677538715336d0 &
    +t2)))
elseif(mc.gt.0.045739d0) then
    t=32.301828283480845016d0*mc-1.4774533238581303702d0
    t2=t*t
    ceis=((1228.30248538045050d0 &
    +t2*(830.59901244101972d0 &
    +t2*8.3445760764027038d0)) &
    +t*(1742.14135275252047d0 &
    +t2*(151.363543340919925d0 &
    +t2*0.0109389121873063296d0)))/ &
    ((1148.72064907022826d0 &
    +t2*(1154.29906916591534d0 &
    +t2*31.7183859972304475d0)) &
    +t*(1922.56327226047254d0 &
    +t2*(299.505582648248342d0 &
    +t2)))
elseif(mc.gt.0.027299d0) then
    t=54.22993492407809111d0*mc-1.4804229934924078091d0
    t2=t*t
    ceis=((1850.32817333933852d0 &
    +t2*(1251.96912463087426d0 &
    +t2*12.7282329449538449d0)) &
    +t*(2624.51676940821043d0 &
    +t2*(228.600812362065591d0 &
    +t2*0.0301997195942178328d0)))/ &
    ((1440.47401213691609d0 &
    +t2*(1388.46074898128627d0 &
    +t2*34.8987556516293270d0)) &
    +t*(2368.67335584121873d0 &
    +t2*(347.987916163266950d0 &
    +t2)))
elseif(mc.gt.0.016280d0) then
    t=90.75233687267447137d0*mc-1.4774480442871403939d0
    t2=t*t
    ceis=((2594.88367393128139d0 &
    +t2*(1777.35557927706750d0 &
    +t2*18.6194756219725130d0)) &
    +t*(3700.86979543655025d0 &
    +t2*(327.656946454522957d0 &
    +t2*0.063865181191751846d0)))/ &
    ((1716.03161920071730d0 &
    +t2*(1612.19856579654867d0 &
    +t2*37.9050752099373943d0)) &
    +t*(2793.09781104785086d0 &
    +t2*(394.278792974606013d0 &
    +t2)))
elseif(mc.gt.0.009692d0) then
    t=151.79113539769277474d0*mc-1.4711596842744383728d0
    t2=t*t
    ceis=((3441.67541433275474d0 &
    +t2*(2401.49730335565496d0 &
    +t2*26.1675237391586565d0)) &
    +t*(4950.25801126609305d0 &
    +t2*(448.806533020458946d0 &
    +t2*0.117111370758512406d0)))/ &
    ((1967.03500815730011d0 &
    +t2*(1820.81267345256858d0 &
    +t2*40.7252758168263498d0)) &
    +t*(3184.72957744623299d0 &
    +t2*(437.723184552594083d0 &
    +t2)))
elseif(mc.gt.0.d0) then
    t=1.d0-103.17787866281469253d0*mc
    ceis=-log(mc*0.0625d0) &
    *(1.94180589378811075d6 &
    +t*(647.22748042899378d0 &
    +t*(-10.1649835519921633d0 &
    +t*0.00183293824579460819d0)))/ &
    (3.79948492107705041d6 &
    +t*(84831.477117968440d0 &
    +t*(568.51804083313882d0 &
    +t))) &
    -(21159.8549313831840d0 &
    +t*(26.4758907608130305d0 &
    +t*0.00284022078125062437d0))/ &
    (10388.4217213415863d0 &
    +t*(203.745109840801902d0 &
    +t))
else
    write(*,"(a20,1pe15.7)") "(ceis) out of range: mc=",mc
endif
!
!write(*,"(1x,a20,0p5f20.15)") "(ceis) mc,K=",mc,ceis
!
return
end