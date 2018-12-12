subroutine rcelbd(mc,elb,eld)
!
!   Single precision general complete elliptic integrals of the second kind
!
!   Reference: T. Fukushima, (2011) Math. Comp., 80, 1725-1743
!      "Precise and Fast Computation of General Complete Elliptic Integral
!       of Second Kind"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
!     Inputs: mc   = complementary parameter 0 <= mc   <= 1
!
!     Output: elb,eld
!
real*4 mc,elk,elb,ele,eld
real*4 m,mx,kkc,nome,eec,kec

real*4 PIQ,PIHALF,PI,PIINV
parameter (PIQ=0.78539816)
parameter (PIHALF=1.57079633)
parameter (PI=3.14159265)
parameter (PIINV=0.318309886)

real*4 mcold,elbold,eldold
save mcold,elbold,eldold

real*4 Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8
parameter (Q1=1.0/16.0,Q2=1.0/32.0,Q3=21.0/1024.0)
parameter (Q4=31.0/2048.0,Q5=6257.0/524288.0)
parameter (Q6=10293.0/1048576.0,Q7=279025.0/33554432.0)
parameter (Q8=483127.0/67108864.0)

real*4 K1,K2,K3,K4
parameter (K1=1.0/4.0)
parameter (K2=9.0/64.0)
parameter (K3=25.0/256.0)
parameter (K4=1225.0/16384.0)

real*4 B1,B2,B3,B4,B5
parameter (B1=1.0/2.0)
parameter (B2=1.0/16.0)
parameter (B3=3.0/128.0)
parameter (B4=25.0/2048.0)
parameter (B5=245.0/32768.0)

real*4 D1,D2,D3,D4,D5
parameter (D1=1.0/2.0)
parameter (D2=3.0/16.0)
parameter (D3=15.0/128.0)
parameter (D4=175.0/2048.0)
parameter (D5=2205.0/32768.0)

real*4 logq2,dkkc,dddc,ddc,dele,delb,elk1

logical first/.TRUE./

if(first) then
    first=.FALSE.
    mcold=1.0
    elbold=PIQ
    eldold=PIQ
endif
m=1.0-mc
if(abs(mc-mcold).lt.1.19e-7*mc) then
    elb=elbold
    eld=eldold
elseif(m.lt.1.19e-7) then
    elb=PIQ
    eld=PIQ
elseif(mc.lt.1.19e-7) then
    elb=1.0
    eld=0.386294361-0.5*log(mc)
elseif(mc.lt.0.1) then
    nome=mc*(Q1+mc*(Q2+mc*(Q3+mc*(Q4+mc*(Q5+mc*(Q6+mc*(Q7+mc*Q8))))))) 
    if(mc.lt.0.01) then
        dkkc=mc*(K1+mc*(K2+mc*(K3+mc*K4)))
        dddc=mc*(D1+mc*(D2+mc*(D3+mc*D4)))
    else
        mx=mc-0.05

! (K'-1)/(pi/2)

        dkkc=0.0128642566+mx*(0.2648342989+mx*(0.1564757379+mx*( &
             0.1142614608+mx*(0.0920272442+mx*(0.0784321883+mx*( &
             0.0693526014))))))

! (K'-E')/(pi/2)

        dddc=0.0254839544+mx*(0.5196738432+mx*(0.2064495111+mx*( &
             0.1361095213+mx*(0.1045801404+mx*(0.0867461292+mx*( &
             0.0753638027))))))
    endif
    kkc=1.0+dkkc
    logq2=-0.5*log(nome)
    elk=kkc*logq2
    dele=-dkkc/kkc+logq2*dddc
    elk1=elk-1.0
    delb=(dele-mc*elk1)/m
    elb=1.0+delb
    eld=elk1-delb
elseif(m.le.0.01) then
    elb=PIHALF*(B1+m*(B2+m*(B3+m*(B4+m*B5))))
    eld=PIHALF*(D1+m*(D2+m*(D3+m*(D4+m*D5))))
elseif(m.le.0.1) then
    mx=0.95-mc
    elb=0.7904014136+mx*(0.1020062662+mx*(0.03987839556+mx*( &
        0.02173713638+mx*(0.01396097977+mx*(0.009892518823)))))
    eld=0.8006020402+mx*(0.3139944778+mx*(0.2059131187+mx*( &
        0.1577443465+mx*(0.1305950773+mx*(0.1133084745)))))
elseif(m.le.0.2) then
    mx=0.85-mc
    elb=0.8010240645+mx*(0.1106953445+mx*(0.0473487467+mx*( &
        0.0284843673+mx*(0.0202778114+mx*(0.0159650059)))))
    eld=0.8342326678+mx*(0.3604952816+mx*(0.2623796641+mx*( &
        0.2237239445+mx*(0.2064478118+mx*(0.1998094409)))))
elseif(m.le.0.3) then
    mx=0.75-mc
    elb=0.8125977729+mx*(0.1211096179+mx*(0.05729337683+mx*( &
        0.03850945160+mx*(0.03078343030+mx*(0.02729056493+mx*( &
        0.02591636929))))))
    eld=0.8731525819+mx*(0.4206222307+mx*(0.3442310616+mx*( &
        0.3311330218+mx*(0.3452772851+mx*(0.3779453222+mx*( &
        0.4273780125))))))
elseif(m.le.0.4) then
    mx=0.65-mc
    elb=0.8253235580+mx*(0.1338621161+mx*(0.07101129360+mx*( &
        0.05417847742+mx*(0.04945174495+mx*(0.05022219622+mx*( &
        0.05474291317))))))
    eld=0.9190270392+mx*(0.5010021593+mx*(0.4688312706+mx*( &
        0.5177142278+mx*(0.6208433913+mx*(0.7823643938+mx*( &
        1.0191145351))))))
elseif(m.le.0.5) then
    mx=0.55-mc
    elb=0.8394795703+mx*(0.1499164403+mx*(0.09083193582+mx*( &
        0.08034703348+mx*(0.08563844050+mx*(0.1019547259+mx*( &
        0.1305748115))))))
    eld=0.9744043665+mx*(0.6132468054+mx*(0.6710966695+mx*( &
        0.8707276202+mx*(1.2295422312+mx*(1.8266059675+mx*( &
        2.8069345310+mx*(4.4187893291)))))))
elseif(m.le.0.6) then
    mx=0.45-mc
    elb=0.8554696152+mx*(0.1708960727+mx*(0.1213352290+mx*( &
        0.1282018836+mx*(0.1646872815+mx*(0.2374189087+mx*( &
        0.3692081047))))))
    eld=1.0434552951+mx*(0.7796257219+mx*(1.0297423609+mx*( &
        1.6220372234+mx*(2.7879895312+mx*(5.0483814874+mx*( &
        9.4632776119+mx*(18.181489949+mx*(35.580980591))))))))
elseif(m.le.0.7) then
    mx=0.35-mc
    elb=0.8739200618+mx*(0.1998140575+mx*(0.1727696159+mx*( &
        0.2281069133+mx*(0.3704681411+mx*(0.6792712529+mx*( &
        1.3480084967+mx*(2.8276709769)))))))
    eld=1.1336783366+mx*(1.0486431737+mx*(1.7534650412+mx*( &
        3.5231827268+mx*(7.7494764138+mx*(17.986450056+mx*( &
        43.255916346+mx*(106.68153445+mx*(268.09848657))))))))
elseif(m.le.0.8) then
    mx=0.25-mc
    elb=0.8959028209+mx*(0.2431400038+mx*(0.2730818756+mx*( &
        0.4862800075+mx*(1.0827474372+mx*(2.7434452910+mx*( &
        7.5558178287+mx*(22.051940825+mx*(67.156406447+mx*( &
        211.27225379)))))))))
    eld=1.2606128266+mx*(1.5486656381+mx*(3.5536694119+mx*( &
        9.9004446761+mx*(30.320566617+mx*(98.180258659+mx*( &
        329.77101043+mx*(1136.6559897+mx*(3993.8343357+mx*( &
        14242.729587+mx*(51394.757292))))))))))
elseif(m.le.0.85) then
    mx=0.175-mc
    elb=0.9159220526+mx*(0.2947142524+mx*(0.4357767093+mx*( &
        1.0673282465+mx*(3.3278441186+mx*(11.904060044+mx*( &
        46.478388202+mx*(192.75560026)))))))
    eld=1.4022005691+mx*(2.3222058979+mx*(7.4621583665+mx*( &
        29.435068908+mx*(128.15909243+mx*(591.08070369+mx*( &
        2830.5462296+mx*(13917.764319+mx*(69786.105252))))))))
else
    mx=0.125-mc
    elb=0.9319060610+mx*(0.3484480295+mx*(0.6668091788+mx*( &
        2.2107691357+mx*(9.4917650489+mx*(47.093047910+mx*( &
        255.92004602+mx*(1480.0295327+mx*(8954.0409047+mx*( &
        56052.482210)))))))))
    eld=1.5416901127+mx*(3.3791762146+mx*(14.940583857+mx*( &
        81.917739292+mx*(497.49005466+mx*(3205.1840102+mx*( &
        21457.322374+mx*(147557.01566+mx*(1.0350452902e6+mx*( &
        7.3719223348e6+mx*(5.3143443951e7))))))))))
endif

mcold=mc
elbold=elb
eldold=eld

return
end
!---------------------------------------------------------------------------
subroutine relcbd(c0,mc,b,dx)

real*4 c0,mc,b,dx
real*4 c,x,y,s,m,d,sy
real*4 yy(11),ss(11)
integer j,i

c=c0
x=c*c
y=1.0-x
s=sqrt(y)
if(x.gt.0.1) then
    call relsbd(s,mc,b,dx)
    return
endif
m=1.0-mc
ss(1)=s
do j=1,10
    d=sqrt(mc+m*x)
    x=(c+d)/(1.0+d)
    y=1.0-x
    yy(j+1)=y
    ss(j+1)=sqrt(y)
    if(x.gt.0.1) then
        goto 1
    endif
    c=sqrt(x)
enddo
write(*,*) "(relcbd) too many iterations: c0,mc=",c0,mc
1 continue
s=ss(j+1)
call relsbd(s,mc,b,dx)
do i=1,j
    sy=ss(j-i+1)*yy(j-i+2)
    b=b+(b-sy)
    dx=dx+(dx+sy)
enddo
return
end
!---------------------------------------------------------------------------
subroutine relsbd(s0,mc,b,d)

real*4 s0,mc,b,d
real*4 m,del,s,y,sy
real*4 yy(11),ss(11)
integer j

m=1.0-mc
del=0.1888-0.0378*m ! F6    Optimum
s=s0
y=s*s
if(y.lt.del) then
    call rserbd(y,m,b,d)
    b=s*b
    d=s*y*d
    return
endif
ss(1)=s
do j=1,10
    y=y/((1.0+sqrt(1.0-y))*(1.0+sqrt(1.0-m*y)))
    yy(j+1)=y
    ss(j+1)=sqrt(y)
    if(y.lt.del) then
        goto 1
    endif
enddo
write(*,*) "(relsbd) too many iterations: s0,m=",s0,m
1 continue
call rserbd(y,m,b,d)
b=ss(j+1)*b
d=ss(j+1)*y*d
do i=1,j
    sy=ss(j-i+1)*yy(j-i+2)
    b=b+(b-sy)
    d=d+(d+sy)
enddo
return
end
!---------------------------------------------------------------------------
subroutine rserbd(y,m,b,d)

real*4 y,m,b,d
real*4 F1,F2,F3,F4,F10,F20,F21,F30,F31,F40,F41,F42
real*4 F5,F50,F51,F52,F6,F60,F61,F62,F63
parameter (F10=1.0/6.0)
parameter (F20=3.0/40.0)
parameter (F21=2.0/40.0)
parameter (F30=5.0/112.0)
parameter (F31=3.0/112.0)
parameter (F40=35.0/1152.0)
parameter (F41=20.0/1152.0)
parameter (F42=18.0/1152.0)
parameter (F50=63.0/2816.0)
parameter (F51=35.0/2816.0)
parameter (F52=30.0/2816.0)
parameter (F60=231.0/13312.0)
parameter (F61=126.0/13312.0)
parameter (F62=105.0/13312.0)
parameter (F63=100.0/13312.0)

real*8 A1,A2,A3,A4,A5,A6
parameter (A1=3.0/5.0)
parameter (A2=5.0/7.0)
parameter (A3=7.0/9.0)
parameter (A4=9.0/11.0)
parameter (A5=11.0/13.0)
parameter (A6=13.0/15.0)

real*8 B1,B2,B3,B4,B5,B6
real*8 D0,D1,D2,D3,D4,D5,D6
parameter (D0=1.0/3.0)

F1=F10+m*F10
F2=F20+m*(F21+m*F20)
F3=F30+m*(F31+m*(F31+m*F30))
F4=F40+m*(F41+m*(F42+m*(F41+m*F40)))
F5=F50+m*(F51+m*(F52+m*(F52+m*(F51+m*F50))))
F6=F60+m*(F61+m*(F62+m*(F63+m*(F62+m*(F61+m*F60)))))

D1=F1*A1
D2=F2*A2
D3=F3*A3
D4=F4*A4
D5=F5*A5
D6=F6*A6

d=D0+y*(D1+y*(D2+y*(D3+y*(D4+y*(D5+y*D6)))))

B1=F1-D0
B2=F2-D1
B3=F3-D2
B4=F4-D3
B5=F5-D4
B6=F6-D5

b=1.0+y*(B1+y*(B2+y*(B3+y*(B4+y*(B5+y*B6)))))

return
end
