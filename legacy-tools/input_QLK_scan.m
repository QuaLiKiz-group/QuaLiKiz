%Creates binary input files based on the content of parameters.m

curdir=pwd;
ind=strfind(curdir,'/'); dir_path=curdir(ind(end)+1:end);

eval(sprintf('%s','!mkdir input'))
eval(sprintf('%s','!mkdir output'))
eval(sprintf('%s','!mkdir output/primitive'))
eval(sprintf('%s',['!ln -s ../../src/QuaLiKiz.exe .']))
eval(sprintf('%s',['!ln -s ../../src/qlk_makeflux.exe .']))

runname=input('Please enter string for batch job name (e.g. ''jsmith_run10''). Hit enter for default of first 8 characters of directory path: ');
if isempty(runname) == 1
	if length(dir_path) > 8
		runname=dir_path(1:8);
	else
		runname=dir_path(1:end);
	end
end

parameters_template; %will read parameter file previously placed by user in newly created run directory

WriteQLKBatchPBS(runname,nprocs) %for Cephee
%Calculate auxilliary quantities for plotting
Epsilonx=Rmin.*x./Ro;
ft=2.*(2.*Epsilonx).^(0.5)./pi; %trapped particle fraction
tau=Tex./Tix(:,1);

%%Quasineutrality (including gradients) test and Zeff calculation
ninormd=ninorm;
ninormd(ion_type==3)=0; %zero out tracers for QN check
ninormd(ion_type==4)=0; %zero out tracers for QN check

if (set_ninorm1 == 1) && (nions > 1) %sets ninorm of 1st species to maintian quasineutrality
   if nions > 2
      ninorm(:,1) = (1 - sum(ninormd(:,2:end)'.*Zi(:,2:end)')')./Zi(:,1);
   else
      ninorm(:,1) = (1 - ninormd(:,2).*Zi(:,2))./Zi(:,1);
   end
end
if (set_Ani1 == 1 && nions > 1) %sets Ani of 1st species to maintian quasineutrality
   if nions > 2
       Ani(:,1) = (Ane - sum(ninormd(:,2:end)'.*Ani(:,2:end)'.*Zi(:,2:end)')')./(Zi(:,1).*ninorm(:,1));
   else
       Ani(:,1) = (Ane - ninormd(:,2).*Ani(:,2).*Zi(:,2))./(Zi(:,1).*ninorm(:,1));
   end
end

quasicheck=-1.*ones(scann,1);
quasicheck_grad=-Ane;
quasitol=1e-5;

for i = 1:scann
  for j = 1:nions
    if ion_type(i,j) < 3
      quasicheck(i)=quasicheck(i) + Zi(i,j).*ninorm(i,j);
      quasicheck_grad(i)=quasicheck_grad(i) + Zi(i,j).*ninorm(i,j).*Ani(i,j);
    end
  end
end

noQNgrad=zeros(scann,1); %set flag for showing QN in gradients failure
for i = 1:scann
  if abs(quasicheck(i))  > quasitol
    error(['Quasineutrality not respected by an absolute value ',num2str(abs(quasicheck(i))),'. Tolerance is ',num2str(quasitol),'. Scan index = ',num2str(i)]);
  end  
  if abs(quasicheck_grad(i)) > quasitol && (set_QN_grad==1)
    noQNgrad(i)=1;
  end  
end
if find(noQNgrad==1)
  error(['Quasineutrality of gradients not respected by an absolute value ',num2str(abs(quasicheck_grad(find(noQNgrad==1))')),'. Tolerance is ',num2str(quasitol),'. Scan indices = ',num2str(find(noQNgrad==1)')]);
end

ninormz=ninorm;
ninormz(ion_type==3)=0; %zero out tracers for Zeff calc
for i=1:scann  	   
   Zeffx(i)=sum(ninormz(i,:).*Zi(i,:).^2); %calculate Zeff
end

figure; 
disp('Figures displayed such that the input data can be verified by eye') 
subplot(221)
set(gca,'FontSize',18)
plot(1:scann,Ati(:,1),'r',1:scann,Tix(:,1),'r--',1:scann,Ate,'b',1:scann,Tex,'b--','LineWidth',2)
xlabel('Scan index')

l1=legend('-R\nabla{T_i}/T_i','Ti','-R\nabla{T_e}/T_e','T_e');
set(l1,'FontSize',12)
legend boxoff
grid on
subplot(222)
set(gca,'FontSize',18)
plot(1:scann,Ane,'g',1:scann,Nex,'g--','LineWidth',2)
xlabel('Scan index')

l2=legend('-R\nabla{n_e}/n_e','n_e');
set(l2,'FontSize',12)
legend boxoff
grid on
subplot(223)
set(gca,'FontSize',18)
plot(1:scann,tau,'c',1:scann,Zeffx,'c--',1:scann,ft,'c-.','LineWidth',2)
xlabel('Scan index')

l3=legend('T_e/T_i','Z_{eff}','f_t');
set(l3,'FontSize',12)
legend boxoff
grid on
subplot(224)
set(gca,'FontSize',18)
plot(1:scann,smag,'m',1:scann,qx,'m--',1:scann,alphax,'m-.s','LineWidth',2)
xlabel('Scan index')
l4=legend('s','q','\alpha');
set(l4,'FontSize',12)
legend boxoff
grid on

%%% End

if(1)

%%% Saving binary files
%% ---------------------------------
%% MANAGING QUALIKIZ INPUT VARIABLES
%% ---------------------------------
kc = 1;
p{kc} = length(x);                  kc=kc+1;%p{1} (-) Number of radial or scan points
p{kc} = length(kthetarhos);         kc=kc+1;%p{2} (-) Number of wavenumbers
p{kc} = nions;                      kc=kc+1;%p{3} (-) Number of ions in system

%Flag input and metadata
p{kc} = phys_meth;	            kc=kc+1;%p{4} (-) Flag for additional calculation (default 0.0)
p{kc} = coll_flag;	            kc=kc+1;%p{5} (-) Flag for collisionality (default 0.0)
p{kc} = rot_flag;	            kc=kc+1;%p{6} (-) Flag  for rotation (default 0.0)
p{kc} = verbose;	            kc=kc+1;%p{7} (-) Flag  for setting level of output verbosity
p{kc} = separateflux;	            kc=kc+1;%p{8} (-) Flag  for separate outputs of ETG, TEM, ITG fluxes
p{kc} = numsols;	            kc=kc+1;%p{9} (-) Number of solutions requested
p{kc} = relacc1;	 	    kc=kc+1;%p{10} (-) Relative accuracy in 1D integrals
p{kc} = relacc2;	 	    kc=kc+1;%p{11} (-) Relative accuracy in 2D integrals
p{kc} = maxruns; 	            kc=kc+1;%p{12} (-) Number of runs jumping directly to Newton between contour checks
p{kc} = maxpts;                     kc=kc+1;%p{13}(-) Number of integrand evaluations done in 2D integral
p{kc} = timeout;                    kc=kc+1;%p{14}(-) Upper time limit (s) beyond which solutions are not sought after at a given wavenumber and radius
p{kc} = ETGmult;                    kc=kc+1;%p{15}(-) ETG multiplier (for testing)
p{kc} = collmult;                   kc=kc+1;%p{16}(-) collisionality multiplier (for testing)

%Geometry input
p{kc} = R0;		            kc=kc+1;%p{17} (m) Geometric major radius used for normalizations
p{kc} = kthetarhos;                 kc=kc+1;%p{18} (-) Wave spectrum input: Vector (dimn)
p{kc} = x;                          kc=kc+1;%p{19} (-) radial normalised coordinate (midplane average)
p{kc} = rho;                        kc=kc+1;%p{20} (-) normalized toroidal flux coordinate
p{kc} = Ro;		            kc=kc+1;%p{21} (m) Major radius. Radial profile due to Shafranov shift
p{kc} = Rmin;	                    kc=kc+1;%p{22} (m) Geometric minor radius. Assumed to be a midplane average at LCFS. Currently a profile but should probably be shifted to a scalar
p{kc} = Bo;		            kc=kc+1;%p{23} (T) Likely not very rigorous to use this sqrt(<B²>) for calculating the Larmor radius % quite close to <Bphi> in practice however 
p{kc} = qx;               	    kc=kc+1;%p{24} (-) Vector (radial grid x(aa))
p{kc} = smag;            	    kc=kc+1;%p{25} (-) Vector (radial grid x(aa))  q is a flux surface quantity --> makes sense to consider s = rho/q dq/drho
p{kc} = alphax;            	    kc=kc+1;%p{26} (-) Vector (radial grid x(aa)) 

%Rotation input
p{kc} = Machtor;            	    kc=kc+1;%p{27} (-) Vector (radial grid x(aa)) 
p{kc} = Autor;            	    kc=kc+1;%p{28} (-) Vector (radial grid x(aa)) 
p{kc} = Machpar;            	    kc=kc+1;%p{29} (-) Vector (radial grid x(aa)) 
p{kc} = Aupar;            	    kc=kc+1;%p{30} (-) Vector (radial grid x(aa)) 
p{kc} = gammaE;            	    kc=kc+1;%p{31} (-) Vector (radial grid x(aa))

%Electron input
p{kc} = Tex;       		    kc=kc+1;%p{32} (keV) Vector (radial grid x(aa))
p{kc} = Nex;    		    kc=kc+1;%p{33} (10^19 m^-3) Vector (radial grid x(aa))
p{kc} = Ate;                        kc=kc+1;%p{34} (-) Vector (radial grid x(aa))
p{kc} = Ane;                	    kc=kc+1;%p{35} (-) Vector (radial grid x(aa))
p{kc} = el_type;                    kc=kc+1;%p{36} Kinetic or adiabatic
p{kc} = anise;                      kc=kc+1;%p{37}  Tperp/Tpar at LFS
p{kc} = danisedr;                   kc=kc+1;%p{38}  d/dr(Tperp/Tpar) at LFS

%Ion inputs (can be for multiple species)
p{kc} = Tix;                kc=kc+1;%p{39} (keV) Vector (radial grid x(aa))
p{kc} = ninorm;             kc=kc+1;%p{40} ni/ne Vector (radial grid x(aa))
p{kc} = Ati;      	    kc=kc+1;%p{41}  (-) Vector (radial grid x(aa))
p{kc} = Ani;                kc=kc+1;%p{42}  (-) Vector (radial grid x(aa))  check calculation w.r.t. Qualikiz electroneutrality assumption
p{kc} = ion_type;           kc=kc+1;%p{43}  Kinetic, adiabatic, tracer
p{kc} = anis;               kc=kc+1;%p{44}  Tperp/Tpar at LFS
p{kc} = danisdr;            kc=kc+1;%p{45}  d/dr(Tperp/Tpar) at LFS
p{kc} = Ai;	            kc=kc+1;%p{46} (-) Ion mass
p{kc} = Zi;     	    kc=kc+1;%p{47} (-) Ion charge

nargu = kc-1;
stringind=[]; %no string indexes. Kept if adding any future string inputs

%Write binary files
for i=1:nargu
  if any(i == stringind)
    eval(sprintf('fid=fopen(\047input/p%d.txt\047,\047wb\047);',i));
    eval(sprintf('fwrite(fid,p{%d});',i));
    fclose(fid);
  else  
    eval(sprintf('fid=fopen(\047input/p%d.bin\047,\047wb\047);',i));
    eval(sprintf('fwrite(fid,p{%d},\047double\047);',i));
    fclose(fid);
  end
end

end


