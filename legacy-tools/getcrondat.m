%THIS ROUTINE EXTRACTS CRONOS DATA FROM CRONOS DATA AND PARAM INPUT FILE. TIME WINDOW FOR AVERAGING IS PROVIDED BY USER. SINGLE TIME SLICES ARE ALSO
%ACCEPTABLE. NOTE, CRONOS MUST BE IN THE BACKGROUND FOR SOME OF THE USED FUNCTIONS TO WORK! OUTPUT DATA IS IN A FORM RELEVANT FOR INPUT TO QUALIKIZ

%NOTE, CAN ALSO BE USED IN A STAND ALONE MANNER FOR ANALYSIS OF CRONOS SIMULATIONS

%nshot is used as a serial number for the data file

function crondat=getcrondat(data,param,nshot)

qe=1.6e-19;
me=9.1e-31;
mp=1.67e-27;
muo=4*pi*1e-7;
Ze = -1;

tps = input(['Please input a time window for averaging (in s), e.g: [5.25 5.75]. For this shot, tmin = ', num2str(min(data.gene.temps)), ' and tmax=',num2str(max(data.gene.temps)) , ': ']);

x=param.gene.x;
dimxr=length(x);
itps=iround(data.gene.temps,tps);

Bo = avg(data.geo.b0(itps(1):itps(2))).*ones(1,dimxr);% Magnetic field array
Ro = avg(data.geo.r0(itps(1):itps(2))).*ones(1,dimxr);% Geometric major radius array

%manipulation of CRONOS R and Z matrices to calculate outboard, midplane average, and flux surface average minor radii
%note that the toroidal flux coordinate axis on R and Z are defined on the data.equi.rhoRZ grid, and are thus interpolated back to the constant
%spaced x grid
R=double(squeeze(avg(data.equi.R(itps(1):itps(2),:,:))));
Z=double(squeeze(avg(data.equi.Z(itps(1):itps(2),:,:))));

xrho=double(avg(data.equi.rhoRZ(itps(1):itps(2),:))); xrho=xrho./xrho(end);

dimtheta=length(R(1,:));

%define arc length on the flux surface which corresponds to each coordinate. Definition takes average of length from each coordinate to 
%its 2 neighbours. This arc length is then used to define the flux surface average minor radius, by weighting each minor radius of each
%coordinate with the arc length when defining the average

ds=(sqrt((circshift(R',1)-R').^2+(circshift(Z',1)-Z').^2)'+sqrt((circshift(R',dimtheta-1)-R').^2+(circshift(Z',dimtheta-1)-Z').^2)')/2;
ds(1,:)=eps; %arc length on axis defined as a small number to avoid division by 0

a=sqrt((R-R(1,1)).^2+(Z-Z(1,1)).^2);

aout=max(R')-R(1,1); %outboard side minor radius
amid=(max(R')-min(R'))/2; %midplane average minor radius

aavg=trapz(linspace(0,1,dimtheta),ds'.*a')./trapz(linspace(0,1,dimtheta),ds');  %flux surface average minor radius

%finally, interpolate the minor radii onto the regular spaced x grid
aout=interp1(xrho,aout,x);
amid=interp1(xrho,amid,x);
aavg=interp1(xrho,aavg,x);

Rmin = aavg(end).*ones(1,dimxr);% <a> array (minor radius)

%calculate profiles
qx=avg(data.prof.q(itps(1):itps(2),:));% q array
smag=x.*gradient(qx,x)./qx; %to keep consistent with the CRONOS definition, smag is defined with respect to the toroidal flux coordinate
Tex=avg(data.prof.te(itps(1):itps(2),:)).*1e-3;% keV
Nex=avg(data.prof.ne(itps(1):itps(2),:)).*1e-19; %1e19m^-3
Zeffx=avg(data.prof.zeff(itps(1):itps(2),:));

hdt = input('Define only 2 bulk species, hydrogenic + impurities [1]. Use all 5 CRONOS species [0]: ');
impur = zeros(length(x),5);
for i=1:5
  impur(:,i) = squeeze(avg(data.impur.impur(itps(1):itps(2),:,i)));
end

clear Zi Ai ninorm Nix Tix

Tix(:,1)=avg(data.prof.ti(itps(1):itps(2),:)).*1e-3;% keV

vtor=Ro.*avg(data.prof.vtor_exp(itps(1):itps(2),:));% m/s

if hdt == 0 %use all CRONOS species
   nions = 5; 
   Zi = param.compo.z;
   Ai = param.compo.a;
   mi=Ai*mp;
   Nix=zeros(length(x),5);
   ninorm=zeros(length(x),5);
   for i=1:nions
      Nix(:,i) = impur(:,i)./1e19;
      Tix(:,i)=Tix(:,1);
   end
   ninorm =  bsxfun(@times,Nix,1./Nex');
   zeffmtest1 = sum(avg(bsxfun(@times,Zi.^2,ninorm)));
end

if hdt == 1 %only 2 "bulk" species. line weighted averages of hydrogenic species for main ions, and all impurities (weighted by Z^2 to maintain zeff)
   nions = 2;
   for i=1:nions
      Tix(:,i)=Tix(:,1);
   end
   Nix=zeros(length(x),2);
   ninorm=zeros(length(x),2);
   indhdt = find(param.compo.z == 1 & any(impur)>0);
   indimp = find(param.compo.z >  1 & any(impur)>0);
   if length(indhdt) > 1
     
      Ai(1) = sum(param.compo.a(indhdt).*trapz(x,impur(:,indhdt)))./sum(trapz(x,impur(:,indhdt))); % atomic mass of the lump species (weighted by line integral density) 
      Zi(1) = 1;	      
      Nix(:,1) = sum(impur(:,indhdt)')/1e19;
      ninorm(:,1) = Nix(:,1)./Nex';
   else
      Ai(1) = param.compo.a(indhdt);
      Zi(1) = 1;
      Nix(:,1) = impur(:,indhdt);
      ninorm(:,1) = Nix(:,1)./Nex';
   end   
   if length(indimp) > 1
      Ai(2) = sum(param.compo.a(indimp).*trapz(x,impur(:,indimp)))./sum(trapz(x,impur(:,indimp))); % atomic mass of the lump species (weighted by line integral density) 
      Zi(2) = sqrt(sum((param.compo.z(indimp).^2).*trapz(x,impur(:,indimp)))./sum(trapz(x,impur(:,indimp)))); % atomic mass of the lump species (weighted by line integral density) 
      Nix(:,2) = sum(impur(:,indimp)')/1e19;
      ninorm(:,2) = Nix(:,2)./Nex';
   else
      %nmain = (squeeze(impur.impur(:,:,1))-prof.xdur)/1e19;
      Ai(2) = param.compo.a(indimp);
      Zi(2) = param.compo.z(indimp);
      Nix(:,2) = impur(:,indimp)';
      ninorm(:,2) = Nix(:,2)./Nex';
   end  
   zeffmtest2 = sum(avg(bsxfun(@times,Zi.^2,ninorm)));
   zeffprof = sum(bsxfun(@times,Zi.^2,ninorm)')';
   ninorm(:,2) = (zeffprof-1)./(Zi(2)^2 - Zi(2));	     
   ninorm(:,1) = 1 - ninorm(:,2)*Zi(2);
   mi=Ai*mp;
end   




Epsilonx=aavg./Ro;
ft=2.*(2.*Epsilonx).^(0.5)./pi;% fraction of trapped electrions
Pix=Nix.*Tix*1.6*1000;
Pex=Nex.*Tex*1.6*1000;

Betax=2*muo*(sum(Pix')+Pex)./(Bo.*Bo);
rhostar=1./Rmin.*sqrt(Tex*1e3*qe/(2*mp))./(qe*Bo/(2*mp));

%calculate gradients
Ate_aout=-gradient(Tex,aout).*Ro./Tex;
Ate_amid=-gradient(Tex,amid).*Ro./Tex;
Ate_aavg=-gradient(Tex,aavg).*Ro./Tex;
Ane_aout=-gradient(Nex,aout).*Ro./Nex;
Ane_amid=-gradient(Nex,amid).*Ro./Nex;
Ane_aavg=-gradient(Nex,aavg).*Ro./Nex;

vtorgrad_aout=-gradient(vtor,aout);
vtorgrad_amid=-gradient(vtor,amid);
vtorgrad_aavg=-gradient(vtor,aavg);

Ati_aout=zeros(length(x),nions); Ati_amid=zeros(length(x),nions); Ati_aavg=zeros(length(x),nions);
Ani_aout=zeros(length(x),nions); Ani_amid=zeros(length(x),nions); Ani_aavg=zeros(length(x),nions);
for i=1:nions
   Ati_aout(:,i)=-gradient(Tix(:,i)',aout).*Ro./Tix(:,i)';
   Ati_amid(:,i)=-gradient(Tix(:,i)',amid).*Ro./Tix(:,i)';
   Ati_aavg(:,i)=-gradient(Tix(:,i)',aavg).*Ro./Tix(:,i)';
   Ani_aout(:,i)=-gradient(Nix(:,i)',aout).*Ro./Nix(:,i)';
   Ani_amid(:,i)=-gradient(Nix(:,i)',amid).*Ro./Nix(:,i)';
   Ani_aavg(:,i)=-gradient(Nix(:,i)',aavg).*Ro./Nix(:,i)';
end
alpionout=Pix.*(Ani_aout+Ati_aout);
alpionmid=Pix.*(Ani_amid+Ati_amid);
alpionavg=Pix.*(Ani_aavg+Ati_aavg);
alphax_aout=2*muo.*(qx).^2./Bo.^2.*(sum(alpionout')+Pex.*(Ate_aout+Ane_aout));
alphax_amid=2*muo.*(qx).^2./Bo.^2.*(sum(alpionmid')+Pex.*(Ate_amid+Ane_amid));
alphax_aavg=2*muo.*(qx).^2./Bo.^2.*(sum(alpionavg')+Pex.*(Ate_aavg+Ane_aavg));

Zeffmtest = avg(bsxfun(@times,ninorm,Zi.^2));
 
%provides user flexibility in the choice of gradient definition. Default is the midplane average definition
choice=input('Press 1/2/3 for using aout/amid/aavg for gradient definitions (hit enter for the amid default): ' );
if isempty(choice)
	choice=2;
end
switch choice
	case 1
		a=aout;
		Ate=Ate_aout;
		Ati=Ati_aout;
		Atz=Atz_aout;
		Ane=Ane_aout;
		Ani=Ani_aout;
		alphax=alphax_aout;
		vtorgrad=vtorgrad_aout;
		graddef='aout';
	case 2
		a=amid;
		Ate=Ate_amid;
		Ati=Ati_amid;
		Ane=Ane_amid;
		Ani=Ani_amid;
		alphax=alphax_amid;
		vtorgrad=vtorgrad_amid;
		graddef='amid';
	case 3
		a=aavg;
		Ate=Ate_aavg;
		Ati=Ati_aavg;
		Ane=Ane_aavg;
		Ani=Ani_aavg;
		alphax=alphax_aavg;
		vtorgrad=vtorgrad_aavg;
		graddef='aavg';
end
		
%save output		
crondat.Ai=Ai;
crondat.Zi=Zi;
crondat.x=x;
crondat.Bo=Bo;
crondat.Ro=Ro;
crondat.a=a;
crondat.Rmin=Rmin;
crondat.qx=qx;
crondat.smag=smag;
crondat.Tex=Tex;
crondat.Tix=Tix;
crondat.Nex=Nex;
crondat.vtor=vtor;
crondat.Zeffmtest=Zeffmtest;
crondat.Nix=Nix;
crondat.ninorm=ninorm;
crondat.Ate=Ate;
crondat.Ati=Ati;
crondat.Ane=Ane;
crondat.Ani=Ani;
crondat.vtorgrad=vtorgrad;
crondat.graddef=graddef;
crondat.Epsilonx=Epsilonx;
crondat.ft=ft;
crondat.mi=mi;
crondat.Pex=Pex;
crondat.Pix=Pix;
crondat.Betax=Betax;
crondat.alphax=alphax;
crondat.rhostar=rhostar;		
crondat.nions = nions;

%save structure in work directory
eval(['save crondat',num2str(nshot),'.mat crondat']);
	
