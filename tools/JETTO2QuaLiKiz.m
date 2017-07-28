% matlab programme to prepare a QuaLiKiz input file from a JETTO catalogued run,
% C Bourdelle, 14/6/17

clear all
display('You are reading JETTO output to make QuaLiKiz input');

phys_meth = 1.0;% Method on additional integrals on particle flux if not 0
coll_flag = 1.0;% Flag for collisionless simulations if 0
rot_flag  = 0.0;% Do not include rotation if 0
verbose   = 1;  % Default verbose output
numsols   = 3;  % Number of solutions in the output
nprocs = 8; % Number of processors in parallel computation
maxruns   = 1; %number of runs between contour checks
maxpts    = 5e5; %number of integrand evaluations done in 2D integral
relacc1   = 1e-3; %relative accuracy in 1D integrals
relacc2   = 2e-2; %relative accuracy in 2D integrals
timeout = 60; %seconds after which further solution searching is stopped, at a given wavenumber and radius
ntheta=64; %resolution of parallel direction
set_ninorm1=1;      %flag for automating main ion concentration for maintaining QN
set_Ani1=1; %flag for automating main ion gradient for maintaining QN
set_QN_grad=1; %flag for maintaining quasineutrality of gradients.
ETGmult=1.0; %ETG saturation rule multiplier (for testing)
collmult=1.0; %Collisionality multiplier (for testing)

%Set the number and range of wave number points
kthetarhos = linspace(0.1,0.8,8); %(-) Wave spectrum input: Vector (dimn)

numn = length(kthetarhos); %number of ktheta in the system

shotn=input('What shot number? ');

uidje=input('enter here the JETTO run user ID:', 's');
seqje=input('enter here the JETTO run sequence number:');

% Btor

   [btje,tje]= ...
   jetreaddata(['ppf/',int2str(shotn),'/jst/btor?uid=',uidje,'+seq=',num2str(seqje)]);

tt = input(['time ? (in s) min = ', num2str(min(tje)), ' max=',num2str(max(tje)) , ' : ']);
   [aa,ind] = min(abs(tje-tt));
   
   Bo = btje(ind); %magnetic field

   [rmije,tje,xrhoje]= ...
   jetreaddata(['ppf/',int2str(shotn),'/jsp/rho?uid=',uidje,'+seq=',num2str(seqje)]);
   [aa,ind] = min(abs(tje-tt));

ro=input(['minimum normalized radius, btw ', num2str(min(xrhoje)),' and ',num2str(max(xrhoje)),',default = 0.2: ']);
% more inside large errors on gradients
if isempty(ro) == 1
	ro = 0.2;
end
r1=input(['minimum normalized radius, btw ', num2str(min(xrhoje)),' and ',num2str(max(xrhoje)),',default = 0.8: ']);
% more outside, radiations + open field lines, not included here.
if isempty(r1) == 1
	r1 = 0.8;
end
dimxr=input(' number of radial points, default = 25 : ');
if isempty(dimxr) == 1
	dimxr = 25;
end
   
   scann = dimxr; % for these cases the length of x determines the number of
   % QuaLiKiz spectra calculated

   srmi=size(rmije);

   x=linspace(ro,r1,scann); % (-) radial normalised coordinate (midplane average)

   Rmin = spline(xrhoje,rmije(ind,:),x); %(m) Minor radius. Radial profile due to Shafranov shift
   rho = spline(xrhoje,rmije(ind,:),x)./Rmin(scann); % radial coordinate outboard midplane
   
   Bo = Bo.*ones(1,scann); %magnetic field, (T) Likely not very rigorous to use this sqrt(<B²>) for calculating the Larmor radius % quite close to <Bphi> in practice however 

   [rmaje,tje,xrje]= ...
   jetreaddata(['ppf/',int2str(shotn),'/jsp/r?uid=',uidje,'+seq=',num2str(seqje)]);
   [aa,ind] = min(abs(tje-tt));
   R0 = rmaje(ind,1); %(m) Geometric major radius used for normalizations
   Ro = spline(xrje,rmaje(ind,:),x); %(m) Major radius. Radial profile due to Shafranov shift
   xx=linspace(x(1),x(scann),100); % xx to calculate gradients with more radial
   %points in the outboard midplane
   
   [qje,tje,xqje]= ...
   jetreaddata(['ppf/',int2str(shotn),'/jsp/q?uid=',uidje,'+seq=',num2str(seqje)]);
   [aa,ind] = min(abs(tje-tt));
   qx = spline(xqje,qje(ind,:),x); % safety factor
   qxx=spline(x,qx,xx);
   rhoxx=spline(x,rho,xx);
   qprimxx=gradient(qxx,rhoxx.*Rmin(scann));
   qprim=spline(xx,qprimxx,x);
   smag=rho.*Rmin(scann)./qx.*qprim; % magnetic shear calculated in OMP

   
%Electrons
   [teje,tje,xteje]= ...
   jetreaddata(['ppf/',int2str(shotn),'/jsp/te?uid=',uidje,'+seq=',num2str(seqje)]);
   [aa,ind] = min(abs(tje-tt));
   Tex = spline(xteje,teje(ind,:),x).*1e-3; %(keV) Vector (radial grid x(aa))
   Texx=spline(x,Tex,xx);
   Atexx=-gradient(Texx,rhoxx.*Rmin(scann)).*R0./Texx;  
   Ate=spline(xx,Atexx,x);
   [neje,tje,xneje]= ...
   jetreaddata(['ppf/',int2str(shotn),'/jsp/ne?uid=',uidje,'+seq=',num2str(seqje)]);
   [aa,ind] = min(abs(tje-tt));
   Nex = spline(xneje,neje(ind,:),x).*1e-19; % (10^19 m^-3) Vector (radial grid x(aa))
   Nexx=spline(x,Nex,xx);
   Anexx=-gradient(Nexx,rhoxx.*Rmin(scann)).*R0./Nexx;  
   Ane=spline(xx,Anexx,x);
   el_type = 1; % 1 for active, 2 for adiabatic, 3 for adiabatic passing at ion scales (kthetarhosmax<2)
   anise = 1.*ones(1,scann); % Tperp/Tpar at LFS
   danisedr = 0.*ones(1,scann); % d/dr(Tperp/Tpar) at LFS
    
%Ion inputs (can be for as many ions as one wants)

   iind=1;
   ion_name{iind}='Main ions';
   [atmje,tje]= ...
   jetreaddata(['ppf/',int2str(shotn),'/jst/atm?uid=',uidje,'+seq=',num2str(seqje)]);
   [aa,ind] = min(abs(tje-tt));
   Ai(iind,:)=atmje(ind).*ones(1,scann);
   if atmje(ind)==4
       Zi(iind,:)=2.*ones(1,scann); % then the main ions are He
   else
       Zi(iind,:)=1.*ones(1,scann); % else they are D, T or H  
   end 
   [tije,tje,xtije]= ...
   jetreaddata(['ppf/',int2str(shotn),'/jsp/ti?uid=',uidje,'+seq=',num2str(seqje)]);
   [aa,ind] = min(abs(tje-tt));   
   Tix(iind,:) = spline(xtije,tije(ind,:),x).*1e-3; %(keV) Vector (radial grid x(aa))
   Tixx=spline(x,Tix(iind,:),xx);
   Atixx=-gradient(Tixx,rhoxx.*Rmin(scann)).*R0./Tixx;  
   Ati(iind,:)=spline(xx,Atixx,x);
   ninorm(iind,:)=0.8.*ones(1,scann);  ; % ni/ne arbitrary for main ions (will be rewritten to ensure quasineutrality)
   Ani(iind,:) = Ane; %arbitrary for main ions, will be rewritten to ensure quasineutrality of gradients
   ion_type(iind,:)=1.*ones(1,scann); %1 for active, 2 for adiabatic, 3 for tracer (also won't be included in quasineutrality checks), 4 for tracer but included in Zeff
   anis(iind,:)=1.*ones(1,scann); %Tperp/Tpar at LFS. Tix defined as Tpar
   danisdr(iind,:)=0.*ones(1,scann); %d/dr(Tperp/Tpar) at LFS. Tix defined as Tpar
   iind=iind+1;


   ion_name{iind}='Impurity';
   % pb only one impurity density available in output of JETTO, so multiple ions
   % cannot be read yet although QuaLiKiz can deal with infinite number of ions
   [nimpje,tje,xnimpje]= ...
   jetreaddata(['ppf/',int2str(shotn),'/jsp/nimp?uid=',uidje,'+seq=',num2str(seqje)]);
   [aa,ind] = min(abs(tje-tt)); 
   % Since JETTO v170517 3 impurities are now catalogued 
   % to be reactivated if atmi and zipi are indeed available in jss, in principle
   % contain up to 3 impurity ions A and Z.
   % and in jss: nim1, nim2, nim3 for the densities, zia1, zia2, zia3 for 
   % the profile of averaged charge
   %[atmije,tje]= ...
   %jetreaddata(['ppf/',int2str(shotn),'/jst/atmi?uid=',uidje,'+seq=',num2str(seqje)]);
   %[aa,ind] = min(abs(tje-tt));
   %Ai(iind,:)=atmije(1,ind).*ones(1,scann);
   %[zipije,tje]= ...
   %jetreaddata(['ppf/',int2str(shotn),'/jst/zipi?uid=',uidje,'+seq=',num2str(seqje)]);
   %[aa,ind] = min(abs(tje-tt));
   %Zi(iind,:)=zipije(1,ind).*ones(1,scann);
   % for now assume Be is the only impurity
   Ai(iind,:)=9.*ones(1,scann);
   Zi(iind,:)=4.*ones(1,scann); 
   Tix(iind,:)=Tix(1,:); % keV, assume all ions at the same temperature
   Nzx = spline(xnimpje,nimpje(ind,:),x).*1e-19; % (10^19 m^-3) Vector (radial grid x(aa))
   ninorm(iind,:)=Nzx./Nex; % ni/ne
   Ati(iind,:)=Ati(1,:);% assume all ions at the same temperature
   Nzxx=spline(x,Nzx,xx);
   Anzxx=-gradient(Nzxx,rhoxx.*Rmin(scann)).*R0./Nzxx;  
   Anz=spline(xx,Anzxx,x);
   Ani(iind,:) = Anz; 
   ion_type(iind,:)=1.*ones(1,scann); %1 for active, 2 for adiabatic, 3 for tracer (also won't be included in quasineutrality checks), 4 for tracer but included in Zeff
   anis(iind,:)=1.*ones(1,scann); %Tperp/Tpar at LFS. Tix defined as Tpar
   danisdr(iind,:)=0.*ones(1,scann); %d/dr(Tperp/Tpar) at LFS. Tix defined as Tpar
   iind=iind+1;

   nions=iind-1; %number of ions in system

%rotation inputs. Intrinsic assumption that all species rotate together
   qe = 1.602176565e-19;
   mp = 1.672621777e-27;
   cref = sqrt(qe*1e3/mp);
   cthi = sqrt(2*qe*Tix(1,:)*1e3/Ai(1)/mp);
   epsilon = x./Ro;

   [vtorje,tje,xvtorje]= ...
   jetreaddata(['ppf/',int2str(shotn),'/jsp/vtor?uid=',uidje,'+seq=',num2str(seqje)]);
   [aa,ind] = min(abs(tje-tt));
   Vtorx = spline(xvtorje,vtorje(ind,:),x); %(m/s) toroidal velocity at Router
   Machtor = Vtorx./cref; % Toroidal mach number: vtor / cref . Can use cthi value above to scale to desired value with respect to main ion
   Vtorxx=spline(x,Vtorx,xx);
   Autorxx=gradient(Vtorxx,rhoxx.*Rmin(scann)).*R0./cref;  % Toroidal velocity shear: dvtor/dr * R/cref
   Autor=spline(xx,Autorxx,x);
   Machpar = Machtor ./ sqrt(1+(epsilon./qx).^2) ;  % Parallel mach number: vpar / cref . This form assumes pure toroidal rotation in the problem, but free to choose any number
   Aupar = Autor./ sqrt(1+(epsilon./qx).^2) ; %Parallel velocity shear: vpar/dr * R/cref . This form assumes pure toroidal rotation in the problem, but free to choose any number

%   gammaE = -epsilon./qx.*Autor; %ExB shear rate normalized by cref/R. This form assumes pure toroidal rotation in the problem, but free to choose any number
   [omebje,tje,xomebje]= ...
   jetreaddata(['ppf/',int2str(shotn),'/jsp/omeb?uid=',uidje,'+seq=',num2str(seqje)]);
   [aa,ind] = min(abs(tje-tt));
   omebx = spline(xomebje,omebje(ind,:),x); %(s^-1) ExB shearing rate
   gammaE =omebx./(cref./R0); %ExB shear rate normalized by cref/R.
   
   % calculating alpha, self-consistently with B, P and grad(P)
   muo=4*pi*1e-7;
   %need ninorm from main ions here... 
   ninorm(1,:)=ones(1,scann);
   if (nions==1)
      ninorm(1,:)=ones(1,scann);
   else
      for ii=1:nions-1
         ninorm(1,:)=ninorm(1,:)-Zi(ii+1,:).*ninorm(ii+1,:);
      end
   end
   Pitot=zeros(1,scann);
   gradPitot=zeros(1,scann);
   for ii=1:nions
      Pix(ii,:)=ninorm(ii,:).*Nex.*Tix(ii,:).*1.6.*1000;
      Pitot=Pitot+Pix(ii,:);
      gradPitot=gradPitot+Pix(ii,:).*(Ani(ii,:)+Ati(ii,:));
   end
   Pex=Nex.*Tex.*1.6.*1000;
   Betax=2.*muo*(Pitot+Pex)./(Bo.*Bo);
   alphax=2.*muo.*qx.^2./Bo.^2.*(gradPitot+Pex.*(Ate+Ane));
