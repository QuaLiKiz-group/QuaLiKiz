phys_meth = 1.0;% Method on additional integrals on particle flux if not 0
coll_flag = 1.0;% Flag for collisionless simulations if 0
rot_flag  = 0.0;% Do not include rotation if 0
verbose   = 1;  % Default verbose output
separateflux = 0;  % Default separateflux output
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

%Set the number and range of wave number points
npoints=8;
kthetarhos = [linspace(0.1,0.8,npoints) linspace(6,48,npoints)]

numn = length(kthetarhos); %number of ktheta in the system
xpoints=3;
scann = xpoints * 3;  %Number of points in parameter scan

%NOTE: for general scans, any of the below can be changed to a vector of size (ones(scann,1))
%e.g. for a q-profile scan: qx=linspace(1,4,scann)'

Bo = 3.*ones(scann,1); %magnetic field
Ro = 3.*ones(scann,1); %major radius
R0 = 3; %geometric major radius
Rmin = 1.*ones(scann,1); %minor radius

x=0.5.*ones(scann,1);	% x is the r/a throughout the scan: impacts the fraction of trapped particles	
rho = x; %assumes circular geometry
qx = 2.00*ones(scann,1); %set q-profile
smag=1.00*ones(scann,1); %set magnetic shear
alphax=0.*ones(scann,1); %MHD alpha

%Electrons
Tex=9.*ones(scann,1); % kev
Nex=5.*ones(scann,1); % ne is in units of 10^19
Ate=[2.*ones(xpoints, 1); linspace(12,2,xpoints)'; 2.*ones(xpoints, 1)]; %logarithmic gradient normalized to R
Ane=[1.*ones(2*xpoints, 1); linspace(4,1,xpoints)']; %logarithmic gradient normalized to R
el_type = 1; % 1 for active, 2 for adiabatic, 3 for adiabatic passing at ion scales (kthetarhosmax<2)
anise = 1.*ones(scann,1); 
danisedr = 0.*ones(scann,1); 
iind=1;
%Ions

ion_name{iind}='Main D';
Ai(:,iind)=2.*ones(scann,1);
Zi(:,iind)=1.*ones(scann,1);
Tix(:,iind)=9.*ones(scann,1); % keV
ninorm(:,iind)=0.8.*ones(scann,1);  ; % ni/ne arbitrary for main ions (will be rewritten to ensure quasineutrality)
Ati(:,iind)=[linspace(12,2,xpoints)'; 2.*ones(xpoints*2, 1)]; %logarithmic gradient normalized to R
Ani(:,iind) = 0.*ones(scann,1); %arbitrary for main ions, will be rewritten to ensure quasineutrality of gradients
ion_type(:,iind)=1.*ones(scann,1); %1 for active, 2 for adiabatic, 3 for tracer (also won't be included in quasineutrality checks), 4 for tracer but included in Zeff
anis(:,iind)=1.*ones(scann,1); %Tperp/Tpar at LFS. Tix defined as Tpar
danisdr(:,iind)=0.*ones(scann,1); %d/dr(Tperp/Tpar) at LFS. Tix defined as Tpar
iind=iind+1;

ion_name{iind}='Be';
Ai(:,iind)=9.*ones(scann,1);
Zi(:,iind)=4.*ones(scann,1);
Tix(:,iind)=9.*ones(scann,1); % keV
ninorm(:,iind)=0.1.*ones(scann,1); % ni/ne
Ati(:,iind)=[linspace(12,2,xpoints)'; 2.*ones(2*xpoints, 1)]; %logarithmic gradient normalized to R
Ani(:,iind) = 0.0.*ones(scann,1);
ion_type(:,iind)=1.*ones(scann,1); %1 for active, 2 for adiabatic, 3 for tracer (also won't be included in quasineutrality checks), 4 for tracer but included in Zeff
anis(:,iind)=1.*ones(scann,1); %Tperp/Tpar at LFS. Tix defined as Tpar
danisdr(:,iind)=0.*ones(scann,1); %d/dr(Tperp/Tpar) at LFS. Tix defined as Tpar
iind=iind+1;

ion_name{iind}='W+42';
Ai(:,iind)=184.*ones(scann,1);
Zi(:,iind)=42.*ones(scann,1);
Tix(:,iind)=9.*ones(scann,1); % keV
ninorm(:,iind)=0.0.*ones(scann,1); % ni/ne
Ati(:,iind)=[linspace(12,2,xpoints)'; 2.*ones(2*xpoints, 1)]; %logarithmic gradient normalized to R
Ani(:,iind) = 0.*ones(scann,1);
ion_type(:,iind)=3.*ones(scann,1); %1 for active, 2 for adiabatic, 3 for tracer, (also won't be included in quasineutrality checks), 4 for tracer but included in Zeff
anis(:,iind)=1.*ones(scann,1); %Tperp/Tpar at LFS. Tix defined as Tpar
danisdr(:,iind)=0.*ones(scann,1); %d/dr(Tperp/Tpar) at LFS. Tix defined as Tpar
iind=iind+1;

nions=iind-1; %number of ions in system

%rotation inputs. Intrinsic assumption that all species rotate together
qe = 1.602176565e-19;
mp = 1.672621777e-27;
cref = sqrt(qe*1e3/mp);
cthi = sqrt(2*qe*Tix(:,1)*1e3/Ai(1)/mp);
epsilon = x./Ro;

Machtor = 0.0.*ones(scann,1); % Toroidal mach number: vtor / cref . Can use cthi value above to scale to desired value with respect to main ion
Autor  = 0.0.*ones(scann,1); % Toroidal velocity shear: dvtor/dr * R/cref 

Machpar = Machtor ./ sqrt(1+(epsilon./qx).^2) ;  % Parallel mach number: vpar / cref . This form assumes pure toroidal rotation in the problem, but free to choose any number
Aupar = Autor./ sqrt(1+(epsilon./qx).^2) ; %Parallel velocity shear: vpar/dr * R/cref . This form assumes pure toroidal rotation in the problem, but free to choose any number

gammaE = -epsilon./qx.*Autor; % - Machtor./qx.*(1-smag./qx)  ; %ExB shear velocity normalized by cref/R. This form assumes pure toroidal rotation in the problem, but free to choose any number


