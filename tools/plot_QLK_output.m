% generic plotting tool for outputs
% C Bourdelle 27/7/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reading inputs

clear all

kthetarhos=load('kthetarhos.dat');
Ane=load('Ane.dat');
Ani=load('Ani.dat');
Ati=load('Ati.dat');
Ate=load('Ate.dat');
qx=load('qx.dat');
smag=load('smag.dat');
Tex=load('Tex.dat');
Tix=load('Tix.dat');
gammaE=load('gammaE.dat');
Aupar=load('Aupar.dat');
Autor=load('Autor.dat');
Machpar=load('Machpar.dat');
Machtor=load('Machtor.dat');
Bo=load('Bo.dat');
Ro=load('Ro.dat');
R0=load('R0.dat');
Nex=load('Nex.dat');
ninorm=load('ninorm.dat');
x=load('x.dat');
Zi=load('Zi.dat');
alphax=load('alphax.dat');
Ai=load('Ai.dat');
Rmin=load('Rmin.dat');
scann=length(x);

Lambe=1-0.078.*log10(Nex.*0.1)+0.15.*log10(Tex);
Nue=1.36e5.*Lambe.*Nex.*0.1./(Tex.^1.5).*Zi(:,1);
q_ele  = 1.6022e-19;
me     = 9.1094e-31;
cthe=sqrt(2*Tex*1e3*q_ele./me);
Athe=cthe./(qx.*Ro);
Epsilonx=Rmin.*x./Ro;
Nuestar = Nue./(Epsilonx.^1.5.*Athe);
ft=2.*(2.*Epsilonx).^(0.5)./pi; %trapped particle fraction
sizions=size(Tix);
Zeffx=ninorm(:,1).*Zi(:,1).^2;
if sizions(2)>1
for ii=1:sizions(2)
    Zeffx=Zeffx+ninorm(:,ii).*Zi(:,ii).^2;
end
end

% a figure summarizing the main input parameters
figure;
subplot(221)
set(gca,'FontSize',18)
plot(x,Ati(:,1),'r',x,Tix(:,1),'r--',x,Ate,'b',x,Tex,'b--','LineWidth',2)
legend('-R\nabla{T_i}/T_i','Ti','-R\nabla{T_e}/T_e','T_e')
xlabel('\rho')
grid on
subplot(222)
set(gca,'FontSize',18)
plot(x,Ane,'g',x,Nex,'g--','LineWidth',2)
legend('-R\nabla{n_e}/n_e','n_e')
xlabel('\rho')
grid on
subplot(223)
set(gca,'FontSize',18)
plot(x,Tex./Tix(:,1),'c',x,Zeffx,'c--',x,ft,'c-.',x,10*Nuestar,'b--','LineWidth',2)
legend('T_e/T_i','Z_{eff}','f_t', '10*\nu^{*}')
xlabel('\rho')
grid on
subplot(224)
set(gca,'FontSize',18)
plot(x,smag,'m',x,qx,'m--',x,alphax,'m-.s','LineWidth',2)
legend('s','q','\alpha')
xlabel('\rho')
grid on

% reading outputs


gam_GB = load('output/gam_GB.dat');
ome_GB = load('output/ome_GB.dat');
gam_SI = load('output/gam_SI.dat');
ome_SI = load('output/ome_SI.dat');

ief_SI = load('output/ief_SI.dat'); 
ief_GB = load('output/ief_GB.dat'); 

eef_SI = load('output/eef_SI.dat'); 
eef_GB = load('output/eef_GB.dat'); 

ipf_SI = load('output/ipf_SI.dat'); 
ipf_GB = load('output/ipf_GB.dat'); 

epf_SI = load('output/epf_SI.dat'); 
epf_GB = load('output/epf_GB.dat'); 

vce_SI = load('output/vce_SI.dat'); 
vte_SI = load('output/vte_SI.dat'); 
vre_SI = load('output/vre_SI.dat'); 

vci_SI = load('output/vci_SI.dat'); 
vti_SI = load('output/vti_SI.dat'); 
vri_SI = load('output/vri_SI.dat'); 

dfi_SI = load('output/dfi_SI.dat'); 
dfe_SI = load('output/dfe_SI.dat'); 

phi = load('output/phi.dat'); 

rmodewidth = load('output/primitive/rmodewidth.dat');
imodewidth = load('output/primitive/imodewidth.dat');
distan = load('output/primitive/distan.dat');


%%%%% Growth rates 3D plots of the 1st growth rate (can have more than one 
%%%%% unstable branch at a time, here the most unstable of the branches only)

figure; 
set(gca,'FontSize',18)
contourf(x,kthetarhos,gam_SI(1:length(x),:)','LineWidth',2);
hold on
colorbar
shading flat
xlabel('\rho');
ylabel('k_\theta\rho');
title('1st root, \gamma in s^{-1}');


Aa=input('do you want to visualize more details on spectra ? all unstable modes, mode frequencies? type 1 for yes, and 2 for no ');
Aa=1;

if isempty(Aa)==1 || Aa==1
jj=1;
   while jj>0
	x1=input(['at which normalized radius btw ', num2str(min(x)), ' and ', num2str(max(x)),' :'])
	[lab,j]=min(abs(x-x1));
	figure
	subplot(211)
	set(gca,'FontSize',18)
	plot(kthetarhos,gam_SI(j,:),kthetarhos,gam_SI(j+length(x),:),kthetarhos,gam_SI(j+2.*length(x),:),'LineWidth',2)
	hold on
	grid on
	legend('1st root','2nd root','3rd root')
	xlabel('k_{\theta}\rho_s')
	ylabel(['\gamma in s^{-1} at \rho=', num2str(x1)])
	subplot(212)
	set(gca,'FontSize',18)
	plot(kthetarhos,ome_SI(j,:),kthetarhos,ome_SI(j+length(x),:),kthetarhos,ome_SI(j+2.*length(x),:),'LineWidth',2)
	hold on
	grid on	
	title('positive: electron drift direction, negative: ion drift direction')
	xlabel('k_{\theta}\rho_s')
	ylabel(['\omega_r in s^{-1}at \rho=', num2str(x1)])
	bb=input('one more radius? to stop and plot the fluxes type 1 ');
	if bb==1
	   jj=-1;
	end
   end
end


%%%%%% Energy fluxes
figure; 
set(gca,'FontSize',18)
plot(x,eef_SI,'b.--',x,ief_SI,'r.--','LineWidth',2);
xlabel('\rho');
title('energy flux W.m^{-2}');
legend('electron', 'main ion')
hold on
grid on

%%%%% Particle fluxes
%%%% checking ambipolarity
% ion flux from D and V
for i=1:sizions(2) 	   
   Nix(:,i)=ninorm(:,i).*Nex; 
   gdni(:,i)=-Ani(:,i).*Nix(:,i)./R0; 
   Fluxi(:,i)=(-dfi_SI(:,i).*gdni(:,i)+(vci_SI(:,i)+vti_SI(:,i)+vri_SI(:,i)).*Nix(:,i)).*1e19;
end

figure; 
set(gca,'FontSize',18)
plot(x,epf_SI,'b.-',x,sum(ipf_SI.*Zi,sizions(2)),'r.-',x,sum(Fluxi.*Zi,sizions(2)),'m.-','LineWidth',2);
xlabel('\rho');
title('particle flux (m^{-2}s^{-1})');
legend('\Gamma_e', 'Sum Z_i \Gamma_i', 'Sum Z_i (-D\nabla n_i+V n_i)')
grid on
hold on

%electron part. flux
figure; 
subplot(1,2,1)
set(gca,'FontSize',18)
plot(x,vce_SI,x,vte_SI,x,vre_SI,x,vce_SI+vte_SI+vre_SI,'m.-','LineWidth',2);
xlabel('\rho');
title('el. convection velocity (m.s^{-1})');
legend('compressibility', 'thermodiffusion', 'rotodiffusion','total')
hold on 
grid on
subplot(1,2,2)
set(gca,'FontSize',18)
plot(x,dfe_SI,'m.-','LineWidth',2);
xlabel('\rho');
title('el. particle diffusion (m^2.s^{-1})');
hold on 
grid on

%ion part. flux

if sizions(2)>1
   for ii=1:sizions(2)
	figure
	subplot(1,2,1)
	set(gca,'FontSize',18)
	plot(x,vci_SI(:,ii),x,vti_SI(:,ii),x,vri_SI(:,ii),x,vci_SI(:,ii)+vti_SI(:,ii)+vri_SI(:,ii),'m.-','LineWidth',2);
	xlabel('\rho');
	title('ion convection velocity (m.s^{-1})' );
	legend('compressibility', 'thermodiffusion', 'rotodiffusion', 'total')
	hold on 
	grid on
	subplot(1,2,2)
	set(gca,'FontSize',18)
	plot(x,dfi_SI(:,ii),'m.-','LineWidth',2);
	xlabel('\rho');
	title(['ion particle diffusion (m^2.s^{-1}). A=',num2str(Ai(1,ii)),' Z=',num2str(Zi(1,ii))]);
	hold on 
	grid on
   end
end
 
