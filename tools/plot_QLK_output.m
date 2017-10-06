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

sizions=size(Tix);
Zeffx=ninorm(:,1).*Zi(:,1).^2;
if sizions(2)>1
for ii=1:sizions(2)
    Zeffx=Zeffx+ninorm(:,ii).*Zi(:,ii).^2;
end
end
Lambe=1-0.078.*log10(Nex.*0.1)+0.15.*log10(Tex);
Nue=1.36e5.*Lambe.*Nex.*0.1./(Tex.^1.5).*Zeffx;
q_ele  = 1.6022e-19;
me     = 9.1094e-31;
cthe=sqrt(2*Tex*1e3*q_ele./me);
Athe=cthe./(qx.*Ro);
Epsilonx=Rmin.*x./Ro;
Nuestar = Nue./(Epsilonx.^1.5.*Athe);
ft=2.*(2.*Epsilonx).^(0.5)./pi; %trapped particle fraction

% a figure summarizing the main input parameters
figure;
subplot(221)
set(gca,'FontSize',18)
plot(x,Ati(:,1),'r',x,Tix(:,1),'r--',x,Ate,'b',x,Tex,'b--','LineWidth',2)
l1=legend('-$R\nabla{T_i}/T_i$','$T_i$','-$R\nabla{T_e}/T_e$','$T_e$')
set(l1,'Interpreter','latex')
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')
grid on
subplot(222)
set(gca,'FontSize',18)
plot(x,Ane,'g',x,Nex,'g--','LineWidth',2)
l1=legend('-$R\nabla{n_e}/n_e$','$n_e$')
set(l1,'Interpreter','latex')
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')
grid on
subplot(223)
set(gca,'FontSize',18)
plot(x,Tex./Tix(:,1),'c',x,Zeffx,'c--',x,ft,'c-.',x,10*Nuestar,'b--','LineWidth',2)
l1=legend('$T_e/T_i$','$Z_{eff}$','$f_t$', '$10\times\nu_e^{*}$')
set(l1,'Interpreter','latex')
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')
grid on
subplot(224)
set(gca,'FontSize',18)
plot(x,smag,'m',x,qx,'m--',x,alphax,'m-.s','LineWidth',2)
l1=legend('s','q','$\alpha$')
set(l1,'Interpreter','latex')
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')
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
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')
l3=ylabel('$k_\theta\rho$');
set(l3,'Interpreter','latex')
l1=title('1st root, $\gamma$ in $s^{-1}$');
set(l1,'Interpreter','latex')


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
	l1=legend('1st root','2nd root','3rd root')
	set(l1,'Interpreter','latex')
	l2=xlabel('$k_{\theta}\rho_s$')
	set(l2,'Interpreter','latex')
	l3=ylabel(['$\gamma$ in $s^{-1}$ at $\rho=$', num2str(x1)])
	set(l3,'Interpreter','latex')
	subplot(212)
	set(gca,'FontSize',18)
	plot(kthetarhos,ome_SI(j,:),kthetarhos,ome_SI(j+length(x),:),kthetarhos,ome_SI(j+2.*length(x),:),'LineWidth',2)
	hold on
	grid on	
	l1=title('positive: electron drift direction, negative: ion drift direction')
	set(l1,'Interpreter','latex')
	l2=xlabel('$k_{\theta}\rho_s$')
	set(l2,'Interpreter','latex')
	l3=ylabel(['$\omega_r$ in $s^{-1}$ at $\rho=$', num2str(x1)])
	set(l3,'Interpreter','latex')
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
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')
l3=title('energy flux $W.m^{-2}$');
set(l3,'Interpreter','latex')
l1=legend('electron', 'main ion')
set(l1,'Interpreter','latex')
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
l2=xlabel('$\rho$');
set(l2,'Interpreter','latex')
l3=title('particle flux $m^{-2}s^{-1}$');
set(l3,'Interpreter','latex')
l1=legend('$ \Gamma_e $', '$ \sum Z_i \Gamma_i $', '$  \sum Z_i (-D\nabla n_i+V n_i)$')
set(l1,'Interpreter','latex')
grid on
hold on

%electron part. flux
figure; 
subplot(1,2,1)
set(gca,'FontSize',18)
plot(x,vce_SI,x,vte_SI,x,vre_SI,x,vce_SI+vte_SI+vre_SI,'m.-','LineWidth',2);
l1=xlabel('$\rho$');
set(l1,'Interpreter','latex')
l2=title('el. convection velocity ($m.s^{-1}$)');
set(l2,'Interpreter','latex')
l1=legend('compressibility', 'thermodiffusion', 'rotodiffusion','total')
set(l1,'Interpreter','latex')
hold on 
grid on
subplot(1,2,2)
set(gca,'FontSize',18)
plot(x,dfe_SI,'m.-','LineWidth',2);
l1=xlabel('$\rho$');
set(l1,'Interpreter','latex')
l2=title('el. particle diffusion ($m^2.s^{-1}$)');
set(l2,'Interpreter','latex')
hold on 
grid on

%ion part. flux

if sizions(2)>1
   for ii=1:sizions(2)
	figure
	subplot(1,2,1)
	set(gca,'FontSize',18)
	plot(x,vci_SI(:,ii),x,vti_SI(:,ii),x,vri_SI(:,ii),x,vci_SI(:,ii)+vti_SI(:,ii)+vri_SI(:,ii),'m.-','LineWidth',2);
	l1=xlabel('$\rho$');
	set(l1,'Interpreter','latex')
	l2=title('ion convection velocity (m.s^{-1})' );
	set(l2,'Interpreter','latex')
	l3=legend('compressibility', 'thermodiffusion', 'rotodiffusion', 'total')
	set(l3,'Interpreter','latex')
	hold on 
	grid on
	subplot(1,2,2)
	set(gca,'FontSize',18)
	plot(x,dfi_SI(:,ii),'m.-','LineWidth',2);
	l1=xlabel('$\rho$');
	set(l1,'Interpreter','latex')
	l2=title(['ion particle diffusion ($m^2.s^{-1}$). A=',num2str(Ai(1,ii)),' Z=',num2str(Zi(1,ii))]);
	set(l2,'Interpreter','latex')
	hold on 
	grid on
   end
end
 
