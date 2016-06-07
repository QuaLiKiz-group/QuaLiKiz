% Script which writes a parameters.m file from CRONOS data, 
% which should be saved in the run directory as crondat.mat
% the presaved mat file needs to be loaded into the workspace

% For now limited to just 1 carbon impurity. Can easily be generalized if needed

nprocs = input('Input number of processors for parallel job (e.g. 4): ');
kthetarhos = input('Input kthetarhos array (e.g. [0.1 0.2 0.4 0.5 0.8 1.5 2.3 4 10 15 20 30]): ');

numn = length(kthetarhos);

radius = input('Input radial range in scan (e.g. [0.2 0.7]): ');
scann = input('Input number of radial points in scan (e.g. 6): ');
xqlk = linspace(radius(1),radius(2),scann);

disp('parameter.m file written. All input variables can also be changed manually');
disp('Run input_QLK_scan.m to create the binary input files and batch script');

k=1;
xcro=linspace(0,1,101);

l{k} = ['phys_meth = 1.0;','\n']; k=k+1;
l{k} = ['coll_flag = 1.0;','\n']; k=k+1;
l{k} = ['rot_flag  = 0.0;','\n']; k=k+1;
l{k} = ['numsols = 3;','\n']; k=k+1;
l{k} = ['nprocs = ',num2str(nprocs),';','\n']; k=k+1;
l{k} = ['maxruns   = 1;','\n']; k=k+1;
l{k} = ['maxpts    = 5e5;','\n']; k=k+1;
l{k} = ['relacc1   = 1e-3;','\n']; k=k+1;
l{k} = ['relacc2   = 2e-2;','\n']; k=k+1;
l{k} = ['timeout   = 60;','\n']; k=k+1;
l{k} = ['ntheta=64;','\n']; k=k+1;
l{k} = ['numecoefs = 13;','\n']; k=k+1;

l{k} = ['set_ninorm1=1;','\n']; k=k+1;
l{k} = ['set_Ani1=1;','\n']; k=k+1;
l{k} = ['set_QN_grad=1;','\n']; k=k+1;

l{k} = '\n'; k=k+1;

l{k} = ['kthetarhos = [',num2str(kthetarhos),'];','\n']; k=k+1;
l{k} = ['numn = ',num2str(numn),';','\n']; k=k+1;
l{k} = ['scann = ',num2str(scann),';','\n']; k=k+1;
l{k} = '\n'; k=k+1;

l{k} = ['Bo = [',num2str(interp1(xcro,crondat.Bo,xqlk)),']'';','\n']; k=k+1;
l{k} = ['Ro = [',num2str(interp1(xcro,crondat.Ro,xqlk)),']'';','\n']; k=k+1;
l{k} = ['Rmin = [',num2str(interp1(xcro,crondat.Rmin,xqlk)),']'';','\n']; k=k+1;
l{k} = '\n'; k=k+1;

l{k} = ['x = [',num2str(interp1(xcro,crondat.x,xqlk)),']'';','\n']; k=k+1;
l{k} = ['qx = [',num2str(interp1(xcro,crondat.qx,xqlk)),']'';','\n']; k=k+1;
l{k} = ['smag = [',num2str(interp1(xcro,crondat.smag,xqlk)),']'';','\n']; k=k+1;
l{k} = ['alphax = [',num2str(interp1(xcro,crondat.alphax,xqlk)),']'';','\n']; k=k+1;
l{k} = '\n'; k=k+1;

l{k} = ['%Electrons','\n']; k=k+1;
l{k} = ['Tex = [',num2str(interp1(xcro,crondat.Tex,xqlk)),']'';','\n']; k=k+1;
l{k} = ['Nex = [',num2str(interp1(xcro,crondat.Nex,xqlk)),']'';','\n']; k=k+1;
l{k} = ['Ate = [',num2str(interp1(xcro,crondat.Ate,xqlk)),']'';','\n']; k=k+1;
l{k} = ['Ane = [',num2str(interp1(xcro,crondat.Ane,xqlk)),']'';','\n']; k=k+1;
l{k} = ['el_type = 1; % 1 for active, 2 for adiabatic, 3 for adiabatic passing at ion scales (kthetarhosmax<2)','\n']; k=k+1;
l{k} = ['anise = 1.*ones(scann,1);','\n']; k=k+1;
l{k} = ['danisedr = 0.*ones(scann,1);','\n']; k=k+1; 
l{k} = ['iind=1;','\n']; k=k+1;
iind=1;

l{k} = '\n'; k=k+1;
l{k} = ['%Ion 1','\n']; k=k+1;
l{k} = ['ion_name{iind}=''Ion 1'';','\n']; k=k+1;
l{k} = ['Ai(iind) = ',num2str(crondat.Ai(iind)),';','\n']; k=k+1;
l{k} = ['Zi(iind) = ',num2str(crondat.Zi(iind)),';','\n']; k=k+1;
l{k} = ['Tix(:,iind) = [',num2str(interp1(xcro,crondat.Tix(:,iind),xqlk)),'];','\n']; k=k+1;
l{k} = ['ninorm(:,iind) = [',num2str(interp1(xcro,crondat.ninorm(:,iind),xqlk)),'];','\n']; k=k+1;
l{k} = ['Ati(:,iind) = [',num2str(interp1(xcro,crondat.Ati(:,iind),xqlk)),'];','\n']; k=k+1;
l{k} = ['Ani(:,iind) = [',num2str(interp1(xcro,crondat.Ani(:,iind),xqlk)),'];','\n']; k=k+1;
if crondat.Zeffmtest(iind) < 0.1
  iontype = 4;
else
  iontype = 1;
end
l{k} = ['ion_type(iind) = ',num2str(iontype),'; % 1 for active, 2 for adiabatic, 3  for tracer (also won''t be included in quasineutrality checks), 4 for tracer but included in Zeff','\n']; k=k+1;
l{k} = ['anis(:,iind) = 1.*ones(scann,1);','\n']; k=k+1;
l{k} = ['danisdr(:,iind) = 0.*ones(scann,1);','\n']; k=k+1; 
l{k} = ['iind=iind+1;','\n']; k=k+1;
iind=iind+1;

l{k} = '\n'; k=k+1;
l{k} = ['%Ion 2','\n']; k=k+1;
l{k} = ['ion_name{iind}=''Ion 2'';','\n']; k=k+1;
l{k} = ['Ai(iind) = ',num2str(crondat.Ai(iind)),';','\n']; k=k+1;
l{k} = ['Zi(iind) = ',num2str(crondat.Zi(iind)),';','\n']; k=k+1;
l{k} = ['Tix(:,iind) = [',num2str(interp1(xcro,crondat.Tix(:,iind),xqlk)),'];','\n']; k=k+1;
l{k} = ['ninorm(:,iind) = [',num2str(interp1(xcro,crondat.ninorm(:,iind),xqlk)),'];','\n']; k=k+1;
l{k} = ['Ati(:,iind) = [',num2str(interp1(xcro,crondat.Ati(:,iind),xqlk)),'];','\n']; k=k+1;
l{k} = ['Ani(:,iind) = [',num2str(interp1(xcro,crondat.Ani(:,iind),xqlk)),'];','\n']; k=k+1;
if crondat.Zeffmtest(iind) < 0.1
  iontype = 4;
else
  iontype = 1;
end
l{k} = ['ion_type(iind) = ',num2str(iontype),'; % 1 for active, 2 for adiabatic, 3  for tracer (also won''t be included in quasineutrality checks), 4 for tracer but included in Zeff','\n']; k=k+1;
l{k} = ['anis(:,iind) = 1.*ones(scann,1);','\n']; k=k+1;
l{k} = ['danisdr(:,iind) = 0.*ones(scann,1);','\n']; k=k+1; 
l{k} = ['iind=iind+1;','\n']; k=k+1;
iind=iind+1;

if crondat.nions > 2
l{k} = '\n'; k=k+1;
l{k} = ['%Ion 3','\n']; k=k+1;
l{k} = ['ion_name{iind}=''Ion 3'';','\n']; k=k+1;
l{k} = ['Ai(iind) = ',num2str(crondat.Ai(iind)),';','\n']; k=k+1;
l{k} = ['Zi(iind) = ',num2str(crondat.Zi(iind)),';','\n']; k=k+1;
l{k} = ['Tix(:,iind) = [',num2str(interp1(xcro,crondat.Tix(:,iind),xqlk)),'];','\n']; k=k+1;
l{k} = ['ninorm(:,iind) = [',num2str(interp1(xcro,crondat.ninorm(:,iind),xqlk)),'];','\n']; k=k+1;
l{k} = ['Ati(:,iind) = [',num2str(interp1(xcro,crondat.Ati(:,iind),xqlk)),'];','\n']; k=k+1;
l{k} = ['Ani(:,iind) = [',num2str(interp1(xcro,crondat.Ani(:,iind),xqlk)),'];','\n']; k=k+1;
if crondat.Zeffmtest(iind) < 0.1
  iontype = 4;
else
  iontype = 1;
end
l{k} = ['ion_type(iind) = ',num2str(iontype),'; % 1 for active, 2 for adiabatic, 3  for tracer (also won''t be included in quasineutrality checks), 4 for tracer but included in Zeff','\n']; k=k+1;
l{k} = ['anis(:,iind) = 1.*ones(scann,1);','\n']; k=k+1;
l{k} = ['danisdr(:,iind) = 0.*ones(scann,1);','\n']; k=k+1; 
l{k} = ['iind=iind+1;','\n']; k=k+1;
iind=iind+1;

l{k} = '\n'; k=k+1;
l{k} = ['%Ion 4','\n']; k=k+1;
l{k} = ['ion_name{iind}=''Ion 4'';','\n']; k=k+1;
l{k} = ['Ai(iind) = ',num2str(crondat.Ai(iind)),';','\n']; k=k+1;
l{k} = ['Zi(iind) = ',num2str(crondat.Zi(iind)),';','\n']; k=k+1;
l{k} = ['Tix(:,iind) = [',num2str(interp1(xcro,crondat.Tix(:,iind),xqlk)),'];','\n']; k=k+1;
l{k} = ['ninorm(:,iind) = [',num2str(interp1(xcro,crondat.ninorm(:,iind),xqlk)),'];','\n']; k=k+1;
l{k} = ['Ati(:,iind) = [',num2str(interp1(xcro,crondat.Ati(:,iind),xqlk)),'];','\n']; k=k+1;
l{k} = ['Ani(:,iind) = [',num2str(interp1(xcro,crondat.Ani(:,iind),xqlk)),'];','\n']; k=k+1;
if crondat.Zeffmtest(iind) < 0.1
  iontype = 4;
else
  iontype = 1;
end
l{k} = ['ion_type(iind) = ',num2str(iontype),'; % 1 for active, 2 for adiabatic, 3  for tracer (also won''t be included in quasineutrality checks), 4 for tracer but included in Zeff','\n']; k=k+1;
l{k} = ['anis(:,iind) = 1.*ones(scann,1);','\n']; k=k+1;
l{k} = ['danisdr(:,iind) = 0.*ones(scann,1);','\n']; k=k+1; 
l{k} = ['iind=iind+1;','\n']; k=k+1;
iind=iind+1;

l{k} = '\n'; k=k+1;
l{k} = ['%Ion 5','\n']; k=k+1;
l{k} = ['ion_name{iind}=''Ion 5'';','\n']; k=k+1;
l{k} = ['Ai(iind) = ',num2str(crondat.Ai(iind)),';','\n']; k=k+1;
l{k} = ['Zi(iind) = ',num2str(crondat.Zi(iind)),';','\n']; k=k+1;
l{k} = ['Tix(:,iind) = [',num2str(interp1(xcro,crondat.Tix(:,iind),xqlk)),'];','\n']; k=k+1;
l{k} = ['ninorm(:,iind) = [',num2str(interp1(xcro,crondat.ninorm(:,iind),xqlk)),'];','\n']; k=k+1;
l{k} = ['Ati(:,iind) = [',num2str(interp1(xcro,crondat.Ati(:,iind),xqlk)),'];','\n']; k=k+1;
l{k} = ['Ani(:,iind) = [',num2str(interp1(xcro,crondat.Ani(:,iind),xqlk)),'];','\n']; k=k+1;
if crondat.Zeffmtest(iind) < 0.1
  iontype = 4;
else
  iontype = 1;
end
l{k} = ['ion_type(iind) = ',num2str(iontype),'; % 1 for active, 2 for adiabatic, 3  for tracer (also won''t be included in quasineutrality checks), 4 for tracer but included in Zeff','\n']; k=k+1;
l{k} = ['anis(:,iind) = 1.*ones(scann,1);','\n']; k=k+1;
l{k} = ['danisdr(:,iind) = 0.*ones(scann,1);','\n']; k=k+1; 
l{k} = ['iind=iind+1;','\n']; k=k+1;
l{k} = '\n'; k=k+1;
iind=iind+1;
end
l{k} = ['nions=iind-1;','\n']; k=k+1; 
l{k} = '\n'; k=k+1;
%rotation inputs. Intrinsic assumption that all species rotate together. For now set to zero. 

l{k} = ['qe = 1.602176565e-19;','\n']; k=k+1; 
l{k} = ['mp = 1.672621777e-27;','\n']; k=k+1; 
l{k} = ['cref = sqrt(qe*1e3/mp);','\n']; k=k+1; 
l{k} = ['cthi = sqrt(2*qe*Tix(:,1)*1e3/Ai(1)/mp);','\n']; k=k+1; 
l{k} = ['epsilon = x./Ro;','\n']; k=k+1; 

l{k} = '\n'; k=k+1;

l{k} = ['Machtor = [',num2str(interp1(xcro,crondat.vtor,xqlk)),']''./cref;','\n']; k=k+1;
l{k} = ['Autor   = [',num2str(interp1(xcro,crondat.vtorgrad,xqlk)),']''.*Ro./cref;','\n']; k=k+1;


l{k} = ['Machpar = Machtor./sqrt(1+(epsilon./qx).^2);','\n']; k=k+1;
l{k} = ['gammaE = -epsilon./qx.*Autor;','\n']; k=k+1;
l{k} = ['Aupar = Autor./ sqrt(1+(epsilon./qx).^2);','\n']; 

fid = fopen('parameters.m','w+');
for i=1:k
	fprintf(fid,l{i});
end

fclose(fid);
