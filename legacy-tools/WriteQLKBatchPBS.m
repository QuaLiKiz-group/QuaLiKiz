% Program which automatically creates the QLK batch scripts for the different machines

function WriteQLKBatchPBS(runname,nprocs)

if nprocs < 12
   nnode=nprocs;
 else
   nnode=12;
end

dire = pwd;
k=1;
clear l

l{k} = ['#PBS -e QuaLiKiz.err\n']; k=k+1;
l{k} = ['#PBS -o QuaLiKiz.out\n']; k=k+1;
l{k} = ['#PBS -N ',runname,'\n']; k=k+1;
l{k} = '#PBS -V\n'; k=k+1;
l{k} = '#PBS -q batch_andromede\n'; k=k+1;
l{k} = ['#PBS -l select=',num2str(ceil(nprocs/12)),':ncpus=',num2str(nnode),'\n']; k=k+1;
l{k} = '\n'; k=k+1;
l{k} = 'cd $PBS_O_WORKDIR'; k=k+1;
l{k} = '\n'; k=k+1;
l{k} = ['/usr/lib64/openmpi/bin/mpiexec -mca btl ^openib -np ',num2str(nprocs),' -machinefile $PBS_NODEFILE QuaLiKiz.exe']; k=k+1;
l{k} = '\n';

fid = fopen('batch_a.pbs','w+');
for i=1:k
	fprintf(fid,l{i});
end

fclose(fid);

if nprocs < 8
   nnode=nprocs;
 else
   nnode=8;
end

dire = pwd;
k=1;
clear l

l{k} = ['#PBS -e QuaLiKiz.err\n']; k=k+1;
l{k} = ['#PBS -o QuaLiKiz.out\n']; k=k+1;
l{k} = ['#PBS -N ',runname,'\n']; k=k+1;
l{k} = '#PBS -V\n'; k=k+1;
l{k} = '#PBS -q batch_cephee\n'; k=k+1;
l{k} = ['#PBS -l select=',num2str(ceil(nprocs/8)),':ncpus=',num2str(nnode),'\n']; k=k+1;
l{k} = '\n'; k=k+1;
l{k} = 'cd $PBS_O_WORKDIR'; k=k+1;
l{k} = '\n'; k=k+1;
l{k} = ['/usr/lib64/openmpi/bin/mpiexec -mca btl ^openib -np ',num2str(nprocs),' -machinefile $PBS_NODEFILE QuaLiKiz.exe']; k=k+1;
l{k} = '\n';


fid = fopen('batch_c.pbs','w+');
for i=1:k
	fprintf(fid,l{i});
end

fclose(fid);


if nprocs < 8
   nnode=nprocs;
 else
   nnode=8;
end

dire = pwd;
k=1;
clear l

l{k} = ['# @ executable = ',dire,'/jacrunscript.sh\n']; k=k+1;
l{k} = ['# @ input = /dev/null\n']; k=k+1;
l{k} = ['# @ output = ',dire,'/QuaLiKiz.out\n']; k=k+1;
l{k} = ['# @ error = ',dire,'/QuaLiKiz.err\n']; k=k+1;
l{k} = ['# @ initialdir = ',dire,'\n']; k=k+1;
l{k} = ['# @ jobtype = openmpi\n']; k=k+1;
l{k} = ['# @ max_processors = ',num2str(nprocs),'\n']; k=k+1;
l{k} = ['# @ min_processors = ',num2str(nprocs),'\n']; k=k+1;
l{k} = ['# @ queue\n']; 

fid = fopen('batch_j','w+');
for i=1:k
	fprintf(fid,l{i});
end

fclose(fid);

dire = pwd;
k=1;
clear l

l{k} = ['mpirun -np ',num2str(nprocs),' ',dire,'/QuaLiKiz.exe']; 

fid = fopen('jacrunscript.sh','w+');
for i=1:k
	fprintf(fid,l{i});
end
unix('chmod 744 jacrunscript.sh');
fclose(fid);

