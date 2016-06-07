% Program which automatically creates the QLK SGE batch script (needed for Zone
% Partenaire)

function WriteQLKBatchSGE(runname)

dire = pwd;

l{1} = '#!/bin/bash\n';  %could also use !/bin/csh\n
l{2} = ['#$ -e ',dire,'/QuaLiKiz.err\n'];
l{3} = ['#$ -o ',dire,'/QuaLiKiz.out\n'];
l{4} = ['#$ -N ',runname,'\n'];
l{5} = '#$ -V\n';
l{6} = '#$ -q batch.q\n';
l{7} = '\n';
l{8} = '\n';
l{9} = 'cd $SGE_O_WORKDIR\n';
l{10} = [dire,'/QuaLiKiz.exe'];
l{11} = '\n';


fid = fopen('batch.sge','w+');
for k=1:11
	fprintf(fid,l{k});
end

fclose(fid);
