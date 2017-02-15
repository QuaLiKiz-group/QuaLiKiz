%THIS ROUTINE SIMPLY CREATES A NEW RUN DIRECTORY, COPIES OVER THE ESSENTIAL M-FILES NEEDED TO SET UP AND VIEW A RUN, AND SETS THE PATH
%TO THE M-FILE TOOLBOX DIRECTORY SUCH THAT VARIOUS FUNCTIONS NECESSARY (E.G. FOR CREATING THE CRONOS PARAMETER FILE) CAN BE ACCESSED

%TYPICALLY AFTER RUNNING NEWRUN.M, THEN AFTER MODIFYING THE PARAMETERS/PARAMETERS_CRONOS INPUT FILES, THEN INPUT_QLK_SCAN/INPUT_QLK_CRONOS CAN BE RUN
%TO PREPARE THE RUN FOR LAUNCHING

runname=input('Please input desired run directory name, e.g. ''run8_JETITB'': ');

%create runs directory if it doesn't already exist 
if isdir('runs') == 0
	mkdir runs
end
runpath=['runs/',runname];

eval(['!mkdir ',runpath]);
eval(['!cp tools/WriteQLKBatchPBS.m ',runpath,'/.']);
eval(['!cp tools/input_QLK_scan.m ',runpath,'/.']);
eval(['!cp tools/make_CRONOS_parameterfile.m ',runpath,'/.']);
eval(['!cp tools/getcrondat.m ',runpath,'/.']);
eval(['!cp tools/avg.m ',runpath,'/.']);
eval(['!cp tools/parameters_template.m ',runpath,'/parameters.m']);

eval(['cd ',runpath]);

