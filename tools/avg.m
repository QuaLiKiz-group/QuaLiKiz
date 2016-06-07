%generalization of the mean function to deal correctly with single timeslices
function [res]=avg(f)

sizearr=size(f);
if sizearr(1) > 1
	res=mean(f);
else
	res=f;
end
