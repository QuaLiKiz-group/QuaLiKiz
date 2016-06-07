#!/bin/bash
for file in *
do
	if [[ $file == *".f90" ]] || [[ $file == *".f" ]]
	then
		echo making $file
		make distclean
		if [[ $file == *".f90" ]]
		then
			make ${file/.f90/.mod}
		else
			make ${file/.f/.mod}
		fi
	fi
done
make distclean
