#!/bin/bash

light=(1000 2000 4000 8000 16000);
medium=(600 1500 3000 5000 10000);
heavy=(400 500 1000 3000 6000);
processes=(1 2 4);

for i in {0..4}; do
	for p in "${processes[@]}"; do
		echo Processes:$p Particles:$((${light[i]}+${medium[i]}+${heavy[i]}))
		mpirun -np $p x.project ${light[i]} ${medium[i]} ${heavy[i]} 10 50 1 1920 1080 "out/temp"
		# echo $p x.project ${light[i]} ${medium[i]} ${heavy[i]} 5 1 1 1920 1080 "out/temp"
    	echo "";
	done
    echo "";
done