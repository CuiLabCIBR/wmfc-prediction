#!/bin/bash

for i_sub in `cat subid_final.txt`

do
    echo "perform of subject: $i_sub"
	sbatch -J ${i_sub} \
		-o log/out.${i_sub}.txt \
		-e log/error.${i_sub}.txt \
		sbatch_get_individual.sh $i_sub
done