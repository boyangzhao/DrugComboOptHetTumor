#!/bin/sh

for subpop1pRange in 50 45 40
	do
	for changeToSubpopOnly in 0
		do
		qsub -N SA_p${subpop1pRange}N${changeToSubpopOnly}_1 ./matjob_batch.sh $subpop1pRange $changeToSubpopOnly
	done
done
