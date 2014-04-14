#!/bin/sh
# request bourne shell as shell for job
#$ -S /bin/sh
#$ -V
#$ -cwd
#$ -m e
#$ -M bozhao@mit.edu

cd ./
matlab -nojvm -nodisplay -r "simulate_batch $1 $2"
#simulate_batch parameters: subpop1pRange changeToSubpopOnly
