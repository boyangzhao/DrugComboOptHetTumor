#~/bin/sh
# request Bourne shell as shell for job
#$ -S /bin/sh
#$ -V
#$ -cwd
#$ -m e
#$ -M bozhao@mit.edu

cd ./
matlab -nojvm -nodisplay < "simulate2.m"
