#!/bin/sh
#$ -cwd
#$ -l mf=50G,h_vmem=70G,h_fsize=2G
#$ -m e -M dengdetian0603@gmail.com

## -t 1-5

Rscript ./AQE_tester_fixSStpr.R 4 1 16 0
