#!/bin/bash
#SBATCH -N1 -c16 -t 12:00:00
out=thetas/$(basename $1 .bamlist)
angsd -minMapQ 10 -anc /tigress/noahr/ref/male_altref_iter3.fna -gl 1 -doSaf 1 -sites norepeat_5_30x_coverage.sites -bam $1 -out $out -nThreads 16
realSFS ${out}.saf.idx -P 16 > ${out}.sfs
angsd -minMapQ 10 -anc /tigress/noahr/ref/male_altref_iter3.fna -gl 1 -pest ${out}.sfs -doThetas 1 -doSaf 1 -sites norepeat_5_30x_coverage.sites -bam $1 -out $out -nThreads 16
thetaStat do_stat ${out}.thetas.idx -type 2 -win 100000 -step 100000 -outnames ${out}


