#!/bin/bash
#SBATCH -N1 -c4 -t 2:00:00

multijoin.sh <(zcat mafs/$1.mafs.gz | awk 'NR > 1{printf("%s:%09d\t%f\n",$1,$2,$6)}' OFS="\t")\
 <(zcat mafs/$2.mafs.gz | awk 'NR > 1{printf("%s:%09d\t%f\n",$1,$2,$6)}' OFS="\t")\
 <(zcat mafs/$3.mafs.gz | awk 'NR > 1{printf("%s:%09d\t%f\n",$1,$2,$6)}' OFS="\t")\
 <(zcat mafs/$4.mafs.gz | awk 'NR > 1{printf("%s:%09d\t%f\n",$1,$2,$6)}' OFS="\t")\
 | tr ':' '\t' | tr ' ' '\t' | awk '{print $1,$2+0,$3,$4,$5,$6}' OFS="\t"\
 | fd.py | awk 'NR > 1{print $1,$2-1,$2,$3,$4,$5}' OFS="\t" | bedtools map -a /tigress/noahr/ref/AaegL5_100kb_intervals.bed -b - -c 4,5,6 -o sum > fst/$1.$2.$3.$4.fd
