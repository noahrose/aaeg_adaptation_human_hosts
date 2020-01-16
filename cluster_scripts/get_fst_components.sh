#!/bin/bash
#SBATCH -N1 -c1 -t 2:00:00

realSFS fst print $1.$2.NC_035107.1.fst.idx | awk '{print $1,$2-1,$2,$3,$4}' OFS="\t" | bedtools map -a /tigress/noahr/ref/AaegL5_100kb_intervals.bed -b - -c 4,5,5 -o sum,sum,count | grep 'NC_035107.1' > ../fst/$1.$2.fst.components.txt
realSFS fst print $1.$2.NC_035108.1.fst.idx | awk '{print $1,$2-1,$2,$3,$4}' OFS="\t" | bedtools map -a /tigress/noahr/ref/AaegL5_100kb_intervals.bed -b - -c 4,5,5 -o sum,sum,count | grep 'NC_035108.1' >> ../fst/$1.$2.fst.components.txt
realSFS fst print $1.$2.NC_035109.1.fst.idx | awk '{print $1,$2-1,$2,$3,$4}' OFS="\t" | bedtools map -a /tigress/noahr/ref/AaegL5_100kb_intervals.bed -b - -c 4,5,5 -o sum,sum,count | grep 'NC_035109.1' >> ../fst/$1.$2.fst.components.txt

