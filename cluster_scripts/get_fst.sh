#!/bin/bash
#SBATCH -N1 -c16 -t 12:00:00

pop1=$1
pop2=$2
while read r; do
	realSFS thetas/$pop1.saf.idx thetas/$pop2.saf.idx -P 16 -r $r > mafs/$pop1.$pop2.$r.ml
	#prepare the fst for easy window analysis etc
	realSFS fst index thetas/$pop1.saf.idx thetas/$pop2.saf.idx -sfs mafs/$pop1.$pop2.$r.ml -r $r -fstout mafs/$pop1.$pop2.$r -P 16
	#get the global estimate
	realSFS fst stats mafs/$pop1.$pop2.$r.fst.idx 
	realSFS fst stats2 mafs/$pop1.$pop2.$r.fst.idx -win 999999 -step 100000 -type 2 > fst/$pop1.$pop2.$r.1mb.fst
        realSFS fst stats2 mafs/$pop1.$pop2.$r.fst.idx -win 99999 -step 100000 -type 2 > fst/$pop1.$pop2.$r.100kb.fst
done < /tigress/noahr/ref/chrs.txt
echo "getting dxy"
getDxy.pl --pop1maf <(zcat mafs/$pop1.mafs.gz ) --pop2maf <(zcat mafs/$pop2.mafs.gz ) --minInd 1 | awk 'NR > 1{print $1,$2-1,$2,$3}' OFS="\t" | bedtools map -a /tigress/noahr/ref/AaegL5_100kb_intervals.bed -b - -c 4 -o sum > fst/$pop1.$pop2.dxy
