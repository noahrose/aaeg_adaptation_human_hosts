#!/usr/bin/env python
import sys

for line in sys.stdin:  
    l = line.strip().split('\t')  
    p1=float(l[2]) 
    p2=float(l[3]) 
    p3=float(l[4]) 
    p4=float(l[5])
    if p1==0 and p2==0 and p3==0 and p4==0:
        continue

    if p4 == 1 or (p4 != 0 and ((p1+p2+p3+p4)/4)>0.5):
        p1=1-p1
        p2=1-p2
        p3=1-p3
        p4=1-p4
    
    abba=(1-p1)*p2*p3*(1-p4)
    baba=p1*(1-p2)*p3*(1-p4)
    dn=abba-baba
    dd=abba+baba
    
    if p2<=p3:
        fdd=(1-p1)*p3*p3*(1-p4) - p1*(1-p3)*p3*(1-p4)
    else:
        fdd=(1-p1)*p2*p2*(1-p4) - p1*(1-p2)*p2*(1-p4)
    
    print l[0]+"\t"+l[1]+"\t"+str(dn)+"\t"+str(dd)+"\t"+str(fdd)
