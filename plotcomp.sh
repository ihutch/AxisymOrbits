#!/bin/bash
# Make various plots for the paper
#
./orbitpoincare -b.9 -w2.5 -E.005 -c -nw39 -r10
ps2png '-r400 plot0001.ps' poin905.png
./swcontour -b.9 -w2.5 -E.005 -c
epstopdf --outfile=phas905.pdf plot0001.ps 

./orbitpoincare -b.9 -w2.5 -E.04 -c -nw39 -r10
ps2png '-r400 plot0001.ps' poin940.png
./swcontour -b.9 -w2.5 -E.04 -c
epstopdf --outfile=phas940.pdf plot0001.ps 

./orbitpoincare -b.9 -w2.5 -E.1 -c -nw39 -r10
ps2png '-r400 plot0001.ps' poin9100.png
./swcontour -b.9 -w2.5 -E.1 -c
epstopdf --outfile=phas9100.pdf plot0001.ps 


./orbitpoincare -b.6 -w2.5 -E.01 -c -nw39 -r10
ps2png '-r400 plot0001.ps' poin610.png
./swcontour -b.6 -w2.5 -E.01 -c
epstopdf --outfile=phas610.pdf plot0001.ps 

./orbitpoincare -b.6 -w2.5 -E.04 -c -nw39 -r10
ps2png '-r400 plot0001.ps' poin640.png
./swcontour -b.6 -w2.5 -E.04 -c
epstopdf --outfile=phas640.pdf plot0001.ps 

./orbitpoincare -b1.8 -w2.5 -E.2 -c -nw39 -r10
ps2png '-r400 plot0001.ps' poin18200.png
./swcontour -b1.8 -w2.5 -E.2 -c
epstopdf --outfile=phas18200.pdf plot0001.ps 

./orbitpoincare -b1.8 -w2.5 -E.04 -c -nw39 -r10
ps2png '-r400 plot0001.ps' poin1840.png

./swcontour -b1.8 -w2.5 -E.04 -c
epstopdf --outfile=phas1840.pdf plot0001.ps 

./trapdomain -c -t1
ps2png '-r400 plot0001.ps' trapdomain.png

./trapdomain -c -t2 trapdinput011.dat trapdinput11.dat trapdinput31.dat trapdinput101.dat
epstopdf --outfile=tdlineouts.pdf plot0001.ps

./trapdomain -b.5 -t3 -L10 -t3 -c
epstopdf --outfile=vspace05.pdf plot0001.ps

./trapdomain -b1.2 -t3 -L10 -t3 -c
epstopdf --outfile=vspace12.pdf plot0001.ps
