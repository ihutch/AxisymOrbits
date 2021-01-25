# This first must be set up with the correct eye.dat file
cp eye1.dat eye.dat
./orbcartpoin helmphiguard.16.dat -r5 -O.9 -nw1 -i -c
mv plot0007.ps gyroholes/3dorbplot1.ps
epstopdf gyroholes/3dorbplot1.ps
mv plot0006.ps  gyroholes/2dorbplot1.ps
epstopdf  gyroholes/2dorbplot1.ps
mv plot0003.ps gyroholes/wpplot1.ps
epstopdf gyroholes/wpplot1.ps
mv plot0001.ps gyroholes/phiofrz1.ps
epstopdf gyroholes/phiofrz1.ps

./orbcartpoin helmphiguard.16.dat -r5 -O.9 -nw1 -f.13 -v3 -c
mv plot0002.ps gyroholes/wpplot2.ps
epstopdf gyroholes/wpplot2.ps

# Poincare plot purple is wp0=-.0120 22 bounces.
./orbcartpoin helmphiguard.16.dat -r5 -O.8939 -nw39 -c
mv plot0001.ps gyroholes/poinplot2.ps
epstopdf  gyroholes/poinplot2.ps

./orbcartpoin helmphiguard.16.dat -r5 -O.5 -nw39 -c
mv plot0001.ps gyroholes/poinplot3.ps
epstopdf  gyroholes/poinplot3.ps

./orbcartpoin helmphiguard.16.dat -r8 -O.5 -c
mv plot0001.ps gyroholes/poinplot4.ps
epstopdf  gyroholes/poinplot4.ps

./orbcartpoin helmphiguard.16.dat -r2 -O.5 -c
mv plot0001.ps gyroholes/poinplot5.ps
epstopdf  gyroholes/poinplot5.ps

# Assuming one has already run helmguardset.sh to make the wpv files.
cp wpv.36.ps gyroholes/
epstopdf gyroholes/wpv.36.ps

cp wpv.04.ps gyroholes/
epstopdf gyroholes/wpv.04.ps

cp wpvwide.36.ps gyroholes/
epstopdf gyroholes/wpvwide.36.ps

cp wpvwide.04.ps gyroholes/
epstopdf gyroholes/wpvwide.04.ps

