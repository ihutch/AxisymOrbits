#!/bin/bash

for PVAL in .36 .25 .16 .09 .04 .01; do
    helmhole/helmhole -p$PVAL -ps-3
    mv helmphiguard.dat helmphiguard${PVAL}.dat
done

for PVAL in .36 .25 .16 .09 .04 .01; do
    OM=`echo "sqrt($PVAL)*3.00" | bc`
    echo $OM
    ./orbcartpoin helmphiguard${PVAL}.dat -O${OM} -mw1.1 -mE5 -c >/dev/null
    mv plot0001.ps wpv${PVAL}.ps
done

#More peaked radial profiles
for PVAL in .36 .25 .16 .09 .04 .01; do
    helmhole/helmhole -rm18 -rt10 -sh.03 -p$PVAL -ps-3
    mv helmphiguard.dat helmphipeak${PVAL}.dat
done

for PVAL in .36 .25 .16 .09 .04 .01; do
    OM=`echo "sqrt($PVAL)*3.00" | bc`
    echo $OM
    ./orbcartpoin helmphipeak${PVAL}.dat -O${OM} -mw1.1 -mE5 -c >/dev/null
    mv plot0001.ps wpvpeak${PVAL}.ps
done

# Even wider profiles with smaller alpha

for PVAL in .36 .25 .16 .09 .04 .01; do
    helmhole/helmhole -rm25 -rt18 -sh.01 -p$PVAL -a.3 -ps-3
    mv helmphiguard.dat helmphiwide${PVAL}.dat
done


for PVAL in .36 .25 .16 .09 .04 .01; do
    OM=`echo "sqrt($PVAL)*3.00" | bc`
    echo $OM
    ./orbcartpoin helmphiwide${PVAL}.dat -O${OM} -mw1.1 -mE5 -c >/dev/null
    mv plot0001.ps wpvwide${PVAL}.ps
done
