#!/bin/bash

# Check of the universal plot by poincare. We want Ervop3=0.2, but
# generally the arguments of orbitpoincare prescribe W not vperp. So
# when Wpar varies, vperp varies. That variation is minimized if
# W>>Wperp, which needs W>>psi and we can compensate by reducing
# Er. So if we take W to be 2 (vperp=sqrt(2W)=2) and psi to be 0.1 the
# variation in vperp will be negligible. Then to make Ervop3=0.2
# implies Eropsi= Ervop3*sqpsi/v =0.2*(0.1)^(1/2)*/2= 0.0316

E=0.00158
#E=0.0158
#E=0.0316
#E=0.0474
#E=0.0632
#E=0.0790
#E=0.158

for B in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 \
1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 ; do
#for B in 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0; do
./orbitpoincare -p0.1 -W2 -E$E -r40 -d -c  -nw99 -a2  -b$B
done
