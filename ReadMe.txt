These codes integrate the orbit in 6 phase-space dimensions of an
electron moving in the solitary potential peak of an electron hole.

orbitpoincare.f90 was used for the analysis reported in the article
"Particle Trapping in Axisymmetric Electron Holes" and

orbcartpoin.f90 for the analysis reported in the article
"Finite gyro-radius multidimensional electron hole equilibria"

They may be used freely, without any warranty whatsoever, provided
those articles are cited in in any scholarly work which uses them. 

They share much code, the main difference being that orbitpoincare
uses cyclindrical coordinates while orbcartpoin uses cartesian
coordinates which performs better in the vicinity of the origin, and
is the current branch of development.

They can be built on a linux system with access to the development
libraries and headers of Xlib, by typing make. Running them with the
flag -h gives some terse help about the command line switches.
They can use the output of the code https://github.com/ihutch/helmhole
to provide the equilibrium. Examples of their usage are in the scripts
helmguardset.sh gyroplots.sh. But broadly they are not intended to be
used except by specialists who can read and understand the source fortran.

They are   Copyright (C) Ian H Hutchinson 2021.
