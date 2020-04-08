#########################################################################
# This makefile should build on linux provided it has access to the X11 library
# which typically requires the development package, a fortran compiler,
# and git. It won't work on MSWindows.
# Decide the FORTRAN compiler and create the accis graphics routines:
include ACCIS.mk
#########################################################################
# This dependency should be included in the application makefile.
$(ACCISX) : $(ACCISHOME)Makefile
	@echo "$(ACCISCHECK)"
	cd $(ACCISHOME); make; cd -
#########################################################################
LIBRARIES := $(LIBRARIES)
LIBDEPS := $(LIBDEPS)
COMPILE-SWITCHES:=$(COMPILE-SWITCHES) -Wno-unused-dummy-argument
#########################################################################
MODULES=acpath.o
#########################################################################
# Patterns for compilation etc.
%.o : %.f makefile ;
	$(FORTRAN) -c $(COMPILE-SWITCHES) $*.f

%.o : %.f90 makefile;
	$(FORTRAN) -c $(COMPILE-SWITCHES) $*.f90

%.o : %.F makefile;
	$(FORTRAN) -c $(COMPILE-SWITCHES) $*.F

% : %.f $(ACCISX) $(LIBDEPS);
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.f $(LIBPATH) $(LIBRARIES)

% : %.f90  makefile $(ACCISX) $(MODULES) $(LIBDEPS);
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.f90 $(MODULES) $(LIBPATH) $(LIBRARIES)

% : %.F $(ACCISX) $(LIBDEPS);
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.F  $(LIBPATH) $(LIBRARIES)
#########################################################################
# A specific library of modified Bessel functions.
libmodbess.a : bessmodIs.f
	$(FORTRAN) -c bessmodIs.f
	ar -crs libmodbess.a bessmodIs.o

clean :
	rm -f *.o *.mod omarray tbedoc verifymain kpsiarray? fcontko dFtdWs bessmodsums omsolve fhgfuncmain libmodbess.a


