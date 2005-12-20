# Top level Makefile to build DPMTA
#
# w. t. rankin
#
# Copyright (c) 1996 Duke University
# All rights reserved
#
# this makefile is set up to compile for a network of sun workstations
# running solaris 2.X.  the compile flags below wil need to be modifed
# for other platforms.
#
#
# RCS Info: $Id: Makefile,v 1.28 1998/05/05 15:44:54 wrankin Exp $
#
################################################################

# these can be overridden on the command line

MAKE   = make
SHELL  = /bin/sh
CC     = gcc
AR     = ar
RM     = rm -f
CP     = cp -f


# specify PVM envronmental variables.  PVM_HOME is the location
# where the executables and libraries will be placed.

#PVM_ROOT = /usr/local/pvm3
PVM_ROOT = /opt/pvm/pvm3

PVM_HOME = $(HOME)/pvm3


# set the compile time flags for the master and slave
# processes.  the available flags are:
# 
#    -DTIME        - include performance timing code
#    -DNOPARCONV   - disable parental conversion during the generation
#                    of the interaction list.
#    -DCOMP_LJ     - include the code to compute Lennard-Jones force and
#                    potentials
#    -DVIRIAL      - include code to compute the virial tensor.  this is
#                    experimental code and should not normally be compiled
#                    into production binaries.
#    -DOLDVIRIAL   - include code to compute the virial tensor.  this is
#                    the older and less accurate virial computation.
#    -DPBC_WRAP    - particles outside the box boundaries get wrapped
#                    around onto the other side of the box.
#                    this does not work with -DPIPED code.
#    -DSORTILIST   - sorts the order of the multipole interaction lists
#                    in attempt to increase accuacy.  best *not* to use.
#    -DMACROSCOPIC - include code that calculates the macroscopic assemblies.
#    -DPMACRO      - include parallel macroscopic generation code (must also
#                    define MACROSCOPIC above for this to work).
#    -DDATADUMP    - test programs dump particle and force
#                    information to output files.
#    -DPVMHOMO     - include optimizations for running on a homogeneous
#                    set of platforms.
#    -DSIZECHECK   - prints out the current size of the DPMTA memory
#                    images.  works for SERIAL case only.
#    -DTESTRESIZE  - test code resizes cell each iteration.
#    -DPIPED       - compile in code to allow parallel-piped simulation
#                    cells.  this code involves a significant performance
#                    hit, so it is not recommended unless you really
#                    need parallel-piped regions.
#    -DDEBUG       - causes each DPMTA process to dump data structure
#                    information.  this is for development purposes
#                    and should not be used by applications.
#    -DHILBERT     - use Hilbert ordering for dividing up cell ownership.
#                    this should result in better data division with
#                    a slight cost in serial performance.
#    -DROWCOL      - use Row/Column ordering for dividing up cell
#                    ownership.  This is provided for comparison and
#                    test purposes only and should not normally be used.
#                    Performance *will* suck.
#

DFLAGS = -DDATADUMP -DPVMHOMO -DTIME

################################################################
#
#  nothing below this line should need to be modified
#
################################################################

#default: solaris
default: linux


solaris:
	(cd mpole; $(MAKE) libmpole.a CC="$(CC)" )
	(cd src; $(MAKE) all install CC="$(CC)" MAKE=$(MAKE) \
           PVM_ARCH=SUN4SOL2 PVM_ROOT=$(PVM_ROOT) PVM_HOME=$(PVM_HOME) \
           DFLAGS="$(DFLAGS)" OSLIBS="-lsocket -lnsl" )

solx86:
	(cd mpole; $(MAKE) libmpole.a CC="$(CC)" )
	(cd src; $(MAKE) all install CC="$(CC)" MAKE=$(MAKE) \
           PVM_ARCH=X86SOL2 PVM_ROOT=$(PVM_ROOT) PVM_HOME=$(PVM_HOME) \
           DFLAGS="$(DFLAGS)" OSLIBS="-lsocket -lnsl" )

#
# for shared memory PVM
# 

sol-mp:
	(cd mpole; $(MAKE) libmpole.a CC="$(CC)" )
	(cd src; $(MAKE)  all install CC="$(CC)" MAKE=$(MAKE) \
           PVM_ARCH=SUNMP PVM_ROOT=$(PVM_ROOT) PVM_HOME=$(PVM_HOME) \
           DFLAGS="$(DFLAGS) -DPVMSHMEM" OSLIBS="-lsocket -lnsl" )

linux:
	(cd mpole; $(MAKE) libmpole.a CC="$(CC)" )
	(cd src; $(MAKE) all install CC="$(CC)" MAKE=$(MAKE) \
           PVM_ARCH=LINUX PVM_ROOT=$(PVM_ROOT) PVM_HOME=$(PVM_HOME) \
           DFLAGS="$(DFLAGS)" OSLIBS="" )

aix4:
	(cd mpole; $(MAKE) libmpole.a CC="$(CC)" \
	   DFLAGS="-mcpu=common" )
	(cd src; $(MAKE) all install CC="$(CC)" MAKE=$(MAKE) \
           PVM_ARCH=RS6K PVM_ROOT=$(PVM_ROOT) PVM_HOME=$(PVM_HOME) \
           DFLAGS="-mcpu=common $(DFLAGS)" OSLIBS="" )

#
#  note the need for gmake here
#

sgi:
	(cd mpole; $(MAKE) libmpole.a CC="$(CC)" \
	   CC=cc DFLAGS="-woff all" )
	(cd src; $(MAKE) all install CC="$(CC)" MAKE=$(MAKE) \
           PVM_ARCH=SGI64 PVM_ROOT=$(PVM_ROOT) PVM_HOME=$(PVM_HOME) \
           DFLAGS="$(DFLAGS) -woff all" OSLIBS="" )


sgimp:
	(cd mpole; $(MAKE) libmpole.a \
	   CC=cc DFLAGS="-woff all" )
	(cd src; $(MAKE) all install CC="$(CC)" MAKE=$(MAKE) \
           PVM_ARCH=SGIMP64 PVM_ROOT=$(PVM_ROOT) PVM_HOME=$(PVM_HOME) \
           DFLAGS="$(DFLAGS) -woff all" OSLIBS="" )

hppa:
	(cd mpole; $(MAKE) libmpole.a \
	   CC=cc DFLAGS="-Ae -Wp,-H500000" )
	(cd src; $(MAKE) all install CC="$(CC)" MAKE=$(MAKE) \
           PVM_ARCH=HPPA PVM_ROOT=$(PVM_ROOT) PVM_HOME=$(PVM_HOME) \
           DFLAGS="$(DFLAGS) -Ae -Wp,-H500000" OSLIBS="" )

#
# the following entry builds the embedded DPMTA library for
# use on the Cray T3D where you cannot have separate DPMTA and
# application executables.
#

t3d:
	(cd mpole; $(MAKE) libmpole.a \
	   CC="cc -T cray-t3d" )
	(cd src; $(MAKE) -f Makefile.mimd all \
           PVM_ARCH=CRAY PVM_ROOT="/usr" PVM_HOME=$(PVM_HOME) \
           CC=cc DFLAGS="-T cray-t3d $(DFLAGS)" OSLIBS="-lsma -lgpvm3" )
	(cd test; $(MAKE) mimd install-mimd \
           PVM_ARCH=CRAY PVM_ROOT="/usr" PVM_HOME=$(PVM_HOME) \
           DFLAGS="-T cray-t3d $(DFLAGS)" OSLIBS="-lsma -lgpvm3" \
	   PMTALIB="-ldpmta2" )


t3e:
	(cd mpole; $(MAKE) libmpole.a \
	   CC="cc -Xm" )
	(cd src; $(MAKE) -f Makefile.mimd all \
           PVM_ARCH=CRAY PVM_ROOT="/usr" PVM_HOME=$(PVM_HOME) \
           CC=cc DFLAGS="-Xm $(DFLAGS)" OSLIBS="-lsma" )
	(cd test; $(MAKE) mimd install-mimd \
           PVM_ARCH=CRAY PVM_ROOT="/usr" PVM_HOME=$(PVM_HOME) \
           DFLAGS="-Xm $(DFLAGS)" OSLIBS="-lsma" \
	   PMTALIB="-ldpmta2" )


#
# the following two build the cray-t3d style embedded library
# on workstation platforms.  these are used mainly for test.
#

sol-mimd:
	(cd mpole; $(MAKE) libmpole.a CC="$(CC)" )
	(cd src; $(MAKE) -f Makefile.mimd all CC="$(CC)" MAKE=$(MAKE) \
           PVM_ARCH=SUN4SOL2 PVM_ROOT=$(PVM_ROOT) PVM_HOME=$(PVM_HOME) \
	   PVM_LDIR=$(PVM_ROOT)/lib/SUN4SOL2 PVM_IDIR=$(PVM_ROOT)/include \
           DFLAGS="$(DFLAGS)" OSLIBS="-lsocket -lnsl -lgpvm3" CC=gcc )
	(cd test; $(MAKE) mimd install-mimd CC="$(CC)" MAKE=$(MAKE) \
           PVM_ARCH=SUN4SOL2 PVM_ROOT=$(PVM_ROOT) PVM_HOME=$(PVM_HOME) \
           DFLAGS="$(DFLAGS)" OSLIBS="-lsocket -lnsl" \
	   PMTALIB="-ldpmta2" )


linux-mimd:
	(cd mpole; $(MAKE) libmpole.a CC="$(CC)" )
	(cd src; $(MAKE) -f Makefile.mimd all CC="$(CC)" MAKE=$(MAKE) \
           PVM_ARCH=LINUX PVM_ROOT=$(PVM_ROOT) PVM_HOME=$(PVM_HOME) \
	   PVM_LDIR=$(PVM_ROOT)/lib/LINUX PVM_IDIR=$(PVM_ROOT)/include \
           DFLAGS="$(DFLAGS)" OSLIBS="-lgpvm3" )
	(cd test; $(MAKE) mimd install-mimd CC="$(CC)" MAKE=$(MAKE) \
           PVM_ARCH=LINUX PVM_ROOT=$(PVM_ROOT) PVM_HOME=$(PVM_HOME) \
           DFLAGS="$(DFLAGS)" OSLIBS="" \
	   PMTALIB="-ldpmta2" )

aix-mimd:
	(cd mpole; $(MAKE) libmpole.a CC="$(CC)" )
	(cd src; $(MAKE) -f Makefile.mimd all CC="$(CC)" MAKE=$(MAKE) \
           PVM_ARCH=RS6K PVM_ROOT=$(PVM_ROOT) PVM_HOME=$(PVM_HOME) \
	   PVM_LDIR=$(PVM_ROOT)/lib/RS6K PVM_IDIR=$(PVM_ROOT)/include \
           DFLAGS="$(DFLAGS)" OSLIBS="-lgpvm3" CC=gcc )
	(cd test; $(MAKE) mimd install-mimd CC="$(CC)" MAKE=$(MAKE) \
           PVM_ARCH=RS6K PVM_ROOT=$(PVM_ROOT) PVM_HOME=$(PVM_HOME) \
           DFLAGS="$(DFLAGS)" OSLIBS="" \
	   PMTALIB="-ldpmta2" )

#
# here is the code which compiles DPMTA as a serial application
#

serial:
	(cd mpole; $(MAKE) libmpole.a CC="$(CC)")
	(cd src; $(MAKE) -f Makefile.serial all CC="$(CC)" MAKE=$(MAKE) \
           DFLAGS="$(DFLAGS)" )

linux-ser:
	(cd mpole; $(MAKE) libmpole.a CC="$(CC)")
	(cd src; $(MAKE) -f Makefile.serial all CC="$(CC)" MAKE=$(MAKE) \
	   DFLAGS="$(DFLAGS) -DLINUX" )

aix-ser:
	(cd mpole; $(MAKE) libmpole.a CC="$(CC)")
	(cd src; $(MAKE) -f Makefile.serial all CC="$(CC)" MAKE=$(MAKE) \
	   DFLAGS="$(DFLAGS) -DRS6K" )

sgi-ser:
	(cd mpole; $(MAKE) libmpole.a CC="$(CC)" DFLAGS="-woff all" )
	(cd src; $(MAKE) -f Makefile.serial all CC="$(CC)" MAKE=$(MAKE) \
	   DFLAGS="$(DFLAGS) -DSGI64 -woff all" )

clean: 
	( cd mpole; $(MAKE) clean )
	( cd src; $(MAKE) clean )

paper: 
	( cd doc; $(MAKE) )

