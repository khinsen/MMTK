# 
# Makefile to build multipole library
#
# w. t. rankin
#
# Copyright (c) 1995 Duke University
# All rights reserved
#
# RCS Info: $Id: Makefile,v 1.11 1998/04/02 21:34:10 wrankin Exp $
#
# Revision History:
# $Log: Makefile,v $
# Revision 1.11  1998/04/02 21:34:10  wrankin
# added Make.dep files to contain all file dependencies
# added sol86 (Solaris X86) architecture to top level Makefile
# cleaned up propagation of CC and MAKE variables in Makefiles
#
# Revision 1.10  1997/11/03 18:46:38  wrankin
# general cleanup/ansi-fication of code.  no new features.
#
# Revision 1.9  1997/05/13 17:52:11  wrankin
# fixed irritation bug seen when math.h is not included
# added a missing prototype to mpe.h
#
# Revision 1.8  1997/05/09  20:14:55  wrankin
# added routines to de-allocate global arrays created in [C,LJ]init()
# added routines to free up multipole expansion matrices
# added LJ prototypes and ansi-fied more procedures
#
# Revision 1.7  1996/10/18  17:00:51  wrankin
# implement new heirarchical make structure
#
# Revision 1.6  1996/09/24  19:05:48  wrankin
# added more test codes.
#
# Revision 1.5  1996/09/19  18:38:02  wrankin
# Makefile now leaves libmpole in current directory
#
# Revision 1.4  1995/09/15  15:08:51  wrankin
# Added Lennard-Jones multipole routines to library.
#
# Revision 1.3  1995/07/19  01:55:14  wrankin
# updates to reflect move of multipole library code to a new directory
#
# Revision 1.2  1995/07/17 01:31:43  wrankin
# Makefile updates to support SGI Power Challenge platform
#
# Revision 1.1.1.1  1995/07/10  13:11:46  wrankin
# Initial release of the Multipole Library
# Based upon W. Elliott's MDMA codes
# Implements Colomb Potentials only (no LJ-potentials yet)
#
#
################################################################


MAKE   = make
SHELL  = /bin/sh
CC     = gcc
AR     = ar
RM     = rm -f
CP     = cp

#  for the HPPA machines, add "-Ae" and "-Wp,-H500000" to CFLAGS
#  for the SGI machines, add "-Woff all" to CFLAGS

CFLAGS = -O3 $(DFLAGS)
LDFLAGS =
LDLIBS = -lm

SRC  = mpe_fft.c mpe_misc.c mpe_mpoleC.c mpe_allocC.c \
       mpe_mpoleLJ.c mpe_allocLJ.c
OBJS = $(SRC:.c=.o)

TSRC = mpe_testC.c mpe_testLJ.c
TOBJS = $(TSRC:.c=.o)

################################################################

default : all

all : libmpole.a

libmpole.a : $(OBJS)
	$(RM) libmpole.a
	$(AR) rs libmpole.a $(OBJS)

install : all

clean:
	$(RM) libmpole.a $(OBJS) $(TOBJS)

test : testC testLJ

testLJ : mpe_testLJ.o libmpole.a
	$(CC) $(LDFLAGS) -o $@ $? $(LDLIBS)

testC : mpe_testC.o libmpole.a
	$(CC) $(LDFLAGS) -o $@ $? $(LDLIBS)


################################################################
#
#  include file dependencies
#

include Make.deps
