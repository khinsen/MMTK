RCS: $Id: README.make,v 1.2 1998/04/02 21:33:57 wrankin Exp $

Warning to HPPA users running the native 'make'.  And also,
apparently, Solaris users, and SGI users.  Well, heck, let's just say
"Anyone not using Gnu make".

If you do not have GNU Make (gmake) installed, then you will need to
edit 'src/Makefile' and 'test/Makefile' and add in explicit
linker/loader steps so that make will link all the '.o' files into an
executable.

An example is given in the Makefile for the 'dpmta_slave' executable.
For any executable which depends upon a bunch of ".o" files - for
example:

	dpmta_test : dpmta_test.o dpmta_testdump.o

there is (should be) an implicit compilation rule:

	$(CC) $(LDFLAGS) -o $@ $? $(LDLIBS)


If you do not like this, then get HP to fix their version of make
to handle object file dependencies correctly.  I tried introducing 
and adding implicit rules, but could not get 'make' to understand it.

If someone comes up with a solution other than hardcoding in the final
compilation steps (which I will *not* do), please send the fix to the
author.

-bill rankin
wrankin@ee.duke.edu


