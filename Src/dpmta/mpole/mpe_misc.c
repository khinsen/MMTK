/*
 * mpe_misc.c - miscellaneous subroutines for multipole library
 *
 * w. elliott
 * w. t. rankin
 *
 * Copyright (c) 1995 Duke University
 * All Rights Reserved.
 *
 */

static char RCSid[] = "$Id: mpe_misc.c,v 1.5 1997/11/03 18:46:46 wrankin Exp $";

/*
 * RSC History:
 *
 * $Log: mpe_misc.c,v $
 * Revision 1.5  1997/11/03 18:46:46  wrankin
 * general cleanup/ansi-fication of code.  no new features.
 *
 * Revision 1.4  1997/09/26 04:28:00  wrankin
 * added routine to dump raw expansion data to file
 *
 * Revision 1.3  1997/05/09 20:15:00  wrankin
 * added routines to de-allocate global arrays created in [C,LJ]init()
 * added routines to free up multipole expansion matrices
 * added LJ prototypes and ansi-fied more procedures
 *
 * Revision 1.2  1996/03/06  21:47:25  wrankin
 * Implemented complete set of F/G multipole calculations
 *   - included evaluation of multipole (not local) potential
 *   - includes evaluation of multipole force.
 *   - converted local potential/force evaluation to use F/G format.
 *
 * Revision 1.1.1.1  1995/07/10  13:11:46  wrankin
 * Initial release of the Multipole Library
 * Based upon W. Elliott's MDMA codes
 * Implements Colomb Potentials only (no LJ-potentials yet)
 *
 *
 */

/*
 * include files
 */

#include <stdio.h>
#include <math.h>
#include "mpe.h"

     
/****************************************************************
 *
 *  Cart2Sph() - converts the cartesian vector 'v' into the
 *  spherical vector 's'
 *
 */
             
void Cart2Sph( Vector v, SphVector *s )
{
   Real x, y, z;
   Real rho, alpha, beta;

   x = v.x; y = v.y; z = v.z;

   rho = sqrt(x*x + y*y + z*z);

   if (rho < fabs(z)) rho = fabs(z);

   if (rho == 0.0) alpha = 0.0;
   else alpha = acos(z/rho);

   if ((x == 0.0) && (y == 0.0)) beta = 0.0;
   else beta = atan2(y, x);

   s->r = rho; s->a = alpha; s->b = beta;

} /* Cart2Sph */


/****************************************************************
 *
 *  dumpY_C() - print out a multipole expansion to standard error.
 *
 */

void dumpY_C( Mtype Y, int p )
{
   int n, m;

   for (n = 0; n < p; n++) {
      for (m = 0; m <= n; m++) {
	 fprintf(stderr, "%.6e %.6e   ", Y[n][m].x, Y[n][m].y);
      } /* for m */
      fprintf(stderr, "\n");
   } /* for n */
   fprintf(stderr, "\n\n");

} /* dumpY_C */


/****************************************************************
 *
 *  dumpYF() - print out an FFT multipole expansion to standard error.
 *
 */

void dumpYF( Real *Ystart, int p )

{
   int n, m, pf;

   pf = 1 << ((int) (log((double) (2 * p - 1)) / log(2.0)));

   for (n = 0; n < 2 * p; n++) {
      for (m = 0; m < pf; m++) {
	 fprintf(stderr, "%e\t%e\n", Ystart[2*n*pf+2*m], Ystart[2*n*pf+2*m+1]);
      } /* for m */
      fprintf(stderr, "\n");
   } /* for n */

} /* dumpYF */


/****************************************************************
 *
 *  MathdumpY_C() - dumps a multipole expansion to a file
 *    specified in s[].
 *
 */

void MathdumpY_C( Mtype Y, int p, char s[] )
{
   int             n, m;
   Real           *scratch;
   FILE           *fpo;

   fpo = fopen(s, "a");

   fprintf(fpo, "multfield = {\n");

   for (n = 0; n < p; n++) {
      fprintf(fpo, "{");
      for (m = 0; m < p; m++) {
	 if (m <= n) {
	    fprintf(fpo, "%.10e ", Y[n][m].x);
	    if (Y[n][m].y >= 0.0)
	       fprintf(fpo, "+ I %.10e ", Y[n][m].y);
	    else
	       fprintf(fpo, "- I %.10e ", -Y[n][m].y);
	 } /* if valid array value */ 
	 else
	    fprintf(fpo, "0");
	 if (m != p - 1)
	    fprintf(fpo, ", ");
      } /* for m */
      fprintf(fpo, "}");
      if (n < p - 1)
	 fprintf(fpo, ",");
      fprintf(fpo, "\n");
   } /* for n */
   fprintf(fpo, "}\n");

   fclose(fpo);

} /* MathdumpY */

/****************************************************************
 *
 *  MDumpRaw_C() - dumps a multipole expansion to a file
 *    specified in *s.
 *
 */

void MDumpRaw_C( Mtype Y, int p, char *s )
{
   int             n, m;
   FILE           *fpo;

   fpo = fopen(s, "w");

   for (n = 0; n < p; n++) {
      for (m = 0; m <=n; m++) {
	 fprintf(fpo, "%20.16lg %20.16lg\n", Y[n][m].x, Y[n][m].y);
      }
   }
   fclose(fpo);

} /* MathdumpY */

