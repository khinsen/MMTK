/*
 *  mpe_mpoleC.c - multipole routines for computing  1/r interactions
 *
 *  w. t. rankin
 *  w. elliott
 *
 *  Copyright (c) 1995 Duke University
 *  All Rights Reserved
 */

static char RCSid[] = "$Id: mpe_mpoleC.c,v 1.13 1997/11/03 18:46:48 wrankin Exp $";

/*
 * RSC History:
 *
 * $Log: mpe_mpoleC.c,v $
 * Revision 1.13  1997/11/03 18:46:48  wrankin
 * general cleanup/ansi-fication of code.  no new features.
 *
 * Revision 1.12  1997/09/26 04:28:01  wrankin
 * added routine to dump raw expansion data to file
 *
 * Revision 1.11  1997/05/09 20:15:02  wrankin
 * added routines to de-allocate global arrays created in [C,LJ]init()
 * added routines to free up multipole expansion matrices
 * added LJ prototypes and ansi-fied more procedures
 *
 * Revision 1.10  1996/11/19  22:13:41  wrankin
 * implemented more efficient MCM routine
 *
 * Revision 1.9  1996/11/14  17:49:53  wrankin
 * performance enhancements to multipole M2L routine.
 * additions to make slaves exit gracefully.
 *
 * Revision 1.8  1996/11/11  20:09:11  wrankin
 * added ANSI-C declarations to mpoleC and prototypes to mpe.h
 *
 * Revision 1.7  1996/03/06  21:47:27  wrankin
 * Implemented complete set of F/G multipole calculations
 *   - included evaluation of multipole (not local) potential
 *   - includes evaluation of multipole force.
 *   - converted local potential/force evaluation to use F/G format.
 *
 * Revision 1.6  1996/02/29  20:54:29  wrankin
 * size of FFT MPE was being set wrong.
 *
 * Revision 1.5  1996/02/12  15:27:30  wrankin
 * Added functions to support Macroscopic Assemblies.
 *
 * Revision 1.4  1996/01/29  21:07:54  wrankin
 * Clean up multipole code.
 *   - removed unused procedures.
 *   - added headers and comments
 *
 * Revision 1.3  1995/12/08  22:58:55  wrankin
 * preliminary release of DPMTA 2.3
 *   - added working Periodic Boundary Conditions (PDC) flag
 *   - added code for Virial Computation (not working yet)
 *   - general cleanup of some modules
 *
 * Revision 1.2  1995/10/01  21:34:57  wrankin
 * added LJsum() function to add multipoles
 * changed name of LJsize() funtion to match Csize()
 * general code cleanup
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpe.h"
#include "mpe_legendre.h"
#include "mpe_fftC.h"


/*
 *  local arrays used in multipole calculations
 */

static Mtype   Y_C;
static Mtype   Hm2l;
static Mtype   L;
static Complex *Yxy;
static Real    **Ycoeff;
static Real    **Fcoeff;
static Real    **Gcoeff;
static Real    **A_C;
static Real    **LegPoly;


/****************************************************************
 *
 *  AddMultipoleC() calculates the Coulomb multipole expansion of a
 *    particle of charge q, at location x,y,z relative to the origin.
 *    The result is accumulated into the blank multipole expansion
 *    array M. The array size is p.
 *
 *    this computation is performed using the simplified taylor series
 *    representation as described by elliott's dissertation.
 *
 *    see TR95-003 Sec. 2.1.3 for more information.
 *
 */

int AddMultipoleC(
   Mtype  M,
   int    p,
   Real   q,
   Vector v )
{
   int             n, m;
   SphVector       sv;

   /*
    * convert input to a spherical vector
    */

   Cart2Sph(v,&sv);

   /*
    * compute the F matrix, leaving the result in Y_C{n][m]
    */

   makeF(p, sv);

   for (n = 0; n < p; n++) {
      for (m = 0; m <= n; m++) {
	 M[n][m].x += q * Y_C[n][m].x;
	 M[n][m].y -= q * Y_C[n][m].y;
      }	/* for m */
   } /* for n */
} /* AddMultipoleC */



/****************************************************************
 *
 *  eval_mpotC() - computes the potential field due to the
 *    multipole, M, at the point <v> from the center of the 
 *    multipole expansion.
 *
 *    this computation is performed by taking the inner product
 *    of the M and G matrices, using the simplified taylor series
 *    representation as described by elliott's dissertation.
 *
 *    see TR95-003 Sec. 2.1.3 for more information.
 *
 */

Real eval_mpotC(
   Mtype   M,
   int     p,
   Vector  v )
{
   int             n, m;
   Real            pot;
   SphVector       sv;


   /*
    * convert input to a spherical vector
    */

   Cart2Sph(v,&sv);

   /*
    * compute the G matrix, leaving the result in Y_C{n][m]
    */

   makeG(p, sv);

   pot = 0;
   for (n = 0; n < p; n++) {
      for (m = 0; m<=n; m++) {
         if (m == 0) {
            pot += Y_C[n][m].x * M[n][m].x - Y_C[n][m].y * M[n][m].y;
	 } /* if */
         else {
            pot += 2.0 * ( Y_C[n][m].x * M[n][m].x - Y_C[n][m].y * M[n][m].y );
	 } /* else */
      } /* for m */
   }/* for n */

   return pot;

} /* eval_mpotC */


/****************************************************************
 *
 *  eval_lpotC() - computes the potential field due to the
 *    local expansion, M, at the point <v> from the center of the 
 *    local expansion.
 *
 *    this computation is performed by taking the inner product
 *    of the L and F* matrices, using the simplified taylor series
 *    representation as described by elliott's dissertation.
 *
 *    see TR95-003 Sec. 2.1.3 for more information.
 *
 */

Real eval_lpotC(
   Mtype   M,
   int     p,
   Vector  v )
{
   int             n, m;
   Real            pot;
   SphVector       sv;


   /*
    * convert input to a spherical vector
    */

   Cart2Sph(v,&sv);

   /*
    * compute the F matrix, leaving the result in Y_C{n][m]
    */

   sv.r = -sv.r;
   makeF(p, sv);

   /*
    * compute the product of L[n,m].F*[n,m]
    */

   pot = 0;
   for (n = 0; n < p; n++) {
      for (m = 0; m<=n; m++) {
         if (m == 0) {
            pot += Y_C[n][m].x * M[n][m].x + Y_C[n][m].y * M[n][m].y;
	 } /* if */
         else {
            pot += 2.0 * ( Y_C[n][m].x * M[n][m].x + Y_C[n][m].y * M[n][m].y );
	 } /* else */
      } /* for m */
   }/* for n */

   return pot;

} /* eval_lpotC */


/****************************************************************
 *
 * M2M_C translates a Coulomb multipole expansion M1 to another
 * expansion M2, both size p, along the vector v
 *
 */

int M2M_C(
   Mtype   M1,
   Mtype   M2,
   int     p,
   Vector  v )
{
   int             n, m, np, mp, startm, endm;
   Real            atemp;
   SphVector       sv;

   /*
    * convert input to a spherical vector
    */

   Cart2Sph(v,&sv);

   /*
    * compute the F matrix, leaving the result in Y_C{n][m]
    */

   sv.r = -sv.r;
   makeF(p, sv);

   for (np = 0; np < p; np++) {
      for (mp = 0; mp <= np; mp++) {
	 for (n = 0; n <= np; n++) {
	    startm = mp - (np - n);
	    startm = (startm < -n ? -n : startm);
	    endm = mp + (np - n);
	    endm = (endm > n ? n : endm);
	    if (startm <= endm) {
	       m = startm;
	       while (m < 0 && m <= endm) {
		  atemp = 1.0 - 2.0 * (Real) (0x0001 & (-m));
		  /* above eqn = -1^m */
		  M2[np][mp].x += atemp *
		     (M1[n][-m].x * Y_C[np - n][mp - m].x -
		      M1[n][-m].y * Y_C[np - n][mp - m].y);
		  M2[np][mp].y -= atemp *
		     (M1[n][-m].x * Y_C[np - n][mp - m].y +
		      M1[n][-m].y * Y_C[np - n][mp - m].x);
		  m++;
	       }		/* while m negative */
	       while (m < mp && m <= endm) {
		  M2[np][mp].x +=
		     (M1[n][m].x * Y_C[np - n][mp - m].x +
		      M1[n][m].y * Y_C[np - n][mp - m].y);
		  M2[np][mp].y +=
		     (M1[n][m].x * -Y_C[np - n][mp - m].y +
		      M1[n][m].y * Y_C[np - n][mp - m].x);
		  m++;
	       } /* while m, positive, less than mp */
	       while (m <= endm) {
		  atemp = 1.0 - 2.0 * (Real) (0x0001 & (mp + m));
		  M2[np][mp].x += atemp *
		     (M1[n][m].x * Y_C[np - n][m - mp].x -
		      M1[n][m].y * Y_C[np - n][m - mp].y);
		  M2[np][mp].y += atemp *
		     (M1[n][m].x * Y_C[np - n][m - mp].y +
		      M1[n][m].y * Y_C[np - n][m - mp].x);
		  m++;
	       } /* while m */
	    } /* if startm <= endm (ie m exists) */
	 } /* for n */
      } /* for mp */
   } /* for np */

   return TRUE;

} /* M2M_C */


/****************************************************************
 *
 * M2M_C translates a Coulomb multipole expansion M1 to another
 * expansion M2, both size p, along the vector v
 *
 */

int M2M_Cshort(
   Mtype   M1,
   Mtype   M2,
   Mtype   H,
   int     p )
{
   int             n, m, np, mp, startm, endm;
   Real            atemp;


   for (np = 0; np < p; np++) {
      for (mp = 0; mp <= np; mp++) {
	 for (n = 0; n <= np; n++) {
	    startm = mp - (np - n);
	    startm = (startm < -n ? -n : startm);
	    endm = mp + (np - n);
	    endm = (endm > n ? n : endm);
	    if (startm <= endm) {
	       m = startm;
	       while (m < 0 && m <= endm) {
		  atemp = 1.0 - 2.0 * (Real) (0x0001 & (-m));
		  /* above eqn = -1^m */
		  M2[np][mp].x += atemp *
		     (M1[n][-m].x * H[np - n][mp - m].x -
		      M1[n][-m].y * H[np - n][mp - m].y);
		  M2[np][mp].y -= atemp *
		     (M1[n][-m].x * H[np - n][mp - m].y +
		      M1[n][-m].y * H[np - n][mp - m].x);
		  m++;
	       }		/* while m negative */
	       while (m < mp && m <= endm) {
		  M2[np][mp].x +=
		     (M1[n][m].x * H[np - n][mp - m].x +
		      M1[n][m].y * H[np - n][mp - m].y);
		  M2[np][mp].y +=
		     (M1[n][m].x * -H[np - n][mp - m].y +
		      M1[n][m].y * H[np - n][mp - m].x);
		  m++;
	       } /* while m, positive, less than mp */
	       while (m <= endm) {
		  atemp = 1.0 - 2.0 * (Real) (0x0001 & (mp + m));
		  M2[np][mp].x += atemp *
		     (M1[n][m].x * H[np - n][m - mp].x -
		      M1[n][m].y * H[np - n][m - mp].y);
		  M2[np][mp].y += atemp *
		     (M1[n][m].x * H[np - n][m - mp].y +
		      M1[n][m].y * H[np - n][m - mp].x);
		  m++;
	       } /* while m */
	    } /* if startm <= endm (ie m exists) */
	 } /* for n */
      } /* for mp */
   } /* for np */

   return TRUE;

} /* M2M_Cshort */




/****************************************************************
 *
 * M2L_Cshort converts a Coulomb multipole expansion M in to a local
 * expansion L, shifting the original expansion using the precomputed
 * transfer matrix H. Both expansions are of size p
 *
 */

int M2L_Cshort(
   Mtype   M,
   Mtype   L,
   Mtype   H,
   int     p )
{
   int             n, m, np, mp;
   Real            atemp, btemp;
   Complex         *Mp, *Lp, *Hp;

   Lp = &(L[0][0]);
   for (np = 0; np < p; np++) {
      atemp = 1.0;
      for (mp = 0; mp <= np; mp++) {
	 for (n = 0; n < p - np; n++) {
	    Mp = &(M[n][n]);
	    Hp = &(H[np+n][n-mp]);
	    for (m = -n; m < -mp; m++) {
	       Lp->x += atemp * (Mp->x * Hp->x - Mp->y * Hp->y);
	       Lp->y -= atemp * (Mp->x * Hp->y + Mp->y * Hp->x);
	       Mp--;
	       Hp--;
	    } /* for m, m neg and |m| > |mp| */
	    Hp = &(H[np+n][mp+m]);
            btemp = 1.0 - 2.0 * (Real) (0x0001 & (-m));
	    for (; m < 0; m++) {
	       Lp->x += btemp * (Mp->x * Hp->x + Mp->y * Hp->y);
	       Lp->y += btemp * (Mp->x * Hp->y - Mp->y * Hp->x);
	       Mp--;
	       Hp++;
	       btemp = -btemp;
            } /* for m, m neg and |m| <= |mp| */
	    for (; m <= n; m++) {
	       Lp->x += (Mp->x * Hp->x - Mp->y * Hp->y);
	       Lp->y += (Mp->x * Hp->y + Mp->y * Hp->x);
	       Mp++;
	       Hp++;
	    } /* for m, m neg and mp-m positive */
	 } /* for n */

	 Lp++;
	 atemp = -atemp;
      } /* for mp */
   } /* for np */

   return TRUE;

} /* M2L_Cshort */



/****************************************************************
 *
 * M2L_Cshort_Old() -
 *
 * converts a Coulomb multipole expansion M in to a local
 * expansion L, shifting the original expansion using the precomputed
 * transfer matrix H. Both expansions are of size p
 *
 * This is an older version of the M2L_Cshort() function which uses
 * less efficient explicit array indexing rather than direct pointer
 * manipulations.  The results should be identical to M2L_Cshort().
 *
 */

int M2L_Cshort_Old(
   Mtype   M,
   Mtype   L,
   Mtype   H,
   int     p )
{
   int             n, m, np, mp;
   Real            atemp;


   for (np = 0; np < p; np++) {
      for (mp = 0; mp <= np; mp++) {
	 for (n = 0; n < p - np; n++) {
	    for (m = -n; m < -mp; m++) {
	       atemp = 1.0 - 2.0 * (Real) (0x0001 & mp);
	       L[np][mp].x += atemp *
		  (M[n][-m].x * H[np + n][-(mp + m)].x -
		   M[n][-m].y * H[np + n][-(mp + m)].y);
	       L[np][mp].y -= atemp *
		  (M[n][-m].x * H[np + n][-(mp + m)].y +
		   M[n][-m].y * H[np + n][-(mp + m)].x);
	    } /* for m, m neg and |m| > |mp| */
	    for (; m < 0; m++) {
	       atemp = 1.0 - 2.0 * (Real) (0x0001 & (-m));
	       L[np][mp].x += atemp *
		  (M[n][-m].x * H[np + n][mp + m].x +
		   M[n][-m].y * H[np + n][mp + m].y);
	       L[np][mp].y += atemp *
		  (M[n][-m].x * H[np + n][mp + m].y -
		   M[n][-m].y * H[np + n][mp + m].x);
            } /* for m, m neg and |m| <= |mp| */
	    for (; m <= n; m++) {
	       L[np][mp].x +=
		  (M[n][m].x * H[np + n][mp + m].x -
		   M[n][m].y * H[np + n][mp + m].y);
	       L[np][mp].y +=
		  (M[n][m].x * H[np + n][mp + m].y +
		   M[n][m].y * H[np + n][mp + m].x);
	    } /* for m, m neg and mp-m positive */
	 } /* for n */
      } /* for mp */
   } /* for np */

   return TRUE;

} /* M2L_Cshort */


/****************************************************************
 *
 * M2L_C() -
 *
 * converts a Coulomb multipole expansion M in to a local
 * expansion L, shifting the original expansion along the vector v.
 * Both expansions are of size p
 *
 */

int M2L_C(
   Mtype   M,
   Mtype   L,
   int     p,
   Vector  v )
{

   int             n, m, np, mp;
   Real            atemp;
   SphVector       sv;


   /*
    * convert input to a spherical vector
    */

   Cart2Sph(v,&sv);

   /*
    * compute the G matrix, leaving the result in Y_C[n][m]
    */

   makeG(p, sv);


   for (np = 0; np < p; np++) {
      for (mp = 0; mp <= np; mp++) {
	 for (n = 0; n < p - np; n++) {
	    for (m = -n; m < -mp; m++) {
	       atemp = 1.0 - 2.0 * (Real) (0x0001 & mp);
	       L[np][mp].x += atemp *
		  (M[n][-m].x * Y_C[np + n][-(mp + m)].x -
		   M[n][-m].y * Y_C[np + n][-(mp + m)].y);
	       L[np][mp].y -= atemp *
		  (M[n][-m].x * Y_C[np + n][-(mp + m)].y +
		   M[n][-m].y * Y_C[np + n][-(mp + m)].x);
	    } /* for m, m neg and |m| > |mp| */
	    for (; m < 0; m++) {
	       atemp = 1.0 - 2.0 * (Real) (0x0001 & (-m));
	       L[np][mp].x += atemp *
		  (M[n][-m].x * Y_C[np + n][mp + m].x +
		   M[n][-m].y * Y_C[np + n][mp + m].y);
	       L[np][mp].y += atemp *
		  (M[n][-m].x * Y_C[np + n][mp + m].y -
		   M[n][-m].y * Y_C[np + n][mp + m].x);
	    } /* for m, m neg and |m| <= |mp| */
	    for (; m <= n; m++) {
	       L[np][mp].x +=
		  (M[n][m].x * Y_C[np + n][mp + m].x -
		   M[n][m].y * Y_C[np + n][mp + m].y);
	       L[np][mp].y +=
		  (M[n][m].x * Y_C[np + n][mp + m].y +
		   M[n][m].y * Y_C[np + n][mp + m].x);
	    } /* for m, m neg and mp-m positive */
	 } /* for n */
      } /* for mp */
   } /* for np */

   return TRUE;

} /* M2L_C */


/****************************************************************
 *
 *  L2L_C translates a local expansion L1 to another local 
 *  expansion L2, both of size p, along vector v.
 *
 */

int L2L_C(
   Mtype   L1,
   Mtype   L2,
   int     p,
   Vector  v )
{
   int             np, mp, n, m;
   int             startm, endm;
   Real            atemp;
   SphVector       sv;

   /*
    * convert input to a spherical vector
    */

   Cart2Sph(v,&sv);

   /*
    * compute the F matrix, leaving the result in Y_C{n][m]
    */

   sv.r = -sv.r;
   makeF(p, sv);


   for (np = 0; np < p; np++) {
      for (mp = 0; mp <= np; mp++) {
	 for (n = np; n < p; n++) {
	    startm = mp - (n - np);
	    startm = (startm < -n ? -n : startm);
	    endm = mp + (n - np);
	    endm = (endm > n ? n : endm);
	    m = startm;
	    if (startm <= endm) {
	       while (m < 0 && m <= endm) {
		  atemp = 1.0 - 2.0 * (Real) (0x0001 & (mp));
		  L2[np][mp].x += atemp *
		     (L1[n][-m].x * Y_C[n - np][-m + mp].x +
		      L1[n][-m].y * Y_C[n - np][-m + mp].y);
		  L2[np][mp].y += atemp *
		     (L1[n][-m].x * Y_C[n - np][-m + mp].y -
		      L1[n][-m].y * Y_C[n - np][-m + mp].x);
		  m++;
	       } /* while m is negative */
	       while (m < mp && m <= endm) {
		  atemp = 1.0 - 2.0 * (Real) (0x0001 & (m + mp));
		  L2[np][mp].x += atemp *
		     (L1[n][m].x * Y_C[n - np][-m + mp].x -
		      L1[n][m].y * Y_C[n - np][-m + mp].y);
		  L2[np][mp].y += atemp *
		     (L1[n][m].x * Y_C[n - np][-m + mp].y +
		      L1[n][m].y * Y_C[n - np][-m + mp].x);
		  m++;
	       } /* while m-mp is negative */
	       while (m <= endm) {
		  L2[np][mp].x +=
		     (L1[n][m].x * Y_C[n - np][m - mp].x +
		      L1[n][m].y * Y_C[n - np][m - mp].y);
		  L2[np][mp].y +=
		     (L1[n][m].x * -Y_C[n - np][m - mp].y +
		      L1[n][m].y * Y_C[n - np][m - mp].x);
		  m++;
	       } /* while mp < endm, m and m-mp positive */
	    } /* if startm <= endm (ie m exists) */
	 } /* for n */
      } /* for mp */
   } /* for np */

   return TRUE;

} /* L2L_C */


/****************************************************************
 *
 * Force_C_Y computes the Coulomb force and potential on the particle
 * with charge q at position x, y, z, due to the local expansion L of
 * size p. The results are put in the array result in the following
 * order: pot., fx, fy, fz.
 *
 * note: this is the old version that used the Ynm matrix.  It has
 * been replaced by the more recent G/F form of the equation.
 *
 */

int Force_C_Y  (
   Mtype   Lin,
   int     p,
   Real    q,
   Vector  pv,
   Real    *potp,
   Vector  *fv )
{
   int             n, m;
   Real            rho, alpha, beta;
   Real            pot, fr, fa, fb;
   Real            cosalpha, sinalpha, cosbeta, sinbeta, rtemp;
   Real            neg_one;
   Vector          vect;
   SphVector       svect;

   rho = sqrt( pv.x*pv.x + pv.y * pv.y + pv.z * pv.z);
   if (rho < fabs(pv.z))
      rho = fabs(pv.z);
   alpha = acos(pv.z / rho);
   if ((pv.x == 0.0) && (pv.y == 0.0)) {
      beta = 0.0;
   } else {
      beta = atan2(pv.y, pv.x);
   }

   cosalpha = cos(alpha);
   sinalpha = sin(alpha);
   if (sinalpha < SMALL_THETA) {
      sinalpha = 0.0;
      if (cosalpha > 0.0)
	 cosalpha = 1.0;
      else
	 cosalpha = -1.0;
   }				/* if a small alpha */
   cosbeta = cos(beta);
   sinbeta = sin(beta);

   pot = 0.0;
   fr = 0.0;
   fa = 0.0;
   fb = 0.0;

   makeYforceC(p, rho, alpha, beta);

   for (n = 0; n < p; n++) {
      for (m = 0; m <= n; m++) {
	 L[n][m].x = A_C[n][m] * Lin[n][m].x;
	 L[n][m].y = -A_C[n][m] * Lin[n][m].y;
      } /* for m */
   } /* for n */

   pot = Y_C[0][0].x * L[0][0].x;
   for (n = 1; n < p; n++) {
      pot += rho * Y_C[n][0].x * L[n][0].x;
      for (m = 1; m <= n; m++) {
	 pot += 2.0 * rho * (Y_C[n][m].x * L[n][m].x 
             - Y_C[n][m].y * L[n][m].y);
      } /* for m */
   } /* for n */

   rtemp = 1.0;
   for (n = 1; n < p; n++) {
      fr += (Real) (-n) * L[n][0].x * Y_C[n][0].x;
      fa += rtemp * Ycoeff[n][0] * L[n][0].x * (-LegPoly[n][1]);
      if (sinalpha != 0.0) {
	 for (m = 1; m < n; m++) {
	    fr += 2.0 * (Real) (-n) * (Y_C[n][m].x * L[n][m].x -
				       Y_C[n][m].y * L[n][m].y);
	    fa += 2.0 * rtemp * Ycoeff[n][m] *
	       (L[n][m].x * Yxy[m].x - L[n][m].y * Yxy[m].y) *
	       ((Real) (-m) * cosalpha / sinalpha * LegPoly[n][m] -
		LegPoly[n][m + 1]);
	    fb += 2.0 / sinalpha * (Real) m *(L[n][m].x * Y_C[n][m].y +
				                   L[n][m].y * Y_C[n][m].x);
	 } /* for m */
	 fr += 2.0 * (Real) (-n) * (Y_C[n][m].x * L[n][m].x -
				    Y_C[n][m].y * L[n][m].y);
	 fa += 2.0 * rtemp * Ycoeff[n][m] *
	    (L[n][m].x * Yxy[m].x - L[n][m].y * Yxy[m].y) *
	    ((Real) (-m) * cosalpha / sinalpha * LegPoly[n][m]);
	 fb += 2.0 / sinalpha * (Real) m *(L[n][m].x * Y_C[n][m].y +
				                   L[n][m].y * Y_C[n][m].x);
      } /* if */ 

      else {
	 neg_one = 1.0 + (cosalpha - 1.0) * (Real) (0x0001 & (n + 1));
	 /* neg one = (-1^(n+1)) if cosalpha = -1, otherwise = 1 */
	 fa += 2.0 * rtemp * Ycoeff[n][1] * neg_one *
	    (L[n][1].x * Yxy[1].x - L[n][1].y * Yxy[1].y) *
	    cosalpha * ((Real) ((n + 1) * n) / 2.0);
	 fb += (L[n][1].x * Yxy[1].y + L[n][1].y * Yxy[1].x) *
	    neg_one * rtemp * Ycoeff[n][1] * (Real) ((n + 1) * -n);
	 for (m = 1; m <= n; m++) {
	    fr += 2.0 * (Real) (-n) * (Y_C[n][m].x * L[n][m].x -
				       Y_C[n][m].y * L[n][m].y);
	 } /* for m */
      } /* else */
      rtemp *= rho;
   } /* for n */


   /*
    *  we had to reverse the signs on the force vectors.
    *  ask bill why we did this?
    */

   *potp = q * pot;
   fv->x = q * (fr*sinalpha*cosbeta + fa*cosalpha*cosbeta - fb*sinbeta);
   fv->y = q * (fr*sinalpha*sinbeta + fa*cosalpha*sinbeta + fb*cosbeta);
   fv->z = q * (fr*cosalpha - fa*sinalpha);

} /* Force_C_Y */


/****************************************************************
 *
 *  Force_C computes the Coulomb force and potential on the
 *  particle with charge q at position x, y, z, due to the local
 *  expansion L of size p.
 *
 */

int Force_C(
   Mtype   Lin,
   int     p,
   Real    q,
   Vector  pv,
   Real    *potp,
   Vector  *fv )
{
   int             n, m;
   Real            pot, fr, fa, fb;
   Real            cosalpha, sinalpha, cosbeta, sinbeta;
   Real            cotalpha;
   Real            rtemp, ntemp, mtemp, ftemp, ptemp;
   SphVector       sv;

   /*
    * convert input to a spherical vector
    */

   Cart2Sph(pv,&sv);

   cosalpha = cos(sv.a);
   sinalpha = sin(sv.a);
   if (sinalpha < SMALL_THETA) {
      sinalpha = 0.0;
      if (cosalpha > 0.0)
	 cosalpha = 1.0;
      else
	 cosalpha = -1.0;
   }				/* if a small alpha */
   cosbeta = cos(sv.b);
   sinbeta = sin(sv.b);
   cotalpha = cosalpha/sinalpha;

   pot = 0.0;
   fr = 0.0;
   fa = 0.0;
   fb = 0.0;

   /*
    * compute the F matrix, leaving the result in Y_C{n][m]
    */

   sv.r = -sv.r;
   makeF(p, sv);


   /*
    * find the potential by computing the product of L[n,m].F*[n,m]
    */

   for (n = 0; n < p; n++) {
      for (m = 0; m<=n; m++) {
         if (m == 0) {
            pot += Y_C[n][m].x*Lin[n][m].x + Y_C[n][m].y*Lin[n][m].y;
	 } /* if */
         else {
            pot += 2.0 * (Y_C[n][m].x*Lin[n][m].x + Y_C[n][m].y*Lin[n][m].y);
	 } /* else */
      } /* for m */
   }/* for n */


   /*
    *  this section computes the special case for the
    *  condition where either alpha=0 (assume beta=0)
    *  or rho=0 (assume alpha=0 and beta=0)
    *
    */

   if ( sinalpha == 0.0 ) {

      if (sv.r == 0.0) {

         fr += Fcoeff[1][0] * LegPoly[1][0] * Lin[1][0].x;
         fr += 2.0 * Fcoeff[1][1] * LegPoly[1][1] * Lin[1][1].x;

         /*
	  * rotate the r-vector to the x-axis and recompute
	  * we can probably get away with only computing the
	  * first couple terms in the Legerdre polynomial
	  */

	 Legendre(LegPoly, p, 0.0);

         fa += Fcoeff[1][0] * LegPoly[1][0] * Lin[1][0].x;
         fa += 2.0 * Fcoeff[1][1] * LegPoly[1][1] * Lin[1][1].x;

         fb += Fcoeff[1][0] * LegPoly[1][0] * Lin[1][0].y;
         fb += 2.0 * Fcoeff[1][1] * LegPoly[1][1] * Lin[1][1].y;

      } /* sv.r == 0 */


      /*
       * otherwise, we can compute the vectors in a more simple
       * fashion.
       */

      else {
         rtemp = 1.0 / sv.r;
	 ftemp = rtemp;
         for (n = 1; n < p; n++) {
            fr += ftemp * (Y_C[n][0].x*Lin[n][0].x + Y_C[n][0].y*Lin[n][0].y);
            for (m = 1; m<=n; m++) {
               fr += 2.0 * ftemp *
		  (Y_C[n][m].x*Lin[n][m].x + Y_C[n][m].y*Lin[n][m].y);
            } /* for m */
	    ftemp += rtemp;
         } /* for n */

         /*
          *  compute the alpha and beta component of the gradient
          */

         rtemp = 1.0;
         ntemp = 0.0;
         for (n = 1; n < p; n++) {
            ntemp -= (double)n;
            ftemp  = 2.0 * Fcoeff[n][1] * rtemp * ntemp;
            fa += ftemp * Lin[n][1].x;
            fb += ftemp * Lin[n][1].y;
            rtemp *= sv.r;
         } /* for n */

      } /* else */

      *potp = q * pot;
      fv->x = q * fa;
      fv->y = q * fb;
      fv->z = q * fr;

   } /* alpha == 0 */


   /*
    * if alpha is not equal to zero, we can compute the vectors
    * in the traditional manner.
    *
    */

   else {

      rtemp = 1.0 / sv.r;
      ftemp = rtemp;
      for (n = 1; n < p; n++) {
         fr += ftemp * (Y_C[n][0].x*Lin[n][0].x + Y_C[n][0].y*Lin[n][0].y);
         for (m = 1; m<=n; m++) {
            fr += (ftemp+ftemp) *
	       (Y_C[n][m].x*Lin[n][m].x + Y_C[n][m].y*Lin[n][m].y);
         } /* for m */
         ftemp += rtemp;
      } /* for n */


      for (n = 0; n < p; n++) {
         ftemp = 2.0;
         for (m = 1; m<=n; m++) {
            fb -= ftemp * (Y_C[n][m].y*Lin[n][m].x - Y_C[n][m].x*Lin[n][m].y);
            ftemp += 2.0;
         } /* for m */
      }/* for n */
      fb /= (sv.r * sinalpha);


      /*
       *  compute the alpha component of the gradient
       *  this is an ugly loop with lots of non-obvious
       *  temp vars, but it does use the least number of
       *  FP multiplies.
       */

      ntemp = 2.0;
      for (n = 1; n < p; n++) {

         fa -= ntemp *
	    ((Y_C[n][1].x*cosbeta + Y_C[n][1].y*sinbeta) * Lin[n][0].x +
             (Y_C[n][1].y*cosbeta - Y_C[n][1].x*sinbeta) * Lin[n][0].y );

         mtemp = 1.0;
         for (m = 1; m < n; m++) {

	    ftemp = 2.0 * mtemp * cotalpha;
	    ptemp = 2.0 * (ntemp+mtemp);

            fa +=  ftemp *
               (Y_C[n][m].x * Lin[n][m].x + Y_C[n][m].y * Lin[n][m].y);

            fa -= ptemp *
	       ((Y_C[n][m+1].x*cosbeta + Y_C[n][m+1].y*sinbeta) * Lin[n][m].x +
                (Y_C[n][m+1].y*cosbeta - Y_C[n][m+1].x*sinbeta) * Lin[n][m].y );

            mtemp += 1.0;
         } /* for m */

         fa += 2.0 * mtemp * cotalpha * 
            (Y_C[n][n].x*Lin[n][n].x + Y_C[n][n].y*Lin[n][n].y); 

         ntemp += 1.0;

      } /* for n */

      fa /= sv.r;

      /*
       *  rotate the spherical vector coordinates to cartesian
       */

      *potp = q * pot;
      fv->x = q * (fr*sinalpha*cosbeta + fa*cosalpha*cosbeta - fb*sinbeta);
      fv->y = q * (fr*sinalpha*sinbeta + fa*cosalpha*sinbeta + fb*cosbeta);
      fv->z = q * (fr*cosalpha - fa*sinalpha);

   } /* else */

} /* Force_C */


/****************************************************************
 *
 *  ForceM_C computes the Coulomb force and potential on the particle
 *  with charge q at position pv.[x,y,z], due to the Multipole
 *  expansion Min of size p.  The results are returned in potp and
 *  fv.[x,y,z]
 *
 */

int ForceM_C(
   Mtype   Min,
   int     p,
   Real    q,
   Vector  pv,
   Real    *potp,
   Vector  *fv )
{
   int             n, m;
   Real            pot, fr, fa, fb;
   Real            cosalpha, sinalpha, cosbeta, sinbeta;
   Real            cotalpha;
   Real            rtemp, ntemp, mtemp, gtemp, ptemp;
   SphVector       sv;

   /*
    * convert input to a spherical vector
    */

   Cart2Sph(pv,&sv);

   cosalpha = cos(sv.a);
   sinalpha = sin(sv.a);
   if (sinalpha < SMALL_THETA) {
      sinalpha = 0.0;
      if (cosalpha > 0.0)
	 cosalpha = 1.0;
      else
	 cosalpha = -1.0;
   }				/* if a small alpha */
   cosbeta = cos(sv.b);
   sinbeta = sin(sv.b);
   cotalpha = cosalpha/sinalpha;

   pot = 0.0;
   fr = 0.0;
   fa = 0.0;
   fb = 0.0;

   /*
    * compute the G matrix, leaving the result in Y_C{n][m]
    */

   makeG(p, sv);


   /*
    * find the potential by computing the product of M[n,m].G[n,m]
    */

   for (n = 0; n < p; n++) {
      pot += Y_C[n][0].x*Min[n][0].x - Y_C[n][0].y*Min[n][0].y;
      for (m = 1; m<=n; m++) {
         pot += 2.0 * (Y_C[n][m].x*Min[n][m].x - Y_C[n][m].y*Min[n][m].y);
      } /* for m */
   }/* for n */


   /*
    *  this section computes the special case for the
    *  condition where alpha=0 (assume beta=0)
    *  rho cannot equal zero because the multipole expansion does 
    *  not converge.
    */

   if ( sinalpha == 0.0 ) {

      rtemp = 1.0 / sv.r;
      gtemp = rtemp;
      for (n = 0; n < p; n++) {
         fr += gtemp * (Y_C[n][0].x*Min[n][0].x - Y_C[n][0].y*Min[n][0].y);
         for (m = 1; m<=n; m++) {
            fr += 2.0 * gtemp *
	       (Y_C[n][m].x*Min[n][m].x - Y_C[n][m].y*Min[n][m].y);
         } /* for m */
         gtemp += rtemp;
      } /* for n */

      /*
       *  compute the alpha and beta component of the gradient
       */

      rtemp = 1.0 / sv.r;
      ntemp = 0.0;
      for (n = 1; n < p; n++) {
         rtemp /= sv.r;
         ntemp += (double)n;
         gtemp  = 2.0 * Gcoeff[n][1] * rtemp * ntemp;
         fa += gtemp * Min[n][1].x;
         fb -= gtemp * Min[n][1].y;
      } /* for n */


      fa /= sv.r;
      fb /= sv.r;

      *potp = q * pot;
      fv->x = q * fa;
      fv->y = q * fb;
      fv->z = q * fr;

   } /* alpha == 0 */


   /*
    * if alpha is not equal to zero, we can compute the vectors
    * in the traditional manner.
    *
    */

   else {

      /* compute the rho component of the gradient */

      rtemp = 1.0 / sv.r;
      gtemp = rtemp;
      for (n = 0; n < p; n++) {
         fr += gtemp * (Y_C[n][0].x*Min[n][0].x - Y_C[n][0].y*Min[n][0].y);
         for (m = 1; m<=n; m++) {
            fr += 2.0 * gtemp *
	       (Y_C[n][m].x*Min[n][m].x - Y_C[n][m].y*Min[n][m].y);
         } /* for m */
         gtemp += rtemp;
      } /* for n */


      /* compute the beta component of the gradient */

      for (n = 1; n < p; n++) {
         gtemp = 2.0;
         for (m = 1; m<=n; m++) {
            fb += gtemp * (Y_C[n][m].y*Min[n][m].x + Y_C[n][m].x*Min[n][m].y);
            gtemp += 2.0;
         } /* for m */
      }/* for n */
      fb /= (sv.r * sinalpha);


      /* compute the alpha component of the gradient */

      ntemp = 1.0;
      for (n = 1; n < p; n++) {
         fa -= ntemp *
	    ((Y_C[n][1].x*cosbeta + Y_C[n][1].y*sinbeta) * Min[n][0].x -
             (Y_C[n][1].y*cosbeta - Y_C[n][1].x*sinbeta) * Min[n][0].y );

         mtemp = 1.0;
         for (m = 1; m < n; m++) {

	    gtemp = 2.0 * mtemp * cotalpha;
	    ptemp = 2.0 * (ntemp-mtemp);

            fa += gtemp * (Y_C[n][m].x*Min[n][m].x - Y_C[n][m].y*Min[n][m].y); 

            fa -= ptemp *
	       ((Y_C[n][m+1].x*cosbeta + Y_C[n][m+1].y*sinbeta) * Min[n][m].x -
                (Y_C[n][m+1].y*cosbeta - Y_C[n][m+1].x*sinbeta) * Min[n][m].y );

            mtemp += 1.0;
         } /* for m */

         fa += 2.0 * mtemp * cotalpha * 
            (Y_C[n][n].x*Min[n][n].x - Y_C[n][n].y*Min[n][n].y); 

      ntemp += 1.0;
      } /* for n */

      fa /= sv.r;
      fa *= -1.0;

      /*
       *  rotate the spherical vector coordinates to cartesian
       */

      *potp = q * pot;
      fv->x = q * (fr*sinalpha*cosbeta + fa*cosalpha*cosbeta - fb*sinbeta);
      fv->y = q * (fr*sinalpha*sinbeta + fa*cosalpha*sinbeta + fb*cosbeta);
      fv->z = q * (fr*cosalpha - fa*sinalpha);

   } /* else */

} /* ForceM_C */


/****************************************************************
 *
 *  makeG() - compute the G transfer matrix, using the equation
 *  from Elliott's dissertation.  the results are stored in the
 *  global array Y_C[n][m].
 *
 */

void makeG(
   int        p,
   SphVector  sv )
{
   int             n, m;
   Real            rinv, rtemp, ytemp;

   Legendre(LegPoly, p, cos(sv.a));
   Fourier_C(p, sv.b);

   rinv = 1.0 / sv.r;
   rtemp = rinv;
   for (n = 0; n < p; n++) {
      for (m = 0; m <= n; m++) {
	 ytemp = rtemp * Gcoeff[n][m] * LegPoly[n][m];
	 Y_C[n][m].x = ytemp * Yxy[m].x;
	 Y_C[n][m].y = ytemp * Yxy[m].y;
      } /* for m */
      rtemp *= rinv;
   } /* for n */

} /* makeG */


/****************************************************************
 *
 *  makeF() - compute the F transfer matrix, as specified in 
 *  Elliott's dissertation.  the results are stored in the global
 *  array Y_C[n][m].
 *
 */

void makeF(
   int        p,
   SphVector  sv )
{
   int     m, n;
   Real    rtemp, ytemp;

   Legendre(LegPoly, p, cos(sv.a));
   Fourier_C(p, sv.b);

   rtemp = 1.0;

   for (n = 0; n < p; n++) {
      for (m = 0; m <= n; m++) {
	 ytemp = rtemp * Fcoeff[n][m] * LegPoly[n][m];
	 Y_C[n][m].x = ytemp * Yxy[m].x;
	 Y_C[n][m].y = ytemp * Yxy[m].y;
      } /* for m */
      rtemp *= sv.r;
   } /* for n */
} /* makeF */


/****************************************************************
 *
 *  makeYforce() - compute the Y matrix used in computing the
 *  forces and potentials due to a local expansion with a cell.
 *  results are left in the global matrix Y_C[n][m].
 *
 *  this is an artifact from the old versions of the code which
 *  did not use the F/G formulation of the multipole computations.
 *
 */

void makeYforceC(
   int     p,
   Real    r,
   Real    a,
   Real    b )
{
   int     m, n;
   Real    rtemp, ytemp;

   Legendre(LegPoly, p, cos(a));
   Fourier_C(p, b);

   ytemp = Ycoeff[0][0] * LegPoly[0][0];
   Y_C[0][0].x = ytemp * Yxy[0].x;
   Y_C[0][0].y = ytemp * Yxy[0].y;

   rtemp = 1.0;
   for (n = 1; n < p; n++) {
      for (m = 0; m <= n; m++) {
	 ytemp = rtemp * Ycoeff[n][m] * LegPoly[n][m];
	 Y_C[n][m].x = ytemp * Yxy[m].x;
	 Y_C[n][m].y = ytemp * Yxy[m].y;
      } /* for m */
      rtemp *= r;
   } /* for n */
} /* makeYforceC */


/****************************************************************
 *
 *  Cinit() - allocation and initialization of global arrays and
 *    data structures used for the multipole calculations.
 *
 */ 

void Cinit( int p )
{
   int      i, l, m, n;
   Real     *scratch, *factorial;
   Complex  *c_scratch;

   scratch = (Real *) malloc(((p * (p + 1)) / 2) * sizeof(Real));
   LegPoly = (Real **) malloc(p * sizeof(Real *));
   for (i = 0; i < p; i++) {
      LegPoly[i] = scratch;
      scratch += i + 1;
   } /* for i */

   scratch = &LegPoly[0][0];
   for (n = 0; n < (p * (p + 1)) / 2; n++) {
      scratch[n] = 0.0;
   }

   c_scratch = (Complex *) malloc(((p * (p + 1)) / 2) * sizeof(Complex));
   Y_C = (Mtype ) malloc(p * sizeof(Complex *));
   for (n = 0; n < p; n++) {
      Y_C[n] = c_scratch;
      c_scratch += n + 1;
   } /* for n */

   scratch = &Y_C[0][0].x;
   for (n = 0; n < (p * (p + 1)); n++)
      scratch[n] = 0.0;

   c_scratch = (Complex *) malloc(((p * (p + 1)) / 2) * sizeof(Complex));
   L = (Mtype ) malloc(p * sizeof(Complex *));
   for (n = 0; n < p; n++) {
      L[n] = c_scratch;
      c_scratch += n + 1;
   } /* for n */

   scratch = &L[0][0].x;
   for (n = 0; n < (p * (p + 1)); n++) {
      scratch[n] = 0.0;
   }

   Yxy = (Complex *) malloc((p + 1) * sizeof(Complex));

   factorial = (Real *) malloc(2 * (p + 1) * sizeof(Real));
   factorial[0] = 1.0;
   for (n = 1; n < 2 * (p + 1); n++) {
      factorial[n] = (Real) n *factorial[n - 1];
   } /* for n */

   scratch = (Real *) malloc(((p * (p + 1)) / 2) * sizeof(Real));
   A_C = (Real **) malloc((p + 1) * sizeof(Real **));
   for (n = 0; n < p; n++) {
      A_C[n] = scratch;
      scratch += n + 1;
   } /* for n */

   for (n = 0; n < p; n++) {
      for (m = 0; m <= n; m++) {
	 A_C[n][m] = 1.0 / sqrt(factorial[n + m] * factorial[n - m]);
      } /* for m */
   } /* for n */

   scratch = (Real *) malloc(((p * (p + 1)) / 2) * sizeof(Real));
   Ycoeff = (Real **) malloc(p * sizeof(Real **));
   for (n = 0; n < p; n++) {
      Ycoeff[n] = scratch;
      scratch += n + 1;
   } /* for n */

   for (n = 0; n < p; n++) {
      for (m = 0; m <= n; m++) {
	 Ycoeff[n][m] = pow(-1.0, (Real) m) *
	    sqrt(factorial[n - m] / factorial[n + m]);
      } /* for m */
   } /* for n */

   scratch = (Real *) malloc(((p * (p + 1)) / 2) * sizeof(Real));
   Fcoeff = (Real **) malloc(p * sizeof(Real **));
   for (n = 0; n < p; n++) {
      Fcoeff[n] = scratch;
      scratch += n + 1;
   } /* for n */

   for (n = 0; n < p; n++) {
      for (m = 0; m <= n; m++) {
	 Fcoeff[n][m] = pow(-1.0, (Real) (n + m)) / factorial[n + m];
      } /* for m */
   } /* for n */

   scratch = (Real *) malloc(((p * (p + 1)) / 2) * sizeof(Real));
   Gcoeff = (Real **) malloc(p * sizeof(Real **));
   for (n = 0; n < p; n++) {
      Gcoeff[n] = scratch;
      scratch += n + 1;
   } /* for n */

   for (n = 0; n < p; n++) {
      for (m = 0; m <= n; m++) {
	 Gcoeff[n][m] = pow(-1.0, (Real) (n + m)) * factorial[n - m];
      } /* for m */
   } /* for n */

   free(factorial);

} /* Cinit */


/****************************************************************
 *
 *  additional initializations needed if you are using the 
 *  FFT and precomputation enhancements.
 *
 */

void CinitFS(int p, int b)
{
   CallocFrevS(&Hm2l, p, b);

} /* CinitFS */


/****************************************************************
 *
 *  additional initializations needed if you are using the 
 *  FFT enhancements.
 *
 */

void CinitF(int p, int b)
{
   CallocFrev(&Hm2l, p, b);

} /* CinitF */


/****************************************************************
 *
 *  Ccleanup() - frees up the space from global arrays and
 *    data structures allocated during Cinit().
 *
 */ 

void Ccleanup( int p )
{

   free(LegPoly[0]);
   free(LegPoly);
   
   free(Y_C[0]);
   free(Y_C);

   free(L[0]);
   free(L);

   free(Yxy);

   free(A_C[0]);
   free(A_C);

   free(Ycoeff[0]);
   free(Ycoeff);

   free(Fcoeff[0]);
   free(Fcoeff);

   free(Gcoeff[0]);
   free(Gcoeff);

} /* Ccleanup */


/****************************************************************
 *
 *  CcleanupFS() - additional cleanup needed if you are
 *  using the FFT with precomputation enhancements and made
 *  an earlier call to CinitFS();
 *
 */

void CcleanupFS(int p, int b)
{
   CfreeFrevS(Hm2l, p, b);

} /* CcleanupFS */


/****************************************************************
 *
 *  CcleanupF() - additional cleanup needed if you are
 *  using the FFT enhancements and made an earlier call
 *  to CinitF();
 *
 */

void CcleanupF(int p, int b)
{
   CfreeFrev(Hm2l, p, b);

} /* CcleanupF */


/****************************************************************
 *
 *  M2L_C_F converts a Coulomb multipole expansion M in to a local
 *   expansion L, shifting the original expansion along the 
 *   vector <v>.  Both expansions are of size p. Uses fourier
 *   transforms with blocking factor b.
 *
 */

int M2L_C_F(
   Mtype  M,
   Mtype  L,
   int    p,
   int    b,
   Vector v )

{
   int     n, m, i, j, blocklen, nblocks, pf;
   Real    rho, alpha, beta, atemp;
   Real   *scratch, *inblock, *hblock, *outblock;
   Real   *inptr, *outptr, *hptr;
   Real    tempr1, tempi1, tempr2, tempi2;
   Real    realtemp, imgtemp;
   Real   *startloop, *endloop;
   Real    negone;

   SphVector   sv;


   /*
    * convert input to a spherical vector
    */

   Cart2Sph(v,&sv);

   /*
    * compute the G matrix, leaving the result in Y_C{n][m]
    */

   makeG(p, sv);


   pf = 1 << ((int) (log((double) (2 * p - 1)) / log(2.0)));
   hblock = &Hm2l[b - 1][0].x;
   for (i = 0; i < 4 * pf * p; i++)
      hblock[i] = 0.0;

   for (n = 0; n < p; n++) {
      negone = 1.0;
      for (m = 0; m <= n; m++) {
	 Hm2l[n][m].x = Y_C[n][m].x * negone;
	 Hm2l[n][m].y = -Y_C[n][m].y * negone;
	 negone *= -1.0;
      } /* for m */
      row_fft(&Hm2l[n][0].x, p);
   } /* for n */
   col_fft(hblock, p, b);


   blocklen = 4 * pf * b;
   nblocks = p / b;
   inblock = &M[0][0].x;
   outblock = &L[b - 1][0].x;
   hblock = &Hm2l[b - 1][0].x;
   for (i = 0; i < nblocks; i++) {
      inptr = inblock;
      hptr = hblock;
      startloop = outblock;
      endloop = outblock;
      endloop += blocklen;
      for (j = i; j < nblocks; j++) {
	 for (outptr = startloop; outptr < endloop;) {
	    tempr1 = *hptr++;
	    tempi1 = *hptr++;
	    tempr2 = *inptr++;
	    tempi2 = *inptr++;
	    realtemp = tempr1 * tempr2 - tempi1 * tempi2;
	    imgtemp = tempr1 * tempi2 + tempi1 * tempr2;
	    *outptr++ += realtemp;
	    *outptr++ += imgtemp;
	 } /* for outptr */
      } /* for j */
      outblock += blocklen;
      hblock += blocklen;
   } /* for i */

   return TRUE;

} /* M2L_C_F */


/****************************************************************
 *
 *  M2L_C_Fshort() - shift a remote multipole into a local expansion,
 *  using the transfer matrix and FFT enhancements.
 *
 */

int M2L_C_Fshort(
   Mtype M,
   Mtype L,
   Mtype H,
   int   p,
   int   b )
{
   int     n, m, i, j;
   int     blocklen, nblocks, pf, fullsize, pblock, sblocklen, step;
   int     pbig;
   Real    scratch, *inblock, *hblock, *outblock;
   Real    *inptr, *outptr, *hptr;
   Real    tempr1, tempi1, tempr2, tempi2;
   Real    realtemp, imgtemp;
   Real    *startloop, *endloop;
   Real    negone;

   pf = 1 << ((int) (log((double) (2 * p - 1)) / log(2.0)));
   blocklen = 4 * pf * b;
   nblocks = p / b;
   inblock = &M[0][0].x;
   outblock = &L[b - 1][0].x;
   hblock = &H[b - 1][0].x;
   for (i = 0; i < nblocks; i++) {
      inptr = inblock;
      hptr = hblock;
      startloop = outblock;
      endloop = outblock;
      endloop += blocklen;
      for (j = i; j < nblocks; j++) {
	 pblock = b * (j + 1);
	 pf = 1 << ((int) (log((double) (2 * pblock - 1)) / log(2.0)));
	 sblocklen = 4 * pf * b;
	 step = blocklen / sblocklen - 1;
	 for (outptr = startloop; outptr < endloop;) {
	    tempr1 = *hptr++;
	    tempi1 = *hptr++;
	    tempr2 = *inptr++;
	    tempi2 = *inptr++;
	    realtemp = tempr1 * tempr2 - tempi1 * tempi2;
	    imgtemp = tempr1 * tempi2 + tempi1 * tempr2;
	    *outptr++ += realtemp;
	    *outptr++ += imgtemp;
	    inptr += (step << 1);
	    outptr += (step << 1);
	 } /* for outptr */
      } /* for j */
      outblock += blocklen;
      pblock = b * (i + 1);
      pf = 1 << ((int) (log((double) (2 * pblock - 1)) / log(2.0)));
      sblocklen = 4 * pf * b;
      hblock += sblocklen;
   } /* for i */

   return TRUE;

} /* M2L_CFshort */



/****************************************************************
 *
 *  Warp_M2L() -take the FFT of the multipole expansion in M1, leaving
 *  the results in M2.
 *
 */

void Warp_M2L(
   Mtype  M1,
   Mtype  M2,
   int    p,
   int    b )
{
   int             n, m;

   for (n = 0; n < p; n++) {
      for (m = 0; m <= n; m++) {
	 M2[n][m].x = M1[n][m].x;
	 M2[n][m].y = M1[n][m].y;
      } /* for m */
      row_fft(&M2[n][0].x, p);
   } /* for n */

   col_fft(&M2[0][0].x, p, b);

} /* Warp_M2L */


/****************************************************************
 *
 *  Unwarp_M2L() - perform the inverst FFT on the local expansion
 *  in M1 and accumulate the results in M2.
 *
 */

void Unwarp_M2L(
   Mtype  M1,
   Mtype  M2,
   int    p,
   int    b )
{
   int             i, j, n, m, pf;
   int             nblocks;
   Real            negone;
   Complex        *row;
   Real           *scratch;

   nblocks = p / b;

   pf = 1 << ((int) (log((double) (2 * p - 1)) / log(2.0)));
   scratch = &M1[b - 1][0].x;
   col_ifft(scratch, p, b);
   for (i = 0; i < 2 * p; i++) {
      row_ifft(scratch, p);
      scratch += 2 * pf;
   } /* for i */

   row = &M1[b - 1][0];
   scratch = &M1[b - 1][0].x;
   for (i = 0; i < nblocks; i++) {
      n = b * (i + 1) - 1;
      for (j = 0; j < 2 * b - 1; j++) {
	 if (n >= 0) {
	    negone = 1.0;
	    for (m = 0; m <= n; m++) {
	       M2[n][m].x += row[m].x * negone / (Real) (4 * pf * b);
	       M2[n][m].y += -row[m].y * negone / (Real) (4 * pf * b);
	       negone *= -1.0;
	    } /* for m */
	 } /* if */
	 n--;
	 row += pf;
      } /* for j */
      row += pf;
   } /* for i */

} /* unwarp_m2l */


/****************************************************************
 *
 *  copyG() - calculate the the G matrix and store the results in
 *  Yout.
 *
 */

void copyG(
   Mtype   Yout,
   int     p,
   Vector  v )
{
   int             n, m;
   SphVector       sv;


   /* convert input to a spherical vector */

   Cart2Sph(v,&sv);

   /* compute the G matrix, leaving the result in Y_C{n][m] */

   makeG(p, sv);

   for (n = 0; n < p; n++) {
      for (m = 0; m <= n; m++) {
	 Yout[n][m].x = Y_C[n][m].x;
	 Yout[n][m].y = Y_C[n][m].y;
      }	/* for m */
   } /* for n */

} /* copyG */


/****************************************************************
 *
 *  addG() - calculate the the G matrix and sum the results in
 *  Yout.
 *
 */

void addG(
   Mtype   Yout,
   int     p,
   Vector  v )
{
   int             n, m;
   SphVector       sv;


   /* convert input to a spherical vector */

   Cart2Sph(v,&sv);

   /* compute the G matrix, leaving the result in Y_C{n][m] */

   makeG(p, sv);

   for (n = 0; n < p; n++) {
      for (m = 0; m <= n; m++) {
	 Yout[n][m].x += Y_C[n][m].x;
	 Yout[n][m].y += Y_C[n][m].y;
      }	/* for m */
   } /* for n */

} /* addG */


/****************************************************************
 *
 *  copyF() - calculate the the F matrix and store the results in
 *  Yout.
 *
 */

void copyF(
   Mtype   Yout,
   int     p,
   Vector  v )
{
   int             n, m;
   SphVector       sv;


   /* convert input to a spherical vector */

   Cart2Sph(v,&sv);

   /* compute the F matrix, leaving the result in Y_C{n][m] */
   sv.r = -sv.r;
   makeF(p, sv);

   for (n = 0; n < p; n++) {
      for (m = 0; m <= n; m++) {
	 Yout[n][m].x = Y_C[n][m].x;
	 Yout[n][m].y = Y_C[n][m].y;
      }	/* for m */
   } /* for n */

} /* copyF */


/****************************************************************
 *
 *  addF() - calculate the the F matrix and sum the results in
 *  Yout.
 *
 */

void addF(
   Mtype   Yout,
   int     p,
   Vector  v )
{
   int             n, m;
   SphVector       sv;


   /* convert input to a spherical vector */

   Cart2Sph(v,&sv);

   /* compute the F matrix, leaving the result in Y_C{n][m] */
   sv.r = -sv.r;
   makeF(p, sv);

   for (n = 0; n < p; n++) {
      for (m = 0; m <= n; m++) {
	 Yout[n][m].x += Y_C[n][m].x;
	 Yout[n][m].y += Y_C[n][m].y;
      }	/* for m */
   } /* for n */

} /* addF */


/****************************************************************
 *
 * Csize() - returns size of multipole expansion in number of 
 *   complex data entries
 *
 */

int Csize( int p )
{
   return (p*(p+1))/2;
}


/****************************************************************
 *
 * CsizeF() - returns size of multipole expansion in number of 
 *   complex data entries for Fourier represented multipoles
 *
 */

int CsizeF( int p )
{
   int pf;
   pf = 1 << ((int) (log((double)(2*p-1))/log(2.0)));
   return (2*pf*p);
}


/****************************************************************
 *
 *  CMclear() - clear a multipole expansion
 *
 */

void CMclear(
   Mtype   M,
   int     p )
{
   int     n, m;

   for (n = 0; n < p; n++) {
      for (m = 0; m <= n; m++) {
	 M[n][m].x = 0.0;
	 M[n][m].y = 0.0;
      } /* for m */
   } /* for n */

} /* CMclear */


/****************************************************************
 *
 *  CMclearF() - clear an FFT multipole expansion.
 *
 */

void CMclearF(
   Mtype   M,
   int     p )
{
   int     i, pf;
   Real    *scratch;

   pf = 1 << ((int) (log((double) (2 * p - 1)) / log(2.0)));

   scratch = &M[0][0].x;
   for (i = 0; i < 4 * pf * p; i++) {
      scratch[i] = 0.0;
   } /* for i */

} /* CMclearF */


/****************************************************************
 *
 *  CMclearFrev() - clear an FFT multipole expansion that was
 *  allocated with CMallocFrev.
 *
 */

void CMclearFrev(
   Mtype   M,
   int     p,
   int     blklen )
{
   int     i, pf;
   Real    *scratch;

   pf = 1 << ((int) (log((double) (2 * p - 1)) / log(2.0)));

   scratch = &M[blklen - 1][0].x;
   for (i = 0; i < 4 * pf * p; i++)
      scratch[i] = 0.0;

} /* CzeroFrev */


/****************************************************************
 *
 * CMclearFshort() - clears the fft version of the multipole
 *   transfer matrix expansion.
 *
 */

void CMclearFshort(
   Mtype   M,
   int     p,
   int     blklen )
{
   int     i, j, nblocks, pblock, pf, fullsize;
   Real    *scratch;

   nblocks = p / blklen;
   fullsize = 0;
   for (i = 0; i < nblocks; i++) {
      pblock = blklen * (i + 1);
      pf = 1 << ((int) (log((double) (2 * pblock - 1)) / log(2.0)));
      fullsize += 4 * pf * blklen;
   } /* for i */

   scratch = &M[blklen - 1][0].x;
   for (i = 0; i < fullsize; i++)
      scratch[i] = 0.0;

} /* CMclearFshort */


/****************************************************************
*
*  CMsum(M1,M2,p) -
*     adds the multipole expansion M1 (of size p) to M2
*
*/

void CMsum(
   Mtype  M1,
   Mtype  M2,
   int    p )
{

   int    i,j;   /* loop counters */

   for ( i=0; i<p; i++) {
      for ( j=0; j<=i; j++ ) {
         M2[i][j].x += M1[i][j].x;
         M2[i][j].y += M1[i][j].y;
      } /* for j */
   } /* for i */
} /* CMsum */


/****************************************************************
*
*  CMsumF(M1,M2,p) -
*     adds the multipole expansion M1 (of size p) to M2
*     the multipole uses FFT representation.
*
*/

void CMsumF(
   Mtype  M1,
   Mtype  M2,
   int    p )
{

   int i;                /* loop counter */
   int pf;               /* length of mpe */
   Complex *scratch1;    /* scratch pointer */
   Complex *scratch2;    /* scratch pointer */

   /* find power of 2 gt or equal to p */
   pf = 1 << ((int) (log((double)(2*p-1))/log(2.0)));

   scratch1 = &(M1[0][0]);
   scratch2 = &(M2[0][0]);

   for ( i=0; i<2*pf*p; i++) {
      scratch2[i].x += scratch1[i].x;
      scratch2[i].y += scratch1[i].y;
   } /* for i */

} /* CMsumF */


/****************************************************************
 *
 *  Warp_Short() - performs an fft on the multipole transfer
 *    matrix expansion.
 *
 */

void Warp_Short(
   Mtype   M,
   int     p,
   int     b )
{
   int             pbig, nblocks, n, m, i, j, pblock, pf;
   Real            warp_val, *hblock;

   pbig = 1 << ((int) (log((double) (2 * p - 1)) / log(2.0)));

   nblocks = p / b;
   n = 0;
   for (i = 0; i < nblocks; i++) {
      pblock = b * (i + 1);
      pf = 1 << ((int) (log((double) (2 * pblock - 1)) / log(2.0)));
      for (j = 0; j < b; j++) {
	 warp_val = (Real) (pbig / pblock);
	 for (m = 0; m <= n; m++) {
	    M[n][m].x *= warp_val;
	    M[n][m].y *= -warp_val;
	    warp_val *= -1.0;
	 } /* for m */
	 row_fft(&M[n][0].x, pf);
	 n++;
      } /* for j */
   } /* for i */
   hblock = &M[b - 1][0].x;
   col_fftS(hblock, p, b);

} /* Warp_Short */



/****************************************************************
 *
 * Compute the Fourier series components of angle b out to p+1 (0->p)
 * terms.  Results are left in the global array Yxy.
 *
 */

void Fourier_C(
   int   p,
   Real  b )
{
   int   m;

   if (Yxy == NULL) {
      fprintf(stderr, "Fourier called with null pointer to array\n");
      exit(0);
   }

   for (m = 0; m <= p; m++) {
      Yxy[m].x = cos((Real) m * b);
      Yxy[m].y = sin((Real) m * b);
   } /* for m */

} /* Fourier_C */



/****************************************************************
 *
 * MCM_C_Orig() - convolves Coulomb multipole expansion M1 with
 *   another expansion M2, both size p, placing the results in M3.
 *   this is shown in the following equation:
 *
 *     M3[n',m'] = Sum{n,mf}( M1[n,m] * M2[n'-n,m'-m] )
 *
 *   when using this to compute multipole sum forces, M1 is the
 *   remote multipole and M2 is the local region.  ForceM_C should
 *   be called with M3 as an argument along with q=1.0 and the
 *   separation vector R pointing from M1 to M2.  this returns
 *   the sum of the effect of M1 on all particles in M2.
 *
 * This is the original version of the code which uses less efficient
 * explicit array indexing instead of the more efficient direct pointer
 * manipulation
 *
 */

int MCM_C_Orig(
   Mtype   M1,
   Mtype   M2,
   Mtype   M3,
   int     p)
{
   int             n, m, np, mp, startm, endm;
   Real            ntemp, mtemp;

   for (np = 0; np < p; np++) {
      for (mp = 0; mp <= np; mp++) {
	 for (n = 0; n <= np; n++) {

	    ntemp = (Real)(1 - 2 * (0x0001 & (np-n)));

	    startm = mp - (np - n);
	    startm = (startm < -n ? -n : startm);
	    endm = mp + (np - n);
	    endm = (endm > n ? n : endm);
	    if (startm <= endm) {
	       m = startm;

               /* m is negative, (mp-m) is positive. */
	       while (m < 0 && m <= endm) {
 	          mtemp = (Real)(1 - 2 * (0x0001 & (-m)));
		  M3[np][mp].x += ntemp * mtemp *
		     (M1[n][-m].x * M2[np - n][mp - m].x +
		      M1[n][-m].y * M2[np - n][mp - m].y);
		  M3[np][mp].y += ntemp * mtemp *
		     (M1[n][-m].x * M2[np - n][mp - m].y -
		      M1[n][-m].y * M2[np - n][mp - m].x);
		  m++;
	       } /* while m negative */

               /* m is positive, (mp-m) is positive. */
	       while (m < mp && m <= endm) {
		  M3[np][mp].x += ntemp *
		     (M1[n][m].x * M2[np - n][mp - m].x -
		      M1[n][m].y * M2[np - n][mp - m].y);
		  M3[np][mp].y += ntemp *
		     (M1[n][m].x * M2[np - n][mp - m].y +
		      M1[n][m].y * M2[np - n][mp - m].x);
		  m++;
	       } /* while m, positive, less than mp */

               /* m is positive, (mp-m) is negative. */
	       while (m <= endm) {
 	          mtemp = (Real)(1 - 2 * (0x0001 & (m-mp)));
		  M3[np][mp].x += ntemp * mtemp *
		     (M1[n][m].x * M2[np - n][m - mp].x +
		      M1[n][m].y * M2[np - n][m - mp].y);
		  M3[np][mp].y += ntemp * mtemp *
		     (M1[n][m].y * M2[np - n][m - mp].x -
		      M1[n][m].x * M2[np - n][m - mp].y);
		  m++;
	       } /* while m */

	    } /* if startm <= endm (ie m exists) */
	 } /* for n */
      } /* for mp */
   } /* for np */

   return TRUE;

} /* MCM_C */


/****************************************************************
 *
 * MCM_C convolves Coulomb multipole expansion M1 with another
 *   expansion M2, both size p, placing the results in M3.  this is
 *   shown in the following equation:
 *
 *     M3[n',m'] = Sum{n,mf}( M1[n,m] * M2[n'-n,m'-m] )
 *
 *   when using this to compute multipole sum forces, M1 is the
 *   remote multipole and M2 is the local region.  ForceM_C should
 *   be called with M3 as an argument along with q=1.0 and the
 *   separation vector R pointing from M1 to M2.  this returns
 *   the sum of the effect of M1 on all particles in M2.
 *
 */

int MCM_C(
   Mtype   M1,
   Mtype   M2,
   Mtype   M3,
   int     p)
{
   int     n, m, np, mp, startm, endm;
   Real    ntemp, mtemp;
   Complex *m1p, *m2p, *m3p;

   m3p = &(M3[0][0]);

   for (np = 0; np < p; np++) {
      for (mp = 0; mp <= np; mp++) {
	 for (n = 0; n <= np; n++) {

	    ntemp = (Real)(1 - 2 * (0x0001 & (np-n)));

	    startm = mp - (np - n);
	    startm = (startm < -n ? -n : startm);
	    endm = mp + (np - n);
	    endm = (endm > n ? n : endm);
	    if (startm <= endm) {
	       m = startm;
               m1p = &(M1[n][-m]);
               m2p = &(M2[np-n][mp-m]);

               /* m is negative, (mp-m) is positive. */
	       while (m < 0 && m <= endm) {
 	          mtemp = (Real)(1 - 2 * (0x0001 & (-m)));
		  m3p->x += ntemp * mtemp * (m1p->x * m2p->x + m1p->y * m2p->y);
		  m3p->y += ntemp * mtemp * (m1p->x * m2p->y - m1p->y * m2p->x);
		  m++;
		  m1p--;
		  m2p--;
	       } /* while m negative */

               m1p = &(M1[n][m]);
               /* m is positive, (mp-m) is positive. */
	       while (m < mp && m <= endm) {
		  m3p->x += ntemp * (m1p->x * m2p->x - m1p->y * m2p->y);
		  m3p->y += ntemp * (m1p->x * m2p->y + m1p->y * m2p->x);
		  m++;
		  m1p++;
		  m2p--;
	       } /* while m, positive, less than mp */

               m2p = &(M2[np-n][m-mp]);
               /* m is positive, (mp-m) is negative. */
	       while (m <= endm) {
 	          mtemp = (Real)(1 - 2 * (0x0001 & (m-mp)));
		  m3p->x += ntemp * mtemp * (m1p->x * m2p->x + m1p->y * m2p->y);
		  m3p->y += ntemp * mtemp * (m1p->y * m2p->x - m1p->x * m2p->y);
		  m++;
		  m1p++;
		  m2p++;
	       } /* while m */

	    } /* if startm <= endm (ie m exists) */
	 } /* for n */
	 m3p++;
      } /* for mp */
   } /* for np */

   return TRUE;

} /* MCM_C */

