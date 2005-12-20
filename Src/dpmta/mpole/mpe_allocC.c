/*
 *  mpe_allocC.c - allocation/deallocation routines for coulomb multipole
 *    expansions.
 *
 *  w. t. rankin
 *  w. elliott
 *
 *  Copyright (c) 1997 Duke University
 *  All Rights Reserved
 */

static char RCSid[] = "$Id: mpe_allocC.c,v 1.3 1997/11/03 18:46:40 wrankin Exp $";

/*
 * RSC History:
 *
 * $Log: mpe_allocC.c,v $
 * Revision 1.3  1997/11/03 18:46:40  wrankin
 * general cleanup/ansi-fication of code.  no new features.
 *
 * Revision 1.2  1997/05/13 17:52:13  wrankin
 * fixed irritation bug seen when math.h is not included
 * added a missing prototype to mpe.h
 *
 * Revision 1.1  1997/05/09  20:14:57  wrankin
 * added routines to de-allocate global arrays created in [C,LJ]init()
 * added routines to free up multipole expansion matrices
 * added LJ prototypes and ansi-fied more procedures
 *
 *
 */

/*
 * include files
 */

#include <stdlib.h>
#include <math.h>
#include "mpe.h"


/****************************************************************
 *
 * Calloc() - allocates a standard multipole expansion array
 *
 */

void Calloc(
   Mtype *Mptr,
   int   p )
{
   int      n, m;
   Real     *scratch;
   Complex  *c_scratch;
   Mtype    M;


   c_scratch = (Complex *) malloc(((p * (p + 1)) / 2) * sizeof(Complex));
   M = (Mtype ) malloc(p * sizeof(Mtype ));
   for (n = 0; n < p; n++) {
      M[n] = c_scratch;
      c_scratch += n + 1;
   } /* for n */
   scratch = &M[0][0].x;
   for (n = 0; n < (p * (p + 1)); n++) {
      scratch[n] = 0.0;
   } /* for n */
   *Mptr = M;

} /* Calloc */


/****************************************************************
 *
 *  CallocF() - allocates a multipole expansion array for use with
 *    FFT enhanced processing.
 *
 */

void CallocF(
   Mtype *Mptr,
   int   p,
   int   b )
{
   int      n, i, j, nblocks, pf;
   Real     *scratch;
   Complex  *c_scratch;
   Mtype    M;

   pf = 1 << ((int) (log((double) (2 * p - 1)) / log(2.0)));
   c_scratch = (Complex *) malloc(2 * pf * p * sizeof(Complex));
   M = (Mtype ) malloc(p * sizeof(Mtype ));
   nblocks = p / b;
   n = 0;
   for (i = 0; i < nblocks; i++) {
      for (j = 0; j < b; j++) {
	 M[n++] = c_scratch;
	 c_scratch += pf;
      } /* for j */
      c_scratch += b * pf;
   } /* for i */
   scratch = &M[0][0].x;
   for (n = 0; n < 4 * pf * p; n++)
      scratch[n] = 0.0;
   *Mptr = M;

} /* CallocF */


/****************************************************************
 *
 *  CallocFrev - allocates a multipole expansion array used for
 *    storing inverse FFT multipole data.
 *
 */

void CallocFrev(
   Mtype *Mptr,
   int   p,
   int   b )
{
   int      n, i, j, nblocks, pf;
   Real     *scratch;
   Complex  *c_scratch;
   Mtype    M;

   nblocks = p / b;

   pf = 1 << ((int) (log((double) (2 * p - 1)) / log(2.0)));
   c_scratch = (Complex *) malloc(2 * pf * p * sizeof(Complex));
   M = (Mtype ) malloc(p * sizeof(Mtype ));
   for (i = 0; i < nblocks; i++) {
      n = (i + 1) * b - 1;
      for (j = 0; j < b; j++) {
	 M[n--] = c_scratch;
	 c_scratch += pf;
      } /* for j */
      c_scratch += b * pf;
   } /* for i */
   scratch = &M[b - 1][0].x;
   for (n = 0; n < 4 * pf * p; n++)
      scratch[n] = 0.0;
   *Mptr = M;

} /* CallocFrev */


/****************************************************************
 *
 *  CallocFrevS - allocates a multipole expansion array used for
 *    storing inverse FFT multipole data.
 *
 *    for use with precomputed multipole transfer matrix
 *    enhancements.
 *
 */

void CallocFrevS(
   Mtype *Mptr,
   int   p,
   int   b )
{
   int      n, i, j, nblocks, pf, fullsize, pblock;
   Real     *scratch;
   Complex  *c_scratch;
   Mtype    M;

   nblocks = p / b;

   fullsize = 0;
   for (i = 0; i < nblocks; i++) {
      pblock = b * (i + 1);
      pf = 1 << ((int) (log((double) (2 * pblock - 1)) / log(2.0)));
      fullsize += 4 * pf * b;
   } /* for i */
   pf = 1 << ((int) (log((double) (2 * p - 1)) / log(2.0)));
   c_scratch = (Complex *) malloc((fullsize / 2) * sizeof(Complex));
   M = (Mtype ) malloc(p * sizeof(Mtype ));
   for (i = 0; i < nblocks; i++) {
      pblock = b * (i + 1);
      pf = 1 << ((int) (log((double) (2 * pblock - 1)) / log(2.0)));
      n = (i + 1) * b - 1;
      for (j = 0; j < b; j++) {
	 M[n--] = c_scratch;
	 c_scratch += pf;
      } /* for j */
      c_scratch += b * pf;
   } /* for i */
   scratch = &M[b - 1][0].x;
   for (n = 0; n < fullsize; n++)
      scratch[n] = 0.0;
   *Mptr = M;

} /* CallocFrevS */



/****************************************************************
 *
 * Cfree() - frees up a standard multipole expansion array
 *
 */

void Cfree(
   Mtype M1,
   int p )
{

   free(M1[0]);
   free(M1);
   
} /* Cfree */


/****************************************************************
 *
 *  CfreeF() - allocates a multipole expansion array for use with
 *    FFT enhanced processing.
 *
 */

void CfreeF(
   Mtype M1,
   int   p,
   int   b )
{

   free(M1[0]);
   free(M1);

} /* CfreeF */


/****************************************************************
 *
 *  CfreeFrev - allocates a multipole expansion array used for
 *    storing inverse FFT multipole data.
 *
 */

void CfreeFrev(
   Mtype M1,
   int   p,
   int   b )
{

   free(M1[b-1]);
   free(M1);
      
} /* CfreeFrev */


/****************************************************************
 *
 *  CfreeFrevS - frees up a multipole expansion array used for
 *    storing inverse FFT multipole data.
 *
 */

void CfreeFrevS(
   Mtype M1,
   int   p,
   int   b )
{

   free(M1[b-1]);
   free(M1);

} /* CfreeFrevS */


