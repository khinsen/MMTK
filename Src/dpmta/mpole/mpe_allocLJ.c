/*
 * mpe_allocLJ.c - allocates and frees LJ multipole expansion matrices
 *
 * w. t. rankin
 *
 * Copyright (c) 1997 Duke University
 * All Rights Reserved.
 *
 */

static char RcsId[] = "$Id: mpe_allocLJ.c,v 1.3 1997/11/03 18:46:41 wrankin Exp $";

/*
 * RCS History:
 *
 * $Log: mpe_allocLJ.c,v $
 * Revision 1.3  1997/11/03 18:46:41  wrankin
 * general cleanup/ansi-fication of code.  no new features.
 *
 * Revision 1.2  1997/05/13 17:52:14  wrankin
 * fixed irritation bug seen when math.h is not included
 * added a missing prototype to mpe.h
 *
 * Revision 1.1  1997/05/09  20:14:59  wrankin
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
 *  LJalloc() - allocate memory to hold an LJ multipole expansion
 *
 */

void LJalloc(
   MtypeLJ *Mptr,
   int     p )
{   
   int n, l, m;
   Real *scratch;
   Complex *scratchC;
   MtypeLJ M;
   
   scratchC = (Complex *) malloc(((p * (p+1) * (p+2))/6) * sizeof(Complex));
   M = (Complex ***) malloc( p * sizeof(Complex **));

   for (n = 0; n < p; n++) {
      M[n] = (Complex **) malloc((n+1) * sizeof(Complex *));
      for (l = 0; l <= n; l++) {
         M[n][l] = scratchC;
         scratchC += n-l+1;
      } /* for l */
   } /* for n */


   /*
    *  initialize array contents to zero
    */

   scratch = &M[0][0][0].x;
   for (n=0; n < (p * (p+1) * (p+2))/3; n++)
      scratch[n] = 0.0;

   /*
    * set return pointer to block
    */

   *Mptr = M;

} /* LJalloc */


/****************************************************************
 *
 *  LJfree() - free up memory from an LJ multipole.
 *
 */

void LJfree( 
   MtypeLJ M1,
   int     p )
{
   int i;

   free(M1[0][0]);
   for ( i=0; i<p; i++ ) {
      free(M1[i]);
   }
   free(M1);

} /* LJfree */
