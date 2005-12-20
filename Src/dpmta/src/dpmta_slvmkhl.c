/**********************************************************************
*
*  dpmta_slvmkhl.c - precompute the mpe transfer matrices and
*     direct interaction vectors
*
*  w.t.rankin
*
*  Copyright (c) 1995 Duke University
*  All rights reserved
*
*/

static char rcsid[] = "$Id: dpmta_slvmkhl.c,v 2.12 1997/11/07 16:49:39 wrankin Exp $";

/*
*  revision history:
* 
*  $Log: dpmta_slvmkhl.c,v $
*  Revision 2.12  1997/11/07 16:49:39  wrankin
*  massive cleanup of code.
*   - ansi-fication and inclusion of prototypes
*   - removed unused variables
*   - all (except the test) code compiles with minimal warnings under gcc.
*
*  Revision 2.11  1997/05/12 18:06:11  wrankin
*  added routines to clean up dynamic memory allocation when
*  PMTAinit() is called
*
 * Revision 2.10  1997/05/07  18:59:39  chumphre
 * major changes to implement ||-piped cells
 *
 * Revision 2.9  1997/05/06  17:05:48  wrankin
 * removed memory leak during x-fer matrix allocation.
 *
 * Revision 2.8  1997/04/08  19:11:20  wrankin
 * Fixed syntax problems in serial and LJ code.
 *
 * Revision 2.7  1997/03/26  20:36:27  wrankin
 * rewrote Ilist and IIlist code to handle reallocation of array sizes
 *   that happen as a resuly of re-scaling a non-cubic unit cell.
 *
 * Revision 2.6  1996/11/18  19:29:36  wrankin
 * update to virial code to use multipole virial calculation
 * old virial code is turned on by defining OLDVIRIAL
 * direct test program now computes absolute virial tensor.
 *
 * Revision 2.5  1996/08/20  17:12:53  wrankin
 * Enhancement to allow any number of processes to be utilized.
 * Formally the number of precesses was limited to a power of 2.
 *
 * Revision 2.4  1996/08/09  15:31:03  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.3  1996/02/29  21:13:48  wrankin
 * New relaease: 2.4 (.1)
 *    - simplified calling structure for initialization
 *    - macroscopic periodic code
 *    - fixed PBC calculation - all particles are now stored in the
 *      cell as positions relative to the cell center.
 *    - virial preasure tensor computed
 *    - fix to allow particles on the outer cube boundary to be
 *      included. (UNC fix)
 *    - fix to order reception of particle data during the distributed
 *      calling sequence (UIUC fix)
 *    - removed M2L code that didn't use transfer matrices
 *    - early hooks in to perform interaction list sorting.
 *    - fixed LJ scaling factor for 1/r^12 potential.
 *    - cleaned up the LJ interface.
 *    - and of course, my continued efforts to ANSI-fy this beast.
 *
 * Revision 2.2  1995/10/02  20:58:43  wrankin
 * changes to support creation and use of MPE x-fer matrix for
 *   LJ multipole calculation
 *
 * Revision 2.1  1995/07/10  02:47:18  wrankin
 * multipole processing code modified to use precomputed mpe transfer
 *   matrices.
 *
 * in addition, all multippole calculation routines have been removed
 *   from this source distribution and will be placed in a separate
 *   multipole library to be supplied externally.
 *
*
*/


/*
*  include files
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "dpmta_cell.h"
#include "dpmta_slvglobals.h"
#include "dpmta_slvmkil.h"


/*****************************************************************
*
*  Init_Hlist() - initialize the data structures that hold the
*    transfer matrices.
*
*/

void Init_Hlist()
{

   int i;

   /* allocate one Hlist for each cell position */
   Dpmta_Hlist = (Hlist *)malloc(8*sizeof(Hlist));
   if ( Dpmta_Hlist == NULL ) {
      fprintf(stderr, "ERROR: Init_Hlist() - malloc() failed\n");
      exit(-1);
   }

   for (i=0; i<8; i++) {
      Dpmta_Hlist[i].psize = 0;
      Dpmta_Hlist[i].ssize = 0;
      Dpmta_Hlist[i].dsize = 0;

      Dpmta_Hlist[i].plist = NULL;
      Dpmta_Hlist[i].slist = NULL;
#ifdef COMP_LJ
      Dpmta_Hlist[i].plist_lj = NULL;
      Dpmta_Hlist[i].slist_lj = NULL;
#endif
#ifdef VIRIAL
      Dpmta_Hlist[i].plist_vec = NULL;
      Dpmta_Hlist[i].slist_vec = NULL;
#endif
      Dpmta_Hlist[i].dlist_vec = NULL;

   } /* for i */

} /* Init_Hlist */


/****************************************************************
* 
*  Compute_Hlist() - construct the set of global mpe transfer
*    matrices that will be used during the M2L phase of the downward
*    mpe pass.
*
*    sice the transfer matrix is computed for a specific level,
*    it must be called once for each level traversed during the
*    downward pass.
*
*/

void Compute_Hlist( int level )
{

   int i,j;
   int sep;
#ifndef PIPED
   Vector clen;
#else
   Vector scaledv1;
   Vector scaledv2;
   Vector scaledv3;
#endif
   Vector hv;
   IntVector iv;

#ifndef PIPED
   /* determine cube edge length for this level */
   clen.x = (Dpmta_CellVector1.x / Dpmta_MaxCellLen) / (double)(1 << level);
   clen.y = (Dpmta_CellVector2.y / Dpmta_MaxCellLen) / (double)(1 << level);
   clen.z = (Dpmta_CellVector3.z / Dpmta_MaxCellLen) / (double)(1 << level);
#else
   scaledv1.x = (Dpmta_CellVector1.x / Dpmta_MaxCellLen) / (double)(1 << level);
   scaledv1.y = (Dpmta_CellVector1.y / Dpmta_MaxCellLen) / (double)(1 << level);
   scaledv1.z = (Dpmta_CellVector1.z / Dpmta_MaxCellLen) / (double)(1 << level);
   scaledv2.x = (Dpmta_CellVector2.x / Dpmta_MaxCellLen) / (double)(1 << level);
   scaledv2.y = (Dpmta_CellVector2.y / Dpmta_MaxCellLen) / (double)(1 << level);
   scaledv2.z = (Dpmta_CellVector2.z / Dpmta_MaxCellLen) / (double)(1 << level);
   scaledv3.x = (Dpmta_CellVector3.x / Dpmta_MaxCellLen) / (double)(1 << level);
   scaledv3.y = (Dpmta_CellVector3.y / Dpmta_MaxCellLen) / (double)(1 << level);
   scaledv3.z = (Dpmta_CellVector3.z / Dpmta_MaxCellLen) / (double)(1 << level);
#endif
   /* zero out the transfer matrices */
   if (Dpmta_FFT) {
      for (i=0; i < 8; i++) {
         for (j=0; j<Dpmta_Intlist[i].pcnt; j++) {
            CMclearFshort(Dpmta_Hlist[i].plist[j], Dpmta_Mp, Dpmta_FftBlock);
         } /* for j */
         for (j=0; j<Dpmta_Intlist[i].scnt; j++) {
            CMclearFshort(Dpmta_Hlist[i].slist[j], Dpmta_Mp, Dpmta_FftBlock);
         } /* for j */
      } /* for i */
   } /* if fft */



   /* compute the transfer matrix and place in array */
   for (i=0; i<8; i++) {
      for (j=0; j<Dpmta_Intlist[i].pcnt; j++) {

         /* compute int vector from parent to remote cell */
 	 /* Sep2Vec(Dpmta_Intlist[i].plist[j], &(iv)); */
         sep = Dpmta_Intlist[i].plist[j];
         SEP2VEC(sep,iv);
#ifndef PIPED
     	 /*  scale parent and add in vector from child to parent center */
         hv.x = ((-2.0 * (double)(iv.x)) - (0.5 - (double)(i&0x01))) * clen.x;
         hv.y = ((-2.0 * (double)(iv.y)) - (0.5 - (double)((i>>1)&0x01))) * clen.y;
         hv.z = ((-2.0 * (double)(iv.z)) - (0.5 - (double)((i>>2)&0x01))) * clen.z;
#else
         hv.x = (((-2.0*(double)(iv.x))-(0.5-(double)(i&0x01)))*scaledv1.x) +
           (((-2.0*(double)(iv.y))-(0.5-(double)((i>>1)&0x01)))*scaledv2.x) +
           (((-2.0*(double)(iv.z))-(0.5-(double)((i>>2)&0x01)))*scaledv3.x);
         hv.y = (((-2.0*(double)(iv.x))-(0.5-(double)(i&0x01)))*scaledv1.y) +
           (((-2.0*(double)(iv.y))-(0.5-(double)((i>>1)&0x01)))*scaledv2.y) +
           (((-2.0*(double)(iv.z))-(0.5-(double)((i>>2)&0x01)))*scaledv3.y);
         hv.z = (((-2.0*(double)(iv.x))-(0.5-(double)(i&0x01)))*scaledv1.z) +
           (((-2.0*(double)(iv.y))-(0.5-(double)((i>>1)&0x01)))*scaledv2.z) +
           (((-2.0*(double)(iv.z))-(0.5-(double)((i>>2)&0x01)))*scaledv3.z);
#endif

	 /* compute the transfer matrix */
         copyG(Dpmta_Hlist[i].plist[j], Dpmta_Mp, hv);

#ifdef COMP_LJ
         copyYI(Dpmta_Hlist[i].plist_lj[j], Dpmta_Mp_LJ, hv);
#endif

#ifdef VIRIAL
         /* store the vectors if we need them */
         Dpmta_Hlist[i].plist_vec[j].x = hv.x;
         Dpmta_Hlist[i].plist_vec[j].y = hv.y;
         Dpmta_Hlist[i].plist_vec[j].z = hv.z;
#endif

      } /* for j */

      for (j=0; j<Dpmta_Intlist[i].scnt; j++) {

         /* compute int vector from local to remote cell */
	 /* Sep2Vec(Dpmta_Intlist[i].slist[j], &(iv)); */
         sep = Dpmta_Intlist[i].slist[j];
         SEP2VEC(sep,iv);
	 
         /* compute vector from local to remote cellt */
#ifndef PIPED
	 hv.x = -1.0 * (double)(iv.x) * clen.x;
	 hv.y = -1.0 * (double)(iv.y) * clen.y;
	 hv.z = -1.0 * (double)(iv.z) * clen.z;
#else
         hv.x = (-1.0 * (double)(iv.x) * scaledv1.x) +
                (-1.0 * (double)(iv.y) * scaledv2.x) +
                (-1.0 * (double)(iv.z) * scaledv3.x);
         hv.y = (-1.0 * (double)(iv.x) * scaledv1.y) +
                (-1.0 * (double)(iv.y) * scaledv2.y) +
                (-1.0 * (double)(iv.z) * scaledv3.y);
         hv.z = (-1.0 * (double)(iv.x) * scaledv1.z) +
                (-1.0 * (double)(iv.y) * scaledv2.z) +
                (-1.0 * (double)(iv.z) * scaledv3.z);
#endif

	 /* compute the transfer matrix */
         copyG(Dpmta_Hlist[i].slist[j], Dpmta_Mp, hv);

#ifdef COMP_LJ
         copyYI(Dpmta_Hlist[i].slist_lj[j], Dpmta_Mp_LJ, hv);
#endif

#ifdef VIRIAL
         /* store the vectors if we need them */
         Dpmta_Hlist[i].slist_vec[j].x = hv.x;
         Dpmta_Hlist[i].slist_vec[j].y = hv.y;
         Dpmta_Hlist[i].slist_vec[j].z = hv.z;
#endif

      } /* for j */
   } /* for i */


   /*
   *  if we are doing FFTs prewarp the transfer matrix
   */
   if (Dpmta_FFT) {
      for (i=0; i < 8; i++) {

         for (j=0; j<Dpmta_Intlist[i].pcnt; j++) {
            Warp_Short(Dpmta_Hlist[i].plist[j], Dpmta_Mp, Dpmta_FftBlock);
         } /* for j */

         for (j=0; j<Dpmta_Intlist[i].scnt; j++) {
            Warp_Short(Dpmta_Hlist[i].slist[j], Dpmta_Mp, Dpmta_FftBlock);
         } /* for j */

      } /* for i */

   } /* if fft */

} /* Compute_Hlist() */


/****************************************************************
* 
*  Make_RelVec() - construct the set of global direct vectors.
*
*  since the direct interactions are only computed at a single level
*  in the tree, this routine only needs to be called once before
*  the beginning of processing, and once whenever the cell is resized.
*/

void Make_RelVec( int level )
{

   int i,j;
   int sep;
#ifndef PIPED
   Vector clen;
#else
   Vector scaledv1;
   Vector scaledv2;
   Vector scaledv3;
#endif
   IntVector iv;

#ifndef PIPED
   /* determine cube edge length for this level */
   clen.x = (Dpmta_CellVector1.x / Dpmta_MaxCellLen) / (double)(1 << level);
   clen.y = (Dpmta_CellVector2.y / Dpmta_MaxCellLen) / (double)(1 << level);
   clen.z = (Dpmta_CellVector3.z / Dpmta_MaxCellLen) / (double)(1 << level);
#else
   scaledv1.x = (Dpmta_CellVector1.x / Dpmta_MaxCellLen) / (double)(1 << level);
   scaledv1.y = (Dpmta_CellVector1.y / Dpmta_MaxCellLen) / (double)(1 << level);
   scaledv1.z = (Dpmta_CellVector1.z / Dpmta_MaxCellLen) / (double)(1 << level);
   scaledv2.x = (Dpmta_CellVector2.x / Dpmta_MaxCellLen) / (double)(1 << level);
   scaledv2.y = (Dpmta_CellVector2.y / Dpmta_MaxCellLen) / (double)(1 << level);
   scaledv2.z = (Dpmta_CellVector2.z / Dpmta_MaxCellLen) / (double)(1 << level);
   scaledv3.x = (Dpmta_CellVector3.x / Dpmta_MaxCellLen) / (double)(1 << level);
   scaledv3.y = (Dpmta_CellVector3.y / Dpmta_MaxCellLen) / (double)(1 << level);
   scaledv3.z = (Dpmta_CellVector3.z / Dpmta_MaxCellLen) / (double)(1 << level);
#endif
   for (i=0; i<8; i++) {
      for (j=0; j<Dpmta_Intlist[i].dcnt; j++) {

         /* compute int vector from local to remote cell */
         /* Sep2Vec(Dpmta_Intlist[i].dlist[j], &(iv)); */
	 sep = Dpmta_Intlist[i].dlist[j];
	 SEP2VEC(sep,iv);

	 /* compute vector from local to remote cellt */
#ifndef PIPED
	 Dpmta_Hlist[i].dlist_vec[j].x = (double)(iv.x) * clen.x;
         Dpmta_Hlist[i].dlist_vec[j].y = (double)(iv.y) * clen.y;
         Dpmta_Hlist[i].dlist_vec[j].z = (double)(iv.z) * clen.z;
#else
         Dpmta_Hlist[i].dlist_vec[j].x = (double)(iv.x) * scaledv1.x +
                                         (double)(iv.y) * scaledv2.x +
                                         (double)(iv.z) * scaledv3.x;
         Dpmta_Hlist[i].dlist_vec[j].y = (double)(iv.x) * scaledv1.y +
                                         (double)(iv.y) * scaledv2.y +
                                         (double)(iv.z) * scaledv3.y;
         Dpmta_Hlist[i].dlist_vec[j].z = (double)(iv.x) * scaledv1.z +
                                         (double)(iv.y) * scaledv2.z +
                                         (double)(iv.z) * scaledv3.z;
#endif

      } /* for j */

   } /* for i */

} /* Make_Relvec() */



/****************************************************************
*
*  Make_Hlist() - allocated strage for the mpe transfer matrix
*    and direct interaction vectors.  Increase storage if more
*    is needed.
* 
*/

void Make_Hlist()
{

   int     i,j;              /* lop counters */
   int     npar, nsib, ndir; /* size of hlist arrays */
   Mtype   *tmp_mpe;         /* temp mpe pointer */
   Vector  *tmp_vec;         /* temp vector pointer */

#ifdef COMP_LJ
   MtypeLJ *tmp_mpe_lj;      /* temp mpe pointer */
#endif
   /*
    * check the size of the h-lists against the alocated size.
    * and increase if needed. 
    */

   
   for (i=0; i<8; i++) {

      npar = Dpmta_Intlist[i].pcnt;
      nsib = Dpmta_Intlist[i].scnt;
      ndir = Dpmta_Intlist[i].dcnt;

      if ( npar > Dpmta_Hlist[i].psize ) {
	 tmp_mpe = Dpmta_Hlist[i].plist;
	 Dpmta_Hlist[i].plist = (Mtype *)realloc(tmp_mpe,npar*sizeof(Mtype));
	 if ( Dpmta_Hlist[i].plist == NULL ) {
	    fprintf(stderr,"ERROR: Make_Hlist() - malloc failed\n");
	    exit(-1);
	 }

	 if (Dpmta_FFT) {
	    for (j=Dpmta_Hlist[i].psize; j < npar; j++) {
	       CallocFrevS(&(Dpmta_Hlist[i].plist[j]), Dpmta_Mp,
			   Dpmta_FftBlock);
	    } /* for j */
	 } /* if fft */
	 else {
	    for (j=Dpmta_Hlist[i].psize; j < npar; j++) {
	       Calloc(&(Dpmta_Hlist[i].plist[j]), Dpmta_Mp);
	    } /* for j */
	 } /* else not fft */

#ifdef COMP_LJ
	 tmp_mpe_lj = Dpmta_Hlist[i].plist_lj;
	 Dpmta_Hlist[i].plist_lj =
	    (MtypeLJ *)realloc(tmp_mpe_lj,npar*sizeof(MtypeLJ));
	 if ( Dpmta_Hlist[i].plist_lj == NULL ) {
	    fprintf(stderr,"ERROR: Make_Hlist() - malloc failed\n");
	    exit(-1);
	 }
	 for (j=Dpmta_Hlist[i].psize; j < npar; j++) {
	    LJalloc(&(Dpmta_Hlist[i].plist_lj[j]), Dpmta_Mp_LJ);
	 } /* for j */
#endif

#ifdef VIRIAL
	 tmp_vec = Dpmta_Hlist[i].plist_vec;
	 Dpmta_Hlist[i].plist_vec =
	    (Vector *)realloc(tmp_vec,npar*sizeof(Vector));
	 if ( Dpmta_Hlist[i].plist_vec == NULL ) {
	    fprintf(stderr,"ERROR: Make_Hlist() - malloc failed\n");
	    exit(-1);
	 }
#endif

	 Dpmta_Hlist[i].psize = npar;

      } /* if npar */


      if ( nsib > Dpmta_Hlist[i].ssize ) {

	 tmp_mpe = Dpmta_Hlist[i].slist;
	 Dpmta_Hlist[i].slist = (Mtype *)realloc(tmp_mpe,nsib*sizeof(Mtype));
	 if ( Dpmta_Hlist[i].slist == NULL ) {
	    fprintf(stderr,"ERROR: Make_Hlist() - malloc failed\n");
	    exit(-1);
	 }

	 if (Dpmta_FFT) {
	    for (j=Dpmta_Hlist[i].ssize; j < nsib; j++) {
	       CallocFrevS(&(Dpmta_Hlist[i].slist[j]), Dpmta_Mp,
			   Dpmta_FftBlock);
	    } /* for j */
	 } /* if fft */
	 else {
	    for (j=Dpmta_Hlist[i].ssize; j < nsib; j++) {
	       Calloc(&(Dpmta_Hlist[i].slist[j]), Dpmta_Mp);
	    } /* for j */
	 } /* else not fft */

#ifdef COMP_LJ
	 tmp_mpe_lj = Dpmta_Hlist[i].slist_lj;
	 Dpmta_Hlist[i].slist_lj =
	    (MtypeLJ *)realloc(tmp_mpe_lj,nsib*sizeof(MtypeLJ));
	 if ( Dpmta_Hlist[i].slist_lj == NULL ) {
	    fprintf(stderr,"ERROR: Make_Hlist() - malloc failed\n");
	    exit(-1);
	 }
	 for (j=Dpmta_Hlist[i].ssize; j < nsib; j++) {
	    LJalloc(&(Dpmta_Hlist[i].slist_lj[j]), Dpmta_Mp_LJ);
	 } /* for j */
#endif

#ifdef VIRIAL
	 tmp_vec = Dpmta_Hlist[i].slist_vec;
	 Dpmta_Hlist[i].slist_vec =
	    (Vector *)realloc(tmp_vec,nsib*sizeof(Vector));
	 if ( Dpmta_Hlist[i].slist_vec == NULL ) {
	    fprintf(stderr,"ERROR: Make_Hlist() - malloc failed\n");
	    exit(-1);
	 }
#endif

	 Dpmta_Hlist[i].ssize = nsib;

      } /* if nsib */


      if ( ndir > Dpmta_Hlist[i].dsize ) {
	 tmp_vec = Dpmta_Hlist[i].dlist_vec;
	 Dpmta_Hlist[i].dlist_vec =
	    (Vector *)realloc(tmp_vec,ndir*sizeof(Vector));
	 if ( Dpmta_Hlist[i].dlist_vec == NULL ) {
	    fprintf(stderr,"ERROR: Make_Hlist() - malloc failed\n");
	    exit(-1);
	 }

	 Dpmta_Hlist[i].dsize = ndir;

      } /* if ndir */

   } /* for i */

} /* Make_Hlist() */


/****************************************************************
 *
 *  Delete_Hlist() -
 *
 *  for completeness, we have added a routine that will free
 *  all allocated data structures that are malloc'd within the
 *  scope of this module.
 *
 *  this will effectively destroy the interaction xfer matrix lists,
 *  which could be possibly reinitialized by another call to
 *  Init_Hlist().
 *
 */

void Delete_Hlist()
{
   int i,j;

   /*
    *  de-alloc the interaction lists
    */

   for ( i=0; i<8; i++ ) {
#ifndef NOPARCONV
      for ( j=0; j<Dpmta_Hlist[i].psize; j++ ) {
         if (Dpmta_FFT) {
	    CfreeFrevS(Dpmta_Hlist[i].plist[j], Dpmta_Mp, Dpmta_FftBlock);
	 }
	 else {
	    Cfree(Dpmta_Hlist[i].plist[j], Dpmta_Mp);
	 }
      } /* for j */

      free( Dpmta_Hlist[i].plist );

#ifdef COMP_LJ
      for ( j=0; j<Dpmta_Hlist[i].psize; j++ ) {
	 LJfree(Dpmta_Hlist[i].plist_lj[j], Dpmta_Mp_LJ);
      }
      free( Dpmta_Hlist[i].plist_lj );
#endif
#ifdef VIRIAL
      free( Dpmta_Hlist[i].plist_vec );
#endif
#endif
      
      for ( j=0; j<Dpmta_Hlist[i].ssize; j++ ) {
         if (Dpmta_FFT) {
	    CfreeFrevS(Dpmta_Hlist[i].slist[j], Dpmta_Mp, Dpmta_FftBlock);
	 }
	 else {
	    Cfree(Dpmta_Hlist[i].slist[j], Dpmta_Mp);
	 }
      }
      free( Dpmta_Hlist[i].slist );
#ifdef COMP_LJ
      for ( j=0; j<Dpmta_Hlist[i].ssize; j++ ) {
	 LJfree(Dpmta_Hlist[i].slist_lj[j], Dpmta_Mp_LJ);
      }
      free( Dpmta_Hlist[i].slist_lj );
#endif
#ifdef VIRIAL
      free( Dpmta_Hlist[i].slist_vec );
#endif
      free( Dpmta_Hlist[i].dlist_vec );

   } /* for i */
   
   free(Dpmta_Intlist);


} /* Delete_Hlist */

