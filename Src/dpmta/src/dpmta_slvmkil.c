/**********************************************************************
*
*  dpmta_slvmkil.c - compute the relative interaction lists.
*
*  w.t.rankin
*
*  Copyright (c) 1995 Duke University
*  All rights reserved
*
*/

static char rcsid[] = "$Id: dpmta_slvmkil.c,v 2.13 1998/04/01 20:08:26 wrankin Exp $";

/*
*  revision history:
* 
*  $Log: dpmta_slvmkil.c,v $
*  Revision 2.13  1998/04/01 20:08:26  wrankin
*  added support for HILBERT ordering
*  general cleanup of code structure (more OOP if you will)
*
*  Revision 2.12  1997/11/07 16:49:46  wrankin
*  massive cleanup of code.
*   - ansi-fication and inclusion of prototypes
*   - removed unused variables
*   - all (except the test) code compiles with minimal warnings under gcc.
*
*  Revision 2.11  1997/05/12 18:06:14  wrankin
*  added routines to clean up dynamic memory allocation when
*  PMTAinit() is called
*
 * Revision 2.10  1997/05/07  18:59:41  chumphre
 * major changes to implement ||-piped cells
 *
 * Revision 2.9  1997/03/26  20:36:31  wrankin
 * rewrote Ilist and IIlist code to handle reallocation of array sizes
 *   that happen as a resuly of re-scaling a non-cubic unit cell.
 *
 * Revision 2.8  1997/02/07  18:13:04  wrankin
 * fixed code to handle -DNOPARCONV flag correctly
 *
 * Revision 2.7  1996/09/24  18:43:18  wrankin
 * many changes for support of version 2.5
 *  - non cubic cells
 *  - non power-of-2 processors
 *  - fixed the resize code
 *  - new virial interface
 *  - changes to macroscopic code (still not working)
 *  - code cleanup
 *  - new test program that reads PDB file
 *  - update code for T3D support
 *
 * Revision 2.6  1996/08/20  17:12:55  wrankin
 * Enhancement to allow any number of processes to be utilized.
 * Formally the number of precesses was limited to a power of 2.
 *
 * Revision 2.5  1996/08/09  15:31:05  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.4  1996/02/29  21:13:51  wrankin
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
 * Revision 2.3  1995/11/29  22:29:30  wrankin
 * addition of periodic boundary code
 * beginning work on virial code.
 *
 * Revision 2.2  1995/07/01  03:27:02  wrankin
 * initial cut at precomputation of mpe transfer functions.
 * this works for vanilla M2L.  it does not work for FFT enhancements.
 *
 * Revision 2.1  1995/06/27  14:20:32  wrankin
 * Implementation of relative interaction list algorithm.
 * Interaction list generation code removed from master library.
 * Cell Table restructured as a 2-D table indexed by level/cell
 * All direct computations are now performed in a single pass.
 * General cleanup and removal of old files.
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


/*
*  prototyping
*/

#include "dpmta_slvmkil.h"
#include "dpmta_slvscale.h"

/* internal prototyping */

void Sort_Ilist( IntVector *, int *, int );

/*
 *  some local globals used in ilist creation
 */

#ifndef NOPARCONV
static IntVector *Tmp_Plist;
#endif
static IntVector *Tmp_Slist;
static IntVector *Tmp_Dlist;
static int Tmp_Size;

#ifdef SORTILIST
static int *Tmp_Sort;
#endif


/****************************************************************
*
*  MAC - multipole acceptance criteria
*
*  this function implementst the multipole acceptance criteria
*  (MAC).  it returns a TRUE(1) or FALSE(0) value depending if the
*  cells are well-separated or not.
*
*  actually, this should probably be put in as a macro, but most
*  decent optimizers will inline the function.
*
*/

int MAC(
   double ra,                   /* max radius of cell of interest */
   double rb,                   /* max radius of remote cell */
   double R,                    /* distance between cells */
   double theta)                /* separation criteria */
{

   if ( (ra+rb) <= (R*theta) ) {
      return(TRUE);
   }
   else {
     return(FALSE);
   }

}


/****************************************************************
*
*  Cell2Cell - translation routines
*
*  this function takes a cell id and a separation index and returns
*  the new cell id pointed to.  the value the function returns is
*  TRUE if the cell is within the boundaries of the simulation space,
*  TRUE if periodic boundary conditions are in effect,
*  and FALSE otherwise.
*
*  the algorithm is based upon section 3.2 of elliotts dissertation
*  with some enhancements.  an important note is that the level mask
*  is stored as a negated value.
*/

int Cell2Cell(
   int level,      /* cell level */
   int cell,       /* local cell id */
   int sep,        /* cell separation index */
   int *outcell,   /* remote cell id */
   int *vflag)     /* overflow flag */
{

   int nlvlmask;   /* a negated level mask - set bits high */
   int tmpout;     /* temp storage for output cell id */

   nlvlmask = -1 << (3*level);
   tmpout =  ((cell | ~ILMASK1) + (ILMASK1 & sep)) & ILMASK1;
   tmpout |= ((cell | ~ILMASK2) + (ILMASK2 & sep)) & ILMASK2;
   tmpout |= ((cell | ~ILMASK3) + (ILMASK3 & sep)) & ILMASK3;
   *vflag = tmpout & nlvlmask;
   *outcell = tmpout & ~nlvlmask;

   if (*vflag == 0)
      return(TRUE);
   else if (Dpmta_PBC)
      return(TRUE);
   else
      return(FALSE);

} /* Cell2Cell */


/****************************************************************
*
*  Vec2Sep - translation routine
*
*  this function takes a  vector offset between two cells and 
*  returns a cell separation index
*
*  the algorithm is based upon section 3.2 of elliotts disertation
*  with some enhancements.
*
*/

int Vec2Sep(
   IntVector offset,    /* remote cell id */
   int *sep)            /* cell separation index */
{

   int i;
   int mask;

   /* pack info from three vectors into a single index */
   *sep = 0;
   offset.y <<= 1;
   offset.z <<= 2;
   mask = 0x1;
   for (i=0; i<LEVELS_MAX; i++) {
      *sep |= offset.x & mask;
      offset.x <<= 2;
      mask <<= 1;
      *sep |= offset.y & mask;
      offset.y <<= 2;
      mask <<= 1;
      *sep |= offset.z & mask;
      offset.z <<= 2;
      mask <<= 1;
   } /* for m */

   return(TRUE);

} /* Vec2Sep */


/****************************************************************
*
*  Sep2Vec - translation routine
*
*  this function takes a cell separation index and returns
*  the offsett to the new cell
*
*  the algorithm is based upon section 3.2 of elliotts disertation
*  with some enhancements.  an important note is that the level mask
*  is stored as a negated value.
*/

int Sep2Vec(
   int sep,            /* cell separation index */
   IntVector *offset)  /* remote cell id */
{

   int i;
   int mask, topmask;

   offset->x = 0;
   offset->y = 0;
   offset->z = 0;

   mask = 0x01;
   for ( i=0; i<LEVELS_MAX; i++ ) {
      offset->x |= mask & sep;
      sep >>= 1;
      offset->y |= mask & sep;
      sep >>= 1;
      offset->z |= mask & sep;
      mask <<= 1;
   }

   topmask = (-1) << LEVELS_MAX;
   mask = 0x01 << (LEVELS_MAX-1);
   if ( offset->x & mask )
      offset->x |= topmask;
   if ( offset->y & mask )
      offset->y |= topmask;
   if ( offset->z & mask )
      offset->z |= topmask;

   return(TRUE);

} /* Sep2Vec */


/****************************************************************
*
*  Make_Ilist() - create the interaction lists
*
*  note that there are two different versions of this code.
*  one version is for use with the old fashion orthogonal cell
*  edges.  the second version supports parallel piped shapes.
*
*  please note that throughout the code, cell indexing is performed
*  in such a way that changes to the lowest bit of the cell index,
*  (bit[0]) represent changes in the direction of the X-axis.  
*  Changes to the second order bit (bit[1]) are changes in the y 
*  direction, and finally, changes to the third-order bit (bit[2])
*  represents changes in the z direction.
*
*  it is important that the code that partitions the particles
*  into the cells use this same convention.
*
*/
void Make_Ilist()
{

#ifndef PIPED

   int i,j,k,l,m,n;          /* loop counters */
   int dsize;                /* size of ilist arrays to allocate */
   int index;                /* temp index into ilist tables */
   int xadj, yadj, zadj;     /* terms used to adjust final vectors */
   void *tmp_ptr;            /* temp pointer */
   double rtmp;              /* temp radius value */
   double xtmp, ytmp, ztmp;  /* temp vector */
   IntVector dist;           /* length of temp cell list array */
   IntVector itmp;           /* temp vector used in calculations */

   int tmp_pcnt, tmp_scnt, tmp_dcnt;  /* length of above lists */

   Vector scaledlen;
   double radius1, radius1_2;

   /*
   *  establish how far from the parent cell we have to go in order
   *  to exceed the MAC.  this will form the limits on our search space
   *  for the child cells.
   */

   scaledlen.x = Dpmta_CellVector1.x / Dpmta_MaxCellLen;
   scaledlen.y = Dpmta_CellVector2.y / Dpmta_MaxCellLen;
   scaledlen.z = Dpmta_CellVector3.z / Dpmta_MaxCellLen;

   radius1 = sqrt( scaledlen.x*scaledlen.x + 
                   scaledlen.y*scaledlen.y + 
                   scaledlen.z*scaledlen.z );
   radius1_2 = radius1 / 2.0;

   dist.x = 1;
   dist.y = 1;
   dist.z = 1;

   while ( MAC(radius1, radius1, scaledlen.x * (double)(dist.x*2), Dpmta_Theta) == FALSE )
      dist.x++;
   while ( MAC(radius1, radius1, scaledlen.y * (double)(dist.y*2), Dpmta_Theta) == FALSE )
      dist.y++;
   while ( MAC(radius1, radius1, scaledlen.z * (double)(dist.z*2), Dpmta_Theta) == FALSE )
      dist.z++;

   /*
   *  now allocate and initialize the temporary ilist, one entry
   *  for every possible cell in the search space.
   *  we allocate entirely too much memory (more than we could 
   *  possibly use) but we'll free it in the end.
   */

   dsize = (dist.x+1)*(dist.y+1)*(dist.z+1);

   if ( dsize > Tmp_Size ) {

#ifndef NOPARCONV
      Tmp_Plist = (IntVector *)realloc(Tmp_Plist,dsize*sizeof(IntVector));
      if ( Tmp_Plist == NULL ) {
	 fprintf(stderr,"ERROR: malloc() failed\n");
	 exit(-1);
      }
#endif

      Tmp_Slist = (IntVector *)realloc(Tmp_Slist,8*dsize*sizeof(IntVector));
      if ( Tmp_Slist == NULL ) {
	 fprintf(stderr,"ERROR: malloc() failed\n");
	 exit(-1);
      }

      Tmp_Dlist = (IntVector *)realloc(Tmp_Dlist,8*dsize*sizeof(IntVector));
      if ( Tmp_Dlist == NULL ) {
	 fprintf(stderr,"ERROR: malloc() failed\n");
	 exit(-1);
      }

#ifdef SORTILIST
      Tmp_Sort = (int *)realloc(Tmp_Sort,8*dsize*sizeof(int));
      if ( Tmp_Sort == NULL ) {
	 fprintf(stderr,"ERROR: malloc() failed\n");
	 exit(-1);
      }
#endif

      Tmp_Size = dsize;

   } /* if dsize */

   tmp_pcnt = 0;
   tmp_scnt = 0;
   tmp_dcnt = 0;


   /*
   *  now we will build the temp interaction list for the child cell
   *  at position <0,0,0> within the center cell.  this corresponds to
   *  a real position of <-0.5, -0.5, -0.5> wrt the center cell.
   *
   *  we can exploit symetry here and only build the interaction list
   *  for a single child cell.  the remaining seven children can be
   *  constructed later by simply reversing the signs of one or more
   *  of the <x,y,z> components of the relative interaction vectors.
   */

   /* cycle through all the parent cells */

   for (i=-dist.x; i<=dist.x; i++) {
      for (j=-dist.y; j<=dist.y; j++) {
         for (k=-dist.z; k<=dist.z; k++) {

	    /*
	    *  if the cell meets the parents MAC, then we do
	    *   not have to consider it.
            */

            xtmp = 2.0 * (double)(i) * scaledlen.x;
            ytmp = 2.0 * (double)(j) * scaledlen.y;
            ztmp = 2.0 * (double)(k) * scaledlen.z;
            rtmp = sqrt( xtmp*xtmp + ytmp*ytmp + ztmp*ztmp );

            if ( MAC(radius1, radius1, rtmp, Dpmta_Theta) == FALSE) {

#ifndef NOPARCONV
               /*
               *  okay, the remote cell does not meet the MAC of the
	       *  parent.
               *
	       *  check to see if the child cell will interact
	       *  with the entire remote cell
	       */

               xtmp = (2.0 * (double)(i) * scaledlen.x) + (0.5 * scaledlen.x);
               ytmp = (2.0 * (double)(j) * scaledlen.y) + (0.5 * scaledlen.y);
               ztmp = (2.0 * (double)(k) * scaledlen.z) + (0.5 * scaledlen.z);
               rtmp = sqrt( xtmp*xtmp + ytmp*ytmp + ztmp*ztmp );

               if ( MAC(radius1_2, radius1, rtmp, Dpmta_Theta) == TRUE ) {
                  index = tmp_pcnt;
                  Tmp_Plist[index].x = i;
                  Tmp_Plist[index].y = j;
                  Tmp_Plist[index].z = k;
                  tmp_pcnt++;
	       } /* if MAC [2] */

               /*
	       *  else the child cannot interact with the remote
	       *  parent cell, so we will have to test each of the
	       *  children of the remote cell.
	       */

	       else {
#endif
                  for (l=0; l<2; l++) {
                     for (m=0; m<2; m++) {
                        for (n=0; n<2; n++) {

                           xtmp = (2.0 * (double)(i) * scaledlen.x) + 
                                  ((double)(l) * scaledlen.x);
                           ytmp = (2.0 * (double)(j) * scaledlen.y) + 
                                  ((double)(m) * scaledlen.y);
                           ztmp = (2.0 * (double)(k) * scaledlen.z) + 
                                  ((double)(n) * scaledlen.z);
                           rtmp = sqrt( xtmp*xtmp + ytmp*ytmp + ztmp*ztmp );

                           /*
			   *  if the child passes the MAC, then we place it
			   *  in the list of sibling interactions.
			   *  otherwise, we put it in the near interaction
			   *  list.
			   */

                           if ( MAC(radius1_2, radius1_2, rtmp, Dpmta_Theta)
				== TRUE ) {
                              index = tmp_scnt;
                              Tmp_Slist[index].x = (2*i)+l;
                              Tmp_Slist[index].y = (2*j)+m;
                              Tmp_Slist[index].z = (2*k)+n;
                              tmp_scnt++;
			   } /* if MAC [3] */

                           else {
                              index = tmp_dcnt;
                              Tmp_Dlist[index].x = (2*i)+l;
                              Tmp_Dlist[index].y = (2*j)+m;
                              Tmp_Dlist[index].z = (2*k)+n;
                              tmp_dcnt++;
			   } /* else MAC [3] */


	                } /* for n */
                     } /* for m */
                  } /* for l */

#ifndef NOPARCONV
	       } /* else MAC [2] */
#endif
	    } /* if MAC [1] */

	 } /* for k */
      } /* for j */
   } /* for i */

#ifdef SORTILIST

   /*
    *  once we have the three vectors, we nees to sort then in
    *  reverse order of their distances.  we do this by creating
    *  an array of the square of the distances to the cells and
    *  use this as a sort key on the interaction lists
    *
    */

#ifndef NOPARCONV
   /* sort the parent array */
   for (i=0; i<tmp_pcnt; i++) {
      Tmp_Sort[i] = Tmp_Plist[i].x * Tmp_Plist[i].x;
      Tmp_Sort[i] += Tmp_Plist[i].y * Tmp_Plist[i].y;
      Tmp_Sort[i] += Tmp_Plist[i].z * Tmp_Plist[i].z;
      Tmp_Sort[i] = 1000 - Tmp_Sort[i];
   } /* for i */
   Sort_Ilist(Tmp_Plist,Tmp_Sort,tmp_pcnt);
#endif

   /* sort the sibling array */
   for (i=0; i<tmp_scnt; i++) {
      Tmp_Sort[i] = Tmp_Slist[i].x * Tmp_Slist[i].x;
      Tmp_Sort[i] += Tmp_Slist[i].y * Tmp_Slist[i].y;
      Tmp_Sort[i] += Tmp_Slist[i].z * Tmp_Slist[i].z;
      Tmp_Sort[i] = 1000 - Tmp_Sort[i];
   } /* for i */
   Sort_Ilist(Tmp_Slist,Tmp_Sort,tmp_scnt);

   /* sort the direct array */
   for (i=0; i<tmp_dcnt; i++) {
      Tmp_Sort[i] = Tmp_Dlist[i].x * Tmp_Dlist[i].x;
      Tmp_Sort[i] += Tmp_Dlist[i].y * Tmp_Dlist[i].y;
      Tmp_Sort[i] += Tmp_Dlist[i].z * Tmp_Dlist[i].z;
   } /* for i */
   Sort_Ilist(Tmp_Dlist,Tmp_Sort,tmp_dcnt);

#endif

   /*
   *  at this point, we have three sorted lists of integer vectors.  these
   *  need to be converted to lists of integer interleaved cell indices.
   *  in addition we need to translate the single computed
   *  vector list into a set of eight lists, one for each child position
   *  within the parent.
   *
   *  the algorithm will proceed by cycling throught a count for each of
   *  the eight children of the parent cell.  for each list (including
   *  the parent list!) the <x,y,z> values are negated based upon the
   *  relative postion of the child within the parent.  then, this
   *  temporary value is used to construct a interleaved cell index
   *  that is then stored in the perminate list.
   *
   */

   for (i=0; i<8; i++) {

#ifndef NOPARCONV
      if ( tmp_pcnt > Dpmta_Intlist[i].psize ) {
	 tmp_ptr = (void *)Dpmta_Intlist[i].plist;
	 Dpmta_Intlist[i].plist = 
	    (int *)realloc(tmp_ptr, tmp_pcnt*sizeof(int));
	 if ( Dpmta_Intlist[i].plist == NULL ) {
	    fprintf(stderr,"ERROR: malloc() failed\n");
	    exit(-1);
	 }
	 Dpmta_Intlist[i].psize = tmp_pcnt;
      }
      Dpmta_Intlist[i].pcnt = tmp_pcnt;
#endif

      if ( tmp_scnt > Dpmta_Intlist[i].ssize ) {
	 tmp_ptr = (void *)Dpmta_Intlist[i].slist;
	 Dpmta_Intlist[i].slist = 
	    (int *)realloc(tmp_ptr, tmp_scnt*sizeof(int));
	 if ( Dpmta_Intlist[i].slist == NULL ) {
	    fprintf(stderr,"ERROR: malloc() failed\n");
	    exit(-1);
	 }
	 Dpmta_Intlist[i].ssize = tmp_scnt;
      }
      Dpmta_Intlist[i].scnt = tmp_scnt;

      if ( tmp_dcnt > Dpmta_Intlist[i].dsize ) {
	 tmp_ptr = (void *)Dpmta_Intlist[i].dlist;
	 Dpmta_Intlist[i].dlist = 
	    (int *)realloc(tmp_ptr, tmp_dcnt*sizeof(int));
	 if ( Dpmta_Intlist[i].dlist == NULL ) {
	    fprintf(stderr,"ERROR: malloc() failed\n");
	    exit(-1);
	 }
	 Dpmta_Intlist[i].dsize = tmp_dcnt;
      }
      Dpmta_Intlist[i].dcnt = tmp_dcnt;

   } /* for i */

   /*
   *  cycle through each interaction list and copy the vector from
   *  the temporary copy to the permanent storage translating
   *  the information from vector format to a single
   *  relative cell separation index.  in addition, the values of the
   *  vectors are changed to adjust for the position of the cell within
   *  the parent.
   */

   xadj = 1;
   yadj = 1;
   zadj = 1;

   for (i=0; i<2; i++) {
      for (j=0; j<2; j++) {
         for (k=0; k<2; k++) {

            /* determine cell index within parent */
            index = (k<<2) + (j<<1) + i;

	    /* process all vectors in parent ilist */
	    for (l=0; l<tmp_pcnt; l++) {

	       /* adjust the vectors depending on cell posn */
	       itmp.x = xadj * Tmp_Plist[l].x;
	       itmp.y = yadj * Tmp_Plist[l].y;
	       itmp.z = zadj * Tmp_Plist[l].z;

	       /* pack info from three vectors into a single index */
	       Vec2Sep( itmp, &(Dpmta_Intlist[index].plist[l]) );

	    } /* for l */


	    /* (repeat) process all vectors in sibling ilist */
	    for (l=0; l<tmp_scnt; l++) {

	       /* adjust the vectors depending on cell posn */
	       itmp.x = xadj * Tmp_Slist[l].x;
	       itmp.y = yadj * Tmp_Slist[l].y;
	       itmp.z = zadj * Tmp_Slist[l].z;

	       /* pack info from three vectors into a single index */
	       Vec2Sep( itmp, &(Dpmta_Intlist[index].slist[l]) );

	    } /* for l */


	    /* (repeat) process all vectors in direct ilist */
	    for (l=0; l<tmp_dcnt; l++) {

	       itmp.x = xadj * Tmp_Dlist[l].x;
	       itmp.y = yadj * Tmp_Dlist[l].y;
	       itmp.z = zadj * Tmp_Dlist[l].z;

	       /* pack info from three vectors into a single index */
	       Vec2Sep( itmp, &(Dpmta_Intlist[index].dlist[l]) );

	    } /* for l */

	    zadj *= -1;
	 } /* for k */
	 yadj *= -1;
      } /* for j */
      xadj *= -1;
   } /* for i */


#else  /* ifdef PIPED */


   int i,j,k,l,m,n;          /* loop counters */
   int dsize;                /* size of ilist arrays to allocate */
   int index, indexp;        /* temp index into ilist tables */
   int xadj, yadj, zadj;     /* terms used to adjust final vectors */
   void *tmp_ptr;            /* temp pointer */
   double rtmp;              /* temp radius value */
   double xtmp, ytmp, ztmp;  /* temp vector */
   double xtmp1, ytmp1, ztmp1; /* temp vector */
   double v2xv3dotv1, v3xv1dotv2, v1xv2dotv3;
   double mag;
   Vector v2xv3, v3xv1, v1xv2;
   Vector dtmp;              /* temp ||piped diagonal vector */
   Vector scaledv1;          /* Scaled ||piped vectors */
   Vector scaledv2;
   Vector scaledv3;
   IntVector dist;           /* length of temp cell list array */
   IntVector itmp;           /* temp vector used in calculations */

   int tmp_pcnt, tmp_scnt, tmp_dcnt;  /* length of above lists */

   Vector scaledlen;
   double radius1, radius1_2;

   /*
   *  establish how far from the parent cell we have to go in order
   *  to exceed the MAC.  this will form the limits on our search space
   *  for the child cells.
   */

   scaledv1.x = Dpmta_CellVector1.x / Dpmta_MaxCellLen;
   scaledv1.y = Dpmta_CellVector1.y / Dpmta_MaxCellLen;
   scaledv1.z = Dpmta_CellVector1.z / Dpmta_MaxCellLen;
   scaledv2.x = Dpmta_CellVector2.x / Dpmta_MaxCellLen;
   scaledv2.y = Dpmta_CellVector2.y / Dpmta_MaxCellLen;
   scaledv2.z = Dpmta_CellVector2.z / Dpmta_MaxCellLen;
   scaledv3.x = Dpmta_CellVector3.x / Dpmta_MaxCellLen;
   scaledv3.y = Dpmta_CellVector3.y / Dpmta_MaxCellLen;
   scaledv3.z = Dpmta_CellVector3.z / Dpmta_MaxCellLen;

   dtmp.x = (scaledv1.x + scaledv2.x + scaledv3.x);
   dtmp.y = (scaledv1.y + scaledv2.y + scaledv3.y);
   dtmp.z = (scaledv1.z + scaledv2.z + scaledv3.z);

   scaledlen.x = Dpmta_CV1Mag / Dpmta_MaxCellLen;
   scaledlen.y = Dpmta_CV2Mag / Dpmta_MaxCellLen;
   scaledlen.z = Dpmta_CV3Mag / Dpmta_MaxCellLen;

   radius1 = sqrt( dtmp.x*dtmp.x + 
                   dtmp.y*dtmp.y + 
                   dtmp.z*dtmp.z );
   radius1_2 = radius1 / 2.0;

   v2xv3.x = scaledv2.y*scaledv3.z - scaledv2.z*scaledv3.y;
   v2xv3.y = scaledv2.z*scaledv3.x - scaledv2.x*scaledv3.z;
   v2xv3.z = scaledv2.x*scaledv3.y - scaledv2.y*scaledv3.x;

   v3xv1.x = scaledv3.y*scaledv1.z - scaledv3.z*scaledv1.y;
   v3xv1.y = scaledv3.z*scaledv1.x - scaledv3.x*scaledv1.z;
   v3xv1.z = scaledv3.x*scaledv1.y - scaledv3.y*scaledv1.x;

   v1xv2.x = scaledv1.y*scaledv2.z - scaledv1.z*scaledv2.y;
   v1xv2.y = scaledv1.z*scaledv2.x - scaledv1.x*scaledv2.z;
   v1xv2.z = scaledv1.x*scaledv2.y - scaledv1.y*scaledv2.x;

   mag=Vec_Mag(&v2xv3);
   v2xv3dotv1 = (v2xv3.x * scaledv1.x + v2xv3.y * scaledv1.y + 
                v2xv3.z * scaledv1.z)/mag;
   mag=Vec_Mag(&v3xv1);
   v3xv1dotv2 = (v3xv1.x * scaledv2.x + v3xv1.y * scaledv2.y + 
                v3xv1.z * scaledv2.z)/mag;
   mag=Vec_Mag(&v1xv2);
   v1xv2dotv3 = (v1xv2.x * scaledv3.x + v1xv2.y * scaledv3.y + 
                v1xv2.z * scaledv3.z)/mag;
   dist.x = 1;
   dist.y = 1;
   dist.z = 1;

   while ( MAC(radius1, radius1, v2xv3dotv1 * (double)(dist.x*2), Dpmta_Theta) == FALSE )
      dist.x++;
   while ( MAC(radius1, radius1, v3xv1dotv2 * (double)(dist.y*2), Dpmta_Theta) == FALSE )
      dist.y++;
   while ( MAC(radius1, radius1, v1xv2dotv3 * (double)(dist.z*2), Dpmta_Theta) == FALSE )
      dist.z++;

   /*
   *  now allocate and initialize the temporary ilist, one entry
   *  for every possible cell in the search space.
   *  we allocate entirely too much memory (more than we could 
   *  possibly use) but we'll free it in the end.
   */

   dsize = (dist.x+1)*(dist.y+1)*(dist.z+1);

   if ( dsize > Tmp_Size ) {

#ifndef NOPARCONV
      Tmp_Plist = (IntVector *)realloc(Tmp_Plist,dsize*sizeof(IntVector));
      if ( Tmp_Plist == NULL ) {
	 fprintf(stderr,"ERROR: malloc() failed\n");
	 exit(-1);
      }
#endif

      Tmp_Slist = (IntVector *)realloc(Tmp_Slist,8*dsize*sizeof(IntVector));
      if ( Tmp_Slist == NULL ) {
	 fprintf(stderr,"ERROR: malloc() failed\n");
	 exit(-1);
      }

      Tmp_Dlist = (IntVector *)realloc(Tmp_Dlist,8*dsize*sizeof(IntVector));
      if ( Tmp_Dlist == NULL ) {
	 fprintf(stderr,"ERROR: malloc() failed\n");
	 exit(-1);
      }

#ifdef SORTILIST
      Tmp_Sort = (int *)realloc(Tmp_Sort,8*dsize*sizeof(int));
      if ( Tmp_Sort == NULL ) {
	 fprintf(stderr,"ERROR: malloc() failed\n");
	 exit(-1);
      }
#endif

      Tmp_Size = dsize;

   } /* if dsize */

   for(index=0; index<4; index++)
   {
     indexp = (~index & 0x07);
     tmp_pcnt = 0;
     tmp_scnt = 0;
     tmp_dcnt = 0;

     /* adjust centers based on child cell index */
     xadj = yadj = zadj = 1;
     switch (index)
     {
       case 1:
         xadj = -1;
         break;
       case 2:
         yadj = -1;
         break;
       case 3:
         xadj = -1;
         yadj = -1;
     }

     /*
     *  now we will build the temp interaction list for the child cell
     *  at position index within the center cell.  index=0 corresponds to
     *  a real position of <-0.5, -0.5, -0.5> wrt the center cell.
     *
     *  we can exploit symetry here and only build the interaction list
     *  for four child cells.  the remaining four children can be
     *  constructed later by simply reversing the signs of all
     *  of the <x,y,z> components of the relative interaction vectors.
     */

     /* cycle through all the parent cells */
  
     for (i=-dist.x; i<=dist.x; i++) {
        for (j=-dist.y; j<=dist.y; j++) {
           for (k=-dist.z; k<=dist.z; k++) {
  
  	    /*
  	    *  if the cell meets the parents MAC, then we do
  	    *   not have to consider it.
              */
  
              xtmp1 = 2.0 * ((double)(i) * scaledv1.x + 
  		(double)(j) * scaledv2.x + 
		(double)(k) * scaledv3.x);
              ytmp1 = 2.0 * ((double)(i) * scaledv1.y + 
  		(double)(j) * scaledv2.y + 
  		(double)(k) * scaledv3.y);
              ztmp1 = 2.0 * ((double)(i) * scaledv1.z + 
  		(double)(j) * scaledv2.z + 
  		(double)(k) * scaledv3.z);
              rtmp = sqrt( xtmp1*xtmp1 + ytmp1*ytmp1 + ztmp1*ztmp1 );
  
              if ( MAC(radius1, radius1, rtmp, Dpmta_Theta) == FALSE) {
  
#ifndef NOPARCONV
                 /*
                 *  okay, the remote cell does not meet the MAC of the
                 *  parent.
                 *
    	         *  check to see if the child cell will interact
	         *  with the entire remote cell
	         */

                 xtmp = xtmp1 + ((double)xadj * 0.5 * dtmp.x);
                 ytmp = ytmp1 + ((double)yadj * 0.5 * dtmp.y);
                 ztmp = ztmp1 + ((double)zadj * 0.5 * dtmp.z);
                 rtmp = sqrt( xtmp*xtmp + ytmp*ytmp + ztmp*ztmp );

                 if ( MAC(radius1_2, radius1, rtmp, Dpmta_Theta) == TRUE ) {
                    Tmp_Plist[tmp_pcnt].x = i;
                    Tmp_Plist[tmp_pcnt].y = j;
                    Tmp_Plist[tmp_pcnt].z = k;
                    tmp_pcnt++;
	         } /* if MAC [2] */

                 /*
	         *  else the child cannot interact with the remote
	         *  parent cell, so we will have to test each of the
	         *  children of the remote cell.
	         */

    	         else {
#endif
                    for (l=0; l<2; l++) {
                       for (m=0; m<2; m++) {
                          for (n=0; n<2; n++) {

                             xtmp = xtmp1 +  
                                ((double)xadj * (double)(l) * scaledv1.x + 
                                (double)yadj * (double)(m) * scaledv2.x + 
                                (double)zadj * (double)(n) * scaledv3.x);
                             ytmp = ytmp1 + 
                                ((double)xadj * (double)(l) * scaledv1.y + 
                                (double)yadj * (double)(m) * scaledv2.y + 
                                (double)zadj * (double)(n) * scaledv3.y);
                             ztmp = ztmp1 + (double)zadj *
                                ((double)xadj * (double)(l) * scaledv1.z + 
                                (double)yadj * (double)(m) * scaledv2.z + 
                                (double)zadj * (double)(n) * scaledv3.z);
                             rtmp = sqrt( xtmp*xtmp + ytmp*ytmp + ztmp*ztmp );

                             /*
			     *  if the child passes the MAC, then we place it
			     *  in the list of sibling interactions.
			     *  otherwise, we put it in the near interaction
			     *  list.
			     */

                             if ( MAC(radius1_2, radius1_2, rtmp, Dpmta_Theta)
		  	  	  == TRUE ) {
                                Tmp_Slist[tmp_scnt].x = (2*i)+ xadj*l;
                                Tmp_Slist[tmp_scnt].y = (2*j)+ yadj*m;
                                Tmp_Slist[tmp_scnt].z = (2*k)+ zadj*n;
                                tmp_scnt++;
			     } /* if MAC [3] */

                             else {
                                Tmp_Dlist[tmp_dcnt].x = (2*i)+ xadj*l;
                                Tmp_Dlist[tmp_dcnt].y = (2*j)+ yadj*m;
                                Tmp_Dlist[tmp_dcnt].z = (2*k)+ zadj*n;
                                tmp_dcnt++;
			     } /* else MAC [3] */


	                  } /* for n */
                       } /* for m */
                    } /* for l */

#ifndef NOPARCONV
	         } /* else MAC [2] */
#endif
	      } /* if MAC [1] */

  	   } /* for k */
        } /* for j */
     } /* for i */

#ifdef SORTILIST

     /*
      *  once we have the three vectors, we need to sort them in
      *  reverse order of their distances.  we do this by creating
      *  an array of the square of the distances to the cells and
      *  use this as a sort key on the interaction lists
      *
      */

#ifndef NOPARCONV
     /* sort the parent array */
     for (i=0; i<tmp_pcnt; i++) {
        Tmp_Sort[i] = Tmp_Plist[i].x * Tmp_Plist[i].x;
        Tmp_Sort[i] += Tmp_Plist[i].y * Tmp_Plist[i].y;
        Tmp_Sort[i] += Tmp_Plist[i].z * Tmp_Plist[i].z;
        Tmp_Sort[i] = 1000 - Tmp_Sort[i];
     } /* for i */
     Sort_Ilist(Tmp_Plist,Tmp_Sort,tmp_pcnt);
#endif

     /* sort the sibling array */
     for (i=0; i<tmp_scnt; i++) {
        Tmp_Sort[i] = Tmp_Slist[i].x * Tmp_Slist[i].x;
        Tmp_Sort[i] += Tmp_Slist[i].y * Tmp_Slist[i].y;
        Tmp_Sort[i] += Tmp_Slist[i].z * Tmp_Slist[i].z;
        Tmp_Sort[i] = 1000 - Tmp_Sort[i];
     } /* for i */
     Sort_Ilist(Tmp_Slist,Tmp_Sort,tmp_scnt);
  
     /* sort the direct array */
     for (i=0; i<tmp_dcnt; i++) {
        Tmp_Sort[i] = Tmp_Dlist[i].x * Tmp_Dlist[i].x;
        Tmp_Sort[i] += Tmp_Dlist[i].y * Tmp_Dlist[i].y;
        Tmp_Sort[i] += Tmp_Dlist[i].z * Tmp_Dlist[i].z;
     } /* for i */
     Sort_Ilist(Tmp_Dlist,Tmp_Sort,tmp_dcnt);

#endif

     /*
     *  at this point, we have three sorted lists of integer vectors.  these
     *  need to be converted to lists of integer interleaved cell indices.
     *  in addition we need to translate the single computed
     *  vector list into a set of eight lists, one for each child position
     *  within the parent.
     *
     *  the algorithm will proceed by cycling throught a count for each of
     *  the eight children of the parent cell.  for each list (including
     *  the parent list!) the <x,y,z> values are negated based upon the
     *  relative postion of the child within the parent.  then, this
     *  temporary value is used to construct a interleaved cell index
     *  that is then stored in the perminate list.
     *
     */


#ifndef NOPARCONV
     if ( tmp_pcnt > Dpmta_Intlist[index].psize ) {
       tmp_ptr = (void *)Dpmta_Intlist[index].plist;
       Dpmta_Intlist[index].plist = 
          (int *)realloc(tmp_ptr, tmp_pcnt*sizeof(int));
       if ( Dpmta_Intlist[index].plist == NULL ) {
         fprintf(stderr,"ERROR: malloc() failed\n");
         exit(-1);
       }
       Dpmta_Intlist[index].psize = tmp_pcnt;
     }
     Dpmta_Intlist[index].pcnt = tmp_pcnt;

     if ( tmp_pcnt > Dpmta_Intlist[indexp].psize ) {
       tmp_ptr = (void *)Dpmta_Intlist[indexp].plist;
       Dpmta_Intlist[indexp].plist = 
          (int *)realloc(tmp_ptr, tmp_pcnt*sizeof(int));
       if ( Dpmta_Intlist[indexp].plist == NULL ) {
         fprintf(stderr,"ERROR: malloc() failed\n");
         exit(-1);
       }
       Dpmta_Intlist[indexp].psize = tmp_pcnt;
     }
     Dpmta_Intlist[indexp].pcnt = tmp_pcnt;
#endif

     if ( tmp_scnt > Dpmta_Intlist[index].ssize ) {
       tmp_ptr = (void *)Dpmta_Intlist[index].slist;
       Dpmta_Intlist[index].slist = 
          (int *)realloc(tmp_ptr, tmp_scnt*sizeof(int));
       if ( Dpmta_Intlist[index].slist == NULL ) {
         fprintf(stderr,"ERROR: malloc() failed\n");
         exit(-1);
       }
       Dpmta_Intlist[index].ssize = tmp_scnt;
     }
     Dpmta_Intlist[index].scnt = tmp_scnt;

     if ( tmp_scnt > Dpmta_Intlist[indexp].ssize ) {
       tmp_ptr = (void *)Dpmta_Intlist[indexp].slist;
       Dpmta_Intlist[indexp].slist = 
          (int *)realloc(tmp_ptr, tmp_scnt*sizeof(int));
       if ( Dpmta_Intlist[indexp].slist == NULL ) {
         fprintf(stderr,"ERROR: malloc() failed\n");
         exit(-1);
       }
       Dpmta_Intlist[indexp].ssize = tmp_scnt;
     }
     Dpmta_Intlist[indexp].scnt = tmp_scnt;

     if ( tmp_dcnt > Dpmta_Intlist[index].dsize ) {
       tmp_ptr = (void *)Dpmta_Intlist[index].dlist;
       Dpmta_Intlist[index].dlist = 
           (int *)realloc(tmp_ptr, tmp_dcnt*sizeof(int));
       if ( Dpmta_Intlist[index].dlist == NULL ) {
         fprintf(stderr,"ERROR: malloc() failed\n");
         exit(-1);
       }
       Dpmta_Intlist[index].dsize = tmp_dcnt;
     }
     Dpmta_Intlist[index].dcnt = tmp_dcnt;

     if ( tmp_dcnt > Dpmta_Intlist[indexp].dsize ) 
     {
       tmp_ptr = (void *)Dpmta_Intlist[indexp].dlist;
       Dpmta_Intlist[indexp].dlist = 
           (int *)realloc(tmp_ptr, tmp_dcnt*sizeof(int));
       if ( Dpmta_Intlist[indexp].dlist == NULL ) 
       {
         fprintf(stderr,"ERROR: malloc() failed\n");
         exit(-1);
       }
       Dpmta_Intlist[indexp].dsize = tmp_dcnt;
     }
     Dpmta_Intlist[indexp].dcnt = tmp_dcnt;


     /*
     *  cycle through each interaction list and copy the vector from
     *  the temporary copy to the permanent storage translating
     *  the information from vector format to a single
     *  relative cell separation index.  in addition, the values of the
     *  vectors are changed to adjust for the position of the cell within
     *  the parent.
     */


    /* process all vectors in parent ilist */
    for (l=0; l<tmp_pcnt; l++) {

       itmp.x = Tmp_Plist[l].x;
       itmp.y = Tmp_Plist[l].y;
       itmp.z = Tmp_Plist[l].z;

       /* pack info from three vectors into a single index */
       Vec2Sep( itmp, &(Dpmta_Intlist[index].plist[l]) );

       itmp.x *= -1;
       itmp.y *= -1;
       itmp.z *= -1;

       /* pack info from three vectors into a single index */
       Vec2Sep( itmp, &(Dpmta_Intlist[indexp].plist[l]) );

    } /* for l */


    /* (repeat) process all vectors in sibling ilist */
    for (l=0; l<tmp_scnt; l++) {

       itmp.x = Tmp_Slist[l].x;
       itmp.y = Tmp_Slist[l].y;
       itmp.z = Tmp_Slist[l].z;

       /* pack info from three vectors into a single index */
       Vec2Sep( itmp, &(Dpmta_Intlist[index].slist[l]) );

       itmp.x *= -1;
       itmp.y *= -1;
       itmp.z *= -1;

       /* pack info from three vectors into a single index */
       Vec2Sep( itmp, &(Dpmta_Intlist[indexp].slist[l]) );

    } /* for l */


    /* (repeat) process all vectors in direct ilist */
    for (l=0; l<tmp_dcnt; l++) {

       itmp.x = Tmp_Dlist[l].x;
       itmp.y = Tmp_Dlist[l].y;
       itmp.z = Tmp_Dlist[l].z;

       /* pack info from three vectors into a single index */
       Vec2Sep( itmp, &(Dpmta_Intlist[index].dlist[l]) );

       itmp.x *= -1;
       itmp.y *= -1;
       itmp.z *= -1;

       /* pack info from three vectors into a single index */
       Vec2Sep( itmp, &(Dpmta_Intlist[indexp].dlist[l]) );

    } /* for l */

  } /* for ic */

#endif  /* ifdef PIPED */

}  /* Make_Ilist() */


/****************************************************************
 *
 *  Init_Ilist() -
 *
 *  this routine performs once only initialization of the global 
 *  data structures and constants used in creating and updating the
 *  interaction lists
 *
 */

void Init_Ilist()
{
   int i;

   /*
   *  first start off by allocating out the final Ilist data structure.
   */

   Dpmta_Intlist = (IlistPtr)malloc(8*sizeof(Ilist));
   if ( Dpmta_Intlist == NULL ) {
      fprintf(stderr,"ERROR: Init_Ilist() - malloc() #1 failed\n");
      exit(-1);
   }

   for (i=0; i<8; i++) {
#ifndef NOPARCONV
      Dpmta_Intlist[i].psize = 0;
      Dpmta_Intlist[i].plist = (int *)NULL;
      Dpmta_Intlist[i].pcnt = 0;
#endif
      Dpmta_Intlist[i].ssize = 0;
      Dpmta_Intlist[i].slist = (int *)NULL;
      Dpmta_Intlist[i].scnt = 0;

      Dpmta_Intlist[i].dsize = 0;
      Dpmta_Intlist[i].dlist = (int *)NULL;
      Dpmta_Intlist[i].dcnt = 0;

   }

#ifndef NOPARCONV
   Tmp_Plist = (IntVector *)NULL;
#endif
   Tmp_Slist = (IntVector *)NULL;
   Tmp_Dlist = (IntVector *)NULL;
   Tmp_Size = 0;

#ifdef SORTILIST
   Tmp_Sort = (int *)NULL;
#endif

} /* Init_Ilist() */


/****************************************************************
 *
 *  Delete_Ilist() -
 *
 *  for completeness, we have added a routine that will free
 *  all allocated data structures that are malloc'd within the
 *  scope of this module.
 *
 *  this will effectively destroy the interaction lists, which
 *  could be possibly reinitialized by another call to Init_Ilist().
 *
 */

void Delete_Ilist()
{
   int i;

   /*
    *  de-alloc the interaction lists
    */

   for (i=0; i<8; i++ ) {
#ifndef NOPARCONV
      free(Dpmta_Intlist[i].plist);
#endif
      free(Dpmta_Intlist[i].slist);
      free(Dpmta_Intlist[i].dlist);
   }
   free(Dpmta_Intlist);

   /*
    *  free up the temporary data structures
    */

#ifndef NOPARCONV
   free(Tmp_Plist);
#endif
   free(Tmp_Slist);
   free(Tmp_Dlist);

#ifdef SORTILIST
   free(Tmp_Sort);
#endif

} /* Delete_Ilist */


#ifdef SORTILIST
/****************************************************************
 *
 *  Sort_Ilist() - 
 *
 *  this routine implements a simple heap sort on the interaction list.
 *  it is based upon the heap sort presented in the first edition of
 *  'Numerical Recipes in C' page 247.  it is slightly modified to 
 *  use C-style indexing (array indexes begin at 0) as well as perform
 *  a reverse sort (highest first).
 *
 */

void Sort_Ilist(
   IntVector *il,       /* array of integer vectors */
   int *sl,             /* integer array used for sorting criteria */
   int n)               /* length of array */
{

   int i, j, l, ir;
   int rra;
   IntVector rrb;


   l = (n>>1);
   ir = n-1;
   
   for (;;) {
      if (l>0) {
	 --l;
	 rra = sl[l];
	 rrb.x = il[l].x;
	 rrb.y = il[l].y;
	 rrb.z = il[l].z;
      }
      else {
	 rra = sl[ir];
	 rrb.x = il[ir].x;
	 rrb.y = il[ir].y;
	 rrb.z = il[ir].z;
	 sl[ir] = sl[0];
	 il[ir].x = il[0].x;
	 il[ir].y = il[0].y;
	 il[ir].z = il[0].z;

	 if (--ir == 0) {
	    sl[0] = rra;
	    il[0].x = rrb.x;
	    il[0].y = rrb.y;
	    il[0].z = rrb.z;
            return;
	 }
      } /* else */

      i = l;
      j = ((l+1)<<1) - 1;

      while ( j <= ir ) {
	 if ( (j < ir) && (sl[j] >= sl[j+1]) ) {
	    ++j;
	 }
	 if (rra >= sl[j]) {
	    sl[i] = sl[j];
	    il[i].x = il[j].x;
	    il[i].y = il[j].y;
	    il[i].z = il[j].z;
            i = j;
	    j += i;
	 }
	 else {
	    j = ir+1;
	 }
      } /* while */

      sl[i] = rra;
      il[i].x = rrb.x;
      il[i].y = rrb.y;
      il[i].z = rrb.z;
   } /* for */

} /* Sort_Ilist() */

#endif


/****************************************************************
*
*  dumps the ilist to a file
*/

void Dump_Ilist() 
{

   int i,j;
   IntVector iv;
   FILE *fp;
   char filename[80];

   sprintf(filename,"/tmp/ilist.pid%d",Dpmta_Pid);
   fp = fopen(filename,"w");

   fprintf(fp," Interaction list, Theta = %f\n\n",Dpmta_Theta);

   for (j=0; j<4; j++) {

#ifndef NOPARCONV
      fprintf(fp,"Parent Ilist %d [%d]\n",j, Dpmta_Intlist[j].pcnt);
      for (i=0; i<Dpmta_Intlist[j].pcnt; i++) {
	 Sep2Vec(Dpmta_Intlist[j].plist[i],&iv);
	 fprintf(fp,"p%d (%d,%d,%d)\n",j,iv.x,iv.y,iv.z);
      }
      fprintf(fp,"================================\n");
#endif

      fprintf(fp,"Sibling Ilist %d [%d]\n", j, Dpmta_Intlist[j].scnt);
      for (i=0; i<Dpmta_Intlist[j].scnt; i++) {
	 Sep2Vec(Dpmta_Intlist[j].slist[i],&iv);
	 fprintf(fp,"s%d (%d,%d,%d)\n",j,iv.x,iv.y,iv.z);
      }
      fprintf(fp,"================================\n");

      fprintf(fp,"Direct Ilist %d [%d]\n",j, Dpmta_Intlist[j].dcnt);
      for (i=0; i<Dpmta_Intlist[j].dcnt; i++) {
	 Sep2Vec(Dpmta_Intlist[j].dlist[i],&iv);
	 fprintf(fp,"d%d (%d,%d,%d)\n",j,iv.x,iv.y,iv.z);
      }
      fprintf(fp,"================================\n");

   }
   fclose(fp);

} /* Dump_Ilist */
