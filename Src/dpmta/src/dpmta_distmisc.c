/* 
*  dpmta_distmisc.c - misc routines to perform mapping between data and
*    distributed processed id.
*
*  these routines are used by both the master and slave processes.  thus
*  they cannot access any global data structures.  most of these
*  functions were taken from dpmta_slvmisc.c
*
*
*  w. t. rankin
*
*  Copyright (c) 1997 Duke University
*  All rights reserved
*
*/

static char rcsid[] = "$Id: dpmta_distmisc.c,v 2.3 1998/04/29 18:36:21 wrankin Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_distmisc.c,v $
 * Revision 2.3  1998/04/29 18:36:21  wrankin
 * fixed code that creates RLcell counter - now works with row/col index
 *
 * Revision 2.2  1998/04/01 20:08:04  wrankin
 * added support for HILBERT ordering
 * general cleanup of code structure (more OOP if you will)
 *
 * Revision 2.1  1997/09/29 20:24:52  wrankin
 * fixed problem with invalid (empty) multipoles during upward pass.
 * cell indexing by processor was inconsistant between master/slave.
 *
 *
 */

/* include files */
#include <stdlib.h>
#include <stdio.h>


/* internal prototypes */

int hil2mort( int, int );
int mort2hil( int, int );

int cart2mort( int*, int );
void mort2cart( int, int, int* );

int rco2mort( int, int );
int mort2rco( int, int );


/* local global data structures */

static int **I2C_Array;
static int **C2I_Array;


/****************************************************************
*
*  Dist_Init() - allocate and initialize internal data structures
*
*/

void Dist_Init( int numlevels )
{

   int i,j;
   int numcells;

   /* allocate arrays for mapping hilbert to morton */
   
   I2C_Array = (int **)malloc(numlevels*sizeof(int *));
   if ( I2C_Array == (int **)NULL ) {
      fprintf(stderr,"ERROR: Dist_Init() - malloc failed\n");
      exit(-1);
   }

   for ( i=0; i<numlevels; i++ ) {
      numcells = 0x1<<(3*i);
      I2C_Array[i] = (int *)malloc(numcells*sizeof(int));
      if ( I2C_Array[i] == (int *)NULL ) {
	 fprintf(stderr,"ERROR: Dist_Init() - malloc failed\n");
	 exit(-1);
      }
   } /* for i */
   
   C2I_Array = (int **)malloc(numlevels*sizeof(int *));
   if ( C2I_Array == (int **)NULL ) {
      fprintf(stderr,"ERROR: Dist_Init() - malloc failed\n");
      exit(-1);
   }

   for ( i=0; i<numlevels; i++ ) {
      numcells = 0x1<<(3*i);
      C2I_Array[i] = (int *)malloc(numcells*sizeof(int));
      if ( C2I_Array[i] == (int *)NULL ) {
	 fprintf(stderr,"ERROR: Dist_Init() - malloc failed\n");
	 exit(-1);
      }
   } /* for i */

   for ( i=0; i<numlevels; i++ ) {
      numcells = 0x1<<(i*3);
      for ( j=0; j<numcells; j++ ) {
#ifdef HILBERT
	 I2C_Array[i][j] = hil2mort(j,i);
	 C2I_Array[i][j] = mort2hil(j,i);
#endif
#ifdef ROWCOL
	 I2C_Array[i][j] = rco2mort(j,i);
	 C2I_Array[i][j] = mort2rco(j,i);
#endif
      } /* for j */
   } /* for i */

} /* Dist_Init() */


/****************************************************************
*
*  Dist_Delete() - free up internal data structures
*
*/

void Dist_Delete( int numlevels )
{
   int i;

   for ( i=0; i<numlevels; i++ ) {
      free(I2C_Array[i]);
      free(C2I_Array[i]);
   } /* for i */
   free(I2C_Array);
   free(C2I_Array);
   
} /* Dist_Delete() */


/****************************************************************
*
*  the following routines are used in indexing the cell table
*  structure.  therefor, all the inputs and outputs are for
*  Morton ordering.
*/

/****************************************************************
*
*  returns the parent cell id for any given cell
*/

int getparent( int cellid )
{
   return (cellid >> 3);
} /* get parent */


/****************************************************************
*
*  returns the lowest indexed child cell id for any given cell
*/

int getfirstchild( int cellid )
{
   return (cellid << 3);
} /* get first child */


/****************************************************************
 *
 *  the following routines map the problem/tree space onto the
 *  processor allocation space.   as far as the calling procesures
 *  are aware, processor allocation is represented by a a mapping
 *  of cells onto a contiguous section of a linear array.
 *
 *  the calling program can make no assumptions as to the
 *  method of mapping for the individual cell structures.
 *
 */


/****************************************************************
*
*  getslvpid() - returns the pid for the slave that owns
*    a given cell
*/

int getslvpid(int np, int level, int cell)
{
   int id;

#ifdef HILBERT
   id = C2I_Array[level][cell];
#else
#ifdef ROWCOL
   id = C2I_Array[level][cell];
#else
   id = cell;
#endif
#endif
   
   return ( (np * id) / (0x01 <<(3*level)) );
}

/****************************************************************
*
*  getslvpid_indx() - returns the pid for the slave that owns
*    a given index
*/

int getslvpid_indx(int np, int level, int id)
{
   return ( (np * id) / (0x01 <<(3*level)) );
}


/*****************************************************************
*
*  returns the first cell that a specific pid owns for a given
*  level.
*/

int getsindex( int np, int pid, int level )
{
   int sc, ec;

   sc = (pid * (0x1 << (3*level)) + (np-1)) / np ;
   ec = ((pid+1) * (0x1 << (3*level)) + (np-1)) / np - 1 ;

   if ( sc > ec )
      return(-1);
   else
      return(sc);
}


/*****************************************************************
*
*  returns the last cell that a specific pid owns for a given
*  level.
*/

int geteindex( int np, int pid, int level )
{
   int sc, ec;

   sc = (pid * (0x1 << (3*level)) + (np-1)) / np ;
   ec = ((pid+1) * (0x1 << (3*level)) + (np-1)) / np - 1 ;

   if ( sc > ec )
      return(-1);
   else
      return(ec);
}


/*****************************************************************
*
*  returns the cell id for a given index
*
*/

int index2cell( int index, int level )
{
   
#ifdef HILBERT
   return I2C_Array[level][index];
   /* hil2mort( index, level+1 ); */
#else
#ifdef ROWCOL
   return I2C_Array[level][index];
#else
   return index;
#endif
#endif

} /* index2cell */


/*****************************************************************
*
*  returns the index for a given cell id.
*
*/

int cell2index( int cellid, int level )
{

#ifdef HILBERT
   return C2I_Array[level][cellid];
   /* mort2hil( index, level+1 ); */
#else
#ifdef ROWCOL
   return C2I_Array[level][cellid];
#else
   return cellid;
#endif
#endif
   
} /* cell2index */


/****************************************************************
*
*  the following routines implement the translation between
*  Hilbert and Morton encoding schemes.  they should never be
*  called by any external routines.
*
*/

/****************************************************************
*
*  mort2hil() - convert a morton coordinate to the corresponsding
*    hilbert index
*
*  parameters -
*    mort - input morton coordinate
*    level - number of levels of decomposition
*
*  returns -
*    hilbert index cooresponding to the morton number
*
*
*/

int mort2hil( int mort, int level )
{
   int i,j;
   int mask1, mask2;
   int hil1, hil2, hil3;
   int curr, cur[3];
   int t[3][3];
   int tmp;
   int width;

   /*
    * initialize xform matrix
    */

   for ( i=0; i<3; i++ ) {
      for ( j=0; j<3; j++ ) {
	 if ( i == j ) {
	    t[i][j] = 1;
	 }
	 else {
	    t[i][j] = 0;
	 }
      } /* for j */
   } /* for i */

   
   /*
    * grab the next top [3] bits
    */

   width = 3 * ( level );
   mask1 = ( 0x1 << 3 ) - 1;
   hil3 = 0;
   
   while ( width >= 0 ) {

      curr = (mort >> width) & mask1;

      /*
      * apply xform matrix to curr (morton), then convert
      * the shifted morton to an inverse grey code and save.
      */

      tmp = curr;
      for ( i=0; i<3; i++ ) {
         cur[i] = tmp & 0x1;
	 tmp >>= 0x1;
      } /* for i */

      hil1 = 0;
      for ( i=2; i>=0; i-- ) {
	 tmp = 0;
	 for ( j=2; j>=0; j-- ) {
  	    tmp |=  ((cur[j] & t[i][j]) ^ ((t[i][j]>>1) & t[i][j]));
	 } /* for j */
	 hil1 = (hil1 << 1) | tmp;
      } /* for i */


      /* convert recent dimension to an inverse grey code */

      hil2 = hil1;
      mask2 = 0x1 << 3;
      while ( mask2 > 1 ) {
	 if ( mask2 & hil2 ) {
	    hil2 ^= (mask2 >> 1);
	 }
	 mask2 >>= 1;
      } /* while mask2 */

      /* accumulate new sub-coordinate in result */

      hil3 = (hil3 << 3) | hil2;
      
      /*
      * based upon position in curr, modify xform
      * matrix by multiplying it by the next incremental
      * transform dependant upon the current relative location
      * within the subcell.
      */

      /*
       * aka. we started with a morton coordinate,
       * and converted it to an inverse grey code to get
       * the corresponding
       * hilbert coordinate within the sub-cell.
       * the resulting subcell is rotated depending upon
       * its position within the larger cell.
       */

      if ( width > 0 ) {
	 switch (hil2) {

	 case 0x00:
	    /* reflect around the x = z plane */
	    tmp = t[0][0];
	    t[0][0] = t[2][0];
	    t[2][0] = tmp;
	    tmp = t[0][1];
	    t[0][1] = t[2][1];
	    t[2][1] = tmp;
	    tmp = t[0][2];
	    t[0][2] = t[2][2];
	    t[2][2] = tmp;
	    break;

	 case 0x01:
	 case 0x02:
	    /* reflect around y=z plane */
	    tmp = t[1][0];
	    t[1][0] = t[0][0];
	    t[0][0] = t[2][0];
	    t[2][0] = tmp;
	    tmp = t[1][1];
	    t[1][1] = t[0][1];
	    t[0][1] = t[2][1];
	    t[2][1] = tmp;
	    tmp = t[1][2];
	    t[1][2] = t[0][2];
	    t[0][2] = t[2][2];
	    t[2][2] = tmp;
	    break;
	 
	 case 0x03:
	 case 0x04:
	    /* reflect around the x=-y plane*/
	    tmp = t[0][0] ^ 0x2;
	    t[0][0] = t[1][0] ^ 0x2;
	    t[1][0] = tmp;
	    tmp = t[0][1] ^ 0x2;
	    t[0][1] = t[1][1] ^ 0x2;
	    t[1][1] = tmp;
	    tmp = t[0][2] ^ 0x2;
	    t[0][2] = t[1][2] ^ 0x2;
	    t[1][2] = tmp;
	    break;

	 case 0x05:
	 case 0x06:
	    /* inversion of 0x03 around y=-y, z=-z */
	    tmp = t[2][0] ^ 0x2;
	    t[2][0] = t[1][0] ^ 0x2;
	    t[1][0] = t[0][0];
	    t[0][0] = tmp;
	    tmp = t[2][1] ^ 0x2;
	    t[2][1] = t[1][1] ^ 0x2;
	    t[1][1] = t[0][1];
	    t[0][1] = tmp;
	    tmp = t[2][2] ^ 0x2;
	    t[2][2] = t[1][2] ^ 0x2;
	    t[1][2] = t[0][2];
	    t[0][2] = tmp;
	    break;
	      
	 case 0x07:
	    /* inversion of 0x00 around x=-x, z=-z */
	    tmp = t[0][0] ^ 0x2;
	    t[0][0] = t[2][0] ^ 0x2;
	    t[2][0] = tmp;
	    tmp = t[0][1] ^ 0x2;
	    t[0][1] = t[2][1] ^ 0x2;
	    t[2][1] = tmp;
	    tmp = t[0][2] ^ 0x2;
	    t[0][2] = t[2][2] ^ 0x2;
	    t[2][2] = tmp;
	    break;

	 } /* switch curr */

      } /* if width */
      
      width -= 3;

   }  /* while */

   return hil3;

} /* mort2hil */


/****************************************************************
*
*  hil2mort() - convert a hilbert index to  the corresponsding
*     morton coordinate.
*
*  parameters -
*    hil - input hilber number
*    level - number of levesl of decomposition
*
*  returns -
*    morton order number
*
*/


int hil2mort( int hil, int level )
{
   
   int i,j;
   int mask1, mask2;
   int mort1, mort2;
   int curr, cur[3];
   int t[3][3];
   int tmp;
   int width;

   /*
    * initialize xform matrix
    */

   for ( i=0; i<3; i++ ) {
      for ( j=0; j<3; j++ ) {
	 if ( i == j ) {
	    t[i][j] = 1;
	 }
	 else {
	    t[i][j] = 0;
	 }
      } /* for j */
   } /* for i */

   
   /*
    * grab the next top [3] bits
    */

   width = 3 * ( level );
   mask1 = ( 0x1 << 3 ) - 1;
   mort2 = 0;

   while ( width >= 0 ) {

      curr = (hil >> width) & mask1;

      /*
       * convert current to grey code
       */

      mort1 = curr;
      mask2 = 0x1 << 3;
      while ( mask2 > 1 ) {
	 if ( mask2 & curr ) {
	    mort1 ^= (mask2 >> 1);
	 }
	 mask2 >>= 1;
      }
      curr = mort1;

      /*
      * apply xform matrix to curr (grey code) and
      * save
      */

      tmp = curr;
      for ( i=0; i<3; i++ ) {
         cur[i] = tmp & 0x1;
	 tmp >>= 0x1;
      }

      for ( i=2; i>=0; i-- ) {
	 tmp = 0;
	 for ( j=2; j>=0; j-- ) {
 	    tmp |=  ((cur[j] & t[i][j]) ^ ((t[i][j]>>1) & t[i][j]));
	 } /* for j */
	 mort2 = (mort2 << 1) | tmp;
      } /* for i */

      /*
      *  based upon position in curr, modify xform
      *  matrix by multiplying it by the next incremental
      *  transform dependant upon the current relative location
      *  within the subcell.
      *
      *  aka. we started with a hilbert coordinate,
      *  and converted it to a grey code to get the corresponding
      *  morton coordinate within the sub-cell.
      *  the resulting subcell is rotated depending upon
      *  its position within the larger cell.
      */

      if ( width > 0 ) {
	 switch (curr) {

	 case 0x00:
	    /* reflect around the x = z plane */
	    tmp = t[0][0];
	    t[0][0] = t[0][2];
	    t[0][2] = tmp;
	    tmp = t[1][0];
	    t[1][0] = t[1][2];
	    t[1][2] = tmp;
	    tmp = t[2][0];
	    t[2][0] = t[2][2];
	    t[2][2] = tmp;
	    break;

	 case 0x01:
	 case 0x03:
	    /*
	     * reflect around x=z plane,
	     * followed by reflection around x=y
	     */
	    tmp = t[0][1];
	    t[0][1] = t[0][0];
	    t[0][0] = t[0][2];
	    t[0][2] = tmp;
	    tmp = t[1][1];
	    t[1][1] = t[1][0];
	    t[1][0] = t[1][2];
	    t[1][2] = tmp;
	    tmp = t[2][1];
	    t[2][1] = t[2][0];
	    t[2][0] = t[2][2];
	    t[2][2] = tmp;
	    break;
	 
	 case 0x02:
	 case 0x06:
	    /* reflect around the x=-y plane*/
	    tmp = t[0][0] ^ 0x2;
	    t[0][0] = t[0][1] ^ 0x2;
	    t[0][1] = tmp;
	    tmp = t[1][0] ^ 0x2;
	    t[1][0] = t[1][1] ^ 0x2;
	    t[1][1] = tmp;
	    tmp = t[2][0] ^ 0x2;
	    t[2][0] = t[2][1] ^ 0x2;
	    t[2][1] = tmp;
	    break;

	 case 0x05:
	 case 0x07:
	    /*
	     * inversion of 0x03 around y=-z
	     * followed by reflection around x=-z
	     */
	    tmp = t[0][2] ^ 0x2;
	    t[0][2] = t[0][1] ^ 0x2;
	    t[0][1] = t[0][0];
	    t[0][0] = tmp;
	    tmp = t[1][2] ^ 0x2;
	    t[1][2] = t[1][1] ^ 0x2;
	    t[1][1] = t[1][0];
	    t[1][0] = tmp;
	    tmp = t[2][2] ^ 0x2;
	    t[2][2] = t[2][1] ^ 0x2;
	    t[2][1] = t[2][0];
	    t[2][0] = tmp;
	    break;
	      
	 case 0x04:
	    /* inversion of 0x00 around x=-x, z=-z */
	    tmp = t[0][0] ^ 0x2;
	    t[0][0] = t[0][2] ^ 0x2;
	    t[0][2] = tmp;
	    tmp = t[1][0] ^ 0x2;
	    t[1][0] = t[1][2] ^ 0x2;
	    t[1][2] = tmp;
	    tmp = t[2][0] ^ 0x2;
	    t[2][0] = t[2][2] ^ 0x2;
	    t[2][2] = tmp;
	    break;

	 } /* switch curr */

      } /* if width */
      
      width -= 3;

   } /* while */

   return mort2;

} /* hil2mort */


/**************************************************************** 
* 
*  mort2cart - translate morton ordered sequence into a 3-d cartesian
*  coordinate
*
*/

void mort2cart( int mort, int width, int* cart )
{

   int i;
   int x,y,z;
   int mask;

   x = 0;
   y = 0;
   z = 0;

   mask = 0x1;

   for ( i=0; i<width; i++ ) {
      x |= (mask & mort);
      mort >>= 1;
      y |= (mask & mort);
      mort >>= 1;
      z |= (mask & mort);
      mask <<= 1;
   } /* for i */

   cart[0] = x;
   cart[1] = y;
   cart[2] = z;

} /* mort2cart() */

/**************************************************************** 
* 
*  cart2mort - translate 3d cartesian coordinates to a morton
*   ordered sequence 
*
*/

int cart2mort( int *cart, int width ) {

   int i;
   int x,y,z;
   int mask;
   int mort;

   x = cart[0];
   y = cart[1];
   z = cart[2];

   mort = 0;
   mask = 0x1;
   y <<= 1;
   z <<= 2;
   for ( i=0; i<width; i++ ) {
      mort |= x & mask;
      mask <<= 1;
      mort |= y & mask;
      mask <<= 1;
      mort |= z & mask;
      mask <<= 1;
      x <<= 2;
      y <<= 2;
      z <<= 2;
   } /* for i */

   return ( mort );

} /* cart2mort() */


/****************************************************************
*
*  rco2mort() - translate from row/col ordering to morton
*    ordering (aka. intereaved bits)
*
*    ie - zzzzyyyyxxxx -> zyxzyxzyxzyx
*
*/

int rco2mort( int rco, int width ) {

   int i;
   int x,y,z;
   int mask;
   int mort;

   mask = ( 0x1 << width ) - 1;
   x = rco & mask;
   rco >>= width;
   y = rco & mask;
   rco >>= width;
   z = rco & mask;

   mask = 0x1;
   mort = 0;

   y <<= 1;
   z <<= 2;
   for ( i=0; i<width; i++ ) {
      mort |= x & mask;
      mask <<= 1;
      mort |= y & mask;
      mask <<= 1;
      mort |= z & mask;
      mask <<= 1;
      x <<= 2;
      y <<= 2;
      z <<= 2;
   } /* for i */

   return ( mort );

} /* rco2mort() */


/****************************************************************
*
*  mort2rco() - maps morton ordered index ro row-colume
*
*    ie - zyxzyxzyxzyx -> zzzzyyyyxxxx
*
*/

int mort2rco( int mort, int width ) {

   int i;
   int x,y,z;
   int mask;
   int rco;

   x = 0;
   y = 0;
   z = 0;

   mask = 0x1;

   for ( i=0; i<width; i++ ) {
      x |= (mask & mort);
      mort >>= 1;
      y |= (mask & mort);
      mort >>= 1;
      z |= (mask & mort);
      mask <<= 1;
   } /* for i */

   mask = (0x1 << width) - 1;
   rco = z & mask;
   rco <<= width;
   rco |= y & mask;
   rco <<= width;
   rco |= x & mask;

   return (rco);

} /* mort2rco() */

