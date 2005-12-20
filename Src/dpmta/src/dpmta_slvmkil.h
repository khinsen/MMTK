/**********************************************************************
*
*  dpmta_slvmkil.h - macro definitions and prototyped for relative
*    interaction lists.
*
*  w.t.rankin
*
*  Copyright (c) 1995 Duke University
*  All rights reserved
*
*/

/*
 *  RCSid: $Id: dpmta_slvmkil.h,v 2.4 1997/11/07 16:49:48 wrankin Exp $
 *
 *  RCS history:
 * 
 *  $Log: dpmta_slvmkil.h,v $
 *  Revision 2.4  1997/11/07 16:49:48  wrankin
 *  massive cleanup of code.
 *   - ansi-fication and inclusion of prototypes
 *   - removed unused variables
 *   - all (except the test) code compiles with minimal warnings under gcc.
 *
 *  Revision 2.3  1995/12/08 23:00:36  wrankin
 *  preliminary release of DPMTA 2.3
 *    - added working Periodic Boundary Conditions (PDC) flag
 *    - added code for Virial Computation (not working yet)
 *    - general cleanup of some modules
 *
 * Revision 2.2  1995/11/29  22:29:33  wrankin
 * addition of periodic boundary code
 * beginning work on virial code.
 *
 * Revision 2.1  1995/07/10  02:47:19  wrankin
 * multipole processing code modified to use precomputed mpe transfer
 *   matrices.
 *
 * in addition, all multippole calculation routines have been removed
 *   from this source distribution and will be placed in a separate
 *   multipole library to be supplied externally.
 *
*/


/*
 *  prototypes 
 */

void Init_Ilist();
void Make_Ilist();
void Delete_Ilist();

int MAC( double, double, double, double );
int Cell2Cell( int, int, int, int *, int *);
int Vec2Sep( IntVector, int * );
int Sep2Vec( int, IntVector * );


/****************************************************************
*
*  Cell2Cell - translation routines macro
*
*/

/* int Level - cell level */
/* int Cell - local cell id */
/* int Sep - cell separation index */
/* int Outcell - remote cell id */
/* int Vflag - overflow flag */

#define CELL2CELL(Level,Cell,Sep,Outcell,Vflag) \
{ \
   int macro_nlvlmask; \
   int macro_tmpout; \
 \
   macro_nlvlmask = -1 << (3*Level); \
   macro_tmpout =  ((Cell | ~ILMASK1) + (ILMASK1 & Sep)) & ILMASK1; \
   macro_tmpout |= ((Cell | ~ILMASK2) + (ILMASK2 & Sep)) & ILMASK2; \
   macro_tmpout |= ((Cell | ~ILMASK3) + (ILMASK3 & Sep)) & ILMASK3; \
   Vflag = macro_tmpout & macro_nlvlmask; \
   if ( Dpmta_PBC == 1 ) \
      Vflag = 0; \
   Outcell = macro_tmpout & ~macro_nlvlmask; \
 \
}



/****************************************************************
*
*  Vec2Sep - translation macro
*
*  this function takes a  vector offset between two cells and 
*  returns a cell separation index
*
*  the algorithm is based upon section 3.2 of elliotts disertation
*  with some enhancements.
*
*/

/* IntVector Offset - remote cell id */
/* int Sep - cell separation index */

#define VEC2SEP(Offset,Sep) \
{ \
 \
   int macro_i; \
   int macro_mask, macro_tmask; \
 \
   (Sep) = 0; \
   Offset.y <<= 1; \
   Offset.z <<= 2; \
   macro_mask = 0x1; \
   for (macro_i=0; macro_i<LEVELS_MAX; macro_i++) { \
      (Sep) |= Offset.x & macro_mask; \
      Offset.x <<= 2; \
      macro_mask <<= 1; \
      (Sep) |= Offset.y & macro_mask; \
      Offset.y <<= 2; \
      macro_mask <<= 1; \
      (Sep) |= Offset.z & macro_mask; \
      Offset.z <<= 2; \
      macro_mask <<= 1; \
   } \
}


/****************************************************************
*
*  Sep2Vec - translation routine
*
*  this function takes a cell separation index and returns
*  the offsett to the new cell
*
*  the algorithm is based upon section 3.2 of elliotts disertation
*  with some enhancements.  an important note is that the level macro_mask
*  is stored as a negated value.
*/

/*
 * int sep - cell separation index
 * IntVector offset - remote cell id
*/

#define SEP2VEC(Sep,Offset) \
{ \
 \
   int macro_i; \
   int macro_sep; \
   int macro_mask, macro_tmask; \
 \
   macro_sep = (Sep); \
   (Offset).x = 0; \
   (Offset).y = 0; \
   (Offset).z = 0; \
   macro_mask = 0x01; \
   for ( macro_i=0; macro_i<LEVELS_MAX; macro_i++ ) { \
      (Offset).x |= macro_mask & macro_sep; \
      macro_sep >>= 1; \
      (Offset).y |= macro_mask & macro_sep; \
      macro_sep >>= 1; \
      (Offset).z |= macro_mask & macro_sep; \
      macro_mask <<= 1; \
   } \
   macro_tmask = (-1) << LEVELS_MAX; \
   macro_mask = 0x01 << (LEVELS_MAX-1); \
   if ( (Offset).x & macro_mask ) \
      (Offset).x |= macro_tmask; \
   if ( (Offset).y & macro_mask ) \
      (Offset).y |= macro_tmask; \
   if ( (Offset).z & macro_mask ) \
      (Offset).z |= macro_tmask; \
}

