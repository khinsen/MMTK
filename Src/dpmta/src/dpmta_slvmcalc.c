/* 
*  dpmta_slvmcalc.c - routines to perform the multipole interactions.
*
*  w. t. rankin
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*  upward pass code taken from dpmta_slvmpexp.c version 2.3.
*
*  downward pass code taken fom dpmta_slvcalc.c version 2.3, since it
*  doesn't make sense to have them in the same module as the particle
*  force interactions.
*
*/

static char rcsid[] = "$Id: dpmta_slvmcalc.c,v 2.25 1998/04/01 20:08:20 wrankin Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_slvmcalc.c,v $
 * Revision 2.25  1998/04/01 20:08:20  wrankin
 * added support for HILBERT ordering
 * general cleanup of code structure (more OOP if you will)
 *
 * Revision 2.24  1997/11/12 13:49:52  wrankin
 * updates to communications routines and general cleanup of code in prep
 *   for introducing load balancing functionality in hte near future
 *
 * Revision 2.23  1997/11/07 16:49:33  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 * Revision 2.22  1997/09/29 20:24:58  wrankin
 * fixed problem with invalid (empty) multipoles during upward pass.
 * cell indexing by processor was inconsistant between master/slave.
 *
 * Revision 2.21  1997/05/12 18:06:07  wrankin
 * added routines to clean up dynamic memory allocation when
 * PMTAinit() is called
 *
 * Revision 2.20  1997/05/09  21:04:30  wrankin
 * added procedures for cleanup of dynamic data structures
 *
 * Revision 2.19  1997/04/11  21:13:38  wrankin
 * fixed error in LJ multipole calculation
 *
 * Revision 2.18  1997/03/26  20:36:24  wrankin
 * rewrote Ilist and IIlist code to handle reallocation of array sizes
 *   that happen as a resuly of re-scaling a non-cubic unit cell.
 *
 * Revision 2.17  1997/03/04  19:18:10  wrankin
 * updates to timing codes
 *
 * Revision 2.16  1997/02/26  16:54:30  wrankin
 * added support for serial non-pvm compilation of DPMTA
 *
 * Revision 2.15  1997/01/13  22:05:20  wrankin
 * added code to perform parallel macroscopic expansion computation
 *
 * Revision 2.14  1996/11/18  19:29:35  wrankin
 * update to virial code to use multipole virial calculation
 * old virial code is turned on by defining OLDVIRIAL
 * direct test program now computes absolute virial tensor.
 *
 * Revision 2.13  1996/09/24  18:42:53  wrankin
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
 * Revision 2.12  1996/08/20  17:12:47  wrankin
 * Enhancement to allow any number of processes to be utilized.
 * Formally the number of precesses was limited to a power of 2.
 *
 * Revision 2.11  1996/08/09  15:30:57  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.10  1996/02/29  21:13:43  wrankin
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
 * Revision 2.9  1996/02/12  15:21:11  wrankin
 * Added Macroscopic Assemblies code.
 *
 * Revision 2.8  1995/12/08  23:00:31  wrankin
 * preliminary release of DPMTA 2.3
 *   - added working Periodic Boundary Conditions (PDC) flag
 *   - added code for Virial Computation (not working yet)
 *   - general cleanup of some modules
 *
 * Revision 2.7  1995/11/29  22:29:26  wrankin
 * addition of periodic boundary code
 * beginning work on virial code.
 *
 * Revision 2.6  1995/10/02  20:58:39  wrankin
 * changes to support creation and use of MPE x-fer matrix for
 *   LJ multipole calculation
 *
 * Revision 2.5  1995/10/01  21:46:10  wrankin
 * implemented LJ multipole calculations - compile time option only
 *
 * Revision 2.4  1995/07/10  02:47:16  wrankin
 * multipole processing code modified to use precomputed mpe transfer
 *   matrices.
 *
 * in addition, all multippole calculation routines have been removed
 *   from this source distribution and will be placed in a separate
 *   multipole library to be supplied externally.
 *
 * Revision 2.3  1995/07/01  03:27:00  wrankin
 * initial cut at precomputation of mpe transfer functions.
 * this works for vanilla M2L.  it does not work for FFT enhancements.
 *
 * Revision 2.2  1995/06/27  14:20:26  wrankin
 * Implementation of relative interaction list algorithm.
 * Interaction list generation code removed from master library.
 * Cell Table restructured as a 2-D table indexed by level/cell
 * All direct computations are now performed in a single pass.
 * General cleanup and removal of old files.
 *
 * Revision 2.1  1995/06/13  04:26:13  wrankin
 * No updates.
 * Sources placed under CVS.
 *
 * Revision 1.4  1995/04/09  22:21:06  wrankin
 * renaming of all global variables to prevent conflicts for
 * port of dpmta library to the T3D.
 *
 * Revision 1.3  1995/02/24  21:20:49  wrankin
 * fixed prototype definition for getnpids()
 * fixed bounds on upward pass so that multipoles for all cells
 *   including the level-0 cell are computed.
 *
 * Revision 1.2  1994/10/19  00:20:15  wrankin
 * added check for valid cell id before calling FFT_Multipole()
 * added check for valid multipole in FFT_Multipole()
 * added Calc_M2L_FFT() routine
 * fixed clearing the local accumulator by calling the
 *   correct routine in Clear_accum()
 *
 * Revision 1.1  1994/10/14  05:06:40  wrankin
 * Initial revision
 *
 *
*/

/* include files */
#include <stdio.h>
#include <math.h>
#include "dpmta_pvm.h"
#include "dpmta_cell.h"
#include "dpmta_slvglobals.h"


/* some required prototyping */

#include "dpmta_distmisc.h"
#include "dpmta_slvcomm.h"
#include "dpmta_slvmkhl.h"
#include "dpmta_slvmkil.h"
#ifdef MACROSCOPIC
#include "dpmta_slvmacro.h"
#endif
#ifdef TIME
#include "dpmta_timer.h"
#endif


/* prototyping internal functions */

void Calc_multipole_exp( int, int );
void Calc_M2M( int, int, int, int );
void Calc_M2L( int, int, int, int, Vector );
void Calc_M2L_FFT( int, int, int, int, Vector );
void Calc_L2L( int, int, int, int );
void Calc_Forces( int, int );
void Clear_mpole( int, int );
void Clear_local( int, int );
void Clear_accum();
void Create_accum();
void Delete_accum();
void FFT_Multipole( int, int );
void IFFT_Local(int, int );

void Calc_M2L_S( int, int, int, int, Mtype, MtypeLJ );
void Calc_M2L_FFT_S( int, int, int, int, Mtype, MtypeLJ );

#ifdef VIRIAL
void Calc_MCM( int, int, int, int, Mtype, Vector );
#endif

#ifdef MACROSCOPIC
void Calc_Macroscopic();
#endif


/****************************************************************
*
*  this procedure performs the one-time mpe initialization of
*  global constants and arrays
*
*/

void MultipoleSetup()
{

   Cinit(Dpmta_Mp);

#ifdef COMP_LJ
   LJinit(Dpmta_Mp_LJ);
#endif

   if (Dpmta_FFT) {
     CinitF(Dpmta_Mp, Dpmta_FftBlock);
     Create_accum();
   }

#ifdef MACROSCOPIC
   if (Dpmta_PBC) {
      MacroInit(Dpmta_K, Dpmta_Mp, Dpmta_FFT, Dpmta_Theta,
		Dpmta_Pid, Dpmta_Nproc, Dpmta_Tids, 0);
   }
#endif

} /* MultipoleSetup */


/****************************************************************
*
*  this procedure performs the cleanup of global structures created
*  during MultipoleSetup()
*
*/

void MultipoleCleanup()
{

   Ccleanup(Dpmta_Mp);

#ifdef COMP_LJ
   LJcleanup(Dpmta_Mp_LJ);
#endif

   if (Dpmta_FFT) {
      CcleanupF(Dpmta_Mp, Dpmta_FftBlock);
      Delete_accum();
   }

#ifdef MACROSCOPIC
   if (Dpmta_PBC) {
      MacroCleanup();
   }
#endif

} /* MultipoleCleanup */



/****************************************************************
*
*  Slave_Mpole_Exp() - this procedure will perform the upward
*  pass of the mulitpole expansion.
*
*  this procedure operates on the cell/tree level.  there are no
*  operations directly on the multipole data structures.  rather,
*  lower level procedures are called.
*
*/

void Slave_Mpole_Exp()
{
   int i,j,k;                  /* loop counters, misc indexes */
   int id;		       /* cell id */
   int level;                  /* cell tree level */
   int done;                   /* lop termination flag */
   
   /*
   *  compute the multipole expansions from the particles 
   *  in the leaf cells 
   */

   level = Dpmta_NumLevels - 1;
   if ( Dpmta_Sindex[level] >= 0 ) {
      for ( i=Dpmta_Sindex[level]; i<=Dpmta_Eindex[level]; i++ ) {
	 id = index2cell(i,level);
	 Clear_mpole(level,id);
	 Calc_multipole_exp(level,id);
      } /* for i */
   } /* if Dpmta_Sindex */
   

   /*
   *  Upward Pass 
   */

   done = FALSE;

   if ( Dpmta_Sindex[level] == -1 ) {
      done = TRUE;
   }

   while (done == FALSE) {

      /*
      *  clear out the multipoles for the parents cells of all
      *  the cells at the current level.  this is in preparation
      *  for doing the M2M translation in the next step
      *
      *  note that we want to clear the parent cells that we 
      *  do not own as well, so a simple traversal of the 
      *  Sindex and Eindex array for the parent level will not
      *  be sufficient.
      */

      k = -1;

      for (i=Dpmta_Sindex[level]; i<=Dpmta_Eindex[level]; i++) {
	 id = index2cell(i,level);
	 j = getparent(id);
	 if ( j != k ) {
	    Clear_mpole(level-1,j);
	    k = j;
	 } /* if (j!=k) */
      } /* for i */


      /*
      *  for each cell in the level, calculate its MPE effect on
      *  its parent.  note that even if we do not own the 
      *  parent of the current cell, a parent has still been
      *  allocated in the cell table.
      *
      *  the previous call to Clear_mpole() verified that the cell
      *  and multipole expansion for the parent actually exists
      */

      for (i=Dpmta_Sindex[level]; i<=Dpmta_Eindex[level]; i++) {
	 id = index2cell(i,level);
         j = getparent(id);
         Calc_M2M(level,id,level-1,j);
      } /* for i */

#ifndef SERIAL
      /*
      *  if we do not own parents of cells at our current level
      *  then distribute the multipoles to those pid who do
      *  own the parent cell(s).
      */

      Send_Mpe_to_Parent(level);
#endif

      /*
      *  pop up to the next level and check if we have any cells
      *  to process at this level.  if so, make sure to pick up any
      *  mpe's that were sent us.
      */

      level--;
      if ( Dpmta_Sindex[level] == -1 )
	 done = TRUE;
#ifndef SERIAL
      else
	 Recv_Mpe_from_Child(level);
#endif

      /*
      *  finally, if we are have gotten to the top level,
      *  then we are through by definition
      */

      if ( level == 0 )
	 done = TRUE;

   } /* while not done */

#ifdef MACROSCOPIC
   /*
    * if we are computing the macroscopic assemblies, then
    * we need to compute the macro expansion for the unit cell.
    * we need to do this before we FFT all the multipoles.
    *
    * the local expansion for the unit cell will be updated by
    * this.
    *
    */

   if ( Dpmta_Sindex[0] != -1 ) {
      if (Dpmta_PBC) {
         Calc_Macroscopic();
      }
   }
#endif

   /*
   *  if we are performing the FFT optimiztions, then at this
   *  point we need to cycle through all the cells that we own
   *  and translate the multipole expansion to the FFT domain.
   *
   */

   if (Dpmta_FFT) {
      for (i=0; i<Dpmta_NumLevels; i++) {
         for (j=Dpmta_Sindex[i]; j<=Dpmta_Eindex[i]; j++) {
            if ( j != -1 ) {
	       id = index2cell(j,i);
               FFT_Multipole(i,id);
	    } /* if j */
         } /* for j */
      } /* for i */
   } /* if Dpmta_FFT */

#ifdef DEBUG2
   /* dump the level 1 multipoles */
   if ( Dpmta_Sindex[2] != -1 ) {
      char filename[80];
      int il;
      for ( il=Dpmta_Sindex[2]; il<=Dpmta_Eindex[2]; il++ ) {
	 id = index2cell(il,2);
	 sprintf(filename,"/tmp/MPD-m%d-p%d-n%d.data",
		 id,Dpmta_Pid,Dpmta_Nproc);
	 MDumpRaw_C(Dpmta_CellTbl[2][id]->m, Dpmta_Mp, filename);
      }
   }
#endif

} /* Slave_Mpole_Exp() */



/****************************************************************
*
*  Slave_MPE_Calc() - this procedure will perform the downward pass
*  of the multipole calculations.
*
*/

void Slave_MPE_Calc()
{

   int i,j,k,l;                  /* loop counters */
   int id;			 /* temp cell id */
   int pcell,plevel;             /* parent cell id and level */
   int rcell,rlevel;             /* remote cell id and level */
   int posn;                     /* cell position within parent */
   int sep;                      /* cell separation index */


   /*
   *  Downward Pass 
   */

#ifdef MACROSCOPIC
   /*
    * if we are computing the macroscopic assemblies, then
    * the local expansion for the unit cell has already been
    * calculated during the upward pass.
    *
    * if there are more than 1 slave processes, then we need to send the
    * expansion to the other processors.
    */

   if ( Dpmta_Sindex[0] != -1 ) {
      if (Dpmta_PBC) {
#ifndef SERIAL
         if (Dpmta_Nproc > 1) {
            Send_Lcl_to_Child(0);
	 }
#endif
      }
      else {
	 Clear_local(0,0);
      }
   }
#endif

   for (i=Dpmta_DownPassStart; i<Dpmta_NumLevels; i++) {

      /* does pid own a cell at this level? */
      if (Dpmta_Sindex[i] != -1) {

#ifndef SERIAL
         /*
         *  if the cell parent is owned by a different slave,
         *  then we need to recieve it.
         */

         Recv_Lcl_from_Parent(i);

#endif
	 /*
	 *  compute the mpe transfer matrix
	 */

         Compute_Hlist(i);

         /*
         *  cycle through each slave in this level and compute the
         *  local multipole expansion
         */

	 for (j=Dpmta_Sindex[i]; j<=Dpmta_Eindex[i]; j++) {
	    id = index2cell(j,i);

            /*
            *  probably want to zero the local expansion here
            */

            Clear_local(i,id);

            /*
            *  if this is not the top level of the expansion, then
            *  shift the local expansion for the parent to the center
            *  of the child cell.
	    *
	    *  if we are doing macroscopic assemblies, then the unit
	    *  cell has a valid local expansion.
            */

#ifdef MACROSCOPIC
            k = getparent(id);
            if ( Dpmta_CellTbl[i-1][k] == NULL ) {
               fprintf(stderr,"ERROR: cell %d/%d not allocated\n",i-1,k);
               exit(-1);
	    }
            Calc_L2L(i-1,k,i,id);

#else
	    if ( i != Dpmta_DownPassStart ) {
	       k = getparent(id);
               if ( Dpmta_CellTbl[i-1][k] == NULL ) {
                  fprintf(stderr,"ERROR: cell %d/%d not allocated\n",i-1,k);
                  exit(-1);
	       }
	       Calc_L2L(i-1,k,i,id);
	    } /* if i */
#endif            


            /*
            *  traverse the cells interaction list and transpose all
            *  the multipole expansion from the near cells to the
            *  local expansion of the current cell.
            *
            *  if we are using the FFT optimizations, then we need to
            *  accumulate the FFT MPE's into a buffer and then IFFT
            *  the entire array in one swell foop.
            */

            if ( Dpmta_FFT ) {

               /* zero out the local accumulator */
               Clear_accum();

	       /* determin cell position within parent cell */
	       posn = id & 0x07;

               /* cycle through the parent level interaction list */
	       plevel = i - 1;
	       pcell = getparent(id);
  	       rlevel = plevel;

               for ( k=0; k<Dpmta_Intlist[posn].pcnt; k++ ) {
                  sep = Dpmta_Intlist[posn].plist[k];
                  CELL2CELL(plevel,pcell,sep,rcell,l);
 		  if ( l == 0 ) {
                     if ( Dpmta_CellTbl[rlevel][rcell] == NULL ) {
                        fprintf(stderr,"ERROR: cell %d/%d not allocated\n",
				rlevel, rcell);
                        exit(-1);
                     }

                     /* accumulate the local expansion */
#ifdef COMP_LJ
		     Calc_M2L_FFT_S(rlevel,rcell,i,id,
				    Dpmta_Hlist[posn].plist[k],
				    Dpmta_Hlist[posn].plist_lj[k]);
#else
		     Calc_M2L_FFT_S(rlevel,rcell,i,id,
				    Dpmta_Hlist[posn].plist[k],
				    (MtypeLJ)NULL );
#endif

		  } /* if Cell2Cell() */
               } /* for k */


               /* cycle through the sibling level interaction list */
  	       rlevel = i;

               for ( k=0; k<Dpmta_Intlist[posn].scnt; k++ ) {
                  sep = Dpmta_Intlist[posn].slist[k];
                  CELL2CELL(i,id,sep,rcell,l);
 		  if ( l == 0 ) {
                     if ( Dpmta_CellTbl[rlevel][rcell] == NULL ) {
                        fprintf(stderr,"ERROR: cell %d/%d not allocated\n",
				rlevel, rcell);
                        exit(-1);
                     }

                     /* accumulate the local expansion */
#ifdef COMP_LJ
		     Calc_M2L_FFT_S(rlevel,rcell,i,id,
				    Dpmta_Hlist[posn].slist[k],
				    Dpmta_Hlist[posn].slist_lj[k]);
#else
		     Calc_M2L_FFT_S(rlevel,rcell,i,id,
				    Dpmta_Hlist[posn].slist[k],
				    (MtypeLJ)NULL );
#endif
		  } /* if Cell2Cell() */
               } /* for k */

               /*
               *  invert the local accum and place in the local
               *  expansion for the cell
               */

               IFFT_Local(i,id);

            } /* FFT */

            /* not FFT */
            else {

	       /* determine cell position within parent cell */
	       posn = id & 0x07;

               /* cycle through the parent level interaction list */
	       plevel = i - 1;
	       pcell = getparent(id);
  	       rlevel = plevel;

               for ( k=0; k<Dpmta_Intlist[posn].pcnt; k++ ) {
                  sep = Dpmta_Intlist[posn].plist[k];
                  CELL2CELL(plevel,pcell,sep,rcell,l);
 		  if ( l == 0 ) {
                     if ( Dpmta_CellTbl[rlevel][rcell] == NULL ) {
                        fprintf(stderr,"ERROR: cell %d/%d not allocated\n",
				rlevel, rcell);
                        exit(-1);
                     }
#ifdef COMP_LJ
                     Calc_M2L_S(rlevel,rcell,i,id,
				Dpmta_Hlist[posn].plist[k],
				Dpmta_Hlist[posn].plist_lj[k]);
#else
                     Calc_M2L_S(rlevel,rcell,i,id,
				Dpmta_Hlist[posn].plist[k],
				(MtypeLJ)NULL );
#endif

#ifdef VIRIAL
                     Calc_MCM(rlevel,rcell,i,id,
			      Dpmta_Hlist[posn].plist[k],
			      Dpmta_Hlist[posn].plist_vec[k]);
#endif

		  } /* if l */
               } /* for k */


               /* cycle through the sibling level interaction list */
  	       rlevel = i;

               for ( k=0; k<Dpmta_Intlist[posn].scnt; k++ ) {
                  sep = Dpmta_Intlist[posn].slist[k];
                  CELL2CELL(i,id,sep,rcell,l);
 		  if ( l == 0 ) {
                     if ( Dpmta_CellTbl[rlevel][rcell] == NULL ) {
                        fprintf(stderr,"ERROR: cell %d/%d not allocated\n",
				rlevel, rcell);
                        exit(-1);
                     }
 
#ifdef COMP_LJ
                     Calc_M2L_S(rlevel,rcell,i,id,
				Dpmta_Hlist[posn].slist[k],
				Dpmta_Hlist[posn].slist_lj[k]);
#else
                     Calc_M2L_S(rlevel,rcell,i,id,
				Dpmta_Hlist[posn].slist[k],
				(MtypeLJ)NULL );
#endif

#ifdef VIRIAL
                     Calc_MCM(rlevel,rcell,i,id,
			      Dpmta_Hlist[posn].slist[k],
			      Dpmta_Hlist[posn].slist_vec[k]);
#endif

		  } /* if l */
               } /* for k */

            } /* not Dpmta_FFT */

         } /* for id */

#ifndef SERIAL
         /*
	 *  check if any children of the cells on the current
	 *  are owned by another pid and if so, send the local
	 *  expansions to that pid.
         */

         Send_Lcl_to_Child(i);
#endif

      } /* if dpmta_scell */
   } /* for i */

} /* end Slave_MPE_Calc */



/*****************************************************************
*
*  we have completed to downward pass to compute the multipole 
*  expansions.  based on this, the forces on the individual particles
*  can be computed from the mpe data
*
*/

void Slave_MPE_Force()
{
   int i;
   int id;
   int level;

   level = Dpmta_NumLevels - 1;

   for ( i=Dpmta_Sindex[level]; i<=Dpmta_Eindex[level]; i++)  {
      id = index2cell(i,level);
      Calc_Forces(level,id);
   }

} /* Slave_MPE_Force() */




/****************************************************************/


/****************************************************************
*
*  all the following are called from within this module.
*
*  they traverse the data structures within the cells and
*  call the appropriate multipole library routines
*
*/

/*
* global multipole acumulator used for FFT processing
*/

static Mtype LocalAccum;


/****************************************************************
*
*  Calc_multipole_exp() - calculate the multipole expansion for all
*     particles within a cell
*
*  w.t.rankin
*
*/

void Calc_multipole_exp( int level, int cell )
{

   int i;                      /* loop counters */
   int num_parts;              /* number of particles in cell */
   ParticlePtr cell_plist;     /* list of particles */
   Mtype mpole;                /* cell multipole expansion array */
   Real q;                     /* particle charge */
   Vector v;                   /* vector from cell center to particle */

#ifdef COMP_LJ
   Real a;                     /* LJ alpha parameter */
   MtypeLJ mpole_lj;           /* cell LJ multipole expansion array */
#endif

   /* 
   *  loop through all particles within the given cell 
   */


   num_parts = Dpmta_CellTbl[level][cell]->n;

   if ( num_parts != 0 ) {

      Dpmta_CellTbl[level][cell]->mvalid = TRUE;

      cell_plist = Dpmta_CellTbl[level][cell]->plist;
      mpole = Dpmta_CellTbl[level][cell]->m;
#ifdef COMP_LJ
      mpole_lj = Dpmta_CellTbl[level][cell]->m_lj;
#endif

      for ( i=0; i<num_parts; i++ ) {

         /* 
         *  add the partilce to the multipole expansion
         */

         v.x = cell_plist[i].p.x;
         v.y = cell_plist[i].p.y;
         v.z = cell_plist[i].p.z;
         q = cell_plist[i].q;
         AddMultipoleC(mpole, Dpmta_Mp, q, v );

#ifdef COMP_LJ
         a = cell_plist[i].a;
         AddMultipoleLJ(mpole_lj, Dpmta_Mp_LJ, a, v );
#endif

      } /* for i */
   } /* if num_parts */

 
  else {
      Dpmta_CellTbl[level][cell]->mvalid = FALSE;
  }

} /* Calc_multipole_exp */


/****************************************************************
*
*   FFT_Multipole - convert multipole expansion in cell to fourier
*   domain
*
*/

void FFT_Multipole( int level, int cell )
{
   Mtype mpole;

   if ( Dpmta_CellTbl[level][cell]->mvalid == TRUE ) {
      mpole = Dpmta_CellTbl[level][cell]->m;
      Warp_M2L(mpole, mpole, Dpmta_Mp, Dpmta_FftBlock);
   } /* if mvalid */

} /* FFT_Multipole() */


/****************************************************************
*
*  Calc_M2M() - shift the multipole expansion for a child cell and add
*    it to the multipole expansion of the parent cell.
*
*  w.t.rankin
*
*/

void Calc_M2M( int clevel, int ccell, int plevel, int pcell )
{

   Mtype cmpole;      /* child multipole expansion */
   Mtype pmpole;      /* parent multipole expansion */
   Vector v;          /* position vector from parent to child center */

#ifdef COMP_LJ
   MtypeLJ cmpole_lj; /* child LJ multipole expansion */
   MtypeLJ pmpole_lj; /* parent LJ multipole expansion */
#endif


   /*
   *  check to see if the child cell has a valid multipole expansion
   */

   if ( Dpmta_CellTbl[clevel][ccell]->mvalid == FALSE )
      return;


   /*
   *  calculate the cartesian vector from the center of the child cell
   *  to the center of the parent cell
   */

   v.x = Dpmta_CellTbl[plevel][pcell]->p.x - Dpmta_CellTbl[clevel][ccell]->p.x;
   v.y = Dpmta_CellTbl[plevel][pcell]->p.y - Dpmta_CellTbl[clevel][ccell]->p.y;
   v.z = Dpmta_CellTbl[plevel][pcell]->p.z - Dpmta_CellTbl[clevel][ccell]->p.z;

   /*
   *  set the parent multipole expansion to valid
   */

   Dpmta_CellTbl[plevel][pcell]->mvalid = TRUE;

   /*
   *  shift the child multipole expansion, placing the result in 
   *  the parent multipole expansion
   */

   cmpole = Dpmta_CellTbl[clevel][ccell]->m;
   pmpole = Dpmta_CellTbl[plevel][pcell]->m;
   M2M_C(cmpole, pmpole, Dpmta_Mp, v);

#ifdef COMP_LJ
   cmpole_lj = Dpmta_CellTbl[clevel][ccell]->m_lj;
   pmpole_lj = Dpmta_CellTbl[plevel][pcell]->m_lj;
   M2M_LJ(cmpole_lj, pmpole_lj, Dpmta_Mp_LJ, v);
#endif

} 


/****************************************************************
*
*  Calc_L2L() - take the local expansion for a parent cell and
*     add it to the local expansion of a child cell.
*
*  w.t.rankin
*
*/

void Calc_L2L( int plevel, int pcell, int clevel, int ccell )
{

   Mtype plocal;      /* parent local expansion */
   Mtype clocal;      /* child local expansion */
   Vector v;          /* position vector from parent to child center */

#ifdef COMP_LJ
   MtypeLJ plocal_lj; /* parent local LJ expansion */
   MtypeLJ clocal_lj; /* child local LJ expansion */
#endif

   /*
   *  check to see if the parent cell has a valid local expansion
   */
   
   if ( Dpmta_CellTbl[plevel][pcell]->mdata->lvalid == TRUE ) {
    

      /*
      *  set the child local expansion to valid
      */

      Dpmta_CellTbl[clevel][ccell]->mdata->lvalid = TRUE;

      /*
      *  calculate the cartesian vector from the center of the parent cell
      *  to the center of the child cell
      */

      v.x = Dpmta_CellTbl[clevel][ccell]->p.x -
	 Dpmta_CellTbl[plevel][pcell]->p.x;
      v.y = Dpmta_CellTbl[clevel][ccell]->p.y -
	 Dpmta_CellTbl[plevel][pcell]->p.y;
      v.z = Dpmta_CellTbl[clevel][ccell]->p.z -
	 Dpmta_CellTbl[plevel][pcell]->p.z;

      /*
      *  shift the remote multipole expansion, placing the result in 
      *  the local multipole expansion
      */

      plocal = Dpmta_CellTbl[plevel][pcell]->mdata->l;
      clocal = Dpmta_CellTbl[clevel][ccell]->mdata->l;
      L2L_C(plocal, clocal, Dpmta_Mp, v);

#ifdef COMP_LJ
      plocal_lj = Dpmta_CellTbl[plevel][pcell]->mdata->l_lj;
      clocal_lj = Dpmta_CellTbl[clevel][ccell]->mdata->l_lj;
      L2L_LJ(plocal_lj, clocal_lj, Dpmta_Mp_LJ, v);
#endif
   } /* if lvalid */

} /* Calc_L2L() */


/****************************************************************
*
*  Calc_M2L() - shift the multipole expansion for a remote cell and add
*    it to the local expansion of the local cell.
*
*  w.t.rankin
*
*/

void Calc_M2L(
   int rlevel,        /* remote level */
   int rcell,         /* remote cell */
   int llevel,        /* local level */
   int lcell,         /* local cell */
   Vector sep )       /* separation vector */
{

   Mtype rmpole;      /* remote multipole expansion */
   Mtype llocal;      /* local cells local expansion */
   Vector v;          /* position vector from parent to child center */

#ifdef COMP_LJ
   MtypeLJ rmpole_lj; /* remote multipole LJ expansion */
   MtypeLJ llocal_lj; /* local cells local LJ expansion */
#endif

   /*
   *  check to see if the remote cell has a valid multipole expansion
   */

   if ( Dpmta_CellTbl[rlevel][rcell]->mvalid == TRUE ) {


      /*
      *  set the local multipole expansion to valid
      */

      Dpmta_CellTbl[llevel][lcell]->mdata->lvalid = TRUE;

      /*
      *  calculate the cartesian vector from the center of the remote cell
      *  to the center of the local cell
      */

      v.x = -sep.x;
      v.y = -sep.y;
      v.z = -sep.z;

      /*
      *  shift the remote multipole expansion, placing the result in 
      *  the local cells local expansion
      */

      rmpole = Dpmta_CellTbl[rlevel][rcell]->m;
      llocal = Dpmta_CellTbl[llevel][lcell]->mdata->l;
      M2L_C(rmpole, llocal, Dpmta_Mp, v);

#ifdef COMP_LJ
      rmpole_lj = Dpmta_CellTbl[rlevel][rcell]->m_lj;
      llocal_lj = Dpmta_CellTbl[llevel][lcell]->mdata->l_lj;
      M2L_LJ(rmpole_lj, llocal_lj, Dpmta_Mp_LJ, v);
#endif

   } /* if lvalid */
   
} /* Calc_M2L() */


/****************************************************************
*
*  Calc_M2L_S() - shift the multipole expansion for a remote cell
*    and add it to the local expansion of the local cell.
*
*    this function uses the mpe transfer matrix to perform this.
**
*/

void Calc_M2L_S(
   int rlevel,       /* remote level */
   int rcell,        /* remote cell */
   int llevel,       /* local level*/
   int lcell,        /* local cell */
   Mtype h,          /* mpe transfer matrix */
   MtypeLJ h_lj )    /* LJ mpe transfer matrix (if used)*/
{

   Mtype rmpole;      /* remote multipole expansion */
   Mtype llocal;      /* local cells local expansion */

#ifdef COMP_LJ
   MtypeLJ rmpole_lj; /* remote LJ multipole */
   MtypeLJ llocal_lj; /* local cells local LJ expansion */
#endif

   /*
   *  check to see if the remote cell has a valid multipole expansion
   */

   if ( Dpmta_CellTbl[rlevel][rcell]->mvalid == TRUE ) {

      /*
      *  set the local multipole expansion to valid
      */

      Dpmta_CellTbl[llevel][lcell]->mdata->lvalid = TRUE;

      /*
      *  shift the remote multipole expansion, placing the result in 
      *  the local cells local expansion
      */

      rmpole = Dpmta_CellTbl[rlevel][rcell]->m;
      llocal = Dpmta_CellTbl[llevel][lcell]->mdata->l;
      M2L_Cshort(rmpole, llocal, h, Dpmta_Mp);

#ifdef COMP_LJ
      rmpole_lj = Dpmta_CellTbl[rlevel][rcell]->m_lj;
      llocal_lj = Dpmta_CellTbl[llevel][lcell]->mdata->l_lj;
      M2L_LJshort(rmpole_lj, llocal_lj, h_lj, Dpmta_Mp_LJ);
#endif

   } /* if mvalid */

} /* Calc_M2L_S() */



/****************************************************************
*
*  Calc_M2L_FFT() - shift the multipole expansion for a remote cell
*    to the center of the local cell and add it to the local 
*    accumulator
*
*/


void Calc_M2L_FFT(
   int rlevel,       /* remote level */
   int rcell,        /* remote cell */
   int llevel,       /* local level*/
   int lcell,        /* local cell */
   Vector sep )      /* separation vector */

{

   Mtype rmpole;      /* remote multipole expansion */
   Vector v;          /* position vector from parent to child center */

#ifdef COMP_LJ
   MtypeLJ rmpole_lj; /* remote multipole LJ expansion */
   MtypeLJ llocal_lj; /* local cells local LJ expansion */
#endif

   /*
   *  check to see if the remote cell has a valid multipole expansion
   */
   
   if ( Dpmta_CellTbl[rlevel][rcell]->mvalid == TRUE ) {

      /*
      *  set the local multipole expansion to valid
      */

      Dpmta_CellTbl[llevel][lcell]->mdata->lvalid = TRUE;

      /*
      *  calculate the cartesian vector from the center of the remote cell
      *  to the center of the local cell
      */

      v.x = -sep.x;
      v.y = -sep.y;
      v.z = -sep.z;

      /*
      *  shift the remote multipole expansion, placing the result in 
      *  the local cells local expansion
      */

      rmpole = Dpmta_CellTbl[rlevel][rcell]->m;
      M2L_C_F(rmpole, LocalAccum, Dpmta_Mp, Dpmta_FftBlock, v);

#ifdef COMP_LJ
      rmpole_lj = Dpmta_CellTbl[rlevel][rcell]->m_lj;
      llocal_lj = Dpmta_CellTbl[llevel][lcell]->mdata->l_lj;
      M2L_LJ(rmpole_lj, llocal_lj, Dpmta_Mp_LJ, v);
#endif

   } /* if mvalid */
   
} /* Calc_M2L_FFT() */


/****************************************************************
*
*  Calc_M2L_FFT_S() - shift the multipole expansion for a remote cell
*    to the center of the local cell and add it to the local 
*    accumulator
*
*    this function uses the mpe transfer matrix to perform this.
*
*/

void Calc_M2L_FFT_S(
   int rlevel,       /* remote level */
   int rcell,        /* remote cell */
   int llevel,       /* local level*/
   int lcell,        /* local cell */
   Mtype h,          /* mpe transfer matrix */
   MtypeLJ h_lj )    /* LJ mpe transfer matrix (if used)*/
{

   Mtype rmpole;      /* remote multipole expansion */

#ifdef COMP_LJ
   MtypeLJ rmpole_lj; /* remote LJ multipole expansion */
   MtypeLJ llocal_lj; /* local cells local LJ expansion */
#endif

   /*
   *  check to see if the remote cell has a valid multipole expansion
   */

   if ( Dpmta_CellTbl[rlevel][rcell]->mvalid == TRUE ) {

      /*
      *  set the local multipole expansion to valid
      */

      Dpmta_CellTbl[llevel][lcell]->mdata->lvalid = TRUE;

      rmpole = Dpmta_CellTbl[rlevel][rcell]->m;

      /*
      *  shift the remote multipole expansion, placing the result in 
      *  the local cells local expansion
      */

      M2L_C_Fshort(rmpole, LocalAccum, h, Dpmta_Mp, Dpmta_FftBlock);

#ifdef COMP_LJ
      rmpole_lj = Dpmta_CellTbl[rlevel][rcell]->m_lj;
      llocal_lj = Dpmta_CellTbl[llevel][lcell]->mdata->l_lj;
      M2L_LJshort(rmpole_lj, llocal_lj, h_lj, Dpmta_Mp_LJ);
#endif

   } /* if mvalid */
   
} /* Calc_M2L_FFT_S() */


/****************************************************************
*
*  IFFT_Local() - inverse FFT the local accumulator and place the
*    result in the local expansion of the local cell.
*
*/

void IFFT_Local(
   int llevel,           /* local level */
   int lcell             /* local cell */
   )
{

   Mtype llocal;      /* local cells local expansion */

   if ( Dpmta_CellTbl[llevel][lcell]->mdata->lvalid == TRUE ) {

      llocal = Dpmta_CellTbl[llevel][lcell]->mdata->l;
      Unwarp_M2L(LocalAccum, llocal, Dpmta_Mp, Dpmta_FftBlock);

   } /* if lvalid */

} /* IFFT_Local */



/****************************************************************
*
*  Calc_Forces() - take the local expansion for a cell and
*     compute the force and potentials for each particle in the
*     cell.
*
*     this routine will also calclate the virial energies, if needed
*
*/

void Calc_Forces(
   int level,       /* cell level */
   int cell         /* cell id */
   )

{

   int i;                      /* loop counters */
   int num_parts;              /* number of particles in cell */
   ParticlePtr cell_plist;     /* list of particles */
   PartInfoPtr cell_flist;     /* list of force result */
   Real q;                     /* particle charge */
   Mtype clocal;               /* cell local expansion */
   Vector v;                   /* position vector from part to center */
   Vector f;                   /* resulting force vector */
   Real pot;                   /* resulting part potential */

#ifdef COMP_LJ
   PartInfoPtr cell_f_lj;      /* list of force result */
   MtypeLJ clocal_lj;          /* cell local LJ expansion */
   Real a;                     /* particle LJ alpha parameter */
#endif

#ifdef OLDVIRIAL
   Real vmag;
   Real fx2, fy2, fz2;
#endif

   /*
   *  check to see if the cell has a valid local expansion
   */
   
   if ( Dpmta_CellTbl[level][cell]->mdata->lvalid == TRUE ) {

   /* 
   *  loop through all particles within the given cell 
   */

      num_parts = Dpmta_CellTbl[level][cell]->n;
      if ( num_parts != 0 ) {

	 cell_plist = Dpmta_CellTbl[level][cell]->plist;
	 cell_flist = Dpmta_CellTbl[level][cell]->mdata->flist;
	 clocal = Dpmta_CellTbl[level][cell]->mdata->l;

#ifdef COMP_LJ
	 cell_f_lj = Dpmta_CellTbl[level][cell]->mdata->f_lj;
	 clocal_lj = Dpmta_CellTbl[level][cell]->mdata->l_lj;
#endif

	 for ( i=0; i<num_parts; i++ ) {

	    /*
	    *  compute the distance from the center point to the particle
	    */

	    v.x = cell_plist[i].p.x;
	    v.y = cell_plist[i].p.y;
	    v.z = cell_plist[i].p.z;

	    /* 
	    *  compute the force on the particle
	    */

	    q = cell_plist[i].q;
	    Force_C(clocal, Dpmta_Mp, q, v, &pot, &f );

	    /* put the forces into the results array */

	    cell_flist[i].f.x += f.x;
	    cell_flist[i].f.y += f.y;
	    cell_flist[i].f.z += f.z;
	    cell_flist[i].v += pot;

#ifdef OLDVIRIAL
	    Dpmta_Vpot += pot * 0.5;
	    fx2 = f.x * f.x;
	    fy2 = f.y * f.y;
	    fz2 = f.z * f.z;
	    vmag = fx2 + fy2 + fz2;
	    if (vmag != 0.0) {
	       vmag = 0.5/vmag;
	       Dpmta_Vf.x -= pot * vmag * fx2;
	       Dpmta_Vf.y -= pot * vmag * fy2;
	       Dpmta_Vf.z -= pot * vmag * fz2;
	    } /* if */
#endif

#ifdef COMP_LJ
	    a = cell_plist[i].a;
	    Force_LJ(clocal_lj, Dpmta_Mp_LJ, a, v, &pot, &f );

	    /* put the forces into the results array */

	    cell_f_lj[i].f.x -= f.x;
	    cell_f_lj[i].f.y -= f.y;
	    cell_f_lj[i].f.z -= f.z;
	    cell_f_lj[i].v -= pot;

#ifdef OLDVIRIAL
	    Dpmta_Vpot_LJ += pot * 0.5;
	    fx2 = f.x * f.x;
	    fy2 = f.y * f.y;
	    fz2 = f.z * f.z;
	    vmag = fx2 + fy2 + fz2;
	    if (vmag != 0.0) {
	       vmag = 0.5/vmag;
	       Dpmta_Vf_LJ.x -= pot * vmag * fx2;
	       Dpmta_Vf_LJ.y -= pot * vmag * fy2;
	       Dpmta_Vf_LJ.z -= pot * vmag * fz2;
	    } /* if */

#endif

#endif

	 } /* for i */
      } /* if num_parts */

   } /* if lvalid */

} /* Calc_Forces() */


#ifdef MACROSCOPIC
/****************************************************************
*
*  Calc_Macroscopic() - take the multipole expansion for the
*    unit cell and compute the resulting local expansion representing
*    the macroscopic expansion of the unit cell.
*
*/

void Calc_Macroscopic()
{
   Mtype min;
   Mtype lout;

   min = Dpmta_CellTbl[0][0]->m;
   lout = Dpmta_CellTbl[0][0]->mdata->l;
   Dpmta_CellTbl[0][0]->mdata->lvalid = TRUE;

#ifdef TIME
   times(&startbuf);
#endif

   MacroCompute(min,lout);

#ifdef TIME
   times(&endbuf);
   times_arr[14] = (double)(endbuf.tms_utime - startbuf.tms_utime) /
      (double)CLK_TCK;
#endif

} /* Calc_Macroscopic() */
#endif


/****************************************************************
*
* misc routines to handle multipoles
*
*/

void Clear_mpole( int level, int cell )
{

   Mtype mpole;
#ifdef COMP_LJ
   MtypeLJ mpole_lj;
#endif

   
   /* 
   *  zero the expansion 
   */
   if ( Dpmta_CellTbl[level][cell] == NULL ) {
      fprintf(stderr,"ERROR: Cell %d/%d not allocated\n",
	      level, cell);
      exit(-1);
   }

   mpole = Dpmta_CellTbl[level][cell]->m;
   if (mpole == NULL) {
      fprintf(stderr,"ERROR: Multipole at cell %d/%d not allocated\n",
	      level, cell);
      exit(-1);
   }

   if (Dpmta_FFT) 
      CMclearF(mpole,Dpmta_Mp);
   else
      CMclear(mpole,Dpmta_Mp);

#ifdef COMP_LJ
   mpole_lj = Dpmta_CellTbl[level][cell]->m_lj;
   if (mpole_lj == NULL) {
      fprintf(stderr,"ERROR: LJ Multipole at cell %d/%d not allocated\n",
	      level, cell);
      exit(-1);
   }

   LJMclear(mpole_lj,Dpmta_Mp_LJ);
#endif

   /* set valid flag to false */

   Dpmta_CellTbl[level][cell]->mvalid = FALSE;
   
} /* Clear_mpole */


/****************************************************************
 *
 *  clear out the local multipole expansions for a given cell.
 *
 */

void Clear_local( int level, int cell )
{
   Mtype local;
#ifdef COMP_LJ
   MtypeLJ local_lj;
#endif

   
   /* 
   *  zero the expansion 
   */

   local = Dpmta_CellTbl[level][cell]->mdata->l;
   if (local == NULL) {
      fprintf(stderr,"ERROR: LocalExp at cell %d not allocated\n", cell);
      exit(-1);
   }

   CMclear(local,Dpmta_Mp);

#ifdef COMP_LJ
   local_lj = Dpmta_CellTbl[level][cell]->mdata->l_lj;
   if (local_lj == NULL) {
      fprintf(stderr,"ERROR: LJ LocalExp at cell %d not allocated\n", cell);
      exit(-1);
   }

   LJMclear(local_lj,Dpmta_Mp_LJ);
#endif

   /* set valid flag to false */

   Dpmta_CellTbl[level][cell]->mdata->lvalid = FALSE;

} /* Clear_local */


/****************************************************************
 *
 *
 */

void Create_accum()
{
   CallocFrev(&LocalAccum,Dpmta_Mp,Dpmta_FftBlock);
}


/****************************************************************
 *
 *
 */

void Clear_accum()
{
   /* 
   *  zero the expansion 
   */
   CMclearFrev(LocalAccum,Dpmta_Mp,Dpmta_FftBlock);

}

/****************************************************************
 *
 *
 */

void Delete_accum()
{
   CfreeFrev(LocalAccum,Dpmta_Mp,Dpmta_FftBlock);
}




#ifdef VIRIAL

/****************************************************************
*
*  Calc_MCM() - compute the virial tensor for the multpoles of 
*    two well separated regions.
*
*  w.t.rankin
*
*/

void Calc_MCM(
   int rlevel,        /* remote level */
   int rcell,         /* remote cell */
   int llevel,        /* local level */
   int lcell,         /* local cell */
   Mtype h,           /* mpe transfer matrix */
   Vector sep )       /* separation vector */
{

   Mtype rmpole;      /* remote multipole expansion */
   Mtype lmpole;      /* local cells local expansion */
   Real pot;          /* potential sum*/
   Vector f;          /* force vector from parent to child center */
   Real req;          /* equivolent radius factor */


   /*
   *  check to see if the remote cell has a valid multipole expansion
   */

   if ( Dpmta_CellTbl[rlevel][rcell]->mvalid == FALSE )
      return;
   if ( Dpmta_CellTbl[llevel][lcell]->mvalid == FALSE )
      return;

   /*
   *  convolve the two multipoles together, leaving the result in the
   *  accumulator
   */
   CMclear(Dpmta_Temp_Mpe, Dpmta_Mp);
   rmpole = Dpmta_CellTbl[rlevel][rcell]->m;
   lmpole = Dpmta_CellTbl[llevel][lcell]->m;
   MCM_C(rmpole, lmpole, Dpmta_Temp_Mpe, Dpmta_Mp);

   ForceM_C(Dpmta_Temp_Mpe, Dpmta_Mp, 1.0, sep, &pot, &f);

   req = (f.x*f.x+f.y*f.y+f.z*f.z);

   if ( req != 0.0 ) {
      req = pot / req;
      Dpmta_Vpot += pot * 0.5;
      Dpmta_Vf.x -= 0.5 * req * f.x * f.x;
      Dpmta_Vf.y -= 0.5 * req * f.y * f.y;
      Dpmta_Vf.z -= 0.5 * req * f.z * f.z;
   }

} /* Calc_M2L() */

#endif
