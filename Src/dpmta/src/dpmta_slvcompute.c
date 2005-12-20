/* 
*  dpmta_slvcompute.c - routine to perform actuall DPMTA processing
*
*  w. t. rankin, w. elliott
*
*  Copyright (c) 1995 Duke University
*  All rights reserved
*
*  this routine implements the actaul DPMTA processing.
*  it was pulled out of the original dpmta_slave.c main() routine
*  for the purpose of t3d integration.
*
*/

static char rcsid[] = "$Id: dpmta_slvcompute.c,v 2.14 1998/04/01 20:08:16 wrankin Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_slvcompute.c,v $
 * Revision 2.14  1998/04/01 20:08:16  wrankin
 * added support for HILBERT ordering
 * general cleanup of code structure (more OOP if you will)
 *
 * Revision 2.13  1997/11/07 16:49:24  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 * Revision 2.12  1997/03/04 19:18:07  wrankin
 * updates to timing codes
 *
 * Revision 2.11  1997/02/26  20:43:31  wrankin
 * updated timing measurements and placed routines in single module
 *
 * Revision 2.10  1997/02/26  16:54:28  wrankin
 * added support for serial non-pvm compilation of DPMTA
 *
 * Revision 2.9  1997/02/24  19:12:28  wrankin
 * fixes to timing code
 *
 * Revision 2.8  1997/01/28  17:10:45  wrankin
 * added new timing codes for macroscopic expansion calculations
 *
 * Revision 2.7  1996/10/18  17:04:59  wrankin
 * improved PVM message passing structures.
 * made current and cleaned up T3D library code
 * added additional test codes
 * implement new heirarchical make structure
 *
 * Revision 2.6  1996/09/24  18:42:38  wrankin
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
 * Revision 2.5  1996/08/20  17:12:43  wrankin
 * Enhancement to allow any number of processes to be utilized.
 * Formally the number of precesses was limited to a power of 2.
 *
 * Revision 2.4  1995/07/17  01:13:35  wrankin
 * updates to support SGI Power Challenge (IRIX 6.0.1)
 * cleanup of Makefile
 * initial work on T3D port (not yet supported in this release)
 *
 * Revision 2.3  1995/07/10  02:47:14  wrankin
 * multipole processing code modified to use precomputed mpe transfer
 *   matrices.
 *
 * in addition, all multippole calculation routines have been removed
 *   from this source distribution and will be placed in a separate
 *   multipole library to be supplied externally.
 *
 * Revision 2.2  1995/06/27  14:20:24  wrankin
 * Implementation of relative interaction list algorithm.
 * Interaction list generation code removed from master library.
 * Cell Table restructured as a 2-D table indexed by level/cell
 * All direct computations are now performed in a single pass.
 * General cleanup and removal of old files.
 *
 * Revision 2.1  1995/06/13  04:26:10  wrankin
 * No updates.
 * Sources placed under CVS.
 *
 * Revision 1.1  1995/04/24  04:13:35  wrankin
 * Initial revision
 *
 *
*/

/* include files */

#include <stdio.h>
#if 0
#include <unistd.h>
#endif
#include "dpmta_cell.h"          /* data type definitions */
#include "dpmta_slvglobals.h"         /* global variable declarations */


/*
*  external prototypes
*/

#include "dpmta_distmisc.h"
#include "dpmta_slvmkcell.h"
#include "dpmta_slvmkil.h"
#include "dpmta_slvmkhl.h"
#include "dpmta_slvmcalc.h"
#include "dpmta_slvpcalc.h"
#include "dpmta_slvscale.h"
#ifndef SERIAL
#include "dpmta_slvcomm.h"
#include "dpmta_slvmkiil.h"
#endif
#ifdef TIME
#include "dpmta_timer.h"
#endif


/****************************************************************
*
*  Slave_Init() - initialize and allocate the global data structures
*    used by DPMTA directly as well as initialization of other
*    procedures.
*
*  this routine should be called once right after the global DPMTA
*  variables have been set.
*
*/

void Slave_Init()
{

   /*
   *  misc setups that need to be done first
   */

   Dist_Init( Dpmta_NumLevels );

#ifndef SERIAL
   Comm_Init();
#endif   
   
   /*
   *  allocate cell table and interaction lists
   */

   Alloc_Cell_Table();

   Init_Ilist();

   Init_Hlist();

#ifndef SERIAL
   Init_Inv_Ilist();
#endif

   /*
   *  compute MPE constants used later
   */

   MultipoleSetup();

  
} /* Slave_Init() */


/****************************************************************
*
*  Slave_Start() - perform any needed initializations needed
*    at the beginning of a DPMTA iteration, including any reallocation
*    of bufffers and data structures.
*
*/

void Slave_Start()
{

#if defined VIRIAL || defined OLDVIRIAL
   Dpmta_Vpot = 0.0;
   Dpmta_Vf.x = 0.0;
   Dpmta_Vf.y = 0.0;
   Dpmta_Vf.z = 0.0;
#ifdef COMP_LJ
   Dpmta_Vpot_LJ = 0.0;
   Dpmta_Vf_LJ.x = 0.0;
   Dpmta_Vf_LJ.y = 0.0;
   Dpmta_Vf_LJ.z = 0.0;
#endif
#endif

   if (Dpmta_Resize==TRUE) {

#ifdef TIME
      /* begin timing resize */
      times(&startbuf);
#endif

      /*
         *  compute interaction lists and place in table.
         *  this procedure will also allocate the cell table entries for
         *  the interaction information.
         * 
         *  allocate all remote cells from the information provided in 
         *  the interaction lists.  also allocate the mpe transfer matrices
         */

      Make_Ilist();

      Make_Hlist();

#ifndef SERIAL      
      Alloc_Ilist_Cells();

      Make_Inv_Ilist();
#endif
      
      MultipoleResize();

#ifdef TIME
      times(&endbuf);
      times_arr[8] = (double)(endbuf.tms_utime - startbuf.tms_utime) /
	 (double)CLK_TCK; 
#endif

      Dpmta_Resize = FALSE;

   } /* if Dpmta_Resize */

} /* Slave_Start */


/****************************************************************
*
*  Slave_Compute() - perform all multipole computations
*
*  this routine performs all multipole comutations on the
*  particle data for a single dpmta iteration.
*
*/

void Slave_Compute()
{
      /*
      *  based on inverse interaction lists, distribute particle
      *  cell information to all other processors
      */

#ifndef SERIAL
      Slave_Send_SDirect();
#endif

      /*
      *  perform the upward pass
      *
      *  each slave will process the multipole expansion for the cells
      *  that it owns.  if the higher level cell is owned by another slave,
      *  (based upon the slaved pid) then this slave will pass the multipole
      *  expansions to the appropriate parent.
      */

#ifdef TIME
      /* begin timing step 1 and 2 */
      times(&startbuf);
#endif

      Slave_Mpole_Exp();

#ifdef TIME
      times(&endbuf);
      times_arr[9] = (double)(endbuf.tms_utime - startbuf.tms_utime) / (double)CLK_TCK;
#endif


      /*
      *  collect the particle information from the other processors
      *  and place the information in the local cell table.  note that
      *  we only receive nproc-1 messages.
      */

#ifndef SERIAL
      Slave_Recv_SDirect();
#endif

      /*
      *  based on the inverse interaction list, distribute mpe
      *  cell information to all other processors
      */

#ifndef SERIAL
      Slave_Send_Multipole();
#endif

#ifdef TIME
      gettimeofday(&runstruct,0);
      times_arr[2] = (double)runstruct.tv_sec +
	 ((double)runstruct.tv_usec/1000000);
      times_arr[2] -= times_arr[0];
#endif


      /*
      *  since we have all the local and remote particle information,
      *  do all the direct particle interaction calculations
      *  while we are waiting for the multipole data to 
      *  arrive.
      */

#ifdef TIME
      /* begin timing direct calculations */
      times(&startbuf);
#endif

      Slave_Direct_Calc();

#ifdef TIME
      /* end time for direct calculations */
      times(&endbuf);
      times_arr[10] = (double)(endbuf.tms_utime - startbuf.tms_utime) / (double)CLK_TCK;
#endif


      /*
      *  collect the mpe data from the other processors
      */

#ifndef SERIAL
      Slave_Recv_Multipole();
#endif

#ifdef TIME
      gettimeofday(&runstruct,0);
      times_arr[3] = (double)runstruct.tv_sec +
	 ((double)runstruct.tv_usec/1000000);
      times_arr[3] -= times_arr[0];
#endif


      /*
      *  compute the mpe-particle interactions.  this is also known as
      *  the downward mpe pass.  after all the local expansions have
      *  been computed, then compute the particle forces based upon
      *  the local expansion.
      */

#ifdef TIME
      /* Time downward pass */
      times(&startbuf);
#endif

      Slave_MPE_Calc();

#ifdef TIME
      /* end time for direct calculations */
      times(&endbuf);
      times_arr[11] = (double)(endbuf.tms_utime - startbuf.tms_utime) / (double)CLK_TCK;
#endif

#ifdef TIME
      /* Time downward pass */
      times(&startbuf);
#endif

      Slave_MPE_Force();

#ifdef TIME
      /* end time for direct calculations */
      times(&endbuf);
      times_arr[12] = (double)(endbuf.tms_utime - startbuf.tms_utime) / (double)CLK_TCK;
#endif

#ifdef TIME
      gettimeofday(&runstruct,0);
      times_arr[4] = (double)runstruct.tv_sec +
	 ((double)runstruct.tv_usec/1000000);
      times_arr[4] -= times_arr[0];
#endif


} /* Slave_Compute() */


/****************************************************************
*
*  Slave_Cleanup() - cleans up the contents of the 
*     the cell table at the end of the force calculations.
*
*  this procedure cleans up the cell table data structures at the 
*  end of one iteration, after the force data has been returned to
*  the calling process and is no longer needed by the slave process.
*  this is done in preparation for the next iteration.
*
*  cleanup is performed by cycling through all cells and setting the
*  particle counts to zero and the multipole valid flags (if any)
*  to false.
*
*  since we use realloc() on the particle/force arrays at each iteration,
*  there is no need to implicitly free() them after each iteration.
*
*/

void Slave_Cleanup()
{

   int i;
   int num_cells;


   num_cells = Dpmta_LevelLocate[Dpmta_NumLevels];

   for (i=0; i<num_cells; i++) {

      if ( Dpmta_CellTbl[0][i] != NULL ) {

         Dpmta_CellTbl[0][i]->mvalid = FALSE;
         Dpmta_CellTbl[0][i]->n = 0;

         if (Dpmta_CellTbl[0][i]->mdata != NULL ) {
            Dpmta_CellTbl[0][i]->mdata->lvalid = FALSE;
	 } /* if mdata != NULL */

      } /* if CellTbl[0][i] != NULL */

   } /* for i */

} /* Slave_Cleanup */


/****************************************************************
*
*  Slave_Delete() - free up all global data structures created in
*    Slave_Init() call.
*
*/

void Slave_Delete()
{

   /* free dynamic structures created by the multipole library */
   MultipoleCleanup();
   
   /* free the cell table structures */
   Delete_Cell_Table();

   /* free up the interaction lists */
   Delete_Ilist();
   Delete_Hlist();

#ifndef SERIAL   
   /* free up the inverse interaction lists */
   Delete_Inv_Ilist();
#endif
   
   /* free up other data structures */
   Dist_Delete( Dpmta_NumLevels );

#ifndef SERIAL
   Comm_Delete();
#endif

} /* Slave_Delete */
