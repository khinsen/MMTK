/* 
*  dpmta_slave.c - pvm slave routine for parallel fma algorithm
*
*  w. t. rankin, w. elliott
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*
*  this routine implements the slave processes for the
*  parallel fma algorithm.
*
*  assumptions/limitations:
*
*    - this algorithm assumes a power of two processors.
*    - no error checking on message passing is done.  this will
*      need to be added in a later release.
*    - this slave performs one pass of the dpmta algorithm.  for multiple 
*      time-steps, the algorithm will need to handle redistributing
*      particles (add/delete) after each pass.
*    - this algorithm may be somewhat memory intensive because of the
*      large amounts of message buffering during particle and mpe
*      redistribution.
*/

static char rcsid[] = "$Id: dpmta_slave.c,v 2.34 1998/03/10 22:22:03 wrankin Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_slave.c,v $
 * Revision 2.34  1998/03/10 22:22:03  wrankin
 * folded start/cleanup functionality into dpmta_slvcompute
 *
 * Revision 2.33  1997/11/07 16:49:15  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 * Revision 2.32  1997/05/07 18:59:27  chumphre
 * major changes to implement ||-piped cells
 *
 * Revision 2.31  1997/03/26  20:36:22  wrankin
 * rewrote Ilist and IIlist code to handle reallocation of array sizes
 *   that happen as a resuly of re-scaling a non-cubic unit cell.
 *
 * Revision 2.30  1997/02/26  20:43:29  wrankin
 * updated timing measurements and placed routines in single module
 *
 * Revision 2.29  1997/02/26  18:52:34  wrankin
 * cleaned up handling of cell resizing
 *
 * Revision 2.28  1997/02/24  19:12:25  wrankin
 * fixes to timing code
 *
 * Revision 2.27  1997/01/28  17:10:42  wrankin
 * added new timing codes for macroscopic expansion calculations
 *
 * Revision 2.26  1996/11/18  19:29:29  wrankin
 * update to virial code to use multipole virial calculation
 * old virial code is turned on by defining OLDVIRIAL
 * direct test program now computes absolute virial tensor.
 *
 * Revision 2.25  1996/11/14  17:50:23  wrankin
 * performance enhancements to multipole M2L routine.
 * additions to make slaves exit gracefully.
 *
 * Revision 2.24  1996/10/18  17:04:50  wrankin
 * improved PVM message passing structures.
 * made current and cleaned up T3D library code
 * added additional test codes
 * implement new heirarchical make structure
 *
 * Revision 2.23  1996/09/24  18:42:04  wrankin
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
 * Revision 2.22  1996/08/20  17:12:36  wrankin
 * Enhancement to allow any number of processes to be utilized.
 * Formally the number of precesses was limited to a power of 2.
 *
 * Revision 2.21  1996/08/09  15:30:46  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.20  1996/02/29  21:13:32  wrankin
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
 * Revision 2.19  1996/02/12  15:20:55  wrankin
 * Added Macroscopic Assemblies code.
 *
 * Revision 2.18  1995/12/08  23:00:21  wrankin
 * preliminary release of DPMTA 2.3
 *   - added working Periodic Boundary Conditions (PDC) flag
 *   - added code for Virial Computation (not working yet)
 *   - general cleanup of some modules
 *
 * Revision 2.17  1995/11/29  22:29:13  wrankin
 * addition of periodic boundary code
 * beginning work on virial code.
 *
 * Revision 2.16  1995/10/01  21:45:58  wrankin
 * implemented LJ multipole calculations - compile time option only
 *
 * Revision 2.15  1995/09/18  19:02:34  wrankin
 * updates to move the scaling of the simulation box to the slave
 * process.
 *
 * Revision 2.14  1995/07/01  03:26:56  wrankin
 * initial cut at precomputation of mpe transfer functions.
 * this works for vanilla M2L.  it does not work for FFT enhancements.
 *
 * Revision 2.13  1995/06/27  14:20:18  wrankin
 * Implementation of relative interaction list algorithm.
 * Interaction list generation code removed from master library.
 * Cell Table restructured as a 2-D table indexed by level/cell
 * All direct computations are now performed in a single pass.
 * General cleanup and removal of old files.
 *
 * Revision 2.12  1995/05/16  04:02:16  wrankin
 * added code to allow resizing the sim cube on the fly
 * the calling process now scales the input to a unit cube
 * and scales the results back to the proper value
 * the slave assumes all data sent to it is in a unit cube
 * with the bottom corner at <0,0,0>
 *
 * Revision 2.11  1995/04/09  22:21:06  wrankin
 * renaming of all global variables to prevent conflicts for
 * port of dpmta library to the T3D.
 *
 * Revision 2.10  1995/03/06  07:57:44  wrankin
 * slave process now receives Theta parameter from master
 * slave process now computes inv. interaction list locally
 *
 * Revision 2.9  1995/01/13  16:07:26  wrankin
 * changed temporary multipole type
 *
 * Revision 2.8  1994/11/30  18:05:26  wrankin
 * added globals for distributed calling sequence
 * timing information now routed correctly for distributed calling sequence
 *
 * Revision 2.7  1994/11/28  04:27:53  wrankin
 * added hooks to handle multiple slave calls to PMTAforce
 *
 * Revision 2.6  1994/11/17  00:27:52  wrankin
 * added code to handled variable cell center
 *
 * Revision 2.5  1994/10/26  02:42:09  wrankin
 * added capability to do multiple iterations
 *
 * Revision 2.4  1994/10/14  04:58:07  wrankin
 * added duke copyright notice
 * added fftblock global variable
 * removed lots of debugging statements
 * removed T3DBUG workarounds
 * general cleanup of commenting.
 *
 * Revision 2.3  1994/08/25  18:53:18  wrankin
 * added T3DBUG message forwarding
 * additional timing information provided to master
 * added call to pvm_joingroup to fix barrier problem
 *
 * Revision 2.2  1994/07/25  23:08:40  wrankin
 * fixed message passing in upward pass portion of code
 * fixed mpe messages that were only sending half the complex data
 *
 * Revision 2.1  1994/06/02  12:52:27  wrankin
 * No change.  Update version number for Release 2
 *
 * Revision 1.21  1994/06/02  03:34:34  wrankin
 * added pvm_barrier() call to synchronize timing runs of slaves.
 *
 * Revision 1.20  1994/05/27  06:31:19  wrankin
 * added timing include files for cray-t3d
 * cleaned up debugging output
 *
 * Revision 1.19  1994/04/19  19:01:30  wrankin
 * added more timing code from uiuc
 * added conditional compilation code from HPPA and LINUX
 *
 * Revision 1.18  1994/04/09  04:16:43  wrankin
 * swapped particle and interaction list processing
 *
 * Revision 1.17  1994/04/09  03:46:08  wrankin
 * added dan gray's performance timing code
 *
 * Revision 1.16  1994/03/31  13:32:58  wrankin
 * added new global CubeLength, passed in from master program
 *
 * Revision 1.15  1994/03/27  20:29:44  wrankin
 * added MultipoleSetup() cell to precompute MPE constants
 * added global array definitions to support MPE
 *
 * Revision 1.14  1994/03/26  16:30:35  wrankin
 * added new global declarations for FFT, Mp, and LclSize
 * enables the cells to the multipole expansion/calculation routines
 *
 * Revision 1.13  1994/03/18  17:40:28  wrankin
 * slave now performs double and single direct interactions
 * debuging output is now provided to a /tmp file
 * slave does not perform MPE
 *
 * Revision 1.12  1994/03/14  15:29:16  wrankin
 * added calls to return force vectors to the master process
 *
 * Revision 1.11  1994/03/11  19:46:28  wrankin
 * added global variables and removed parameter passing to subroutines
 *
 * Revision 1.10  1994/03/10  04:49:17  wrankin
 * contains calls to all routines to calculate particle interactions
 *
 * Revision 1.9  1994/02/22  03:58:04  wrankin
 * added calls to distribute particle info between slaves
 *
 * Revision 1.8  1994/02/16  19:24:30  wrankin
 * fixed interaction list message processing
 *
 * Revision 1.7  1994/02/16  14:09:32  wrankin
 * added call to Alloc_Ilist_Cells() - untested.
 *
 * Revision 1.6  1994/02/16  01:04:27  wrankin
 * interaction list processing now works
 *
 * Revision 1.5  1994/02/12  15:48:56  wrankin
 * updated subroutine interface to use sinfo structure
 *
 * Revision 1.4  1994/02/01  20:04:01  wrankin
 * added sleep() stubs to simulate fma processing between message sending
 *
 * Revision 1.3  1994/01/31  11:03:11  wrankin
 * fixed compilation errors - code compiles under gcc.
 *
 * Revision 1.2  1994/01/31  10:19:24  wrankin
 * filled in message passing sections with dummy calls
 * code is still untested
 *
 * Revision 1.1  1994/01/30  17:29:03  wrankin
 * Initial revision
 *
*/

/* include files */

#include <stdio.h>
#include <unistd.h>
#include "pvm3.h"
#include "dpmta_pvm.h"           /* pvm messaging declarations */
#include "dpmta_cell.h"          /* data type definitions */
#include "dpmta_slvglobals.h"    /* global variables */

/*
 * prototypes - lots of them here, because this is
 *   the top level routine and calls everything else.
 */

#include "dpmta_slvcompute.h"
#include "dpmta_slvcomm.h"
#include "dpmta_slvscale.h"

#ifdef TIME
#include "dpmta_timer.h"
#endif


/****************************************************************/

int main( int argc, char *argv[] )
{

   /* enroll in pvm */
   if ((Dpmta_MyTid = pvm_mytid()) < 0) {
        fprintf(stderr,"error: exiting\n");
        fflush(stderr);
        exit(0);
   }


   /* determine master id */
   Dpmta_MasterTid = pvm_parent();


   /****************************************************************
   *
   *  receive initialization message from master
   */

   Recv_Master_Info();

#ifdef TIME
   pvm_joingroup("DpmtaSlave");
#endif

   /****************************************************************
   *
   *  initialize buffers and data structures
   *
   */

   Slave_Init();


   /****************************************************************
   *
   *  This is the beginning of the particle processing.  We will 
   *  loop continuously waiting for new sets of particles, computing
   *  the forces, sending the results to the master process, cleaning
   *  up the tempoary data structures, and going back to wait for
   *  more particles.
   *
   *  Recv_Particles() - returns TRUE as long as it receives valid
   *  particle data.
   *
   */

   while ( Recv_Particles() ) {

#ifdef TIME
      {
	 int i;
	 for (i=0; i<TIMINGSIZE; i++) {
	    times_arr[i] = 0.0;
	 }
      }
#endif

#ifdef TIME
      /* synchronized start */
      pvm_barrier("DpmtaSlave",Dpmta_Nproc);
      gettimeofday(&runstruct,0);
      times_arr[0] = (double)runstruct.tv_sec +
	 ((double)runstruct.tv_usec/1000000);
#endif

      Slave_Start();

      Rescale_Particles();


      /*
      *  compute all forces and potentials
      */

      Slave_Compute();

      
      /*
      *  send completed data to parent
      */

      Rescale_Results();

      Send_Results();

#ifdef TIME
      gettimeofday(&runstruct,0);
      times_arr[1] = (double)runstruct.tv_sec +
	 ((double)runstruct.tv_usec/1000000);
      times_arr[1] -= times_arr[0];
#endif

#ifdef TIME
   /*
   *  send timing results to master
   */
   
   if ( Dpmta_CallingNum == 0 )
      Send_Slave_Times(Dpmta_Pid, times_arr, Dpmta_MasterTid);
   else
      Send_Slave_Times(Dpmta_Pid, times_arr, Dpmta_CallingTids[0]);

#endif

      /*
      *  clean and free up the particle allocations and other data structures
      *  that are re-allocated on each time step.  note that this can be done
      *  after normal processing has been complete so it does not show up
      *  as part of the timings.
      */

      Slave_Cleanup();

      /*
      *  end of single iteration.  go back and wait for new particles
      *  if necessary
      */

   } /* while Recv_Particles() */

   /*
   *  we have recived a termination message -
   *  perform any needed final cleanup and exit
   */

   Slave_Delete();

   pvm_exit();
   exit(0);
   
} /* Main */


