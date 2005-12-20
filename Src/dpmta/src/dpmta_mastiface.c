/* 
*  dpmta_mastiface.c - Distributed PMTA interface routines
*
*  w. t. rankin
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*  this program provides the interface routines between application 
*  programs and the distributed PMTA routines.
*
*/

static char rcsid[]="$Id: dpmta_mastiface.c,v 2.21 1998/04/01 20:08:08 wrankin Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_mastiface.c,v $
 * Revision 2.21  1998/04/01 20:08:08  wrankin
 * added support for HILBERT ordering
 * general cleanup of code structure (more OOP if you will)
 *
 * Revision 2.20  1998/01/05 15:35:31  wrankin
 * cleaned up processing of PMTAregister() and PMTAexit().
 *
 * Revision 2.19  1997/11/07 16:49:04  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 * Revision 2.18  1997/05/12 18:06:01  wrankin
 * added routines to clean up dynamic memory allocation when
 * PMTAinit() is called
 *
 * Revision 2.17  1997/05/07  21:27:46  chumphre
 * cubectr changed to cellctr and other minor changes related to the interface
 *
 * Revision 2.16  1997/05/07  19:31:44  chumphre
 * implement a uniform interface for rectangular and ||-piped cells
 *
 * Revision 2.15  1997/05/07  18:59:19  chumphre
 * major changes to implement ||-piped cells
 *
 * Revision 2.14  1997/04/28  16:07:10  wrankin
 * number of calling processes is now passed in CallingNum, regardless
 *   of if there is a single master calling routine or not.
 *
 * Revision 2.13  1997/02/26  18:52:30  wrankin
 * cleaned up handling of cell resizing
 *
 * Revision 2.12  1997/02/24  19:12:17  wrankin
 * fixes to timing code
 *
 * Revision 2.11  1996/11/18  19:29:26  wrankin
 * update to virial code to use multipole virial calculation
 * old virial code is turned on by defining OLDVIRIAL
 * direct test program now computes absolute virial tensor.
 *
 * Revision 2.10  1996/11/14  17:50:12  wrankin
 * performance enhancements to multipole M2L routine.
 * additions to make slaves exit gracefully.
 *
 * Revision 2.9  1996/10/18  17:04:37  wrankin
 * improved PVM message passing structures.
 * made current and cleaned up T3D library code
 * added additional test codes
 * implement new heirarchical make structure
 *
 * Revision 2.8  1996/09/24  18:41:45  wrankin
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
 * Revision 2.7  1996/08/09  15:30:38  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.6  1996/02/29  21:13:24  wrankin
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
 * Revision 2.5  1995/11/29  22:29:03  wrankin
 * addition of periodic boundary code
 * beginning work on virial code.
 *
 * Revision 2.4  1995/10/01  21:45:48  wrankin
 * implemented LJ multipole calculations - compile time option only
 *
 * Revision 2.3  1995/09/18  19:02:27  wrankin
 * updates to move the scaling of the simulation box to the slave
 * process.
 *
 * Revision 2.2  1995/06/27  14:20:17  wrankin
 * Implementation of relative interaction list algorithm.
 * Interaction list generation code removed from master library.
 * Cell Table restructured as a 2-D table indexed by level/cell
 * All direct computations are now performed in a single pass.
 * General cleanup and removal of old files.
 *
 * Revision 2.1  1995/05/16  04:02:16  wrankin
 * added code to allow resizing the sim cube on the fly
 * the calling process now scales the input to a unit cube
 * and scales the results back to the proper value
 * the slave assumes all data sent to it is in a unit cube
 * with the bottom corner at <0,0,0>
 *
 * Revision 1.7  1995/03/06  07:45:08  wrankin
 * master now sends Theta parameter to all slaves
 *
 * Revision 1.6  1995/02/24  21:19:18  wrankin
 * added new theta paramenter to PMTAinit() for ilist generation
 *
 * Revision 1.5  1994/12/06  19:04:05  wrankin
 * added code to handle cell center offsets
 *
 * Revision 1.4  1994/11/30  17:59:57  wrankin
 * added code to handle distributed calling sequence
 * interface now uses types from dpmta.h
 * all procedures now return values
 *
 * Revision 1.3  1994/11/28  04:27:53  wrankin
 * added hooks to handle multiple slave calls to PMTAforce
 *
 * Revision 1.2  1994/11/17  00:27:52  wrankin
 * added code to handled variable cell center
 *
 * Revision 1.1  1994/10/26  02:37:29  wrankin
 * Initial revision
 *
 *
*/


/* include files */
#include <stdio.h>
#include "pvm3.h"
#include "dpmta.h"
#include "dpmta_pvm.h"

/*
 * external prototypes
 */

#include "dpmta_mastiter.h"
#include "dpmta_distmisc.h"
#ifdef TIME
#include "dpmta_timer.h"
#endif

/* some global defines */

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/*
*  globals used to communicate between the different routines
*/

static int MyProcNum;              /* slaves index into CallingTids */
static int Tids[MAXPROC];          /* tids for the DPMTA slave processes */
static int Nprocs;                 /* number of DPMTA slave processes */
static int CallingTids[MAXPROC];   /* tids for the application slaves */
static int CallingNum;             /* number of application slaves*/
static int Nlevels;                /* number of levels of spactial decomp */
static PmtaVector CellVector1;        /* vectors defining edges of cell */
static PmtaVector CellVector2;
static PmtaVector CellVector3;
static PmtaVector CellCenter;      /* center coordinates of simulation space */
static int Resize;                 /* flag to indicate a cell resize */

/****************************************************************
*
*  PMTAinit() - initialize slave processes.
*
*  this routine creates the PVM slave processes and sends them
*  initializing data.  it makes the assumption that the celling 
*  process is already enrolled in PVM.
*
*  the tids of the slave processes are returned in the array
*  rtn_tids[].  the user is responsible for allocating enough
*  space to hold the tids.
*
*/

int PMTAinit(PmtaInitDataPtr initdata, int *rtn_tids)
{

   int i;
   int p;

   /* if this is called from a master, the proc id is 0 */
   MyProcNum = 0;

   /* set globals */
   Nprocs  = initdata->nprocs; 
   Nlevels = initdata->nlevels;

   CallingNum = initdata->calling_num;
   if ( CallingNum > 0 ) {
      for ( i=0; i<CallingNum; i++ ) {
         CallingTids[i] = initdata->calling_tids[i];
      }
   }
   else {
      CallingNum = 1;
      CallingTids[0] = pvm_mytid();
   }

   CellVector1.x = initdata->v1.x;
   CellVector1.y = initdata->v1.y;
   CellVector1.z = initdata->v1.z;
   CellVector2.x = initdata->v2.x;
   CellVector2.y = initdata->v2.y;
   CellVector2.z = initdata->v2.z;
   CellVector3.x = initdata->v3.x;
   CellVector3.y = initdata->v3.y;
   CellVector3.z = initdata->v3.z;

   CellCenter.x = initdata->cellctr.x;
   CellCenter.y = initdata->cellctr.y;
   CellCenter.z = initdata->cellctr.z;

   Resize = TRUE;

   /* start up slave tasks */
   p = pvm_spawn("dpmta_slave",NULL,0,"",Nprocs,Tids);
   if (Nprocs != p) {
      fprintf(stderr,"Error: only spawned %d of %d processes, exiting\n",
         p, Nprocs);
      return(-1);
   }

   /*
   *  load tids into return array.  note that we do not just return the
   *  address to our global Tids array, because we dont want to run the 
   *  risk of having it modified by the user process.  so we must copy
   *  the data.
   */

   for ( i=0; i<Nprocs; i++ )
      rtn_tids[i] = Tids[i];

   /*
    * initialize particle distributions arrays
    */

   Dist_Init( Nlevels );


   /*
   *  send initial data to slaves
   */

   Send_Slave_Info(Nprocs,Tids,Nlevels,initdata->fft,initdata->pbc,
      initdata->mp,initdata->mp_lj,initdata->kterm,initdata->theta,
      initdata->fftblock,CallingNum,CallingTids);

   /*
   *  send global data to other registered calling slaves if needed
   *
   *  i probably need to move this off into a subroutine.
   *  next release, i promise!
   *
   *  on the other hand, if there is only one calling process, then
   *  we need to take care of a couple initializations that would
   *  otherwise be called in PMTAregister().
   *
   *  on the other hand, we will just make the limitation that
   *  only in the case where CallingNum is set to zero, we will
   *  allow the application to *not* call PMTAregister.  so
   *  the check will be for CallingNum > 0.
   *
   */

   for ( i=0; i<CallingNum; i++ ) {
      pvm_initsend(DATA_NORMAL_PVM);
      pvm_pkint(&i,1,1);
      pvm_pkint(&Nprocs,1,1);
      pvm_pkint(Tids,Nprocs,1);
      pvm_pkint(&CallingNum,1,1);
      pvm_pkint(CallingTids,CallingNum,1);
      pvm_pkint(&Nlevels,1,1);
      pvm_pkint(&Resize,1,1);
      pvm_pkdouble(&(CellVector1.x),3,1);
      pvm_pkdouble(&(CellVector2.x),3,1);
      pvm_pkdouble(&(CellVector3.x),3,1);
      pvm_pkdouble(&(CellCenter.x),3,1);
      pvm_send(CallingTids[i],MSG_INIT2);
   } /* for i */

   return(0);
   
} /* PMTAinit */



/****************************************************************
*
*  PMTAregister() - register application slave processes.
*
*  this routine downloads the initializing data from the master process.
*  this data will be needed in order to correctly distribute the particle
*  data to the dpmta slaves.
*
*  if this process has called DPMTAinit() then it still need to call 
*  this routine if it wants to get rid of the outstanding initialization
*  message to itself.
*
*  message format -
*
*    int - my processor index (*not* the tid!)
*    int - number of dpmta slaves
*    int[] - tids of dpmta slaves
*    int - number of application slaves
*    int[] - tids of application slaves
*    int - number of decomp levels
*    double[3] - cell length
*    double[3] - cell center
*
*/


int PMTAregister()
{

   pvm_recv(-1,MSG_INIT2);

   pvm_upkint(&MyProcNum,1,1);
   pvm_upkint(&Nprocs,1,1);
   pvm_upkint(Tids,Nprocs,1);
   pvm_upkint(&CallingNum,1,1);
   pvm_upkint(CallingTids,CallingNum,1);
   pvm_upkint(&Nlevels,1,1);
   pvm_upkint(&Resize,1,1);
   pvm_upkdouble(&(CellVector1.x),3,1);
   pvm_upkdouble(&(CellVector2.x),3,1);
   pvm_upkdouble(&(CellVector3.x),3,1);
   pvm_upkdouble(&(CellCenter.x),3,1);
   Resize = TRUE;

   Init_Local_Buffers();

   return(0);
   
} /* PMTAregister */



/****************************************************************
*
*  PMTAforce - calculate colomb forces
*
*  send particles to slaves
*  wait for returned forces and potential energies
*  prints out timing information if needed
* 
*/

int PMTAforce( 
   int nparts,
   PmtaParticlePtr particles,
   PmtaPartInfoPtr results,
   PmtaPartInfoPtr results_lj )
{


   Send_Slave_Particles(nparts, particles, Nprocs, Tids,
      Nlevels, Resize, &(CellVector1), &(CellVector2), &(CellVector3),
      &(CellCenter), MyProcNum );

   Resize = FALSE;

   Recv_Slave_Results(nparts, results, results_lj,
		      Nprocs );

#ifdef TIME
   Recv_Slave_Times(MyProcNum,Nprocs);
#endif

   return(0);

} /* PMTAforce */


/****************************************************************
*
*  PMTAresize() - resizes the simulation cube
*
*  this routine changes the values of the CellCenter and CubeLength
*  globals, which are used to scale the particle results.
*
*  note that for the distributed case, everybody need to make this call.
*
*
*/

int PMTAresize(
   PmtaVector *pvector1,
   PmtaVector *pvector2,
   PmtaVector *pvector3,
   PmtaVector *cellctr )  /* center of simulation cube */
{

   CellVector1.x = pvector1->x;
   CellVector1.y = pvector1->y;
   CellVector1.z = pvector1->z;
   CellVector2.x = pvector2->x;
   CellVector2.y = pvector2->y;
   CellVector2.z = pvector2->z;
   CellVector3.x = pvector3->x;
   CellVector3.y = pvector3->y;
   CellVector3.z = pvector3->z;

   CellCenter.x = cellctr->x;
   CellCenter.y = cellctr->y;
   CellCenter.z = cellctr->z;

   Resize = TRUE;

   return(0);
   
} /* PMTAresize */


/****************************************************************
 *
 *  PMTAvirial() - return the virial sums from the previous
 *
 */

int PMTAvirial(
   double *vp,
   PmtaVector *vf,
   double *vp_lj,
   PmtaVector *vf_lj )
{

#if defined VIRIAL || defined OLDVIRIAL
   Return_Virial(vp,vf,vp_lj,vf_lj);
   return(0);
#else
   return(-1);
#endif

} /* PMTAvirial */


/****************************************************************
*
*  PMTAexit() - shuts down slave processes
*
*  kills off all slaves processes.  all processes that called
*  PMTAregister() must call this routine.
*
*  to allow graceful termination of the slave processes, this routine
*  sends a message to the slaves that will cause the slave to exit
*  on its own, as opposed to using pvm_kill().
*
*  this procedure does not exit PVM.  it assumes that the calling
*  process will take care of this.
*
*  in addition, this process will free up all the dynamic data
*  structures created by the master side of DPMTA.
*
*/

int PMTAexit()
{

   int i;
   int msg_type;

   msg_type = MSG_EXIT;

   for (i=0; i<Nprocs; i++) {
      pvm_initsend(DATA_NORMAL_PVM);
      pvm_pkint(&msg_type,1,1);
      pvm_send(Tids[i],MSG_PART1);
   }

   /* clean up local allocations */
   Delete_Local_Buffers();

   Dist_Delete( Nlevels );
   return(0);

} /* PMTAexit */
