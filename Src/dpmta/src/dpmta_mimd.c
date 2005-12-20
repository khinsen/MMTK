/*
*  dpmta_t3dlib - Cray T3D library interface.
*
*  w. t. rankin
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*  these routines merge the functionality of dpmta_slave and the
*  master interface routines into a single module for execution
*  on the Cray T3D.
*
*
*/

static char rcsid[] = "$Id: dpmta_mimd.c,v 2.3 1998/03/10 22:21:59 wrankin Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_mimd.c,v $
 * Revision 2.3  1998/03/10 22:21:59  wrankin
 * folded start/cleanup functionality into dpmta_slvcompute
 *
 * Revision 2.2  1998/01/05 13:42:00  wrankin
 * added more timing code.
 *
 * Revision 2.1  1997/11/07 16:49:11  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 * Revision 2.20  1997/05/14 17:06:05  chumphre
 * Implemented ||-Piped changes for t3d/e
 *
 * Revision 2.19  1997/05/12  18:06:16  wrankin
 * added routines to clean up dynamic memory allocation when
 * PMTAinit() is called
 *
 * Revision 2.18  1997/05/07  21:27:55  chumphre
 * cubectr changed to cellctr and other minor changes related to the interface
 *
 * Revision 2.17  1997/04/28  16:01:45  wrankin
 * fix to timer code to make consistant with distributed impl.
 * fit to handling of cell length to make consistant with distributed impl.
 *
 * Revision 2.16  1997/03/26  20:36:35  wrankin
 * rewrote Ilist and IIlist code to handle reallocation of array sizes
 *   that happen as a resuly of re-scaling a non-cubic unit cell.
 *
 * Revision 2.15  1997/02/26  20:43:32  wrankin
 * updated timing measurements and placed routines in single module
 *
 * Revision 2.14  1997/02/26  18:52:39  wrankin
 * cleaned up handling of cell resizing
 *
 * Revision 2.13  1997/02/24  19:12:30  wrankin
 * fixes to timing code
 *
 * Revision 2.12  1996/11/18  19:29:42  wrankin
 * update to virial code to use multipole virial calculation
 * old virial code is turned on by defining OLDVIRIAL
 * direct test program now computes absolute virial tensor.
 *
 * Revision 2.11  1996/11/01  02:26:48  wrankin
 * modifications to support cray t3d compilation
 * version update for 2.5 release
 *
 * Revision 2.10  1996/10/18  17:05:06  wrankin
 * improved PVM message passing structures.
 * made current and cleaned up T3D library code
 * added additional test codes
 * implement new heirarchical make structure
 *
 * Revision 2.9  1996/09/24  18:43:55  wrankin
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
 * Revision 2.8  1996/08/09  15:31:14  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.7  1996/02/29  21:13:59  wrankin
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
 * Revision 2.6  1995/10/01  21:46:25  wrankin
 * implemented LJ multipole calculations - compile time option only
 *
 * Revision 2.5  1995/09/18  19:02:48  wrankin
 * updates to move the scaling of the simulation box to the slave
 * process.
 *
 * Revision 2.4  1995/07/26  01:52:58  wrankin
 * updated Makefiles
 * slight performance improvement for direct calculations
 * clean up test code example
 * fixed t3dlib for new multipole code
 *
 * Revision 2.3  1995/07/17 01:13:36  wrankin
 * updates to support SGI Power Challenge (IRIX 6.0.1)
 * cleanup of Makefile
 * initial work on T3D port (not yet supported in this release)
 *
 * Revision 2.2  1995/05/18  04:28:34  wrankin
 * added section to receive and print timing information from
 * other slaves (all the other ecode was already there).
 *
 * Revision 2.1  1995/05/16  04:02:16  wrankin
 * added code to allow resizing the sim cube on the fly
 * the calling process now scales the input to a unit cube
 * and scales the results back to the proper value
 * the slave assumes all data sent to it is in a unit cube
 * with the bottom corner at <0,0,0>
 *
 * Revision 1.1  1995/04/09  22:15:33  wrankin
 * Initial revision
 *
 *
 */


/* include files */
#include <stdio.h>
#include "pvm3.h"
#include "dpmta.h"
#include "dpmta_pvm.h"
#include "dpmta_cell.h"
#include "dpmta_slvglobals.h"


/*
*  some prototypes
*/

#include "dpmta_mastiter.h"
#include "dpmta_slvcompute.h"
#include "dpmta_slvcomm.h"
#include "dpmta_slvscale.h"

#ifdef TIME
#include "dpmta_timer.h"
#endif



/****************************************************************
*
*  PMTAinit() - initialize slave processes.
*
*  this routine creates the PVM slave processes and sends them
*  initializing data.  it makes the assumption that the celling 
*  process is already enrolled in PVM.
*
*  for the T3D interface, it is assumed that only processor #0
*  will call this routine.  the application is responsible for
*  assuring this.
*
*  the tids of the slave processes are returned in the array
*  rtn_tids[].  the user is responsible for allocating enough
*  space to hold the tids.  since on the t3d, the slaves have
*  already been started and are actually the calling processes,
*  "rtn_tids" is simply the contents of "calling_tids".
*
*  note also that we should probably verify that 
*  calling_num == nprocs.
*
*/

int PMTAinit( PmtaInitDataPtr initdata, int *rtn_tids )
{
   int i;              /* loop counters */

   /*
    * in the case of the cray, the calling tids are same
    * as the return tids (they are the same processes.
    * therefor, we need to copy the values for the
    * calling tids to the rtn_tid array
    */

   for (i=0; i<(initdata->nprocs); i++) {
      rtn_tids[i] = initdata->calling_tids[i];
   }

   /*
   *  send global data to other registered slaves
   *  (this is received in the PMTA register routine...)
   *
   */

   Send_Slave_Info(initdata->nprocs, rtn_tids,
      initdata->nlevels,initdata->fft,initdata->pbc,
      initdata->mp,initdata->mp_lj,initdata->kterm,
      initdata->theta,initdata->fftblock,
      initdata->calling_num,initdata->calling_tids);

   /*
    * each of the other slaves needs to know the size of the
    * simulation cell so that this information can be included
    * in the header of the particle message
    */

   if ( initdata->calling_num > 1 ) {
      for ( i=0; i<(initdata->calling_num); i++ ) {
	 pvm_initsend(DATA_INPLACE_PVM);
	 pvm_pkdouble(&(initdata->v1.x),3,1);
	 pvm_pkdouble(&(initdata->v2.x),3,1);
	 pvm_pkdouble(&(initdata->v3.x),3,1);
	 pvm_pkdouble(&(initdata->cellctr.x),3,1);
	 pvm_send(initdata->calling_tids[i],MSG_INIT2);
      } /* for i */
   }
   else {
      Dpmta_CellVector1.x = initdata->v1.x;
      Dpmta_CellVector1.y = initdata->v1.y;
      Dpmta_CellVector1.z = initdata->v1.z;
      Dpmta_CellVector2.x = initdata->v2.x;
      Dpmta_CellVector2.y = initdata->v2.y;
      Dpmta_CellVector2.z = initdata->v2.z;
      Dpmta_CellVector3.x = initdata->v3.x;
      Dpmta_CellVector3.y = initdata->v3.y;
      Dpmta_CellVector3.z = initdata->v3.z;
#ifdef PIPED
      Dpmta_CV1Mag = Vec_Mag(&Dpmta_CellVector1);
      Dpmta_CV2Mag = Vec_Mag(&Dpmta_CellVector2);
      Dpmta_CV3Mag = Vec_Mag(&Dpmta_CellVector3);
#endif
      Dpmta_CellCenter.x = initdata->cellctr.x;
      Dpmta_CellCenter.y = initdata->cellctr.y;
      Dpmta_CellCenter.z = initdata->cellctr.z;
      Dpmta_MaxCellLen = Max_CellLength();
      Dpmta_Resize = TRUE;
   }

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
*  for the T3D version, this process also performs all the initialzation
*  steps formerly run by the dpmta_slave slave process upon its creation.
*
*/

int PMTAregister()
{

   /****************************************************************
   *
   *  receive initialization message from master
   */

   Recv_Master_Info();

   if ( Dpmta_CallingNum > 1 ) {
      pvm_recv(-1,MSG_INIT2);
      pvm_upkdouble(&(Dpmta_CellVector1.x),3,1);
      pvm_upkdouble(&(Dpmta_CellVector2.x),3,1);
      pvm_upkdouble(&(Dpmta_CellVector3.x),3,1);
#ifdef PIPED
      Dpmta_CV1Mag = Vec_Mag(&Dpmta_CellVector1);
      Dpmta_CV2Mag = Vec_Mag(&Dpmta_CellVector2);
      Dpmta_CV3Mag = Vec_Mag(&Dpmta_CellVector3);
#endif
      pvm_upkdouble(&(Dpmta_CellCenter.x),3,1);
      Dpmta_MaxCellLen = Max_CellLength();
      Dpmta_Resize = TRUE;
   }

#ifdef TIME
   pvm_joingroup("DpmtaSlave");
#endif

   /*
   *  now each slave needs to process this information and
   *  allocate the appropriate data structures and constants 
   */

   Slave_Init();
   
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

#ifdef TIME
   {
      int i;
      for (i=0; i<TIMINGSIZE; i++)
	  times_arr[i] = 0.0;
   }
#endif

#ifdef TIME
   /* synchronized start */
   pvm_barrier("DpmtaSlave",Dpmta_Nproc);
   gettimeofday(&runstruct,0);
   times_arr[0] = (double)runstruct.tv_sec +
      ((double)runstruct.tv_usec/1000000);

#endif


   /*
   *  redistribute all particles to the other processors.
   *
   *  this routine also performs particle scaling to map the simulation
   *  space of the calling process to the unit cube of the slave.
   *
   *  we should probably combine this routine with the receive routine
   *  inorder to interleave the sending and receiving of messages so 
   *  that there is only one message out at any one time (saving 
   *  message buffers).
   */

   Send_Slave_Particles(nparts, particles,
      Dpmta_Nproc, Dpmta_Tids, Dpmta_NumLevels, Dpmta_Resize,
      (PmtaVector *)(&Dpmta_CellVector1), 
      (PmtaVector *)(&Dpmta_CellVector2), 
      (PmtaVector *)(&Dpmta_CellVector3),
      (PmtaVector *)(&Dpmta_CellCenter), Dpmta_Pid);

#ifdef TIME
   gettimeofday(&runstruct,0);
   times_arr[5] = (double)runstruct.tv_sec +
      ((double)runstruct.tv_usec/1000000);
   times_arr[5] -= times_arr[0];
#endif


   /*
   *  receive particle lists and place in cell table.
   */

   Recv_Particles();

#ifdef TIME
   gettimeofday(&runstruct,0);
   times_arr[6] = (double)runstruct.tv_sec +
      ((double)runstruct.tv_usec/1000000);
   times_arr[6] -= times_arr[0];
#endif


#ifdef TIME
   /* begin timing resize */
   times(&startbuf);
#endif

   Slave_Start();

#ifdef TIME
   times(&endbuf);
   times_arr[8] = (double)(endbuf.tms_utime - startbuf.tms_utime) / (double)CLK_TCK;
#endif

   Rescale_Particles();


   /*
   *  compute all forces and potentials
   */

   Slave_Compute();


   /*
   *  send completed data to owning processor.
   */

   Rescale_Results();

   Send_Results();


   /*
   *  clean and free up the particle allocations and other data
   *  structures that are re-allocated on each time step.  perform
   *  this while the particle information is in transit.
   */

   Slave_Cleanup();


   /*
   *  collect the results from each other process
   */

   Recv_Slave_Results(nparts, results, results_lj, Dpmta_Nproc);

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

   /*
   *  if needed, receive the timing data from all the other processes
   */

   Recv_Slave_Times(Dpmta_Pid, Dpmta_Nproc);

#endif

   return(0);

} /* PMTAforce */



/****************************************************************
*
*  PMTAresize() - resizes the simulation cube
*
*  this routine changes the values of the CubeCenter and CubeLength
*  globals, which are used to scale the particle results.
*
*  note that for the distributed case, everybody need to make this call.
*
*
*/

int PMTAresize(
   PmtaVector *pvector1,
   PmtaVector *pvector2,
   PmtaVector *pvector3,   /* length of cell edges */
   PmtaVector *cellctr )   /* center of simulation cell */

{

   Dpmta_CellVector1.x = pvector1->x;
   Dpmta_CellVector1.y = pvector1->y;
   Dpmta_CellVector1.z = pvector1->z;
   Dpmta_CellVector2.x = pvector2->x;
   Dpmta_CellVector2.y = pvector2->y;
   Dpmta_CellVector2.z = pvector2->z;
   Dpmta_CellVector3.x = pvector3->x;
   Dpmta_CellVector3.y = pvector3->y;
   Dpmta_CellVector3.z = pvector3->z;
#ifdef PIPED
   Dpmta_CV1Mag = Vec_Mag(&Dpmta_CellVector1);
   Dpmta_CV2Mag = Vec_Mag(&Dpmta_CellVector2);
   Dpmta_CV3Mag = Vec_Mag(&Dpmta_CellVector3);
#endif
   Dpmta_CellCenter.x = cellctr->x;
   Dpmta_CellCenter.y = cellctr->y;
   Dpmta_CellCenter.z = cellctr->z;
   Dpmta_MaxCellLen = Max_CellLength();
   Dpmta_Resize = TRUE;

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
*  for the t3d version this routine calls routines to free up
*  all the dynamic memory buffers.
*
*/

int PMTAexit()
{

   /* free dynamic structures created by the multipole library */
   Slave_Delete();
   
   /* finally, free up the MIMD local send buffers */
   Delete_Local_Buffers();
   
   return(0);

}
