/*
*  dpmta_serial - implement a serial version of DPMTA
*
*  w. t. rankin
*
*  Copyright (c) 1997 Duke University
*  All rights reserved
*
*  these routines merge the functionality of dpmta_slave and the
*  master interface routines into a single module for execution
*  on a serial workstation.
*
*  they are based upon the dpmta_t3dlib module veriosn 2.13
*
*/

static char rcsid[] = "$Id: dpmta_serial.c,v 2.17 1998/03/10 22:22:01 wrankin Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_serial.c,v $
 * Revision 2.17  1998/03/10 22:22:01  wrankin
 * folded start/cleanup functionality into dpmta_slvcompute
 *
 * Revision 2.16  1997/11/12 13:49:46  wrankin
 * updates to communications routines and general cleanup of code in prep
 *   for introducing load balancing functionality in hte near future
 *
 * Revision 2.15  1997/11/07 16:49:13  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 * Revision 2.14  1997/05/12 18:06:05  wrankin
 * added routines to clean up dynamic memory allocation when
 * PMTAinit() is called
 *
 * Revision 2.13  1997/05/07  21:27:52  chumphre
 * cubectr changed to cellctr and other minor changes related to the interface
 *
 * Revision 2.12  1997/05/07  19:31:46  chumphre
 * implement a uniform interface for rectangular and ||-piped cells
 *
 * Revision 2.11  1997/05/07  18:59:24  chumphre
 * major changes to implement ||-piped cells
 *
 * Revision 2.10  1997/04/30  20:40:13  wrankin
 * added -DPBC_WRAP flag to serial code
 *
 * Revision 2.9  1997/04/28  16:03:00  wrankin
 * removed silly comment
 *
 * Revision 2.8  1997/04/11  18:55:17  wrankin
 * fixed the sort process timing information
 *
 * Revision 2.7  1997/04/10  18:30:26  wrankin
 * fixed new sorting routine to correct bug in LJ code
 *
 * Revision 2.6  1997/04/10  17:26:17  wrankin
 * improved particle sorting algorithm from serial runs.
 * removed extraneous timing code frm within sort proc.
 *
 * Revision 2.5  1997/04/08  19:11:18  wrankin
 * Fixed syntax problems in serial and LJ code.
 *
 * Revision 2.4  1997/03/26  20:36:20  wrankin
 * rewrote Ilist and IIlist code to handle reallocation of array sizes
 *   that happen as a resuly of re-scaling a non-cubic unit cell.
 *
 * Revision 2.3  1997/03/12  20:15:45  wrankin
 * new more efficient serial sorting
 *
 * Revision 2.2  1997/02/26  20:43:27  wrankin
 * updated timing measurements and placed routines in single module
 *
 * Revision 2.1  1997/02/26  16:54:26  wrankin
 * added support for serial non-pvm compilation of DPMTA
 *
 *
 */


/* to handle lack of pvm settings */
#define MAXPROC 1


/* include files */

#include <stdlib.h>
#include "dpmta.h"
#include "dpmta_cell.h"
#include "dpmta_slvglobals.h"

/*
 * prototypes - lots of them here, because this is
 *   the top level routine and calls everything else.
 */

#include "dpmta_slvcompute.h"
#include "dpmta_slvscale.h"


#ifdef TIME
#include "dpmta_timer.h"
#endif


/****************************************************************
*
*  some prototypes
*/


/* prototypes for internal functions */
void Sort_Particles();
void Init_Local_Buffers();
void Delete_Local_Buffers();
void Return_Results( int, PmtaPartInfoPtr, PmtaPartInfoPtr );


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

   /*
    * if we are running a serial process, we do not return
    * any tids since there are none.
    */

   rtn_tids = (int *)NULL;

   Dpmta_Pid = 0;
   Dpmta_Nproc = 1;
   Dpmta_NumLevels = initdata->nlevels;
   Dpmta_FFT = initdata->fft;
   Dpmta_PBC = initdata->pbc;
   Dpmta_Mp = initdata->mp;
#ifdef COMP_LJ
   Dpmta_Mp_LJ = initdata->mp_lj;
#endif
   Dpmta_Theta = initdata->theta;
#ifdef MACROSCOPIC
   Dpmta_K = initdata->kterm;
#endif
   Dpmta_FftBlock = initdata->fftblock;
   Dpmta_CallingNum = 0;

   if (Dpmta_PBC != 0)
      Dpmta_DownPassStart = 1;
   else
      Dpmta_DownPassStart = 2;

   
#ifndef PIPED
   Dpmta_CellVector1.x = initdata->v1.x;
   Dpmta_CellVector2.y = initdata->v2.y;
   Dpmta_CellVector3.z = initdata->v3.z;
#else
   Dpmta_CellVector1.x = initdata->v1.x;
   Dpmta_CellVector1.y = initdata->v1.y;
   Dpmta_CellVector1.z = initdata->v1.z;
   Dpmta_CellVector2.x = initdata->v2.x;
   Dpmta_CellVector2.y = initdata->v2.y;
   Dpmta_CellVector2.z = initdata->v2.z;
   Dpmta_CellVector3.x = initdata->v3.x;
   Dpmta_CellVector3.y = initdata->v3.y;
   Dpmta_CellVector3.z = initdata->v3.z;
   Dpmta_CV1Mag = Vec_Mag(&Dpmta_CellVector1);
   Dpmta_CV2Mag = Vec_Mag(&Dpmta_CellVector2);
   Dpmta_CV3Mag = Vec_Mag(&Dpmta_CellVector3);
#endif
   Dpmta_CellCenter.x = initdata->cellctr.x;
   Dpmta_CellCenter.y = initdata->cellctr.y;
   Dpmta_CellCenter.z = initdata->cellctr.z;
   Dpmta_MaxCellLen = Max_CellLength();

   /*
    * force an initial resizing
    */

   Dpmta_Resize = TRUE;

   /*
    * initialize local data arrays
    */

   Init_Local_Buffers();

  
   /*
   *  now each slave needs to process this information and
   *  allocate teh appropriate data structures and constants 
   */

   Slave_Init();

   return(0);

} /* PMTAinit */


/****************************************************************
*
*  PMTAregister() - register application slave processes.
*
*  in the serial version, all processing done in this module
*  has been moved to PMTAinit().
*
*/

int PMTAregister()
{
   return(0);
}



/****************************************************************
*
*  PMTAforce - calculate colomb forces
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
   gettimeofday(&runstruct,0);
   times_arr[0] = (double)runstruct.tv_sec +
      ((double)runstruct.tv_usec/1000000);
#endif

#ifdef TIME
   /* begin timing sort and resize  */
   times(&startbuf);
#endif

   Slave_Start();

   /*
    * Load particles into cell table
   */

   Sort_Particles(nparts, particles);

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

   /*
   *  collect the results from each cell
   */

   Return_Results(nparts, results, results_lj);

   /*
   *  clean and free up the particle allocations and other data
   *  structures that are re-allocated on each time step.
   */

   Slave_Cleanup();

#ifdef TIME
   gettimeofday(&runstruct,0);
   times_arr[1] = (double)runstruct.tv_sec +
      ((double)runstruct.tv_usec/1000000);
   times_arr[1] -= times_arr[0];

   Print_Times(times_arr);
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

int PMTAresize( PmtaVector *v1, PmtaVector *v2, PmtaVector *v3,
                PmtaVector *cellctr)   /* center of simulation cell */

{
#ifndef PIPED
   Dpmta_CellVector1.x = v1->x;
   Dpmta_CellVector2.y = v2->y;
   Dpmta_CellVector3.z = v3->z;
#else
   Dpmta_CellVector1.x = v1->x;
   Dpmta_CellVector1.y = v1->y;
   Dpmta_CellVector1.z = v1->z;
   Dpmta_CellVector2.x = v2->x;
   Dpmta_CellVector2.y = v2->y;
   Dpmta_CellVector2.z = v2->z;
   Dpmta_CellVector3.x = v3->x;
   Dpmta_CellVector3.y = v3->y;
   Dpmta_CellVector3.z = v3->z;
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
   *vp = Dpmta_Vpot;
   vf->x = Dpmta_Vf.x;
   vf->y = Dpmta_Vf.y;
   vf->z = Dpmta_Vf.z;

#ifdef COMP_LJ
   *vp_lj = Dpmta_Vpot_LJ;
   vf_lj->x = Dpmta_Vf_LJ.x;
   vf_lj->y = Dpmta_Vf_LJ.y;
   vf_lj->z = Dpmta_Vf_LJ.z;
#endif
   return(0);
#else
   return(-1);
#endif

} /* PMTAvirial */


/****************************************************************
*
*  PMTAexit() - shuts down slave processes
*
*  for the serial version this routine needs to free up
*  all the dynamic data structures created by DPMTA
*
*/

int PMTAexit()
{

   /* free dynamic structures  */
   Slave_Delete();
   
   /* free up local buffers (serial version specific) */
   Delete_Local_Buffers();
   
   return(0);

}


/****************************************************************
 *
 *  beyond here are support functions not externally referenced
 *
 */

/*
*  static buffers for sorting particle data.
*  Make them global so that we have them all
*  in one place in case we need to free() them later
*/

static int         *SendPartCnt;
static int         SendPartCntSz;
static int         *SendCellId;
static int         SendCellIdSz;

static ParticlePtr *CellPart;
static PartInfoPtr *CellInfo;
#ifdef COMP_LJ
static PartInfoPtr *CellInfoLJ;
#endif 
static int         **CellPartId;
static int         CellPartSz;


/****************************************************************
*
* Init_Local_Buffers() - initializes local datastructures
*
*/

void Init_Local_Buffers()
{

   SendPartCnt = NULL;
   SendCellId = NULL;
   CellPart = NULL;
   CellInfo = NULL;
   CellPartId = NULL;
#ifdef COMP_LJ
   CellInfoLJ = NULL;
#endif
   
   SendPartCntSz = 0;
   SendCellIdSz = 0;
   CellPartSz = 0;
   
}
   

/****************************************************************
*
*  Sort_Particles() - sort and distribute particle lists
*     to the slave processes.
*
*/

void Sort_Particles(
   int num_parts,               /* number of particles */
   PmtaParticlePtr part_table)   /* particle table */

{
   int    i,j;                  /* loop counters */
   int    level;                /* points to bottom level on cell table */
   int    num_cells;            /* number of cells */
   int    cell_x, cell_y, cell_z; /* cartesian cell offsets */
   int    cell_id;              /* cell index */
   int    nparts;               /* particle count */
   int    wrap_mask;            /* mask to wrap cells around boundary */
   double cell_edge;            /* number of cells per edge */
#ifdef PIPED
   double v2xv3dotp, v3xv1dotp, v1xv2dotp; /* Dot products */
   double shiftc;               /* Proportional shift of center */
   PmtaVector v1xv2, v2xv3, v3xv1;
   PmtaVector tmpc;             /* Temporary Parallelepiped center */
#else
   int    cell_mask;            /* cell id mask */
#endif

   ParticlePtr part_tmp;        /* temporary pointer to list of particles */
   PartInfoPtr flist_tmp;       /* temporary pointer to force vector */
#ifdef COMP_LJ
   PartInfoPtr f_lj_tmp;        /* temporary pointer to force vector */
#endif
   int *idlist_tmp;             /* temporary pointer to id vector */

   ParticlePtr *part_p_tmp;     /* temporary pointer to list of particles */
   PartInfoPtr *flist_p_tmp;    /* temporary pointer to force vector */
#ifdef COMP_LJ
   PartInfoPtr *f_lj_p_tmp;     /* temporary pointer to force vector */
#endif
   int **idlist_p_tmp;          /* temporary pointer to id vector */

   int *cellid_tmp;             /* temp cell ID pointer */
   CellPtrPtr cell_tmp;         /* temporary cell pointer */
   PmtaParticlePtr part_in;     /* input part. array ptr from application */

   level = Dpmta_NumLevels - 1;
   num_cells = Dpmta_Power8[level];

   /*
   *  check to see if we have enough buffer space to hold the 
   *  particle and cell buffers and allocate extra if needed
   */

   if ( num_parts > SendCellIdSz ) {
      SendCellId = (int *)realloc((void *)SendCellId,num_parts*sizeof(int));
      SendCellIdSz = num_parts;
   }

   /*
   *  allocate cell indices -  note that this should never have
   *  to be reallocated until we get into some funky reallocation
   *  scheme,  but is always pays to be prepared.
   */

   if ( num_cells > SendPartCntSz ) {
      SendPartCnt = (int *)realloc((void *)SendPartCnt,num_cells*sizeof(int));
      SendPartCntSz = num_cells;
   }

   if ( num_cells > CellPartSz ) {
      CellPart = (ParticlePtr *)realloc((void *)CellPart,
				num_cells*sizeof(ParticlePtr));
      CellInfo = (PartInfoPtr *)realloc((void *)CellInfo,
				num_cells*sizeof(PartInfoPtr));
      CellPartId = (int **)realloc((void *)CellPartId,
				num_cells*sizeof(int *));
#ifdef COMP_LJ
      CellInfoLJ = (PartInfoPtr *)realloc((void *)CellInfoLJ,
				num_cells*sizeof(PartInfoPtr));
#endif
      CellPartSz = num_cells;
   }


   /* initailize counters */
   for (i=0; i<num_cells; i++) {
      SendPartCnt[i]=0;
   }

   /* compute cell numbers and size */

   cell_edge = (double)(0x1 << (level));
   wrap_mask = (0x1 << level) - 1;

   /* cycle through all particles */
   for (i=0; i<num_parts; i++) {

#ifndef PIPED
      /* compute integer cell coordinates */
      cell_x = (int)(((part_table[i].p.x - Dpmta_CellCenter.x) /
		      Dpmta_CellVector1.x + 0.5) * cell_edge);
      cell_y = (int)(((part_table[i].p.y - Dpmta_CellCenter.y) /
		      Dpmta_CellVector2.y + 0.5) * cell_edge);
      cell_z = (int)(((part_table[i].p.z - Dpmta_CellCenter.z) /
		      Dpmta_CellVector3.z + 0.5) * cell_edge);

      /*
       * when running PBCs, it is a good thing if particles
       * do not get wrapped, but rather if they wander outside
       * the box, they are kept in the outter-most cell nearest
       * them.  otherwise it would look (to the simulation) that
       * a buch of particles were disappearing off of one side of
       * simulation space and appearing "magically" on the other.
       * this could play hell with constant energy sims.
       *
       * if the user wants wrapping, then simply define the PBC_WRAP
       * flag.
       */

#ifdef PBC_WRAP
      cell_x &= wrap_mask;
      cell_y &= wrap_mask;
      cell_z &= wrap_mask;
#else
      if ( cell_x < 0 )
	 cell_x = 0;
      if ( cell_x > wrap_mask )      
	 cell_x = wrap_mask;
      if ( cell_y < 0 )      
	 cell_y = 0;
      if ( cell_y > wrap_mask )
	 cell_y = wrap_mask;
      if ( cell_z < 0 )
	 cell_z = 0;
      if ( cell_z > wrap_mask )
	 cell_z = wrap_mask;
#endif


      /*
      *  build cell id from absolute coordinates
      *  there may be a more efficient way to do this
      */

      cell_mask = 0x1;
      cell_id = 0;
      cell_y = cell_y << 1;
      cell_z = cell_z << 2;
      for (j=0; j<level; j++) {
         cell_id |= cell_x & cell_mask;
         cell_mask = cell_mask << 1;
         cell_x = cell_x << 2;
         cell_id |= cell_y & cell_mask;
         cell_mask = cell_mask << 1;
         cell_y = cell_y << 2;
         cell_id |= cell_z & cell_mask;
         cell_mask = cell_mask << 1;
         cell_z = cell_z << 2;
      } /* for j */

#else
      tmpc.x = Dpmta_CellCenter.x;
      tmpc.y = Dpmta_CellCenter.y;
      tmpc.z = Dpmta_CellCenter.z;
      cell_id = 0;
      shiftc = 2.0;

      v2xv3.x = Dpmta_CellVector2.y*Dpmta_CellVector3.z - 
		Dpmta_CellVector2.z*Dpmta_CellVector3.y;
      v2xv3.y = Dpmta_CellVector2.z*Dpmta_CellVector3.x - 
		Dpmta_CellVector2.x*Dpmta_CellVector3.z;
      v2xv3.z = Dpmta_CellVector2.x*Dpmta_CellVector3.y - 
		Dpmta_CellVector2.y*Dpmta_CellVector3.x;

      v3xv1.x = Dpmta_CellVector3.y*Dpmta_CellVector1.z - 
		Dpmta_CellVector3.z*Dpmta_CellVector1.y;
      v3xv1.y = Dpmta_CellVector3.z*Dpmta_CellVector1.x - 
		Dpmta_CellVector3.x*Dpmta_CellVector1.z;
      v3xv1.z = Dpmta_CellVector3.x*Dpmta_CellVector1.y - 
		Dpmta_CellVector3.y*Dpmta_CellVector1.x;

      v1xv2.x = Dpmta_CellVector1.y*Dpmta_CellVector2.z - 
		Dpmta_CellVector1.z*Dpmta_CellVector2.y;
      v1xv2.y = Dpmta_CellVector1.z*Dpmta_CellVector2.x - 
		Dpmta_CellVector1.x*Dpmta_CellVector2.z;
      v1xv2.z = Dpmta_CellVector1.x*Dpmta_CellVector2.y - 
		Dpmta_CellVector1.y*Dpmta_CellVector2.x;

      for(j=0;j<level;j++)
      {
        cell_id = cell_id << 3;

        v2xv3dotp = ((part_table[i].p.x - tmpc.x)*v2xv3.x) +
                 ((part_table[i].p.y - tmpc.y)*v2xv3.y) +
                 ((part_table[i].p.z - tmpc.z)*v2xv3.z);
        v3xv1dotp = ((part_table[i].p.x - tmpc.x)*v3xv1.x) +
                 ((part_table[i].p.y - tmpc.y)*v3xv1.y) +
                 ((part_table[i].p.z - tmpc.z)*v3xv1.z);
        v1xv2dotp = ((part_table[i].p.x - tmpc.x)*v1xv2.x) +
                 ((part_table[i].p.y - tmpc.y)*v1xv2.y) +
                 ((part_table[i].p.z - tmpc.z)*v1xv2.z);

        if(v2xv3dotp<0) cell_x = 0;
        else cell_x = 1;

        if(v3xv1dotp<0) cell_y = 0;
        else cell_y = 1;

        if(v1xv2dotp<0) cell_z = 0;
        else cell_z = 1;

        cell_y = cell_y << 1;
        cell_z = cell_z << 2;

        cell_id |= (cell_x | cell_y | cell_z);

        shiftc *= 2.0;

        if(cell_x)
        {
          tmpc.x += Dpmta_CellVector1.x/shiftc;
          tmpc.y += Dpmta_CellVector1.y/shiftc;
          tmpc.z += Dpmta_CellVector1.z/shiftc;
        }
        else
        {
          tmpc.x -= Dpmta_CellVector1.x/shiftc;
          tmpc.y -= Dpmta_CellVector1.y/shiftc;
          tmpc.z -= Dpmta_CellVector1.z/shiftc;
        }

        if(cell_y)
        {
          tmpc.x += Dpmta_CellVector2.x/shiftc;
          tmpc.y += Dpmta_CellVector2.y/shiftc;
          tmpc.z += Dpmta_CellVector2.z/shiftc;
        }
        else
        {
          tmpc.x -= Dpmta_CellVector2.x/shiftc;
          tmpc.y -= Dpmta_CellVector2.y/shiftc;
          tmpc.z -= Dpmta_CellVector2.z/shiftc;
        }

        if(cell_z)
        {
          tmpc.x += Dpmta_CellVector3.x/shiftc;
          tmpc.y += Dpmta_CellVector3.y/shiftc;
          tmpc.z += Dpmta_CellVector3.z/shiftc;
        }
        else
        {
          tmpc.x -= Dpmta_CellVector3.x/shiftc;
          tmpc.y -= Dpmta_CellVector3.y/shiftc;
          tmpc.z -= Dpmta_CellVector3.z/shiftc;
        }
      }
#endif

      /* store the cell index for the particle */
      SendCellId[i] = cell_id;

      /* increment the partical counter for that cell */
      SendPartCnt[cell_id] += 1;

   } /* for i */


   cell_tmp = Dpmta_CellTbl[level];
   part_p_tmp = CellPart;
   flist_p_tmp =  CellInfo;
   idlist_p_tmp = CellPartId;
#ifdef COMP_LJ
   f_lj_p_tmp = CellInfoLJ;
#endif

   for (i=0; i<num_cells; i++) {

      nparts = SendPartCnt[i];
      (*cell_tmp)->n = nparts;

      /*
      *  here is where we need to check psize against the number 
      *  of particles in the list and realloc is it is not large
      *  enough.  since we are receiving all of out particles from
      *  a single source (the master) we don't have to worry about
      *  preserving any existing data.
      */

      if ( nparts > (*cell_tmp)->psize ) {

         part_tmp = (*cell_tmp)->plist;
         part_tmp = (ParticlePtr)realloc(part_tmp,
                       nparts*sizeof(Particle));
         (*cell_tmp)->plist = part_tmp;

         idlist_tmp = (*cell_tmp)->mdata->part_id;
         idlist_tmp = (int *)realloc(idlist_tmp,nparts*sizeof(int));
         (*cell_tmp)->mdata->part_id = idlist_tmp;

         flist_tmp = (*cell_tmp)->mdata->flist;
         flist_tmp = (PartInfoPtr)realloc(flist_tmp,
                       nparts*sizeof(PartInfo));
         (*cell_tmp)->mdata->flist = flist_tmp;

#ifdef COMP_LJ
         f_lj_tmp = (*cell_tmp)->mdata->f_lj;
         f_lj_tmp = (PartInfoPtr)realloc(f_lj_tmp,
                       nparts*sizeof(PartInfo));
         (*cell_tmp)->mdata->f_lj = f_lj_tmp;
#endif 

         (*cell_tmp)->psize = nparts;
      } /* if nparts */

      /*
       * load pointers that reference the particle and force
       * arrays within each cell.  we will use these to keep
       * track of the next particle position within each
       * cell list
       */

      (*part_p_tmp) = (*cell_tmp)->plist;
      (*flist_p_tmp) = (*cell_tmp)->mdata->flist;
      (*idlist_p_tmp) = (*cell_tmp)->mdata->part_id;
#ifdef COMP_LJ
      (*f_lj_p_tmp) = (*cell_tmp)->mdata->f_lj;
#endif

      part_p_tmp++;
      flist_p_tmp++;
      idlist_p_tmp++;
#ifdef COMP_LJ
      f_lj_p_tmp++;
#endif

      cell_tmp++;

   } /* for i */

   /*
    * cycle through all the particles and place each
    * one at a time into the appropriate cell structure.
    */

   cellid_tmp = SendCellId;
   part_in = part_table;

   for ( i=0; i<num_parts; i++ ) {

      cell_id = *cellid_tmp;

      part_tmp = CellPart[cell_id];
      flist_tmp = CellInfo[cell_id];
      idlist_tmp = CellPartId[cell_id];
#ifdef COMP_LJ
      f_lj_tmp = CellInfoLJ[cell_id];
#endif

      part_tmp->p.x = part_in->p.x;
      part_tmp->p.y = part_in->p.y;
      part_tmp->p.z = part_in->p.z;
      part_tmp->q = part_in->q;
#ifdef COMP_LJ
      part_tmp->a = part_in->a;
      part_tmp->b = part_in->b;
#endif

      flist_tmp->f.x = 0.0;
      flist_tmp->f.y = 0.0;
      flist_tmp->f.z = 0.0;
      flist_tmp->v = 0.0;

      *idlist_tmp = i;

#ifdef COMP_LJ
      f_lj_tmp->f.x = 0.0;
      f_lj_tmp->f.y = 0.0;
      f_lj_tmp->f.z = 0.0;
      f_lj_tmp->v = 0.0;
#endif

      part_in++;
      cellid_tmp++;

      CellPart[cell_id]++;
      CellInfo[cell_id]++;
      CellPartId[cell_id]++;
#ifdef COMP_LJ
      CellInfoLJ[cell_id]++;
#endif

   } /* for i */

} /* Sort_Particles */


/****************************************************************
*
*   Return_Results() - put result data from cell table back into
*     source array
*/

void Return_Results(
   int nparts,                   /* number of particles to generate */
   PmtaPartInfoPtr results,      /* array of force results */
   PmtaPartInfoPtr results_lj)   /* array of LJ force results */
{
   int i,j;                      /* loop counters */
   int partid;                   /* sorting counters */
   int num_cells, num_parts;     /* sorting counters */
   int level;

   int *idlist_tmp;              /* temporary pointer to id vector */
   PartInfoPtr flist_tmp;        /* temporary pointer to force vector */
#ifdef COMP_LJ
   PartInfoPtr f_lj_tmp;         /* temporary pointer to force vector */
#endif

   /*
   *  cycle through cell table and store the resulting force
   *  in the result array indexed by the particles id.
   */

   level = Dpmta_NumLevels - 1;
   num_cells = Dpmta_Power8[level];

   for (i=0; i<num_cells; i++) {
      num_parts = Dpmta_CellTbl[level][i]->n;
      idlist_tmp = Dpmta_CellTbl[level][i]->mdata->part_id;
      flist_tmp = Dpmta_CellTbl[level][i]->mdata->flist;
#ifdef COMP_LJ
      f_lj_tmp = Dpmta_CellTbl[level][i]->mdata->f_lj;
#endif

      for (j=0; j<num_parts; j++) {
	 partid = idlist_tmp[j];
	 results[partid].f.x = flist_tmp[j].f.x;
	 results[partid].f.y = flist_tmp[j].f.y;
	 results[partid].f.z = flist_tmp[j].f.z;
	 results[partid].v = flist_tmp[j].v;
#ifdef COMP_LJ
	 results_lj[partid].f.x = f_lj_tmp[j].f.x;
	 results_lj[partid].f.y = f_lj_tmp[j].f.y;
	 results_lj[partid].f.z = f_lj_tmp[j].f.z;
	 results_lj[partid].v = f_lj_tmp[j].v;
#endif

	 /* need to add rescaling code here */

      } /* for j */
   } /* for i */
} /* Return_Results */


/****************************************************************
*
* Delete_Local_Buffers() - free up locally allocated dynamic data
*   structures.
*
*/

void Delete_Local_Buffers()
{

   free( SendPartCnt );
   free( SendCellId );
   free( CellPart );
   free( CellInfo );
   free( CellPartId );
#ifdef COMP_LJ
   free( CellInfoLJ );
#endif

} /* Delete_Local_Buffers */
