/* 
*  dpmta_slvscale.c - misc routines to perform particle and result
*     rescaling.  also routines to support dynamic cell resizing.
*
*  w. t. rankin
*
*  Copyright (c) 1996 Duke University
*  All rights reserved
*
*/

static char rcsid[] = "$Id: dpmta_slvscale.c,v 2.10 1998/04/01 20:08:29 wrankin Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_slvscale.c,v $
 * Revision 2.10  1998/04/01 20:08:29  wrankin
 * added support for HILBERT ordering
 * general cleanup of code structure (more OOP if you will)
 *
 * Revision 2.9  1997/11/07 16:49:52  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 * Revision 2.8  1997/05/09 18:32:30  chumphre
 * Implementation of macroscopic assemblies for ||-pipeds
 *
 * Revision 2.7  1997/05/07  18:59:44  chumphre
 * major changes to implement ||-piped cells
 *
 * Revision 2.6  1997/03/04  19:18:14  wrankin
 * updates to timing codes
 *
 * Revision 2.5  1997/01/13  22:05:18  wrankin
 * added code to perform parallel macroscopic expansion computation
 *
 * Revision 2.4  1996/11/18  19:29:39  wrankin
 * update to virial code to use multipole virial calculation
 * old virial code is turned on by defining OLDVIRIAL
 * direct test program now computes absolute virial tensor.
 *
 * Revision 2.3  1996/09/24  18:43:44  wrankin
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
 * Revision 2.2  1996/08/09  15:31:10  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.1  1996/02/29  21:13:56  wrankin
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
 *
 */

/* include files */
#include <math.h>
#include "dpmta_cell.h"
#include "dpmta_slvglobals.h"

/*
 * external prototypes
 */

#include "dpmta_slvmkcell.h"
#include "dpmta_slvmkhl.h"
#include "dpmta_distmisc.h"
#ifdef MACROSCOPIC
#include "dpmta_slvmacro.h"
#endif
#ifdef TIME
#include "dpmta_timer.h"
#endif

/****************************************************************
 *
 * Vec_Mag - Calculate the magnitude of a vector
 *
 */

double Vec_Mag(Vector *vec1)
{
  double mag1;

  mag1 = (vec1->x * vec1->x) + (vec1->y * vec1->y) + (vec1->z * vec1->z);
  mag1 = sqrt(mag1);

  return(mag1);

}

/****************************************************************
 *
 *  Max_CellLength - Calculate maximum length of cell edge
 *
 */

double Max_CellLength()
{
  double maxlen;

#ifdef PIPED
  double l1, l2, l3;

  l1 = (Dpmta_CellVector1.x * Dpmta_CellVector1.x) + 
       (Dpmta_CellVector1.y * Dpmta_CellVector1.y) + 
       (Dpmta_CellVector1.z * Dpmta_CellVector1.z);
  l2 = (Dpmta_CellVector2.x * Dpmta_CellVector2.x) + 
       (Dpmta_CellVector2.y * Dpmta_CellVector2.y) + 
       (Dpmta_CellVector2.z * Dpmta_CellVector2.z);
  l3 = (Dpmta_CellVector3.x * Dpmta_CellVector3.x) + 
       (Dpmta_CellVector3.y * Dpmta_CellVector3.y) + 
       (Dpmta_CellVector3.z * Dpmta_CellVector3.z);

  if(l1 > l2)
  {
    if(l1 > l3) maxlen = sqrt(l1);
    else maxlen = sqrt(l3);
  }
  else
  {
    if(l2 > l3) maxlen = sqrt(l2);
    else maxlen = sqrt(l3);
  }

#else

  if ( Dpmta_CellVector1.x > Dpmta_CellVector2.y )
     if ( Dpmta_CellVector1.x > Dpmta_CellVector3.z )
	maxlen = Dpmta_CellVector1.x;
     else
	maxlen = Dpmta_CellVector3.z;
  else
     if ( Dpmta_CellVector2.y > Dpmta_CellVector3.z )
	maxlen = Dpmta_CellVector2.y;
     else 
	maxlen = Dpmta_CellVector3.z;

#endif

  return maxlen;
}

/****************************************************************
 *
 *  Rescale_Particles - scale particle data to fit in a unit simulation
 *    cell.  all particles positions are then converted to relative
 *    positions w.r.t. the individual cell center.
 *
 */

void Rescale_Particles()
{

  int i,j;
  int id;
  int level;
  int num_parts;
  Vector *part;
  Vector *center;

  level = Dpmta_NumLevels - 1;

  for ( i=Dpmta_Sindex[level]; i<=Dpmta_Eindex[level]; i++ ) {
     id = index2cell(i,level);
     num_parts = Dpmta_CellTbl[level][id]->n;
     for ( j=0; j<num_parts; j++ ) {
	part = &(Dpmta_CellTbl[level][id]->plist[j].p);
        center = &(Dpmta_CellTbl[level][id]->p);
	part->x = (part->x - Dpmta_CellCenter.x)/Dpmta_MaxCellLen + 0.5;
        part->x -= center->x;
	part->y = (part->y - Dpmta_CellCenter.y)/Dpmta_MaxCellLen + 0.5;
        part->y -= center->y;
	part->z = (part->z - Dpmta_CellCenter.z)/Dpmta_MaxCellLen + 0.5;
        part->z -= center->z;
     } /* for j */
  } /* for i */

} /* Rescale_Particle */

      
/****************************************************************
 *
 *  Rescale_Results - rescale force data from a unit simulation
 *    cell back to the original cell size.
 *
 */

void Rescale_Results()
{

  int i,j;
  int id;
  int level;
  int num_parts;
  double pscale, fscale;
#ifdef COMP_LJ
  double pscale_lj, fscale_lj;
#endif

  level = Dpmta_NumLevels - 1;

  pscale = 1.0 / Dpmta_MaxCellLen;
  fscale = pscale * pscale; 

#ifdef COMP_LJ
  pscale_lj = fscale * fscale * fscale; 
  fscale_lj = pscale_lj * pscale; 
#endif

  for ( i=Dpmta_Sindex[level]; i<=Dpmta_Eindex[level]; i++ ) {
     id = index2cell(i,level);
     num_parts = Dpmta_CellTbl[level][id]->n;
     for ( j=0; j<num_parts; j++ ) {

        Dpmta_CellTbl[level][id]->mdata->flist[j].v *= pscale;
        Dpmta_CellTbl[level][id]->mdata->flist[j].f.x *= fscale;
        Dpmta_CellTbl[level][id]->mdata->flist[j].f.y *= fscale;
        Dpmta_CellTbl[level][id]->mdata->flist[j].f.z *= fscale;

#ifdef COMP_LJ
        Dpmta_CellTbl[level][id]->mdata->f_lj[j].v *= pscale_lj;
        Dpmta_CellTbl[level][id]->mdata->f_lj[j].f.x *= fscale_lj;
        Dpmta_CellTbl[level][id]->mdata->f_lj[j].f.y *= fscale_lj;
        Dpmta_CellTbl[level][id]->mdata->f_lj[j].f.z *= fscale_lj;
#endif


     } /* for j */
  } /* for i */

#if defined VIRIAL || defined OLDVIRIAL
  Dpmta_Vpot *= pscale;
  Dpmta_Vf.x *= pscale;
  Dpmta_Vf.y *= pscale;
  Dpmta_Vf.z *= pscale;
#ifdef COMP_LJ
  Dpmta_Vpot_LJ *= pscale_lj;
  Dpmta_Vf_LJ.x *= pscale_lj;
  Dpmta_Vf_LJ.y *= pscale_lj;
  Dpmta_Vf_LJ.z *= pscale_lj;
#endif
#endif

} /* Rescale_Results */

      
/****************************************************************
 *
 *  this procedure will perform any processing needed as a result
 *  of cell resizing.  it will also be called prior to the first
 *  iteration.
 *
 */

void MultipoleResize()
{

#ifdef TIME
   /* Time macroscopic */
   times(&startbuf);
#endif

#ifdef MACROSCOPIC
   if (Dpmta_PBC) {
      MacroPreComp(Dpmta_CellVector1, Dpmta_CellVector2,
          Dpmta_CellVector3, Dpmta_MaxCellLen);
   }
#endif

#ifdef TIME
   /* end time for macroscopic calculations */
   times(&endbuf);
   times_arr[13] = (double)(endbuf.tms_utime - startbuf.tms_utime) / (double)CLK_TCK;
#endif

   /* compute and store relative separation vectors */
   Make_RelVec(Dpmta_NumLevels-1);

   /* reset the centers of all the cells */
   Make_Cell_Centers();

} /* MultipoleResize */


