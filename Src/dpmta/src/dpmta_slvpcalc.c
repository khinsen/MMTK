/* 
*  dpmta_slvpcalc.c - routines to perform the single and double direct
*     particle interactions.
*
*  w. t. rankin
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*
*  originally, these were part of the dpmta_slvcalc.c (version 2.3) module
*  but were pulled out separately since the particle-particle computations
*  are independant of whatever multipole calculations we are running.
*
*/

static char rcsid[] = "$Id: dpmta_slvpcalc.c,v 2.14 1998/04/01 20:08:28 wrankin Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_slvpcalc.c,v $
 * Revision 2.14  1998/04/01 20:08:28  wrankin
 * added support for HILBERT ordering
 * general cleanup of code structure (more OOP if you will)
 *
 * Revision 2.13  1997/11/07 16:49:49  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 * Revision 2.12  1997/03/26 20:36:33  wrankin
 * rewrote Ilist and IIlist code to handle reallocation of array sizes
 *   that happen as a resuly of re-scaling a non-cubic unit cell.
 *
 * Revision 2.11  1996/11/18  19:29:38  wrankin
 * update to virial code to use multipole virial calculation
 * old virial code is turned on by defining OLDVIRIAL
 * direct test program now computes absolute virial tensor.
 *
 * Revision 2.10  1996/09/24  18:43:31  wrankin
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
 * Revision 2.9  1996/08/14  16:06:59  wrankin
 * Fixed LJ computations
 *
 * Revision 2.8  1996/08/09  15:31:08  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.7  1996/02/29  21:13:54  wrankin
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
 * Revision 2.6  1995/12/08  23:00:38  wrankin
 * preliminary release of DPMTA 2.3
 *   - added working Periodic Boundary Conditions (PDC) flag
 *   - added code for Virial Computation (not working yet)
 *   - general cleanup of some modules
 *
 * Revision 2.5  1995/11/29  22:29:36  wrankin
 * addition of periodic boundary code
 * beginning work on virial code.
 *
 * Revision 2.4  1995/10/01  21:46:21  wrankin
 * implemented LJ multipole calculations - compile time option only
 *
 * Revision 2.3  1995/07/26  01:52:57  wrankin
 * updated Makefiles
 * slight performance improvement for direct calculations
 * clean up test code example
 * fixed t3dlib for new multipole code
 *
 * Revision 2.2  1995/06/27  14:20:33  wrankin
 * Implementation of relative interaction list algorithm.
 * Interaction list generation code removed from master library.
 * Cell Table restructured as a 2-D table indexed by level/cell
 * All direct computations are now performed in a single pass.
 * General cleanup and removal of old files.
 *
 * Revision 2.1  1995/06/13  04:26:20  wrankin
 * No updates.
 * Sources placed under CVS.
 *
 * Revision 1.2  1995/04/09  22:21:06  wrankin
 * renaming of all global variables to prevent conflicts for
 * port of dpmta library to the T3D.
 *
 * Revision 1.1  1994/10/14  05:13:34  wrankin
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


/*
 * external prototyping
 */

#include "dpmta_slvmkil.h"
#include "dpmta_distmisc.h"

/*
 * prototyping of functions internal to module
 *
 */

void Cell_Calc_Self( int, int );
void Cell_Calc_DDirect( int, int, int, int, Vector * );
void Cell_Calc_SDirect( int, int, int, int, Vector * );

/*
 * globals used for LJ scaling, if needed
 */

#ifdef COMP_LJ
static double LJscale; 
#endif

/****************************************************************
*
*  this procedure will perform the double direct interactions
*  for each of the cells that the processor owns.  it cycles through 
*  each of the cells owned by the processor and for each cell, it
*  first processes the interactions between all the particles in a cell.
*  then  computes the direct interactions between the particles in that
*  cell and the particles in the cells listed in the cells direct
*  interaction list.
*
*  this procedure will perform the single direct interactions
*  for each of the cells that the processor owns.  it cycles through 
*  each of the cells owned by the processor and for each cell, it
*  computes the direct interactions between the particles in that
*  cell and the particles in the cells listed in the cells single
*  direct interaction list.
*
*/

void Slave_Direct_Calc()
{
   int i,j;         /* loop counters */
   int id;
   int level;
   int rcell;       /* remote cell id */
   int posn;
   int sep;         /* cell separation vector */
   int ovfl;        /* overflow boundary flag */


#ifdef COMP_LJ
   /* update LJ scaling factor */
   LJscale = 1.0 / Dpmta_MaxCellLen;
   LJscale = LJscale * LJscale * LJscale;
   LJscale *= LJscale;
#endif

   level = Dpmta_NumLevels - 1;

   /* cycle through all cells at bottom level */
   for ( i=Dpmta_Sindex[level]; i<=Dpmta_Eindex[level]; i++) {
      id = index2cell(i,level);
      
      /* determin cell position within parent cell */
      posn = id & 0x07;

      for ( j=0; j<Dpmta_Intlist[posn].dcnt; j++ ) {
         sep = Dpmta_Intlist[posn].dlist[j];
         if ( Cell2Cell(level,id,sep,&rcell,&ovfl) ) {

	    /* check for remote cell not allocated */
            if ( Dpmta_CellTbl[level][rcell] == NULL ) {
               fprintf(stderr,"ERROR: cell %d/%d not allocated\n",
		       level,rcell);
               exit(-1);
            } /* if Dpmta_CellTbl */

	    /*
	    *  check to see if we own this cell and the cell is located
	    *  within the boundary of the unit cell.  if so, then
	    *  we can perform double direct interactions.
	    *  we only do so when the remote cell is at a higher cell
	    *  index. this is to prevent doing the same pair of
	    *  particles twice.
	    */

            if ( Dpmta_CellTbl[level][rcell]->pid == Dpmta_Pid ) {
 	       if ( rcell > id ) {
                  Cell_Calc_DDirect(level,id,level,rcell,
				    &(Dpmta_Hlist[posn].dlist_vec[j]));
	       } /* if rcell and ovfl */
	    } /* if pid */

	    /*
	    *  otherwise, perform single-direct calculations
	    */

	    else {
	       Cell_Calc_SDirect(level,id,level,rcell,
				 &(Dpmta_Hlist[posn].dlist_vec[j]));
	    } /* else not pid */

	 } /* if Cell2cell */
      } /* for j */
   } /* for i */

   /*
   * in a futile attempt to improve accuracy, do not compute inter-cell
   * potentials until after all extra-cell calculations are done.
   */

   for ( i=Dpmta_Sindex[level]; i<=Dpmta_Eindex[level]; i++) {
      id = index2cell(i,level);
      Cell_Calc_Self(level,id);
   } /* for i */

} /* Slave_Direct_Calc() */


/****************************************************************
*
*  this procedure will perform the direct force calculations for
*  each particle within a cell interacting with all the other
*  particles within that same cell.
*
*/

void Cell_Calc_Self( int level, int cell) 
{
   int i,j;
   double dx,dy,dz;
   double wtq;
   double ir, ir2;
   double ir_q, ir3_q;
   double f_x, f_y, f_z;
   ParticlePtr cell_particles;
   PartInfoPtr cell_forces;

#ifdef COMP_LJ
   double wta, wtb;
   double ir6, ir6_a, ir8_a;
   double ir12_b, ir14_b;
   double pot, ir_lj;
   PartInfoPtr cell_f_lj;
#endif

   cell_particles = Dpmta_CellTbl[level][cell]->plist;
   cell_forces = Dpmta_CellTbl[level][cell]->mdata->flist;

#ifdef COMP_LJ
   cell_f_lj = Dpmta_CellTbl[level][cell]->mdata->f_lj;
#endif

   for ( i=0; i < (Dpmta_CellTbl[level][cell]->n)-1; i++) {
      for ( j=i+1; j < Dpmta_CellTbl[level][cell]->n; j++) {

         wtq = cell_particles[i].q * cell_particles[j].q;

         dx = cell_particles[j].p.x - cell_particles[i].p.x;
         dy = cell_particles[j].p.y - cell_particles[i].p.y;
         dz = cell_particles[j].p.z - cell_particles[i].p.z;

         ir2 = 1.0/(dx*dx + dy*dy + dz*dz);
         ir = sqrt(ir2);

         ir_q = wtq * ir;
         cell_forces[j].v += ir_q;
         cell_forces[i].v += ir_q;

         ir3_q = ir_q * ir2;
	 f_x = ir3_q * dx;
	 f_y = ir3_q * dy;
	 f_z = ir3_q * dz;
         cell_forces[j].f.x += f_x;
         cell_forces[i].f.x -= f_x;
         cell_forces[j].f.y += f_y;
         cell_forces[i].f.y -= f_y;
         cell_forces[j].f.z += f_z;
         cell_forces[i].f.z -= f_z;

#if defined VIRIAL || defined OLDVIRIAL
         Dpmta_Vpot += ir_q;
         Dpmta_Vf.x -= f_x * dx;
         Dpmta_Vf.y -= f_y * dy;
         Dpmta_Vf.z -= f_z * dz;
#endif



#ifdef COMP_LJ

	 /*
	  *  note that the 1/r^6 is always an attractive force
	  *  so the signs are reversed from the colomb above.
	  *  the 1/r^12 is repulsive, so the signs are the same.
	  */

         wta = cell_particles[i].a * cell_particles[j].a;
         wtb = cell_particles[i].b * cell_particles[j].b;

	 ir6 = ir2 * ir2 * ir2;
         ir6_a = ir6 * wta;
	 ir8_a = ir6_a * ir2;
	 ir12_b = ir6 * ir6 * wtb;
	 ir14_b = ir12_b * ir2;

	 pot = (ir12_b * LJscale) - ir6_a;
         cell_f_lj[j].v += pot;
         cell_f_lj[i].v += pot;

	 ir_lj = 6.0 * ( (2.0 * ir14_b * LJscale) - ir8_a);
	 f_x = ir_lj * dx;
	 f_y = ir_lj * dy;
	 f_z = ir_lj * dz;
         cell_f_lj[j].f.x += f_x;
         cell_f_lj[i].f.x -= f_x;
         cell_f_lj[j].f.y += f_y;
         cell_f_lj[i].f.y -= f_y;
         cell_f_lj[j].f.z += f_z;
         cell_f_lj[i].f.z -= f_z;

#if defined VIRIAL || defined OLDVIRIAL
         Dpmta_Vpot_LJ += pot;
         Dpmta_Vf_LJ.x -= f_x * dx;
         Dpmta_Vf_LJ.y -= f_y * dy;
         Dpmta_Vf_LJ.z -= f_z * dz;
#endif

#endif
	 
      } /* for j */
   } /* for i */
   
} /* Cell_Calc_Self */


/****************************************************************
*
*  Cell_Calc_DDirect() -
*
*  this procedure will perform the direct force calculations for
*  each particle within the target cell interacting with all the
*  other particles within the remote cell.  the contents of the
*  remote cell are updated.
*
*/

void Cell_Calc_DDirect(
   int tlevel,     /* target level */
   int tcell,      /* target cell */
   int rlevel,     /* remote level */
   int rcell,      /* remote cell */
   Vector *sep)    /* cell separation vector */

{
   int i,j;
   double dx,dy,dz;
   double wtq;
   double ir, ir2;
   double ir_q, ir3_q;
   double f_x, f_y, f_z;
   ParticlePtr tcell_parts, rcell_parts;
   PartInfoPtr tcell_forces, rcell_forces;

#ifdef COMP_LJ
   double wta, wtb;
   double ir6, ir6_a, ir8_a;
   double ir12_b, ir14_b;
   double pot, ir_lj;
   PartInfoPtr tcell_f_lj, rcell_f_lj;
#endif

   tcell_parts = Dpmta_CellTbl[tlevel][tcell]->plist;
   rcell_parts = Dpmta_CellTbl[rlevel][rcell]->plist;
   tcell_forces = Dpmta_CellTbl[tlevel][tcell]->mdata->flist;
   rcell_forces = Dpmta_CellTbl[rlevel][rcell]->mdata->flist;

#ifdef COMP_LJ
   tcell_f_lj = Dpmta_CellTbl[tlevel][tcell]->mdata->f_lj;
   rcell_f_lj = Dpmta_CellTbl[rlevel][rcell]->mdata->f_lj;
#endif

   for ( i=0; i < (Dpmta_CellTbl[tlevel][tcell]->n); i++) {
      for ( j=0; j < Dpmta_CellTbl[rlevel][rcell]->n; j++) {

         wtq = tcell_parts[i].q * rcell_parts[j].q;

         dx = rcell_parts[j].p.x - tcell_parts[i].p.x + sep->x;
         dy = rcell_parts[j].p.y - tcell_parts[i].p.y + sep->y;
         dz = rcell_parts[j].p.z - tcell_parts[i].p.z + sep->z;

         ir2 = 1.0/(dx*dx + dy*dy + dz*dz);
         ir = sqrt(ir2);

         ir_q = wtq * ir;
         tcell_forces[i].v += ir_q;
         rcell_forces[j].v += ir_q;

         ir3_q = ir_q * ir2;
	 f_x = ir3_q * dx;
	 f_y = ir3_q * dy;
	 f_z = ir3_q * dz;
         tcell_forces[i].f.x -= f_x;
         tcell_forces[i].f.y -= f_y;
         tcell_forces[i].f.z -= f_z;
         rcell_forces[j].f.x += f_x;
         rcell_forces[j].f.y += f_y;
         rcell_forces[j].f.z += f_z;

#if defined VIRIAL || defined OLDVIRIAL
         Dpmta_Vpot += ir_q;
         Dpmta_Vf.x -= f_x * dx;
         Dpmta_Vf.y -= f_y * dy;
         Dpmta_Vf.z -= f_z * dz;
#endif

#ifdef COMP_LJ

	 /*
	  *  note that the 1/r^6 is always an attractive force
	  *  so the signs are reversed from the colomb above.
	  *  the 1/r^12 is repulsive, so the signs are the same.
	  */

         wta = tcell_parts[i].a * rcell_parts[j].a;
         wtb = tcell_parts[i].b * rcell_parts[j].b;

	 ir6 = ir2 * ir2 * ir2;
         ir6_a = ir6 * wta;
	 ir8_a = ir6_a * ir2;
	 ir12_b = ir6 * ir6 * wtb;
	 ir14_b = ir12_b * ir2;

	 pot = (ir12_b * LJscale) - ir6_a;
         rcell_f_lj[j].v += pot;
         tcell_f_lj[i].v += pot;

	 ir_lj = 6.0 * ( (2.0 * ir14_b * LJscale) - ir8_a);
	 f_x = ir_lj * dx;
	 f_y = ir_lj * dy;
	 f_z = ir_lj * dz;
         tcell_f_lj[i].f.x -= f_x;
         tcell_f_lj[i].f.y -= f_y;
         tcell_f_lj[i].f.z -= f_z;
         rcell_f_lj[j].f.y += f_y;
         rcell_f_lj[j].f.x += f_x;
         rcell_f_lj[j].f.z += f_z;

#if defined VIRIAL || defined OLDVIRIAL
         Dpmta_Vpot_LJ += pot;
         Dpmta_Vf_LJ.x -= f_x * dx;
         Dpmta_Vf_LJ.y -= f_y * dy;
         Dpmta_Vf_LJ.z -= f_z * dz;
#endif

#endif

      } /* for j */
   } /* for i */
} /* Cell_Calc_DDirect */



/****************************************************************
*
*  Cell_Calc_SDirect( ) -
*
*  this procedure will perform the direct force calculations for each
*  particle within the target cell interacting with all the other
*  particles within the remote cell.  the contents of the remote cell
*  are not updated.
* */

void Cell_Calc_SDirect(
   int tlevel,     /* target level */
   int tcell,      /* target cell */
   int rlevel,     /* remote level */
   int rcell,      /* remote cell */
   Vector *sep)    /* cell separation vector */

{
   int i,j;
   double dx,dy,dz;
   double wtq;
   double ir, ir2;
   double ir_q, ir3_q;
   double f_x, f_y, f_z;
   ParticlePtr tcell_parts, rcell_parts;
   PartInfoPtr tcell_forces;

#ifdef COMP_LJ
   double wta, wtb;
   double ir6, ir6_a, ir8_a;
   double ir12_b, ir14_b;
   double pot, ir_lj;
   PartInfoPtr tcell_f_lj;
#endif

   tcell_parts = Dpmta_CellTbl[tlevel][tcell]->plist;
   rcell_parts = Dpmta_CellTbl[rlevel][rcell]->plist;
   tcell_forces = Dpmta_CellTbl[tlevel][tcell]->mdata->flist;

#ifdef COMP_LJ
   tcell_f_lj = Dpmta_CellTbl[tlevel][tcell]->mdata->f_lj;
#endif

   for ( i=0; i < (Dpmta_CellTbl[tlevel][tcell]->n); i++) {
      for ( j=0; j < Dpmta_CellTbl[rlevel][rcell]->n; j++) {

         wtq = tcell_parts[i].q * rcell_parts[j].q;

         dx = rcell_parts[j].p.x - tcell_parts[i].p.x + sep->x;
         dy = rcell_parts[j].p.y - tcell_parts[i].p.y + sep->y;
         dz = rcell_parts[j].p.z - tcell_parts[i].p.z + sep->z;

         ir2 = 1.0/(dx*dx + dy*dy + dz*dz);
	 ir = sqrt(ir2);

         ir_q = wtq * ir;
         tcell_forces[i].v += ir_q;

         ir3_q = ir_q * ir2;
	 f_x = ir3_q * dx;
	 f_y = ir3_q * dy;
	 f_z = ir3_q * dz;
         tcell_forces[i].f.x -= f_x;
         tcell_forces[i].f.y -= f_y;
         tcell_forces[i].f.z -= f_z;

#if defined VIRIAL || defined OLDVIRIAL
         Dpmta_Vpot += ir_q * 0.5;
         Dpmta_Vf.x -= f_x * dx * 0.5;
         Dpmta_Vf.y -= f_y * dy * 0.5;
         Dpmta_Vf.z -= f_z * dz * 0.5;
#endif

#ifdef COMP_LJ

	 /*
	  *  note that the 1/r^6 is always an attractive force
	  *  so the signs are reversed from the colomb above.
	  *  the 1/r^12 is repulsive, so the signs are the same.
	  */

         wta = tcell_parts[i].a * rcell_parts[j].a;
         wtb = tcell_parts[i].b * rcell_parts[j].b;

	 ir6 = ir2 * ir2 * ir2;
         ir6_a = ir6 * wta;
	 ir8_a = ir6_a * ir2;
	 ir12_b = ir6 * ir6 * wtb;
	 ir14_b = ir12_b * ir2;

	 pot = (ir12_b * LJscale) - ir6_a;
         tcell_f_lj[i].v += pot;

	 ir_lj = 6.0 * ( (2.0 * ir14_b * LJscale) - ir8_a);
	 f_x = ir_lj * dx;
	 f_y = ir_lj * dy;
	 f_z = ir_lj * dz;
         tcell_f_lj[i].f.x -= f_x;
         tcell_f_lj[i].f.y -= f_y;
         tcell_f_lj[i].f.z -= f_z;

#if defined VIRIAL || defined OLDVIRIAL
         Dpmta_Vpot_LJ += pot * 0.5;
         Dpmta_Vf_LJ.x -= f_x * dx * 0.5;
         Dpmta_Vf_LJ.y -= f_y * dy * 0.5;
         Dpmta_Vf_LJ.z -= f_z * dz * 0.5;
#endif
#endif

      } /* for j */
   } /* for i */
   
} /* Cell_Calc_SDirect() */


