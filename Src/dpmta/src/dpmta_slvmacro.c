/* 
*  dpmta_slvmacro.c - routines to compute the local expansion representataion
*    of Chris Lamberts macro assemblies of particles.
*
*  w. t. rankin
*
*  Copyright (c) 1996 Duke University
*  All rights reserved
*
*  The procedures within this module attempt to implement the "Macro
*  Assembley" concept as described in the document:
*
*     "A Multipole-Based Algorithm for Efficient Calculation of
*      Forces and Potentials in Macroscopic Periodic Assemblies
*      of Particles"
*     Christophe Lambert
*     Duke University Dept. of Comp. Sci.
*     TR 95-001a
*     Revised October 1995.
*
*/

static char rcsid[] = "$Id: dpmta_slvmacro.c,v 2.13 1997/11/12 13:47:55 wrankin Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_slvmacro.c,v $
 * Revision 2.13  1997/11/12 13:47:55  wrankin
 * results are now collected at 'MastPid' instead of automatically going to
 *  pid-0.  we should probably just have MastPid set using getslvpid().
 *
 * Revision 2.12  1997/11/07 16:49:29  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 * Revision 2.11  1997/05/09 21:04:28  wrankin
 * added procedures for cleanup of dynamic data structures
 *
 * Revision 2.10  1997/05/09  18:32:27  chumphre
 * Implementation of macroscopic assemblies for ||-pipeds
 *
 * Revision 2.9  1997/03/04  19:18:08  wrankin
 * updates to timing codes
 *
 * Revision 2.8  1997/01/13  22:05:16  wrankin
 * added code to perform parallel macroscopic expansion computation
 *
 * Revision 2.7  1996/10/28  23:01:58  wrankin
 * additions to test routines to provide processing paramenters as
 *   part of output file.
 * dpmta_direct will now perform macroscopic computations, reading
 *   in processing parameters from particle position file.
 *
 * Revision 2.6  1996/09/24  18:42:45  wrankin
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
 * Revision 2.5  1996/09/10  18:07:03  wrankin
 * fixed macroscopic processing for non-cubic cells.
 * fixed passing of local expansions for non-2^n processors
 *
 * Revision 2.4  1996/08/20  17:12:45  wrankin
 * Enhancement to allow any number of processes to be utilized.
 * Formally the number of precesses was limited to a power of 2.
 *
 * Revision 2.3  1996/08/09  15:30:54  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.2  1996/02/29  21:13:40  wrankin
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
 * Revision 2.1  1996/02/12  15:21:07  wrankin
 * Added Macroscopic Assemblies code.
 *
 *
*/

/* include files */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "dpmta_cell.h"
#include "dpmta_slvmkil.h"
#include "dpmta_slvscale.h"

#ifdef PMACRO
#include "pvm3.h"
#include "dpmta_pvm.h"
#endif


/*
 * local prototypes for functions used within this module
 */

void MacroUpward( Mtype );
void MacroDownward( Mtype );
void MacroDump( Mtype *, Mtype );
#ifdef PMACRO
void MacroSend();
void MacroReceive();
#endif


/* global arrays and constants */
/* declared static to limit the scope */

static Mtype  *MacroMpole;      /* array of higher level multipoles */
static Mtype  *MacroM2M;        /* array of macro M2M xfer functions */
static Mtype  *MacroM2L;        /* array of macro M2L xfer functions */
static int    Kterm;            /* number of levels in the macroscopic */
static int    Mp;               /* number of terms in the multipole */
static int    Fft;              /* fft flag */
static Real   Theta;            /* multipole separation criteria */
static int    MyPid;            /* my processor id (0 for serial case) */
static int    MastPid;          /* proc id of whomever owns Cell [0,0] */

#ifdef PMACRO
static Mtype  MacroTemp;        /* temporary macro expansions */
static int    MpeSize;          /* size of multipole expansion (in doubles) */
static int    Nprocs;           /* number of slave processes */
static int    MastTid;          /* Pvm Tid of whomever owns Cell [0,0] */
static int    PreCompFlag;      /* flag to indicate rescale */
#endif

/*****************************************************************
 *
 *  MacroInit() - initialize data structures used in macroscopic 
 *    processing.  This should be called once by every process.
 *
 */

void MacroInit( int level, int mp, int fft, Real theta,
	   int mypid, int nprocs, int *tids, int mastpid )
{

   int i;                /* loop variables */

   /* set global variables */

   Kterm = level;
   Mp = mp;
   Fft = fft;
   Theta = theta;
   MyPid = mypid;
   MastPid = mastpid;

#ifdef PMACRO
   Nprocs = nprocs;
   MastTid = tids[mastpid];
   MpeSize = Csize(mp);
   PreCompFlag = FALSE;
#endif

   /* if we specify no levels, return */
   if ( Kterm == 0 ) {
      return;
   }

   /* allocate storage for everything */
   /* put the malloc checks in later */

   MacroMpole = (Mtype *)malloc(level*sizeof(Mtype));
   for ( i=0; i<level; i++ )
      Calloc(&(MacroMpole[i]), mp);
	 
   MacroM2M = (Mtype *)malloc(level*sizeof(Mtype));
   for ( i=0; i<level; i++ )
      Calloc(&(MacroM2M[i]), mp);
	 
   MacroM2L = (Mtype *)malloc(level*sizeof(Mtype));
   for ( i=0; i<level; i++ )
      Calloc(&(MacroM2L[i]), mp);
	 
#ifdef PMACRO
   Calloc(&(MacroTemp), mp);
#endif

} /* MacroInit */



/****************************************************************
 *
 *  MacroPreComp() - precompute the transfer matrices.  This should
 *  be called once initially and then once everytime the cell is 
 *  resized.
 *
 */

void MacroPreComp( Vector v1, Vector v2, Vector v3, Real edge_scale )
{   

   int i,j;                  /* loop variables */
   int x,y,z;                /* loop variables */
   int x2,y2,z2;             /* loop variables */

#ifndef PIPED
   Vector sedge;             /* scaled edge amount */
   Vector sep;               /* cell separation distance */

#else
   Vector sv1, sv2, sv3;     /* scaled ||-piped vectors */
   Vector dtmp;              /* dtmp = sv1+sv2+sv3 */
   Vector v1xv2, v2xv3, v3xv1; /* cross products */
   Real mag;                 /* vector magnitude */
   Real v2xv3dotv1,          /* triple product */
        v3xv1dotv2,
        v1xv2dotv3;
   Real xtmp1, ytmp1, ztmp1; /* temp distance vector */
   int iv1, iv2, iv3;
#endif

   Vector hv;                /* cell separation vector */
   IntVector dist;           /* cell region bounding vector */
   Real xtmp, ytmp, ztmp;    /* temp distance vector */
   Real rtmp;                /* temp distance */
   Real rad1, rad2;          /* temp distance */

#ifdef PMACRO
   int count;                /* counter to partition parallel exec */
#endif

#ifndef PMACRO
   /*
    * unless we are running parallel macro, only master pid needs to
    * precompute the x-fer matrices
    */
   if ( MyPid != MastPid ) {
      return;
   }
#endif

   /* if we specified no levels, return */
   if ( Kterm == 0 ) {
      return;
   }

   /*
    * we need to determin the size of the edges of the
    * cube scaled down to a unit cell.  we should really
    * be passed this value directly from the calling routine.
    */

#ifndef PIPED
   sedge.x = v1.x / edge_scale;
   sedge.y = v2.y / edge_scale;
   sedge.z = v3.z / edge_scale;
#else
   sv1.x = v1.x / edge_scale;
   sv1.y = v1.y / edge_scale;
   sv1.z = v1.z / edge_scale;
   sv2.x = v2.x / edge_scale;
   sv2.y = v2.y / edge_scale;
   sv2.z = v2.z / edge_scale;
   sv3.x = v3.x / edge_scale;
   sv3.y = v3.y / edge_scale;
   sv3.z = v3.z / edge_scale;
#endif

   /*
    * compute M2M transfer arrays:
    *
    *  for each level
    *   for each cell in the x,y,z direction
    *      add M2M xfer matrix into MpoleXfer[level];
    *   done x,y,z
    *  done level
    */

#ifndef PIPED
   sep.x = sedge.x;
   sep.y = sedge.y;
   sep.z = sedge.z;
#endif

#ifdef PMACRO
   count = 0;
#endif
   
   for ( i=1; i<Kterm; i++ ) {
      CMclear(MacroM2M[i],Mp);
#ifndef PIPED
      for ( x=-1; x<=1; x++ ) {
	 for ( y=-1; y<=1; y++ ) {
	    for ( z=-1; z<=1; z++ ) {

#ifdef PMACRO
	       if ( (count % Nprocs) == MyPid ) {
#endif
		  /* compute vector between center */
		  hv.x = (double)(x) * sep.x;
		  hv.y = (double)(y) * sep.y;
		  hv.z = (double)(z) * sep.z;
		  /* calculate transfer matrix */
		  addF(MacroM2M[i], Mp, hv);
#ifdef PMACRO
	       }
	       count++;
#endif
	    } /* for z */
	 } /* for y */
      } /* for x */

      /* increase separation distance for next level */
      sep.x *= 3.0;
      sep.y *= 3.0;
      sep.z *= 3.0;

#else /* PIPED */

      for ( iv1=-1; iv1<=1; iv1++ ) {
	 for ( iv2=-1; iv2<=1; iv2++ ) {
	    for ( iv3=-1; iv3<=1; iv3++ ) {

#ifdef PMACRO
	       if ( (count % Nprocs) == MyPid ) {
#endif
		  /* compute vector between center */

                  hv.x = (double)(iv1) * sv1.x + (double)(iv2) * sv2.x +
			 (double)(iv3) * sv3.x;
                  hv.y = (double)(iv1) * sv1.y + (double)(iv2) * sv2.y +
			 (double)(iv3) * sv3.y;
                  hv.z = (double)(iv1) * sv1.z + (double)(iv2) * sv2.z +
			 (double)(iv3) * sv3.z;
		  /* calculate transfer matrix */
		  addF(MacroM2M[i], Mp, hv);
#ifdef PMACRO
	       }
	       count++;
#endif
	    } /* for iv3 */
	 } /* for iv2 */
      } /* for iv1 */

      /* increase separation distance for next level */
      sv1.x *= 3.0;
      sv1.y *= 3.0;
      sv1.z *= 3.0;
      sv2.x *= 3.0;
      sv2.y *= 3.0;
      sv2.z *= 3.0;
      sv3.x *= 3.0;
      sv3.y *= 3.0;
      sv3.z *= 3.0;

#endif /* PIPED */



   } /* for i */


   /*
   * compute M2L transfer arrays:
   *
   * determine the maximum x/y/z directions for the cells
   * at the current level as goverened by theta.  this is
   * basically equivalent to the process used in generating
   * the relative interaction list.
   *
   * for each level bottom to top
   *   for each cell in the x,y,z direction
   *     if MAC(cell,theta) and not MAC(cell_parent,theta)
   *       rho = separation vector
   *       add xfer matrix into MpoleXfer[level];
   *     endif
   *   done x,y,z
   *  done level
   */

   /*
    * first, clear out the M2L arrays
    */

   for ( i=0; i<Kterm; i++ ) {
      CMclear(MacroM2L[i],Mp);
   }


   /*
   *  establish how far from the parent cell we have to go in order
   *  to exceed the MAC.  this will form the limits on our search space
   *  for the child cells.
   */

#ifndef PIPED
   rad1 = sqrt(sedge.x*sedge.x + sedge.y*sedge.y + sedge.z*sedge.z) * 0.5;
   rad2 = rad1 * 3.0;

   dist.x = 1;
   while ( MAC(rad2, rad2, (double)(dist.x*3)*sedge.x, Theta) == FALSE )
      dist.x++;

   dist.y = 1;
   while ( MAC(rad2, rad2, (double)(dist.y*3)*sedge.y, Theta) == FALSE )
      dist.y++;

   dist.z = 1;
   while ( MAC(rad2, rad2, (double)(dist.z*3)*sedge.z, Theta) == FALSE )
      dist.z++;

   for ( x=-dist.x; x<=dist.x; x++ ) {
      for ( y=-dist.y; y<=dist.y; y++ ) {
	 for ( z=-dist.z; z<=dist.z; z++ ) {

  	    /*
	     * if the cell meets the parents MAC, then we do
	     * do not have to consider it.
             */

            xtmp = (double)(x*3) * sedge.x;
            ytmp = (double)(y*3) * sedge.y;
	    ztmp = (double)(z*3) * sedge.z;
	    rtmp = sqrt( xtmp*xtmp + ytmp*ytmp + ztmp*ztmp );

	    if ( MAC(rad2,rad2,rtmp,Theta) == FALSE) {

               /*
		* okay, for each child in the 3x3 parent cell,
	        * lets check to see if it should be included
		*/
               for ( x2=-1; x2<=1; x2++ ) {
	          for ( y2=-1; y2<=1; y2++ ) {
	             for ( z2=-1; z2<=1; z2++ ) {

                        xtmp = (double)(3*x+x2) * sedge.x;
                        ytmp = (double)(3*y+y2) * sedge.y;
	                ztmp = (double)(3*z+z2) * sedge.z;
	                rtmp = sqrt( xtmp*xtmp + ytmp*ytmp + ztmp*ztmp );

	                if ( MAC(rad1,rad1,rtmp,Theta) == TRUE) {

			   j = 1;
			   for ( i=0; i<Kterm; i++ ) {
#ifdef PMACRO
			      if ( (count % Nprocs) == MyPid ) {
#endif
				 hv.x = xtmp*(double)(j);
				 hv.y = ytmp*(double)(j);
				 hv.z = ztmp*(double)(j);
				 addG(MacroM2L[i],Mp,hv);
#ifdef PMACRO
			      }
			      count++;
#endif
			      j *= 3;
			   } /* for i */

			} /* if MAC(2) */
		     } /* for z2 */
		  } /* for y2 */
	       } /* for x2 */
	    } /* if MAC */
	 } /* for z */
      } /* for y */
   } /* for x */
#else
   sv1.x = v1.x / edge_scale;
   sv1.y = v1.y / edge_scale;
   sv1.z = v1.z / edge_scale;
   sv2.x = v2.x / edge_scale;
   sv2.y = v2.y / edge_scale;
   sv2.z = v2.z / edge_scale;
   sv3.x = v3.x / edge_scale;
   sv3.y = v3.y / edge_scale;
   sv3.z = v3.z / edge_scale;

   dtmp.x = sv1.x + sv2.x + sv3.x;
   dtmp.y = sv1.y + sv2.y + sv3.y;
   dtmp.z = sv1.z + sv2.z + sv3.z;

   rad1 = sqrt(dtmp.x * dtmp.x + dtmp.y * dtmp.y + dtmp.z * dtmp.z);
   rad2 = rad1 * 3.0; 

   v2xv3.x = sv2.y*sv3.z - sv2.z*sv3.y;
   v2xv3.y = sv2.z*sv3.x - sv2.x*sv3.z;
   v2xv3.z = sv2.x*sv3.y - sv2.y*sv3.x;

   v3xv1.x = sv3.y*sv1.z - sv3.z*sv1.y;
   v3xv1.y = sv3.z*sv1.x - sv3.x*sv1.z;
   v3xv1.z = sv3.x*sv1.y - sv3.y*sv1.x;

   v1xv2.x = sv1.y*sv2.z - sv1.z*sv2.y;
   v1xv2.y = sv1.z*sv2.x - sv1.x*sv2.z;
   v1xv2.z = sv1.x*sv2.y - sv1.y*sv2.x;

   mag=Vec_Mag(&v2xv3);
   v2xv3dotv1 = (v2xv3.x * sv1.x + v2xv3.y * sv1.y + 
                v2xv3.z * sv1.z)/mag;
   mag=Vec_Mag(&v3xv1);
   v3xv1dotv2 = (v3xv1.x * sv2.x + v3xv1.y * sv2.y + 
                v3xv1.z * sv2.z)/mag;
   mag=Vec_Mag(&v1xv2);
   v1xv2dotv3 = (v1xv2.x * sv3.x + v1xv2.y * sv3.y + 
                v1xv2.z * sv3.z)/mag;
   dist.x = 1;
   dist.y = 1;
   dist.z = 1;

   while ( MAC(rad2, rad2, (double)(dist.x*3) * v2xv3dotv1, Theta) == FALSE)
       dist.x++;
   while ( MAC(rad2, rad2, (double)(dist.y*3) * v3xv1dotv2, Theta) == FALSE)
       dist.y++;
   while ( MAC(rad2, rad2, (double)(dist.z*3) * v1xv2dotv3, Theta) == FALSE)
       dist.z++;

   for ( x=-dist.x; x<=dist.x; x++ ) {
      for ( y=-dist.y; y<=dist.y; y++ ) {
	 for ( z=-dist.z; z<=dist.z; z++ ) {

  	    /*
	     * if the cell meets the parents MAC, then we do
	     * do not have to consider it.
             */

            xtmp1 = 3.0 * ((double)x * sv1.x + (double)y * sv2.x +
                   (double)z * sv3.x);
            ytmp1 = 3.0 * ((double)x * sv1.y + (double)y * sv2.y +
                   (double)z * sv3.y);
            ztmp1 = 3.0 * ((double)x * sv1.z + (double)y * sv2.z +
                   (double)z * sv3.z);
	    rtmp = sqrt( xtmp1*xtmp1 + ytmp1*ytmp1 + ztmp1*ztmp1 );

	    if ( MAC(rad2,rad2,rtmp,Theta) == FALSE) {

               /*
		* okay, for each child in the 3x3 parent cell,
	        * lets check to see if it should be included
		*/
   
               for ( x2=-1; x2<=1; x2++ ) {
	          for ( y2=-1; y2<=1; y2++ ) {
	             for ( z2=-1; z2<=1; z2++ ) {

                        xtmp = xtmp1 + ((double)(x2) * sv1.x +
                               (double)(y2) * sv2.x +
                               (double)(z2) * sv3.x);
                        ytmp = ytmp1 + ((double)(x2) * sv1.y +
                               (double)(y2) * sv2.y +
                               (double)(z2) * sv3.y);
                        ztmp = ztmp1 + ((double)(x2) * sv1.z +
                               (double)(y2) * sv2.z +
                               (double)(z2) * sv3.z);
	                rtmp = sqrt( xtmp*xtmp + ytmp*ytmp + ztmp*ztmp );

	                if ( MAC(rad1,rad1,rtmp,Theta) == TRUE) {

			   j = 1;
			   for ( i=0; i<Kterm; i++ ) {
#ifdef PMACRO
			      if ( (count % Nprocs) == MyPid ) {
#endif
				 hv.x = xtmp*(double)(j);
				 hv.y = ytmp*(double)(j);
				 hv.z = ztmp*(double)(j);
				 addG(MacroM2L[i],Mp,hv);
#ifdef PMACRO
			      }
			      count++;
#endif
			      j *= 3;
			   } /* for i */

			} /* if MAC(2) */
		     } /* for z2 */
		  } /* for y2 */
	       } /* for x2 */
	    } /* if MAC */
	 } /* for z */
      } /* for y */
   } /* for x */
#endif

#ifdef PMACRO
   /*
    * if we are not the zero pid, then we need to send the
    * results to be collected.
    */

   if ( MyPid != MastPid ) {
      MacroSend();
   } /* if MyPid */

   PreCompFlag = TRUE;


#endif

} /* MacroPreComp() */



/****************************************************************
 *
 *  MacroCompute() - perform the upward and downward pass of
 *    the macroscopic expansion code.
 *
 */

void MacroCompute( Mtype Min, Mtype Lout )
{

   /* if we specified no levels, return */
   if ( Kterm == 0 ) {
      return;
   }

#ifdef PMACRO
   /*
    * if we just recomputed the xfer matrices,
    * then receive them.
    */
   if ( PreCompFlag == TRUE ) {
      PreCompFlag = FALSE;
      if ( MyPid == MastPid ) {
	 MacroReceive();
      }
   }
#endif

   MacroUpward( Min );

   MacroDownward( Lout );

#ifdef DEBUG
   MacroDump( MacroMpole, Lout );
#endif
   
} /* MacroCompute */


/****************************************************************
 *
 *  MacroUpward() - given the multipole expansion of the unit
 *    cell, compute the expansions for the macroscopic assembly
 *
 */

void MacroUpward( Mtype Min )
{

   int i;     /* loop counter */

   /*
    * upward pass -
    *
    * for each level in the macro expansion
    *   M[level] = M2M(M[level-1],Xfer[level])
    *   level++
    * done
    */

   CMclear(MacroMpole[0],Mp);
   CMsum(Min,MacroMpole[0],Mp);

   for (i=1; i<Kterm; i++) {
      CMclear(MacroMpole[i],Mp);
      M2M_Cshort(MacroMpole[i-1],MacroMpole[i],MacroM2M[i],Mp);
   } /* for i */

} /* MacroUpward() */


/****************************************************************
 *
 *  MacroDownward() - given the multipole macro expansion of the unit
 *    cell, compute the local expansion for the macroscopic assembly
 *
 */

void MacroDownward( Mtype Lout )
{

   int i;     /* loop counter */

   /*
    * downward pass -
    * for each level in the macro expansion
    *   Local += M2L(M[level],Xfer[level])
    *   level--
    * done
    */

   CMclear(Lout,Mp);

   for (i=Kterm-1; i>=0; i--) {
      M2L_Cshort(MacroMpole[i],Lout,MacroM2L[i],Mp);
   } /* for i */

} /* MacroDownward() */

#ifdef PMACRO

/****************************************************************
 *
 *  MacroSend() - send macro xfer matrices to processor that
 *    owns the unit cell.
 *
 */

void MacroSend()
{

   int i;

   pvm_initsend(DATA_INPLACE_PVM);
   for ( i=0; i<Kterm; i++ ) {
      pvm_pkdouble(&(MacroM2M[i][0][0].x),MpeSize*2,1);
   }
   for ( i=0; i<Kterm; i++ ) {
      pvm_pkdouble(&(MacroM2L[i][0][0].x),MpeSize*2,1);
   }
   pvm_send(MastTid,MSG_MACRO);

} /* MacroSend() */


/****************************************************************
 *
 *  MacroReceive() - received and sum up transfer arrays from
 *    other processes.
 *
 */

void MacroReceive()
{

   int i,j;

   for ( i=1; i<Nprocs; i++ ) {
      pvm_recv(-1,MSG_MACRO);
      for ( j=0; j<Kterm; j++ ) {
	 pvm_upkdouble(&(MacroTemp[0][0].x),MpeSize*2,1);
	 CMsum(MacroTemp,MacroM2M[j],Mp);
      }

      for ( j=0; j<Kterm; j++ ) {
	 pvm_upkdouble(&(MacroTemp[0][0].x),MpeSize*2,1);
	 CMsum(MacroTemp,MacroM2L[j],Mp);
      }
   }

} /* MacroReceive() */

#endif

/****************************************************************
 *
 *  MacroCleanup() - free up any allocated storage
 *
 */

void MacroCleanup()
{
   int i;       /* loop variables */

   /* if we specify no Kterms, return */
   if ( Kterm == 0 ) {
      return;
   }

   /* allocate storage for everything */
   /* put the malloc checks in later */

   for ( i=0; i<Kterm; i++ )
      Cfree(MacroMpole[i], Mp);
   free(MacroMpole);
	 
   for ( i=0; i<Kterm; i++ )
      Cfree(MacroM2M[i], Mp);
   free(MacroM2M);
	 
   for ( i=0; i<Kterm; i++ )
      Cfree(MacroM2L[i], Mp);
   free(MacroM2L);
	 
#ifdef PMACRO
   Cfree(MacroTemp, Mp);
#endif

} /* MacroCleanup() */


/****************************************************************
 *
 *  MacroDump() - dump the contents of the macro expansion to
 *   a file for inspection.
 *
 */

void MacroDump( Mtype *Min, Mtype Lout)
{
   int  i,n,m;
   FILE *fp;

   fp = fopen("/tmp/DpmtaMacro.out","w");

   for (n = 0; n < Mp; n++) {
      for (m = 0; m <= n; m++) {
	 fprintf(fp, "%.6e %.6e   ", Lout[n][m].x, Lout[n][m].y);
      } /* for m */
      fprintf(fp, "\n");
   } /* for n */

   for ( i=0; i<Kterm; i++ ) {
      for (n = 0; n < Mp; n++) {
	 for (m = 0; m <= n; m++) {
	    fprintf(fp, "%.6e %.6e   ", Min[i][n][m].x, Min[i][n][m].y);
	 } /* for m */
	 fprintf(fp, "\n");
      } /* for n */
      fprintf(fp, "\n\n");
   } /* for i */

   fclose(fp);

} /* MacroDump() */
