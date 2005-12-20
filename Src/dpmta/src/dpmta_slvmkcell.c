/* 
*  dpmta_slvmkcell.c - routines to identify and allocate cell table 
*
*  w. t. rankin
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*  these routines are a repackaging of the dpmta_rcvpart.c module,
*  version 2.4, and dpmta_rcvilist.c version 2.3
*/

static char rcsid[] = "$Id: dpmta_slvmkcell.c,v 2.16 1998/04/29 18:36:23 wrankin Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_slvmkcell.c,v $
 * Revision 2.16  1998/04/29 18:36:23  wrankin
 * fixed code that creates RLcell counter - now works with row/col index
 *
 * Revision 2.15  1998/04/01 20:08:22  wrankin
 * added support for HILBERT ordering
 * general cleanup of code structure (more OOP if you will)
 *
 * Revision 2.14  1997/11/12 13:49:54  wrankin
 * updates to communications routines and general cleanup of code in prep
 *   for introducing load balancing functionality in hte near future
 *
 * Revision 2.13  1997/11/07 16:49:36  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 * Revision 2.12  1997/09/29 20:25:02  wrankin
 * fixed problem with invalid (empty) multipoles during upward pass.
 * cell indexing by processor was inconsistant between master/slave.
 *
 * Revision 2.11  1997/05/12 18:06:09  wrankin
 * added routines to clean up dynamic memory allocation when
 * PMTAinit() is called
 *
 * Revision 2.10  1997/05/07  18:59:36  chumphre
 * major changes to implement ||-piped cells
 *
 * Revision 2.9  1997/02/26  16:54:32  wrankin
 * added support for serial non-pvm compilation of DPMTA
 *
 * Revision 2.8  1996/11/01  02:26:46  wrankin
 * modifications to support cray t3d compilation
 * version update for 2.5 release
 *
 * Revision 2.7  1996/09/24  18:43:01  wrankin
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
 * Revision 2.6  1996/09/10  18:07:06  wrankin
 * fixed macroscopic processing for non-cubic cells.
 * fixed passing of local expansions for non-2^n processors
 *
 * Revision 2.5  1996/08/20  17:12:51  wrankin
 * Enhancement to allow any number of processes to be utilized.
 * Formally the number of precesses was limited to a power of 2.
 *
 * Revision 2.4  1996/08/09  15:31:00  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.3  1995/10/01  21:46:18  wrankin
 * implemented LJ multipole calculations - compile time option only
 *
 * Revision 2.2  1995/06/27  14:20:29  wrankin
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
 * Revision 1.7  1995/04/09  22:21:06  wrankin
 * renaming of all global variables to prevent conflicts for
 * port of dpmta library to the T3D.
 *
 * Revision 1.6  1995/01/13  16:46:45  wrankin
 * particle and force arrays are no longer allocated during
 *   cell table creation.
 * removed some extranious procedures
 *
 * Revision 1.5  1994/11/28  04:27:53  wrankin
 * added hooks to handle multiple slave calls to PMTAforce
 * added static allocation of particle data structures
 *
 * Revision 1.4  1994/11/24  13:54:24  wrankin
 * updates to use static particle structure allocation
 * in preparation for having multiple slave entry points
 *
 * Revision 1.3  1994/11/17  00:27:52  wrankin
 * added code to handled variable cell center
 *
 * Revision 1.2  1994/10/19  00:08:16  wrankin
 * allocation of cells is now done through bill elliotts routines
 * from dpmta_slvmultipole.c
 *
 * Revision 1.1  1994/10/14  05:10:51  wrankin
 * Initial revision
 *
 *
*/

/* include files */
#include <stdlib.h>
#include <stdio.h>
#include "dpmta_cell.h"
#include "dpmta_slvglobals.h"

/* a little prototyping */

#include "dpmta_distmisc.h"
#include "dpmta_slvmkil.h"

/*
 * prototyping for functions internal to module
 */

void cell_identify();
void cell_center( int, int );
void alloc_local_cell( CellPtrPtr );
void alloc_remote_cell( CellPtrPtr );
void free_cell( CellPtr );


/****************************************************************
*
*  this funtion will allocate and initialize the cell table for
*  a specific slave processor.  the individual cell and particle
*  structures will be allocated for the cells owned by that process.
*
*/

void Alloc_Cell_Table()
{

   int  num_cells;         /* total number of cells */
   int  i, j, k;           /* loop counters */
   int  id;                /* cell id */

   /*
   *  initialize the cell indices
   *  these indices identify which cells (for each level) that
   *  this processor owns.
   */

   cell_identify();

   /* allocate array of cell pointers */

   Dpmta_CellTbl = (CellPtrPtrPtr)malloc(Dpmta_NumLevels*sizeof(CellPtrPtr));
   if ( Dpmta_CellTbl == NULL ) {
      fprintf(stderr,"Alloc_Cell_Table(): malloc failed [1]\n");
      exit(-1);
   }
   num_cells = Dpmta_LevelLocate[Dpmta_NumLevels];
   Dpmta_CellTbl[0] = (CellPtrPtr)malloc(num_cells*sizeof(CellPtr));
   if ( Dpmta_CellTbl[0] == NULL ) {
      fprintf(stderr,"Alloc_Cell_Table(): malloc failed [2]\n");
      exit(-1);
   }
   for ( i=1; i<Dpmta_NumLevels; i++ )
      Dpmta_CellTbl[i] = &(Dpmta_CellTbl[0][Dpmta_LevelLocate[i]]);
   for ( i=0; i<num_cells; i++ )
      Dpmta_CellTbl[0][i] = (CellPtr) NULL;

   /*
   *  allocate the cell table entries owned by the process
   */

   /* cycle through each level */
   for ( i=0; i<Dpmta_NumLevels; i++ ) {
      if ( (Dpmta_Sindex[i]) != -1 ) {

	 /*
         *  all allocation is pretty much done.
         *  cycle through and setup pointers and initial values
         */

         for ( j=Dpmta_Sindex[i]; j<=Dpmta_Eindex[i]; j++ ) {
	    id = index2cell(j,i);

	    alloc_local_cell( &(Dpmta_CellTbl[i][id]) );

            Dpmta_CellTbl[i][id]->pid = Dpmta_Pid;
            Dpmta_CellTbl[i][id]->id = Dpmta_LevelLocate[i] + id;
            Dpmta_CellTbl[i][id]->n = 0;
            Dpmta_CellTbl[i][id]->mvalid = FALSE;
            Dpmta_CellTbl[i][id]->mdata->lvalid = FALSE;

         } /* for j */

         /*
         *  now, if we are not at the top level, we need to make sure
         *  that a parent has been allocated for each cell on the level
         *  that we just allocated.  we must have a local parent of our
         *  cells to hold the local mpe expansion.
         *
         *  the below implementation allocates only the parent cell
         *  one level above topmost cell the processor owns (if this
         *  cell exists - it doesn't for pid 0)
         */

         if ( i > 0 ) {

 	    /* cycle through all children and allocate parents */
 	    for ( j=Dpmta_Sindex[i]; j<=Dpmta_Eindex[i]; j++ ) {
	       id = index2cell(j,i);
               k = getparent(id);
	       alloc_local_cell( &(Dpmta_CellTbl[i-1][k]) );
	       Dpmta_CellTbl[i-1][k]->pid = getslvpid(Dpmta_Nproc, i-1, k);;
	       Dpmta_CellTbl[i-1][k]->id = Dpmta_LevelLocate[i-1] + k;
	       Dpmta_CellTbl[i-1][k]->n = 0;
	       Dpmta_CellTbl[i-1][k]->mvalid = FALSE;
	       Dpmta_CellTbl[i-1][k]->mdata->lvalid = FALSE;

	    } /* for j */
         } /* if i */
      } /* if Dpmta_Sindex */
   } /* for i */

   /*
   *  get size of the multipole and local expansions.  these are
   *  used later for sending and receiving MPE and local expansions.
   *  also allocate temporary buffers used to send cell MPEs
   *  between parent and child pids.
   *
   *  actually this should probably be moved to MultipoleSetup()
   */

   if (Dpmta_FFT) {
      Dpmta_MpeSize = CsizeF(Dpmta_Mp);
      CallocF(&(Dpmta_Temp_Mpe), Dpmta_Mp, Dpmta_FftBlock);
   }
   else {
      Dpmta_MpeSize = Csize(Dpmta_Mp);
      Calloc(&(Dpmta_Temp_Mpe), Dpmta_Mp);
   }

   Dpmta_LclSize = Csize(Dpmta_Mp);


#ifdef COMP_LJ
   Dpmta_MpeSize_LJ = LJsize(Dpmta_Mp_LJ);
   Dpmta_LclSize_LJ = LJsize(Dpmta_Mp_LJ);
   LJalloc(&(Dpmta_Temp_Mpe_LJ), Dpmta_Mp_LJ);
#endif

} /* Alloc_Cell_TAble */



/****************************************************************
*
*  the following procedure will compute an array of indices
*  into the cell table which identify, for a given level, 
*  which cells that processor owns.
*
*  the algorithm cycles through each level.  for a given level, 
*  if there are more processors than there are cells at that level,
*  then this processor may own at most a single cell.
*
*  otherwise, the process owns one or more cells at that level,
*  and we can compute the range of cells it owns.  note that
*  this makes the assumption that the cells owned by a process are
*  all adjacent within that level.
*
*  in addition, the procedure computes an array of the number of
*  multipole expansions that are received from other pids as a result
*  of the upward M2M pass of the algorithm.
*/

void cell_identify()
{

   int i,j,k,l,m,n;  	  /* loop counters, what else? */
   int id;		  /* temp cell id */

   for (i=0; i<Dpmta_NumLevels; i++) {

      Dpmta_Sindex[i] = getsindex(Dpmta_Nproc, Dpmta_Pid, i);
      Dpmta_Eindex[i] = geteindex(Dpmta_Nproc, Dpmta_Pid, i);

   } /* for i */

   /*
   *  compute multipole receive index -
   *
   *  cycle through all the levels (except the bottom).
   *  for each of the children of each cell, count the number of
   *  other pids which own those children.  this will be the number
   *  of processes that will send multipole expansions over the course
   *  of the upward pass.
   *
   *  note that this *assumes* that id another processor owns children
   *  of the current cell, then all these children are contiguous with
   *  respect to the index space.
   *
   *  note that like the start and end cell arrays computed above, this
   *  array will need to be recomputed anytime cells are reassigned to
   *  new processors.
   *
   *  if we are not doing PBC or Macroscopic Expansions, then
   *  there is no reason to propagate the MPEs up past one
   *  level above our max.
   */

   for (i=0; i<Dpmta_NumLevels-1; i++) {
      Dpmta_RMcell[i] = 0;
   }

   for (i=(Dpmta_DownPassStart-1); i<Dpmta_NumLevels-1; i++) {
      if ( Dpmta_Sindex[i] != -1 ) {
         for (j=Dpmta_Sindex[i]; j<=Dpmta_Eindex[i]; j++) {
	   /* id = index2cell(j,i); */
            k = getfirstchild(j);
            m = -1;
            for (l=k; l<(k+8); l++) {
               n = getslvpid_indx(Dpmta_Nproc,i+1,l);
               if ( m != n ) {
                  m = n;
                  if ( n != Dpmta_Pid ) {
	             Dpmta_RMcell[i]++;
                  } /* if n */
               } /* if m */
            } /* for l */
         } /* for j */
      } /* if Dpmta_Sindex */
   } /* for i */


   /*
   *  compute local exp receive index -
   *
   *  cycle through all the levels (except the top).  for each of
   *  the unique parents, count the number of parents which are
   *  owned by other pids.  this will be the number of processes that
   *  will send local expansions over the course of the downward
   *  pass.
   *
   *  note that like the start and end cell arrays computed above, this
   *  array will need to be recomputed anytime cells are reassigned to
   *  new processors.
   */

   for (i=0; i<Dpmta_NumLevels; i++) {
      Dpmta_RLcell[i] = 0;
   }

#ifdef OLDRLCOMP
   for (i=Dpmta_DownPassStart; i<Dpmta_NumLevels; i++) {
      if ( Dpmta_Sindex[i] != -1 ) {
         l = -1;
         for (j=Dpmta_Sindex[i]; j<=Dpmta_Eindex[i]; j++) {
	    id = index2cell(j,i);
            k = getparent(id);
            if ( l != k ) {
	       l = k;
               m = getslvpid(Dpmta_Nproc,i-1,k);
               if ( m != Dpmta_Pid ) {
                  Dpmta_RLcell[i]++;
               } /* if m */
            } /* if l */
         } /* for j */
      } /* if Dpmta_Sindex */
   } /* for i */
#else
   /* note downpassstart had better be greater than zero */
   for ( i=Dpmta_DownPassStart-1; i<Dpmta_NumLevels-1; i++ ) {
      for ( j=0; j<Dpmta_Power8[i]; j++ ) {
	 if ( (j<Dpmta_Sindex[i]) || (j>Dpmta_Eindex[i]) ) {
	    id = index2cell(j,i);
	    k = getfirstchild(id);
	    for ( l=k; l<k+8; l++ ) {
	       m = getslvpid(Dpmta_Nproc,i+1,l);
	       if ( m == Dpmta_Pid ) {
		  Dpmta_RLcell[i+1]++;
		  l = k+8;
	       } /* if m */
	    } /* for l */
	 } /* if j */
      } /* for j */
   } /* for i */
#endif

   /*
    * there is no need ot pass locals at the upper levels
    * if we are not doing PBC or Macroscopic Expansions, so
    * filter these out.
    */

#ifdef MACROSCOPIC
   if ( Dpmta_K == 0 ) {
      Dpmta_RLcell[1] = 0;
   }
#else
   Dpmta_RLcell[1] = 0;
#endif

   if ( Dpmta_PBC == 0 ) {
      Dpmta_RLcell[2] = 0;
   }

} /* cell_identify */


/****************************************************************
*
*  cell_center - set the center of the cell table
*
*/

void cell_center( int level, int cell )
{

   int i;                  /* loop counter */
   int mul;                /* number of cells per edge */
   int index;              /* cell id */
#ifndef PIPED
   int x,y,z;              /* cell coordinates */
   Vector cube;            /* size of the smallest cube */
#else
   int v1, v2, v3;
   double imul;
   Vector scaledv1, scaledv2, scaledv3;
   Vector tmpc;
#endif

#ifndef PIPED
   x = y = z = 0;

   mul = 1 << level;
   index = cell;
   for ( i=0; i<level; i++ ) {
      x |= (index & 01) << i;
      index >>= 1;
      y |= (index & 01) << i;
      index >>= 1;
      z |= (index & 01) << i;
      index >>= 1;
   }


   /* cube the length of a cube edge at that level. */

   cube.x = (1.0/((double) mul))* (Dpmta_CellVector1.x / Dpmta_MaxCellLen);
   cube.y = (1.0/((double) mul))* (Dpmta_CellVector2.y / Dpmta_MaxCellLen);
   cube.z = (1.0/((double) mul))* (Dpmta_CellVector3.z / Dpmta_MaxCellLen);

   Dpmta_CellTbl[level][cell]->p.x = cube.x*((double)x + 0.5);
   Dpmta_CellTbl[level][cell]->p.y = cube.y*((double)y + 0.5);
   Dpmta_CellTbl[level][cell]->p.z = cube.z*((double)z + 0.5);

#else
   v1 = v2 = v3 = 0;

   mul = 1 << level;
   index = cell;
   for ( i=0; i<level; i++ ) {
      v1 |= (index & 01) << i;
      index >>= 1;
      v2 |= (index & 01) << i;
      index >>= 1;
      v3 |= (index & 01) << i;
      index >>= 1;
   }

   scaledv1.x = Dpmta_CellVector1.x / Dpmta_MaxCellLen;
   scaledv1.y = Dpmta_CellVector1.y / Dpmta_MaxCellLen;
   scaledv1.z = Dpmta_CellVector1.z / Dpmta_MaxCellLen;
   scaledv2.x = Dpmta_CellVector2.x / Dpmta_MaxCellLen;
   scaledv2.y = Dpmta_CellVector2.y / Dpmta_MaxCellLen;
   scaledv2.z = Dpmta_CellVector2.z / Dpmta_MaxCellLen;
   scaledv3.x = Dpmta_CellVector3.x / Dpmta_MaxCellLen;
   scaledv3.y = Dpmta_CellVector3.y / Dpmta_MaxCellLen;
   scaledv3.z = Dpmta_CellVector3.z / Dpmta_MaxCellLen;

   tmpc.x = (scaledv1.x * ((double)v1 + 0.5)) +
            (scaledv2.x * ((double)v2 + 0.5)) +
            (scaledv3.x * ((double)v3 + 0.5));
   tmpc.y = (scaledv1.y * ((double)v1 + 0.5)) +
            (scaledv2.y * ((double)v2 + 0.5)) +
            (scaledv3.y * ((double)v3 + 0.5));
   tmpc.z = (scaledv1.z * ((double)v1 + 0.5)) +
            (scaledv2.z * ((double)v2 + 0.5)) +
            (scaledv3.z * ((double)v3 + 0.5));
   imul = 1.0 / (double)mul;

   Dpmta_CellTbl[level][cell]->p.x = tmpc.x * imul;
   Dpmta_CellTbl[level][cell]->p.y = tmpc.y * imul;
   Dpmta_CellTbl[level][cell]->p.z = tmpc.z * imul;
#endif
} /* cell_center */



/****************************************************************
*
*  Alloc_Ilist_Cells() - allocate the cell table entries for the
*     cells pointed to in the interaction lists.
*
*  this routine will traverse the set of 'private' cells in the
*  cell table and for each interaction list, check to make sure
*  that the complete set of cell table entries are allocated
*  for the remote cells that will be accessed by this processor.
*
*  note that we do not look at the double direct interaction list
*  since all the cells in this table are defined as being located
*  local to the current processor, and are already allocated.
*/

void Alloc_Ilist_Cells()
{
   int i,j,k;               /* loop counters */
   int id;	  	    /* cell id */
   int temp1, temp2;        /* temp cell ids */
   int sep;                 /* cell separation index */
   int pcell, plevel;       /* parent cell id */
   int ovfl;                /* overflow cell boundary flag */


   /*  cycle through all the levels */
   for ( i=Dpmta_DownPassStart; i<Dpmta_NumLevels; i++ ) {

      /*
      *  see if we own any cells at this level. if we do
      *  then go through their interaction lists and allocate any
      *  needed cells
      */

      if ( Dpmta_Sindex[i] != -1 ) {
         for ( j=Dpmta_Sindex[i]; j<=Dpmta_Eindex[i]; j++ ) {
	    id = index2cell(j,i);

	    /*
	    *  determine relative cell position
	    */
	    temp1 = id & 0x07;

            /*
            *  check direct interaction lists
            */
            for ( k=0; k<Dpmta_Intlist[temp1].dcnt; k++ ) {
	       sep = Dpmta_Intlist[temp1].dlist[k];
               if ( Cell2Cell(i,id,sep,&temp2,&ovfl) ) {

                  alloc_remote_cell( &(Dpmta_CellTbl[i][temp2]) );

		  Dpmta_CellTbl[i][temp2]->pid =
		     getslvpid(Dpmta_Nproc, i, temp2);
		  Dpmta_CellTbl[i][temp2]->id = Dpmta_LevelLocate[i] + temp2;

	       } /* if Cell2Cell() */
            } /* for k */

            /*
            *  check multipole sibling interaction lists
            */
            for ( k=0; k<Dpmta_Intlist[temp1].scnt; k++ ) {
	       sep = Dpmta_Intlist[temp1].slist[k];
               if ( Cell2Cell(i,id,sep,&temp2,&ovfl) ) {

                  alloc_remote_cell( &(Dpmta_CellTbl[i][temp2]) );

		  Dpmta_CellTbl[i][temp2]->pid =
		     getslvpid(Dpmta_Nproc, i, temp2);
		  Dpmta_CellTbl[i][temp2]->id = Dpmta_LevelLocate[i] + temp2;

	       } /* if Cell2Cell */
	    } /* for k */

            /*
            *  check multipole parental interaction lists
	    *
	    *  NOTE: ONLY IF (i>0) !!!!!!
            */
            for ( k=0; k<Dpmta_Intlist[temp1].pcnt; k++ ) {
	       sep = Dpmta_Intlist[temp1].plist[k];
	       pcell = getparent(id);
	       plevel = i-1;
               if ( Cell2Cell(plevel,pcell,sep,&temp2,&ovfl) ) {

                  alloc_remote_cell( &(Dpmta_CellTbl[plevel][temp2]) );

		  Dpmta_CellTbl[plevel][temp2]->pid =
		     getslvpid(Dpmta_Nproc, plevel, temp2);
		  Dpmta_CellTbl[plevel][temp2]->id =
		     Dpmta_LevelLocate[plevel] + temp2;
		  Dpmta_CellTbl[plevel][temp2]->n = 0;


	       } /* if Cell2Cell */
	    } /* for k */

	 } /* for j */
      } /* if Dpmta_Sindex[] */
   } /* for i */
 
} /* Alloc_Ilist_Cells */


/****************************************************************
*
*  Make_Cell_Centers() - cycle through all cells and set centers
*    of the cells.
*/

void Make_Cell_Centers()
{

   int i,j;

   for (i=0; i<Dpmta_NumLevels; i++) {
      for (j=0; j<Dpmta_Power8[i]; j++) {
         if ( Dpmta_CellTbl[i][j] != NULL ) {
 	    cell_center(i,j);
	 } /* if */
      } /* for j */
   }/* for i */

} /* Make_Cell_Centers */



/****************************************************************
*
*  Delete_Cell_Table() - free all dynamic data structures
*
*  this funtion will delete the cell table and all associated
*  structures.  it traverses the cell table, free()-ing each data
*  structure as it finds them.
*
*  it's slow and plodding, but at this stage of the endgame, who cares?
*
*/

void Delete_Cell_Table()
{

   int i;          /* loop counters */
   int num_cells;  /* size of cell table */

   num_cells = Dpmta_LevelLocate[Dpmta_NumLevels];

   for ( i=0; i<num_cells; i++ ) {
      free_cell( Dpmta_CellTbl[0][i] );
   } /* for i */
   
   /* free the cell table array */
   free( Dpmta_CellTbl[0] );
   free( Dpmta_CellTbl );

   /*
    *  here are some other data arrays that are allocated
    *  somewhere in this module and thus need to be free
    */

   if (Dpmta_FFT) {
      CfreeF(Dpmta_Temp_Mpe, Dpmta_Mp, Dpmta_FftBlock);
   }
   else {
      Cfree(Dpmta_Temp_Mpe, Dpmta_Mp);
   }
   
#ifdef COMP_LJ
   LJfree(Dpmta_Temp_Mpe_LJ, Dpmta_Mp_LJ);
#endif

} /* Delete_Cell_Table */


/****************************************************************
*
*  this funtion will travers the cell table and verify the
*  contents against the current cell allocation.  any missing
*  date structures will be added.
*
*  note that we assume that the basic cell table pointer array
*  has already been allocated via an initial call to the
*  Alloc_Cell_Table() function.
*
*/

void Realloc_Cell_Table()
{

   int  i, j, k;           /* loop counters */
   int  id;		   /* temp cell pointer */
   CellPtrPtr celltemp;    /* temporary cell pointer */

   /*
   *  allocate the cell table entries owned by the process
   */

   /* cycle through each level */
   for ( i=0; i<Dpmta_NumLevels; i++ ) {
      if ( (Dpmta_Sindex[i]) != -1 ) {

         /*
         *  all allocation is pretty much done.
         *  cycle through and setup pointers and initial values
         */

         for ( j=Dpmta_Sindex[i]; j<=Dpmta_Eindex[i]; j++ ) {
	    id = index2cell(j,i);
	    
	    /*
	     *  check to see if cell is allocated
	     */
	    celltemp = &(Dpmta_CellTbl[i][id]);
	    alloc_local_cell( celltemp );
	    
	    (*celltemp)->pid = Dpmta_Pid;
	    (*celltemp)->id = Dpmta_LevelLocate[i] + id;
	    (*celltemp)->n = 0;
	    (*celltemp)->mvalid = FALSE;
	    (*celltemp)->mdata->lvalid = FALSE;
	 } /* for j */
   
	 /*  now, if we are not at the top level, we need to
	  *  make sure that a parent has been allocated for each
	  *  cell on the level that we just allocated.  we must
	  *  have a local parent of our cells to hold the local
	  *  mpe expansion.
	  *
	  *  the below implementation allocates only the parent
	  *  cell one level above topmost cell the processor
	  *  owns (if this cell exists - it doesn't for pid 0)
	  *
	  */

         if ( i > 0 ) {

 	    /* cycle through all children and allocate parents */
 	    for ( j=Dpmta_Sindex[i]; j<=Dpmta_Eindex[i]; j++ ) {
	       id = index2cell(j,i);
               k = getparent(id);
	       celltemp = &(Dpmta_CellTbl[i-1][k]);
	       alloc_local_cell( celltemp );

	       (*celltemp)->pid = getslvpid(Dpmta_Nproc, i-1, k);
               (*celltemp)->id = Dpmta_LevelLocate[i-1] + k;
	       (*celltemp)->mvalid = FALSE;
	       (*celltemp)->n = 0;
	       (*celltemp)->mdata->lvalid = FALSE;
	    } /* for j */
         } /* if i */
      } /* if Dpmta_Sindex */
   } /* for i */

} /* Realloc_Cell_TAble */



/****************************************************************
*
*  alloc_local_cell();
*
*  this routine will check and allocate any missing datastructures for
*  a local cell
*
*/

void alloc_local_cell( CellPtrPtr lcell )
{

   MdataPtr mdatatemp;     /* temporary mdata pointer */
   
   /*
    *  check to see if there is any cell allocated
    */

   if ( (*lcell) == (CellPtr)NULL ) {
      (*lcell) = (CellPtr)malloc(sizeof(Cell));
      if ((*lcell) == NULL ) {
	 fprintf(stderr,"alloc_local_cell(): malloc() failed\n");
	 exit(-1);
      }
      (*lcell)->plist = (ParticlePtr)NULL;
      (*lcell)->psize = 0;
      (*lcell)->mdata = (MdataPtr)NULL;
      (*lcell)->m = (Mtype)NULL;
#ifdef COMP_LJ
      (*lcell)->m_lj = (MtypeLJ)NULL;
#endif
   } /* if (*lcell) */
   
   if ( (*lcell)->m == (Mtype)NULL ) {
      if (Dpmta_FFT)
	 CallocF(&((*lcell)->m), Dpmta_Mp, Dpmta_FftBlock);
      else
	 Calloc(&((*lcell)->m), Dpmta_Mp);
   } /* if (*lcell)->m */

#ifdef COMP_LJ
   if ( (*lcell)->m_lj == (MtypeLJ)NULL ) {
      LJalloc(&((*lcell)->m_lj), Dpmta_Mp_LJ);
   } /* if m_lj */
#endif

    
   mdatatemp = (*lcell)->mdata;
   if ( mdatatemp == (MdataPtr)NULL ) {
      mdatatemp = (MdataPtr)malloc(sizeof(Mdata));
      if ( mdatatemp == NULL ) {
	 fprintf(stderr,"alloc_local_cell(): malloc() failed\n");
	 exit(-1);
      }
      mdatatemp->flist = NULL;
      mdatatemp->part_id = NULL;
      mdatatemp->proc_id = NULL;
      
      Calloc(&(mdatatemp->l), Dpmta_Mp);
#ifdef COMP_LJ
      LJalloc(&(mdatatemp->l_lj), Dpmta_Mp_LJ);
      mdatatemp->f_lj = NULL;
#endif
      
      (*lcell)->mdata = mdatatemp;
   } /* if mdatatemp */

} /* alloc_local_cell */       


/****************************************************************
*
*  alloc_remote_cell();
*
*  this routine will check and allocate any missing
*  datastructures for a remote cell.  it is the same as the
*  local cell allocation routine, except for the alloclation of
*  the mdate strucure.
*
*/

void alloc_remote_cell( CellPtrPtr rcell )
{

   /*
    *  check to see if there is any cell allocated
    */

   if ( (*rcell) == (CellPtr)NULL ) {
      (*rcell) = (CellPtr)malloc(sizeof(Cell));
      if ((*rcell) == NULL ) {
	 fprintf(stderr,"alloc_remote_cell(): malloc() failed\n");
	 exit(-1);
      }
      (*rcell)->plist = (ParticlePtr)NULL;
      (*rcell)->psize = 0;
      (*rcell)->mdata = (MdataPtr)NULL;
      (*rcell)->m = (Mtype)NULL;
#ifdef COMP_LJ
      (*rcell)->m_lj = (MtypeLJ)NULL;
#endif
   } /* if (*rcell) */
   
   
   if ( (*rcell)->m == (Mtype)NULL ) {
      if (Dpmta_FFT)
	 CallocF(&((*rcell)->m), Dpmta_Mp, Dpmta_FftBlock);
      else
	 Calloc(&((*rcell)->m), Dpmta_Mp);
   } /* if (*rcell)-m */

#ifdef COMP_LJ
   if ( (*rcell)->m_lj == (MtypeLJ)NULL ) {
      LJalloc(&((*rcell)->m_lj), Dpmta_Mp_LJ);
   } /* if m_lj */
#endif

} /* alloc_remote_cell() */       


/****************************************************************
*
*  free_cell() - free up cell data structure
*
*  this routine will free up all the dynamically allocated
*  data structures representing one cell
*/

void free_cell( CellPtr cellid )
{

   if ( cellid != (CellPtr) NULL ) {

      /* check if we have an mdata structure for this cell */

      if ( cellid->mdata != (MdataPtr) NULL ) {

	 if ( cellid->mdata->flist != NULL ) {
	    free( cellid->mdata->flist );
	 }

#ifdef COMP_LJ
	 if ( cellid->mdata->f_lj != NULL ) {
	    free( cellid->mdata->f_lj );
	 }
#endif	    
	 if ( cellid->mdata->part_id != NULL ) {
	    free( cellid->mdata->part_id );
	 }
	    
	 if ( cellid->mdata->proc_id != NULL ) {
	    free( cellid->mdata->proc_id );
	 }
	    
	 if ( cellid->mdata->l != NULL ) {
	    Cfree( cellid->mdata->l, Dpmta_Mp );
	 }

#ifdef COMP_LJ
	 if ( cellid->mdata->l_lj != NULL ) {
	    LJfree( cellid->mdata->l_lj, Dpmta_Mp_LJ );
	 }
#endif

	 /* free up the mdata structure */
	 free( cellid->mdata );

      } /* if mdata */

      /* check if we have a particle list and free it */
      if ( cellid->plist != NULL ) {
	 free( cellid->plist );
      }

      /* free up multipoles */
      if ( cellid->m != NULL ) {
	 if (Dpmta_FFT) {
	    CfreeF( cellid->m, Dpmta_Mp, Dpmta_FftBlock );
	 }
	 else {
	    Cfree( cellid->m, Dpmta_Mp );
	 }
      }
      
#ifdef COMP_LJ
      if ( cellid->m_lj != NULL ) {
	 LJfree( cellid->m_lj, Dpmta_Mp_LJ );
      }
#endif	    

      /* finally, free the cell */
      free( cellid );

   } /* if cellid */

} /* free_cell() */
