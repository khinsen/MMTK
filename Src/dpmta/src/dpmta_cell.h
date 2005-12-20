/*
*  dpmta_cell.h - include file for DPMTA data structures
*
*  w. t. rankin
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*/

/* $Id: dpmta_cell.h,v 2.19 1997/05/13 17:53:21 wrankin Exp $
 *
 * revision history:
 *
 * $Log: dpmta_cell.h,v $
 * Revision 2.19  1997/05/13 17:53:21  wrankin
 * use basic types defined in mpe.h instead of redeclaring them
 *
 * Revision 2.18  1997/03/26  20:36:18  wrankin
 * rewrote Ilist and IIlist code to handle reallocation of array sizes
 *   that happen as a resuly of re-scaling a non-cubic unit cell.
 *
 * Revision 2.17  1996/11/18  19:29:22  wrankin
 * update to virial code to use multipole virial calculation
 * old virial code is turned on by defining OLDVIRIAL
 * direct test program now computes absolute virial tensor.
 *
 * Revision 2.16  1996/10/18  17:04:30  wrankin
 * improved PVM message passing structures.
 * made current and cleaned up T3D library code
 * added additional test codes
 * implement new heirarchical make structure
 *
 * Revision 2.15  1996/09/24  18:41:35  wrankin
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
 * Revision 2.14  1996/08/09  15:30:36  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.13  1996/02/29  21:13:18  wrankin
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
 * Revision 2.12  1995/09/18  19:02:23  wrankin
 * updates to move the scaling of the simulation box to the slave
 * process.
 *
 * Revision 2.11  1995/07/01  03:26:53  wrankin
 * initial cut at precomputation of mpe transfer functions.
 * this works for vanilla M2L.  it does not work for FFT enhancements.
 *
 * Revision 2.10  1995/06/27  14:20:15  wrankin
 * Implementation of relative interaction list algorithm.
 * Interaction list generation code removed from master library.
 * Cell Table restructured as a 2-D table indexed by level/cell
 * All direct computations are now performed in a single pass.
 * General cleanup and removal of old files.
 *
 * Revision 2.9  1995/04/09  22:21:06  wrankin
 * renaming of all global variables to prevent conflicts for
 * port of dpmta library to the T3D.
 *
 * Revision 2.8  1995/03/06  07:55:35  wrankin
 * added counters to IIsort data structure
 * removed unused structures
 *
 * Revision 2.7  1995/02/24  21:18:14  wrankin
 * added new datatypes for interaction list generation
 *
 * Revision 2.6  1995/01/13  17:06:36  wrankin
 * added 'psize' field to cell table
 *
 * Revision 2.5  1994/11/30  18:08:19  wrankin
 * added constant definitions for particle id encoding
 *
 * Revision 2.4  1994/11/24  13:54:24  wrankin
 * updates to use static particle structure allocation
 * in preparation for having multiple slave entry points
 *
 * Revision 2.3  1994/11/17  00:27:52  wrankin
 * added code to handled variable cell center
 *
 * Revision 2.2  1994/10/14  04:44:21  wrankin
 * added typedef for Real to support new multipole library
 *
 * Revision 2.1  1994/06/02  12:52:27  wrankin
 * No change.  Update version number for Release 2
 *
 * Revision 1.13  1994/06/02  03:26:50  wrankin
 * moved id field from particle structure and place in separate array
 *
 * Revision 1.12  1994/03/27  23:27:36  wrankin
 * added spherical vector data type
 *
 * Revision 1.11  1994/03/27  19:51:40  wrankin
 * added global definitions for support of Mps_Setup
 *
 * Revision 1.10  1994/03/11  19:38:56  wrankin
 * removed sinfo and minfo structures in order to make the data global
 *
 * Revision 1.9  1994/02/16  01:04:27  wrankin
 * interaction list processing now works
 *
 * Revision 1.8  1994/02/12  15:48:56  wrankin
 * added slave information structure
 *
 * Revision 1.7  1994/02/11  05:07:42  wrankin
 * added slaveinfo data type to pass processing parameters around
 *
 * Revision 1.6  1994/02/11  00:19:32  wrankin
 * moved power8[] and levellocate[] arrays to here, and made
 * them static so there are no compile conflicts.
 *
 * Revision 1.5  1994/02/09  19:17:16  wrankin
 * moved particle id to particle definition
 *
 * Revision 1.4  1994/02/09  04:19:44  wrankin
 * modifications to support bill elliots initialization code
 *
 * Revision 1.3  1994/02/08  22:44:47  wrankin
 * changes to support interaction list code
 *
 * Revision 1.2  1994/02/08  19:49:18  wrankin
 * syntax corrections for type definitions
 *
 * Revision 1.1  1994/02/03  17:04:54  wrankin
 * Initial revision
 *
 *
*/

/*
*  Constant definitions 
*/

/* do NOT raise this above 10.  things will break */
#define LEVELS_MAX  10

#define DIR_X  1
#define DIR_Y  2
#define DIR_Z  4

#define TRUE   1
#define FALSE  0

#define SQRT3    1.73205080757
#define SQRT3_2  0.86602540379

/*
*  constant definitions for particle id encoding
*
*  for the case where we have multiple slaves sending particle
*  data to the PMTA processing slaves, it is much easier if we
*  store both the particle ID and the source processor id in
*  the same field.  to this end, we will declare the bottom
*  20 bits of the particle id field to be the index of the
*  particle as originally sent to the processor.  the top
*  12 bit will hold the process id (0-4k) of the originating
*  process.
*
*  note this limits us to 4096 processors and 1 Meg particles
*  per process, but i think we can live within these constraints
*  for a few years, until 64 bit ints become more common ;-)
*
*  the following define shift distances and logic masks to encode
*  and extract the particle codes from the id fields.
*/

#define PID_SHIFT 20
#define PID_MASK  0x0FFF
#define PART_SHIFT 0
#define PART_MASK 0x0FFFFF

/*
*  constants used in decoding the interaction list separation
*  indices.  note the dependance on LEVELS_MAX.
*/

/* #define ILMASK1 0x09249249UL */
/* #define ILMASK2 0x12492492UL */
/* #define ILMASK3 0x24924924UL */

#define ILMASK1 (0x09249249UL>>(30-LEVELS_MAX*3))
#define ILMASK2 (0x12492492UL>>(30-LEVELS_MAX*3))
#define ILMASK3 (0x24924924UL>>(30-LEVELS_MAX*3))



/* data structures */

/*
*  cell table entry.
*
*  the cell table entry contains the minimum amount of data
*  for a cell.  this includes the cell position, the multipole
*  expansion for that cell, and a particle list, if needed.
*
*  if the cell is local to the processor an additional date structure
*  is created which contains the additional fields needed for the
*  multipole/force calculations.
*/

/* certain basic data types are defines in the multipole include file */

#include "mpe.h"

/* integer vectors for cell separation processing */

typedef struct intvector {
   int x;
   int y;
   int z;
   } IntVector;


/*
*  this is the minimum amount of data required to describe a particle
*/

typedef struct particle {
   Vector p;                   /* particle position */
   double q;                   /* particle charge */
#ifdef COMP_LJ
   double a,b;                 /* alpha and beta paramenters for LJ */
#endif
   } Particle;

typedef Particle *ParticlePtr;


/*
*  this is the additional data needed if the particle is local to a
*  processor and will have the multipole calculations done on it.
*/

typedef struct partinfo {
   Vector f;                   /* force vector */
   double v;                   /* potential */
   } PartInfo;

typedef PartInfo *PartInfoPtr;


/*
*  the following information is only needed if the cell is local to
*  the processor and will have the multipole calculations performed
*  on its particles
*/

typedef struct mdata {
   Mtype       l;              /* local expansion for Colomb interaction */
   int         lvalid;         /* local valid flag */
   int         *part_id;       /* array of particle idents */
   int         *proc_id;       /* array of processor idents */
   PartInfoPtr flist;          /* particle force/interaction list */
   Vector      h;              /* length of cell */
#ifdef COMP_LJ
   MtypeLJ     l_lj;           /* local expansion for LJ interaction */
   PartInfoPtr f_lj;           /* LJ force/interaction list */
#endif
   } Mdata;

typedef Mdata *MdataPtr;
   

/*
*  this is the structure of the cell table entry.  the mlist
*  pointer is set to NULL if the cell is not owned by the 
*  processor.
*/

typedef struct cell {
   int         id;             /* cell id */
   int         pid;            /* processor that owns cell */
   Vector      p;              /* cell position */
   Mtype       m;              /* multipole expansion */
   int         mvalid;         /* multipole valid flag */
   int         n;              /* number of particles in cell */
   int         psize;          /* length of alloc'd plist array */
   ParticlePtr plist;          /* list of particles in cell */
   MdataPtr    mdata;          /* list of multipole information */
#ifdef COMP_LJ
   MtypeLJ     m_lj;           /* LJ multipole expansion */
#endif
   } Cell;

typedef Cell *CellPtr;         /* we always declare pointers this way */
typedef CellPtr *CellPtrPtr;   /* it's an ugly declaration, but true */
typedef CellPtrPtr *CellPtrPtrPtr;   /* even worse */


/*
*
*  structure definitions for processing interaction lists
*
*/

typedef struct ilist {
   int *plist;                 /* list of parent cell interactions */
   int pcnt;                   /* number of parent cell interactions */
   int psize;                  /* size of plist array */
   int *slist;                 /* list of sibling cell interactions */
   int scnt;                   /* number of sibling cell interactions */
   int ssize;                  /* size of slist array */
   int *dlist;                 /* list of near (direct) cells */
   int dcnt;                   /* number of near (direct) cells */
   int dsize;                  /* size of dlist array */
} Ilist;

typedef Ilist *IlistPtr;

typedef struct iidata {
   int   *mlen;         /* lengths of mexp lists 1 for ea rcv processor */
   int   *dlen;         /* lengths of single direct lists 1 for ea rcv proc */
   int   *msize;        /* size of mlist[] arrays */
   int   *dsize;        /* size of dlist[] arrays */
   int   **mlist;       /* mexp lists for each rcv processor */
   int   **dlist;       /* single direct lists for each rcv processor */
   } IIdata;            /* data block with all info for a sending processor */

typedef IIdata *IIdataPtr;

typedef struct iisort {
   int mcnt;                   /* length of mexp bitfield array */
   int dcnt;                   /* length of dcnt bitfield array */
   int *mexp;                  /* bitfield array used in sort algorith */
   int *direct;                /* bitfield array used in sort algorith */
   } IIsort;


/*
*  structure definitions for processing mpe transfer matrices
*  and reletive vector lists.
*
*/

typedef struct hlist {
   Mtype  *plist;              /* list of parent x-fer matrices */
   Mtype  *slist;              /* list of sibling x-fer matrices */
#ifdef COMP_LJ
   MtypeLJ  *plist_lj;         /* list of parent LJ x-fer matrices */
   MtypeLJ  *slist_lj;         /* list of sibling LJ x-fer matrices */
#endif
#ifdef VIRIAL
   Vector *plist_vec;          /* corresponding direction vectors */
   Vector *slist_vec;          /* corresponding direction vectors */
#endif
   Vector *dlist_vec;          /* list of direct interaction vectors */

   int    psize;               /* size of the above arrays */
   int    ssize;               /* size of the above arrays */
   int    dsize;               /* size of the above arrays */
} Hlist;

typedef Hlist *HlistPtr;


