/*
*  dpmta.h - include file for DPMTA data structures and prototype
*    declarations
*
*  w. t. rankin
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*/

/* $Id: dpmta.h,v 2.9 1997/05/07 21:27:44 chumphre Exp $
 *
 * revision history:
 *
 * $Log: dpmta.h,v $
 * Revision 2.9  1997/05/07 21:27:44  chumphre
 * cubectr changed to cellctr and other minor changes related to the interface
 *
 * Revision 2.8  1997/05/07  19:31:41  chumphre
 * implement a uniform interface for rectangular and ||-piped cells
 *
 * Revision 2.7  1997/05/07  18:59:13  chumphre
 * major changes to implement ||-piped cells
 *
 * Revision 2.6  1996/09/24  18:41:30  wrankin
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
 * Revision 2.5  1996/08/09  15:30:33  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.4  1996/02/29  21:13:16  wrankin
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
 * Revision 2.3  1995/11/29  22:29:01  wrankin
 * addition of periodic boundary code
 * beginning work on virial code.
 *
 * Revision 2.2  1995/10/01  21:45:42  wrankin
 * implemented LJ multipole calculations - compile time option only
 *
 * Revision 2.1  1995/05/16  04:02:16  wrankin
 * added code to allow resizing the sim cube on the fly
 * the calling process now scales the input to a unit cube
 * and scales the results back to the proper value
 * the slave assumes all data sent to it is in a unit cube
 * with the bottom corner at <0,0,0>
 *
 * Revision 1.2  1995/02/24  21:17:07  wrankin
 * update PMTAinit prototype to include theta
 *
 * Revision 1.1  1994/11/30  17:48:54  wrankin
 * Initial revision
 *
 *
 *
*/

/* in case we're interfacing with C++ */
#ifdef __cplusplus
extern "C" {
#endif


typedef struct pmtavector {
   double x;
   double y;
   double z;
   }  PmtaVector;

/*
*  this is the minimum amount of data required to describe a particle
*  if the LJ processing has not been compiled into DPMTA, then the
*  alpha and beta parameters may be left unset.
*/

typedef struct pmtaparticle {
   PmtaVector p;               /* particle position */
   double q;                   /* particle charge */
   double a,b;                 /* alpha and beta paramenters for LJ */
   } PmtaParticle;

typedef PmtaParticle *PmtaParticlePtr;


/*
*  this structure holds the force and potential energy results
*/

typedef struct pmtapartinfo {
   PmtaVector f;               /* force vector */
   double v;                   /* potential */
   } PmtaPartInfo;

typedef PmtaPartInfo *PmtaPartInfoPtr;


/*
*  this structure holds the execution/setup parameters
*/

typedef struct pmtainitdata {
   int nprocs;
   int nlevels;
   int mp;
   int mp_lj;
   int fft;
   int fftblock;
   int pbc;
   int kterm;
   double theta;
   PmtaVector v1;
   PmtaVector v2;
   PmtaVector v3;
   PmtaVector cellctr;
   int calling_num;
   int *calling_tids;
} PmtaInitData;

typedef PmtaInitData *PmtaInitDataPtr;


/*
*  prototypes for interface routines
*/

int PMTAinit( PmtaInitDataPtr, int* );

int PMTAregister();

int PMTAforce( int, PmtaParticlePtr, PmtaPartInfoPtr, PmtaPartInfoPtr );

int PMTAresize( PmtaVector *, PmtaVector *, PmtaVector* , PmtaVector* );

int PMTAvirial( double*, PmtaVector*, double*, PmtaVector* );

int PMTAexit();

#ifdef __cplusplus
}
#endif
