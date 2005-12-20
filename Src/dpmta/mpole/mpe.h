/*
 *  mpe.h - type and constant definitions for multipole library
 *
 *  w. t. rankin
 *
 *  Copyright (C)1995 Duke University
 *  All Rights Reserved
 */

/*
 * RCS Id:
 *
 * $Id: mpe.h,v 1.7 1997/11/07 16:45:16 wrankin Exp $
 *
 * RSC History:
 *
 * $Log: mpe.h,v $
 * Revision 1.7  1997/11/07 16:45:16  wrankin
 * added missing prototype
 *
 * Revision 1.6  1997/09/26 04:27:59  wrankin
 * added routine to dump raw expansion data to file
 *
 * Revision 1.5  1997/05/13 17:52:12  wrankin
 * fixed irritation bug seen when math.h is not included
 * added a missing prototype to mpe.h
 *
 * Revision 1.4  1997/05/09  20:14:56  wrankin
 * added routines to de-allocate global arrays created in [C,LJ]init()
 * added routines to free up multipole expansion matrices
 * added LJ prototypes and ansi-fied more procedures
 *
 * Revision 1.3  1996/11/11  20:09:09  wrankin
 * added ANSI-C declarations to mpoleC and prototypes to mpe.h
 *
 * Revision 1.2  1995/09/15  15:08:53  wrankin
 * Added Lennard-Jones multipole routines to library.
 *
 * Revision 1.1.1.1  1995/07/10  13:11:46  wrankin
 * Initial release of the Multipole Library
 * Based upon W. Elliott's MDMA codes
 * Implements Colomb Potentials only (no LJ-potentials yet)
 *
 *
 */


#ifndef TRUE
   #define TRUE 1
#endif
#ifndef FALSE
   #define FALSE 0
#endif

#define SMALL_THETA 1.0e-10


typedef double Real;

typedef struct {
   Real x, y;
} Complex;

typedef struct {
   Real x, y, z;
} Vector;

typedef struct {
   Real r, a, b;
} SphVector;

typedef Complex **Mtype;
typedef Complex ***MtypeLJ;


/*
 *  prototypes for mpe_mpoleC.c
 */

int AddMultipoleC( Mtype, int, Real, Vector );
int Csize( int );
int CsizeF( int );
int ForceM_C( Mtype, int, Real, Vector, Real *, Vector * );
int Force_C( Mtype, int, Real, Vector, Real *, Vector * );
int Force_C_Y( Mtype, int, Real, Vector, Real *, Vector * );
int L2L_C( Mtype, Mtype, int, Vector );
int M2L_C( Mtype, Mtype, int, Vector );
int M2L_Cshort( Mtype, Mtype, Mtype, int );
int M2L_C_F( Mtype, Mtype, int, int, Vector );
int M2L_C_Fshort( Mtype, Mtype, Mtype, int, int );
int M2M_C( Mtype, Mtype, int, Vector );
int M2M_Cshort( Mtype, Mtype, Mtype, int );
int MCM_C( Mtype, Mtype, Mtype, int );

void CMclear( Mtype, int );
void CMclearF( Mtype, int );
void CMclearFrev( Mtype, int, int );
void CMclearFshort( Mtype, int, int );
void CMsum( Mtype, Mtype, int );
void CMsumF( Mtype, Mtype, int );
void Calloc( Mtype *, int );
void CallocF( Mtype *, int, int );
void CallocFrev( Mtype *, int, int );
void CallocFrevS( Mtype *, int, int );
void Ccleanup( int );
void CcleanupF(int, int );
void CcleanupFS( int, int );
void Cfree( Mtype, int );
void CfreeF( Mtype, int, int );
void CfreeFrev( Mtype, int, int );
void CfreeFrevS( Mtype, int, int );
void Cinit( int );
void CinitF( int, int );
void CinitFS( int, int );
void Fourier_C( int, Real );
void MathdumpY_C( Mtype, int, char * );
void MDumpRaw_C( Mtype, int, char * );
void Unwarp_M2L( Mtype, Mtype, int, int );
void Warp_M2L( Mtype, Mtype, int, int );
void Warp_Short( Mtype, int, int );
void addF( Mtype, int, Vector );
void addG( Mtype, int, Vector );
void copyF( Mtype, int, Vector );
void copyG( Mtype, int, Vector );
void dumpYF( Real *, int );
void dumpY_C( Mtype, int );
void makeF( int, SphVector );
void makeG( int, SphVector );
void makeYforceC( int, Real, Real, Real );

Real eval_lpotC( Mtype, int, Vector );
Real eval_mpotC( Mtype, int, Vector );


/*
 * prototypes for mpe_mpoleLJ.c and mpe_allocC.c
 */

void AddMultipoleLJ( MtypeLJ, int, Real, Vector );
void Force_LJ( MtypeLJ, int, Real, Vector, Real *, Vector * );
void Fourier_LJ( int, Real );
void Gegenbauer( int, Real );
void L2L_LJ( MtypeLJ, MtypeLJ, int, Vector );
void LJMclear( MtypeLJ, int );
void LJMsum( MtypeLJ, MtypeLJ, int );
void LJalloc( MtypeLJ *, int );
void LJcleanup( int );
void LJfree( MtypeLJ, int );
void LJinit( int );
void M2L_LJ( MtypeLJ, MtypeLJ, int, Vector );
void M2L_LJshort( MtypeLJ, MtypeLJ, MtypeLJ, int );
void M2M_LJ( MtypeLJ, MtypeLJ, int, Vector );
void copyYI( MtypeLJ, int, Vector );
void makeYI( int, SphVector );
void makeYII( int, SphVector );
void makeYIIforce( int, SphVector );
void makeYIIforce0( int, SphVector );

int LJsize( int );

Real eval_mpotLJ( MtypeLJ, int, Vector );

