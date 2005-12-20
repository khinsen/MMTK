/****************************************************************
*
*  dpmta_slvglobals.h - global variables for dpmta pvm implementation
*
*  w. t. rankin
* 
*  these are the external references to the  global variables
*  used by the slave process during computation
*
*  Copyright (c) 1997 Duke University
*  All rights reserved
*
*  note: dpmta_cell.h need to be included prior to this file in
*  all source modules.
*
*/

/*
 *  RCS info: $Id: dpmta_slvglobals.h,v 2.2 1998/04/01 20:08:18 wrankin Exp $
 *
 *  revision history:
 *
 *  $Log: dpmta_slvglobals.h,v $
 *  Revision 2.2  1998/04/01 20:08:18  wrankin
 *  added support for HILBERT ordering
 *  general cleanup of code structure (more OOP if you will)
 *
 *  Revision 2.1  1997/11/07 16:49:28  wrankin
 *  massive cleanup of code.
 *   - ansi-fication and inclusion of prototypes
 *   - removed unused variables
 *   - all (except the test) code compiles with minimal warnings under gcc.
 *
 *
 *
*/


extern CellPtrPtrPtr Dpmta_CellTbl;    /* pointer to cell table list */
extern IlistPtr      Dpmta_Intlist;    /* interaction list table */
extern HlistPtr      Dpmta_Hlist;      /* mpe transfer matrix table */
extern IIdata        Dpmta_IIlist;     /* inverse interaction data table */

/* PVM processing parameters */

extern int Dpmta_Pid;                  /* process id of slave */
extern int Dpmta_Nproc;                /* total number of processors */
extern int Dpmta_MyTid;                /* my task id */
extern int Dpmta_Tids[];               /* task ids of other slave processes */
extern int Dpmta_MasterTid;            /* task id of master process */

extern int Dpmta_CallingNum;           /* interface for slave particles */
extern int Dpmta_CallingTids[];        /* task ids of slave procs */

/* MPE processing parameters */

extern int Dpmta_NumLevels;            /* number of levels in spatial decomp */
extern int Dpmta_DownPassStart;        /* level where downward pass starts */
extern int Dpmta_FFT;                  /* FFT processing flag */
extern int Dpmta_PBC;                  /* periodic boundary cond. flag */
extern int Dpmta_Mp;                   /* # terms in the multipole exp (p) */
extern int Dpmta_FftBlock;             /* FFT Blocking size (b) */
extern int Dpmta_MpeSize;              /* # Complex in the multipole exp */
extern int Dpmta_LclSize;              /* # Complex in the local exp */

extern double Dpmta_Theta;             /* multipole acceptance parameter */
extern Vector Dpmta_CellVector1;       /* ||-piped vectors and magnitudes */
extern Vector Dpmta_CellVector2;
extern Vector Dpmta_CellVector3;
extern double Dpmta_MaxCellLen;	       /* length of longest side */
extern Vector Dpmta_CellCenter;        /* position of simulation cube center */
extern int Dpmta_Resize;               /* flag to initiate cell resizing */

#ifdef PIPED
extern double Dpmta_CV1Mag;
extern double Dpmta_CV2Mag;
extern double Dpmta_CV3Mag;
#endif

#ifdef COMP_LJ
extern int Dpmta_Mp_LJ;                /* # terms in mpe for LJ potential */
extern int Dpmta_MpeSize_LJ;           /* # Complex in the LJ multipole exp */
extern int Dpmta_LclSize_LJ;           /* # Complex in the LJ local exp */

extern MtypeLJ Dpmta_Temp_Mpe_LJ;      /* temporary LJ multipole exp buffer */
#endif

/* misc processing variables */

extern int Dpmta_Sindex[];              /* index to first owned cell */
extern int Dpmta_Eindex[];              /* index to last owned cells */
extern int Dpmta_RMcell[];             /* index to number of mpe's xfer'd */
extern int Dpmta_RLcell[];             /* index to number of loc's xfer'd */

extern Mtype Dpmta_Temp_Mpe;           /* temporary multipole exp buffer */

#if defined VIRIAL || defined OLDVIRIAL
extern Real Dpmta_Vpot;                /* virital potential */
extern Vector Dpmta_Vf;                /* virital force summation */
#ifdef COMP_LJ
extern Real Dpmta_Vpot_LJ;             /* virital potential */
extern Vector Dpmta_Vf_LJ;             /* virital force summation */
#endif
#endif

extern int Dpmta_Power8[];             /* arrays to aid cell table indexing */
extern int Dpmta_LevelLocate[];        /* arrays to aid cell table indexing */

#ifdef MACROSCOPIC
extern int Dpmta_K;                    /* # of levels in macro expansion */
#endif
