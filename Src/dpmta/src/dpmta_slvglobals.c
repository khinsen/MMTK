/* 
*  dpmta_slvglobals.c - global variables used by DPMTA
*
*  w. t. rankin
*
*  Copyright (c) 1997 Duke University
*  All rights reserved
*
*  this file contains the declarations for all the global
*  data structures accessed by the DPMTA routines.
*
*  these were originally contained in the dpmta_slave.c file,
*  but were moved here for clarity.
*
*/

static char rcsid[] = "$Id";

/*
 * revision history:
 *
 * $Log: dpmta_slvglobals.c,v $
 * Revision 2.2  1998/04/01 20:08:17  wrankin
 * added support for HILBERT ordering
 * general cleanup of code structure (more OOP if you will)
 *
 * Revision 2.1  1997/11/07 16:49:26  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 *
*/

/* include files */

#include "dpmta_cell.h"          /* data type definitions */
#include "dpmta_pvm.h"           /* processor max definiton */


/*
*  global variables - exported in dpmta_slvglobals.h
*/

CellPtrPtrPtr Dpmta_CellTbl;          /* pointer to cell table list */
IlistPtr      Dpmta_Intlist;          /* interaction list table */
HlistPtr      Dpmta_Hlist;            /* mpe transfer matrix table */
IIdata        Dpmta_IIlist;           /* inverse interaction data table */

int Dpmta_Pid;                        /* process id of slave */
int Dpmta_Nproc;                      /* total number of processors */
int Dpmta_MyTid;                      /* my task id */
int Dpmta_Tids[MAXPROC];              /* task ids of other slave processes */
int Dpmta_MasterTid;                  /* task id of master process */

int Dpmta_CallingNum;                 /* interface for slave particles */
int Dpmta_CallingTids[MAXPROC];       /* task ids of slave procs */

int Dpmta_NumLevels;                  /* number of levels in spatial decomp */
int Dpmta_DownPassStart;              /* level where downward pass starts */
int Dpmta_FFT;                        /* FFT processing flag */
int Dpmta_PBC;                        /* periodic boundary cond. flag */
int Dpmta_Mp;                         /* # terms in the multipole exp (p) */
int Dpmta_FftBlock;                   /* FFT Blocking size (b) */
int Dpmta_MpeSize;                    /* # Complex in the multipole exp */
int Dpmta_LclSize;                    /* # Complex in the local exp */

double Dpmta_Theta;                   /* multipole acceptance parameter */

Vector Dpmta_CellVector1;             /* vectors describing the three */
Vector Dpmta_CellVector2;             /*   unit cell edges. */
Vector Dpmta_CellVector3;

double Dpmta_MaxCellLen;	      /* lenth of longest side */
Vector Dpmta_CellCenter;              /* position of simulation cube center */

int Dpmta_Resize;                     /* flag to initiate cell resizing */

#ifdef PIPED
double Dpmta_CV1Mag;
double Dpmta_CV2Mag;
double Dpmta_CV3Mag;
#endif

int Dpmta_Sindex[LEVELS_MAX];         /* index to first owned cell */
int Dpmta_Eindex[LEVELS_MAX];         /* index to last owned cells */
int Dpmta_RMcell[LEVELS_MAX];         /* index to number of mpe's xfer'd */
int Dpmta_RLcell[LEVELS_MAX];         /* index to number of loc's xfer'd */

Mtype Dpmta_Temp_Mpe;                 /* temporary multipole exp buffer */

#ifdef COMP_LJ
int Dpmta_Mp_LJ;                      /* # terms in mpe for LJ potential */
int Dpmta_MpeSize_LJ;                 /* # Complex in the LJ multipole exp */
int Dpmta_LclSize_LJ;                 /* # Complex in the LJ local exp */

MtypeLJ Dpmta_Temp_Mpe_LJ;            /* temporary LJ multipole exp buffer */
#endif

#if defined VIRIAL || defined OLDVIRIAL
double Dpmta_Vpot;                    /* virital potential */
Vector Dpmta_Vf;                      /* virital force summation */
#ifdef COMP_LJ
double Dpmta_Vpot_LJ;                 /* virital potential */
Vector Dpmta_Vf_LJ;                   /* virital force summation */
#endif
#endif

#ifdef MACROSCOPIC
int Dpmta_K;                          /* # of levels in macro expansion */
#endif


/****************************************************************
*
*  some other global definitions
*
*  note, these two arrays only have eleven entries simply because 
*  a thirty two bit integer cannot represent larger numbers.
*
*/

int Dpmta_Power8[] = { 1, 8, 64, 512, 4096, 32768, 262144, 2097152,
                        16777216, 134217728, 1073741824 };

int Dpmta_LevelLocate[] = { 0, 1, 9, 73, 585, 4681, 37449, 299593,
                        2396745, 19173961, 153391689 };



