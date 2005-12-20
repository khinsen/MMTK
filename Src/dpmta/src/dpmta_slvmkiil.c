/**********************************************************************
*
*  dpmta_slvmkiil.c - compute and send the inverse interaction
*     lists to the other slave processors.
*
*  w.t.rankin
*
*  Copyright (c) 1995 Duke University
*  All rights reserved
*
*/

static char rcsid[] = "$Id: dpmta_slvmkiil.c,v 2.9 1998/04/01 20:08:24 wrankin Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_slvmkiil.c,v $
 * Revision 2.9  1998/04/01 20:08:24  wrankin
 * added support for HILBERT ordering
 * general cleanup of code structure (more OOP if you will)
 *
 * Revision 2.8  1997/11/07 16:49:42  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 * Revision 2.7  1997/09/29 20:25:03  wrankin
 * fixed problem with invalid (empty) multipoles during upward pass.
 * cell indexing by processor was inconsistant between master/slave.
 *
 * Revision 2.6  1997/05/12 18:06:12  wrankin
 * added routines to clean up dynamic memory allocation when
 * PMTAinit() is called
 *
 * Revision 2.5  1997/03/26  20:36:28  wrankin
 * rewrote Ilist and IIlist code to handle reallocation of array sizes
 *   that happen as a resuly of re-scaling a non-cubic unit cell.
 *
 * Revision 2.4  1996/10/18  17:05:03  wrankin
 * improved PVM message passing structures.
 * made current and cleaned up T3D library code
 * added additional test codes
 * implement new heirarchical make structure
 *
 * Revision 2.3  1996/09/24  18:43:11  wrankin
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
 * Revision 2.2  1995/06/27  14:20:31  wrankin
 * Implementation of relative interaction list algorithm.
 * Interaction list generation code removed from master library.
 * Cell Table restructured as a 2-D table indexed by level/cell
 * All direct computations are now performed in a single pass.
 * General cleanup and removal of old files.
 *
 * Revision 2.1  1995/06/13  04:26:15  wrankin
 * No updates.
 * Sources placed under CVS.
 *
 * Revision 1.3  1995/04/09  22:21:06  wrankin
 * renaming of all global variables to prevent conflicts for
 * port of dpmta library to the T3D.
 *
 * Revision 1.2  1995/03/11  03:07:27  wrankin
 * fixed loop error in inv-ilist generation
 *
 * Revision 1.1  1995/03/06  07:37:03  wrankin
 * Initial revision
 *
 *
*/

/* include files */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pvm3.h"
#include "dpmta_pvm.h"
#include "dpmta_cell.h"
#include "dpmta_slvglobals.h"


/*
*  prototypes for external routines
*/

#include "dpmta_distmisc.h"
#include "dpmta_slvmkil.h"


/*
*  prototypes used internal to this module
*/

void init_ii_table();
void fill_ii_table( int );
void sort_ii_list();
void send_ii_list( int );
void recv_ii_list();
void dump_ii_list();
void dump_ii_count();


/*
*  global variables local to these procedures
*/

static IIsort     *SortTable;    /* int ilist sorting table */

static int        *MpSndBuf;     /* buffer for packing mp list */
static int        MpSndBufSz;    /* size of MpSndBuf */
static int        MpSndCnt;      /* number cells in MpSndBuf */
static int        *DirSndBuf;    /* buffer for packing mp list */
static int        DirSndBufSz;   /* size of MpSndBuf */
static int        DirSndCnt;     /* number cells in DirSndBuf */


/****************************************************************
 *
 * Init_Inv_Ilist() - inititalize data structures for the
 *  inverse interaction lists processing
 *
 */

void Init_Inv_Ilist()
{
   int i;
   int ncells, nlongs;

   /* allocate the inverse interaction list */
   Dpmta_IIlist.mlen = (int *)malloc(Dpmta_Nproc*sizeof(int));
   if ( Dpmta_IIlist.mlen == NULL) {
      fprintf(stderr,"ERROR: Init_Inv_Ilist - malloc() failed\n");
      exit(-1);
   }

   Dpmta_IIlist.dlen = (int *)malloc(Dpmta_Nproc*sizeof(int));
   if ( Dpmta_IIlist.dlen == NULL ) {
      fprintf(stderr,"ERROR: Init_Inv_Ilist - malloc() failed\n");
      exit(-1);
   }

   Dpmta_IIlist.msize = (int*)malloc(Dpmta_Nproc*sizeof(int));
   if ( Dpmta_IIlist.msize == NULL ) {
      fprintf(stderr,"ERROR: Init_Inv_Ilist - malloc() failed\n");
      exit(-1);
   }

   Dpmta_IIlist.dsize = (int*)malloc(Dpmta_Nproc*sizeof(int));
   if ( Dpmta_IIlist.dsize == NULL ) {
      fprintf(stderr,"ERROR: Init_Inv_Ilist - malloc() failed\n");
      exit(-1);
   }

   Dpmta_IIlist.mlist = (int**)malloc(Dpmta_Nproc*sizeof(int*));
   if ( Dpmta_IIlist.mlist == NULL ) {
      fprintf(stderr,"ERROR: Init_Inv_Ilist - malloc() failed\n");
      exit(-1);
   }

   Dpmta_IIlist.dlist = (int**)malloc(Dpmta_Nproc*sizeof(int*));
   if ( Dpmta_IIlist.dlist == NULL ) {
      fprintf(stderr,"ERROR: Init_Inv_Ilist - malloc() failed\n");
      exit(-1);
   }

   /* initialize mlen and dlen to prevent traversing problems later */
   for ( i=0; i<(Dpmta_Nproc); i++ ) {
      Dpmta_IIlist.dlen[i] = 0;
      Dpmta_IIlist.mlen[i] = 0;
      Dpmta_IIlist.dsize[i] = 0;
      Dpmta_IIlist.msize[i] = 0;
      Dpmta_IIlist.dlist[i] = (int *)NULL;
      Dpmta_IIlist.mlist[i] = (int *)NULL;
   }

   /* 
   *  Build a table of ulongs with enough bits for each cell
   *  and one row for each processor. zero out the bit fields
   */

   ncells = Dpmta_LevelLocate[Dpmta_NumLevels];
   nlongs = ncells/32 + 1;

   SortTable = (IIsort *)malloc(Dpmta_Nproc * sizeof(IIsort));
   if ( SortTable == NULL) {
      fprintf(stderr,"ERROR: Init_Inv_Ilist - malloc() failed\n");
      exit(-1);
   }

   for ( i=0; i<Dpmta_Nproc; i++) {

      SortTable[i].mexp = (int *)malloc(nlongs*sizeof(int));
      if ( SortTable[i].mexp == NULL) {
	 fprintf(stderr,"ERROR: Init_Inv_Ilist - malloc() failed\n");
	 exit(-1);
      }

      SortTable[i].direct = (int *)malloc(nlongs*sizeof(int));
      if ( SortTable[i].mexp == NULL) {
	 fprintf(stderr,"ERROR: Init_Inv_Ilist - malloc() failed\n");
	 exit(-1);
      }

   } /* for i */

   /*
   *  initialize other structures
   */

   MpSndBuf = NULL;
   MpSndBufSz = 0;
   DirSndBuf = NULL;
   DirSndBufSz = 0;


} /* Init_Inv_Ilist */


/****************************************************************
*
*  generate the inverse lists and send to each processor
*  then receive our own ilist from each processor.  In order to
*  keep the number of buffered messages down, we interleave the
*  sending and receiving of the interaction list so that we never
*  have more than one outstanding message at a time.
*  
*/

void Make_Inv_Ilist()
{
   int i;
   int procid;

   /* make sure there is more than one slave */
   if ( Dpmta_Nproc > 1 ) {

      /* allocate the perminate and temp buffers used */
      init_ii_table();

      /* traverse the cell table and create list */
      sort_ii_list();

      /*
      *  for each other processor, 
      *  unpack the appropriate bit fields, 
      *  pack and send a messages,
      *  then recv a message to keep buffering down
      */

      procid = (Dpmta_Pid + 1) % Dpmta_Nproc;
      fill_ii_table(procid);
      send_ii_list(procid);

      for ( i=2; i<Dpmta_Nproc; i++ ) {
         procid = (Dpmta_Pid + i) % Dpmta_Nproc;
         fill_ii_table(procid);
         recv_ii_list();
         send_ii_list(procid);
      } /* for i */

      recv_ii_list();

   } /* if Dpmta_Nproc */

#ifdef DEBUG
   dump_ii_list();
   dump_ii_count();
#endif

} /* Make_Inv_Ilist */   



/****************************************************************
*
*  Delete_Inv_List() - dealloc the contents of the inverse
*     interaction lists and other structures.
*
*/

void Delete_Inv_Ilist()
{
   int i;

   for ( i=0; i<Dpmta_Nproc; i++ ) {
      free(SortTable[i].mexp);
      free(SortTable[i].direct);
   } /* for i */
   free(SortTable);

   free(MpSndBuf);
   free(DirSndBuf);

   for ( i=0; i<Dpmta_Nproc; i++ ) {
      free(Dpmta_IIlist.mlist[i]);
      free(Dpmta_IIlist.dlist[i]);
   } /* for i */
   
   free(Dpmta_IIlist.mlist);
   free(Dpmta_IIlist.dlist);
   free(Dpmta_IIlist.mlen);
   free(Dpmta_IIlist.dlen);
   free(Dpmta_IIlist.msize);
   free(Dpmta_IIlist.dsize);

} /* Delete_Inv_Ilist() */



/****************************************************************
*
*  Inverse Interaction List Processing
*
*  the following routines generate the inverse interaction list
*  from the contents of the interaction list.
*
*  the algorith creates two zeroed out bitmaps for each processor.
*  the length of the bitmap is the number of cells.  the algorithm
*  then traverses the cell table and for takes the contents of each
*  interaction list and sets the bit in the approprite processors
*  bitmap, indicating that the cell is needed by this processor.
*
*/

void init_ii_table()
{
   int   i, j;
   int   ncells;
   int   nlongs;


   /* 
   *  Build a table of ulongs with enough bits for each cell
   *  and one row for each processor. zero out the bit fields
   */

   ncells = Dpmta_LevelLocate[Dpmta_NumLevels];
   nlongs = ncells/32 + 1;

   for (i=0; i<Dpmta_Nproc; i++) {
      SortTable[i].mcnt = 0;
      SortTable[i].dcnt = 0;
      for (j=0; j<nlongs; j++) {
         SortTable[i].mexp[j] = 0;
         SortTable[i].direct[j] = 0;
      } /* for j */
   } /* for i */

   /*
   *  initialize other structures
   */

   MpSndCnt = 0;
   DirSndCnt = 0;

} /* init_ii_table */


/****************************************************************
*
*  traverse the cell table and add the contents of the interaction
*  lists to the inverse interaction list
*/

void sort_ii_list()
{
   int i,j;
   int level, cell, proc;
   int tcell, pcell, index, mask;
   int temp1, temp2;
   int ovfl;

   /* cycle through all cells at all levels */
   for (level=Dpmta_DownPassStart; level<Dpmta_NumLevels; level++) {
      if ( Dpmta_Sindex[level] != -1 ) {
         for (j=Dpmta_Sindex[level]; j<=Dpmta_Eindex[level]; j++) {
	    cell = index2cell(j,level);

	    /* identify cell posn wrt parent */
            temp1 = cell & 0x07;

            /* cycle through all parent interaction cells */
            pcell = getparent(cell);
            for ( i=0; i<Dpmta_Intlist[temp1].pcnt; i++ ) {
               temp2 = Dpmta_Intlist[temp1].plist[i];

               if ( Cell2Cell(level-1,pcell,temp2,&tcell,&ovfl) ) {
                  proc = getslvpid(Dpmta_Nproc,level-1,tcell);
		  tcell += Dpmta_LevelLocate[level-1];
                  index = tcell / 32;
                  mask = 0x01 << (tcell % 32);
                  if ( (SortTable[proc].mexp[index] & mask) == 0 ) {
                     SortTable[proc].mcnt++;
                     SortTable[proc].mexp[index] |= mask;
		  } /* if sorttbl & mask */
	       } /* if Cell2Cell */

            } /* for i */

            /* cycle through all sibling interaction cells */
            for ( i=0; i<Dpmta_Intlist[temp1].scnt; i++ ) {
               temp2 = Dpmta_Intlist[temp1].slist[i];

               if ( Cell2Cell(level,cell,temp2,&tcell,&ovfl) ) {
                  proc = getslvpid(Dpmta_Nproc,level,tcell);
		  tcell += Dpmta_LevelLocate[level];
                  index = tcell / 32;
                  mask = 0x01 << (tcell % 32);
                  if ( (SortTable[proc].mexp[index] & mask) == 0 ) {
                     SortTable[proc].mcnt++;
                     SortTable[proc].mexp[index] |= mask;
		  } /* if sorttbl & mask */
	       } /* if Cell2Cell */

            } /* for i */

            /* cycle through all direct interaction cells */
            for ( i=0; i<Dpmta_Intlist[temp1].dcnt; i++ ) {
               temp2 = Dpmta_Intlist[temp1].dlist[i];

               if ( Cell2Cell(level,cell,temp2,&tcell,&ovfl) ) {
                  proc = getslvpid(Dpmta_Nproc,level,tcell);
		  tcell += Dpmta_LevelLocate[level];
                  index = tcell / 32;
                  mask = 0x01 << (tcell % 32);
                  if ( (SortTable[proc].direct[index] & mask) == 0 ) {
                     SortTable[proc].dcnt++;
                     SortTable[proc].direct[index] |= mask;
		  } /* if sorttbl & mask */
	       } /* if Cell2Cell */

            } /* for i */

	 } /* for cell */
      } /* if Dpmta_Sindex */
   } /* for level */

} /* sort_ii_list() */


/****************************************************************
*
*  this routine takes the data contained in the sorted table and
*  shifts out the cell ids and places them into an array for sending
*  to the other processors.
*/

void fill_ii_table( int dproc )
{

   int i,j;
   int cellid;
   int ncells, nlongs;
   int temp;

   /* realloc buffers if needed */
   if ( MpSndBufSz < SortTable[dproc].mcnt ) {
      MpSndBufSz = SortTable[dproc].mcnt;
      MpSndBuf = (int *)realloc( MpSndBuf, MpSndBufSz * sizeof(int));
      if ( MpSndBuf == NULL ) {
         fprintf(stderr,"Make_Inv_Ilist(): malloc() failed\n");
         exit(-1);
      }
   }

   if ( DirSndBufSz < SortTable[dproc].dcnt ) {
      DirSndBufSz = SortTable[dproc].dcnt;
      DirSndBuf = (int *)realloc( DirSndBuf, DirSndBufSz * sizeof(int));
      if ( DirSndBuf == NULL ) {
         fprintf(stderr,"Make_Inv_Ilist(): malloc() failed\n");
         exit(-1);
      }
   }

   MpSndCnt = 0;
   DirSndCnt = 0;

   ncells = Dpmta_LevelLocate[Dpmta_NumLevels];
   nlongs = ncells/32 + 1;

   for ( i=0; i<nlongs; i++ ) {
      temp = SortTable[dproc].mexp[i];
      for ( j=0; j<32; j++ ) {
         if ( temp & 0x01 ) {
            cellid = i * 32 + j;
            MpSndBuf[MpSndCnt] = cellid;
            MpSndCnt++;
	 } /* if temp */
         temp = temp >> 1;
      } /* for j */
   } /* for i */

   for ( i=0; i<nlongs; i++ ) {
      temp = SortTable[dproc].direct[i];
      for ( j=0; j<32; j++ ) {
         if ( temp & 0x01 ) {
            cellid = i * 32 + j;
            DirSndBuf[DirSndCnt] = cellid;
            DirSndCnt++;
	 } /* if temp */
         temp = temp >> 1;
      } /* for j */
   } /* for i */

} /* fill_ii_table */


/****************************************************************
*
*  send_ii_list() - sends inverse interaction lists to a processor
*
*  for a processr pointd to by index, send that processor the inverse
*  interaction list.
*  
*  message format:
*    1) int - sending processor id
*    2) int - number of elements in inverse mpe list
*    3) int[] - array of inv. inter elements, length spec. in (2)
*    4) int - number of elements in inverse particle list
*    5) int[] - array of inv. inter elements, length spec. in (4)
*
*/

void send_ii_list( int dproc )
{

   pvm_initsend(DATA_NORMAL_PVM);

   pvm_pkint(&(Dpmta_Pid),1,1);  
   pvm_pkint(&(MpSndCnt),1,1);
   if ( MpSndCnt != 0 )
      pvm_pkint(MpSndBuf,MpSndCnt,1);
   pvm_pkint(&(DirSndCnt),1,1);
   if ( DirSndCnt != 0 )
      pvm_pkint(DirSndBuf,DirSndCnt,1);

   pvm_send(Dpmta_Tids[dproc],MSG_INTR2);

} /* send_ii_list */


/****************************************************************
*
*  recv_ii_list() - receive inverse interaction list.
*
*  this procedure will receive an inverse interaction list from
*  another slave process and store the values in the inverse
*  interaction list structure.
*
*  message format:
*    1) int - sending processor id
*    2) int - number of elements in inverse mpe list
*    3) int[] - array of inv. inter elements, length spec. in (2)
*    4) int - number of elements in inverse particle list
*    5) int[] - array of inv. inter elements, length spec. in (4)
*
*/

void recv_ii_list()
{

   int sproc;               /* id of processor for inv interaction */
   int num_cells;           /* number of cells in inv interaction list */
   int *tmp_ptr;            /* temp integer pointer for allocation */

   /* wait for message */
   pvm_recv(-1,MSG_INTR2);

   pvm_upkint(&sproc,1,1);

   /*
   *  unpack multipole Dpmta_IIlist
   */

   pvm_upkint(&num_cells,1,1);
   Dpmta_IIlist.mlen[sproc] = num_cells;

   if ( num_cells > Dpmta_IIlist.msize[sproc] ) {
      tmp_ptr = Dpmta_IIlist.mlist[sproc];
      Dpmta_IIlist.mlist[sproc] = (int *)realloc(tmp_ptr,num_cells*sizeof(int));
      if ( Dpmta_IIlist.mlist[sproc] == NULL ) {
         fprintf(stderr,"Make_Inv_Ilist(): malloc() failed\n");
         exit(-1);
      }
      Dpmta_IIlist.msize[sproc] = num_cells;
   } /* if num_cells */

   if ( num_cells > 0 ) {
      pvm_upkint((Dpmta_IIlist.mlist[sproc]),num_cells,1);
   }

   /*
   *  unpack direct iilist
   */

   pvm_upkint(&num_cells,1,1);
   Dpmta_IIlist.dlen[sproc] = num_cells;

   if ( num_cells > Dpmta_IIlist.dsize[sproc] ) {
      tmp_ptr = Dpmta_IIlist.dlist[sproc];
      Dpmta_IIlist.dlist[sproc] = (int *)realloc(tmp_ptr,num_cells*sizeof(int));
      if ( Dpmta_IIlist.dlist[sproc] == NULL ) {
         fprintf(stderr,"Make_Inv_Ilist(): malloc() failed\n");
         exit(-1);
      }
      Dpmta_IIlist.dsize[sproc] = num_cells;
   } /* if num_cells */

   if ( num_cells > 0 ) {
      pvm_upkint((Dpmta_IIlist.dlist[sproc]),num_cells,1);
   }

} /* recv_ii_list */


/****************************************************************
*
*  dump_ii_list() - dump the contents of the current inverse
*     interaction list to a file
*
*/

void dump_ii_list()
{
   int i,j,k;
   int nlongs;
   int ncells;
   int temp;
   int cellid;
   FILE *fp;
   char filename[80];

   sprintf(filename,"/tmp/iilist.pid%d",Dpmta_Pid);
   fp = fopen(filename,"w");

   for ( i=0; i<Dpmta_Nproc; i++ ) {
      fprintf(fp,"Inv Ilist for Proc #%d\n",i);
      fprintf(fp,"----------------------\n");

      fprintf(fp,"processor %d:\n",i);
      fprintf(fp,"   mpe[%d]=",Dpmta_IIlist.mlen[i]);
      for (j=0; j<(Dpmta_IIlist.mlen[i]); j++) {
         fprintf(fp,"%d,",Dpmta_IIlist.mlist[i][j]);
      }
      fprintf(fp,"\n");

      fprintf(fp,"   dir[%d]=",Dpmta_IIlist.dlen[i]);
      for (j=0; j<(Dpmta_IIlist.dlen[i]); j++) {
         fprintf(fp,"%d,",Dpmta_IIlist.dlist[i][j]);
      }
      fprintf(fp,"\n");

      fprintf(fp,"====================================\n\n");

   } /* for i */

   fprintf(fp,"SortTable Listing:\n");

   ncells = Dpmta_LevelLocate[Dpmta_NumLevels];
   nlongs = ncells/32 + 1;

   for (i=0; i < Dpmta_Nproc; i++) {
      fprintf(fp,"Proc %d:\n",i);
      fprintf(fp," MPE[%d]: ",SortTable[i].mcnt);
      for (j = 0; j < nlongs; j++)
         fprintf(fp,"%08x ", SortTable[i].mexp[j]);
      fprintf(fp,"\n");
      fprintf(fp," Dir[%d]: ",SortTable[i].dcnt);
      for (j = 0; j < nlongs; j++)
         fprintf(fp,"%08x ", SortTable[i].direct[j]);
      fprintf(fp,"\n");
   } /* for i */
   fprintf(fp,"====================================\n\n");

   fprintf(fp,"SortTable Dump Listing:\n");
   for (i=0; i < Dpmta_Nproc; i++) {
      fprintf(fp,"Proc %d:\n",i);
      fprintf(fp," MPE[%d]: ",SortTable[i].mcnt);

      for ( j=0; j<nlongs; j++ ) {
	 temp = SortTable[i].mexp[j];
	 for ( k=0; k<32; k++ ) {
	    if ( temp & 0x01 ) {
	       cellid = j * 32 + k;
	       fprintf(fp,"%d ",cellid);
	    } /* if temp */
	    temp = temp >> 1;
	 } /* for k */
      } /* for j */
      fprintf(fp,"\n");

      fprintf(fp," Dir[%d]: ",SortTable[i].dcnt);
      for ( j=0; j<nlongs; j++ ) {
	 temp = SortTable[i].direct[j];
	 for ( k=0; k<32; k++ ) {
	    if ( temp & 0x01 ) {
	       cellid = j * 32 + k;
	       fprintf(fp,"%d ",cellid);
	    } /* if temp */
	    temp = temp >> 1;
	 } /* for k */
      } /* for j */
      fprintf(fp,"\n");
      
   } /* for i */
   fprintf(fp,"====================================\n\n");
   
   fclose(fp);

} /* dump_ii_list */


/****************************************************************
 *
 *  dump_ii_count() - dump counts of the ilist sizes
 *
 */

void dump_ii_count() {

   int i;
   int mtot, dtot;
   FILE *fp;
   char filename[80];

   sprintf(filename,"/home/scicomp0/wrankin/research/hilbert/count/tmp/iicount.pid%d",Dpmta_Pid);
   /*   sprintf(filename,"/tmp/iicount.pid%d",Dpmta_Pid); */
   fp = fopen(filename,"w");

   mtot = 0;
   dtot = 0;
   
   for (i=0; i < Dpmta_Nproc; i++) {
      if ( i != Dpmta_Pid ) {
	 mtot += SortTable[i].mcnt;
	 dtot += SortTable[i].dcnt;
      } /* if */
   } /* for i */


   fprintf(fp,"nproc pid lvl theta    pbc  m_total  d_total\n");
   fprintf(fp," %2d   %2d  %2d  %2f  %1d     %3d     %3d\n",
	   Dpmta_Nproc, Dpmta_Pid, Dpmta_NumLevels, Dpmta_Theta,
	   Dpmta_PBC, mtot, dtot);

   fprintf(fp,"  MPE:");
   for (i=0; i < Dpmta_Nproc; i++) {
      if ( i == Dpmta_Pid ) {
	 fprintf(fp," %3d",0);
      }
      else {
	 fprintf(fp," %3d",SortTable[i].mcnt);
      }
   }
   fprintf(fp,"\n");


   fprintf(fp,"  DIR:");
   for (i=0; i < Dpmta_Nproc; i++) {
      if ( i == Dpmta_Pid ) {
	 fprintf(fp," %3d",0);
      }
      else {
	 fprintf(fp," %3d",SortTable[i].dcnt);
      }
   } /* for i */
   fprintf(fp,"\n");

   fclose(fp);
      
} /* dump_ii_count() */
