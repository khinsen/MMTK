/*
 * fft.c
 *
 * Miscellaneous fft subroutines for Coulomb computations
 *
 * Copyright (c) 1995 Duke University
 * All Rights Reserved.
 *
 * Bill Elliott welliott@ee.duke.edu
 * W. T. Rankin <wrankin@ee.duke.edu>
 *
 */

static char RCSid[] = "$Id: mpe_fft.c,v 1.2 1997/11/03 18:46:42 wrankin Exp $";

/*
 *  RCS History:
 *
 * $Log: mpe_fft.c,v $
 * Revision 1.2  1997/11/03 18:46:42  wrankin
 * general cleanup/ansi-fication of code.  no new features.
 *
 * Revision 1.1.1.1  1995/07/10 13:11:47  wrankin
 * Initial release of the Multipole Library
 * Based upon W. Elliott's MDMA codes
 * Implements Colomb Potentials only (no LJ-potentials yet)
 *
 *
 */

#include <stdio.h>
#include <math.h>
#include "mpe.h"
#include "mpe_fftC.h"

row_fft( Real *row, int p )
{

   /* ffth(row, p);  */
   switch(p) {
      case 4:
         FFTH8(row);
         break;
      case 8:
         FFTH16(row);
         break;
      case 12:
      case 16:
         FFTH32(row);
         break;
      case 20:
      case 24:
      case 28:
      case 32:
         ffth(row, 32);
         break;
      default:
         fprintf(stderr, "wrong FFT size\n");
         exit(0);
   } /* switch */
} /* row_fft */
 

row_ifft( Real *row, int p )
{
   /*   iffth(row, p); */
   switch(p) {
      case 4:
         IFFTH8(row);
         break;
      case 8:
         IFFTH16(row);
         break;
      case 12:
      case 16:
         IFFTH32(row);
         break;
      case 20:
      case 24:
      case 28:
      case 32:
         iffth(row, 32);
         break;
      default:
         fprintf(stderr, "wrong FFT size\n");
         exit(0);
   } /* switch */
} /* row_ifft */


col_fft( Real *block, int p, int b )
{
   int i, j, blocklen, nblocks, pf;
   Real *tempptr;

   pf = 1 << ((int) (log((double)(2*p-1))/log(2.0)));

   blocklen = 4 * pf * b;
   nblocks = p/b;

   switch(b) {
      case 4:
         for (i=0; i < nblocks; i++)  {
            tempptr = &block[i*blocklen];
            FFTV8(tempptr, pf);
         }
         break;     
      default:
         for (i=0; i < nblocks; i++) 
            for (j=0; j < pf; j++)  {
               fftv(&block[i*blocklen+2*j], 2*b, pf, 1);
            }
         break;
   } /* switch */

} /* col_fft */


col_ifft( Real *block, int p, int b )
{
   int i, j, blocklen, nblocks, pf;
   Real *tempptr;

   pf = 1 << ((int) (log((double)(2*p-1))/log(2.0)));

   blocklen = 4 * pf * b;
   nblocks = p/b;

   switch(b) {
      case 4:
         for (i=0; i < nblocks; i++)  {
            tempptr = &block[i*blocklen];
            IFFTV8(tempptr, pf);
         }
         break;  
      default:
         for (i=0; i < nblocks; i++)
            for (j=0; j < pf; j++)  {
               fftv(&block[i*blocklen+2*j], 2*b, pf, -1);
            }
         break;
   } /* switch */

} /* col_ifft */

col_fftS( Real *block, int p, int b )
{
   int i, j, blocklen, nblocks, pf, pblock;
   Real *tempptr;

   if (b != 4) {
      fprintf(stderr, "Block length must be 4 for short FFT's\n");
      exit(0);
   } /* if */


/*   blocklen = 4 * pf * b; dead code? */
   nblocks = p/b;

   tempptr = &block[0];
   for (i=0; i < nblocks; i++)  {
      pblock = b * (i+1);
      pf = 1 << ((int) (log((double)(2*pblock-1))/log(2.0)));
      FFTV8(tempptr, pf);
      blocklen = 4 * pf * b;
      tempptr += blocklen;
   } /* for i */

} /* col_fftS */

col_ifftS( Real *block, int p, int b )
{
   int i, j, blocklen, nblocks, pf, pblock;
   Real *tempptr;

/*   if (b != 4) {
      fprintf(stderr, "Block length must be 4 for short FFT's\n");
      exit(0);
   }*/ /* if */


/*   blocklen = 4 * pf * b; dead code */
   nblocks = p/b;

   tempptr = &block[0];
   for (i=0; i < nblocks; i++)  {
      pblock = b * (i+1);
      pf = 1 << ((int) (log((double)(2*pblock-1))/log(2.0)));
      IFFTV8(tempptr, pf);
      blocklen = 4 * pf * b;
      tempptr += blocklen;
   } /* for i */

} /* col_ifftS */


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

fftv( Real *data, unsigned long nn, int sc, int isign )
{
   unsigned long n,mmax,m,j,istep,i;
   double wtemp,wr,wpr,wpi,wi,theta;
   double tempr,tempi;

   n=nn << 1;
   j=1;
   for (i=1;i<n;i+=2) {
      if (j > i) {
         SWAP(data[(j-1)*sc],data[(i-1)*sc]);
         SWAP(data[(j-1)*sc+1],data[(i-1)*sc+1]);
      }
      m=n >> 1;
      while (m >= 2 && j > m) {
         j -= m;
         m >>= 1;
      }
      j += m;
   }
   mmax=2;
   while (n > mmax) {
      istep=mmax << 1;
      theta=isign*(6.28318530717959/mmax);
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (m=1;m<mmax;m+=2) {
         for (i=m;i<=n;i+=istep) {
            j=i+mmax;
            tempr=wr*data[(j-1)*sc]-wi*data[(j-1)*sc+1];
            tempi=wr*data[(j-1)*sc+1]+wi*data[(j-1)*sc];
            data[(j-1)*sc]=data[(i-1)*sc]-tempr;
            data[(j-1)*sc+1]=data[(i-1)*sc+1]-tempi;
            data[(i-1)*sc] += tempr;
            data[(i-1)*sc+1] += tempi;
         }
         wr=(wtemp=wr)*wpr-wi*wpi+wr;
         wi=wi*wpr+wtemp*wpi+wi;
      }
      mmax=istep;
   }
}


int ffth( Real *row, int len )
{
   Real newrow[2*MAXP], negone;
   int m, i;

   for (m=0; m < 4*len; m++)
      newrow[m] = 0.0;

   newrow[0] = row[0];
   newrow[1] = 0.0;
   for (m=1; m < len; m++) {
      newrow[2*m] = row[2*m];
      newrow[2*m+1] = row[2*m+1];
      negone = 1.0 - 2.0 * (Real) (0x0001 & m);
      newrow[4*len-2*m] = row[2*m] * negone;
      newrow[4*len-2*m+1] = -row[2*m+1] * negone;
   } /* for m */

   four1(newrow-1, 2*len, 1);

   for (m=0; m < len; m++) {
      row[2*m] = newrow[2*m];
      row[2*m+1] = newrow[2*m+1];
   } /* for m */

} /* ffth */
     
int iffth( Real *row, int len )
{
   Real newrow[2*MAXP], negone;
   int m, i;

   for (m=0; m < 4*len; m++)
      newrow[m] = 0.0;

   for (m=0; m < len; m++) {
      newrow[2*m] = row[2*m];
      newrow[2*m+1] = row[2*m+1];
      newrow[2*(m+len)] = row[2*m];
      newrow[2*(m+len)+1] = -row[2*m+1];
   } /* for m */

   four1(newrow-1, 2*len, -1);

   for (m=0; m < len; m++) {
      row[2*m] = newrow[2*m];
      row[2*m+1] = newrow[2*m+1];
   } /* for m */

} /* iffth */
 


int four1( Real *data, unsigned long nn, int isign )
{
   unsigned long n,mmax,m,j,istep,i;
   double wtemp,wr,wpr,wpi,wi,theta;
   double tempr,tempi;

   n=nn << 1;
   j=1;
   for (i=1;i<n;i+=2) {
      if (j > i) {
         SWAP(data[j],data[i]);
         SWAP(data[j+1],data[i+1]);
      }
      m=n >> 1;
      while (m >= 2 && j > m) {
         j -= m;
         m >>= 1;
      }
      j += m;
   }
   mmax=2;
   while (n > mmax) {
      istep=mmax << 1;
      theta=isign*(6.28318530717959/mmax);
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (m=1;m<mmax;m+=2) {
         for (i=m;i<=n;i+=istep) {
            j=i+mmax;
            tempr=wr*data[j]-wi*data[j+1];
            tempi=wr*data[j+1]+wi*data[j];
            data[j]=data[i]-tempr;
            data[j+1]=data[i+1]-tempi;
            data[i] += tempr;
            data[i+1] += tempi;
         }
         wr=(wtemp=wr)*wpr-wi*wpi+wr;
         wi=wi*wpr+wtemp*wpi+wi;
      }
      mmax=istep;
   }
}
#undef SWAP
