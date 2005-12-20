/*
 * mpe_mpoleLJ.c - multipole expansion code for 1/r^6 LJ potential.
 *
 * bill elliott
 * w. t. rankin
 *
 * Copyright (c) 1995 Duke University
 * All Rights Reserved.
 *
 */

static char RcsId[] = "$Id: mpe_mpoleLJ.c,v 1.5 1997/11/03 18:46:51 wrankin Exp $";

/*
 * RCS History:
 *
 * $Log: mpe_mpoleLJ.c,v $
 * Revision 1.5  1997/11/03 18:46:51  wrankin
 * general cleanup/ansi-fication of code.  no new features.
 *
 * Revision 1.4  1997/05/09 20:15:05  wrankin
 * added routines to de-allocate global arrays created in [C,LJ]init()
 * added routines to free up multipole expansion matrices
 * added LJ prototypes and ansi-fied more procedures
 *
 * Revision 1.3  1995/11/22  16:44:21  wrankin
 * fix line in potential computation code.
 *
 * Revision 1.2  1995/10/01  21:35:02  wrankin
 * added LJsum() function to add multipoles
 * changed name of LJsize() funtion to match Csize()
 * general code cleanup
 *
 * Revision 1.1  1995/09/15  15:08:56  wrankin
 * Added Lennard-Jones multipole routines to library.
 *
 *
 */


/*
 * include files
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpe.h"


/*
 *  local arrays used in multipole calculations
 */

static MtypeLJ Y_LJ;
static MtypeLJ Y_LJf;
static MtypeLJ L;
static Complex *Yxy; 
static Real    ***YIcoeff;
static Real    ***A_LJ;
static Real    ***A_LJi;
static Real    **GegPoly;


/*****************************************************************
 *
 * AddMultipoleLJ() - calculates the 1/r^6 multipole expansion
 * of a particle of weight b_lj, at location x,y,z relative to
 * the origin. The result is accumulated into the blank
 * multipole expansion array M. The array size is p.
 *
 */

void AddMultipoleLJ(
   MtypeLJ    M,
   int        p,
   Real       b_lj,
   Vector     v )
{
   int        n, l, m;
   SphVector  sv;


   Cart2Sph(v, &sv);

   makeYII(p, sv);

   for (n=0; n < p; n++) {
      for (l=0; l <= n; l++)  {
         for (m = (0x0001 & n+l); m <= (n-l); m += 2) {
            M[n][l][m].x += b_lj * Y_LJ[n][l][m].x;
            M[n][l][m].y += b_lj * Y_LJ[n][l][m].y;
         } /* for m */
      } /* for l */
   } /* for n */


} /* MakeMultipoleLJ */


/*****************************************************************
 *
 * M2M_LJ() - translates a Lennard-Jones multipole expansion M1
 * to another expansion M2, both size p, along the vector x,y,z
 *
 */

void M2M_LJ(
   MtypeLJ M1,
   MtypeLJ M2,
   int     p,
   Vector  v )
{
   int       n, l, m, np, lp, mp;
   int       startm, endm;
   Real      atemp;
   SphVector sv;


   Cart2Sph(v, &sv);

   sv.r = -sv.r; /* don't ask me why... */
   makeYII(p, sv);

   atemp = 1.0;
   for (np=0; np < p; np++)  {
      for (lp=0; lp <= np; lp++)  {
         for (mp = 0x0001 & (np+lp); mp <= (np-lp); mp += 2)  {
            for (n=0; n <= np; n++) {
               for (l=0; l <= lp; l++)  {
                  startm = mp - (np-n) + (lp-l);
                  startm = (startm < -(n-l) ? -(n-l) : startm);
                  endm = mp + (np-n) - (lp-l);
                  endm = (endm > n-l ? n-l : endm);
                  if (startm <= endm) {
                     m = startm;
                     while (m < 0 && m <= endm) { 
                        M2[np][lp][mp].x += 
                           (M1[n][l][-m].x * Y_LJ[np-n][lp-l][mp-m].x +
                           M1[n][l][-m].y * Y_LJ[np-n][lp-l][mp-m].y);
                        M2[np][lp][mp].y += /* atemp *  */
                           (M1[n][l][-m].x * Y_LJ[np-n][lp-l][mp-m].y -
                           M1[n][l][-m].y * Y_LJ[np-n][lp-l][mp-m].x);
                        m += 2;
                     } /* while m */
                     while (m < mp && m <= endm) {
                        M2[np][lp][mp].x += 
                           (M1[n][l][m].x * Y_LJ[np-n][lp-l][mp-m].x - 
                           M1[n][l][m].y * Y_LJ[np-n][lp-l][mp-m].y);
                        M2[np][lp][mp].y += /* atemp *  */
                           (M1[n][l][m].x * Y_LJ[np-n][lp-l][mp-m].y + 
                           M1[n][l][m].y * Y_LJ[np-n][lp-l][mp-m].x);
                        m += 2;
                     } /* while m */
                     while (m <= endm) {
                        M2[np][lp][mp].x += 
                           (M1[n][l][m].x * Y_LJ[np-n][lp-l][m-mp].x + 
                           M1[n][l][m].y * Y_LJ[np-n][lp-l][m-mp].y);
                        M2[np][lp][mp].y += /* atemp * */
                           (M1[n][l][m].x * -Y_LJ[np-n][lp-l][m-mp].y +
                           M1[n][l][m].y * Y_LJ[np-n][lp-l][m-mp].x);
                        m += 2;
                     } /* while m */
                  } /* if startm <= endm (ie m exists) */
               } /* for l */
            } /* for n */
         } /* for mp */
      } /* for lp */
   } /* for np */
  
} /* M2M_LJ */


/*****************************************************************
 *
 * M2L_LJ() - converts a Lennard-Jones multipole expansion M in
 * to a local expansion L, shifting the original expansion along
 * the vector x, y, z.  Both expansions are of size p.
 *
 */

void M2L_LJ(
   MtypeLJ M,
   MtypeLJ L,
   int     p,
   Vector  v )
{
   int       np, lp, mp, n, l, m;
   int       startm, endm;
   Real      atemp;
   SphVector sv;


   /*
    * convert cartesean to spherical vector
    */

   Cart2Sph(v, &sv);

   makeYI(p, sv);

   atemp = 1.0;

   for (np=0; np < p; np++)  {
      for (lp=0; lp <= np; lp++)  {
         for (mp = 0x0001 & (np+lp); mp <= (np-lp); mp += 2)  {
            for (n=0; n < p-np; n++) {
               for (l=0; l <= n; l++)  {
                  for (m = 0x0001 & (n+l); m <= (n-l); m += 2)  {
                     L[np][lp][mp].x += /* atemp * */
                        (M[n][l][m].x * Y_LJ[np+n][lp+l][mp+m].x -
                        M[n][l][m].y * Y_LJ[np+n][lp+l][mp+m].y);
                     L[np][lp][mp].y += /* atemp * */
                        (M[n][l][m].x * Y_LJ[np+n][lp+l][mp+m].y +
                        M[n][l][m].y * Y_LJ[np+n][lp+l][mp+m].x);
                  } /* for m, m and mp positive */

                  m = (0x0001 & (n+l+1)) + 1; /* find first neg m */
                  if (m <= n-l) {
                     while (m <= mp && m <= n-l)  {
                        L[np][lp][mp].x += /* atemp * */
                           (M[n][l][m].x * Y_LJ[np+n][lp+l][mp-m].x +
                           M[n][l][m].y * Y_LJ[np+n][lp+l][mp-m].y);
                        L[np][lp][mp].y += /* atemp * */
                           (M[n][l][m].x * Y_LJ[np+n][lp+l][mp-m].y -
                           M[n][l][m].y * Y_LJ[np+n][lp+l][mp-m].x);
                        m += 2;
                     } /* while  m <= mp */
                  } /* if m */

                  while (m <= n-l)  {
                     L[np][lp][mp].x += /* atemp * */
                        (M[n][l][m].x * Y_LJ[np+n][lp+l][m-mp].x -
                        M[n][l][m].y * Y_LJ[np+n][lp+l][m-mp].y);
                     L[np][lp][mp].y += /* -atemp * */ -1.0 *
                        (M[n][l][m].x * Y_LJ[np+n][lp+l][m-mp].y +
                        M[n][l][m].y * Y_LJ[np+n][lp+l][m-mp].x);
                     m += 2;
                  } /* while n-l > m > mp */
               } /* for l */
            } /* for n */
         } /* for mp */
      } /* for lp */
   } /* for np */

}  /* M2L_LJ */


/*****************************************************************
 *
 * M2L_LJshort() - M2L using transfer matrix Y
 *
 */

void M2L_LJshort(
   MtypeLJ M,
   MtypeLJ L,
   MtypeLJ Y,
   int     p )

{
   int np, lp, mp, n, l, m;
   int startm, endm;
   Real atemp;
   Real *Yrowstart, *Mrowstart, *Mend, *Mptr, *Yptr, *Lxptr, *Lyptr;
   Real Mtempr, Mtempi, Ytempr, Ytempi;


   atemp = 1.0;
   for (np=0; np < p; np++)  {
      for (lp=0; lp <= np; lp++)  {
         for (mp = 0x0001 & (np+lp); mp <= (np-lp); mp += 2)  {
            Lxptr = &L[np][lp][mp].x; Lyptr = &L[np][lp][mp].y;
            for (n=0; n < p-np; n++) {
               for (l=0; l <= n; l++)  {

                  Yrowstart = &Y[np+n][lp+l][0].x;
		  Mrowstart = &M[n][l][0].x;
                  m = 0x0001 & (n+l);
                  Mptr = Mrowstart;
		  Mptr += (m << 1); 
                  Yptr = Yrowstart;
		  Yptr += ((mp+m) << 1);
                  for (; m <= (n-l); m += 2) { 
                     Mtempr = *Mptr++;
		     Mtempi = *Mptr;
                     Ytempr = *Yptr++;
		     Ytempi = *Yptr;
                     *Lxptr += Mtempr * Ytempr - Mtempi * Ytempi;
                     *Lyptr += Mtempr * Ytempi + Mtempi * Ytempr;
                     Mptr += 3;
		     Yptr += 3;
                  } /* for m, m and mp positive */ 

                  m = (0x0001 & (n+l+1)) + 1; /* find first neg m */
                  Mptr = Mrowstart;
		  Mptr += (m << 1);
                  Yptr = Yrowstart;
		  Yptr += ((mp-m) << 1);
                  if (m <= n-l) {
                     Mptr = Mrowstart;
		     Mptr += (m << 1);
                     Yptr = Yrowstart;
		     Yptr += ((mp-m) << 1);
                     while (m <= mp && m <= n-l) {
                        Mtempr = *Mptr++;
			Mtempi = *Mptr;
                        Ytempr = *Yptr++;
			Ytempi = *Yptr;
                        *Lxptr += Mtempr * Ytempr + Mtempi * Ytempi;
                        *Lyptr += Mtempr * Ytempi - Mtempi * Ytempr;
                        Mptr += 3;
			Yptr -= 5;
                        m += 2; 
                     } /* while  m <= mp */

		     Yptr = Yrowstart;
		     Yptr += ((m-mp) << 1); 
                     while (m <= n-l)  {
                        Mtempr = *Mptr++;
			Mtempi = *Mptr;
                        Ytempr = *Yptr++;
			Ytempi = *Yptr;
                        *Lxptr += Mtempr * Ytempr - Mtempi * Ytempi;
                        *Lyptr -= Mtempr * Ytempi + Mtempi * Ytempr;
                        Mptr += 3;
			Yptr += 3;
                        m += 2; 
                     } /* while n-l > m > mp */
                  } /* if m */
               } /* for l */
            } /* for n */
         } /* for mp */
      } /* for lp */
   } /* for np */

} /* M2L_LJshort */

/*****************************************************************
 *
 * L2L_LJ() - translates a local expansion L1 to another local
 * expansion L2, both of size p, along vector x, y, z
 *
 */

void L2L_LJ(
   MtypeLJ L1,
   MtypeLJ L2,
   int p,
   Vector v )
{
   int       np, lp, mp, n, l, m;
   int       startm, endm;
   Real      atemp;
   SphVector sv;


   Cart2Sph(v, &sv);

   makeYII(p, sv);

   atemp = 1.0;

   for (np=0; np <= p; np++)  {
      for (lp=0; lp <= np; lp++)  {
         for (mp = 0x0001 & (np+lp); mp <= (np-lp); mp += 2)  {
            for (n=np; n < p; n++) {
               for (l=lp; l <= n; l++)  {
                  startm = mp - (n-np) + (l-lp);
                  startm = (startm < -(n-l) ? -(n-l): startm);
                  endm = mp + (n-np) - (l-lp);
                  endm = (endm > n-l ? n-l : endm);
                  m = startm;
                  if (startm <= endm) {
                     while (m < 0 && m <= endm)  {
                        atemp = (1.0 - 2.0 * (Real) (0x0001 & (n-np)));
                        L2[np][lp][mp].x += atemp *
                           (L1[n][l][-m].x * Y_LJ[n-np][l-lp][-m+mp].x -
                           L1[n][l][-m].y * Y_LJ[n-np][l-lp][-m+mp].y);
                        L2[np][lp][mp].y += -atemp *
                           (L1[n][l][-m].x * Y_LJ[n-np][l-lp][-m+mp].y +
                           L1[n][l][-m].y * Y_LJ[n-np][l-lp][-m+mp].x);
                        m += 2;
                     } /* while m is negative */

                     while (m < mp && m <= endm)  {
                        atemp = (1.0 - 2.0 * (Real) (0x0001 & (n-np)));
                        L2[np][lp][mp].x += atemp *
                           (L1[n][l][m].x * Y_LJ[n-np][l-lp][-m+mp].x +
                           L1[n][l][m].y * Y_LJ[n-np][l-lp][-m+mp].y);
                        L2[np][lp][mp].y += atemp *
                           (L1[n][l][m].x * -Y_LJ[n-np][l-lp][-m+mp].y +
                           L1[n][l][m].y * Y_LJ[n-np][l-lp][-m+mp].x);
                        m += 2;
                     } /* while m-mp is negative */

                     while (m <= endm)  {
                        atemp = (1.0 - 2.0 * (Real) (0x0001 & (n-np)));
                        L2[np][lp][mp].x += atemp *
                           (L1[n][l][m].x * Y_LJ[n-np][l-lp][m-mp].x - 
                           L1[n][l][m].y * Y_LJ[n-np][l-lp][m-mp].y);
                        L2[np][lp][mp].y += atemp *
                           (L1[n][l][m].x * Y_LJ[n-np][l-lp][m-mp].y +
                           L1[n][l][m].y * Y_LJ[n-np][l-lp][m-mp].x);
                        m += 2;
                     } /* while mp < endm, m and m-mp positive */

                  } /* if startm <= endm (ie m exists) */
               } /* for l */
            } /* for n */
         } /* for mp */
      } /* for lp */
   } /* for np */

} /* L2L_LJ */


/*****************************************************************
 *
 * Force_LJ() - computes the Lennard-Jones force and potential
 * on the particle with charge q at position vector pv, due to
 * the local expansion L of size p.
 *
 * The results are returned as a potential, potp, and force
 * vector, fv.
 *
 */

void Force_LJ(
   MtypeLJ    Lin,
   int        p,
   Real       b_lj,
   Vector     pv,
   Real      *potp,
   Vector    *fv )
{

   int        n, l, m;
   Real       pot, fr, fa, fb;
   Real       realtemp, imgtemp, cosalpha, cosbeta, sinalpha, sinbeta;
   SphVector  sv;


   Cart2Sph(pv, &sv);

   pot = 0.0;
   fr = 0.0;
   fa = 0.0;
   fb = 0.0;

   for (n=0; n < p; n++) {
      for (l=0; l <= n; l++)  {
         for (m = (0x0001 & n+l); m <= (n-l); m += 2) {
            L[n][l][m].x = Lin[n][l][m].x  * A_LJ[n][l][m];
            L[n][l][m].y = Lin[n][l][m].y  * A_LJ[n][l][m];
         } /* for m */
      } /* for l */
   } /* for n */

   sv.r = -sv.r;
   makeYIIforce(p, sv);

   pot = Y_LJf[0][0][0].x * L[0][0][0].x;   
   for (n=1; n < p; n++) {
      for (l=0; l <= n; l++) {
         for (m = (0x0001 & (n+l)); m <= (n-l); m += 2) {

            if (m == 0) {
               pot += sv.r * Y_LJf[n][l][0].x * L[n][l][0].x;
               fr += (Real) n * Y_LJf[n][l][0].x * L[n][l][0].x;
               if (l > 0)
                  fa +=  L[n][l][0].x * (Real) (-l) * Y_LJf[n][l-1][0].x;
               if (l < n)
                  fa += L[n][l][0].x * (Real) (n-l) * Y_LJf[n][l+1][0].x;
            } /* if m == 0 */

            else {
               pot += 2.0 * sv.r * (Y_LJf[n][l][m].x * L[n][l][m].x -
                  Y_LJf[n][l][m].y * L[n][l][m].y);
               fr += 2.0 * (Real) n *  
                  (Y_LJf[n][l][m].x * L[n][l][m].x - 
                  Y_LJf[n][l][m].y * L[n][l][m].y);
               realtemp = 0.0;
	       imgtemp = 0.0;
               if (l > 0) {
                  realtemp = (Real) (-l) * Y_LJf[n][l-1][m].x;
                  imgtemp = (Real) (-l) * Y_LJf[n][l-1][m].y;
               }
               if (l < n) {
                  realtemp += (Real) (n-l) * Y_LJf[n][l+1][m].x;
                  imgtemp += (Real) (n-l) * Y_LJf[n][l+1][m].y; 
               }
               fa += 2.0 * (L[n][l][m].x * realtemp - 
                  L[n][l][m].y * imgtemp);
               fb += 2.0 * (Real) m * (L[n][l][m].x * Y_LJf[n][l][m].y +
                  L[n][l][m].y * Y_LJf[n][l][m].x);   
            } /* else (m != 0) */

         } /* for m */
      } /* for l */
   } /* for n */

   cosalpha = cos(sv.a);
   sinalpha = sin(sv.a);
   cosbeta = cos(sv.b);
   sinbeta = sin(sv.b);

   if (sinalpha != 0.0) {
      fb *= 1.0/sinalpha;
   } /* if */

   else {
      fb = 0.0;
      makeYIIforce0(p, sv);
      for (n=1; n < p; n++) {
         for (l=0; l <= n; l++) {
            for (m = (0x0001 & (n+l)); m <= (n-l); m += 2) {
               if (m != 0) {
                  fb += 2.0 * (Real) m * (L[n][l][m].x * Y_LJf[n][l][m].y +
                     L[n][l][m].y * Y_LJf[n][l][m].x);
               } /* if */
            } /* for m */
         } /* for l */
      } /* for n */
   } /* else */


   *potp = b_lj * pot;
   fv->x = b_lj * (fr*sinalpha*cosbeta + fa*cosalpha*cosbeta - fb*sinbeta); 
   fv->y = b_lj * (fr*sinalpha*sinbeta + fa*cosalpha*sinbeta + fb*cosbeta);
   fv->z = b_lj * (fr*cosalpha - fa*sinalpha);

} /* Force_LJ */


/*****************************************************************
 *
 *  makeYI() - computes and updates the Y_LJ[][][] array
 *
 */

void makeYI(
   int       p,
   SphVector sv )
{
   int a, b, i, l, m, n, mstart;
   Real rinv, rtemp;
   Real stheta, stheta2, sthetaM1, sthetaM2;


   Gegenbauer(p, cos(sv.a));
   Fourier_LJ(p, sv.b);

   if (sv.r == 0.0) {
      fprintf(stderr, "zero radius passed to MakeYI\n");
      exit(0);
   }

   rinv = 1.0/sv.r;
   rtemp = rinv*rinv*rinv;
   rtemp *= rtemp;
   stheta = sin(sv.a);
   stheta2 = stheta*stheta;

   for (n=0; n < p; n++) {
      for (l=0; l <= n; l++) {

         mstart = 0x0001 & (n+l);
         if (mstart)
            sthetaM1 = stheta;
         else
            sthetaM1 = 1.0;

         for (m = mstart; m <= (n-l); m += 2) {
            a = (n-l-m)/2;
            b = (n-l+m)/2;
            Y_LJ[n][l][m].x = 0.0;
            sthetaM2 = 1.0;

            for (i=0; i <= a; i++) {
               Y_LJ[n][l][m].x +=  sthetaM2 *
                  GegPoly[l][b+i] * YIcoeff[b][a][i];
               sthetaM2 *= stheta2;
            } /* for i */

            Y_LJ[n][l][m].x *= sthetaM1 * rtemp / A_LJ[n][l][m];
            Y_LJ[n][l][m].y = Y_LJ[n][l][m].x * Yxy[m].y;
            Y_LJ[n][l][m].x *= Yxy[m].x;
            sthetaM1 *= stheta2;
         } /* for m */
      } /* for l */

      rtemp *= rinv;

   } /* for n */
} /* make YI */


/*****************************************************************
 *
 *  makeYII() - computes and updates the Y_LJ[][][] array
 *
 */

void makeYII(
   int       p,
   SphVector sv )
{
   int l, m, n;
   Real rtemp;
   Real sina, cosa;
   Real sinaM, cosaM, sinainv;


   Fourier_LJ(p, sv.b);

   sina = sin(sv.a);
   if (sina != 0.0) {
      sinainv = 1.0/sina;
      cosa = cos(sv.a);

      rtemp = 1.0;
   
      for (n=0; n < p; n++) {
         sinaM = pow(sina, (Real) n);
         cosaM = 1.0;
         for (l=0; l <= n; l++) {
            for (m = (0x0001 & n+l); m <= (n-l); m += 2) {
               Y_LJ[n][l][m].x = cosaM * sinaM * rtemp * A_LJ[n][l][m];
               Y_LJ[n][l][m].y = Y_LJ[n][l][m].x * -Yxy[m].y;
               Y_LJ[n][l][m].x *= Yxy[m].x;
            } /* for m */
            sinaM *= sinainv;
            cosaM *= cosa;
         } /* for l */
         rtemp *= sv.r;   
      } /* for n */
   } /* if sina != 0.0 */

   else {
      cosa = cos(sv.a);
      rtemp = 1.0;
  
      for (n=0; n < p; n++) {
         cosaM = 1.0;
         for (l=0; l <= n; l++) {
            for (m = (0x0001 & n+l); m <= (n-l); m += 2) {

               if (n == l) {
                  Y_LJ[n][l][m].x = cosaM * rtemp * A_LJ[n][l][m];
                  Y_LJ[n][l][m].y = Y_LJ[n][l][m].x * -Yxy[m].y;
                  Y_LJ[n][l][m].x *= Yxy[m].x;
               } /* if */

               else {
                  Y_LJ[n][l][m].x = 0.0;
                  Y_LJ[n][l][m].y = 0.0;
               } /* else */   

            } /* for m */

            cosaM *= cosa;

         } /* for l */

         rtemp *= sv.r;

      } /* for n */
   } /* else */
} /* makeYII */


/****************************************************************
 *
 *  makeYIIforce() - computes and updates the Y_LJ[][][] array
 *
 */


void makeYIIforce(
   int       p,
   SphVector sv )
{
   int l, m, n;
   Real rtemp;
   Real sina, cosa;
   Real sinaM, cosaM, sinainv;


   Fourier_LJ(p, sv.b);

   sina = sin(sv.a);
   cosa = cos(sv.a);

   if (sina == 0.0) {
      sina = 0.0;
      if (cosa > 0.0)
	 cosa = 1.0;
      else
	 cosa = -1.0;
   } /* if sina */

   if (sina != 0.0) {

      sinainv = 1.0/sina;

      Y_LJf[0][0][0].x = Yxy[0].x;
      Y_LJf[0][0][0].y = 0.0;

      rtemp = 1.0;

      for (n=1; n < p; n++) {
         sinaM = pow(sina, (Real) n);
         cosaM = 1.0;
         for (l=0; l <= n; l++) {
            for (m = 0; m <= n-l+1; m++) {
               Y_LJf[n][l][m].x = cosaM * sinaM * rtemp;
               Y_LJf[n][l][m].y = Y_LJf[n][l][m].x * -Yxy[m].y;
               Y_LJf[n][l][m].x *= Yxy[m].x;
            } /* for m */
            cosaM *= cosa;
            sinaM *= sinainv;
         } /* for l */
         rtemp *= sv.r;
      } /* for n */
   } /* if sina */

   else {
      Y_LJf[0][0][0].x = Yxy[0].x;
      Y_LJf[0][0][0].y = 0.0;

      rtemp = 1.0;

      for (n=1; n < p; n++) {
         cosaM = 1.0;
         for (l=0; l <= n; l++) {
            for (m = 0; m <= n-l+1; m++) {
               if (n == l) {
                  Y_LJf[n][l][m].x = cosaM * rtemp;
                  Y_LJf[n][l][m].y = Y_LJf[n][l][m].x * -Yxy[m].y;
                  Y_LJf[n][l][m].x *= Yxy[m].x;
               } /* if n */
               else {
                  Y_LJf[n][l][m].x = 0.0;
                  Y_LJf[n][l][m].y = 0.0;
               } /* else */
            } /* for m */
            cosaM *= cosa;
         } /* for l */
         rtemp *= sv.r;
      } /* for n */
   } /* else */
   
} /* makeYIIforce */


/*****************************************************************
 *
 *  makeYIIforce0() - computes and updates the Y_LJ[][][] array
 *
 */

void makeYIIforce0(
   int       p,
   SphVector sv )
{
   int l, m, n;
   Real rtemp;
   Real sina, cosa;
   Real sinaM, cosaM, sinainv;


   Fourier_LJ(p, sv.b);

   sina = sin(sv.a);
   cosa = cos(sv.a);

   if (sina == 0.0 ) {
      sina = 0.0;
      if (cosa > 0.0)
	cosa = 1.0;
      else
	cosa = -1.0;
   } /* if sina */

   if (sina != 0.0) {

      sinainv = 1.0/sina;

      Y_LJf[0][0][0].x = Yxy[0].x;
      Y_LJf[0][0][0].y = 0.0;

      rtemp = 1.0;

      for (n=1; n < p; n++) {
         sinaM = pow(sina, (Real) n);
         cosaM = 1.0;
         for (l=0; l <= n; l++) {
            for (m = 0; m <= n-l+1; m++) {
               Y_LJf[n][l][m].x = cosaM * sinaM * rtemp;
               Y_LJf[n][l][m].y = Y_LJf[n][l][m].x * -Yxy[m].y;
               Y_LJf[n][l][m].x *= Yxy[m].x;
            } /* for m */
            cosaM *= cosa;
            sinaM *= sinainv;
         } /* for l */
         rtemp *= sv.r;
      } /* for n */
   } /* if sina */

   else {

      Y_LJf[0][0][0].x = Yxy[0].x;
      Y_LJf[0][0][0].y = 0.0;

      rtemp = 1.0;

      for (n=1; n < p; n++) {
         cosaM = 1.0;
         for (l=0; l <= n; l++) {
            for (m = 0; m <= n-l+1; m++) {
               if ((n-l) == 1) {
                  Y_LJf[n][l][m].x = cosaM * rtemp;
                  Y_LJf[n][l][m].y = Y_LJf[n][l][m].x * -Yxy[m].y;
                  Y_LJf[n][l][m].x *= Yxy[m].x;
               } /* if */
               else {
                  Y_LJf[n][l][m].x = 0.0;
                  Y_LJf[n][l][m].y = 0.0;
               } /* else */
            } /* for m */
            cosaM *= cosa;
         } /* for l */
         rtemp *= sv.r;
      } /* for n */
   } /* else */

} /* makeYIIforce0 */


/*****************************************************************
 *
 * Gegenbauer() - Makes a Gegenbauer array shifted in the lambda 
 *   index by 3
 *
 */

void Gegenbauer(
   int    p,
   Real theta )
{
   int l, n;


   if (GegPoly == NULL) {
      fprintf(stderr,"Null pointer passed to Gegenbauer subroutine\n");
      exit(0);
   }

   for (l=0; l < p; l++)
      GegPoly[0][l] = 1.0;

   for (l=0; l < p-1; l++)
      GegPoly[1][l] = 2.0 * (Real) (l+3) * theta;
  
   for (l=0; l < p-2; l++)
      GegPoly[2][l] =  GegPoly[1][l] * (Real) (l+4) * theta - (Real) (l+3);

   for (n=3; n < p; n++) {
      for (l=0; l < p-n; l++) {
         GegPoly[n][l] = (2.0 * (Real)(l+3)/(Real) n) *
            (theta * GegPoly[n-1][l+1] - GegPoly[n-2][l+1]);
      } /* for l */
   } /* for n */

} /* Gegenbauer */


/****************************************************************
 *
 * LJinit() - initialize constants used in LJ computations
 *
 */

void LJinit( int p )
{
   int i, l, m, n, a, b;
   Real *factorial, *scratch;
   Complex *scratchC;


   /*
    * allocate and initialize space for Gegenbauer polynomial
    */

   scratch = (Real *) malloc(((p * (p+1))/2) * sizeof(Real));
   GegPoly = (Real **) malloc(p * sizeof(Real *));
   GegPoly[0] = scratch;
   scratch += p;
   for (i=1; i < p; i++) {
      GegPoly[i] = scratch;
      scratch += (p-i);
   }
   scratch = &GegPoly[0][0];
   for (n=0; n < (p * (p+1))/2; n++)
      scratch[n] = 0.0;

   /*
    * allocate and initialize Y_LJ array
    */

   scratchC = (Complex *) malloc(((p * (p+1) * (p+2))/6) * sizeof(Complex));
   Y_LJ = (Complex ***) malloc( p * sizeof(Complex **));
   for (n = 0; n < p; n++) {
      Y_LJ[n] = (Complex **) malloc((n+1) * sizeof(Complex *));
      for (l = 0; l <= n; l++) {
         Y_LJ[n][l] = scratchC;
         scratchC += n-l+1;
      } /* for l */
   } /* for n */
   scratch = &Y_LJ[0][0][0].x;
   for (n=0; n < (p * (p+1) * (p+2))/3; n++)
      scratch[n] = 0.0;

   /*
    * allocate and initialize L array
    */

   scratchC = (Complex *) malloc(((p * (p+1) * (p+2))/6) * sizeof(Complex));
   L = (Complex ***) malloc( p * sizeof(Complex **));
   for (n = 0; n < p; n++) {
      L[n] = (Complex **) malloc((n+1) * sizeof(Complex *));
      for (l = 0; l <= n; l++) {
         L[n][l] = scratchC;
         scratchC += n-l+1;
      } /* for l */
   } /* for n */
   scratch = &L[0][0][0].x;
   for (n=0; n < (p * (p+1) * (p+2))/3; n++)
      scratch[n] = 0.0;

   /*
    * allocate and initialize Y_LJF array
    */

   scratchC = (Complex *) malloc(p * p * p * sizeof(Complex));
   Y_LJf = (Complex ***) malloc( p * sizeof(Complex **));
   for (n = 0; n < p; n++) {
      Y_LJf[n] = (Complex **) malloc(p * sizeof(Complex *));
      for (l = 0; l < p; l++) {
         Y_LJf[n][l] = scratchC;
         scratchC += p;
      } /* for l */
   } /* for n */
   scratch = &Y_LJf[0][0][0].x;
   for (n=0; n < p * p * p * 2; n++)
      scratch[n] = 0.0;

   /*
    * allocate and initialize Yxy array
    */

   Yxy = (Complex *) malloc((p+1) * sizeof(Complex));


   /*
    * allocate and compute factorial array
    */

   factorial = (Real *) malloc(2*(p+1)*sizeof(Real));
   factorial[0] = 1.0;
   for(n=1;n < 2*(p+1); n++) {
      factorial[n] = (Real) n * factorial[n-1];
   }


   /*
    * allocate and compute A_LJ array
    */

   scratch = (Real *) malloc( ((p * (p+1) * (p+2))/6) * sizeof(Real));
   A_LJ = (Real ***) malloc( (p+1) * sizeof(Real **));
   for (n = 0; n < p; n++) {
      A_LJ[n] = (Real **) malloc((n+1) * sizeof(Real *));
      for (l = 0; l <= n; l++) {
         A_LJ[n][l] = scratch;
         scratch += n-l+1;
      } /* for l */
   } /* for n */

   for (n=0; n < p; n++) {
      for (l=0; l <= n; l++) {
         for (m = 0x0001 & (n+l); m <= (n-l); m += 2) {
            A_LJ[n][l][m] = pow(-1.0, (Real) n) /
               (factorial[(n-l-m)/2] * factorial[(n-l+m)/2] * factorial[l]);
         } /* for m */
      } /* for l */
   } /* for n */

   scratch = (Real *) malloc( ((p * (p+1) * (p+2))/6) * sizeof(Real));
   A_LJi = (Real ***) malloc( (p+1) * sizeof(Real **));
   for (n = 0; n < p; n++) {
      A_LJi[n] = (Real **) malloc((n+1) * sizeof(Real *));
      for (l = 0; l <= n; l++) {
         A_LJi[n][l] = scratch;
         scratch += n-l+1;
      } /* for l */
   } /* for n */


   /*
    * allocate and compute A_LJi (inverse) array
    */

   for (n=0; n < p; n++) {
      for (l=0; l <= n; l++) {
         for (m = 0x0001 & (n+l); m <= (n-l); m += 2) {
            A_LJi[n][l][m] = 1.0/A_LJ[n][l][m];
         } /* for m */
      } /* for l */
   } /* for n */


   /*
    * allocate and compute VIcoeff array
    */

   scratch = (Real *) malloc( ((p * (p+1) * (p+2))/6) * sizeof(Real));
   YIcoeff = (Real ***) malloc( p * sizeof(Real **));
   for (n = 0; n < p; n++) {
      YIcoeff[n] = (Real **) malloc((n+1) * sizeof(Real *));
      for (l = 0; l <= n; l++) {
         YIcoeff[n][l] = scratch;
         scratch += l+1;
      } /* for m */
   } /* for n */

   for (b=0; b < p; b++)  {
      for (a=0; a <= b; a++) {
         for (i=0; i <= a; i++)  {
            YIcoeff[b][a][i] = pow(-1.0, (Real) (a+i)) * factorial[b+i+2] /
               (2.0 * factorial[i] * factorial[a-i] * factorial[b+i-a]);
         } /* for i */
      } /* for a */
   } /* for b */

   /*
    * free up the factorial array.  we no longer need it.
    */

   free(factorial);

} /* LJinit */


/****************************************************************
 *
 * LJcleanup() - deallocate constants used in LJ computations
 *
 */

void LJcleanup( int p )
{
   int i;

   free( GegPoly[0] );
   free( GegPoly );


   free( Y_LJ[0][0] );
   for ( i=0; i<p; i++ ) {
      free( Y_LJ[i] );
   }
   free( Y_LJ );

   free( L[0][0] );
   for ( i=0; i<p; i++ ) {
      free( L[i] );
   }
   free( L );

   free( Y_LJf[0][0] );
   for ( i=0; i<p; i++ ) {
      free( Y_LJf[i] );
   }
   free( Y_LJf );

   free( Yxy );

   free( A_LJ[0][0] );
   for ( i=0; i<p; i++ ) {
      free( A_LJ[i] );
   }
   free( A_LJ );

   free( A_LJi[0][0] );
   for ( i=0; i<p; i++ ) {
      free( A_LJi[i] );
   }
   free( A_LJi );


   free( YIcoeff[0][0] );
   for ( i=0; i<p; i++ ) {
      free( YIcoeff[i] );
   }
   free( YIcoeff );


} /* LJcleanup */


/****************************************************************
 *
 *  LJsize() - returns size of LJ multipole array in number of
 *    complex data elements.
 *
 */

int LJsize( int p )
{
   return (p * (p+1) * (p+2))/6;
}


/****************************************************************
 *
 *  eval_mpoleLJ - returns potential from an LJ multipole 
 *    evaluated at the point specified by the vector v.
 *
 */

Real eval_mpotLJ(
   MtypeLJ M,
   int     p,
   Vector  v )
{
   int       l, m, n;
   Real      pot;
   SphVector sv;


   Cart2Sph(v, &sv);

   makeYI(p, sv);

   pot = Y_LJ[0][0][0].x * M[0][0][0].x;

   for (n=1; n < p; n++) {
      for (l=0; l <= n; l++) {
         for (m = (0x0001 & (n+l)); m <= (n-l); m += 2) {
            if (m == 0) {
               pot += Y_LJ[n][l][0].x * M[n][l][0].x;
            } /* if m == 0 */
            else {
               pot += 2.0 * (Y_LJ[n][l][m].x * M[n][l][m].x -
                  Y_LJ[n][l][m].y * M[n][l][m].y);
            } /* else */
         } /* for m */
      } /* for l */
   } /* for n */

   return pot;

} /* eval_mpotC */


/****************************************************************
 *
 * LJMclear - clear the LJ multipole expansion
 *
 */

void LJMclear(
   MtypeLJ M,
   int     p )
{   
   int n;
   Real *scratch;


   scratch = &M[0][0][0].x;
   for (n=0; n < (p * (p+1) * (p+2))/3; n++)
      scratch[n] = 0.0;

} /* LJMclear */


/****************************************************************
*
*  LJMsum(M1,M2,p) -
*     adds the multipole expansion M1 (of size p) to M2
*
*/

void LJMsum(
   MtypeLJ M1,
   MtypeLJ M2,
   int p )
{

   int i;                /* loop counters */
   Complex *scratch1;    /* temp pointer */
   Complex *scratch2;    /* temp pointer */


   scratch1 = &(M1[0][0][0]);
   scratch2 = &(M2[0][0][0]);

   for ( i=0; i < (p*(p+1)*(p+2))/6; i++) {
      scratch2[i].x += scratch1[i].x;
      scratch2[i].y += scratch1[i].y;
   } /* for i */

} /* LJMsum */


/****************************************************************
 *
 *  copyYi() - precompute the Y_LJ transfer matrix for a given
 *  vector, v, and expansion size, p, returning the matrix in Yout.
 *
 */

void copyYI(
   MtypeLJ Yout,
   int p,
   Vector v )
{
   int n, l, m;
   SphVector sv;


   Cart2Sph(v, &sv);

   makeYI(p, sv);

   for (n=0; n < p; n++) {
      for (l=0; l <= n; l++)  {
         for (m = (0x0001 & n+l); m <= (n-l); m += 2) {
            Yout[n][l][m].x = Y_LJ[n][l][m].x;
            Yout[n][l][m].y = Y_LJ[n][l][m].y;
         } /* for m */
      } /* for l */
   } /* for n */
} /* copy YI */


/****************************************************************
 *
 * Fourier_LJ() - Fourier series components of angle b out to p terms
 *
 */

void Fourier_LJ(
   int    p,
   Real   b )
{
   int m;


   if (Yxy == NULL) {
      fprintf(stderr, "Fourier called with null pointer to array\n");
      exit(0);
   }

   for (m=0; m <= p; m++) {
      Yxy[m].x = cos((Real) m * b);
      Yxy[m].y = sin((Real) m * b);
   } /* for m */

} /* Fourier */



