/* f2c.h  --  Standard Fortran to C header file */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."

	- From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */

#include <math.h>

#ifdef _LONG64_
typedef int integer;
typedef int logical;
#else
typedef long int integer;
typedef long int logical;
#endif
typedef char *address;
typedef short int shortint;
typedef struct { float r, i; } complex;
typedef struct { double r, i; } doublecomplex;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;

#define TRUE_ (1)
#define FALSE_ (0)

/* Extern is for use with -E */
#ifndef Extern
#define Extern static
#endif
#define extern static

/* I/O stuff */

#ifdef f2c_i2
/* for -i2 */
typedef short flag;
typedef short ftnlen;
typedef short ftnint;
#else
#ifdef _LONG64_
typedef int flag;
typedef int ftnlen;
typedef int ftnint;
#else
typedef long flag;
typedef long ftnlen;
typedef long ftnint;
#endif
#endif

/*external read, write*/
typedef struct
{	flag cierr;
	ftnint ciunit;
	flag ciend;
	char *cifmt;
	ftnint cirec;
} cilist;

/*internal read, write*/
typedef struct
{	flag icierr;
	char *iciunit;
	flag iciend;
	char *icifmt;
	ftnint icirlen;
	ftnint icirnum;
} icilist;

/*open*/
typedef struct
{	flag oerr;
	ftnint ounit;
	char *ofnm;
	ftnlen ofnmlen;
	char *osta;
	char *oacc;
	char *ofm;
	ftnint orl;
	char *oblnk;
} olist;

/*close*/
typedef struct
{	flag cerr;
	ftnint cunit;
	char *csta;
} cllist;

/*rewind, backspace, endfile*/
typedef struct
{	flag aerr;
	ftnint aunit;
} alist;

/* inquire */
typedef struct
{	flag inerr;
	ftnint inunit;
	char *infile;
	ftnlen infilen;
	ftnint	*inex;	/*parameters in standard's order*/
	ftnint	*inopen;
	ftnint	*innum;
	ftnint	*innamed;
	char	*inname;
	ftnlen	innamlen;
	char	*inacc;
	ftnlen	inacclen;
	char	*inseq;
	ftnlen	inseqlen;
	char 	*indir;
	ftnlen	indirlen;
	char	*infmt;
	ftnlen	infmtlen;
	char	*inform;
	ftnint	informlen;
	char	*inunf;
	ftnlen	inunflen;
	ftnint	*inrecl;
	ftnint	*innrec;
	char	*inblank;
	ftnlen	inblanklen;
} inlist;

#define VOID void

union Multitype {	/* for multiple entry points */
	shortint h;
	integer i;
	float r;
	double d;
	complex c;
	doublecomplex z;
	};

typedef union Multitype Multitype;

struct Vardesc {	/* for Namelist */
	char *name;
	char *addr;
	ftnlen *dims;
	int  type;
	};
typedef struct Vardesc Vardesc;

struct Namelist {
	char *name;
	Vardesc **vars;
	int nvars;
	};
typedef struct Namelist Namelist;

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (double)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (double)min(a,b)
#define dmax(a,b) (double)max(a,b)

/* procedure parameter types for -A and -C++ */

#define F2C_proc_par_types 1
#ifdef __cplusplus
typedef int /* Unknown procedure type */ (*U_fp)(...);
typedef shortint (*J_fp)(...);
typedef integer (*I_fp)(...);
typedef float (*R_fp)(...);
typedef double (*D_fp)(...), (*E_fp)(...);
typedef /* Complex */ VOID (*C_fp)(...);
typedef /* Double Complex */ VOID (*Z_fp)(...);
typedef logical (*L_fp)(...);
typedef shortlogical (*K_fp)(...);
typedef /* Character */ VOID (*H_fp)(...);
typedef /* Subroutine */ int (*S_fp)(...);
#else
typedef int /* Unknown procedure type */ (*U_fp)();
typedef shortint (*J_fp)();
typedef integer (*I_fp)();
typedef float (*R_fp)();
typedef double (*D_fp)(), (*E_fp)();
typedef /* Complex */ VOID (*C_fp)();
typedef /* Double Complex */ VOID (*Z_fp)();
typedef logical (*L_fp)();
typedef shortlogical (*K_fp)();
typedef /* Character */ VOID (*H_fp)();
typedef /* Subroutine */ int (*S_fp)();
#endif
/* E_fp is for float functions when -R is not specified */
typedef VOID C_f;	/* complex function */
typedef VOID H_f;	/* character function */
typedef VOID Z_f;	/* double complex function */
typedef double E_f;	/* float function with -R not specified */

/* undef any lower-case symbols that your C compiler predefines, e.g.: */

#ifndef Skip_f2c_Undefs
#undef cray
#undef gcos
#undef mc68010
#undef mc68020
#undef mips
#undef pdp11
#undef sgi
#undef sparc
#undef sun
#undef sun2
#undef sun3
#undef sun4
#undef u370
#undef u3b
#undef u3b2
#undef u3b5
#undef unix
#undef vax
#endif

/* End of f2c.h */



/* Subroutine */ int dsyev_(char *jobz, char *uplo, integer *n, double *a,
	 integer *lda, double *w, double *work, integer *lwork, 
	integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSYEV computes all eigenvalues and, optionally, eigenvectors of a   
    float symmetric matrix A.   

    Arguments   
    =========   

    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the   
            leading N-by-N upper triangular part of A contains the   
            upper triangular part of the matrix A.  If UPLO = 'L',   
            the leading N-by-N lower triangular part of A contains   
            the lower triangular part of the matrix A.   
            On exit, if JOBZ = 'V', then if INFO = 0, A contains the   
            orthonormal eigenvectors of the matrix A.   
            If JOBZ = 'N', then on exit the lower triangle (if UPLO='L') 
  
            or the upper triangle (if UPLO='U') of A, including the   
            diagonal, is destroyed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    W       (output) DOUBLE PRECISION array, dimension (N)   
            If INFO = 0, the eigenvalues in ascending order.   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of the array WORK.  LWORK >= max(1,3*N-1).   
            For optimal efficiency, LWORK >= (NB+2)*N,   
            where NB is the blocksize for DSYTRD returned by ILAENV.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the algorithm failed to converge; i   
                  off-diagonal elements of an intermediate tridiagonal   
                  form did not converge to zero.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__0 = 0;
    static double c_b12 = 1.;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    double d__1;
    /* Local variables */
    static integer inde;
    static double anrm;
    static integer imax;
    static double rmin, rmax;
    static integer lopt;
    extern /* Subroutine */ int dscal_(integer *, double *, double *, 
	    integer *);
    static double sigma;
    extern logical lsame_(char *, char *);
    static integer iinfo;
    static logical lower, wantz;
    extern double dlamch_(char *);
    static integer iscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    double *, double *, integer *, integer *, double *, 
	    integer *, integer *);
    static double safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static double bignum;
    static integer indtau;
    extern /* Subroutine */ int dsterf_(integer *, double *, double *,
	     integer *);
    extern double dlansy_(char *, char *, integer *, double *, 
	    integer *, double *);
    static integer indwrk;
    extern /* Subroutine */ int dorgtr_(char *, integer *, double *, 
	    integer *, double *, double *, integer *, integer *), dsteqr_(char *, integer *, double *, double *, 
	    double *, integer *, double *, integer *), 
	    dsytrd_(char *, integer *, double *, integer *, double *, 
	    double *, double *, double *, integer *, integer *);
    static integer llwork;
    static double smlnum, eps;



#define W(I) w[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    wantz = lsame_(jobz, "V");
    lower = lsame_(uplo, "L");

    *info = 0;
    if (! (wantz || lsame_(jobz, "N"))) {
	*info = -1;
    } else if (! (lower || lsame_(uplo, "U"))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n * 3 - 1;
	if (*lwork < max(i__1,i__2)) {
	    *info = -8;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYEV ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	WORK(1) = 1.;
	return 0;
    }

    if (*n == 1) {
	W(1) = A(1,1);
	WORK(1) = 3.;
	if (wantz) {
	    A(1,1) = 1.;
	}
	return 0;
    }

/*     Get machine constants. */

    safmin = dlamch_("Safe minimum");
    eps = dlamch_("Precision");
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

    anrm = dlansy_("M", uplo, n, &A(1,1), lda, &WORK(1));
    iscale = 0;
    if (anrm > 0. && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	dlascl_(uplo, &c__0, &c__0, &c_b12, &sigma, n, n, &A(1,1), lda, 
		info);
    }

/*     Call DSYTRD to reduce symmetric matrix to tridiagonal form. */

    inde = 1;
    indtau = inde + *n;
    indwrk = indtau + *n;
    llwork = *lwork - indwrk + 1;
    dsytrd_(uplo, n, &A(1,1), lda, &W(1), &WORK(inde), &WORK(indtau), &
	    WORK(indwrk), &llwork, &iinfo);
    lopt = (integer) ((*n << 1) + WORK(indwrk));

/*     For eigenvalues only, call DSTERF.  For eigenvectors, first call   
       DORGTR to generate the orthogonal matrix, then call DSTEQR. */

    if (! wantz) {
	dsterf_(n, &W(1), &WORK(inde), info);
    } else {
	dorgtr_(uplo, n, &A(1,1), lda, &WORK(indtau), &WORK(indwrk), &
		llwork, &iinfo);
	dsteqr_(jobz, n, &W(1), &WORK(inde), &A(1,1), lda, &WORK(indtau),
		 info);
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

    if (iscale == 1) {
	if (*info == 0) {
	    imax = *n;
	} else {
	    imax = *info - 1;
	}
	d__1 = 1. / sigma;
	dscal_(&imax, &d__1, &W(1), &c__1);
    }

/*     Set WORK(1) to optimal workspace size.   

   Computing MAX */
    i__1 = *n * 3 - 1;
    WORK(1) = (double) max(i__1,lopt);

    return 0;

/*     End of DSYEV */

} /* dsyev_ */


/* Subroutine */
static  int dlae2_(double *a, double *b, double *c, 
	double *rt1, double *rt2)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix   
       [  A   B  ]   
       [  B   C  ].   
    On return, RT1 is the eigenvalue of larger absolute value, and RT2   
    is the eigenvalue of smaller absolute value.   

    Arguments   
    =========   

    A       (input) DOUBLE PRECISION   
            The (1,1) element of the 2-by-2 matrix.   

    B       (input) DOUBLE PRECISION   
            The (1,2) and (2,1) elements of the 2-by-2 matrix.   

    C       (input) DOUBLE PRECISION   
            The (2,2) element of the 2-by-2 matrix.   

    RT1     (output) DOUBLE PRECISION   
            The eigenvalue of larger absolute value.   

    RT2     (output) DOUBLE PRECISION   
            The eigenvalue of smaller absolute value.   

    Further Details   
    ===============   

    RT1 is accurate to a few ulps barring over/underflow.   

    RT2 may be inaccurate if there is massive cancellation in the   
    determinant A*C-B*B; higher precision or correctly rounded or   
    correctly truncated arithmetic would be needed to compute RT2   
    accurately in all cases.   

    Overflow is possible only if RT1 is within a factor of 5 of overflow. 
  
    Underflow is harmless if the input data is 0 or exceeds   
       underflow_threshold / macheps.   

   ===================================================================== 
  


       Compute the eigenvalues */
    /* System generated locals */
    double d__1;
    /* Local variables */
    static double acmn, acmx, ab, df, tb, sm, rt, adf;


    sm = *a + *c;
    df = *a - *c;
    adf = abs(df);
    tb = *b + *b;
    ab = abs(tb);
    if (abs(*a) > abs(*c)) {
	acmx = *a;
	acmn = *c;
    } else {
	acmx = *c;
	acmn = *a;
    }
    if (adf > ab) {
/* Computing 2nd power */
	d__1 = ab / adf;
	rt = adf * sqrt(d__1 * d__1 + 1.);
    } else if (adf < ab) {
/* Computing 2nd power */
	d__1 = adf / ab;
	rt = ab * sqrt(d__1 * d__1 + 1.);
    } else {

/*        Includes case AB=ADF=0 */

	rt = ab * sqrt(2.);
    }
    if (sm < 0.) {
	*rt1 = (sm - rt) * .5;

/*        Order of execution important.   
          To get fully accurate smaller eigenvalue,   
          next line needs to be executed in higher precision. */

	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else if (sm > 0.) {
	*rt1 = (sm + rt) * .5;

/*        Order of execution important.   
          To get fully accurate smaller eigenvalue,   
          next line needs to be executed in higher precision. */

	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else {

/*        Includes case RT1 = RT2 = 0 */

	*rt1 = rt * .5;
	*rt2 = rt * -.5;
    }
    return 0;

/*     End of DLAE2 */

} /* dlae2_ */


/* Subroutine */
static  int dlaev2_(double *a, double *b, double *c, 
	double *rt1, double *rt2, double *cs1, double *sn1)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix   
       [  A   B  ]   
       [  B   C  ].   
    On return, RT1 is the eigenvalue of larger absolute value, RT2 is the 
  
    eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right 
  
    eigenvector for RT1, giving the decomposition   

       [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]   
       [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].   

    Arguments   
    =========   

    A       (input) DOUBLE PRECISION   
            The (1,1) element of the 2-by-2 matrix.   

    B       (input) DOUBLE PRECISION   
            The (1,2) element and the conjugate of the (2,1) element of   
            the 2-by-2 matrix.   

    C       (input) DOUBLE PRECISION   
            The (2,2) element of the 2-by-2 matrix.   

    RT1     (output) DOUBLE PRECISION   
            The eigenvalue of larger absolute value.   

    RT2     (output) DOUBLE PRECISION   
            The eigenvalue of smaller absolute value.   

    CS1     (output) DOUBLE PRECISION   
    SN1     (output) DOUBLE PRECISION   
            The vector (CS1, SN1) is a unit right eigenvector for RT1.   

    Further Details   
    ===============   

    RT1 is accurate to a few ulps barring over/underflow.   

    RT2 may be inaccurate if there is massive cancellation in the   
    determinant A*C-B*B; higher precision or correctly rounded or   
    correctly truncated arithmetic would be needed to compute RT2   
    accurately in all cases.   

    CS1 and SN1 are accurate to a few ulps barring over/underflow.   

    Overflow is possible only if RT1 is within a factor of 5 of overflow. 
  
    Underflow is harmless if the input data is 0 or exceeds   
       underflow_threshold / macheps.   

   ===================================================================== 
  


       Compute the eigenvalues */
    /* System generated locals */
    double d__1;
    /* Local variables */
    static double acmn, acmx, ab, df, cs, ct, tb, sm, tn, rt, adf, acs;
    static integer sgn1, sgn2;


    sm = *a + *c;
    df = *a - *c;
    adf = abs(df);
    tb = *b + *b;
    ab = abs(tb);
    if (abs(*a) > abs(*c)) {
	acmx = *a;
	acmn = *c;
    } else {
	acmx = *c;
	acmn = *a;
    }
    if (adf > ab) {
/* Computing 2nd power */
	d__1 = ab / adf;
	rt = adf * sqrt(d__1 * d__1 + 1.);
    } else if (adf < ab) {
/* Computing 2nd power */
	d__1 = adf / ab;
	rt = ab * sqrt(d__1 * d__1 + 1.);
    } else {

/*        Includes case AB=ADF=0 */

	rt = ab * sqrt(2.);
    }
    if (sm < 0.) {
	*rt1 = (sm - rt) * .5;
	sgn1 = -1;

/*        Order of execution important.   
          To get fully accurate smaller eigenvalue,   
          next line needs to be executed in higher precision. */

	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else if (sm > 0.) {
	*rt1 = (sm + rt) * .5;
	sgn1 = 1;

/*        Order of execution important.   
          To get fully accurate smaller eigenvalue,   
          next line needs to be executed in higher precision. */

	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else {

/*        Includes case RT1 = RT2 = 0 */

	*rt1 = rt * .5;
	*rt2 = rt * -.5;
	sgn1 = 1;
    }

/*     Compute the eigenvector */

    if (df >= 0.) {
	cs = df + rt;
	sgn2 = 1;
    } else {
	cs = df - rt;
	sgn2 = -1;
    }
    acs = abs(cs);
    if (acs > ab) {
	ct = -tb / cs;
	*sn1 = 1. / sqrt(ct * ct + 1.);
	*cs1 = ct * *sn1;
    } else {
	if (ab == 0.) {
	    *cs1 = 1.;
	    *sn1 = 0.;
	} else {
	    tn = -cs / tb;
	    *cs1 = 1. / sqrt(tn * tn + 1.);
	    *sn1 = tn * *cs1;
	}
    }
    if (sgn1 == sgn2) {
	tn = *cs1;
	*cs1 = -(*sn1);
	*sn1 = tn;
    }
    return 0;

/*     End of DLAEV2 */

} /* dlaev2_ */


static double dlamch_(char *cmach)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMCH determines double precision machine parameters.   

    Arguments   
    =========   

    CMACH   (input) CHARACTER*1   
            Specifies the value to be returned by DLAMCH:   
            = 'E' or 'e',   DLAMCH := eps   
            = 'S' or 's ,   DLAMCH := sfmin   
            = 'B' or 'b',   DLAMCH := base   
            = 'P' or 'p',   DLAMCH := eps*base   
            = 'N' or 'n',   DLAMCH := t   
            = 'R' or 'r',   DLAMCH := rnd   
            = 'M' or 'm',   DLAMCH := emin   
            = 'U' or 'u',   DLAMCH := rmin   
            = 'L' or 'l',   DLAMCH := emax   
            = 'O' or 'o',   DLAMCH := rmax   

            where   

            eps   = relative machine precision   
            sfmin = safe minimum, such that 1/sfmin does not overflow   
            base  = base of the machine   
            prec  = eps*base   
            t     = number of (base) digits in the mantissa   
            rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise   
            emin  = minimum exponent before (gradual) underflow   
            rmin  = underflow threshold - base**(emin-1)   
            emax  = largest exponent before overflow   
            rmax  = overflow threshold  - (base**emax)*(1-eps)   

   ===================================================================== 
*/
/* >>Start of File<<   
       Initialized data */
    static logical first = TRUE_;
    /* System generated locals */
    integer i__1;
    double ret_val;
    /* Builtin functions */
    extern double pow_di(double *, integer *);
    /* Local variables */
    static double base;
    static integer beta;
    static double emin, prec, emax;
    static integer imin, imax;
    static logical lrnd;
    static double rmin, rmax, t, rmach;
    extern logical lsame_(char *, char *);
    static double small, sfmin;
    extern /* Subroutine */ int dlamc2_(integer *, integer *, logical *, 
	    double *, integer *, double *, integer *, double *);
    static integer it;
    static double rnd, eps;



    if (first) {
	first = FALSE_;
	dlamc2_(&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
	base = (double) beta;
	t = (double) it;
	if (lrnd) {
	    rnd = 1.;
	    i__1 = 1 - it;
	    eps = pow_di(&base, &i__1) / 2;
	} else {
	    rnd = 0.;
	    i__1 = 1 - it;
	    eps = pow_di(&base, &i__1);
	}
	prec = eps * base;
	emin = (double) imin;
	emax = (double) imax;
	sfmin = rmin;
	small = 1. / rmax;
	if (small >= sfmin) {

/*           Use SMALL plus a bit, to avoid the possibility of rou
nding   
             causing overflow when computing  1/sfmin. */

	    sfmin = small * (eps + 1.);
	}
    }

    if (lsame_(cmach, "E")) {
	rmach = eps;
    } else if (lsame_(cmach, "S")) {
	rmach = sfmin;
    } else if (lsame_(cmach, "B")) {
	rmach = base;
    } else if (lsame_(cmach, "P")) {
	rmach = prec;
    } else if (lsame_(cmach, "N")) {
	rmach = t;
    } else if (lsame_(cmach, "R")) {
	rmach = rnd;
    } else if (lsame_(cmach, "M")) {
	rmach = emin;
    } else if (lsame_(cmach, "U")) {
	rmach = rmin;
    } else if (lsame_(cmach, "L")) {
	rmach = emax;
    } else if (lsame_(cmach, "O")) {
	rmach = rmax;
    }

    ret_val = rmach;
    return ret_val;

/*     End of DLAMCH */

} /* dlamch_ */


/* Subroutine */
static  int dlamc1_(integer *beta, integer *t, logical *rnd, logical 
	*ieee1)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC1 determines the machine parameters given by BETA, T, RND, and   
    IEEE1.   

    Arguments   
    =========   

    BETA    (output) INTEGER   
            The base of the machine.   

    T       (output) INTEGER   
            The number of ( BETA ) digits in the mantissa.   

    RND     (output) LOGICAL   
            Specifies whether proper rounding  ( RND = .TRUE. )  or   
            chopping  ( RND = .FALSE. )  occurs in addition. This may not 
  
            be a reliable guide to the way in which the machine performs 
  
            its arithmetic.   

    IEEE1   (output) LOGICAL   
            Specifies whether rounding appears to be done in the IEEE   
            'round to nearest' style.   

    Further Details   
    ===============   

    The routine is based on the routine  ENVRON  by Malcolm and   
    incorporates suggestions by Gentleman and Marovich. See   

       Malcolm M. A. (1972) Algorithms to reveal properties of   
          floating-point arithmetic. Comms. of the ACM, 15, 949-951.   

       Gentleman W. M. and Marovich S. B. (1974) More on algorithms   
          that reveal properties of floating point arithmetic units.   
          Comms. of the ACM, 17, 276-277.   

   ===================================================================== 
*/
    /* Initialized data */
    static logical first = TRUE_;
    /* System generated locals */
    double d__1, d__2;
    /* Local variables */
    static logical lrnd;
    static double a, b, c, f;
    static integer lbeta;
    static double savec;
    extern double dlamc3_(double *, double *);
    static logical lieee1;
    static double t1, t2;
    static integer lt;
    static double one, qtr;



    if (first) {
	first = FALSE_;
	one = 1.;

/*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BE
TA,   
          IEEE1, T and RND.   

          Throughout this routine  we use the function  DLAMC3  to ens
ure   
          that relevant values are  stored and not held in registers, 
 or   
          are not affected by optimizers.   

          Compute  a = 2.0**m  with the  smallest positive integer m s
uch   
          that   

             fl( a + 1.0 ) = a. */

	a = 1.;
	c = 1.;

/* +       WHILE( C.EQ.ONE )LOOP */
L10:
	if (c == one) {
	    a *= 2;
	    c = dlamc3_(&a, &one);
	    d__1 = -a;
	    c = dlamc3_(&c, &d__1);
	    goto L10;
	}
/* +       END WHILE   

          Now compute  b = 2.0**m  with the smallest positive integer 
m   
          such that   

             fl( a + b ) .gt. a. */

	b = 1.;
	c = dlamc3_(&a, &b);

/* +       WHILE( C.EQ.A )LOOP */
L20:
	if (c == a) {
	    b *= 2;
	    c = dlamc3_(&a, &b);
	    goto L20;
	}
/* +       END WHILE   

          Now compute the base.  a and c  are neighbouring floating po
int   
          numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and
 so   
          their difference is beta. Adding 0.25 to c is to ensure that
 it   
          is truncated to beta and not ( beta - 1 ). */

	qtr = one / 4;
	savec = c;
	d__1 = -a;
	c = dlamc3_(&c, &d__1);
	lbeta = (integer) (c + qtr);

/*        Now determine whether rounding or chopping occurs,  by addin
g a   
          bit  less  than  beta/2  and a  bit  more  than  beta/2  to 
 a. */

	b = (double) lbeta;
	d__1 = b / 2;
	d__2 = -b / 100;
	f = dlamc3_(&d__1, &d__2);
	c = dlamc3_(&f, &a);
	if (c == a) {
	    lrnd = TRUE_;
	} else {
	    lrnd = FALSE_;
	}
	d__1 = b / 2;
	d__2 = b / 100;
	f = dlamc3_(&d__1, &d__2);
	c = dlamc3_(&f, &a);
	if (lrnd && c == a) {
	    lrnd = FALSE_;
	}

/*        Try and decide whether rounding is done in the  IEEE  'round
 to   
          nearest' style. B/2 is half a unit in the last place of the 
two   
          numbers A and SAVEC. Furthermore, A is even, i.e. has last  
bit   
          zero, and SAVEC is odd. Thus adding B/2 to A should not  cha
nge   
          A, but adding B/2 to SAVEC should change SAVEC. */

	d__1 = b / 2;
	t1 = dlamc3_(&d__1, &a);
	d__1 = b / 2;
	t2 = dlamc3_(&d__1, &savec);
	lieee1 = t1 == a && t2 > savec && lrnd;

/*        Now find  the  mantissa, t.  It should  be the  integer part
 of   
          log to the base beta of a,  however it is safer to determine
  t   
          by powering.  So we find t as the smallest positive integer 
for   
          which   

             fl( beta**t + 1.0 ) = 1.0. */

	lt = 0;
	a = 1.;
	c = 1.;

/* +       WHILE( C.EQ.ONE )LOOP */
L30:
	if (c == one) {
	    ++lt;
	    a *= lbeta;
	    c = dlamc3_(&a, &one);
	    d__1 = -a;
	    c = dlamc3_(&c, &d__1);
	    goto L30;
	}
/* +       END WHILE */

    }

    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *ieee1 = lieee1;
    return 0;

/*     End of DLAMC1 */

} /* dlamc1_ */


/* Subroutine */
static  int dlamc2_(integer *beta, integer *t, logical *rnd, 
	double *eps, integer *emin, double *rmin, integer *emax, 
	double *rmax)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC2 determines the machine parameters specified in its argument   
    list.   

    Arguments   
    =========   

    BETA    (output) INTEGER   
            The base of the machine.   

    T       (output) INTEGER   
            The number of ( BETA ) digits in the mantissa.   

    RND     (output) LOGICAL   
            Specifies whether proper rounding  ( RND = .TRUE. )  or   
            chopping  ( RND = .FALSE. )  occurs in addition. This may not 
  
            be a reliable guide to the way in which the machine performs 
  
            its arithmetic.   

    EPS     (output) DOUBLE PRECISION   
            The smallest positive number such that   

               fl( 1.0 - EPS ) .LT. 1.0,   

            where fl denotes the computed value.   

    EMIN    (output) INTEGER   
            The minimum exponent before (gradual) underflow occurs.   

    RMIN    (output) DOUBLE PRECISION   
            The smallest normalized number for the machine, given by   
            BASE**( EMIN - 1 ), where  BASE  is the floating point value 
  
            of BETA.   

    EMAX    (output) INTEGER   
            The maximum exponent before overflow occurs.   

    RMAX    (output) DOUBLE PRECISION   
            The largest positive number for the machine, given by   
            BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point 
  
            value of BETA.   

    Further Details   
    ===============   

    The computation of  EPS  is based on a routine PARANOIA by   
    W. Kahan of the University of California at Berkeley.   

   ===================================================================== 
*/
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* Initialized data */
    static logical first = TRUE_;
    static logical iwarn = FALSE_;
    /* System generated locals */
    integer i__1;
    double d__1, d__2, d__3, d__4, d__5;
    /* Builtin functions */
    extern double pow_di(double *, integer *);
    /* Local variables */
    static logical ieee;
    static double half;
    static logical lrnd;
    static double leps, zero, a, b, c;
    static integer i, lbeta;
    static double rbase;
    static integer lemin, lemax, gnmin;
    static double small;
    static integer gpmin;
    static double third, lrmin, lrmax, sixth;
    extern /* Subroutine */ int dlamc1_(integer *, integer *, logical *, 
	    logical *);
    extern double dlamc3_(double *, double *);
    static logical lieee1;
    extern /* Subroutine */ int dlamc4_(integer *, double *, integer *), 
	    dlamc5_(integer *, integer *, integer *, logical *, integer *, 
	    double *);
    static integer lt, ngnmin, ngpmin;
    static double one, two;



    if (first) {
	first = FALSE_;
	zero = 0.;
	one = 1.;
	two = 2.;

/*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values
 of   
          BETA, T, RND, EPS, EMIN and RMIN.   

          Throughout this routine  we use the function  DLAMC3  to ens
ure   
          that relevant values are stored  and not held in registers, 
 or   
          are not affected by optimizers.   

          DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1. 
*/

	dlamc1_(&lbeta, &lt, &lrnd, &lieee1);

/*        Start to find EPS. */

	b = (double) lbeta;
	i__1 = -lt;
	a = pow_di(&b, &i__1);
	leps = a;

/*        Try some tricks to see whether or not this is the correct  E
PS. */

	b = two / 3;
	half = one / 2;
	d__1 = -half;
	sixth = dlamc3_(&b, &d__1);
	third = dlamc3_(&sixth, &sixth);
	d__1 = -half;
	b = dlamc3_(&third, &d__1);
	b = dlamc3_(&b, &sixth);
	b = abs(b);
	if (b < leps) {
	    b = leps;
	}

	leps = 1.;

/* +       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP */
L10:
	if (leps > b && b > zero) {
	    leps = b;
	    d__1 = half * leps;
/* Computing 5th power */
	    d__3 = two, d__4 = d__3, d__3 *= d__3;
/* Computing 2nd power */
	    d__5 = leps;
	    d__2 = d__4 * (d__3 * d__3) * (d__5 * d__5);
	    c = dlamc3_(&d__1, &d__2);
	    d__1 = -c;
	    c = dlamc3_(&half, &d__1);
	    b = dlamc3_(&half, &c);
	    d__1 = -b;
	    c = dlamc3_(&half, &d__1);
	    b = dlamc3_(&half, &c);
	    goto L10;
	}
/* +       END WHILE */

	if (a < leps) {
	    leps = a;
	}

/*        Computation of EPS complete.   

          Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3
)).   
          Keep dividing  A by BETA until (gradual) underflow occurs. T
his   
          is detected when we cannot recover the previous A. */

	rbase = one / lbeta;
	small = one;
	for (i = 1; i <= 3; ++i) {
	    d__1 = small * rbase;
	    small = dlamc3_(&d__1, &zero);
/* L20: */
	}
	a = dlamc3_(&one, &small);
	dlamc4_(&ngpmin, &one, &lbeta);
	d__1 = -one;
	dlamc4_(&ngnmin, &d__1, &lbeta);
	dlamc4_(&gpmin, &a, &lbeta);
	d__1 = -a;
	dlamc4_(&gnmin, &d__1, &lbeta);
	ieee = FALSE_;

	if (ngpmin == ngnmin && gpmin == gnmin) {
	    if (ngpmin == gpmin) {
		lemin = ngpmin;
/*            ( Non twos-complement machines, no gradual under
flow;   
                e.g.,  VAX ) */
	    } else if (gpmin - ngpmin == 3) {
		lemin = ngpmin - 1 + lt;
		ieee = TRUE_;
/*            ( Non twos-complement machines, with gradual und
erflow;   
                e.g., IEEE standard followers ) */
	    } else {
		lemin = min(ngpmin,gpmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else if (ngpmin == gpmin && ngnmin == gnmin) {
	    if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1) {
		lemin = max(ngpmin,ngnmin);
/*            ( Twos-complement machines, no gradual underflow
;   
                e.g., CYBER 205 ) */
	    } else {
		lemin = min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1 && gpmin == gnmin)
		 {
	    if (gpmin - min(ngpmin,ngnmin) == 3) {
		lemin = max(ngpmin,ngnmin) - 1 + lt;
/*            ( Twos-complement machines with gradual underflo
w;   
                no known machine ) */
	    } else {
		lemin = min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else {
/* Computing MIN */
	    i__1 = min(ngpmin,ngnmin), i__1 = min(i__1,gpmin);
	    lemin = min(i__1,gnmin);
/*         ( A guess; no known machine ) */
	    iwarn = TRUE_;
	}
/* **   
   Comment out this if block if EMIN is ok */
	if (iwarn) {
	    first = TRUE_;
	    printf("\n\n WARNING. The value EMIN may be incorrect:- ");
	    printf("EMIN = %8i\n",lemin);
	    printf("If, after inspection, the value EMIN looks acceptable");
            printf("please comment out \n the IF block as marked within the"); 
            printf("code of routine DLAMC2, \n otherwise supply EMIN"); 
            printf("explicitly.\n");
	}
/* **   

          Assume IEEE arithmetic if we found denormalised  numbers abo
ve,   
          or if arithmetic seems to round in the  IEEE style,  determi
ned   
          in routine DLAMC1. A true IEEE machine should have both  thi
ngs   
          true; however, faulty machines may have one or the other. */

	ieee = ieee || lieee1;

/*        Compute  RMIN by successive division by  BETA. We could comp
ute   
          RMIN as BASE**( EMIN - 1 ),  but some machines underflow dur
ing   
          this computation. */

	lrmin = 1.;
	i__1 = 1 - lemin;
	for (i = 1; i <= 1-lemin; ++i) {
	    d__1 = lrmin * rbase;
	    lrmin = dlamc3_(&d__1, &zero);
/* L30: */
	}

/*        Finally, call DLAMC5 to compute EMAX and RMAX. */

	dlamc5_(&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
    }

    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *eps = leps;
    *emin = lemin;
    *rmin = lrmin;
    *emax = lemax;
    *rmax = lrmax;

    return 0;


/*     End of DLAMC2 */

} /* dlamc2_ */


static double dlamc3_(double *a, double *b)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC3  is intended to force  A  and  B  to be stored prior to doing 
  
    the addition of  A  and  B ,  for use in situations where optimizers 
  
    might hold one of these in a register.   

    Arguments   
    =========   

    A, B    (input) DOUBLE PRECISION   
            The values A and B.   

   ===================================================================== 
*/
/* >>Start of File<<   
       System generated locals */
    double ret_val;



    ret_val = *a + *b;

    return ret_val;

/*     End of DLAMC3 */

} /* dlamc3_ */


/* Subroutine */
static  int dlamc4_(integer *emin, double *start, integer *base)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC4 is a service routine for DLAMC2.   

    Arguments   
    =========   

    EMIN    (output) EMIN   
            The minimum exponent before (gradual) underflow, computed by 
  
            setting A = START and dividing by BASE until the previous A   
            can not be recovered.   

    START   (input) DOUBLE PRECISION   
            The starting point for determining EMIN.   

    BASE    (input) INTEGER   
            The base of the machine.   

   ===================================================================== 
*/
    /* System generated locals */
    integer i__1;
    double d__1;
    /* Local variables */
    static double zero, a;
    static integer i;
    static double rbase, b1, b2, c1, c2, d1, d2;
    extern double dlamc3_(double *, double *);
    static double one;



    a = *start;
    one = 1.;
    rbase = one / *base;
    zero = 0.;
    *emin = 1;
    d__1 = a * rbase;
    b1 = dlamc3_(&d__1, &zero);
    c1 = a;
    c2 = a;
    d1 = a;
    d2 = a;
/* +    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.   
      $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP */
L10:
    if (c1 == a && c2 == a && d1 == a && d2 == a) {
	--(*emin);
	a = b1;
	d__1 = a / *base;
	b1 = dlamc3_(&d__1, &zero);
	d__1 = b1 * *base;
	c1 = dlamc3_(&d__1, &zero);
	d1 = zero;
	i__1 = *base;
	for (i = 1; i <= *base; ++i) {
	    d1 += b1;
/* L20: */
	}
	d__1 = a * rbase;
	b2 = dlamc3_(&d__1, &zero);
	d__1 = b2 / rbase;
	c2 = dlamc3_(&d__1, &zero);
	d2 = zero;
	i__1 = *base;
	for (i = 1; i <= *base; ++i) {
	    d2 += b2;
/* L30: */
	}
	goto L10;
    }
/* +    END WHILE */

    return 0;

/*     End of DLAMC4 */

} /* dlamc4_ */


/* Subroutine */
static  int dlamc5_(integer *beta, integer *p, integer *emin, 
	logical *ieee, integer *emax, double *rmax)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC5 attempts to compute RMAX, the largest machine floating-point   
    number, without overflow.  It assumes that EMAX + abs(EMIN) sum   
    approximately to a power of 2.  It will fail on machines where this   
    assumption does not hold, for example, the Cyber 205 (EMIN = -28625, 
  
    EMAX = 28718).  It will also fail if the value supplied for EMIN is   
    too large (i.e. too close to zero), probably with overflow.   

    Arguments   
    =========   

    BETA    (input) INTEGER   
            The base of floating-point arithmetic.   

    P       (input) INTEGER   
            The number of base BETA digits in the mantissa of a   
            floating-point value.   

    EMIN    (input) INTEGER   
            The minimum exponent before (gradual) underflow.   

    IEEE    (input) LOGICAL   
            A logical flag specifying whether or not the arithmetic   
            system is thought to comply with the IEEE standard.   

    EMAX    (output) INTEGER   
            The largest exponent before overflow   

    RMAX    (output) DOUBLE PRECISION   
            The largest machine floating-point number.   

   ===================================================================== 
  


       First compute LEXP and UEXP, two powers of 2 that bound   
       abs(EMIN). We then assume that EMAX + abs(EMIN) will sum   
       approximately to the bound that is closest to abs(EMIN).   
       (EMAX is the exponent of the required number RMAX). */
    /* Table of constant values */
    static double c_b5 = 0.;
    
    /* System generated locals */
    integer i__1;
    double d__1;
    /* Local variables */
    static integer lexp;
    static double oldy;
    static integer uexp, i;
    static double y, z;
    static integer nbits;
    extern double dlamc3_(double *, double *);
    static double recbas;
    static integer exbits, expsum, try__;



    lexp = 1;
    exbits = 1;
L10:
    try__ = lexp << 1;
    if (try__ <= -(*emin)) {
	lexp = try__;
	++exbits;
	goto L10;
    }
    if (lexp == -(*emin)) {
	uexp = lexp;
    } else {
	uexp = try__;
	++exbits;
    }

/*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater   
       than or equal to EMIN. EXBITS is the number of bits needed to   
       store the exponent. */

    if (uexp + *emin > -lexp - *emin) {
	expsum = lexp << 1;
    } else {
	expsum = uexp << 1;
    }

/*     EXPSUM is the exponent range, approximately equal to   
       EMAX - EMIN + 1 . */

    *emax = expsum + *emin - 1;
    nbits = exbits + 1 + *p;

/*     NBITS is the total number of bits needed to store a   
       floating-point number. */

    if (nbits % 2 == 1 && *beta == 2) {

/*        Either there are an odd number of bits used to store a   
          floating-point number, which is unlikely, or some bits are 
  
          not used in the representation of numbers, which is possible
,   
          (e.g. Cray machines) or the mantissa has an implicit bit,   
          (e.g. IEEE machines, Dec Vax machines), which is perhaps the
   
          most likely. We have to assume the last alternative.   
          If this is true, then we need to reduce EMAX by one because 
  
          there must be some way of representing zero in an implicit-b
it   
          system. On machines like Cray, we are reducing EMAX by one 
  
          unnecessarily. */

	--(*emax);
    }

    if (*ieee) {

/*        Assume we are on an IEEE machine which reserves one exponent
   
          for infinity and NaN. */

	--(*emax);
    }

/*     Now create RMAX, the largest machine number, which should   
       be equal to (1.0 - BETA**(-P)) * BETA**EMAX .   

       First compute 1.0 - BETA**(-P), being careful that the   
       result is less than 1.0 . */

    recbas = 1. / *beta;
    z = *beta - 1.;
    y = 0.;
    i__1 = *p;
    for (i = 1; i <= *p; ++i) {
	z *= recbas;
	if (y < 1.) {
	    oldy = y;
	}
	y = dlamc3_(&y, &z);
/* L20: */
    }
    if (y >= 1.) {
	y = oldy;
    }

/*     Now multiply by BETA**EMAX to get RMAX. */

    i__1 = *emax;
    for (i = 1; i <= *emax; ++i) {
	d__1 = y * *beta;
	y = dlamc3_(&d__1, &c_b5);
/* L30: */
    }

    *rmax = y;
    return 0;

/*     End of DLAMC5 */

} /* dlamc5_ */


static double dlanst_(char *norm, integer *n, double *d, double *e)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLANST  returns the value of the one norm,  or the Frobenius norm, or 
  
    the  infinity norm,  or the  element of  largest absolute value  of a 
  
    float symmetric tridiagonal matrix A.   

    Description   
    ===========   

    DLANST returns the value   

       DLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
                (   
                ( norm1(A),         NORM = '1', 'O' or 'o'   
                (   
                ( normI(A),         NORM = 'I' or 'i'   
                (   
                ( normF(A),         NORM = 'F', 'f', 'E' or 'e'   

    where  norm1  denotes the  one norm of a matrix (maximum column sum), 
  
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and 
  
    normF  denotes the  Frobenius norm of a matrix (square root of sum of 
  
    squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies the value to be returned in DLANST as described   
            above.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.  When N = 0, DLANST is   
            set to zero.   

    D       (input) DOUBLE PRECISION array, dimension (N)   
            The diagonal elements of A.   

    E       (input) DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) sub-diagonal or super-diagonal elements of A.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer i__1;
    double ret_val, d__1, d__2, d__3, d__4, d__5;
    /* Local variables */
    static integer i;
    static double scale;
    extern logical lsame_(char *, char *);
    static double anorm;
    extern /* Subroutine */ int dlassq_(integer *, double *, integer *, 
	    double *, double *);
    static double sum;



#define E(I) e[(I)-1]
#define D(I) d[(I)-1]


    if (*n <= 0) {
	anorm = 0.;
    } else if (lsame_(norm, "M")) {

/*        Find max(abs(A(i,j))). */

	anorm = (d__1 = D(*n), abs(d__1));
	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {
/* Computing MAX */
	    d__2 = anorm, d__3 = (d__1 = D(i), abs(d__1));
	    anorm = max(d__2,d__3);
/* Computing MAX */
	    d__2 = anorm, d__3 = (d__1 = E(i), abs(d__1));
	    anorm = max(d__2,d__3);
/* L10: */
	}
    } else if (lsame_(norm, "O") || *(unsigned char *)norm == '1' || 
	    lsame_(norm, "I")) {

/*        Find norm1(A). */

	if (*n == 1) {
	    anorm = abs(D(1));
	} else {
/* Computing MAX */
	    d__3 = abs(D(1)) + abs(E(1)), d__4 = (d__1 = E(*n - 1), abs(d__1))
		     + (d__2 = D(*n), abs(d__2));
	    anorm = max(d__3,d__4);
	    i__1 = *n - 1;
	    for (i = 2; i <= *n-1; ++i) {
/* Computing MAX */
		d__4 = anorm, d__5 = (d__1 = D(i), abs(d__1)) + (d__2 = E(i), 
			abs(d__2)) + (d__3 = E(i - 1), abs(d__3));
		anorm = max(d__4,d__5);
/* L20: */
	    }
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.;
	sum = 1.;
	if (*n > 1) {
	    i__1 = *n - 1;
	    dlassq_(&i__1, &E(1), &c__1, &scale, &sum);
	    sum *= 2;
	}
	dlassq_(n, &D(1), &c__1, &scale, &sum);
	anorm = scale * sqrt(sum);
    }

    ret_val = anorm;
    return ret_val;

/*     End of DLANST */

} /* dlanst_ */


static double dlansy_(char *norm, char *uplo, integer *n, double *a, integer 
	*lda, double *work)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLANSY  returns the value of the one norm,  or the Frobenius norm, or 
  
    the  infinity norm,  or the  element of  largest absolute value  of a 
  
    float symmetric matrix A.   

    Description   
    ===========   

    DLANSY returns the value   

       DLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
                (   
                ( norm1(A),         NORM = '1', 'O' or 'o'   
                (   
                ( normI(A),         NORM = 'I' or 'i'   
                (   
                ( normF(A),         NORM = 'F', 'f', 'E' or 'e'   

    where  norm1  denotes the  one norm of a matrix (maximum column sum), 
  
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and 
  
    normF  denotes the  Frobenius norm of a matrix (square root of sum of 
  
    squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies the value to be returned in DLANSY as described   
            above.   

    UPLO    (input) CHARACTER*1   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is to be referenced.   
            = 'U':  Upper triangular part of A is referenced   
            = 'L':  Lower triangular part of A is referenced   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.  When N = 0, DLANSY is   
            set to zero.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
            The symmetric matrix A.  If UPLO = 'U', the leading n by n   
            upper triangular part of A contains the upper triangular part 
  
            of the matrix A, and the strictly lower triangular part of A 
  
            is not referenced.  If UPLO = 'L', the leading n by n lower   
            triangular part of A contains the lower triangular part of   
            the matrix A, and the strictly upper triangular part of A is 
  
            not referenced.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(N,1).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),   
            where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,   
            WORK is not referenced.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    double ret_val, d__1, d__2, d__3;
    /* Local variables */
    static double absa;
    static integer i, j;
    static double scale;
    extern logical lsame_(char *, char *);
    static double value;
    extern /* Subroutine */ int dlassq_(integer *, double *, integer *, 
	    double *, double *);
    static double sum;



#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    if (*n == 0) {
	value = 0.;
    } else if (lsame_(norm, "M")) {

/*        Find max(abs(A(i,j))). */

	value = 0.;
	if (lsame_(uplo, "U")) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j;
		for (i = 1; i <= j; ++i) {
/* Computing MAX */
		    d__2 = value, d__3 = (d__1 = A(i,j), abs(d__1))
			    ;
		    value = max(d__2,d__3);
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = *n;
		for (i = j; i <= *n; ++i) {
/* Computing MAX */
		    d__2 = value, d__3 = (d__1 = A(i,j), abs(d__1))
			    ;
		    value = max(d__2,d__3);
/* L30: */
		}
/* L40: */
	    }
	}
    } else if (lsame_(norm, "I") || lsame_(norm, "O") || *(
	    unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

	value = 0.;
	if (lsame_(uplo, "U")) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		sum = 0.;
		i__2 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    absa = (d__1 = A(i,j), abs(d__1));
		    sum += absa;
		    WORK(i) += absa;
/* L50: */
		}
		WORK(j) = sum + (d__1 = A(j,j), abs(d__1));
/* L60: */
	    }
	    i__1 = *n;
	    for (i = 1; i <= *n; ++i) {
/* Computing MAX */
		d__1 = value, d__2 = WORK(i);
		value = max(d__1,d__2);
/* L70: */
	    }
	} else {
	    i__1 = *n;
	    for (i = 1; i <= *n; ++i) {
		WORK(i) = 0.;
/* L80: */
	    }
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		sum = WORK(j) + (d__1 = A(j,j), abs(d__1));
		i__2 = *n;
		for (i = j + 1; i <= *n; ++i) {
		    absa = (d__1 = A(i,j), abs(d__1));
		    sum += absa;
		    WORK(i) += absa;
/* L90: */
		}
		value = max(value,sum);
/* L100: */
	    }
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.;
	sum = 1.;
	if (lsame_(uplo, "U")) {
	    i__1 = *n;
	    for (j = 2; j <= *n; ++j) {
		i__2 = j - 1;
		dlassq_(&i__2, &A(1,j), &c__1, &scale, &sum);
/* L110: */
	    }
	} else {
	    i__1 = *n - 1;
	    for (j = 1; j <= *n-1; ++j) {
		i__2 = *n - j;
		dlassq_(&i__2, &A(j+1,j), &c__1, &scale, &sum);
/* L120: */
	    }
	}
	sum *= 2;
	i__1 = *lda + 1;
	dlassq_(n, &A(1,1), &i__1, &scale, &sum);
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of DLANSY */

} /* dlansy_ */


static double dlapy2_(double *x, double *y)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary 
  
    overflow.   

    Arguments   
    =========   

    X       (input) DOUBLE PRECISION   
    Y       (input) DOUBLE PRECISION   
            X and Y specify the values x and y.   

    ===================================================================== 
*/
/* >>Start of File<<   
       System generated locals */
    double ret_val, d__1;
    /* Local variables */
    static double xabs, yabs, w, z;



    xabs = abs(*x);
    yabs = abs(*y);
    w = max(xabs,yabs);
    z = min(xabs,yabs);
    if (z == 0.) {
	ret_val = w;
    } else {
/* Computing 2nd power */
	d__1 = z / w;
	ret_val = w * sqrt(d__1 * d__1 + 1.);
    }
    return ret_val;

/*     End of DLAPY2 */

} /* dlapy2_ */


/* Subroutine */
static  int dlarf_(char *side, integer *m, integer *n, double *v,
	 integer *incv, double *tau, double *c, integer *ldc, 
	double *work)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLARF applies a float elementary reflector H to a float m by n matrix   
    C, from either the left or the right. H is represented in the form   

          H = I - tau * v * v'   

    where tau is a float scalar and v is a float vector.   

    If tau = 0, then H is taken to be the unit matrix.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': form  H * C   
            = 'R': form  C * H   

    M       (input) INTEGER   
            The number of rows of the matrix C.   

    N       (input) INTEGER   
            The number of columns of the matrix C.   

    V       (input) DOUBLE PRECISION array, dimension   
                       (1 + (M-1)*abs(INCV)) if SIDE = 'L'   
                    or (1 + (N-1)*abs(INCV)) if SIDE = 'R'   
            The vector v in the representation of H. V is not used if   
            TAU = 0.   

    INCV    (input) INTEGER   
            The increment between elements of v. INCV <> 0.   

    TAU     (input) DOUBLE PRECISION   
            The value tau in the representation of H.   

    C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the m by n matrix C.   
            On exit, C is overwritten by the matrix H * C if SIDE = 'L', 
  
            or C * H if SIDE = 'R'.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= max(1,M).   

    WORK    (workspace) DOUBLE PRECISION array, dimension   
                           (N) if SIDE = 'L'   
                        or (M) if SIDE = 'R'   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static double c_b4 = 1.;
    static double c_b5 = 0.;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer c_dim1, c_offset;
    double d__1;
    /* Local variables */
    extern /* Subroutine */ int dger_(integer *, integer *, double *, 
	    double *, integer *, double *, integer *, double *, 
	    integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    double *, double *, integer *, double *, integer *, 
	    double *, double *, integer *);



#define V(I) v[(I)-1]
#define WORK(I) work[(I)-1]

#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    if (lsame_(side, "L")) {

/*        Form  H * C */

	if (*tau != 0.) {

/*           w := C' * v */

	    dgemv_("Transpose", m, n, &c_b4, &C(1,1), ldc, &V(1), incv, &
		    c_b5, &WORK(1), &c__1);

/*           C := C - v * w' */

	    d__1 = -(*tau);
	    dger_(m, n, &d__1, &V(1), incv, &WORK(1), &c__1, &C(1,1), 
		    ldc);
	}
    } else {

/*        Form  C * H */

	if (*tau != 0.) {

/*           w := C * v */

	    dgemv_("No transpose", m, n, &c_b4, &C(1,1), ldc, &V(1), 
		    incv, &c_b5, &WORK(1), &c__1);

/*           C := C - w * v' */

	    d__1 = -(*tau);
	    dger_(m, n, &d__1, &WORK(1), &c__1, &V(1), incv, &C(1,1), 
		    ldc);
	}
    }
    return 0;

/*     End of DLARF */

} /* dlarf_ */


/* Subroutine */
static  int dlarfb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, double *v, integer *
	ldv, double *t, integer *ldt, double *c, integer *ldc, 
	double *work, integer *ldwork)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLARFB applies a float block reflector H or its transpose H' to a   
    float m by n matrix C, from either the left or the right.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': apply H or H' from the Left   
            = 'R': apply H or H' from the Right   

    TRANS   (input) CHARACTER*1   
            = 'N': apply H (No transpose)   
            = 'T': apply H' (Transpose)   

    DIRECT  (input) CHARACTER*1   
            Indicates how H is formed from a product of elementary   
            reflectors   
            = 'F': H = H(1) H(2) . . . H(k) (Forward)   
            = 'B': H = H(k) . . . H(2) H(1) (Backward)   

    STOREV  (input) CHARACTER*1   
            Indicates how the vectors which define the elementary   
            reflectors are stored:   
            = 'C': Columnwise   
            = 'R': Rowwise   

    M       (input) INTEGER   
            The number of rows of the matrix C.   

    N       (input) INTEGER   
            The number of columns of the matrix C.   

    K       (input) INTEGER   
            The order of the matrix T (= the number of elementary   
            reflectors whose product defines the block reflector).   

    V       (input) DOUBLE PRECISION array, dimension   
                                  (LDV,K) if STOREV = 'C'   
                                  (LDV,M) if STOREV = 'R' and SIDE = 'L' 
  
                                  (LDV,N) if STOREV = 'R' and SIDE = 'R' 
  
            The matrix V. See further details.   

    LDV     (input) INTEGER   
            The leading dimension of the array V.   
            If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);   
            if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);   
            if STOREV = 'R', LDV >= K.   

    T       (input) DOUBLE PRECISION array, dimension (LDT,K)   
            The triangular k by k matrix T in the representation of the   
            block reflector.   

    LDT     (input) INTEGER   
            The leading dimension of the array T. LDT >= K.   

    C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the m by n matrix C.   
            On exit, C is overwritten by H*C or H'*C or C*H or C*H'.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDA >= max(1,M).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)   

    LDWORK  (input) INTEGER   
            The leading dimension of the array WORK.   
            If SIDE = 'L', LDWORK >= max(1,N);   
            if SIDE = 'R', LDWORK >= max(1,M).   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static double c_b14 = 1.;
    static double c_b25 = -1.;
    
    /* System generated locals */
    integer c_dim1, c_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1, 
	    work_offset, i__1, i__2;
    /* Local variables */
    static integer i, j;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, double *, double *, integer *, double *, 
	    integer *, double *, double *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dcopy_(integer *, double *, integer *, 
	    double *, integer *);
    extern int dtrmm_(char *side, char *uplo, char *transa, char *diag,
		      integer *m, integer *n,
		      double *alpha, double *a, integer *lda,
		      double *b, integer *ldb);
    static char transt[1];



#undef V
#define V(I,J) v[(I)-1 + ((J)-1)* ( *ldv)]
#define T(I,J) t[(I)-1 + ((J)-1)* ( *ldt)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]
#undef WORK
#define WORK(I,J) work[(I)-1 + ((J)-1)* ( *ldwork)]

    if (*m <= 0 || *n <= 0) {
	return 0;
    }

    if (lsame_(trans, "N")) {
	*(unsigned char *)transt = 'T';
    } else {
	*(unsigned char *)transt = 'N';
    }

    if (lsame_(storev, "C")) {

	if (lsame_(direct, "F")) {

/*           Let  V =  ( V1 )    (first K rows)   
                       ( V2 )   
             where  V1  is unit lower triangular. */

	    if (lsame_(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in 
WORK)   

                W := C1' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    dcopy_(n, &C(j,1), ldc, &WORK(1,j), &
			    c__1);
/* L10: */
		}

/*              W := W * V1 */

		dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b14,
			 &V(1,1), ldv, &WORK(1,1), ldwork);
		if (*m > *k) {

/*                 W := W + C2'*V2 */

		    i__1 = *m - *k;
		    dgemm_("Transpose", "No transpose", n, k, &i__1, &c_b14, &
			    C(*k+1,1), ldc, &V(*k+1,1), ldv,
			     &c_b14, &WORK(1,1), ldwork);
		}

/*              W := W * T'  or  W * T */

		dtrmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);

/*              C := C - V * W' */

		if (*m > *k) {

/*                 C2 := C2 - V2 * W' */

		    i__1 = *m - *k;
		    dgemm_("No transpose", "Transpose", &i__1, n, k, &c_b25, &
			    V(*k+1,1), ldv, &WORK(1,1), 
			    ldwork, &c_b14, &C(*k+1,1), ldc)
			    ;
		}

/*              W := W * V1' */

		dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
			V(1,1), ldv, &WORK(1,1), ldwork);

/*              C1 := C1 - W' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    i__2 = *n;
		    for (i = 1; i <= *n; ++i) {
			C(j,i) -= WORK(i,j);
/* L20: */
		    }
/* L30: */
		}

	    } else if (lsame_(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
   

                W := C * V  =  (C1*V1 + C2*V2)  (stored in WOR
K)   

                W := C1 */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    dcopy_(m, &C(1,j), &c__1, &WORK(1,j), &c__1);
/* L40: */
		}

/*              W := W * V1 */

		dtrmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b14,
			 &V(1,1), ldv, &WORK(1,1), ldwork);
		if (*n > *k) {

/*                 W := W + C2 * V2 */

		    i__1 = *n - *k;
		    dgemm_("No transpose", "No transpose", m, k, &i__1, &
			    c_b14, &C(1,*k+1), ldc, &V(*k+1,1), ldv, &c_b14, &WORK(1,1), 
			    ldwork);
		}

/*              W := W * T  or  W * T' */

		dtrmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);

/*              C := C - W * V' */

		if (*n > *k) {

/*                 C2 := C2 - W * V2' */

		    i__1 = *n - *k;
		    dgemm_("No transpose", "Transpose", m, &i__1, k, &c_b25, &
			    WORK(1,1), ldwork, &V(*k+1,1), 
			    ldv, &c_b14, &C(1,*k+1), ldc);
		}

/*              W := W * V1' */

		dtrmm_("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
			V(1,1), ldv, &WORK(1,1), ldwork);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			C(i,j) -= WORK(i,j);
/* L50: */
		    }
/* L60: */
		}
	    }

	} else {

/*           Let  V =  ( V1 )   
                       ( V2 )    (last K rows)   
             where  V2  is unit upper triangular. */

	    if (lsame_(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in 
WORK)   

                W := C2' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    dcopy_(n, &C(*m-*k+j,1), ldc, &WORK(1,j), &c__1);
/* L70: */
		}

/*              W := W * V2 */

		dtrmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b14,
			 &V(*m-*k+1,1), ldv, &WORK(1,1), 
			ldwork);
		if (*m > *k) {

/*                 W := W + C1'*V1 */

		    i__1 = *m - *k;
		    dgemm_("Transpose", "No transpose", n, k, &i__1, &c_b14, &
			    C(1,1), ldc, &V(1,1), ldv, &c_b14, &
			    WORK(1,1), ldwork);
		}

/*              W := W * T'  or  W * T */

		dtrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);

/*              C := C - V * W' */

		if (*m > *k) {

/*                 C1 := C1 - V1 * W' */

		    i__1 = *m - *k;
		    dgemm_("No transpose", "Transpose", &i__1, n, k, &c_b25, &
			    V(1,1), ldv, &WORK(1,1), ldwork, &
			    c_b14, &C(1,1), ldc);
		}

/*              W := W * V2' */

		dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
			V(*m-*k+1,1), ldv, &WORK(1,1), 
			ldwork);

/*              C2 := C2 - W' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    i__2 = *n;
		    for (i = 1; i <= *n; ++i) {
			C(*m-*k+j,i) -= WORK(i,j)
				;
/* L80: */
		    }
/* L90: */
		}

	    } else if (lsame_(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
   

                W := C * V  =  (C1*V1 + C2*V2)  (stored in WOR
K)   

                W := C2 */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    dcopy_(m, &C(1,*n-*k+j), &c__1, &WORK(1,j), &c__1);
/* L100: */
		}

/*              W := W * V2 */

		dtrmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b14,
			 &V(*n-*k+1,1), ldv, &WORK(1,1), 
			ldwork);
		if (*n > *k) {

/*                 W := W + C1 * V1 */

		    i__1 = *n - *k;
		    dgemm_("No transpose", "No transpose", m, k, &i__1, &
			    c_b14, &C(1,1), ldc, &V(1,1), ldv, &
			    c_b14, &WORK(1,1), ldwork);
		}

/*              W := W * T  or  W * T' */

		dtrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);

/*              C := C - W * V' */

		if (*n > *k) {

/*                 C1 := C1 - W * V1' */

		    i__1 = *n - *k;
		    dgemm_("No transpose", "Transpose", m, &i__1, k, &c_b25, &
			    WORK(1,1), ldwork, &V(1,1), ldv, &
			    c_b14, &C(1,1), ldc);
		}

/*              W := W * V2' */

		dtrmm_("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
			V(*n-*k+1,1), ldv, &WORK(1,1), 
			ldwork);

/*              C2 := C2 - W */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			C(i,*n-*k+j) -= WORK(i,j);
/* L110: */
		    }
/* L120: */
		}
	    }
	}

    } else if (lsame_(storev, "R")) {

	if (lsame_(direct, "F")) {

/*           Let  V =  ( V1  V2 )    (V1: first K columns)   
             where  V1  is unit upper triangular. */

	    if (lsame_(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V'  =  (C1'*V1' + C2'*V2') (stored i
n WORK)   

                W := C1' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    dcopy_(n, &C(j,1), ldc, &WORK(1,j), &
			    c__1);
/* L130: */
		}

/*              W := W * V1' */

		dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
			V(1,1), ldv, &WORK(1,1), ldwork);
		if (*m > *k) {

/*                 W := W + C2'*V2' */

		    i__1 = *m - *k;
		    dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b14, &C(*k+1,1), ldc, &V(1,*k+1), 
			    ldv, &c_b14, &WORK(1,1), ldwork);
		}

/*              W := W * T'  or  W * T */

		dtrmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);

/*              C := C - V' * W' */

		if (*m > *k) {

/*                 C2 := C2 - V2' * W' */

		    i__1 = *m - *k;
		    dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b25, &V(1,*k+1), ldv, &WORK(1,1), 
			    ldwork, &c_b14, &C(*k+1,1), ldc);
		}

/*              W := W * V1 */

		dtrmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b14,
			 &V(1,1), ldv, &WORK(1,1), ldwork);

/*              C1 := C1 - W' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    i__2 = *n;
		    for (i = 1; i <= *n; ++i) {
			C(j,i) -= WORK(i,j);
/* L140: */
		    }
/* L150: */
		}

	    } else if (lsame_(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
   

                W := C * V'  =  (C1*V1' + C2*V2')  (stored in 
WORK)   

                W := C1 */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    dcopy_(m, &C(1,j), &c__1, &WORK(1,j), &c__1);
/* L160: */
		}

/*              W := W * V1' */

		dtrmm_("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
			V(1,1), ldv, &WORK(1,1), ldwork);
		if (*n > *k) {

/*                 W := W + C2 * V2' */

		    i__1 = *n - *k;
		    dgemm_("No transpose", "Transpose", m, k, &i__1, &c_b14, &
			    C(1,*k+1), ldc, &V(1,*k+1), ldv, &c_b14, &WORK(1,1), 
			    ldwork);
		}

/*              W := W * T  or  W * T' */

		dtrmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);

/*              C := C - W * V */

		if (*n > *k) {

/*                 C2 := C2 - W * V2 */

		    i__1 = *n - *k;
		    dgemm_("No transpose", "No transpose", m, &i__1, k, &
			    c_b25, &WORK(1,1), ldwork, &V(1,*k+1), ldv, &c_b14, &C(1,*k+1), ldc);
		}

/*              W := W * V1 */

		dtrmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b14,
			 &V(1,1), ldv, &WORK(1,1), ldwork);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			C(i,j) -= WORK(i,j);
/* L170: */
		    }
/* L180: */
		}

	    }

	} else {

/*           Let  V =  ( V1  V2 )    (V2: last K columns)   
             where  V2  is unit lower triangular. */

	    if (lsame_(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V'  =  (C1'*V1' + C2'*V2') (stored i
n WORK)   

                W := C2' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    dcopy_(n, &C(*m-*k+j,1), ldc, &WORK(1,j), &c__1);
/* L190: */
		}

/*              W := W * V2' */

		dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
			V(1,*m-*k+1), ldv, &WORK(1,1)
			, ldwork);
		if (*m > *k) {

/*                 W := W + C1'*V1' */

		    i__1 = *m - *k;
		    dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b14, &C(1,1), ldc, &V(1,1), ldv, &c_b14, &WORK(1,1), ldwork);
		}

/*              W := W * T'  or  W * T */

		dtrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);

/*              C := C - V' * W' */

		if (*m > *k) {

/*                 C1 := C1 - V1' * W' */

		    i__1 = *m - *k;
		    dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b25, &V(1,1), ldv, &WORK(1,1), ldwork, &
			    c_b14, &C(1,1), ldc);
		}

/*              W := W * V2 */

		dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b14,
			 &V(1,*m-*k+1), ldv, &WORK(1,1), ldwork);

/*              C2 := C2 - W' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    i__2 = *n;
		    for (i = 1; i <= *n; ++i) {
			C(*m-*k+j,i) -= WORK(i,j)
				;
/* L200: */
		    }
/* L210: */
		}

	    } else if (lsame_(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
   

                W := C * V'  =  (C1*V1' + C2*V2')  (stored in 
WORK)   

                W := C2 */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    dcopy_(m, &C(1,*n-*k+j), &c__1, &WORK(1,j), &c__1);
/* L220: */
		}

/*              W := W * V2' */

		dtrmm_("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
			V(1,*n-*k+1), ldv, &WORK(1,1)
			, ldwork);
		if (*n > *k) {

/*                 W := W + C1 * V1' */

		    i__1 = *n - *k;
		    dgemm_("No transpose", "Transpose", m, k, &i__1, &c_b14, &
			    C(1,1), ldc, &V(1,1), ldv, &c_b14, &
			    WORK(1,1), ldwork);
		}

/*              W := W * T  or  W * T' */

		dtrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);

/*              C := C - W * V */

		if (*n > *k) {

/*                 C1 := C1 - W * V1 */

		    i__1 = *n - *k;
		    dgemm_("No transpose", "No transpose", m, &i__1, k, &
			    c_b25, &WORK(1,1), ldwork, &V(1,1), 
			    ldv, &c_b14, &C(1,1), ldc);
		}

/*              W := W * V2 */

		dtrmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b14,
			 &V(1,*n-*k+1), ldv, &WORK(1,1), ldwork);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			C(i,*n-*k+j) -= WORK(i,j);
/* L230: */
		    }
/* L240: */
		}

	    }

	}
    }

    return 0;

/*     End of DLARFB */

} /* dlarfb_ */


/* Subroutine */
static  int dlarfg_(integer *n, double *alpha, double *x, 
	integer *incx, double *tau)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLARFG generates a float elementary reflector H of order n, such   
    that   

          H * ( alpha ) = ( beta ),   H' * H = I.   
              (   x   )   (   0  )   

    where alpha and beta are scalars, and x is an (n-1)-element float   
    vector. H is represented in the form   

          H = I - tau * ( 1 ) * ( 1 v' ) ,   
                        ( v )   

    where tau is a float scalar and v is a float (n-1)-element   
    vector.   

    If the elements of x are all zero, then tau = 0 and H is taken to be 
  
    the unit matrix.   

    Otherwise  1 <= tau <= 2.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the elementary reflector.   

    ALPHA   (input/output) DOUBLE PRECISION   
            On entry, the value alpha.   
            On exit, it is overwritten with the value beta.   

    X       (input/output) DOUBLE PRECISION array, dimension   
                           (1+(N-2)*abs(INCX))   
            On entry, the vector x.   
            On exit, it is overwritten with the vector v.   

    INCX    (input) INTEGER   
            The increment between elements of X. INCX > 0.   

    TAU     (output) DOUBLE PRECISION   
            The value tau.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1;
    double d__1;
    /* Builtin functions */
    extern double d_sign(double *, double *);
    /* Local variables */
    static double beta;
    extern double dnrm2_(integer *, double *, integer *);
    static integer j;
    extern /* Subroutine */ int dscal_(integer *, double *, double *, 
	    integer *);
    static double xnorm;
    extern double dlapy2_(double *, double *), dlamch_(char *);
    static double safmin, rsafmn;
    static integer knt;


#define X(I) x[(I)-1]


    if (*n <= 1) {
	*tau = 0.;
	return 0;
    }

    i__1 = *n - 1;
    xnorm = dnrm2_(&i__1, &X(1), incx);

    if (xnorm == 0.) {

/*        H  =  I */

	*tau = 0.;
    } else {

/*        general case */

	d__1 = dlapy2_(alpha, &xnorm);
	beta = -d_sign(&d__1, alpha);
	safmin = dlamch_("S") / dlamch_("E");
	if (abs(beta) < safmin) {

/*           XNORM, BETA may be inaccurate; scale X and recompute 
them */

	    rsafmn = 1. / safmin;
	    knt = 0;
L10:
	    ++knt;
	    i__1 = *n - 1;
	    dscal_(&i__1, &rsafmn, &X(1), incx);
	    beta *= rsafmn;
	    *alpha *= rsafmn;
	    if (abs(beta) < safmin) {
		goto L10;
	    }

/*           New BETA is at most 1, at least SAFMIN */

	    i__1 = *n - 1;
	    xnorm = dnrm2_(&i__1, &X(1), incx);
	    d__1 = dlapy2_(alpha, &xnorm);
	    beta = -d_sign(&d__1, alpha);
	    *tau = (beta - *alpha) / beta;
	    i__1 = *n - 1;
	    d__1 = 1. / (*alpha - beta);
	    dscal_(&i__1, &d__1, &X(1), incx);

/*           If ALPHA is subnormal, it may lose relative accuracy 
*/

	    *alpha = beta;
	    i__1 = knt;
	    for (j = 1; j <= knt; ++j) {
		*alpha *= safmin;
/* L20: */
	    }
	} else {
	    *tau = (beta - *alpha) / beta;
	    i__1 = *n - 1;
	    d__1 = 1. / (*alpha - beta);
	    dscal_(&i__1, &d__1, &X(1), incx);
	    *alpha = beta;
	}
    }

    return 0;

/*     End of DLARFG */

} /* dlarfg_ */


/* Subroutine */
static  int dlarft_(char *direct, char *storev, integer *n, integer *
	k, double *v, integer *ldv, double *tau, double *t, 
	integer *ldt)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLARFT forms the triangular factor T of a float block reflector H   
    of order n, which is defined as a product of k elementary reflectors. 
  

    If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular; 
  

    If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular. 
  

    If STOREV = 'C', the vector which defines the elementary reflector   
    H(i) is stored in the i-th column of the array V, and   

       H  =  I - V * T * V'   

    If STOREV = 'R', the vector which defines the elementary reflector   
    H(i) is stored in the i-th row of the array V, and   

       H  =  I - V' * T * V   

    Arguments   
    =========   

    DIRECT  (input) CHARACTER*1   
            Specifies the order in which the elementary reflectors are   
            multiplied to form the block reflector:   
            = 'F': H = H(1) H(2) . . . H(k) (Forward)   
            = 'B': H = H(k) . . . H(2) H(1) (Backward)   

    STOREV  (input) CHARACTER*1   
            Specifies how the vectors which define the elementary   
            reflectors are stored (see also Further Details):   
            = 'C': columnwise   
            = 'R': rowwise   

    N       (input) INTEGER   
            The order of the block reflector H. N >= 0.   

    K       (input) INTEGER   
            The order of the triangular factor T (= the number of   
            elementary reflectors). K >= 1.   

    V       (input/output) DOUBLE PRECISION array, dimension   
                                 (LDV,K) if STOREV = 'C'   
                                 (LDV,N) if STOREV = 'R'   
            The matrix V. See further details.   

    LDV     (input) INTEGER   
            The leading dimension of the array V.   
            If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K. 
  

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i).   

    T       (output) DOUBLE PRECISION array, dimension (LDT,K)   
            The k by k triangular factor T of the block reflector.   
            If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is 
  
            lower triangular. The rest of the array is not used.   

    LDT     (input) INTEGER   
            The leading dimension of the array T. LDT >= K.   

    Further Details   
    ===============   

    The shape of the matrix V and the storage of the vectors which define 
  
    the H(i) is best illustrated by the following example with n = 5 and 
  
    k = 3. The elements equal to 1 are not stored; the corresponding   
    array elements are modified but restored on exit. The rest of the   
    array is not used.   

    DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R': 
  

                 V = (  1       )                 V = (  1 v1 v1 v1 v1 ) 
  
                     ( v1  1    )                     (     1 v2 v2 v2 ) 
  
                     ( v1 v2  1 )                     (        1 v3 v3 ) 
  
                     ( v1 v2 v3 )   
                     ( v1 v2 v3 )   

    DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R': 
  

                 V = ( v1 v2 v3 )                 V = ( v1 v1  1       ) 
  
                     ( v1 v2 v3 )                     ( v2 v2 v2  1    ) 
  
                     (  1 v2 v3 )                     ( v3 v3 v3 v3  1 ) 
  
                     (     1 v3 )   
                     (        1 )   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static double c_b8 = 0.;
    
    /* System generated locals */
    integer t_dim1, t_offset, v_dim1, v_offset, i__1, i__2, i__3;
    double d__1;
    /* Local variables */
    static integer i, j;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    double *, double *, integer *, double *, integer *, 
	    double *, double *, integer *), dtrmv_(char *, 
	    char *, char *, integer *, double *, integer *, double *, 
	    integer *);
    static double vii;



#define TAU(I) tau[(I)-1]

#define V(I,J) v[(I)-1 + ((J)-1)* ( *ldv)]
#define T(I,J) t[(I)-1 + ((J)-1)* ( *ldt)]

    if (*n == 0) {
	return 0;
    }

    if (lsame_(direct, "F")) {
	i__1 = *k;
	for (i = 1; i <= *k; ++i) {
	    if (TAU(i) == 0.) {

/*              H(i)  =  I */

		i__2 = i;
		for (j = 1; j <= i; ++j) {
		    T(j,i) = 0.;
/* L10: */
		}
	    } else {

/*              general case */

		vii = V(i,i);
		V(i,i) = 1.;
		if (lsame_(storev, "C")) {

/*                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' 
* V(i:n,i) */

		    i__2 = *n - i + 1;
		    i__3 = i - 1;
		    d__1 = -TAU(i);
		    dgemv_("Transpose", &i__2, &i__3, &d__1, &V(i,1), 
			    ldv, &V(i,i), &c__1, &c_b8, &T(1,i), &c__1);
		} else {

/*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) *
 V(i,i:n)' */

		    i__2 = i - 1;
		    i__3 = *n - i + 1;
		    d__1 = -TAU(i);
		    dgemv_("No transpose", &i__2, &i__3, &d__1, &V(1,i), ldv, &V(i,i), ldv, &c_b8, &T(1,i), &c__1);
		}
		V(i,i) = vii;

/*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i) */

		i__2 = i - 1;
		dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &T(1,1), ldt, &T(1,i), &c__1);
		T(i,i) = TAU(i);
	    }
/* L20: */
	}
    } else {
	for (i = *k; i >= 1; --i) {
	    if (TAU(i) == 0.) {

/*              H(i)  =  I */

		i__1 = *k;
		for (j = i; j <= *k; ++j) {
		    T(j,i) = 0.;
/* L30: */
		}
	    } else {

/*              general case */

		if (i < *k) {
		    if (lsame_(storev, "C")) {
			vii = V(*n-*k+i,i);
			V(*n-*k+i,i) = 1.;

/*                    T(i+1:k,i) :=   
                              - tau(i) * V(1:n-k+i,i+1
:k)' * V(1:n-k+i,i) */

			i__1 = *n - *k + i;
			i__2 = *k - i;
			d__1 = -TAU(i);
			dgemv_("Transpose", &i__1, &i__2, &d__1, &V(1,i+1), ldv, &V(1,i), &c__1, &
				c_b8, &T(i+1,i), &c__1);
			V(*n-*k+i,i) = vii;
		    } else {
			vii = V(i,*n-*k+i);
			V(i,*n-*k+i) = 1.;

/*                    T(i+1:k,i) :=   
                              - tau(i) * V(i+1:k,1:n-k
+i) * V(i,1:n-k+i)' */

			i__1 = *k - i;
			i__2 = *n - *k + i;
			d__1 = -TAU(i);
			dgemv_("No transpose", &i__1, &i__2, &d__1, &V(i+1,1), ldv, &V(i,1), ldv, &c_b8, &
				T(i+1,i), &c__1);
			V(i,*n-*k+i) = vii;
		    }

/*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,
i) */

		    i__1 = *k - i;
		    dtrmv_("Lower", "No transpose", "Non-unit", &i__1, &T(i+1,i+1), ldt, &T(i+1,i)
			    , &c__1);
		}
		T(i,i) = TAU(i);
	    }
/* L40: */
	}
    }
    return 0;

/*     End of DLARFT */

} /* dlarft_ */


/* Subroutine */
static  int dlartg_(double *f, double *g, double *cs, 
	double *sn, double *r)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLARTG generate a plane rotation so that   

       [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.   
       [ -SN  CS  ]     [ G ]     [ 0 ]   

    This is a slower, more accurate version of the BLAS1 routine DROTG,   
    with the following other differences:   
       F and G are unchanged on return.   
       If G=0, then CS=1 and SN=0.   
       If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any   
          floating point operations (saves work in DBDSQR when   
          there are zeros on the diagonal).   

    If F exceeds G in magnitude, CS will be positive.   

    Arguments   
    =========   

    F       (input) DOUBLE PRECISION   
            The first component of vector to be rotated.   

    G       (input) DOUBLE PRECISION   
            The second component of vector to be rotated.   

    CS      (output) DOUBLE PRECISION   
            The cosine of the rotation.   

    SN      (output) DOUBLE PRECISION   
            The sine of the rotation.   

    R       (output) DOUBLE PRECISION   
            The nonzero component of the rotated vector.   

    ===================================================================== 
*/
    /* Initialized data */
    static logical first = TRUE_;
    /* System generated locals */
    integer i__1;
    double d__1, d__2;
    /* Builtin functions */
    extern double pow_di(double *, integer *);
    /* Local variables */
    static integer i;
    static double scale;
    static integer count;
    static double f1, g1, safmn2, safmx2;
    extern double dlamch_(char *);
    static double safmin, eps;



    if (first) {
	first = FALSE_;
	safmin = dlamch_("S");
	eps = dlamch_("E");
	d__1 = dlamch_("B");
	i__1 = (integer) (log(safmin / eps) / log(dlamch_("B")) / 2.);
	safmn2 = pow_di(&d__1, &i__1);
	safmx2 = 1. / safmn2;
    }
    if (*g == 0.) {
	*cs = 1.;
	*sn = 0.;
	*r = *f;
    } else if (*f == 0.) {
	*cs = 0.;
	*sn = 1.;
	*r = *g;
    } else {
	f1 = *f;
	g1 = *g;
/* Computing MAX */
	d__1 = abs(f1), d__2 = abs(g1);
	scale = max(d__1,d__2);
	if (scale >= safmx2) {
	    count = 0;
L10:
	    ++count;
	    f1 *= safmn2;
	    g1 *= safmn2;
/* Computing MAX */
	    d__1 = abs(f1), d__2 = abs(g1);
	    scale = max(d__1,d__2);
	    if (scale >= safmx2) {
		goto L10;
	    }
/* Computing 2nd power */
	    d__1 = f1;
/* Computing 2nd power */
	    d__2 = g1;
	    *r = sqrt(d__1 * d__1 + d__2 * d__2);
	    *cs = f1 / *r;
	    *sn = g1 / *r;
	    i__1 = count;
	    for (i = 1; i <= count; ++i) {
		*r *= safmx2;
/* L20: */
	    }
	} else if (scale <= safmn2) {
	    count = 0;
L30:
	    ++count;
	    f1 *= safmx2;
	    g1 *= safmx2;
/* Computing MAX */
	    d__1 = abs(f1), d__2 = abs(g1);
	    scale = max(d__1,d__2);
	    if (scale <= safmn2) {
		goto L30;
	    }
/* Computing 2nd power */
	    d__1 = f1;
/* Computing 2nd power */
	    d__2 = g1;
	    *r = sqrt(d__1 * d__1 + d__2 * d__2);
	    *cs = f1 / *r;
	    *sn = g1 / *r;
	    i__1 = count;
	    for (i = 1; i <= count; ++i) {
		*r *= safmn2;
/* L40: */
	    }
	} else {
/* Computing 2nd power */
	    d__1 = f1;
/* Computing 2nd power */
	    d__2 = g1;
	    *r = sqrt(d__1 * d__1 + d__2 * d__2);
	    *cs = f1 / *r;
	    *sn = g1 / *r;
	}
	if (abs(*f) > abs(*g) && *cs < 0.) {
	    *cs = -(*cs);
	    *sn = -(*sn);
	    *r = -(*r);
	}
    }
    return 0;

/*     End of DLARTG */

} /* dlartg_ */


/* Subroutine */
static  int dlascl_(char *type, integer *kl, integer *ku, double 
	*cfrom, double *cto, integer *m, integer *n, double *a, 
	integer *lda, integer *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLASCL multiplies the M by N float matrix A by the float scalar   
    CTO/CFROM.  This is done without over/underflow as long as the final 
  
    result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that 
  
    A may be full, upper triangular, lower triangular, upper Hessenberg, 
  
    or banded.   

    Arguments   
    =========   

    TYPE    (input) CHARACTER*1   
            TYPE indices the storage type of the input matrix.   
            = 'G':  A is a full matrix.   
            = 'L':  A is a lower triangular matrix.   
            = 'U':  A is an upper triangular matrix.   
            = 'H':  A is an upper Hessenberg matrix.   
            = 'B':  A is a symmetric band matrix with lower bandwidth KL 
  
                    and upper bandwidth KU and with the only the lower   
                    half stored.   
            = 'Q':  A is a symmetric band matrix with lower bandwidth KL 
  
                    and upper bandwidth KU and with the only the upper   
                    half stored.   
            = 'Z':  A is a band matrix with lower bandwidth KL and upper 
  
                    bandwidth KU.   

    KL      (input) INTEGER   
            The lower bandwidth of A.  Referenced only if TYPE = 'B',   
            'Q' or 'Z'.   

    KU      (input) INTEGER   
            The upper bandwidth of A.  Referenced only if TYPE = 'B',   
            'Q' or 'Z'.   

    CFROM   (input) DOUBLE PRECISION   
    CTO     (input) DOUBLE PRECISION   
            The matrix A is multiplied by CTO/CFROM. A(I,J) is computed   
            without over/underflow if the final result CTO*A(I,J)/CFROM   
            can be represented without over/underflow.  CFROM must be   
            nonzero.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,M)   
            The matrix to be multiplied by CTO/CFROM.  See TYPE for the   
            storage type.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    INFO    (output) INTEGER   
            0  - successful exit   
            <0 - if INFO = -i, the i-th argument had an illegal value.   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    /* Local variables */
    static logical done;
    static double ctoc;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer itype, k1, k2, k3, k4;
    static double cfrom1;
    extern double dlamch_(char *);
    static double cfromc;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static double bignum, smlnum, mul, cto1;



#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;

    if (lsame_(type, "G")) {
	itype = 0;
    } else if (lsame_(type, "L")) {
	itype = 1;
    } else if (lsame_(type, "U")) {
	itype = 2;
    } else if (lsame_(type, "H")) {
	itype = 3;
    } else if (lsame_(type, "B")) {
	itype = 4;
    } else if (lsame_(type, "Q")) {
	itype = 5;
    } else if (lsame_(type, "Z")) {
	itype = 6;
    } else {
	itype = -1;
    }

    if (itype == -1) {
	*info = -1;
    } else if (*cfrom == 0.) {
	*info = -4;
    } else if (*m < 0) {
	*info = -6;
    } else if (*n < 0 || itype == 4 && *n != *m || itype == 5 && *n != *m) {
	*info = -7;
    } else if (itype <= 3 && *lda < max(1,*m)) {
	*info = -9;
    } else if (itype >= 4) {
/* Computing MAX */
	i__1 = *m - 1;
	if (*kl < 0 || *kl > max(i__1,0)) {
	    *info = -2;
	} else /* if(complicated condition) */ {
/* Computing MAX */
	    i__1 = *n - 1;
	    if (*ku < 0 || *ku > max(i__1,0) || (itype == 4 || itype == 5) && 
		    *kl != *ku) {
		*info = -3;
	    } else if (itype == 4 && *lda < *kl + 1 || itype == 5 && *lda < *
		    ku + 1 || itype == 6 && *lda < (*kl << 1) + *ku + 1) {
		*info = -9;
	    }
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLASCL", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *m == 0) {
	return 0;
    }

/*     Get machine parameters */

    smlnum = dlamch_("S");
    bignum = 1. / smlnum;

    cfromc = *cfrom;
    ctoc = *cto;

L10:
    cfrom1 = cfromc * smlnum;
    cto1 = ctoc / bignum;
    if (abs(cfrom1) > abs(ctoc) && ctoc != 0.) {
	mul = smlnum;
	done = FALSE_;
	cfromc = cfrom1;
    } else if (abs(cto1) > abs(cfromc)) {
	mul = bignum;
	done = FALSE_;
	ctoc = cto1;
    } else {
	mul = ctoc / cfromc;
	done = TRUE_;
    }

    if (itype == 0) {

/*        Full matrix */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = *m;
	    for (i = 1; i <= *m; ++i) {
		A(i,j) *= mul;
/* L20: */
	    }
/* L30: */
	}

    } else if (itype == 1) {

/*        Lower triangular matrix */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = *m;
	    for (i = j; i <= *m; ++i) {
		A(i,j) *= mul;
/* L40: */
	    }
/* L50: */
	}

    } else if (itype == 2) {

/*        Upper triangular matrix */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = min(j,*m);
	    for (i = 1; i <= min(j,*m); ++i) {
		A(i,j) *= mul;
/* L60: */
	    }
/* L70: */
	}

    } else if (itype == 3) {

/*        Upper Hessenberg matrix */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MIN */
	    i__3 = j + 1;
	    i__2 = min(i__3,*m);
	    for (i = 1; i <= min(j+1,*m); ++i) {
		A(i,j) *= mul;
/* L80: */
	    }
/* L90: */
	}

    } else if (itype == 4) {

/*        Lower half of a symmetric band matrix */

	k3 = *kl + 1;
	k4 = *n + 1;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MIN */
	    i__3 = k3, i__4 = k4 - j;
	    i__2 = min(i__3,i__4);
	    for (i = 1; i <= min(k3,k4-j); ++i) {
		A(i,j) *= mul;
/* L100: */
	    }
/* L110: */
	}

    } else if (itype == 5) {

/*        Upper half of a symmetric band matrix */

	k1 = *ku + 2;
	k3 = *ku + 1;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MAX */
	    i__2 = k1 - j;
	    i__3 = k3;
	    for (i = max(k1-j,1); i <= k3; ++i) {
		A(i,j) *= mul;
/* L120: */
	    }
/* L130: */
	}

    } else if (itype == 6) {

/*        Band matrix */

	k1 = *kl + *ku + 2;
	k2 = *kl + 1;
	k3 = (*kl << 1) + *ku + 1;
	k4 = *kl + *ku + 1 + *m;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MAX */
	    i__3 = k1 - j;
/* Computing MIN */
	    i__4 = k3, i__5 = k4 - j;
	    i__2 = min(i__4,i__5);
	    for (i = max(k1-j,k2); i <= min(k3,k4-j); ++i) {
		A(i,j) *= mul;
/* L140: */
	    }
/* L150: */
	}

    }

    if (! done) {
	goto L10;
    }

    return 0;

/*     End of DLASCL */

} /* dlascl_ */


/* Subroutine */
static  int dlaset_(char *uplo, integer *m, integer *n, double *
	alpha, double *beta, double *a, integer *lda)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLASET initializes an m-by-n matrix A to BETA on the diagonal and   
    ALPHA on the offdiagonals.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies the part of the matrix A to be set.   
            = 'U':      Upper triangular part is set; the strictly lower 
  
                        triangular part of A is not changed.   
            = 'L':      Lower triangular part is set; the strictly upper 
  
                        triangular part of A is not changed.   
            Otherwise:  All of the matrix A is set.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    ALPHA   (input) DOUBLE PRECISION   
            The constant to which the offdiagonal elements are to be set. 
  

    BETA    (input) DOUBLE PRECISION   
            The constant to which the diagonal elements are to be set.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On exit, the leading m-by-n submatrix of A is set as follows: 
  

            if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,   
            if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,   
            otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,   

            and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer i, j;
    extern logical lsame_(char *, char *);



#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    if (lsame_(uplo, "U")) {

/*        Set the strictly upper triangular or trapezoidal part of the
   
          array to ALPHA. */

	i__1 = *n;
	for (j = 2; j <= *n; ++j) {
/* Computing MIN */
	    i__3 = j - 1;
	    i__2 = min(i__3,*m);
	    for (i = 1; i <= min(j-1,*m); ++i) {
		A(i,j) = *alpha;
/* L10: */
	    }
/* L20: */
	}

    } else if (lsame_(uplo, "L")) {

/*        Set the strictly lower triangular or trapezoidal part of the
   
          array to ALPHA. */

	i__1 = min(*m,*n);
	for (j = 1; j <= min(*m,*n); ++j) {
	    i__2 = *m;
	    for (i = j + 1; i <= *m; ++i) {
		A(i,j) = *alpha;
/* L30: */
	    }
/* L40: */
	}

    } else {

/*        Set the leading m-by-n submatrix to ALPHA. */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = *m;
	    for (i = 1; i <= *m; ++i) {
		A(i,j) = *alpha;
/* L50: */
	    }
/* L60: */
	}
    }

/*     Set the first min(M,N) diagonal elements to BETA. */

    i__1 = min(*m,*n);
    for (i = 1; i <= min(*m,*n); ++i) {
	A(i,i) = *beta;
/* L70: */
    }

    return 0;

/*     End of DLASET */

} /* dlaset_ */


/* Subroutine */
static  int dlasr_(char *side, char *pivot, char *direct, integer *m,
	 integer *n, double *c, double *s, double *a, integer *
	lda)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLASR   performs the transformation   

       A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )   

       A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )   

    where A is an m by n float matrix and P is an orthogonal matrix,   
    consisting of a sequence of plane rotations determined by the   
    parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l' 
  
    and z = n when SIDE = 'R' or 'r' ):   

    When  DIRECT = 'F' or 'f'  ( Forward sequence ) then   

       P = P( z - 1 )*...*P( 2 )*P( 1 ),   

    and when DIRECT = 'B' or 'b'  ( Backward sequence ) then   

       P = P( 1 )*P( 2 )*...*P( z - 1 ),   

    where  P( k ) is a plane rotation matrix for the following planes:   

       when  PIVOT = 'V' or 'v'  ( Variable pivot ),   
          the plane ( k, k + 1 )   

       when  PIVOT = 'T' or 't'  ( Top pivot ),   
          the plane ( 1, k + 1 )   

       when  PIVOT = 'B' or 'b'  ( Bottom pivot ),   
          the plane ( k, z )   

    c( k ) and s( k )  must contain the  cosine and sine that define the 
  
    matrix  P( k ).  The two by two plane rotation part of the matrix   
    P( k ), R( k ), is assumed to be of the form   

       R( k ) = (  c( k )  s( k ) ).   
                ( -s( k )  c( k ) )   

    This version vectorises across rows of the array A when SIDE = 'L'.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            Specifies whether the plane rotation matrix P is applied to   
            A on the left or the right.   
            = 'L':  Left, compute A := P*A   
            = 'R':  Right, compute A:= A*P'   

    DIRECT  (input) CHARACTER*1   
            Specifies whether P is a forward or backward sequence of   
            plane rotations.   
            = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )   
            = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )   

    PIVOT   (input) CHARACTER*1   
            Specifies the plane for which P(k) is a plane rotation   
            matrix.   
            = 'V':  Variable pivot, the plane (k,k+1)   
            = 'T':  Top pivot, the plane (1,k+1)   
            = 'B':  Bottom pivot, the plane (k,z)   

    M       (input) INTEGER   
            The number of rows of the matrix A.  If m <= 1, an immediate 
  
            return is effected.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  If n <= 1, an   
            immediate return is effected.   

    C, S    (input) DOUBLE PRECISION arrays, dimension   
                    (M-1) if SIDE = 'L'   
                    (N-1) if SIDE = 'R'   
            c(k) and s(k) contain the cosine and sine that define the   
            matrix P(k).  The two by two plane rotation part of the   
            matrix P(k), R(k), is assumed to be of the form   
            R( k ) = (  c( k )  s( k ) ).   
                     ( -s( k )  c( k ) )   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            The m by n matrix A.  On exit, A is overwritten by P*A if   
            SIDE = 'R' or by A*P' if SIDE = 'L'.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    static integer info;
    static double temp;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static double ctemp, stemp;
    extern /* Subroutine */ int xerbla_(char *, integer *);


#undef C
#define C(I) c[(I)-1]
#define S(I) s[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! (lsame_(side, "L") || lsame_(side, "R"))) {
	info = 1;
    } else if (! (lsame_(pivot, "V") || lsame_(pivot, "T") || 
	    lsame_(pivot, "B"))) {
	info = 2;
    } else if (! (lsame_(direct, "F") || lsame_(direct, "B")))
	     {
	info = 3;
    } else if (*m < 0) {
	info = 4;
    } else if (*n < 0) {
	info = 5;
    } else if (*lda < max(1,*m)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("DLASR ", &info);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }
    if (lsame_(side, "L")) {

/*        Form  P * A */

	if (lsame_(pivot, "V")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m - 1;
		for (j = 1; j <= *m-1; ++j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *n;
			for (i = 1; i <= *n; ++i) {
			    temp = A(j+1,i);
			    A(j+1,i) = ctemp * temp - stemp * A(j,i);
			    A(j,i) = stemp * temp + ctemp * A(j,i);
/* L10: */
			}
		    }
/* L20: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *n;
			for (i = 1; i <= *n; ++i) {
			    temp = A(j+1,i);
			    A(j+1,i) = ctemp * temp - stemp * A(j,i);
			    A(j,i) = stemp * temp + ctemp * A(j,i);
/* L30: */
			}
		    }
/* L40: */
		}
	    }
	} else if (lsame_(pivot, "T")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m;
		for (j = 2; j <= *m; ++j) {
		    ctemp = C(j - 1);
		    stemp = S(j - 1);
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *n;
			for (i = 1; i <= *n; ++i) {
			    temp = A(j,i);
			    A(j,i) = ctemp * temp - stemp * A(1,i);
			    A(1,i) = stemp * temp + ctemp * A(1,i);
/* L50: */
			}
		    }
/* L60: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m; j >= 2; --j) {
		    ctemp = C(j - 1);
		    stemp = S(j - 1);
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *n;
			for (i = 1; i <= *n; ++i) {
			    temp = A(j,i);
			    A(j,i) = ctemp * temp - stemp * A(1,i);
			    A(1,i) = stemp * temp + ctemp * A(1,i);
/* L70: */
			}
		    }
/* L80: */
		}
	    }
	} else if (lsame_(pivot, "B")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m - 1;
		for (j = 1; j <= *m-1; ++j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *n;
			for (i = 1; i <= *n; ++i) {
			    temp = A(j,i);
			    A(j,i) = stemp * A(*m,i) + 
				    ctemp * temp;
			    A(*m,i) = ctemp * A(*m,i) - 
				    stemp * temp;
/* L90: */
			}
		    }
/* L100: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *n;
			for (i = 1; i <= *n; ++i) {
			    temp = A(j,i);
			    A(j,i) = stemp * A(*m,i) + 
				    ctemp * temp;
			    A(*m,i) = ctemp * A(*m,i) - 
				    stemp * temp;
/* L110: */
			}
		    }
/* L120: */
		}
	    }
	}
    } else if (lsame_(side, "R")) {

/*        Form A * P' */

	if (lsame_(pivot, "V")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n - 1;
		for (j = 1; j <= *n-1; ++j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *m;
			for (i = 1; i <= *m; ++i) {
			    temp = A(i,j+1);
			    A(i,j+1) = ctemp * temp - stemp * 
				    A(i,j);
			    A(i,j) = stemp * temp + ctemp * A(i,j);
/* L130: */
			}
		    }
/* L140: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *m;
			for (i = 1; i <= *m; ++i) {
			    temp = A(i,j+1);
			    A(i,j+1) = ctemp * temp - stemp * 
				    A(i,j);
			    A(i,j) = stemp * temp + ctemp * A(i,j);
/* L150: */
			}
		    }
/* L160: */
		}
	    }
	} else if (lsame_(pivot, "T")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n;
		for (j = 2; j <= *n; ++j) {
		    ctemp = C(j - 1);
		    stemp = S(j - 1);
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *m;
			for (i = 1; i <= *m; ++i) {
			    temp = A(i,j);
			    A(i,j) = ctemp * temp - stemp * A(i,1);
			    A(i,1) = stemp * temp + ctemp * A(i,1);
/* L170: */
			}
		    }
/* L180: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n; j >= 2; --j) {
		    ctemp = C(j - 1);
		    stemp = S(j - 1);
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *m;
			for (i = 1; i <= *m; ++i) {
			    temp = A(i,j);
			    A(i,j) = ctemp * temp - stemp * A(i,1);
			    A(i,1) = stemp * temp + ctemp * A(i,1);
/* L190: */
			}
		    }
/* L200: */
		}
	    }
	} else if (lsame_(pivot, "B")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n - 1;
		for (j = 1; j <= *n-1; ++j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *m;
			for (i = 1; i <= *m; ++i) {
			    temp = A(i,j);
			    A(i,j) = stemp * A(i,*n) + 
				    ctemp * temp;
			    A(i,*n) = ctemp * A(i,*n) - 
				    stemp * temp;
/* L210: */
			}
		    }
/* L220: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *m;
			for (i = 1; i <= *m; ++i) {
			    temp = A(i,j);
			    A(i,j) = stemp * A(i,*n) + 
				    ctemp * temp;
			    A(i,*n) = ctemp * A(i,*n) - 
				    stemp * temp;
/* L230: */
			}
		    }
/* L240: */
		}
	    }
	}
    }

    return 0;

/*     End of DLASR */

} /* dlasr_ */


/* Subroutine */
static  int dlasrt_(char *id, integer *n, double *d, integer *
	info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    Sort the numbers in D in increasing order (if ID = 'I') or   
    in decreasing order (if ID = 'D' ).   

    Use Quick Sort, reverting to Insertion sort on arrays of   
    size <= 20. Dimension of STACK limits N to about 2**32.   

    Arguments   
    =========   

    ID      (input) CHARACTER*1   
            = 'I': sort D in increasing order;   
            = 'D': sort D in decreasing order.   

    N       (input) INTEGER   
            The length of the array D.   

    D       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, the array to be sorted.   
            On exit, D has been sorted into increasing order   
            (D(1) <= ... <= D(N) ) or into decreasing order   
            (D(1) >= ... >= D(N) ), depending on ID.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input paramters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1, i__2;
    /* Local variables */
    static integer endd, i, j;
    extern logical lsame_(char *, char *);
    static integer stack[64]	/* was [2][32] */;
    static double dmnmx, d1, d2, d3;
    static integer start;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static integer stkpnt, dir;
    static double tmp;


#define STACK(I) stack[(I)]
#define WAS(I) was[(I)]
#define D(I) d[(I)-1]


    *info = 0;
    dir = -1;
    if (lsame_(id, "D")) {
	dir = 0;
    } else if (lsame_(id, "I")) {
	dir = 1;
    }
    if (dir == -1) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLASRT", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 1) {
	return 0;
    }

    stkpnt = 1;
    STACK(0) = 1;
    STACK(1) = *n;
L10:
    start = STACK((stkpnt << 1) - 2);
    endd = STACK((stkpnt << 1) - 1);
    --stkpnt;
    if (endd - start <= 20 && endd - start > 0) {

/*        Do Insertion sort on D( START:ENDD ) */

	if (dir == 0) {

/*           Sort into decreasing order */

	    i__1 = endd;
	    for (i = start + 1; i <= endd; ++i) {
		i__2 = start + 1;
		for (j = i; j >= start+1; --j) {
		    if (D(j) > D(j - 1)) {
			dmnmx = D(j);
			D(j) = D(j - 1);
			D(j - 1) = dmnmx;
		    } else {
			goto L30;
		    }
/* L20: */
		}
L30:
		;
	    }

	} else {

/*           Sort into increasing order */

	    i__1 = endd;
	    for (i = start + 1; i <= endd; ++i) {
		i__2 = start + 1;
		for (j = i; j >= start+1; --j) {
		    if (D(j) < D(j - 1)) {
			dmnmx = D(j);
			D(j) = D(j - 1);
			D(j - 1) = dmnmx;
		    } else {
			goto L50;
		    }
/* L40: */
		}
L50:
		;
	    }

	}

    } else if (endd - start > 20) {

/*        Partition D( START:ENDD ) and stack parts, largest one first
   

          Choose partition entry as median of 3 */

	d1 = D(start);
	d2 = D(endd);
	i = (start + endd) / 2;
	d3 = D(i);
	if (d1 < d2) {
	    if (d3 < d1) {
		dmnmx = d1;
	    } else if (d3 < d2) {
		dmnmx = d3;
	    } else {
		dmnmx = d2;
	    }
	} else {
	    if (d3 < d2) {
		dmnmx = d2;
	    } else if (d3 < d1) {
		dmnmx = d3;
	    } else {
		dmnmx = d1;
	    }
	}

	if (dir == 0) {

/*           Sort into decreasing order */

	    i = start - 1;
	    j = endd + 1;
L60:
L70:
	    --j;
	    if (D(j) < dmnmx) {
		goto L70;
	    }
L80:
	    ++i;
	    if (D(i) > dmnmx) {
		goto L80;
	    }
	    if (i < j) {
		tmp = D(i);
		D(i) = D(j);
		D(j) = tmp;
		goto L60;
	    }
	    if (j - start > endd - j - 1) {
		++stkpnt;
		STACK((stkpnt << 1) - 2) = start;
		STACK((stkpnt << 1) - 1) = j;
		++stkpnt;
		STACK((stkpnt << 1) - 2) = j + 1;
		STACK((stkpnt << 1) - 1) = endd;
	    } else {
		++stkpnt;
		STACK((stkpnt << 1) - 2) = j + 1;
		STACK((stkpnt << 1) - 1) = endd;
		++stkpnt;
		STACK((stkpnt << 1) - 2) = start;
		STACK((stkpnt << 1) - 1) = j;
	    }
	} else {

/*           Sort into increasing order */

	    i = start - 1;
	    j = endd + 1;
L90:
L100:
	    --j;
	    if (D(j) > dmnmx) {
		goto L100;
	    }
L110:
	    ++i;
	    if (D(i) < dmnmx) {
		goto L110;
	    }
	    if (i < j) {
		tmp = D(i);
		D(i) = D(j);
		D(j) = tmp;
		goto L90;
	    }
	    if (j - start > endd - j - 1) {
		++stkpnt;
		STACK((stkpnt << 1) - 2) = start;
		STACK((stkpnt << 1) - 1) = j;
		++stkpnt;
		STACK((stkpnt << 1) - 2) = j + 1;
		STACK((stkpnt << 1) - 1) = endd;
	    } else {
		++stkpnt;
		STACK((stkpnt << 1) - 2) = j + 1;
		STACK((stkpnt << 1) - 1) = endd;
		++stkpnt;
		STACK((stkpnt << 1) - 2) = start;
		STACK((stkpnt << 1) - 1) = j;
	    }
	}
    }
    if (stkpnt > 0) {
	goto L10;
    }
    return 0;

/*     End of DLASRT */

} /* dlasrt_ */


/* Subroutine */
static  int dlassq_(integer *n, double *x, integer *incx, 
	double *scale, double *sumsq)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLASSQ  returns the values  scl  and  smsq  such that   

       ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq, 
  

    where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is   
    assumed to be non-negative and  scl  returns the value   

       scl = max( scale, abs( x( i ) ) ).   

    scale and sumsq must be supplied in SCALE and SUMSQ and   
    scl and smsq are overwritten on SCALE and SUMSQ respectively.   

    The routine makes only one pass through the vector x.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of elements to be used from the vector X.   

    X       (input) DOUBLE PRECISION   
            The vector for which a scaled sum of squares is computed.   
               x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.   

    INCX    (input) INTEGER   
            The increment between successive values of the vector X.   
            INCX > 0.   

    SCALE   (input/output) DOUBLE PRECISION   
            On entry, the value  scale  in the equation above.   
            On exit, SCALE is overwritten with  scl , the scaling factor 
  
            for the sum of squares.   

    SUMSQ   (input/output) DOUBLE PRECISION   
            On entry, the value  sumsq  in the equation above.   
            On exit, SUMSQ is overwritten with  smsq , the basic sum of   
            squares from which  scl  has been factored out.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1, i__2;
    double d__1;
    /* Local variables */
    static double absxi;
    static integer ix;


#define X(I) x[(I)-1]


    if (*n > 0) {
	i__1 = (*n - 1) * *incx + 1;
	i__2 = *incx;
	for (ix = 1; *incx < 0 ? ix >= (*n-1)**incx+1 : ix <= (*n-1)**incx+1; ix += *incx) {
	    if (X(ix) != 0.) {
		absxi = (d__1 = X(ix), abs(d__1));
		if (*scale < absxi) {
/* Computing 2nd power */
		    d__1 = *scale / absxi;
		    *sumsq = *sumsq * (d__1 * d__1) + 1;
		    *scale = absxi;
		} else {
/* Computing 2nd power */
		    d__1 = absxi / *scale;
		    *sumsq += d__1 * d__1;
		}
	    }
/* L10: */
	}
    }
    return 0;

/*     End of DLASSQ */

} /* dlassq_ */


/* Subroutine */
static  int dlatrd_(char *uplo, integer *n, integer *nb, double *
	a, integer *lda, double *e, double *tau, double *w, 
	integer *ldw)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLATRD reduces NB rows and columns of a float symmetric matrix A to   
    symmetric tridiagonal form by an orthogonal similarity   
    transformation Q' * A * Q, and returns the matrices V and W which are 
  
    needed to apply the transformation to the unreduced part of A.   

    If UPLO = 'U', DLATRD reduces the last NB rows and columns of a   
    matrix, of which the upper triangle is supplied;   
    if UPLO = 'L', DLATRD reduces the first NB rows and columns of a   
    matrix, of which the lower triangle is supplied.   

    This is an auxiliary routine called by DSYTRD.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is stored:   
            = 'U': Upper triangular   
            = 'L': Lower triangular   

    N       (input) INTEGER   
            The order of the matrix A.   

    NB      (input) INTEGER   
            The number of rows and columns to be reduced.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading 
  
            n-by-n upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n-by-n lower triangular part of A contains the lower 
  
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   
            On exit:   
            if UPLO = 'U', the last NB columns have been reduced to   
              tridiagonal form, with the diagonal elements overwriting   
              the diagonal elements of A; the elements above the diagonal 
  
              with the array TAU, represent the orthogonal matrix Q as a 
  
              product of elementary reflectors;   
            if UPLO = 'L', the first NB columns have been reduced to   
              tridiagonal form, with the diagonal elements overwriting   
              the diagonal elements of A; the elements below the diagonal 
  
              with the array TAU, represent the  orthogonal matrix Q as a 
  
              product of elementary reflectors.   
            See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= (1,N).   

    E       (output) DOUBLE PRECISION array, dimension (N-1)   
            If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal   
            elements of the last NB columns of the reduced matrix;   
            if UPLO = 'L', E(1:nb) contains the subdiagonal elements of   
            the first NB columns of the reduced matrix.   

    TAU     (output) DOUBLE PRECISION array, dimension (N-1)   
            The scalar factors of the elementary reflectors, stored in   
            TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'. 
  
            See Further Details.   

    W       (output) DOUBLE PRECISION array, dimension (LDW,NB)   
            The n-by-nb matrix W required to update the unreduced part   
            of A.   

    LDW     (input) INTEGER   
            The leading dimension of the array W. LDW >= max(1,N).   

    Further Details   
    ===============   

    If UPLO = 'U', the matrix Q is represented as a product of elementary 
  
    reflectors   

       Q = H(n) H(n-1) . . . H(n-nb+1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a float scalar, and v is a float vector with   
    v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i), 
  
    and tau in TAU(i-1).   

    If UPLO = 'L', the matrix Q is represented as a product of elementary 
  
    reflectors   

       Q = H(1) H(2) . . . H(nb).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a float scalar, and v is a float vector with   
    v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i), 
  
    and tau in TAU(i).   

    The elements of the vectors v together form the n-by-nb matrix V   
    which is needed, with W, to apply the transformation to the unreduced 
  
    part of the matrix, using a symmetric rank-2k update of the form:   
    A := A - V*W' - W*V'.   

    The contents of A on exit are illustrated by the following examples   
    with n = 5 and nb = 2:   

    if UPLO = 'U':                       if UPLO = 'L':   

      (  a   a   a   v4  v5 )              (  d                  )   
      (      a   a   v4  v5 )              (  1   d              )   
      (          a   1   v5 )              (  v1  1   a          )   
      (              d   1  )              (  v1  v2  a   a      )   
      (                  d  )              (  v1  v2  a   a   a  )   

    where d denotes a diagonal element of the reduced matrix, a denotes   
    an element of the original matrix that is unchanged, and vi denotes   
    an element of the vector defining H(i).   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static double c_b5 = -1.;
    static double c_b6 = 1.;
    static integer c__1 = 1;
    static double c_b16 = 0.;
    
    /* System generated locals */
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3;
    /* Local variables */
    extern double ddot_(integer *, double *, integer *, double *, 
	    integer *);
    static integer i;
    static double alpha;
    extern /* Subroutine */ int dscal_(integer *, double *, double *, 
	    integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    double *, double *, integer *, double *, integer *, 
	    double *, double *, integer *), daxpy_(integer *, 
	    double *, double *, integer *, double *, integer *), 
	    dsymv_(char *, integer *, double *, double *, integer *, 
	    double *, integer *, double *, double *, integer *), dlarfg_(integer *, double *, double *, integer *,
	     double *);
    static integer iw;



#define E(I) e[(I)-1]
#define TAU(I) tau[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#undef W
#define W(I,J) w[(I)-1 + ((J)-1)* ( *ldw)]

    if (*n <= 0) {
	return 0;
    }

    if (lsame_(uplo, "U")) {

/*        Reduce last NB columns of upper triangle */

	i__1 = *n - *nb + 1;
	for (i = *n; i >= *n-*nb+1; --i) {
	    iw = i - *n + *nb;
	    if (i < *n) {

/*              Update A(1:i,i) */

		i__2 = *n - i;
		dgemv_("No transpose", &i, &i__2, &c_b5, &A(1,i+1), lda, &W(i,iw+1), ldw, &c_b6, &A(1,i), &c__1);
		i__2 = *n - i;
		dgemv_("No transpose", &i, &i__2, &c_b5, &W(1,iw+1), ldw, &A(i,i+1), lda, &c_b6, &A(1,i), &c__1);
	    }
	    if (i > 1) {

/*              Generate elementary reflector H(i) to annihila
te   
                A(1:i-2,i) */

		i__2 = i - 1;
		dlarfg_(&i__2, &A(i-1,i), &A(1,i), &
			c__1, &TAU(i - 1));
		E(i - 1) = A(i-1,i);
		A(i-1,i) = 1.;

/*              Compute W(1:i-1,i) */

		i__2 = i - 1;
		dsymv_("Upper", &i__2, &c_b6, &A(1,1), lda, &A(1,i), &c__1, &c_b16, &W(1,iw), &
			c__1);
		if (i < *n) {
		    i__2 = i - 1;
		    i__3 = *n - i;
		    dgemv_("Transpose", &i__2, &i__3, &c_b6, &W(1,iw+1), ldw, &A(1,i), &c__1, &
			    c_b16, &W(i+1,iw), &c__1);
		    i__2 = i - 1;
		    i__3 = *n - i;
		    dgemv_("No transpose", &i__2, &i__3, &c_b5, &A(1,i+1), lda, &W(i+1,iw), &c__1, 
			    &c_b6, &W(1,iw), &c__1);
		    i__2 = i - 1;
		    i__3 = *n - i;
		    dgemv_("Transpose", &i__2, &i__3, &c_b6, &A(1,i+1), lda, &A(1,i), &c__1, &
			    c_b16, &W(i+1,iw), &c__1);
		    i__2 = i - 1;
		    i__3 = *n - i;
		    dgemv_("No transpose", &i__2, &i__3, &c_b5, &W(1,iw+1), ldw, &W(i+1,iw), &c__1, 
			    &c_b6, &W(1,iw), &c__1);
		}
		i__2 = i - 1;
		dscal_(&i__2, &TAU(i - 1), &W(1,iw), &c__1);
		i__2 = i - 1;
		alpha = TAU(i - 1) * -.5 * ddot_(&i__2, &W(1,iw), &
			c__1, &A(1,i), &c__1);
		i__2 = i - 1;
		daxpy_(&i__2, &alpha, &A(1,i), &c__1, &W(1,iw), &c__1);
	    }

/* L10: */
	}
    } else {

/*        Reduce first NB columns of lower triangle */

	i__1 = *nb;
	for (i = 1; i <= *nb; ++i) {

/*           Update A(i:n,i) */

	    i__2 = *n - i + 1;
	    i__3 = i - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b5, &A(i,1), lda, &
		    W(i,1), ldw, &c_b6, &A(i,i), &c__1)
		    ;
	    i__2 = *n - i + 1;
	    i__3 = i - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b5, &W(i,1), ldw, &
		    A(i,1), lda, &c_b6, &A(i,i), &c__1)
		    ;
	    if (i < *n) {

/*              Generate elementary reflector H(i) to annihila
te   
                A(i+2:n,i) */

		i__2 = *n - i;
/* Computing MIN */
		i__3 = i + 2;
		dlarfg_(&i__2, &A(i+1,i), &A(min(i+2,*n),i), &c__1, &TAU(i));
		E(i) = A(i+1,i);
		A(i+1,i) = 1.;

/*              Compute W(i+1:n,i) */

		i__2 = *n - i;
		dsymv_("Lower", &i__2, &c_b6, &A(i+1,i+1), 
			lda, &A(i+1,i), &c__1, &c_b16, &W(i+1,i), &c__1);
		i__2 = *n - i;
		i__3 = i - 1;
		dgemv_("Transpose", &i__2, &i__3, &c_b6, &W(i+1,1), 
			ldw, &A(i+1,i), &c__1, &c_b16, &W(1,i), &c__1);
		i__2 = *n - i;
		i__3 = i - 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b5, &A(i+1,1)
			, lda, &W(1,i), &c__1, &c_b6, &W(i+1,i), &c__1);
		i__2 = *n - i;
		i__3 = i - 1;
		dgemv_("Transpose", &i__2, &i__3, &c_b6, &A(i+1,1), 
			lda, &A(i+1,i), &c__1, &c_b16, &W(1,i), &c__1);
		i__2 = *n - i;
		i__3 = i - 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b5, &W(i+1,1)
			, ldw, &W(1,i), &c__1, &c_b6, &W(i+1,i), &c__1);
		i__2 = *n - i;
		dscal_(&i__2, &TAU(i), &W(i+1,i), &c__1);
		i__2 = *n - i;
		alpha = TAU(i) * -.5 * ddot_(&i__2, &W(i+1,i), &
			c__1, &A(i+1,i), &c__1);
		i__2 = *n - i;
		daxpy_(&i__2, &alpha, &A(i+1,i), &c__1, &W(i+1,i), &c__1);
	    }

/* L20: */
	}
    }

    return 0;

/*     End of DLATRD */

} /* dlatrd_ */


/* Subroutine */
static  int dorg2l_(integer *m, integer *n, integer *k, double *
	a, integer *lda, double *tau, double *work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DORG2L generates an m by n float matrix Q with orthonormal columns,   
    which is defined as the last n columns of a product of k elementary   
    reflectors of order m   

          Q  =  H(k) . . . H(2) H(1)   

    as returned by DGEQLF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. M >= N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the 
  
            matrix Q. N >= K >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the (n-k+i)-th column must contain the vector which 
  
            defines the elementary reflector H(i), for i = 1,2,...,k, as 
  
            returned by DGEQLF in the last k columns of its array   
            argument A.   
            On exit, the m by n matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= max(1,M).   

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEQLF.   

    WORK    (workspace) DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument has an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    double d__1;
    /* Local variables */
    static integer i, j, l;
    extern /* Subroutine */ int dscal_(integer *, double *, double *, 
	    integer *), dlarf_(char *, integer *, integer *, double *, 
	    integer *, double *, double *, integer *, double *);
    static integer ii;
    extern /* Subroutine */ int xerbla_(char *, integer *);



#define TAU(I) tau[(I)-1]
#undef WORK
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORG2L", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	return 0;
    }

/*     Initialise columns 1:n-k to columns of the unit matrix */

    i__1 = *n - *k;
    for (j = 1; j <= *n-*k; ++j) {
	i__2 = *m;
	for (l = 1; l <= *m; ++l) {
	    A(l,j) = 0.;
/* L10: */
	}
	A(*m-*n+j,j) = 1.;
/* L20: */
    }

    i__1 = *k;
    for (i = 1; i <= *k; ++i) {
	ii = *n - *k + i;

/*        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left */

	A(*m-*n+ii,ii) = 1.;
	i__2 = *m - *n + ii;
	i__3 = ii - 1;
	dlarf_("Left", &i__2, &i__3, &A(1,ii), &c__1, &TAU(i), &A(1,1), lda, &WORK(1));
	i__2 = *m - *n + ii - 1;
	d__1 = -TAU(i);
	dscal_(&i__2, &d__1, &A(1,ii), &c__1);
	A(*m-*n+ii,ii) = 1. - TAU(i);

/*        Set A(m-k+i+1:m,n-k+i) to zero */

	i__2 = *m;
	for (l = *m - *n + ii + 1; l <= *m; ++l) {
	    A(l,ii) = 0.;
/* L30: */
	}
/* L40: */
    }
    return 0;

/*     End of DORG2L */

} /* dorg2l_ */


/* Subroutine */
static  int dorg2r_(integer *m, integer *n, integer *k, double *
	a, integer *lda, double *tau, double *work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DORG2R generates an m by n float matrix Q with orthonormal columns,   
    which is defined as the first n columns of a product of k elementary 
  
    reflectors of order m   

          Q  =  H(1) H(2) . . . H(k)   

    as returned by DGEQRF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. M >= N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the 
  
            matrix Q. N >= K >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the i-th column must contain the vector which   
            defines the elementary reflector H(i), for i = 1,2,...,k, as 
  
            returned by DGEQRF in the first k columns of its array   
            argument A.   
            On exit, the m-by-n matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= max(1,M).   

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEQRF.   

    WORK    (workspace) DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument has an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    double d__1;
    /* Local variables */
    static integer i, j, l;
    extern /* Subroutine */ int dscal_(integer *, double *, double *, 
	    integer *), dlarf_(char *, integer *, integer *, double *, 
	    integer *, double *, double *, integer *, double *), xerbla_(char *, integer *);



#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORG2R", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	return 0;
    }

/*     Initialise columns k+1:n to columns of the unit matrix */

    i__1 = *n;
    for (j = *k + 1; j <= *n; ++j) {
	i__2 = *m;
	for (l = 1; l <= *m; ++l) {
	    A(l,j) = 0.;
/* L10: */
	}
	A(j,j) = 1.;
/* L20: */
    }

    for (i = *k; i >= 1; --i) {

/*        Apply H(i) to A(i:m,i:n) from the left */

	if (i < *n) {
	    A(i,i) = 1.;
	    i__1 = *m - i + 1;
	    i__2 = *n - i;
	    dlarf_("Left", &i__1, &i__2, &A(i,i), &c__1, &TAU(i), &
		    A(i,i+1), lda, &WORK(1));
	}
	if (i < *m) {
	    i__1 = *m - i;
	    d__1 = -TAU(i);
	    dscal_(&i__1, &d__1, &A(i+1,i), &c__1);
	}
	A(i,i) = 1. - TAU(i);

/*        Set A(1:i-1,i) to zero */

	i__1 = i - 1;
	for (l = 1; l <= i-1; ++l) {
	    A(l,i) = 0.;
/* L30: */
	}
/* L40: */
    }
    return 0;

/*     End of DORG2R */

} /* dorg2r_ */


/* Subroutine */
static  int dorgql_(integer *m, integer *n, integer *k, double *
	a, integer *lda, double *tau, double *work, integer *lwork, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DORGQL generates an M-by-N float matrix Q with orthonormal columns,   
    which is defined as the last N columns of a product of K elementary   
    reflectors of order M   

          Q  =  H(k) . . . H(2) H(1)   

    as returned by DGEQLF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. M >= N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the 
  
            matrix Q. N >= K >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the (n-k+i)-th column must contain the vector which 
  
            defines the elementary reflector H(i), for i = 1,2,...,k, as 
  
            returned by DGEQLF in the last k columns of its array   
            argument A.   
            On exit, the M-by-N matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= max(1,M).   

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEQLF.   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= max(1,N).   
            For optimum performance LWORK >= N*NB, where NB is the   
            optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument has an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__3 = 3;
    static integer c__2 = 2;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    static integer i, j, l, nbmin, iinfo;
    extern /* Subroutine */ int dorg2l_(integer *, integer *, integer *, 
	    double *, integer *, double *, double *, integer *);
    static integer ib, nb, kk;
    extern /* Subroutine */ int dlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, double *, integer *, 
	    double *, integer *, double *, integer *, double *, 
	    integer *);
    static integer nx;
    extern /* Subroutine */ int dlarft_(char *, char *, integer *, integer *, 
	    double *, integer *, double *, double *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ldwork, iws;



#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    } else if (*lwork < max(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORGQL", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	WORK(1) = 1.;
	return 0;
    }

/*     Determine the block size. */

    nb = ilaenv_(&c__1, "DORGQL", " ", m, n, k, &c_n1, 6L, 1L);
    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code.
   

   Computing MAX */
	i__1 = 0, i__2 = ilaenv_(&c__3, "DORGQL", " ", m, n, k, &c_n1, 6L, 1L)
		;
	nx = max(i__1,i__2);
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked co
de. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduc
e NB and   
                determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = ilaenv_(&c__2, "DORGQL", " ", m, n, k, &c_n1,
			 6L, 1L);
		nbmin = max(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the first block.   
          The last kk columns are handled by the block method.   

   Computing MIN */
	i__1 = *k, i__2 = (*k - nx + nb - 1) / nb * nb;
	kk = min(i__1,i__2);

/*        Set A(m-kk+1:m,1:n-kk) to zero. */

	i__1 = *n - kk;
	for (j = 1; j <= *n-kk; ++j) {
	    i__2 = *m;
	    for (i = *m - kk + 1; i <= *m; ++i) {
		A(i,j) = 0.;
/* L10: */
	    }
/* L20: */
	}
    } else {
	kk = 0;
    }

/*     Use unblocked code for the first or only block. */

    i__1 = *m - kk;
    i__2 = *n - kk;
    i__3 = *k - kk;
    dorg2l_(&i__1, &i__2, &i__3, &A(1,1), lda, &TAU(1), &WORK(1), &iinfo)
	    ;

    if (kk > 0) {

/*        Use blocked code */

	i__1 = *k;
	i__2 = nb;
	for (i = *k - kk + 1; nb < 0 ? i >= *k : i <= *k; i += nb) {
/* Computing MIN */
	    i__3 = nb, i__4 = *k - i + 1;
	    ib = min(i__3,i__4);
	    if (*n - *k + i > 1) {

/*              Form the triangular factor of the block reflec
tor   
                H = H(i+ib-1) . . . H(i+1) H(i) */

		i__3 = *m - *k + i + ib - 1;
		dlarft_("Backward", "Columnwise", &i__3, &ib, &A(1,*n-*k+i), lda, &TAU(i), &WORK(1), &ldwork);

/*              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the 
left */

		i__3 = *m - *k + i + ib - 1;
		i__4 = *n - *k + i - 1;
		dlarfb_("Left", "No transpose", "Backward", "Columnwise", &
			i__3, &i__4, &ib, &A(1,*n-*k+i), lda,
			 &WORK(1), &ldwork, &A(1,1), lda, &WORK(ib + 1), 
			&ldwork);
	    }

/*           Apply H to rows 1:m-k+i+ib-1 of current block */

	    i__3 = *m - *k + i + ib - 1;
	    dorg2l_(&i__3, &ib, &ib, &A(1,*n-*k+i), lda, &
		    TAU(i), &WORK(1), &iinfo);

/*           Set rows m-k+i+ib:m of current block to zero */

	    i__3 = *n - *k + i + ib - 1;
	    for (j = *n - *k + i; j <= *n-*k+i+ib-1; ++j) {
		i__4 = *m;
		for (l = *m - *k + i + ib; l <= *m; ++l) {
		    A(l,j) = 0.;
/* L30: */
		}
/* L40: */
	    }
/* L50: */
	}
    }

    WORK(1) = (double) iws;
    return 0;

/*     End of DORGQL */

} /* dorgql_ */


/* Subroutine */
static  int dorgqr_(integer *m, integer *n, integer *k, double *
	a, integer *lda, double *tau, double *work, integer *lwork, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DORGQR generates an M-by-N float matrix Q with orthonormal columns,   
    which is defined as the first N columns of a product of K elementary 
  
    reflectors of order M   

          Q  =  H(1) H(2) . . . H(k)   

    as returned by DGEQRF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. M >= N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the 
  
            matrix Q. N >= K >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the i-th column must contain the vector which   
            defines the elementary reflector H(i), for i = 1,2,...,k, as 
  
            returned by DGEQRF in the first k columns of its array   
            argument A.   
            On exit, the M-by-N matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= max(1,M).   

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEQRF.   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= max(1,N).   
            For optimum performance LWORK >= N*NB, where NB is the   
            optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument has an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__3 = 3;
    static integer c__2 = 2;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer i, j, l, nbmin, iinfo;
    extern /* Subroutine */ int dorg2r_(integer *, integer *, integer *, 
	    double *, integer *, double *, double *, integer *);
    static integer ib, nb, ki, kk;
    extern /* Subroutine */ int dlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, double *, integer *, 
	    double *, integer *, double *, integer *, double *, 
	    integer *);
    static integer nx;
    extern /* Subroutine */ int dlarft_(char *, char *, integer *, integer *, 
	    double *, integer *, double *, double *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ldwork, iws;



#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    } else if (*lwork < max(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORGQR", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	WORK(1) = 1.;
	return 0;
    }

/*     Determine the block size. */

    nb = ilaenv_(&c__1, "DORGQR", " ", m, n, k, &c_n1, 6L, 1L);
    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code.
   

   Computing MAX */
	i__1 = 0, i__2 = ilaenv_(&c__3, "DORGQR", " ", m, n, k, &c_n1, 6L, 1L)
		;
	nx = max(i__1,i__2);
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked co
de. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduc
e NB and   
                determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = ilaenv_(&c__2, "DORGQR", " ", m, n, k, &c_n1,
			 6L, 1L);
		nbmin = max(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the last block.   
          The first kk columns are handled by the block method. */

	ki = (*k - nx - 1) / nb * nb;
/* Computing MIN */
	i__1 = *k, i__2 = ki + nb;
	kk = min(i__1,i__2);

/*        Set A(1:kk,kk+1:n) to zero. */

	i__1 = *n;
	for (j = kk + 1; j <= *n; ++j) {
	    i__2 = kk;
	    for (i = 1; i <= kk; ++i) {
		A(i,j) = 0.;
/* L10: */
	    }
/* L20: */
	}
    } else {
	kk = 0;
    }

/*     Use unblocked code for the last or only block. */

    if (kk < *n) {
	i__1 = *m - kk;
	i__2 = *n - kk;
	i__3 = *k - kk;
	dorg2r_(&i__1, &i__2, &i__3, &A(kk+1,kk+1), lda, &
		TAU(kk + 1), &WORK(1), &iinfo);
    }

    if (kk > 0) {

/*        Use blocked code */

	i__1 = -nb;
	for (i = ki + 1; -nb < 0 ? i >= 1 : i <= 1; i += -nb) {
/* Computing MIN */
	    i__2 = nb, i__3 = *k - i + 1;
	    ib = min(i__2,i__3);
	    if (i + ib <= *n) {

/*              Form the triangular factor of the block reflec
tor   
                H = H(i) H(i+1) . . . H(i+ib-1) */

		i__2 = *m - i + 1;
		dlarft_("Forward", "Columnwise", &i__2, &ib, &A(i,i), lda, &TAU(i), &WORK(1), &ldwork);

/*              Apply H to A(i:m,i+ib:n) from the left */

		i__2 = *m - i + 1;
		i__3 = *n - i - ib + 1;
		dlarfb_("Left", "No transpose", "Forward", "Columnwise", &
			i__2, &i__3, &ib, &A(i,i), lda, &WORK(1), &
			ldwork, &A(i,i+ib), lda, &WORK(ib + 1),
			 &ldwork);
	    }

/*           Apply H to rows i:m of current block */

	    i__2 = *m - i + 1;
	    dorg2r_(&i__2, &ib, &ib, &A(i,i), lda, &TAU(i), &WORK(
		    1), &iinfo);

/*           Set rows 1:i-1 of current block to zero */

	    i__2 = i + ib - 1;
	    for (j = i; j <= i+ib-1; ++j) {
		i__3 = i - 1;
		for (l = 1; l <= i-1; ++l) {
		    A(l,j) = 0.;
/* L30: */
		}
/* L40: */
	    }
/* L50: */
	}
    }

    WORK(1) = (double) iws;
    return 0;

/*     End of DORGQR */

} /* dorgqr_ */


/* Subroutine */
static  int dorgtr_(char *uplo, integer *n, double *a, integer *
	lda, double *tau, double *work, integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DORGTR generates a float orthogonal matrix Q which is defined as the   
    product of n-1 elementary reflectors of order N, as returned by   
    DSYTRD:   

    if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),   

    if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U': Upper triangle of A contains elementary reflectors   
                   from DSYTRD;   
            = 'L': Lower triangle of A contains elementary reflectors   
                   from DSYTRD.   

    N       (input) INTEGER   
            The order of the matrix Q. N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the vectors which define the elementary reflectors, 
  
            as returned by DSYTRD.   
            On exit, the N-by-N orthogonal matrix Q.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= max(1,N).   

    TAU     (input) DOUBLE PRECISION array, dimension (N-1)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DSYTRD.   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= max(1,N-1).   
            For optimum performance LWORK >= (N-1)*NB, where NB is   
            the optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer iinfo;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *), dorgql_(
	    integer *, integer *, integer *, double *, integer *, 
	    double *, double *, integer *, integer *), dorgqr_(
	    integer *, integer *, integer *, double *, integer *, 
	    double *, double *, integer *, integer *);


#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n - 1;
	if (*lwork < max(i__1,i__2)) {
	    *info = -7;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORGTR", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	WORK(1) = 1.;
	return 0;
    }

    if (upper) {

/*        Q was determined by a call to DSYTRD with UPLO = 'U'   

          Shift the vectors which define the elementary reflectors one
   
          column to the left, and set the last row and column of Q to 
  
          those of the unit matrix */

	i__1 = *n - 1;
	for (j = 1; j <= *n-1; ++j) {
	    i__2 = j - 1;
	    for (i = 1; i <= j-1; ++i) {
		A(i,j) = A(i,j+1);
/* L10: */
	    }
	    A(*n,j) = 0.;
/* L20: */
	}
	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {
	    A(i,*n) = 0.;
/* L30: */
	}
	A(*n,*n) = 1.;

/*        Generate Q(1:n-1,1:n-1) */

	i__1 = *n - 1;
	i__2 = *n - 1;
	i__3 = *n - 1;
	dorgql_(&i__1, &i__2, &i__3, &A(1,1), lda, &TAU(1), &WORK(1), 
		lwork, &iinfo);

    } else {

/*        Q was determined by a call to DSYTRD with UPLO = 'L'.   

          Shift the vectors which define the elementary reflectors one
   
          column to the right, and set the first row and column of Q t
o   
          those of the unit matrix */

	for (j = *n; j >= 2; --j) {
	    A(1,j) = 0.;
	    i__1 = *n;
	    for (i = j + 1; i <= *n; ++i) {
		A(i,j) = A(i,j-1);
/* L40: */
	    }
/* L50: */
	}
	A(1,1) = 1.;
	i__1 = *n;
	for (i = 2; i <= *n; ++i) {
	    A(i,1) = 0.;
/* L60: */
	}
	if (*n > 1) {

/*           Generate Q(2:n,2:n) */

	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    i__3 = *n - 1;
	    dorgqr_(&i__1, &i__2, &i__3, &A(2,2), lda, &TAU(1), 
		    &WORK(1), lwork, &iinfo);
	}
    }
    return 0;

/*     End of DORGTR */

} /* dorgtr_ */


/* Subroutine */
static  int dsteqr_(char *compz, integer *n, double *d, 
	double *e, double *z, integer *ldz, double *work, integer 
	*info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSTEQR computes all eigenvalues and, optionally, eigenvectors of a   
    symmetric tridiagonal matrix using the implicit QL or QR method.   
    The eigenvectors of a full or band symmetric matrix can also be found 
  
    if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to 
  
    tridiagonal form.   

    Arguments   
    =========   

    COMPZ   (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only.   
            = 'V':  Compute eigenvalues and eigenvectors of the original 
  
                    symmetric matrix.  On entry, Z must contain the   
                    orthogonal matrix used to reduce the original matrix 
  
                    to tridiagonal form.   
            = 'I':  Compute eigenvalues and eigenvectors of the   
                    tridiagonal matrix.  Z is initialized to the identity 
  
                    matrix.   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, the diagonal elements of the tridiagonal matrix.   
            On exit, if INFO = 0, the eigenvalues in ascending order.   

    E       (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, the (n-1) subdiagonal elements of the tridiagonal   
            matrix.   
            On exit, E has been destroyed.   

    Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)   
            On entry, if  COMPZ = 'V', then Z contains the orthogonal   
            matrix used in the reduction to tridiagonal form.   
            On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the   
            orthonormal eigenvectors of the original symmetric matrix,   
            and if COMPZ = 'I', Z contains the orthonormal eigenvectors   
            of the symmetric tridiagonal matrix.   
            If COMPZ = 'N', then Z is not referenced.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1, and if   
            eigenvectors are desired, then  LDZ >= max(1,N).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2)) 
  
            If COMPZ = 'N', then WORK is not referenced.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  the algorithm has failed to find all the eigenvalues in 
  
                  a total of 30*N iterations; if INFO = i, then i   
                  elements of E have not converged to zero; on exit, D   
                  and E contain the elements of a symmetric tridiagonal   
                  matrix which is orthogonally similar to the original   
                  matrix.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static double c_b9 = 0.;
    static double c_b10 = 1.;
    static integer c__0 = 0;
    static integer c__1 = 1;
    static integer c__2 = 2;
    
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
    double d__1, d__2;
    /* Builtin functions */
    extern double d_sign(double *, double *);
    /* Local variables */
    static integer lend, jtot;
    extern /* Subroutine */ int dlae2_(double *, double *, double 
	    *, double *, double *);
    static double b, c, f, g;
    static integer i, j, k, l, m;
    static double p, r, s;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dlasr_(char *, char *, char *, integer *, 
	    integer *, double *, double *, double *, integer *);
    static double anorm;
    extern /* Subroutine */ int dswap_(integer *, double *, integer *, 
	    double *, integer *);
    static integer l1;
    extern /* Subroutine */ int dlaev2_(double *, double *, 
	    double *, double *, double *, double *, 
	    double *);
    static integer lendm1, lendp1;
    extern double dlapy2_(double *, double *);
    static integer ii;
    extern double dlamch_(char *);
    static integer mm, iscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    double *, double *, integer *, integer *, double *, 
	    integer *, integer *), dlaset_(char *, integer *, integer 
	    *, double *, double *, double *, integer *);
    static double safmin;
    extern /* Subroutine */ int dlartg_(double *, double *, 
	    double *, double *, double *);
    static double safmax;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern double dlanst_(char *, integer *, double *, double *);
    extern /* Subroutine */ int dlasrt_(char *, integer *, double *, 
	    integer *);
    static integer lendsv;
    static double ssfmin;
    static integer nmaxit, icompz;
    static double ssfmax;
    static integer lm1, mm1, nm1;
    static double rt1, rt2, eps;
    static integer lsv;
    static double tst, eps2;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define WORK(I) work[(I)-1]

#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    *info = 0;

    if (lsame_(compz, "N")) {
	icompz = 0;
    } else if (lsame_(compz, "V")) {
	icompz = 1;
    } else if (lsame_(compz, "I")) {
	icompz = 2;
    } else {
	icompz = -1;
    }
    if (icompz < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n)) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSTEQR", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (*n == 1) {
	if (icompz == 2) {
	    Z(1,1) = 1.;
	}
	return 0;
    }

/*     Determine the unit roundoff and over/underflow thresholds. */

    eps = dlamch_("E");
/* Computing 2nd power */
    d__1 = eps;
    eps2 = d__1 * d__1;
    safmin = dlamch_("S");
    safmax = 1. / safmin;
    ssfmax = sqrt(safmax) / 3.;
    ssfmin = sqrt(safmin) / eps2;

/*     Compute the eigenvalues and eigenvectors of the tridiagonal   
       matrix. */

    if (icompz == 2) {
	dlaset_("Full", n, n, &c_b9, &c_b10, &Z(1,1), ldz);
    }

    nmaxit = *n * 30;
    jtot = 0;

/*     Determine where the matrix splits and choose QL or QR iteration   
       for each block, according to whether top or bottom diagonal   
       element is smaller. */

    l1 = 1;
    nm1 = *n - 1;

L10:
    if (l1 > *n) {
	goto L160;
    }
    if (l1 > 1) {
	E(l1 - 1) = 0.;
    }
    if (l1 <= nm1) {
	i__1 = nm1;
	for (m = l1; m <= nm1; ++m) {
	    tst = (d__1 = E(m), abs(d__1));
	    if (tst == 0.) {
		goto L30;
	    }
	    if (tst <= sqrt((d__1 = D(m), abs(d__1))) * sqrt((d__2 = D(m + 1),
		     abs(d__2))) * eps) {
		E(m) = 0.;
		goto L30;
	    }
/* L20: */
	}
    }
    m = *n;

L30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
	goto L10;
    }

/*     Scale submatrix in rows and columns L to LEND */

    i__1 = lend - l + 1;
    anorm = dlanst_("I", &i__1, &D(l), &E(l));
    iscale = 0;
    if (anorm == 0.) {
	goto L10;
    }
    if (anorm > ssfmax) {
	iscale = 1;
	i__1 = lend - l + 1;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &D(l), n, 
		info);
	i__1 = lend - l;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &E(l), n, 
		info);
    } else if (anorm < ssfmin) {
	iscale = 2;
	i__1 = lend - l + 1;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &D(l), n, 
		info);
	i__1 = lend - l;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &E(l), n, 
		info);
    }

/*     Choose between QL and QR iteration */

    if ((d__1 = D(lend), abs(d__1)) < (d__2 = D(l), abs(d__2))) {
	lend = lsv;
	l = lendsv;
    }

    if (lend > l) {

/*        QL Iteration   

          Look for small subdiagonal element. */

L40:
	if (l != lend) {
	    lendm1 = lend - 1;
	    i__1 = lendm1;
	    for (m = l; m <= lendm1; ++m) {
/* Computing 2nd power */
		d__2 = (d__1 = E(m), abs(d__1));
		tst = d__2 * d__2;
		if (tst <= eps2 * (d__1 = D(m), abs(d__1)) * (d__2 = D(m + 1),
			 abs(d__2)) + safmin) {
		    goto L60;
		}
/* L50: */
	    }
	}

	m = lend;

L60:
	if (m < lend) {
	    E(m) = 0.;
	}
	p = D(l);
	if (m == l) {
	    goto L80;
	}

/*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2   
          to compute its eigensystem. */

	if (m == l + 1) {
	    if (icompz > 0) {
		dlaev2_(&D(l), &E(l), &D(l + 1), &rt1, &rt2, &c, &s);
		WORK(l) = c;
		WORK(*n - 1 + l) = s;
		dlasr_("R", "V", "B", n, &c__2, &WORK(l), &WORK(*n - 1 + l), &
			Z(1,l), ldz);
	    } else {
		dlae2_(&D(l), &E(l), &D(l + 1), &rt1, &rt2);
	    }
	    D(l) = rt1;
	    D(l + 1) = rt2;
	    E(l) = 0.;
	    l += 2;
	    if (l <= lend) {
		goto L40;
	    }
	    goto L140;
	}

	if (jtot == nmaxit) {
	    goto L140;
	}
	++jtot;

/*        Form shift. */

	g = (D(l + 1) - p) / (E(l) * 2.);
	r = dlapy2_(&g, &c_b10);
	g = D(m) - p + E(l) / (g + d_sign(&r, &g));

	s = 1.;
	c = 1.;
	p = 0.;

/*        Inner loop */

	mm1 = m - 1;
	i__1 = l;
	for (i = mm1; i >= l; --i) {
	    f = s * E(i);
	    b = c * E(i);
	    dlartg_(&g, &f, &c, &s, &r);
	    if (i != m - 1) {
		E(i + 1) = r;
	    }
	    g = D(i + 1) - p;
	    r = (D(i) - g) * s + c * 2. * b;
	    p = s * r;
	    D(i + 1) = g + p;
	    g = c * r - b;

/*           If eigenvectors are desired, then save rotations. */

	    if (icompz > 0) {
		WORK(i) = c;
		WORK(*n - 1 + i) = -s;
	    }

/* L70: */
	}

/*        If eigenvectors are desired, then apply saved rotations. */

	if (icompz > 0) {
	    mm = m - l + 1;
	    dlasr_("R", "V", "B", n, &mm, &WORK(l), &WORK(*n - 1 + l), &Z(1,l), ldz);
	}

	D(l) -= p;
	E(l) = g;
	goto L40;

/*        Eigenvalue found. */

L80:
	D(l) = p;

	++l;
	if (l <= lend) {
	    goto L40;
	}
	goto L140;

    } else {

/*        QR Iteration   

          Look for small superdiagonal element. */

L90:
	if (l != lend) {
	    lendp1 = lend + 1;
	    i__1 = lendp1;
	    for (m = l; m >= lendp1; --m) {
/* Computing 2nd power */
		d__2 = (d__1 = E(m - 1), abs(d__1));
		tst = d__2 * d__2;
		if (tst <= eps2 * (d__1 = D(m), abs(d__1)) * (d__2 = D(m - 1),
			 abs(d__2)) + safmin) {
		    goto L110;
		}
/* L100: */
	    }
	}

	m = lend;

L110:
	if (m > lend) {
	    E(m - 1) = 0.;
	}
	p = D(l);
	if (m == l) {
	    goto L130;
	}

/*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2   
          to compute its eigensystem. */

	if (m == l - 1) {
	    if (icompz > 0) {
		dlaev2_(&D(l - 1), &E(l - 1), &D(l), &rt1, &rt2, &c, &s);
		WORK(m) = c;
		WORK(*n - 1 + m) = s;
		dlasr_("R", "V", "F", n, &c__2, &WORK(m), &WORK(*n - 1 + m), &
			Z(1,l-1), ldz);
	    } else {
		dlae2_(&D(l - 1), &E(l - 1), &D(l), &rt1, &rt2);
	    }
	    D(l - 1) = rt1;
	    D(l) = rt2;
	    E(l - 1) = 0.;
	    l += -2;
	    if (l >= lend) {
		goto L90;
	    }
	    goto L140;
	}

	if (jtot == nmaxit) {
	    goto L140;
	}
	++jtot;

/*        Form shift. */

	g = (D(l - 1) - p) / (E(l - 1) * 2.);
	r = dlapy2_(&g, &c_b10);
	g = D(m) - p + E(l - 1) / (g + d_sign(&r, &g));

	s = 1.;
	c = 1.;
	p = 0.;

/*        Inner loop */

	lm1 = l - 1;
	i__1 = lm1;
	for (i = m; i <= lm1; ++i) {
	    f = s * E(i);
	    b = c * E(i);
	    dlartg_(&g, &f, &c, &s, &r);
	    if (i != m) {
		E(i - 1) = r;
	    }
	    g = D(i) - p;
	    r = (D(i + 1) - g) * s + c * 2. * b;
	    p = s * r;
	    D(i) = g + p;
	    g = c * r - b;

/*           If eigenvectors are desired, then save rotations. */

	    if (icompz > 0) {
		WORK(i) = c;
		WORK(*n - 1 + i) = s;
	    }

/* L120: */
	}

/*        If eigenvectors are desired, then apply saved rotations. */

	if (icompz > 0) {
	    mm = l - m + 1;
	    dlasr_("R", "V", "F", n, &mm, &WORK(m), &WORK(*n - 1 + m), &Z(1,m), ldz);
	}

	D(l) -= p;
	E(lm1) = g;
	goto L90;

/*        Eigenvalue found. */

L130:
	D(l) = p;

	--l;
	if (l >= lend) {
	    goto L90;
	}
	goto L140;

    }

/*     Undo scaling if necessary */

L140:
    if (iscale == 1) {
	i__1 = lendsv - lsv + 1;
	dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &D(lsv), n, 
		info);
	i__1 = lendsv - lsv;
	dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &E(lsv), n, 
		info);
    } else if (iscale == 2) {
	i__1 = lendsv - lsv + 1;
	dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &D(lsv), n, 
		info);
	i__1 = lendsv - lsv;
	dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &E(lsv), n, 
		info);
    }

/*     Check for no convergence to an eigenvalue after a total   
       of N*MAXIT iterations. */

    if (jtot < nmaxit) {
	goto L10;
    }
    i__1 = *n - 1;
    for (i = 1; i <= *n-1; ++i) {
	if (E(i) != 0.) {
	    ++(*info);
	}
/* L150: */
    }
    goto L190;

/*     Order eigenvalues and eigenvectors. */

L160:
    if (icompz == 0) {

/*        Use Quick Sort */

	dlasrt_("I", n, &D(1), info);

    } else {

/*        Use Selection Sort to minimize swaps of eigenvectors */

	i__1 = *n;
	for (ii = 2; ii <= *n; ++ii) {
	    i = ii - 1;
	    k = i;
	    p = D(i);
	    i__2 = *n;
	    for (j = ii; j <= *n; ++j) {
		if (D(j) < p) {
		    k = j;
		    p = D(j);
		}
/* L170: */
	    }
	    if (k != i) {
		D(k) = D(i);
		D(i) = p;
		dswap_(n, &Z(1,i), &c__1, &Z(1,k), &
			c__1);
	    }
/* L180: */
	}
    }

L190:
    return 0;

/*     End of DSTEQR */

} /* dsteqr_ */


/* Subroutine */
static  int dsterf_(integer *n, double *d, double *e, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSTERF computes all eigenvalues of a symmetric tridiagonal matrix   
    using the Pal-Walker-Kahan variant of the QL or QR algorithm.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, the n diagonal elements of the tridiagonal matrix. 
  
            On exit, if INFO = 0, the eigenvalues in ascending order.   

    E       (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, the (n-1) subdiagonal elements of the tridiagonal   
            matrix.   
            On exit, E has been destroyed.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  the algorithm failed to find all of the eigenvalues in 
  
                  a total of 30*N iterations; if INFO = i, then i   
                  elements of E have not converged to zero.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__0 = 0;
    static integer c__1 = 1;
    static double c_b32 = 1.;
    
    /* System generated locals */
    integer i__1;
    double d__1, d__2;
    /* Builtin functions */
    extern double d_sign(double *, double *);
    /* Local variables */
    static double oldc;
    static integer lend, jtot;
    extern /* Subroutine */ int dlae2_(double *, double *, double 
	    *, double *, double *);
    static double c;
    static integer i, l, m;
    static double p, gamma, r, s, alpha, sigma, anorm;
    static integer l1, lendm1, lendp1;
    extern double dlapy2_(double *, double *);
    static double bb;
    extern double dlamch_(char *);
    static integer iscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    double *, double *, integer *, integer *, double *, 
	    integer *, integer *);
    static double oldgam, safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static double safmax;
    extern double dlanst_(char *, integer *, double *, double *);
    extern /* Subroutine */ int dlasrt_(char *, integer *, double *, 
	    integer *);
    static integer lendsv;
    static double ssfmin;
    static integer nmaxit;
    static double ssfmax;
    static integer lm1, mm1, nm1;
    static double rt1, rt2, eps, rte;
    static integer lsv;
    static double tst, eps2;



#define E(I) e[(I)-1]
#define D(I) d[(I)-1]


    *info = 0;

/*     Quick return if possible */

    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	xerbla_("DSTERF", &i__1);
	return 0;
    }
    if (*n <= 1) {
	return 0;
    }

/*     Determine the unit roundoff for this environment. */

    eps = dlamch_("E");
/* Computing 2nd power */
    d__1 = eps;
    eps2 = d__1 * d__1;
    safmin = dlamch_("S");
    safmax = 1. / safmin;
    ssfmax = sqrt(safmax) / 3.;
    ssfmin = sqrt(safmin) / eps2;

/*     Compute the eigenvalues of the tridiagonal matrix. */

    nmaxit = *n * 30;
    sigma = 0.;
    jtot = 0;

/*     Determine where the matrix splits and choose QL or QR iteration   
       for each block, according to whether top or bottom diagonal   
       element is smaller. */

    l1 = 1;
    nm1 = *n - 1;

L10:
    if (l1 > *n) {
	goto L170;
    }
    if (l1 > 1) {
	E(l1 - 1) = 0.;
    }
    if (l1 <= nm1) {
	i__1 = nm1;
	for (m = l1; m <= nm1; ++m) {
	    tst = (d__1 = E(m), abs(d__1));
	    if (tst == 0.) {
		goto L30;
	    }
	    if (tst <= sqrt((d__1 = D(m), abs(d__1))) * sqrt((d__2 = D(m + 1),
		     abs(d__2))) * eps) {
		E(m) = 0.;
		goto L30;
	    }
/* L20: */
	}
    }
    m = *n;

L30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
	goto L10;
    }

/*     Scale submatrix in rows and columns L to LEND */

    i__1 = lend - l + 1;
    anorm = dlanst_("I", &i__1, &D(l), &E(l));
    iscale = 0;
    if (anorm > ssfmax) {
	iscale = 1;
	i__1 = lend - l + 1;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &D(l), n, 
		info);
	i__1 = lend - l;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &E(l), n, 
		info);
    } else if (anorm < ssfmin) {
	iscale = 2;
	i__1 = lend - l + 1;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &D(l), n, 
		info);
	i__1 = lend - l;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &E(l), n, 
		info);
    }

    i__1 = lend - 1;
    for (i = l; i <= lend-1; ++i) {
/* Computing 2nd power */
	d__1 = E(i);
	E(i) = d__1 * d__1;
/* L40: */
    }

/*     Choose between QL and QR iteration */

    if ((d__1 = D(lend), abs(d__1)) < (d__2 = D(l), abs(d__2))) {
	lend = lsv;
	l = lendsv;
    }

    if (lend >= l) {

/*        QL Iteration   

          Look for small subdiagonal element. */

L50:
	if (l != lend) {
	    lendm1 = lend - 1;
	    i__1 = lendm1;
	    for (m = l; m <= lendm1; ++m) {
		tst = (d__1 = E(m), abs(d__1));
		if (tst <= eps2 * (d__1 = D(m) * D(m + 1), abs(d__1))) {
		    goto L70;
		}
/* L60: */
	    }
	}

	m = lend;

L70:
	if (m < lend) {
	    E(m) = 0.;
	}
	p = D(l);
	if (m == l) {
	    goto L90;
	}

/*        If remaining matrix is 2 by 2, use DLAE2 to compute its   
          eigenvalues. */

	if (m == l + 1) {
	    rte = sqrt(E(l));
	    dlae2_(&D(l), &rte, &D(l + 1), &rt1, &rt2);
	    D(l) = rt1;
	    D(l + 1) = rt2;
	    E(l) = 0.;
	    l += 2;
	    if (l <= lend) {
		goto L50;
	    }
	    goto L150;
	}

	if (jtot == nmaxit) {
	    goto L150;
	}
	++jtot;

/*        Form shift. */

	rte = sqrt(E(l));
	sigma = (D(l + 1) - p) / (rte * 2.);
	r = dlapy2_(&sigma, &c_b32);
	sigma = p - rte / (sigma + d_sign(&r, &sigma));

	c = 1.;
	s = 0.;
	gamma = D(m) - sigma;
	p = gamma * gamma;

/*        Inner loop */

	mm1 = m - 1;
	i__1 = l;
	for (i = mm1; i >= l; --i) {
	    bb = E(i);
	    r = p + bb;
	    if (i != m - 1) {
		E(i + 1) = s * r;
	    }
	    oldc = c;
	    c = p / r;
	    s = bb / r;
	    oldgam = gamma;
	    alpha = D(i);
	    gamma = c * (alpha - sigma) - s * oldgam;
	    D(i + 1) = oldgam + (alpha - gamma);
	    if (c != 0.) {
		p = gamma * gamma / c;
	    } else {
		p = oldc * bb;
	    }
/* L80: */
	}

	E(l) = s * p;
	D(l) = sigma + gamma;
	goto L50;

/*        Eigenvalue found. */

L90:
	D(l) = p;

	++l;
	if (l <= lend) {
	    goto L50;
	}
	goto L150;

    } else {

/*        QR Iteration   

          Look for small superdiagonal element. */

L100:
	if (l != lend) {
	    lendp1 = lend + 1;
	    i__1 = lendp1;
	    for (m = l; m >= lendp1; --m) {
		tst = (d__1 = E(m - 1), abs(d__1));
		if (tst <= eps2 * (d__1 = D(m) * D(m - 1), abs(d__1))) {
		    goto L120;
		}
/* L110: */
	    }
	}

	m = lend;

L120:
	if (m > lend) {
	    E(m - 1) = 0.;
	}
	p = D(l);
	if (m == l) {
	    goto L140;
	}

/*        If remaining matrix is 2 by 2, use DLAE2 to compute its   
          eigenvalues. */

	if (m == l - 1) {
	    rte = sqrt(E(l - 1));
	    dlae2_(&D(l), &rte, &D(l - 1), &rt1, &rt2);
	    D(l) = rt1;
	    D(l - 1) = rt2;
	    E(l - 1) = 0.;
	    l += -2;
	    if (l >= lend) {
		goto L100;
	    }
	    goto L150;
	}

	if (jtot == nmaxit) {
	    goto L150;
	}
	++jtot;

/*        Form shift. */

	rte = sqrt(E(l - 1));
	sigma = (D(l - 1) - p) / (rte * 2.);
	r = dlapy2_(&sigma, &c_b32);
	sigma = p - rte / (sigma + d_sign(&r, &sigma));

	c = 1.;
	s = 0.;
	gamma = D(m) - sigma;
	p = gamma * gamma;

/*        Inner loop */

	lm1 = l - 1;
	i__1 = lm1;
	for (i = m; i <= lm1; ++i) {
	    bb = E(i);
	    r = p + bb;
	    if (i != m) {
		E(i - 1) = s * r;
	    }
	    oldc = c;
	    c = p / r;
	    s = bb / r;
	    oldgam = gamma;
	    alpha = D(i + 1);
	    gamma = c * (alpha - sigma) - s * oldgam;
	    D(i) = oldgam + (alpha - gamma);
	    if (c != 0.) {
		p = gamma * gamma / c;
	    } else {
		p = oldc * bb;
	    }
/* L130: */
	}

	E(lm1) = s * p;
	D(l) = sigma + gamma;
	goto L100;

/*        Eigenvalue found. */

L140:
	D(l) = p;

	--l;
	if (l >= lend) {
	    goto L100;
	}
	goto L150;

    }

/*     Undo scaling if necessary */

L150:
    if (iscale == 1) {
	i__1 = lendsv - lsv + 1;
	dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &D(lsv), n, 
		info);
    }
    if (iscale == 2) {
	i__1 = lendsv - lsv + 1;
	dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &D(lsv), n, 
		info);
    }

/*     Check for no convergence to an eigenvalue after a total   
       of N*MAXIT iterations. */

    if (jtot == nmaxit) {
	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {
	    if (E(i) != 0.) {
		++(*info);
	    }
/* L160: */
	}
	return 0;
    }
    goto L10;

/*     Sort eigenvalues in increasing order. */

L170:
    dlasrt_("I", n, &D(1), info);

    return 0;

/*     End of DSTERF */

} /* dsterf_ */


/* Subroutine */
static  int dsytd2_(char *uplo, integer *n, double *a, integer *
	lda, double *d, double *e, double *tau, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DSYTD2 reduces a float symmetric matrix A to symmetric tridiagonal   
    form T by an orthogonal similarity transformation: Q' * A * Q = T.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is stored:   
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading 
  
            n-by-n upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n-by-n lower triangular part of A contains the lower 
  
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   
            On exit, if UPLO = 'U', the diagonal and first superdiagonal 
  
            of A are overwritten by the corresponding elements of the   
            tridiagonal matrix T, and the elements above the first   
            superdiagonal, with the array TAU, represent the orthogonal   
            matrix Q as a product of elementary reflectors; if UPLO   
            = 'L', the diagonal and first subdiagonal of A are over-   
            written by the corresponding elements of the tridiagonal   
            matrix T, and the elements below the first subdiagonal, with 
  
            the array TAU, represent the orthogonal matrix Q as a product 
  
            of elementary reflectors. See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    D       (output) DOUBLE PRECISION array, dimension (N)   
            The diagonal elements of the tridiagonal matrix T:   
            D(i) = A(i,i).   

    E       (output) DOUBLE PRECISION array, dimension (N-1)   
            The off-diagonal elements of the tridiagonal matrix T:   
            E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'. 
  

    TAU     (output) DOUBLE PRECISION array, dimension (N-1)   
            The scalar factors of the elementary reflectors (see Further 
  
            Details).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    If UPLO = 'U', the matrix Q is represented as a product of elementary 
  
    reflectors   

       Q = H(n-1) . . . H(2) H(1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a float scalar, and v is a float vector with   
    v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in   
    A(1:i-1,i+1), and tau in TAU(i).   

    If UPLO = 'L', the matrix Q is represented as a product of elementary 
  
    reflectors   

       Q = H(1) H(2) . . . H(n-1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a float scalar, and v is a float vector with   
    v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i), 
  
    and tau in TAU(i).   

    The contents of A on exit are illustrated by the following examples   
    with n = 5:   

    if UPLO = 'U':                       if UPLO = 'L':   

      (  d   e   v2  v3  v4 )              (  d                  )   
      (      d   e   v3  v4 )              (  e   d              )   
      (          d   e   v4 )              (  v1  e   d          )   
      (              d   e  )              (  v1  v2  e   d      )   
      (                  d  )              (  v1  v2  v3  e   d  )   

    where d and e denote diagonal and off-diagonal elements of T, and vi 
  
    denotes an element of the vector defining H(i).   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static double c_b8 = 0.;
    static double c_b14 = -1.;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    extern double ddot_(integer *, double *, integer *, double *, 
	    integer *);
    static double taui;
    extern /* Subroutine */ int dsyr2_(char *, integer *, double *, 
	    double *, integer *, double *, integer *, double *, 
	    integer *);
    static integer i;
    static double alpha;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int daxpy_(integer *, double *, double *, 
	    integer *, double *, integer *);
    static logical upper;
    extern /* Subroutine */ int dsymv_(char *, integer *, double *, 
	    double *, integer *, double *, integer *, double *, 
	    double *, integer *), dlarfg_(integer *, double *,
	     double *, integer *, double *), xerbla_(char *, integer *
	    );



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define TAU(I) tau[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYTD2", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	return 0;
    }

    if (upper) {

/*        Reduce the upper triangle of A */

	for (i = *n - 1; i >= 1; --i) {

/*           Generate elementary reflector H(i) = I - tau * v * v'
   
             to annihilate A(1:i-1,i+1) */

	    dlarfg_(&i, &A(i,i+1), &A(1,i+1), &
		    c__1, &taui);
	    E(i) = A(i,i+1);

	    if (taui != 0.) {

/*              Apply H(i) from both sides to A(1:i,1:i) */

		A(i,i+1) = 1.;

/*              Compute  x := tau * A * v  storing x in TAU(1:
i) */

		dsymv_(uplo, &i, &taui, &A(1,1), lda, &A(1,i+1), &c__1, &c_b8, &TAU(1), &c__1);

/*              Compute  w := x - 1/2 * tau * (x'*v) * v */

		alpha = taui * -.5 * ddot_(&i, &TAU(1), &c__1, &A(1,i+1), &c__1);
		daxpy_(&i, &alpha, &A(1,i+1), &c__1, &TAU(1), &
			c__1);

/*              Apply the transformation as a rank-2 update: 
  
                   A := A - v * w' - w * v' */

		dsyr2_(uplo, &i, &c_b14, &A(1,i+1), &c__1, &
			TAU(1), &c__1, &A(1,1), lda);

		A(i,i+1) = E(i);
	    }
	    D(i + 1) = A(i+1,i+1);
	    TAU(i) = taui;
/* L10: */
	}
	D(1) = A(1,1);
    } else {

/*        Reduce the lower triangle of A */

	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {

/*           Generate elementary reflector H(i) = I - tau * v * v'
   
             to annihilate A(i+2:n,i) */

	    i__2 = *n - i;
/* Computing MIN */
	    i__3 = i + 2;
	    dlarfg_(&i__2, &A(i+1,i), &A(min(i+2,*n),i), &c__1, &taui);
	    E(i) = A(i+1,i);

	    if (taui != 0.) {

/*              Apply H(i) from both sides to A(i+1:n,i+1:n) 
*/

		A(i+1,i) = 1.;

/*              Compute  x := tau * A * v  storing y in TAU(i:
n-1) */

		i__2 = *n - i;
		dsymv_(uplo, &i__2, &taui, &A(i+1,i+1), lda, 
			&A(i+1,i), &c__1, &c_b8, &TAU(i), &c__1);

/*              Compute  w := x - 1/2 * tau * (x'*v) * v */

		i__2 = *n - i;
		alpha = taui * -.5 * ddot_(&i__2, &TAU(i), &c__1, &A(i+1,i), &c__1);
		i__2 = *n - i;
		daxpy_(&i__2, &alpha, &A(i+1,i), &c__1, &TAU(i), 
			&c__1);

/*              Apply the transformation as a rank-2 update: 
  
                   A := A - v * w' - w * v' */

		i__2 = *n - i;
		dsyr2_(uplo, &i__2, &c_b14, &A(i+1,i), &c__1, &
			TAU(i), &c__1, &A(i+1,i+1), lda);

		A(i+1,i) = E(i);
	    }
	    D(i) = A(i,i);
	    TAU(i) = taui;
/* L20: */
	}
	D(*n) = A(*n,*n);
    }

    return 0;

/*     End of DSYTD2 */

} /* dsytd2_ */


/* Subroutine */
static  int dsytrd_(char *uplo, integer *n, double *a, integer *
	lda, double *d, double *e, double *tau, double *work, 
	integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSYTRD reduces a float symmetric matrix A to float symmetric   
    tridiagonal form T by an orthogonal similarity transformation:   
    Q**T * A * Q = T.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading 
  
            N-by-N upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of A contains the lower 
  
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   
            On exit, if UPLO = 'U', the diagonal and first superdiagonal 
  
            of A are overwritten by the corresponding elements of the   
            tridiagonal matrix T, and the elements above the first   
            superdiagonal, with the array TAU, represent the orthogonal   
            matrix Q as a product of elementary reflectors; if UPLO   
            = 'L', the diagonal and first subdiagonal of A are over-   
            written by the corresponding elements of the tridiagonal   
            matrix T, and the elements below the first subdiagonal, with 
  
            the array TAU, represent the orthogonal matrix Q as a product 
  
            of elementary reflectors. See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    D       (output) DOUBLE PRECISION array, dimension (N)   
            The diagonal elements of the tridiagonal matrix T:   
            D(i) = A(i,i).   

    E       (output) DOUBLE PRECISION array, dimension (N-1)   
            The off-diagonal elements of the tridiagonal matrix T:   
            E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'. 
  

    TAU     (output) DOUBLE PRECISION array, dimension (N-1)   
            The scalar factors of the elementary reflectors (see Further 
  
            Details).   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= 1.   
            For optimum performance LWORK >= N*NB, where NB is the   
            optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    If UPLO = 'U', the matrix Q is represented as a product of elementary 
  
    reflectors   

       Q = H(n-1) . . . H(2) H(1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a float scalar, and v is a float vector with   
    v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in   
    A(1:i-1,i+1), and tau in TAU(i).   

    If UPLO = 'L', the matrix Q is represented as a product of elementary 
  
    reflectors   

       Q = H(1) H(2) . . . H(n-1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a float scalar, and v is a float vector with   
    v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i), 
  
    and tau in TAU(i).   

    The contents of A on exit are illustrated by the following examples   
    with n = 5:   

    if UPLO = 'U':                       if UPLO = 'L':   

      (  d   e   v2  v3  v4 )              (  d                  )   
      (      d   e   v3  v4 )              (  e   d              )   
      (          d   e   v4 )              (  v1  e   d          )   
      (              d   e  )              (  v1  v2  e   d      )   
      (                  d  )              (  v1  v2  v3  e   d  )   

    where d and e denote diagonal and off-diagonal elements of T, and vi 
  
    denotes an element of the vector defining H(i).   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__3 = 3;
    static integer c__2 = 2;
    static double c_b22 = -1.;
    static double c_b23 = 1.;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer nbmin, iinfo;
    static logical upper;
    extern /* Subroutine */ int dsytd2_(char *, integer *, double *, 
	    integer *, double *, double *, double *, integer *), dsyr2k_(char *, char *, integer *, integer *, double 
	    *, double *, integer *, double *, integer *, double *,
	     double *, integer *);
    static integer nb, kk, nx;
    extern /* Subroutine */ int dlatrd_(char *, integer *, integer *, 
	    double *, integer *, double *, double *, double *,
	     integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ldwork, iws;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else if (*lwork < 1) {
	*info = -9;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYTRD", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	WORK(1) = 1.;
	return 0;
    }

/*     Determine the block size. */

    nb = ilaenv_(&c__1, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, 6L, 1L);
    nx = *n;
    iws = 1;
    if (nb > 1 && nb < *n) {

/*        Determine when to cross over from blocked to unblocked code 
  
          (last block is always handled by unblocked code).   

   Computing MAX */
	i__1 = nb, i__2 = ilaenv_(&c__3, "DSYTRD", uplo, n, &c_n1, &c_n1, &
		c_n1, 6L, 1L);
	nx = max(i__1,i__2);
	if (nx < *n) {

/*           Determine if workspace is large enough for blocked co
de. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  deter
mine the   
                minimum value of NB, and reduce NB or force us
e of   
                unblocked code by setting NX = N.   

   Computing MAX */
		i__1 = *lwork / ldwork;
		nb = max(i__1,1);
		nbmin = ilaenv_(&c__2, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1,
			 6L, 1L);
		if (nb < nbmin) {
		    nx = *n;
		}
	    }
	} else {
	    nx = *n;
	}
    } else {
	nb = 1;
    }

    if (upper) {

/*        Reduce the upper triangle of A.   
          Columns 1:kk are handled by the unblocked method. */

	kk = *n - (*n - nx + nb - 1) / nb * nb;
	i__1 = kk + 1;
	i__2 = -nb;
	for (i = *n - nb + 1; -nb < 0 ? i >= kk+1 : i <= kk+1; i += -nb) {

/*           Reduce columns i:i+nb-1 to tridiagonal form and form 
the   
             matrix W which is needed to update the unreduced part
 of   
             the matrix */

	    i__3 = i + nb - 1;
	    dlatrd_(uplo, &i__3, &nb, &A(1,1), lda, &E(1), &TAU(1), &
		    WORK(1), &ldwork);

/*           Update the unreduced submatrix A(1:i-1,1:i-1), using 
an   
             update of the form:  A := A - V*W' - W*V' */

	    i__3 = i - 1;
	    dsyr2k_(uplo, "No transpose", &i__3, &nb, &c_b22, &A(1,i), lda, &WORK(1), &ldwork, &c_b23, &A(1,1), lda);

/*           Copy superdiagonal elements back into A, and diagonal
   
             elements into D */

	    i__3 = i + nb - 1;
	    for (j = i; j <= i+nb-1; ++j) {
		A(j-1,j) = E(j - 1);
		D(j) = A(j,j);
/* L10: */
	    }
/* L20: */
	}

/*        Use unblocked code to reduce the last or only block */

	dsytd2_(uplo, &kk, &A(1,1), lda, &D(1), &E(1), &TAU(1), &iinfo);
    } else {

/*        Reduce the lower triangle of A */

	i__2 = *n - nx;
	i__1 = nb;
	for (i = 1; nb < 0 ? i >= *n-nx : i <= *n-nx; i += nb) {

/*           Reduce columns i:i+nb-1 to tridiagonal form and form 
the   
             matrix W which is needed to update the unreduced part
 of   
             the matrix */

	    i__3 = *n - i + 1;
	    dlatrd_(uplo, &i__3, &nb, &A(i,i), lda, &E(i), &TAU(i),
		     &WORK(1), &ldwork);

/*           Update the unreduced submatrix A(i+ib:n,i+ib:n), usin
g   
             an update of the form:  A := A - V*W' - W*V' */

	    i__3 = *n - i - nb + 1;
	    dsyr2k_(uplo, "No transpose", &i__3, &nb, &c_b22, &A(i+nb,i), lda, &WORK(nb + 1), &ldwork, &c_b23, &A(i+nb,i+nb), lda);

/*           Copy subdiagonal elements back into A, and diagonal 
  
             elements into D */

	    i__3 = i + nb - 1;
	    for (j = i; j <= i+nb-1; ++j) {
		A(j+1,j) = E(j);
		D(j) = A(j,j);
/* L30: */
	    }
/* L40: */
	}

/*        Use unblocked code to reduce the last or only block */

	i__1 = *n - i + 1;
	dsytd2_(uplo, &i__1, &A(i,i), lda, &D(i), &E(i), &TAU(i), &
		iinfo);
    }

    WORK(1) = (double) iws;
    return 0;

/*     End of DSYTRD */

} /* dsytrd_ */


static integer ilaenv_(integer *ispec, char *name, char *opts, integer *n1, integer *
	n2, integer *n3, integer *n4, ftnlen name_len, ftnlen opts_len)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ILAENV is called from the LAPACK routines to choose problem-dependent 
  
    parameters for the local environment.  See ISPEC for a description of 
  
    the parameters.   

    This version provides a set of parameters which should give good,   
    but not optimal, performance on many of the currently available   
    computers.  Users are encouraged to modify this subroutine to set   
    the tuning parameters for their particular machine using the option   
    and problem size information in the arguments.   

    This routine will not function correctly if it is converted to all   
    lower case.  Converting it to all upper case is allowed.   

    Arguments   
    =========   

    ISPEC   (input) INTEGER   
            Specifies the parameter to be returned as the value of   
            ILAENV.   
            = 1: the optimal blocksize; if this value is 1, an unblocked 
  
                 algorithm will give the best performance.   
            = 2: the minimum block size for which the block routine   
                 should be used; if the usable block size is less than   
                 this value, an unblocked routine should be used.   
            = 3: the crossover point (in a block routine, for N less   
                 than this value, an unblocked routine should be used)   
            = 4: the number of shifts, used in the nonsymmetric   
                 eigenvalue routines   
            = 5: the minimum column dimension for blocking to be used;   
                 rectangular blocks must have dimension at least k by m, 
  
                 where k is given by ILAENV(2,...) and m by ILAENV(5,...) 
  
            = 6: the crossover point for the SVD (when reducing an m by n 
  
                 matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds 
  
                 this value, a QR factorization is used first to reduce   
                 the matrix to a triangular form.)   
            = 7: the number of processors   
            = 8: the crossover point for the multishift QR and QZ methods 
  
                 for nonsymmetric eigenvalue problems.   

    NAME    (input) CHARACTER*(*)   
            The name of the calling subroutine, in either upper case or   
            lower case.   

    OPTS    (input) CHARACTER*(*)   
            The character options to the subroutine NAME, concatenated   
            into a single character string.  For example, UPLO = 'U',   
            TRANS = 'T', and DIAG = 'N' for a triangular routine would   
            be specified as OPTS = 'UTN'.   

    N1      (input) INTEGER   
    N2      (input) INTEGER   
    N3      (input) INTEGER   
    N4      (input) INTEGER   
            Problem dimensions for the subroutine NAME; these may not all 
  
            be required.   

   (ILAENV) (output) INTEGER   
            >= 0: the value of the parameter specified by ISPEC   
            < 0:  if ILAENV = -k, the k-th argument had an illegal value. 
  

    Further Details   
    ===============   

    The following conventions have been used when calling ILAENV from the 
  
    LAPACK routines:   
    1)  OPTS is a concatenation of all of the character options to   
        subroutine NAME, in the same order that they appear in the   
        argument list for NAME, even if they are not used in determining 
  
        the value of the parameter specified by ISPEC.   
    2)  The problem dimensions N1, N2, N3, N4 are specified in the order 
  
        that they appear in the argument list for NAME.  N1 is used   
        first, N2 second, and so on, and unused problem dimensions are   
        passed a value of -1.   
    3)  The parameter value returned by ILAENV is checked for validity in 
  
        the calling subroutine.  For example, ILAENV is used to retrieve 
  
        the optimal blocksize for STRTRI as follows:   

        NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )   
        IF( NB.LE.1 ) NB = MAX( 1, N )   

    ===================================================================== 
*/
/* >>Start of File<<   
       System generated locals */
    integer ret_val;
    /* Builtin functions   
       Subroutine */ extern int s_copy(char *, char *, ftnlen, ftnlen);
    extern integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Local variables */
    static integer i;
    static logical cname, sname;
    static integer nbmin;
    static char c1[1], c2[2], c3[3], c4[2];
    static integer ic, nb, iz, nx;
    static char subnam[6];



    switch (*ispec) {
	case 1:  goto L100;
	case 2:  goto L100;
	case 3:  goto L100;
	case 4:  goto L400;
	case 5:  goto L500;
	case 6:  goto L600;
	case 7:  goto L700;
	case 8:  goto L800;
    }

/*     Invalid value for ISPEC */

    ret_val = -1;
    return ret_val;

L100:

/*     Convert NAME to upper case if the first character is lower case. */

    ret_val = 1;
    s_copy(subnam, name, 6L, name_len);
    ic = *(unsigned char *)subnam;
    iz = 'Z';
    if (iz == 90 || iz == 122) {

/*        ASCII character set */

	if (ic >= 97 && ic <= 122) {
	    *(unsigned char *)subnam = (char) (ic - 32);
	    for (i = 2; i <= 6; ++i) {
		ic = *(unsigned char *)&subnam[i - 1];
		if (ic >= 97 && ic <= 122) {
		    *(unsigned char *)&subnam[i - 1] = (char) (ic - 32);
		}
/* L10: */
	    }
	}

    } else if (iz == 233 || iz == 169) {

/*        EBCDIC character set */

	if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 && 
		ic <= 169) {
	    *(unsigned char *)subnam = (char) (ic + 64);
	    for (i = 2; i <= 6; ++i) {
		ic = *(unsigned char *)&subnam[i - 1];
		if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 
			162 && ic <= 169) {
		    *(unsigned char *)&subnam[i - 1] = (char) (ic + 64);
		}
/* L20: */
	    }
	}

    } else if (iz == 218 || iz == 250) {

/*        Prime machines:  ASCII+128 */

	if (ic >= 225 && ic <= 250) {
	    *(unsigned char *)subnam = (char) (ic - 32);
	    for (i = 2; i <= 6; ++i) {
		ic = *(unsigned char *)&subnam[i - 1];
		if (ic >= 225 && ic <= 250) {
		    *(unsigned char *)&subnam[i - 1] = (char) (ic - 32);
		}
/* L30: */
	    }
	}
    }

    *(unsigned char *)c1 = *(unsigned char *)subnam;
    sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
    cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
    if (! (cname || sname)) {
	return ret_val;
    }
    s_copy(c2, subnam + 1, 2L, 2L);
    s_copy(c3, subnam + 3, 3L, 3L);
    s_copy(c4, c3 + 1, 2L, 2L);

    switch (*ispec) {
	case 1:  goto L110;
	case 2:  goto L200;
	case 3:  goto L300;
    }

L110:

/*     ISPEC = 1:  block size   

       In these examples, separate code is provided for setting NB for   
       float and complex.  We assume that NB will take the same value in   
       single or double precision. */

    nb = 1;

    if (s_cmp(c2, "GE", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	} else if (s_cmp(c3, "QRF", 3L, 3L) == 0 || s_cmp(c3, "RQF", 3L, 3L) 
		== 0 || s_cmp(c3, "LQF", 3L, 3L) == 0 || s_cmp(c3, "QLF", 3L, 
		3L) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "HRD", 3L, 3L) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "BRD", 3L, 3L) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "TRI", 3L, 3L) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "PO", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "SY", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	} else if (sname && s_cmp(c3, "TRD", 3L, 3L) == 0) {
	    nb = 1;
	} else if (sname && s_cmp(c3, "GST", 3L, 3L) == 0) {
	    nb = 64;
	}
    } else if (cname && s_cmp(c2, "HE", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    nb = 64;
	} else if (s_cmp(c3, "TRD", 3L, 3L) == 0) {
	    nb = 1;
	} else if (s_cmp(c3, "GST", 3L, 3L) == 0) {
	    nb = 64;
	}
    } else if (sname && s_cmp(c2, "OR", 2L, 2L) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nb = 32;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nb = 32;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", 2L, 2L) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nb = 32;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nb = 32;
	    }
	}
    } else if (s_cmp(c2, "GB", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    if (sname) {
		if (*n4 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    } else {
		if (*n4 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    }
	}
    } else if (s_cmp(c2, "PB", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    if (sname) {
		if (*n2 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    } else {
		if (*n2 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    }
	}
    } else if (s_cmp(c2, "TR", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRI", 3L, 3L) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "LA", 2L, 2L) == 0) {
	if (s_cmp(c3, "UUM", 3L, 3L) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (sname && s_cmp(c2, "ST", 2L, 2L) == 0) {
	if (s_cmp(c3, "EBZ", 3L, 3L) == 0) {
	    nb = 1;
	}
    }
    ret_val = nb;
    return ret_val;

L200:

/*     ISPEC = 2:  minimum block size */

    nbmin = 2;
    if (s_cmp(c2, "GE", 2L, 2L) == 0) {
	if (s_cmp(c3, "QRF", 3L, 3L) == 0 || s_cmp(c3, "RQF", 3L, 3L) == 0 || 
		s_cmp(c3, "LQF", 3L, 3L) == 0 || s_cmp(c3, "QLF", 3L, 3L) == 
		0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "HRD", 3L, 3L) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "BRD", 3L, 3L) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "TRI", 3L, 3L) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	}
    } else if (s_cmp(c2, "SY", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    if (sname) {
		nbmin = 8;
	    } else {
		nbmin = 8;
	    }
	} else if (sname && s_cmp(c3, "TRD", 3L, 3L) == 0) {
	    nbmin = 2;
	}
    } else if (cname && s_cmp(c2, "HE", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRD", 3L, 3L) == 0) {
	    nbmin = 2;
	}
    } else if (sname && s_cmp(c2, "OR", 2L, 2L) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nbmin = 2;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nbmin = 2;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", 2L, 2L) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nbmin = 2;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nbmin = 2;
	    }
	}
    }
    ret_val = nbmin;
    return ret_val;

L300:

/*     ISPEC = 3:  crossover point */

    nx = 0;
    if (s_cmp(c2, "GE", 2L, 2L) == 0) {
	if (s_cmp(c3, "QRF", 3L, 3L) == 0 || s_cmp(c3, "RQF", 3L, 3L) == 0 || 
		s_cmp(c3, "LQF", 3L, 3L) == 0 || s_cmp(c3, "QLF", 3L, 3L) == 
		0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	} else if (s_cmp(c3, "HRD", 3L, 3L) == 0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	} else if (s_cmp(c3, "BRD", 3L, 3L) == 0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	}
    } else if (s_cmp(c2, "SY", 2L, 2L) == 0) {
	if (sname && s_cmp(c3, "TRD", 3L, 3L) == 0) {
	    nx = 1;
	}
    } else if (cname && s_cmp(c2, "HE", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRD", 3L, 3L) == 0) {
	    nx = 1;
	}
    } else if (sname && s_cmp(c2, "OR", 2L, 2L) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nx = 128;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", 2L, 2L) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nx = 128;
	    }
	}
    }
    ret_val = nx;
    return ret_val;

L400:

/*     ISPEC = 4:  number of shifts (used by xHSEQR) */

    ret_val = 6;
    return ret_val;

L500:

/*     ISPEC = 5:  minimum column dimension (not used) */

    ret_val = 2;
    return ret_val;

L600:

/*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */

    ret_val = (integer) ((float) min(*n1,*n2) * 1.6f);
    return ret_val;

L700:

/*     ISPEC = 7:  number of processors (not used) */

    ret_val = 1;
    return ret_val;

L800:

/*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

    ret_val = 50;
    return ret_val;

/*     End of ILAENV */

} /* ilaenv_ */


static logical lsame_(char *ca, char *cb)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    LSAME returns .TRUE. if CA is the same letter as CB regardless of   
    case.   

    Arguments   
    =========   

    CA      (input) CHARACTER*1   
    CB      (input) CHARACTER*1   
            CA and CB specify the single characters to be compared.   

   ===================================================================== 
  


       Test if the characters are equal */
    /* System generated locals */
    logical ret_val;
    /* Local variables */
    static integer inta, intb, zcode;


    ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
    if (ret_val) {
	return ret_val;
    }

/*     Now test for equivalence if both characters are alphabetic. */

    zcode = 'Z';

/*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime   
       machines, on which ICHAR returns a value with bit 8 set.   
       ICHAR('A') on Prime machines returns 193 which is the same as   
       ICHAR('A') on an EBCDIC machine. */

    inta = *(unsigned char *)ca;
    intb = *(unsigned char *)cb;

    if (zcode == 90 || zcode == 122) {

/*        ASCII is assumed - ZCODE is the ASCII code of either lower o
r   
          upper case 'Z'. */

	if (inta >= 97 && inta <= 122) {
	    inta += -32;
	}
	if (intb >= 97 && intb <= 122) {
	    intb += -32;
	}

    } else if (zcode == 233 || zcode == 169) {

/*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower
 or   
          upper case 'Z'. */

	if (inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta 
		>= 162 && inta <= 169) {
	    inta += 64;
	}
	if (intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb 
		>= 162 && intb <= 169) {
	    intb += 64;
	}

    } else if (zcode == 218 || zcode == 250) {

/*        ASCII is assumed, on Prime machines - ZCODE is the ASCII cod
e   
          plus 128 of either lower or upper case 'Z'. */

	if (inta >= 225 && inta <= 250) {
	    inta += -32;
	}
	if (intb >= 225 && intb <= 250) {
	    intb += -32;
	}
    }
    ret_val = inta == intb;

/*     RETURN   

       End of LSAME */

    return ret_val;
} /* lsame_ */


/* Subroutine */
#ifndef NOBLAS
static
#endif
int xerbla_(char *srname, integer *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    XERBLA  is an error handler for the LAPACK routines.   
    It is called by an LAPACK routine if an input parameter has an   
    invalid value.  A message is printed and execution stops.   

    Installers may consider modifying the STOP statement in order to   
    call system-specific exception-handling facilities.   

    Arguments   
    =========   

    SRNAME  (input) CHARACTER*6   
            The name of the routine which called XERBLA.   

    INFO    (input) INTEGER   
            The position of the invalid parameter in the parameter list   

            of the calling routine.   

   ===================================================================== 
*/

    printf("** On entry to %6s, parameter number %2i had an illegal value\n",
		srname, *info);

/*     End of XERBLA */

    return 0;
} /* xerbla_ */


#ifndef NOBLAS

/* Subroutine */
static  int dsyr2_(char *uplo, integer *n, double *alpha, 
	double *x, integer *incx, double *y, integer *incy, 
	double *a, integer *lda)
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static double temp1, temp2;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*  Purpose   
    =======   

    DSYR2  performs the symmetric rank 2 operation   

       A := alpha*x*y' + alpha*y*x' + A,   

    where alpha is a scalar, x and y are n element vectors and A is an n 
  
    by n symmetric matrix.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the array A is to be referenced as   
             follows:   

                UPLO = 'U' or 'u'   Only the upper triangular part of A   
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the lower triangular part of A   
                                    is to be referenced.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    X      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    Y      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y.   
             Unchanged on exit.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular part of the symmetric matrix and the strictly   
             lower triangular part of A is not referenced. On exit, the   
             upper triangular part of the array A is overwritten by the   
             upper triangular part of the updated matrix.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular part of the symmetric matrix and the strictly   
             upper triangular part of A is not referenced. On exit, the   
             lower triangular part of the array A is overwritten by the   
             lower triangular part of the updated matrix.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   


    Level 2 Blas routine.   

    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*n)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("DSYR2 ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0.) {
	return 0;
    }

/*     Set up the start points in X and Y if the increments are not both 
  
       unity. */

    if (*incx != 1 || *incy != 1) {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*n - 1) * *incx;
	}
	if (*incy > 0) {
	    ky = 1;
	} else {
	    ky = 1 - (*n - 1) * *incy;
	}
	jx = kx;
	jy = ky;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A. */

    if (lsame_(uplo, "U")) {

/*        Form  A  when A is stored in the upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(j) != 0. || Y(j) != 0.) {
		    temp1 = *alpha * Y(j);
		    temp2 = *alpha * X(j);
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			A(i,j) = A(i,j) + X(i) * temp1 
				+ Y(i) * temp2;
/* L10: */
		    }
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(jx) != 0. || Y(jy) != 0.) {
		    temp1 = *alpha * Y(jy);
		    temp2 = *alpha * X(jx);
		    ix = kx;
		    iy = ky;
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			A(i,j) = A(i,j) + X(ix) * temp1 
				+ Y(iy) * temp2;
			ix += *incx;
			iy += *incy;
/* L30: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when A is stored in the lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(j) != 0. || Y(j) != 0.) {
		    temp1 = *alpha * Y(j);
		    temp2 = *alpha * X(j);
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			A(i,j) = A(i,j) + X(i) * temp1 
				+ Y(i) * temp2;
/* L50: */
		    }
		}
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(jx) != 0. || Y(jy) != 0.) {
		    temp1 = *alpha * Y(jy);
		    temp2 = *alpha * X(jx);
		    ix = jx;
		    iy = jy;
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			A(i,j) = A(i,j) + X(ix) * temp1 
				+ Y(iy) * temp2;
			ix += *incx;
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    }

    return 0;

/*     End of DSYR2 . */

} /* dsyr2_ */


/* Subroutine */
static  int dsymv_(char *uplo, integer *n, double *alpha, 
	double *a, integer *lda, double *x, integer *incx, double 
	*beta, double *y, integer *incy)
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static double temp1, temp2;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*  Purpose   
    =======   

    DSYMV  performs the matrix-vector  operation   

       y := alpha*A*x + beta*y,   

    where alpha and beta are scalars, x and y are n element vectors and   
    A is an n by n symmetric matrix.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the array A is to be referenced as   
             follows:   

                UPLO = 'U' or 'u'   Only the upper triangular part of A   
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the lower triangular part of A   
                                    is to be referenced.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular part of the symmetric matrix and the strictly   
             lower triangular part of A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular part of the symmetric matrix and the strictly   
             upper triangular part of A is not referenced.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   

    X      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    BETA   - DOUBLE PRECISION.   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   

    Y      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y. On exit, Y is overwritten by the updated   
             vector y.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   


    Level 2 Blas routine.   

    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*lda < max(1,*n)) {
	info = 5;
    } else if (*incx == 0) {
	info = 7;
    } else if (*incy == 0) {
	info = 10;
    }
    if (info != 0) {
	xerbla_("DSYMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0. && *beta == 1.) {
	return 0;
    }

/*     Set up the start points in  X  and  Y. */

    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (*n - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A.   

       First form  y := beta*y. */

    if (*beta != 1.) {
	if (*incy == 1) {
	    if (*beta == 0.) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    Y(i) = 0.;
/* L10: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    Y(i) = *beta * Y(i);
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    Y(iy) = 0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    Y(iy) = *beta * Y(iy);
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == 0.) {
	return 0;
    }
    if (lsame_(uplo, "U")) {

/*        Form  y  when A is stored in upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(j);
		temp2 = 0.;
		i__2 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    Y(i) += temp1 * A(i,j);
		    temp2 += A(i,j) * X(i);
/* L50: */
		}
		Y(j) = Y(j) + temp1 * A(j,j) + *alpha * temp2;
/* L60: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(jx);
		temp2 = 0.;
		ix = kx;
		iy = ky;
		i__2 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    Y(iy) += temp1 * A(i,j);
		    temp2 += A(i,j) * X(ix);
		    ix += *incx;
		    iy += *incy;
/* L70: */
		}
		Y(jy) = Y(jy) + temp1 * A(j,j) + *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    } else {

/*        Form  y  when A is stored in lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(j);
		temp2 = 0.;
		Y(j) += temp1 * A(j,j);
		i__2 = *n;
		for (i = j + 1; i <= *n; ++i) {
		    Y(i) += temp1 * A(i,j);
		    temp2 += A(i,j) * X(i);
/* L90: */
		}
		Y(j) += *alpha * temp2;
/* L100: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(jx);
		temp2 = 0.;
		Y(jy) += temp1 * A(j,j);
		ix = jx;
		iy = jy;
		i__2 = *n;
		for (i = j + 1; i <= *n; ++i) {
		    ix += *incx;
		    iy += *incy;
		    Y(iy) += temp1 * A(i,j);
		    temp2 += A(i,j) * X(ix);
/* L110: */
		}
		Y(jy) += *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;

/*     End of DSYMV . */

} /* dsymv_ */


/* Subroutine */
static  int dsyr2k_(char *uplo, char *trans, integer *n, integer *k, 
	double *alpha, double *a, integer *lda, double *b, 
	integer *ldb, double *beta, double *c, integer *ldc)
{


    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer info;
    static double temp1, temp2;
    static integer i, j, l;
    extern logical lsame_(char *, char *);
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*  Purpose   
    =======   

    DSYR2K  performs one of the symmetric rank 2k operations   

       C := alpha*A*B' + alpha*B*A' + beta*C,   

    or   

       C := alpha*A'*B + alpha*B'*A + beta*C,   

    where  alpha and beta  are scalars, C is an  n by n  symmetric matrix 
  
    and  A and B  are  n by k  matrices  in the  first  case  and  k by n 
  
    matrices in the second case.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On  entry,   UPLO  specifies  whether  the  upper  or  lower 
  
             triangular  part  of the  array  C  is to be  referenced  as 
  
             follows:   

                UPLO = 'U' or 'u'   Only the  upper triangular part of  C 
  
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the  lower triangular part of  C 
  
                                    is to be referenced.   

             Unchanged on exit.   

    TRANS  - CHARACTER*1.   
             On entry,  TRANS  specifies the operation to be performed as 
  
             follows:   

                TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +   
                                          beta*C.   

                TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +   
                                          beta*C.   

                TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +   
                                          beta*C.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry,  N specifies the order of the matrix C.  N must be 
  
             at least zero.   
             Unchanged on exit.   

    K      - INTEGER.   
             On entry with  TRANS = 'N' or 'n',  K  specifies  the number 
  
             of  columns  of the  matrices  A and B,  and on  entry  with 
  
             TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number 
  
             of rows of the matrices  A and B.  K must be at least  zero. 
  
             Unchanged on exit.   

    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is 
  
             k  when  TRANS = 'N' or 'n',  and is  n  otherwise.   
             Before entry with  TRANS = 'N' or 'n',  the  leading  n by k 
  
             part of the array  A  must contain the matrix  A,  otherwise 
  
             the leading  k by n  part of the array  A  must contain  the 
  
             matrix A.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' 
  
             then  LDA must be at least  max( 1, n ), otherwise  LDA must 
  
             be at least  max( 1, k ).   
             Unchanged on exit.   

    B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is 
  
             k  when  TRANS = 'N' or 'n',  and is  n  otherwise.   
             Before entry with  TRANS = 'N' or 'n',  the  leading  n by k 
  
             part of the array  B  must contain the matrix  B,  otherwise 
  
             the leading  k by n  part of the array  B  must contain  the 
  
             matrix B.   
             Unchanged on exit.   

    LDB    - INTEGER.   
             On entry, LDB specifies the first dimension of B as declared 
  
             in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' 
  
             then  LDB must be at least  max( 1, n ), otherwise  LDB must 
  
             be at least  max( 1, k ).   
             Unchanged on exit.   

    BETA   - DOUBLE PRECISION.   
             On entry, BETA specifies the scalar beta.   
             Unchanged on exit.   

    C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).   
             Before entry  with  UPLO = 'U' or 'u',  the leading  n by n 
  
             upper triangular part of the array C must contain the upper 
  
             triangular part  of the  symmetric matrix  and the strictly 
  
             lower triangular part of C is not referenced.  On exit, the 
  
             upper triangular part of the array  C is overwritten by the 
  
             upper triangular part of the updated matrix.   
             Before entry  with  UPLO = 'L' or 'l',  the leading  n by n 
  
             lower triangular part of the array C must contain the lower 
  
             triangular part  of the  symmetric matrix  and the strictly 
  
             upper triangular part of C is not referenced.  On exit, the 
  
             lower triangular part of the array  C is overwritten by the 
  
             lower triangular part of the updated matrix.   

    LDC    - INTEGER.   
             On entry, LDC specifies the first dimension of C as declared 
  
             in  the  calling  (sub)  program.   LDC  must  be  at  least 
  
             max( 1, n ).   
             Unchanged on exit.   


    Level 3 Blas routine.   


    -- Written on 8-February-1989.   
       Jack Dongarra, Argonne National Laboratory.   
       Iain Duff, AERE Harwell.   
       Jeremy Du Croz, Numerical Algorithms Group Ltd.   
       Sven Hammarling, Numerical Algorithms Group Ltd.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#undef C
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    if (lsame_(trans, "N")) {
	nrowa = *n;
    } else {
	nrowa = *k;
    }
    upper = lsame_(uplo, "U");

    info = 0;
    if (! upper && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (! lsame_(trans, "N") && ! lsame_(trans, "T") &&
	     ! lsame_(trans, "C")) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*k < 0) {
	info = 4;
    } else if (*lda < max(1,nrowa)) {
	info = 7;
    } else if (*ldb < max(1,nrowa)) {
	info = 9;
    } else if (*ldc < max(1,*n)) {
	info = 12;
    }
    if (info != 0) {
	xerbla_("DSYR2K", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == 0.) {
	if (upper) {
	    if (*beta == 0.) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			C(i,j) = 0.;
/* L10: */
		    }
/* L20: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			C(i,j) = *beta * C(i,j);
/* L30: */
		    }
/* L40: */
		}
	    }
	} else {
	    if (*beta == 0.) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			C(i,j) = 0.;
/* L50: */
		    }
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			C(i,j) = *beta * C(i,j);
/* L70: */
		    }
/* L80: */
		}
	    }
	}
	return 0;
    }

/*     Start the operations. */

    if (lsame_(trans, "N")) {

/*        Form  C := alpha*A*B' + alpha*B*A' + C. */

	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (*beta == 0.) {
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			C(i,j) = 0.;
/* L90: */
		    }
		} else if (*beta != 1.) {
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			C(i,j) = *beta * C(i,j);
/* L100: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= *k; ++l) {
		    if (A(j,l) != 0. || B(j,l) != 0.) {
			temp1 = *alpha * B(j,l);
			temp2 = *alpha * A(j,l);
			i__3 = j;
			for (i = 1; i <= j; ++i) {
			    C(i,j) = C(i,j) + A(i,l) * temp1 + B(i,l) * 
				    temp2;
/* L110: */
			}
		    }
/* L120: */
		}
/* L130: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (*beta == 0.) {
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			C(i,j) = 0.;
/* L140: */
		    }
		} else if (*beta != 1.) {
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			C(i,j) = *beta * C(i,j);
/* L150: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= *k; ++l) {
		    if (A(j,l) != 0. || B(j,l) != 0.) {
			temp1 = *alpha * B(j,l);
			temp2 = *alpha * A(j,l);
			i__3 = *n;
			for (i = j; i <= *n; ++i) {
			    C(i,j) = C(i,j) + A(i,l) * temp1 + B(i,l) * 
				    temp2;
/* L160: */
			}
		    }
/* L170: */
		}
/* L180: */
	    }
	}
    } else {

/*        Form  C := alpha*A'*B + alpha*B'*A + C. */

	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j;
		for (i = 1; i <= j; ++i) {
		    temp1 = 0.;
		    temp2 = 0.;
		    i__3 = *k;
		    for (l = 1; l <= *k; ++l) {
			temp1 += A(l,i) * B(l,j);
			temp2 += B(l,i) * A(l,j);
/* L190: */
		    }
		    if (*beta == 0.) {
			C(i,j) = *alpha * temp1 + *alpha * temp2;
		    } else {
			C(i,j) = *beta * C(i,j) + *
				alpha * temp1 + *alpha * temp2;
		    }
/* L200: */
		}
/* L210: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = *n;
		for (i = j; i <= *n; ++i) {
		    temp1 = 0.;
		    temp2 = 0.;
		    i__3 = *k;
		    for (l = 1; l <= *k; ++l) {
			temp1 += A(l,i) * B(l,j);
			temp2 += B(l,i) * A(l,j);
/* L220: */
		    }
		    if (*beta == 0.) {
			C(i,j) = *alpha * temp1 + *alpha * temp2;
		    } else {
			C(i,j) = *beta * C(i,j) + *
				alpha * temp1 + *alpha * temp2;
		    }
/* L230: */
		}
/* L240: */
	    }
	}
    }

    return 0;

/*     End of DSYR2K. */

} /* dsyr2k_ */


/* Subroutine */
static  int dtrmm_(char *side, char *uplo, char *transa, char *diag,
			    integer *m, integer *n,
			    double *alpha, double *a,
			    integer *lda, double *b, integer *ldb)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer info;
    static double temp;
    static integer i, j, k;
    static logical lside;
    extern logical blas_lsame_();
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int xerbla_();
    static logical nounit;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DTRMM  performs one of the matrix-matrix operations */

/*     B := alpha*op( A )*B,   or   B := alpha*B*op( A ), */

/*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or 
*/
/*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of 
*/

/*     op( A ) = A   or   op( A ) = A'. */

/*  Parameters */
/*  ========== */

/*  SIDE   - CHARACTER*1. */
/*           On entry,  SIDE specifies whether  op( A ) multiplies B from 
*/
/*           the left or right as follows: */

/*              SIDE = 'L' or 'l'   B := alpha*op( A )*B. */

/*              SIDE = 'R' or 'r'   B := alpha*B*op( A ). */

/*           Unchanged on exit. */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix A is an upper or 
*/
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANSA - CHARACTER*1. */
/*           On entry, TRANSA specifies the form of op( A ) to be used in 
*/
/*           the matrix multiplication as follows: */

/*              TRANSA = 'N' or 'n'   op( A ) = A. */

/*              TRANSA = 'T' or 't'   op( A ) = A'. */

/*              TRANSA = 'C' or 'c'   op( A ) = A'. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit triangular 
*/
/*           as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  M      - INTEGER. */
/*           On entry, M specifies the number of rows of B. M must be at 
*/
/*           least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the number of columns of B.  N must be 
*/
/*           at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is 
*/
/*           zero then  A is not referenced and  B need not be set before 
*/
/*           entry. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m 
*/
/*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'. 
*/
/*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k 
*/
/*           upper triangular part of the array  A must contain the upper 
*/
/*           triangular matrix  and the strictly lower triangular part of 
*/
/*           A is not referenced. */
/*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k 
*/
/*           lower triangular part of the array  A must contain the lower 
*/
/*           triangular matrix  and the strictly upper triangular part of 
*/
/*           A is not referenced. */
/*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of 
*/
/*           A  are not referenced either,  but are assumed to be  unity. 
*/
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared 
*/
/*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then 
*/
/*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r' 
*/
/*           then LDA must be at least max( 1, n ). */
/*           Unchanged on exit. */

/*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ). */
/*           Before entry,  the leading  m by n part of the array  B must 
*/
/*           contain the matrix  B,  and  on exit  is overwritten  by the 
*/
/*           transformed matrix. */

/*  LDB    - INTEGER. */
/*           On entry, LDB specifies the first dimension of B as declared 
*/
/*           in  the  calling  (sub)  program.   LDB  must  be  at  least 
*/
/*           max( 1, m ). */
/*           Unchanged on exit. */


/*  Level 3 Blas routine. */

/*  -- Written on 8-February-1989. */
/*     Jack Dongarra, Argonne National Laboratory. */
/*     Iain Duff, AERE Harwell. */
/*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/*     Sven Hammarling, Numerical Algorithms Group Ltd. */


/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Parameters .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    b_dim1 = *ldb;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    lside = blas_lsame_(side, "L", 1L, 1L);
    if (lside) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    nounit = blas_lsame_(diag, "N", 1L, 1L);
    upper = blas_lsame_(uplo, "U", 1L, 1L);

    info = 0;
    if (! lside && ! blas_lsame_(side, "R", 1L, 1L)) {
	info = 1;
    } else if (! upper && ! blas_lsame_(uplo, "L", 1L, 1L)) {
	info = 2;
    } else if (! blas_lsame_(transa, "N", 1L, 1L) && ! blas_lsame_(transa, "T", 1L, 1L) 
	    && ! blas_lsame_(transa, "C", 1L, 1L)) {
	info = 3;
    } else if (! blas_lsame_(diag, "U", 1L, 1L) && ! blas_lsame_(diag, "N", 1L, 1L)) {
	info = 4;
    } else if (*m < 0) {
	info = 5;
    } else if (*n < 0) {
	info = 6;
    } else if (*lda < max(1,nrowa)) {
	info = 9;
    } else if (*ldb < max(1,*m)) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("DTRMM ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == 0.) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i = 1; i <= i__2; ++i) {
		b[i + j * b_dim1] = 0.;
/* L10: */
	    }
/* L20: */
	}
	return 0;
    }

/*     Start the operations. */

    if (lside) {
	if (blas_lsame_(transa, "N", 1L, 1L)) {

/*           Form  B := alpha*A*B. */

	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (k = 1; k <= i__2; ++k) {
			if (b[k + j * b_dim1] != 0.) {
			    temp = *alpha * b[k + j * b_dim1];
			    i__3 = k - 1;
			    for (i = 1; i <= i__3; ++i) {
				b[i + j * b_dim1] += temp * a[i + k * a_dim1];
/* L30: */
			    }
			    if (nounit) {
				temp *= a[k + k * a_dim1];
			    }
			    b[k + j * b_dim1] = temp;
			}
/* L40: */
		    }
/* L50: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    for (k = *m; k >= 1; --k) {
			if (b[k + j * b_dim1] != 0.) {
			    temp = *alpha * b[k + j * b_dim1];
			    b[k + j * b_dim1] = temp;
			    if (nounit) {
				b[k + j * b_dim1] *= a[k + k * a_dim1];
			    }
			    i__2 = *m;
			    for (i = k + 1; i <= i__2; ++i) {
				b[i + j * b_dim1] += temp * a[i + k * a_dim1];
/* L60: */
			    }
			}
/* L70: */
		    }
/* L80: */
		}
	    }
	} else {

/*           Form  B := alpha*B*A'. */

	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    for (i = *m; i >= 1; --i) {
			temp = b[i + j * b_dim1];
			if (nounit) {
			    temp *= a[i + i * a_dim1];
			}
			i__2 = i - 1;
			for (k = 1; k <= i__2; ++k) {
			    temp += a[k + i * a_dim1] * b[k + j * b_dim1];
/* L90: */
			}
			b[i + j * b_dim1] = *alpha * temp;
/* L100: */
		    }
/* L110: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i = 1; i <= i__2; ++i) {
			temp = b[i + j * b_dim1];
			if (nounit) {
			    temp *= a[i + i * a_dim1];
			}
			i__3 = *m;
			for (k = i + 1; k <= i__3; ++k) {
			    temp += a[k + i * a_dim1] * b[k + j * b_dim1];
/* L120: */
			}
			b[i + j * b_dim1] = *alpha * temp;
/* L130: */
		    }
/* L140: */
		}
	    }
	}
    } else {
	if (blas_lsame_(transa, "N", 1L, 1L)) {

/*           Form  B := alpha*B*A. */

	    if (upper) {
		for (j = *n; j >= 1; --j) {
		    temp = *alpha;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__1 = *m;
		    for (i = 1; i <= i__1; ++i) {
			b[i + j * b_dim1] = temp * b[i + j * b_dim1];
/* L150: */
		    }
		    i__1 = j - 1;
		    for (k = 1; k <= i__1; ++k) {
			if (a[k + j * a_dim1] != 0.) {
			    temp = *alpha * a[k + j * a_dim1];
			    i__2 = *m;
			    for (i = 1; i <= i__2; ++i) {
				b[i + j * b_dim1] += temp * b[i + k * b_dim1];
/* L160: */
			    }
			}
/* L170: */
		    }
/* L180: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = *alpha;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = *m;
		    for (i = 1; i <= i__2; ++i) {
			b[i + j * b_dim1] = temp * b[i + j * b_dim1];
/* L190: */
		    }
		    i__2 = *n;
		    for (k = j + 1; k <= i__2; ++k) {
			if (a[k + j * a_dim1] != 0.) {
			    temp = *alpha * a[k + j * a_dim1];
			    i__3 = *m;
			    for (i = 1; i <= i__3; ++i) {
				b[i + j * b_dim1] += temp * b[i + k * b_dim1];
/* L200: */
			    }
			}
/* L210: */
		    }
/* L220: */
		}
	    }
	} else {

/*           Form  B := alpha*B*A'. */

	    if (upper) {
		i__1 = *n;
		for (k = 1; k <= i__1; ++k) {
		    i__2 = k - 1;
		    for (j = 1; j <= i__2; ++j) {
			if (a[j + k * a_dim1] != 0.) {
			    temp = *alpha * a[j + k * a_dim1];
			    i__3 = *m;
			    for (i = 1; i <= i__3; ++i) {
				b[i + j * b_dim1] += temp * b[i + k * b_dim1];
/* L230: */
			    }
			}
/* L240: */
		    }
		    temp = *alpha;
		    if (nounit) {
			temp *= a[k + k * a_dim1];
		    }
		    if (temp != 1.) {
			i__2 = *m;
			for (i = 1; i <= i__2; ++i) {
			    b[i + k * b_dim1] = temp * b[i + k * b_dim1];
/* L250: */
			}
		    }
/* L260: */
		}
	    } else {
		for (k = *n; k >= 1; --k) {
		    i__1 = *n;
		    for (j = k + 1; j <= i__1; ++j) {
			if (a[j + k * a_dim1] != 0.) {
			    temp = *alpha * a[j + k * a_dim1];
			    i__2 = *m;
			    for (i = 1; i <= i__2; ++i) {
				b[i + j * b_dim1] += temp * b[i + k * b_dim1];
/* L270: */
			    }
			}
/* L280: */
		    }
		    temp = *alpha;
		    if (nounit) {
			temp *= a[k + k * a_dim1];
		    }
		    if (temp != 1.) {
			i__1 = *m;
			for (i = 1; i <= i__1; ++i) {
			    b[i + k * b_dim1] = temp * b[i + k * b_dim1];
/* L290: */
			}
		    }
/* L300: */
		}
	    }
	}
    }

    return 0;

/*     End of DTRMM . */

} /* dtrmm_ */


/* Subroutine */
static  int dscal_(n, da, dx, incx)
integer *n;
double *da, *dx;
integer *incx;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, mp1;


/*     scales a vector by a constant. */
/*     uses unrolled loops for increment equal to one. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified to correct problem with negative increment, 8/21/90. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dx[ix] = *da * dx[ix];
	ix += *incx;
/* L10: */
    }
    return 0;

/*        code for increment equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	dx[i] = *da * dx[i];
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 5) {
	dx[i] = *da * dx[i];
	dx[i + 1] = *da * dx[i + 1];
	dx[i + 2] = *da * dx[i + 2];
	dx[i + 3] = *da * dx[i + 3];
	dx[i + 4] = *da * dx[i + 4];
/* L50: */
    }
    return 0;
} /* dscal_ */


/* Subroutine */
static  int dgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
char *trans;
integer *m, *n;
double *alpha, *a;
integer *lda;
double *x;
integer *incx;
double *beta, *y;
integer *incy;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static double temp;
    static integer lenx, leny, i, j;
    extern logical blas_lsame_();
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_();

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGEMV  performs one of the matrix-vector operations */

/*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y, */

/*  where alpha and beta are scalars, x and y are vectors and A is an */
/*  m by n matrix. */

/*  Parameters */
/*  ========== */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the operation to be performed as */
/*           follows: */

/*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y. */

/*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y. */

/*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y. */

/*           Unchanged on exit. */

/*  M      - INTEGER. */
/*           On entry, M specifies the number of rows of the matrix A. */
/*           M must be at least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the number of columns of the matrix A. 
*/
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/*           Before entry, the leading m by n part of the array A must */
/*           contain the matrix of coefficients. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared 
*/
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, m ). */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of DIMENSION at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' */
/*           and at least */
/*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. */
/*           Before entry, the incremented array X must contain the */
/*           vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  BETA   - DOUBLE PRECISION. */
/*           On entry, BETA specifies the scalar beta. When BETA is */
/*           supplied as zero then Y need not be set on input. */
/*           Unchanged on exit. */

/*  Y      - DOUBLE PRECISION array of DIMENSION at least */
/*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' */
/*           and at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. */
/*           Before entry with BETA non-zero, the incremented array Y */
/*           must contain the vector y. On exit, Y is overwritten by the 
*/
/*           updated vector y. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --y;
    --x;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    info = 0;
    if (! blas_lsame_(trans, "N", 1L, 1L) && ! blas_lsame_(trans, "T", 1L, 1L) && ! 
	    blas_lsame_(trans, "C", 1L, 1L)) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*lda < max(1,*m)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("DGEMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
	return 0;
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set 
*/
/*     up the start points in  X  and  Y. */

    if (blas_lsame_(trans, "N", 1L, 1L)) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

    if (*beta != 1.) {
	if (*incy == 1) {
	    if (*beta == 0.) {
		i__1 = leny;
		for (i = 1; i <= i__1; ++i) {
		    y[i] = 0.;
/* L10: */
		}
	    } else {
		i__1 = leny;
		for (i = 1; i <= i__1; ++i) {
		    y[i] = *beta * y[i];
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.) {
		i__1 = leny;
		for (i = 1; i <= i__1; ++i) {
		    y[iy] = 0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = leny;
		for (i = 1; i <= i__1; ++i) {
		    y[iy] = *beta * y[iy];
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == 0.) {
	return 0;
    }
    if (blas_lsame_(trans, "N", 1L, 1L)) {

/*        Form  y := alpha*A*x + y. */

	jx = kx;
	if (*incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.) {
		    temp = *alpha * x[jx];
		    i__2 = *m;
		    for (i = 1; i <= i__2; ++i) {
			y[i] += temp * a[i + j * a_dim1];
/* L50: */
		    }
		}
		jx += *incx;
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.) {
		    temp = *alpha * x[jx];
		    iy = ky;
		    i__2 = *m;
		    for (i = 1; i <= i__2; ++i) {
			y[iy] += temp * a[i + j * a_dim1];
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
/* L80: */
	    }
	}
    } else {

/*        Form  y := alpha*A'*x + y. */

	jy = ky;
	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp = 0.;
		i__2 = *m;
		for (i = 1; i <= i__2; ++i) {
		    temp += a[i + j * a_dim1] * x[i];
/* L90: */
		}
		y[jy] += *alpha * temp;
		jy += *incy;
/* L100: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp = 0.;
		ix = kx;
		i__2 = *m;
		for (i = 1; i <= i__2; ++i) {
		    temp += a[i + j * a_dim1] * x[ix];
		    ix += *incx;
/* L110: */
		}
		y[jy] += *alpha * temp;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;

/*     End of DGEMV . */

} /* dgemv_ */


/* Subroutine */
static  int dswap_(n, dx, incx, dy, incy)
integer *n;
double *dx;
integer *incx;
double *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m;
    static double dtemp;
    static integer ix, iy, mp1;


/*     interchanges two vectors. */
/*     uses unrolled loops for increments equal one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dtemp = dx[ix];
	dx[ix] = dy[iy];
	dy[iy] = dtemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       code for both increments equal to 1 */


/*       clean-up loop */

L20:
    m = *n % 3;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	dtemp = dx[i];
	dx[i] = dy[i];
	dy[i] = dtemp;
/* L30: */
    }
    if (*n < 3) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 3) {
	dtemp = dx[i];
	dx[i] = dy[i];
	dy[i] = dtemp;
	dtemp = dx[i + 1];
	dx[i + 1] = dy[i + 1];
	dy[i + 1] = dtemp;
	dtemp = dx[i + 2];
	dx[i + 2] = dy[i + 2];
	dy[i + 2] = dtemp;
/* L50: */
    }
    return 0;
} /* dswap_ */


/* Subroutine */
static  int dcopy_(n, dx, incx, dy, incy)
integer *n;
double *dx;
integer *incx;
double *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     copies a vector, x, to a vector, y. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dy[iy] = dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	dy[i] = dx[i];
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 7) {
	dy[i] = dx[i];
	dy[i + 1] = dx[i + 1];
	dy[i + 2] = dx[i + 2];
	dy[i + 3] = dx[i + 3];
	dy[i + 4] = dx[i + 4];
	dy[i + 5] = dx[i + 5];
	dy[i + 6] = dx[i + 6];
/* L50: */
    }
    return 0;
} /* dcopy_ */


static double ddot_(n, dx, incx, dy, incy)
integer *n;
double *dx;
integer *incx;
double *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1;
    double ret_val;

    /* Local variables */
    static integer i, m;
    static double dtemp;
    static integer ix, iy, mp1;


/*     forms the dot product of two vectors. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dtemp += dx[ix] * dy[iy];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	dtemp += dx[i] * dy[i];
/* L30: */
    }
    if (*n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 5) {
	dtemp = dtemp + dx[i] * dy[i] + dx[i + 1] * dy[i + 1] + dx[i + 2] * 
		dy[i + 2] + dx[i + 3] * dy[i + 3] + dx[i + 4] * dy[i + 4];
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* ddot_ */


/* Subroutine */
static  int dgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, 
	beta, c, ldc)
char *transa, *transb;
integer *m, *n, *k;
double *alpha, *a;
integer *lda;
double *b;
integer *ldb;
double *beta, *c;
integer *ldc;
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer info;
    static logical nota, notb;
    static double temp;
    static integer i, j, l, ncola;
    extern logical blas_lsame_();
    static integer nrowa, nrowb;
    extern /* Subroutine */ int xerbla_();

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGEMM  performs one of the matrix-matrix operations */

/*     C := alpha*op( A )*op( B ) + beta*C, */

/*  where  op( X ) is one of */

/*     op( X ) = X   or   op( X ) = X', */

/*  alpha and beta are scalars, and A, B and C are matrices, with op( A ) 
*/
/*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix. 
*/

/*  Parameters */
/*  ========== */

/*  TRANSA - CHARACTER*1. */
/*           On entry, TRANSA specifies the form of op( A ) to be used in 
*/
/*           the matrix multiplication as follows: */

/*              TRANSA = 'N' or 'n',  op( A ) = A. */

/*              TRANSA = 'T' or 't',  op( A ) = A'. */

/*              TRANSA = 'C' or 'c',  op( A ) = A'. */

/*           Unchanged on exit. */

/*  TRANSB - CHARACTER*1. */
/*           On entry, TRANSB specifies the form of op( B ) to be used in 
*/
/*           the matrix multiplication as follows: */

/*              TRANSB = 'N' or 'n',  op( B ) = B. */

/*              TRANSB = 'T' or 't',  op( B ) = B'. */

/*              TRANSB = 'C' or 'c',  op( B ) = B'. */

/*           Unchanged on exit. */

/*  M      - INTEGER. */
/*           On entry,  M  specifies  the number  of rows  of the  matrix 
*/
/*           op( A )  and of the  matrix  C.  M  must  be at least  zero. 
*/
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry,  N  specifies the number  of columns of the matrix 
*/
/*           op( B ) and the number of columns of the matrix C. N must be 
*/
/*           at least zero. */
/*           Unchanged on exit. */

/*  K      - INTEGER. */
/*           On entry,  K  specifies  the number of columns of the matrix 
*/
/*           op( A ) and the number of rows of the matrix op( B ). K must 
*/
/*           be at least  zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is 
*/
/*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise. */
/*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k 
*/
/*           part of the array  A  must contain the matrix  A,  otherwise 
*/
/*           the leading  k by m  part of the array  A  must contain  the 
*/
/*           matrix A. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared 
*/
/*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then 
*/
/*           LDA must be at least  max( 1, m ), otherwise  LDA must be at 
*/
/*           least  max( 1, k ). */
/*           Unchanged on exit. */

/*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is 
*/
/*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise. */
/*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n 
*/
/*           part of the array  B  must contain the matrix  B,  otherwise 
*/
/*           the leading  n by k  part of the array  B  must contain  the 
*/
/*           matrix B. */
/*           Unchanged on exit. */

/*  LDB    - INTEGER. */
/*           On entry, LDB specifies the first dimension of B as declared 
*/
/*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then 
*/
/*           LDB must be at least  max( 1, k ), otherwise  LDB must be at 
*/
/*           least  max( 1, n ). */
/*           Unchanged on exit. */

/*  BETA   - DOUBLE PRECISION. */
/*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is 
*/
/*           supplied as zero then C need not be set on input. */
/*           Unchanged on exit. */

/*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ). */
/*           Before entry, the leading  m by n  part of the array  C must 
*/
/*           contain the matrix  C,  except when  beta  is zero, in which 
*/
/*           case C need not be set on entry. */
/*           On exit, the array  C  is overwritten by the  m by n  matrix 
*/
/*           ( alpha*op( A )*op( B ) + beta*C ). */

/*  LDC    - INTEGER. */
/*           On entry, LDC specifies the first dimension of C as declared 
*/
/*           in  the  calling  (sub)  program.   LDC  must  be  at  least 
*/
/*           max( 1, m ). */
/*           Unchanged on exit. */


/*  Level 3 Blas routine. */

/*  -- Written on 8-February-1989. */
/*     Jack Dongarra, Argonne National Laboratory. */
/*     Iain Duff, AERE Harwell. */
/*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/*     Sven Hammarling, Numerical Algorithms Group Ltd. */


/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Parameters .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not 
*/
/*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows 
*/
/*     and  columns of  A  and the  number of  rows  of  B  respectively. 
*/

    /* Parameter adjustments */
    c_dim1 = *ldc;
    c_offset = c_dim1 + 1;
    c -= c_offset;
    b_dim1 = *ldb;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    nota = blas_lsame_(transa, "N", 1L, 1L);
    notb = blas_lsame_(transb, "N", 1L, 1L);
    if (nota) {
	nrowa = *m;
	ncola = *k;
    } else {
	nrowa = *k;
	ncola = *m;
    }
    if (notb) {
	nrowb = *k;
    } else {
	nrowb = *n;
    }

/*     Test the input parameters. */

    info = 0;
    if (! nota && ! blas_lsame_(transa, "C", 1L, 1L) && ! blas_lsame_(transa, "T", 1L, 
	    1L)) {
	info = 1;
    } else if (! notb && ! blas_lsame_(transb, "C", 1L, 1L) && ! blas_lsame_(transb, 
	    "T", 1L, 1L)) {
	info = 2;
    } else if (*m < 0) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*k < 0) {
	info = 5;
    } else if (*lda < max(1,nrowa)) {
	info = 8;
    } else if (*ldb < max(1,nrowb)) {
	info = 10;
    } else if (*ldc < max(1,*m)) {
	info = 13;
    }
    if (info != 0) {
	xerbla_("DGEMM ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
	return 0;
    }

/*     And if  alpha.eq.zero. */

    if (*alpha == 0.) {
	if (*beta == 0.) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i = 1; i <= i__2; ++i) {
		    c[i + j * c_dim1] = 0.;
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i = 1; i <= i__2; ++i) {
		    c[i + j * c_dim1] = *beta * c[i + j * c_dim1];
/* L30: */
		}
/* L40: */
	    }
	}
	return 0;
    }

/*     Start the operations. */

    if (notb) {
	if (nota) {

/*           Form  C := alpha*A*B + beta*C. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.) {
		    i__2 = *m;
		    for (i = 1; i <= i__2; ++i) {
			c[i + j * c_dim1] = 0.;
/* L50: */
		    }
		} else if (*beta != 1.) {
		    i__2 = *m;
		    for (i = 1; i <= i__2; ++i) {
			c[i + j * c_dim1] = *beta * c[i + j * c_dim1];
/* L60: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (b[l + j * b_dim1] != 0.) {
			temp = *alpha * b[l + j * b_dim1];
			i__3 = *m;
			for (i = 1; i <= i__3; ++i) {
			    c[i + j * c_dim1] += temp * a[i + l * a_dim1];
/* L70: */
			}
		    }
/* L80: */
		}
/* L90: */
	    }
	} else {

/*           Form  C := alpha*A'*B + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i = 1; i <= i__2; ++i) {
		    temp = 0.;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp += a[l + i * a_dim1] * b[l + j * b_dim1];
/* L100: */
		    }
		    if (*beta == 0.) {
			c[i + j * c_dim1] = *alpha * temp;
		    } else {
			c[i + j * c_dim1] = *alpha * temp + *beta * c[i + j * 
				c_dim1];
		    }
/* L110: */
		}
/* L120: */
	    }
	}
    } else {
	if (nota) {

/*           Form  C := alpha*A*B' + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.) {
		    i__2 = *m;
		    for (i = 1; i <= i__2; ++i) {
			c[i + j * c_dim1] = 0.;
/* L130: */
		    }
		} else if (*beta != 1.) {
		    i__2 = *m;
		    for (i = 1; i <= i__2; ++i) {
			c[i + j * c_dim1] = *beta * c[i + j * c_dim1];
/* L140: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (b[j + l * b_dim1] != 0.) {
			temp = *alpha * b[j + l * b_dim1];
			i__3 = *m;
			for (i = 1; i <= i__3; ++i) {
			    c[i + j * c_dim1] += temp * a[i + l * a_dim1];
/* L150: */
			}
		    }
/* L160: */
		}
/* L170: */
	    }
	} else {

/*           Form  C := alpha*A'*B' + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i = 1; i <= i__2; ++i) {
		    temp = 0.;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp += a[l + i * a_dim1] * b[j + l * b_dim1];
/* L180: */
		    }
		    if (*beta == 0.) {
			c[i + j * c_dim1] = *alpha * temp;
		    } else {
			c[i + j * c_dim1] = *alpha * temp + *beta * c[i + j * 
				c_dim1];
		    }
/* L190: */
		}
/* L200: */
	    }
	}
    }

    return 0;

/*     End of DGEMM . */

} /* dgemm_ */


/* Subroutine */
static  int dger_(m, n, alpha, x, incx, y, incy, a, lda)
integer *m, *n;
double *alpha, *x;
integer *incx;
double *y;
integer *incy;
double *a;
integer *lda;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static double temp;
    static integer i, j, ix, jy, kx;
    extern /* Subroutine */ int xerbla_();

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGER   performs the rank 1 operation */

/*     A := alpha*x*y' + A, */

/*  where alpha is a scalar, x is an m element vector, y is an n element 
*/
/*  vector and A is an m by n matrix. */

/*  Parameters */
/*  ========== */

/*  M      - INTEGER. */
/*           On entry, M specifies the number of rows of the matrix A. */
/*           M must be at least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the number of columns of the matrix A. 
*/
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( m - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the m */
/*           element vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  Y      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ). */
/*           Before entry, the incremented array Y must contain the n */
/*           element vector y. */
/*           Unchanged on exit. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/*           Before entry, the leading m by n part of the array A must */
/*           contain the matrix of coefficients. On exit, A is */
/*           overwritten by the updated matrix. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared 
*/
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, m ). */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --y;
    --x;

    /* Function Body */
    info = 0;
    if (*m < 0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*m)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("DGER  ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == 0.) {
	return 0;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

    if (*incy > 0) {
	jy = 1;
    } else {
	jy = 1 - (*n - 1) * *incy;
    }
    if (*incx == 1) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (y[jy] != 0.) {
		temp = *alpha * y[jy];
		i__2 = *m;
		for (i = 1; i <= i__2; ++i) {
		    a[i + j * a_dim1] += x[i] * temp;
/* L10: */
		}
	    }
	    jy += *incy;
/* L20: */
	}
    } else {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*m - 1) * *incx;
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (y[jy] != 0.) {
		temp = *alpha * y[jy];
		ix = kx;
		i__2 = *m;
		for (i = 1; i <= i__2; ++i) {
		    a[i + j * a_dim1] += x[ix] * temp;
		    ix += *incx;
/* L30: */
		}
	    }
	    jy += *incy;
/* L40: */
	}
    }

    return 0;

/*     End of DGER  . */

} /* dger_ */


/* Subroutine */
static  int daxpy_(n, da, dx, incx, dy, incy)
integer *n;
double *da, *dx;
integer *incx;
double *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     constant times a vector plus a vector. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*da == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dy[iy] += *da * dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
	dy[i] += *da * dx[i];
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 4) {
	dy[i] += *da * dx[i];
	dy[i + 1] += *da * dx[i + 1];
	dy[i + 2] += *da * dx[i + 2];
	dy[i + 3] += *da * dx[i + 3];
/* L50: */
    }
    return 0;
} /* daxpy_ */


static double dnrm2_(n, dx, incx)
integer *n;
double *dx;
integer *incx;
{
    /* Initialized data */

    static double zero = 0.;
    static double one = 1.;
    static double cutlo = 8.232e-11;
    static double cuthi = 1.304e19;

    /* Format strings */
    static char fmt_30[] = "";
    static char fmt_50[] = "";
    static char fmt_70[] = "";
    static char fmt_110[] = "";

    /* System generated locals */
    integer i__1;
    double ret_val, d__1;

    /* Local variables */
    static double xmax;
    static integer next, i, j, ix;
    static double hitest, sum;

    /* Assigned format variables */
    char *next_fmt;

    /* Parameter adjustments */
    --dx;

    /* Function Body */

/*     euclidean norm of the n-vector stored in dx() with storage */
/*     increment incx . */
/*     if    n .le. 0 return with result = 0. */
/*     if n .ge. 1 then incx must be .ge. 1 */

/*           c.l.lawson, 1978 jan 08 */
/*     modified to correct problem with negative increment, 8/21/90. */
/*     modified to correct failure to update ix, 1/25/92. */

/*     four phase method     using two built-in constants that are */
/*     hopefully applicable to all machines. */
/*         cutlo = maximum of  dsqrt(u/eps)  over all known machines. */
/*         cuthi = minimum of  dsqrt(v)      over all known machines. */
/*     where */
/*         eps = smallest no. such that eps + 1. .gt. 1. */
/*         u   = smallest positive no.   (underflow limit) */
/*         v   = largest  no.            (overflow  limit) */

/*     brief outline of algorithm.. */

/*     phase 1    scans zero components. */
/*     move to phase 2 when a component is nonzero and .le. cutlo */
/*     move to phase 3 when a component is .gt. cutlo */
/*     move to phase 4 when a component is .ge. cuthi/m */
/*     where m = n for x() float and m = 2*n for complex. */

/*     values for cutlo and cuthi.. */
/*     from the environmental parameters listed in the imsl converter */
/*     document the limiting values are as follows.. */
/*     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are 
*/
/*                   univac and dec at 2**(-103) */
/*                   thus cutlo = 2**(-51) = 4.44089e-16 */
/*     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec. */
/*                   thus cuthi = 2**(63.5) = 1.30438e19 */
/*     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec. */
/*                   thus cutlo = 2**(-33.5) = 8.23181d-11 */
/*     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19 */
/*     data cutlo, cuthi / 8.232d-11,  1.304d19 / */
/*     data cutlo, cuthi / 4.441e-16,  1.304e19 / */

    if (*n > 0) {
	goto L10;
    }
    ret_val = zero;
    goto L300;

L10:
    next = 0;
    next_fmt = fmt_30;
    sum = zero;
    i = 1;
    if (*incx < 0) {
	i = (-(*n) + 1) * *incx + 1;
    }
    ix = 1;
/*                                                 begin main loop */
L20:
    switch ((int)next) {
	case 0: goto L30;
	case 1: goto L50;
	case 2: goto L70;
	case 3: goto L110;
    }
L30:
    if ((d__1 = dx[i], abs(d__1)) > cutlo) {
	goto L85;
    }
    next = 1;
    next_fmt = fmt_50;
    xmax = zero;

/*                        phase 1.  sum is zero */

L50:
    if (dx[i] == zero) {
	goto L200;
    }
    if ((d__1 = dx[i], abs(d__1)) > cutlo) {
	goto L85;
    }

/*                                prepare for phase 2. */
    next = 2;
    next_fmt = fmt_70;
    goto L105;

/*                                prepare for phase 4. */

L100:
    ix = j;
    next = 3;
    next_fmt = fmt_110;
    sum = sum / dx[i] / dx[i];
L105:
    xmax = (d__1 = dx[i], abs(d__1));
    goto L115;

/*                   phase 2.  sum is small. */
/*                             scale to avoid destructive underflow. */

L70:
    if ((d__1 = dx[i], abs(d__1)) > cutlo) {
	goto L75;
    }

/*                     common code for phases 2 and 4. */
/*                     in phase 4 sum is large.  scale to avoid overflow. 
*/

L110:
    if ((d__1 = dx[i], abs(d__1)) <= xmax) {
	goto L115;
    }
/* Computing 2nd power */
    d__1 = xmax / dx[i];
    sum = one + sum * (d__1 * d__1);
    xmax = (d__1 = dx[i], abs(d__1));
    goto L200;

L115:
/* Computing 2nd power */
    d__1 = dx[i] / xmax;
    sum += d__1 * d__1;
    goto L200;


/*                  prepare for phase 3. */

L75:
    sum = sum * xmax * xmax;


/*     for float or d.p. set hitest = cuthi/n */
/*     for complex      set hitest = cuthi/(2*n) */

L85:
    hitest = cuthi / (float) (*n);

/*                   phase 3.  sum is mid-range.  no scaling. */

    i__1 = *n;
    for (j = ix; j <= i__1; ++j) {
	if ((d__1 = dx[i], abs(d__1)) >= hitest) {
	    goto L100;
	}
/* Computing 2nd power */
	d__1 = dx[i];
	sum += d__1 * d__1;
	i += *incx;
/* L95: */
    }
    ret_val = sqrt(sum);
    goto L300;

L200:
    ++ix;
    i += *incx;
    if (ix <= *n) {
	goto L20;
    }

/*              end of main loop. */

/*              compute square root and adjust for scaling. */

    ret_val = xmax * sqrt(sum);
L300:
    return ret_val;
} /* dnrm2_ */


/* Subroutine */
static  int dtrmv_(uplo, trans, diag, n, a, lda, x, incx)
char *uplo, *trans, *diag;
integer *n;
double *a;
integer *lda;
double *x;
integer *incx;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static double temp;
    static integer i, j;
    extern logical blas_lsame_();
    static integer ix, jx, kx;
    extern /* Subroutine */ int xerbla_();
    static logical nounit;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DTRMV  performs one of the matrix-vector operations */

/*     x := A*x,   or   x := A'*x, */

/*  where x is an n element vector and  A is an n by n unit, or non-unit, 
*/
/*  upper or lower triangular matrix. */

/*  Parameters */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix is an upper or */
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the operation to be performed as */
/*           follows: */

/*              TRANS = 'N' or 'n'   x := A*x. */

/*              TRANS = 'T' or 't'   x := A'*x. */

/*              TRANS = 'C' or 'c'   x := A'*x. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit */
/*           triangular as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/*           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/*           upper triangular part of the array A must contain the upper 
*/
/*           triangular matrix and the strictly lower triangular part of 
*/
/*           A is not referenced. */
/*           Before entry with UPLO = 'L' or 'l', the leading n by n */
/*           lower triangular part of the array A must contain the lower 
*/
/*           triangular matrix and the strictly upper triangular part of 
*/
/*           A is not referenced. */
/*           Note that when  DIAG = 'U' or 'u', the diagonal elements of 
*/
/*           A are not referenced either, but are assumed to be unity. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared 
*/
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, n ). */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element vector x. On exit, X is overwritten with the */
/*           tranformed vector x. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --x;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    info = 0;
    if (! blas_lsame_(uplo, "U", 1L, 1L) && ! blas_lsame_(uplo, "L", 1L, 1L)) {
	info = 1;
    } else if (! blas_lsame_(trans, "N", 1L, 1L) && ! blas_lsame_(trans, "T", 1L, 1L) &&
	     ! blas_lsame_(trans, "C", 1L, 1L)) {
	info = 2;
    } else if (! blas_lsame_(diag, "U", 1L, 1L) && ! blas_lsame_(diag, "N", 1L, 1L)) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < max(1,*n)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    }
    if (info != 0) {
	xerbla_("DTRMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = blas_lsame_(diag, "N", 1L, 1L);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

    if (blas_lsame_(trans, "N", 1L, 1L)) {

/*        Form  x := A*x. */

	if (blas_lsame_(uplo, "U", 1L, 1L)) {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[j] != 0.) {
			temp = x[j];
			i__2 = j - 1;
			for (i = 1; i <= i__2; ++i) {
			    x[i] += temp * a[i + j * a_dim1];
/* L10: */
			}
			if (nounit) {
			    x[j] *= a[j + j * a_dim1];
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[jx] != 0.) {
			temp = x[jx];
			ix = kx;
			i__2 = j - 1;
			for (i = 1; i <= i__2; ++i) {
			    x[ix] += temp * a[i + j * a_dim1];
			    ix += *incx;
/* L30: */
			}
			if (nounit) {
			    x[jx] *= a[j + j * a_dim1];
			}
		    }
		    jx += *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (x[j] != 0.) {
			temp = x[j];
			i__1 = j + 1;
			for (i = *n; i >= i__1; --i) {
			    x[i] += temp * a[i + j * a_dim1];
/* L50: */
			}
			if (nounit) {
			    x[j] *= a[j + j * a_dim1];
			}
		    }
/* L60: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    if (x[jx] != 0.) {
			temp = x[jx];
			ix = kx;
			i__1 = j + 1;
			for (i = *n; i >= i__1; --i) {
			    x[ix] += temp * a[i + j * a_dim1];
			    ix -= *incx;
/* L70: */
			}
			if (nounit) {
			    x[jx] *= a[j + j * a_dim1];
			}
		    }
		    jx -= *incx;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := A'*x. */

	if (blas_lsame_(uplo, "U", 1L, 1L)) {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = x[j];
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    for (i = j - 1; i >= 1; --i) {
			temp += a[i + j * a_dim1] * x[i];
/* L90: */
		    }
		    x[j] = temp;
/* L100: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    for (i = j - 1; i >= 1; --i) {
			ix -= *incx;
			temp += a[i + j * a_dim1] * x[ix];
/* L110: */
		    }
		    x[jx] = temp;
		    jx -= *incx;
/* L120: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[j];
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = *n;
		    for (i = j + 1; i <= i__2; ++i) {
			temp += a[i + j * a_dim1] * x[i];
/* L130: */
		    }
		    x[j] = temp;
/* L140: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = *n;
		    for (i = j + 1; i <= i__2; ++i) {
			ix += *incx;
			temp += a[i + j * a_dim1] * x[ix];
/* L150: */
		    }
		    x[jx] = temp;
		    jx += *incx;
/* L160: */
		}
	    }
	}
    }

    return 0;

/*     End of DTRMV . */

} /* dtrmv_ */


static logical blas_lsame_(ca, cb, ca_len, cb_len)
char *ca, *cb;
ftnlen ca_len;
ftnlen cb_len;
{
    /* System generated locals */
    logical ret_val;

    /* Local variables */
    static integer inta, intb, zcode;


/*  -- LAPACK auxiliary routine (version 1.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 29, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  LSAME returns .TRUE. if CA is the same letter as CB regardless of */
/*  case. */

/*  Arguments */
/*  ========= */

/*  CA      (input) CHARACTER*1 */
/*  CB      (input) CHARACTER*1 */
/*          CA and CB specify the single characters to be compared. */

/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test if the characters are equal */

    ret_val = *ca == *cb;
    if (ret_val) {
	return ret_val;
    }

/*     Now test for equivalence if both characters are alphabetic. */

    zcode = 'Z';

/*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime */
/*     machines, on which ICHAR returns a value with bit 8 set. */
/*     ICHAR('A') on Prime machines returns 193 which is the same as */
/*     ICHAR('A') on an EBCDIC machine. */

    inta = *ca;
    intb = *cb;

    if (zcode == 90 || zcode == 122) {

/*        ASCII is assumed - ZCODE is the ASCII code of either lower o
r */
/*        upper case 'Z'. */

	if (inta >= 97 && inta <= 122) {
	    inta += -32;
	}
	if (intb >= 97 && intb <= 122) {
	    intb += -32;
	}

    } else if (zcode == 233 || zcode == 169) {

/*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower
 or */
/*        upper case 'Z'. */

	if (inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta 
		>= 162 && inta <= 169) {
	    inta += 64;
	}
	if (intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb 
		>= 162 && intb <= 169) {
	    intb += 64;
	}

    } else if (zcode == 218 || zcode == 250) {

/*        ASCII is assumed, on Prime machines - ZCODE is the ASCII cod
e */
/*        plus 128 of either lower or upper case 'Z'. */

	if (inta >= 225 && inta <= 250) {
	    inta += -32;
	}
	if (intb >= 225 && intb <= 250) {
	    intb += -32;
	}
    }
    ret_val = inta == intb;

/*     RETURN */

/*     End of LSAME */

    return ret_val;
} /* lsame_ */

#endif

#ifdef KR_headers
static double d_sign(a,b) double *a, *b;
#else
static double d_sign(double *a, double *b)
#endif
{
double x;
x = (*a >= 0 ? *a : - *a);
return( *b >= 0 ? x : -x);
}


/* compare two strings */

#ifdef KR_headers
static integer s_cmp(a0, b0, la, lb) char *a0, *b0; ftnlen la, lb;
#else
static integer s_cmp(char *a0, char *b0, ftnlen la, ftnlen lb)
#endif
{
register unsigned char *a, *aend, *b, *bend;
a = (unsigned char *)a0;
b = (unsigned char *)b0;
aend = a + la;
bend = b + lb;

if(la <= lb)
	{
	while(a < aend)
		if(*a != *b)
			return( *a - *b );
		else
			{ ++a; ++b; }

	while(b < bend)
		if(*b != ' ')
			return( ' ' - *b );
		else	++b;
	}

else
	{
	while(b < bend)
		if(*a == *b)
			{ ++a; ++b; }
		else
			return( *a - *b );
	while(a < aend)
		if(*a != ' ')
			return(*a - ' ');
		else	++a;
	}
return(0);
}


#ifdef KR_headers
static double pow_di(ap, bp) double *ap; integer *bp;
#else
static double pow_di(double *ap, integer *bp)
#endif
{
double pow, x;
integer n;
unsigned long u;

pow = 1;
x = *ap;
n = *bp;

if(n != 0)
	{
	if(n < 0)
		{
		n = -n;
		x = 1/x;
		}
	for(u = n; ; )
		{
		if(u & 01)
			pow *= x;
		if(u >>= 1)
			x *= x;
		else
			break;
		}
	}
return(pow);
}


/* assign strings:  a = b */

#ifdef KR_headers
static int s_copy(a, b) char *a, *b; ftnlen la, lb;
#else
static int s_copy(char *a, char *b, ftnlen la, ftnlen lb)
#endif
{
	register char *aend, *bend;

	aend = a + la;

	if(la <= lb)
#ifndef NO_OVERWRITE
		if (a <= b || a >= b + la)
#endif
			while(a < aend)
				*a++ = *b++;
#ifndef NO_OVERWRITE
		else
			for(b += la; a < aend; )
				*--aend = *--b;
#endif

	else {
		bend = b + lb;
#ifndef NO_OVERWRITE
		if (a <= b || a >= bend)
#endif
			while(b < bend)
				*a++ = *b++;
#ifndef NO_OVERWRITE
		else {
			a += lb;
			while(b < bend)
				*--a = *--bend;
			a += lb;
			}
#endif
		while(a < aend)
			*a++ = ' ';
		}
        return 0;
	}
