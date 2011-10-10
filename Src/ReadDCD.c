/***************************************************************************
 *cr                                                                       
 *cr            (C) Copyright 1995 The Board of Trustees of the           
 *cr                        University of Illinois                       
 *cr                         All Rights Reserved                        
 *cr                                                                   
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: ReadDCD.C,v $
 *      $Author: dalke $        $Locker:  $                $State: Exp $
 *      $Revision: 1.6 $      $Date: 1997/03/13 17:38:56 $
 *
 ***************************************************************************
 *
 * Modifications for use in MMTK:
 *
 * - Argument DELTA in read_dcdheader changed to FLOAT.
 * - Check for file existence in open_dcd_write removed.
 * - Include file name changed to mmtk_readdcd.h
 * - open/read/write/close replaced by fopen/fread/fwrite/fclose
 *   (for portability to Windows)
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 * C routines to read and write binary DCD files (which use the goofy
 * FORTRAN UNFORMATTED format).  These routines are courtesy of
 * Mark Nelson (autographs available upon request and $1 tip).
 *
 ***************************************************************************/

#include <string.h>
#include "MMTK/readdcd.h"



void pad(char *s, int len)
{
	int curlen;
	int i;

	curlen=strlen(s);

	if (curlen>len)
	{
		s[len]=0;
		return;
	}

	for (i=curlen; i<len; i++)
	{
		s[i]=' ';
	}

	s[i]=0;
}

/************************************************************************/
/*									*/
/*			FUNCTION open_dcd_read				*/
/*									*/
/*   INPUTS:								*/
/*	filename - filename to read in					*/
/*									*/
/*   OUTPUTS:								*/
/*	an open filedesriptor (integer) is returned if the call is      */
/*   successful, otherwise a negative number is returned		*/
/*									*/
/*	open_dcd_read opens a dcd file for input.  It first does a check*/
/*   to see if the file really exists.  If the file does not exist,     */
/*   a DCD_DNE is returned, if the file exists but can' be opened for   */
/*   some reason, a DCD_OPENFAILED is returned.  If the open is         */
/*   successful, the file descriptor is returned.			*/
/*									*/
/************************************************************************/

FILE *open_dcd_read(const char *filename)

{
	FILE *dcdfd;		/*  file descriptor for dcd file	*/

	dcdfd = fopen(filename, "rb");
	return(dcdfd);
}

/****************************************************************/
/*								*/
/*			FUNCTION read_dcdheader			*/
/*								*/
/*   INPUTS:							*/
/*	fd - file descriptor for the dcd file			*/
/*	N - Number of atoms					*/
/*	NSET - Number of sets of coordinates			*/
/*	ISTART - Starting timestep of DCD file			*/
/*	NSAVC - Timesteps between DCD saves			*/
/*	DELTA - length of a timestep				*/
/*								*/
/*   OUTPUTS:							*/
/*	N, NSET, ISTART, NSAVC, and DELTA are returned populated*/
/*   a 0 is returned if the call is successful, or a negative   */
/*   number if errors are detected				*/
/*								*/
/*	read_header reads in the header from a DCD file and     */
/*   returns the timestep size and the number of atoms in the   */
/*   system.  A 0 is returned if the header is successfully     */
/*   read, or a negative error code is returned otherwise.      */
/*								*/
/****************************************************************/

int read_dcdheader(FILE *fd, int *N, int *NSET, int *ISTART, 
		   int *NSAVC, float *DELTA, int *NAMNF, 
		   int **FREEINDEXES)

{
	int input_integer;	/*  Integer buffer space	*/
	char bigbuf[256];	/*  A large string buffer	*/
	int ret_val;		/*  Return value from read	*/
	int i;			/*  Loop counter		*/
	char HDR[5];		/*  Header = "CORD"		*/
	int I;
	int NTITLE;

	/*  First thing in the file should be an 84		*/
	ret_val = fread(&input_integer, sizeof(int), 1, fd);
	if (ret_val != 1)
	  return(DCD_BADREAD);

	if (input_integer != 84)
	  return(DCD_BADFORMAT);

	/*  Read in the string "CORD"			*/
	ret_val = fread(HDR, sizeof(char), 4, fd);
	if (ret_val != 4)
	  return(DCD_BADREAD);
	HDR[4]=0;


	/*  Read in the number of Sets of coordinates, NSET  */
	ret_val = fread(NSET, sizeof(int), 1, fd);
	if (ret_val != 1)
	  return(DCD_BADREAD);

	/*  Read in ISTART, the starting timestep	     */
	ret_val = fread(ISTART, sizeof(int), 1, fd);
	if (ret_val != 1)
	  return(DCD_BADREAD);

	/*  Read in NSAVC, the number of timesteps between   */
	/*  dcd saves					     */
	ret_val = fread(NSAVC, sizeof(int), 1, fd);
	if (ret_val != 1)
	  return(DCD_BADREAD);


	/*  Skip blank integers				     */
	for (i=0; i<5; i++)
	{
		ret_val = fread(&input_integer, sizeof(int), 1, fd);
		if (ret_val != 1)
		  return(DCD_BADREAD);
	}

	/*  Read NAMNF, the number of free atoms	     */
	ret_val = fread(NAMNF, sizeof(int), 1, fd);
	if (ret_val != 1)
	  return(DCD_BADREAD);

	/*  Read in the timestep, DELTA				*/
	ret_val = fread(DELTA, sizeof(float), 1, fd);
	if (ret_val != 1)
	  return(DCD_BADREAD);

	/*  Skip blank integers					*/
	for (i=0; i<10; i++)
	{
		ret_val = fread(&I, sizeof(int), 1, fd);
		if (ret_val != 1)
		  return(DCD_BADREAD);
	}

	/*  Get the end size of the first block			*/
	ret_val = fread(&input_integer, sizeof(int), 1, fd);
	if (ret_val != 1)
	  return(DCD_BADREAD);

	if (input_integer != 84)
	{
		return(DCD_BADFORMAT);
	}

	/*  Read in the size of the next block			*/
	ret_val = fread(&input_integer, sizeof(int), 1, fd);
	if (ret_val != 1)
	  return(DCD_BADREAD);

	if ( ((input_integer-4)%80) == 0)
	{
		/*  Read NTITLE, the number of 80 characeter    */
		/*  title strings there are			*/
		ret_val = fread(&NTITLE, sizeof(int), 1, fd);
		if (ret_val != 1)
		  return(DCD_BADREAD);

		for (i=0; i<NTITLE; i++)
		{
			ret_val = fread(bigbuf, sizeof(char), 80, fd);
			if (ret_val != 80)
			  return(DCD_BADREAD);
		}

		/*  Get the ending size for this block		*/
		ret_val = fread(&input_integer, sizeof(int), 1, fd);
		if (ret_val != 1)
		  return(DCD_BADREAD);
	}
	else
	{
		return(DCD_BADFORMAT);
	}

	/*  Read in a 4				*/
	ret_val = fread(&input_integer, sizeof(int), 1, fd);
	if (ret_val != 1)
	  return(DCD_BADREAD);

	if (input_integer != 4)
	{
		return(DCD_BADFORMAT);
	}

	/*  Read in the number of atoms			*/
	ret_val = fread(N, sizeof(int), 1, fd);
	if (ret_val != 1)
	  return(DCD_BADREAD);

	/*  Read in an 4				*/
	ret_val = fread(&input_integer, sizeof(int), 1, fd);
	if (ret_val != 1)
	  return(DCD_BADREAD);

	if (input_integer != 4)
	{
		return(DCD_BADFORMAT);
	}

	if (*NAMNF != 0)
	{
		(*FREEINDEXES) = (int *) calloc(((*N)-(*NAMNF)), sizeof(int));

		if (*FREEINDEXES == NULL)
			return(DCD_BADMALLOC);
	
		/*  Read in an size */
		ret_val = fread(&input_integer, sizeof(int), 1, fd);
		if (ret_val != 1)
		  return(DCD_BADREAD);

		if (input_integer != ((*N)-(*NAMNF))*4)
		{
			return(DCD_BADFORMAT);
		}
		
		ret_val = fread((*FREEINDEXES), sizeof(int), (*N)-(*NAMNF), fd);
		if (ret_val != (*N)-(*NAMNF))
		  return(DCD_BADREAD);

		ret_val = fread(&input_integer, sizeof(int), 1, fd);
		if (ret_val != 1)
		  return(DCD_BADREAD);

		if (input_integer != ((*N)-(*NAMNF))*4)
		{
			return(DCD_BADFORMAT);
		}
	}

	return(0);
}

/************************************************************************/
/*									*/
/*			FUNCTION read_dcdstep				*/
/*									*/
/*   INPUTS:								*/
/*	fd - file descriptor to use					*/
/*	N - Number of atoms						*/
/*	X - array of X coordinates					*/
/*	Y - array of Y coordinates					*/
/*	Z - array of Z coordinates					*/
/*	num_fixed - Number of fixed atoms				*/
/*	first - flag 1->first time called				*/
/*	indexes - indexes for free atoms				*/
/*									*/
/*   OUTPUTS:								*/
/*	read step populates the arrays X, Y, Z and returns a 0 if the   */
/*   read is successful.  If the read fails, a negative error code is   */
/*   returned.								*/
/*									*/
/*	read_step reads in the coordinates from one step.  It places    */
/*   these coordinates into the arrays X, Y, and Z.			*/
/*									*/
/************************************************************************/

int read_dcdstep(FILE *fd, int N, float *X, float *Y, float *Z, int num_fixed,
		 int first, int *indexes)

{
	int ret_val;		/*  Return value from read		*/
	int input_integer;	/*  Integer buffer space		*/
	int i;			/*  Loop counter			*/
	static float *tmpX;

	if (first && num_fixed)
	{
		tmpX = (float *) calloc(N-num_fixed, sizeof(float));

		if (tmpX==NULL)
		{
			return(DCD_BADMALLOC);
		}
	}

	/*  Get the first size from the file				*/
	ret_val = fread(&input_integer, sizeof(int), 1, fd);

	/*  See if we've reached the end of the file			*/
	if (ret_val == 0)
	{
		free(tmpX);

		return(-1);
	}

	if ( (num_fixed==0) || first)
	{
		if (input_integer != 4*N)
		{
			return(DCD_BADFORMAT);
		}

		ret_val = fread(X, sizeof(float), N, fd);
		if (ret_val != N)
		  return(DCD_BADREAD);
	
		ret_val = fread(&input_integer, sizeof(int), 1, fd);
		if (ret_val != 1)
		  return(DCD_BADREAD);

		if (input_integer != 4*N)
		{
			return(DCD_BADFORMAT);
		}

		ret_val = fread(&input_integer, sizeof(int), 1, fd);
		if (ret_val != 1)
		  return(DCD_BADREAD);

		if (input_integer != 4*N)
		{
			return(DCD_BADFORMAT);
		}

		ret_val = fread(Y, sizeof(float), N, fd);
		if (ret_val != N)
		  return(DCD_BADREAD);

		ret_val = fread(&input_integer, sizeof(int), 1, fd);
		if (ret_val != 1)
		  return(DCD_BADREAD);

		if (input_integer != 4*N)
		{
			return(DCD_BADFORMAT);
		}

		ret_val = fread(&input_integer, sizeof(int), 1, fd);
		if (ret_val != 1)
		  return(DCD_BADREAD);

		if (input_integer != 4*N)
		{
			return(DCD_BADFORMAT);
		}

		ret_val = fread(Z, sizeof(float), N, fd);
		if (ret_val != N)
		  return(DCD_BADREAD);

		ret_val = fread(&input_integer, sizeof(int), 1, fd);
		if (ret_val != 1)
		  return(DCD_BADREAD);

		if (input_integer != 4*N)
		{
			return(DCD_BADFORMAT);
		}
	}
	else
	{
		if (input_integer != 4*(N-num_fixed))
		{
			return(DCD_BADFORMAT);
		}

		ret_val = fread(tmpX, sizeof(float), N-num_fixed, fd);
		if (ret_val != N-num_fixed)
		  return(DCD_BADREAD);
	
		for (i=0; i<N-num_fixed; i++)
		{
			X[indexes[i]-1]=tmpX[i];
		}

		ret_val = fread(&input_integer, sizeof(int), 1, fd);
		if (ret_val != 1)
		  return(DCD_BADREAD);

		if (input_integer != 4*(N-num_fixed))
		{
			return(DCD_BADFORMAT);
		}

		ret_val = fread(&input_integer, sizeof(int), 1, fd);
		if (ret_val != 1)
		  return(DCD_BADREAD);

		if (input_integer != 4*(N-num_fixed))
		{
			return(DCD_BADFORMAT);
		}

		ret_val = fread(tmpX, sizeof(float), N-num_fixed, fd);
		if (ret_val != N-num_fixed)
		  return(DCD_BADREAD);
	
		for (i=0; i<N-num_fixed; i++)
		{
			Y[indexes[i]-1]=tmpX[i];
		}

		ret_val = fread(&input_integer, sizeof(int), 1, fd);
		if (ret_val != 1)
		  return(DCD_BADREAD);

		if (input_integer != 4*(N-num_fixed))
		{
			return(DCD_BADFORMAT);
		}

		ret_val = fread(&input_integer, sizeof(int), 1, fd);
		if (ret_val != 1)
		  return(DCD_BADREAD);

		if (input_integer != 4*(N-num_fixed))
		{
			return(DCD_BADFORMAT);
		}

		ret_val = fread(tmpX, sizeof(float), N-num_fixed, fd);
		if (ret_val != N-num_fixed)
		  return(DCD_BADREAD);

		for (i=0; i<N-num_fixed; i++)
		{
			Z[indexes[i]-1]=tmpX[i];
		}

		ret_val = fread(&input_integer, sizeof(int), 1, fd);
		if (ret_val != 1)
		  return(DCD_BADREAD);

		if (input_integer != 4*(N-num_fixed))
		{
			return(DCD_BADFORMAT);
		}
	}

	return(0);
}


/*********************************************************************/
/*								     */
/*			FUNCTION open_dcd_write			     */
/*								     */
/*   INPUTS:							     */
/*	dcdfile - Name of the dcd file				     */
/*								     */
/*   OUTPUTS:							     */
/*	returns an open file descriptor for writing		     */
/*								     */
/*	This function will open a dcd file for writing.  It takes    */
/*   the filename to open as its only argument.	 It will return a    */
/*   valid file descriptor if successful or DCD_OPENFAILED if the    */
/*   open fails for some reason.                                     */
/*								     */
/*********************************************************************/

FILE *open_dcd_write(char *dcdname)

{
	return fopen(dcdname, "wb");
}

/************************************************************************/
/*									*/
/*				FUNCTION write_dcdstep			*/
/*									*/
/*   INPUTS:								*/
/*	fd - file descriptor for the DCD file to write to		*/
/*	N - Number of atoms						*/
/*	X - X coordinates						*/
/*	Y - Y coordinates						*/
/*	Z - Z coordinates						*/
/*									*/
/*   OUTPUTS:								*/
/*	none								*/
/*									*/
/*	write_dcdstep writes the coordinates out for a given timestep   */
/*   to the specified DCD file.						*/
/*									*/
/************************************************************************/

int write_dcdstep(FILE *fd, int N, float *X, float *Y, float *Z)

{
	int out_integer;

	out_integer = N*4;
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) X, 4, N, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) Y, 4, N, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) Z, 4, N, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	
	return 1;
}

/*****************************************************************************/
/*									     */
/*				FUNCTION write_dcdheader		     */
/*									     */
/*   INPUTS:								     */
/*	fd - file descriptor for the dcd file				     */
/*	N - Number of atoms						     */
/*	NSET - Number of sets of coordinates				     */
/*	ISTART - Starting timestep of DCD file				     */
/*	NSAVC - Timesteps between DCD saves				     */
/*	DELTA - length of a timestep					     */
/*									     */
/*   OUTPUTS:								     */
/*	none								     */
/*									     */
/*	This function prints the "header" information to the DCD file.  Since*/
/*   this is duplicating an unformatted binary output from FORTRAN, its ugly.*/
/*   So if you're squimish, don't look.					     */
/*									     */
/*****************************************************************************/

int write_dcdheader(FILE *fd, char *filename, int N, int NSET, int ISTART, 
		   int NSAVC, double DELTA)
{
	int	out_integer;
	char	title_string[200];
	time_t 	cur_time;
	struct  tm *tmbuf;
	char    time_str[11];
	float   DELTA_float = (float)DELTA;

	out_integer = 84;
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	strcpy(title_string, "CORD");
	fwrite(title_string, sizeof(char), 4, fd);
	fwrite((char *) &NSET, sizeof(int), 1, fd);
	fwrite((char *) &ISTART, sizeof(int), 1, fd);
	fwrite((char *) &NSAVC, sizeof(int), 1, fd);
	out_integer=0;
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) &DELTA_float, sizeof(float), 1, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	fwrite((char *) &out_integer, sizeof(int), 1, fd);
	out_integer = 84;
	fwrite((char *) & out_integer, sizeof(int), 1, fd);

	out_integer = 164;
	fwrite((char *) & out_integer, sizeof(int), 1, fd);
	out_integer = 2;
	fwrite((char *) & out_integer, sizeof(int), 1, fd);

	snprintf(title_string, sizeof(title_string),
                 "REMARKS FILENAME=%s CREATED BY VMD", filename);
	pad(title_string, 80);
	fwrite(title_string, sizeof(char), 80, fd);

	cur_time=time(NULL);
	tmbuf=localtime(&cur_time);
	strftime(time_str, 10, "%m/%d/%y", tmbuf);

	snprintf(title_string, sizeof(title_string),
                 "REMARKS DATE: %s CREATED BY MMTK.", time_str);
	pad(title_string, 80);
	fwrite(title_string, sizeof(char), 80, fd);
	out_integer = 164;
	fwrite((char *) & out_integer, sizeof(int), 1, fd);
	out_integer = 4;
	fwrite((char *) & out_integer, sizeof(int), 1, fd);
	out_integer = N;
	fwrite((char *) & out_integer, sizeof(int), 1, fd);
	out_integer = 4;
	fwrite((char *) & out_integer, sizeof(int), 1, fd);
	
	return 1;
}

/****************************************************************/
/*								*/
/*			FUNCTION close_dcd_read			*/
/*								*/
/*   INPUTS:							*/
/*	fd - file descriptor to close				*/
/*	num_fixed - the number of fixed atoms			*/
/*	indexes - Array of indexes to be deallocated		*/
/*								*/
/*   OUTPUTS:							*/
/*	the file pointed to by fd is closed and the memory      */
/*   pointed to by indexes is freed if any was allocated	*/
/*								*/
/*	close_dcd_read closes a dcd file that was opened for    */
/*   reading.  It also deallocated memory used for the indexes  */
/*   used for the free atom list, if there were fixed atoms.    */
/*								*/
/****************************************************************/

void close_dcd_read(FILE *fd, int num_fixed, int *indexes)

{
	fclose(fd);

	if (num_fixed)
	{
		free(indexes);
	}
}

/****************************************************************/
/*								*/
/*			FUNCTION close_dcd_write		*/
/*								*/
/*   INPUTS:							*/
/*	fd - file descriptor to close				*/
/*								*/
/*   OUTPUTS:							*/
/*	the file pointed to by fd				*/
/*								*/
/*	close_dcd_write close a dcd file that was opened for    */
/*   writing							*/
/*								*/
/****************************************************************/

void close_dcd_write(FILE *fd)

{
	fclose(fd);
}

