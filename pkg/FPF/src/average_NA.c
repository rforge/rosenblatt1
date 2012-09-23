/* 
Author Shlomi Lifshits (lifshits_il@yahoo.com)
Date: September 2012 
*/

#include <stdio.h>
#include <stdlib.h>

#include <R.h>
#include <Rinternals.h>

#include "average_NA.h"

static void get_neighbors_vals(SEXP vec_p, int i, int j, int k, int dim_x, int dim_y, neig_val_arr_p neig_p)
{
  int ii,jj,kk;
  int linear_ind;
  
  neig_p->do_average=FALSE;
  neig_p->num_in_vals=0;
  
  for (ii=(i-1);ii<=(i+1);ii++)
    {
      for (jj=(j-1);jj<=(j+1);jj++)
	{
	  for (kk=(k-1);kk<=(k+1);kk++)
	    {
	      if (ii==i && jj==j && kk==k)
		continue;
	      
	      linear_ind=SUB2IND(ii,jj,kk,dim_x,dim_y);
	      
	      if (!ISNA(REAL(vec_p)[linear_ind]))
		{
		  neig_p->vals[neig_p->num_in_vals]=REAL(vec_p)[linear_ind];
		  (neig_p->num_in_vals)++;
		}
	    }
	}
    }
  
  if ((neig_p->num_in_vals)>NA_NUM)
    {
      neig_p->do_average=TRUE;
    }
  
  return;
}

static double calc_average_over_neig(neig_val_arr_p neig_p)
{
  int i;
  double sum=0.0;
  int num_in_vals=neig_p->num_in_vals;

  for (i=0;i<num_in_vals;i++)
    {
      sum+=neig_p->vals[i];
    }

  return(sum/num_in_vals);
}

SEXP fill_NA_voxels_with_neig_mean(SEXP vec_p, SEXP dim_x_p, SEXP dim_y_p, SEXP dim_z_p)
{
  int i,j,k,l,my_linear_index,len;
  neig_val_arr neig;
  SEXP out_vec;
  int dim_x,dim_y,dim_z;
  int num_protected=0;

  dim_x=INTEGER(dim_x_p)[0];
  dim_y=INTEGER(dim_y_p)[0];
  dim_z=INTEGER(dim_z_p)[0];
  len=dim_x*dim_y*dim_z;

  PROTECT(out_vec = allocVector(REALSXP, len));
  num_protected++;

  /*copy input vector*/
  for (l=0;l<len;l++)
    {
      REAL(out_vec)[l]=REAL(vec_p)[l];
    }

  for (k=2;k<=(dim_z-1);k++)
    {
      for (j=2;j<=(dim_y-1);j++)
	{
	  for (i=2;i<=(dim_x-1);i++)
	    {
              my_linear_index= SUB2IND(i,j,k,dim_x,dim_y);

	      if (ISNA(REAL(out_vec)[my_linear_index]))
		{
		  get_neighbors_vals(out_vec,i,j,k,dim_x,dim_y,&neig);

		  if (neig.do_average==TRUE)
		    {
		      REAL(out_vec)[my_linear_index]=calc_average_over_neig(&neig);
		    }
		}
	    }
	}
    }

  UNPROTECT(num_protected);
  return(out_vec);
}
