/*
 Author Shlomi Lifshits (lifshits_il@yahoo.com)
 Date: September 2012 
*/

#ifndef AVERAGE_NA
#define AVERAGE_NA

#define NA_NUM 20
#define NEIG_SIZE 26

#define SUB2IND(i,j,k,dim_x,dim_y) (((k)-1)*(dim_x)*(dim_y)+((j)-1)*(dim_x)+(i)-1)

typedef struct
{
  double vals[NEIG_SIZE];
  int num_in_vals;
  int do_average;
} neig_val_arr;

typedef neig_val_arr* neig_val_arr_p;

static void get_neighbors_vals(SEXP vec_p,int i, int j, int k, int dim_x, int dim_y, neig_val_arr_p neig_p);
static double calc_average_over_neig(neig_val_arr_p neig_p);
SEXP fill_NA_voxels_with_neig_mean(SEXP vec_p, SEXP dim_x, SEXP dim_y, SEXP dim_z);

#endif
