
#include <stdio.h>
#include <R.h>
#include <math.h>


//int main(int argc, char **argv)
//{
//
//	return 0;
//}

double phi(
		double mu,
		double sigma_sq,
		double x){
	return(1/sqrt(2.0*PI*sigma_sq)*exp(-pow((x-mu),2)  / (2.0*sigma_sq)));
}

void likelihood(
		double *p1,
		double *p2,
		double *p3,
		double *mu,
		double *A,
		double *B,
		double *C,
		double *rawData,
		int *n,
		double *result) {
	*result=0;
	double temp=0.0;
	for(int i=0; i < *n; ++i)  {
		temp= *p1*phi(0.0,*A,rawData[i]) + *p2*phi(0.0,*B,rawData[i]) + *p3*phi(*mu,*C,rawData[i]);
		*result+= log(temp);
	}
}
