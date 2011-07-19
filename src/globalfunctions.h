/*
 * globalfunctions.h
 *
 *  Created on: Jul 19, 2011
 *      Author: fenrisulfr
 */

#ifndef GLOBALFUNCTIONS_H_
#define GLOBALFUNCTIONS_H_
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <fstream>
#include <iostream>
#include <cstring>

double Samplemax (double size, int pkPts, double wv);
void interp (double *xs,int nxs,double *inx, double *iny, int ninvalues, double *outvalues);
void splinterp (double *xs,int nxs,double *inx, double *iny, int nin, double *outvalues);
void linterp (double *xs,int nxs,double *inx, double *iny, int nin, double *outvalues);
bool addelement(double**myarray, double newelement, unsigned int &elements);

double Samplemax (double size, double wv,  int pkPts)
/*
 * For a cubic particle, gives smallest peak size, well smallest
 */
{
	double ib = .83*wv/size;
	double sigma = sqrt(ib/(2*M_PI));
	return 6*sigma/pkPts;
}
void interp (double *xs,int nxs,double *inx, double *iny, int ninvalues, double **outvalues)
{
	unsigned int nixes = 0;
	for (int i = 0; i<nxs; i++)
	{
		int adv = 0;
		for ( ; inx[adv]<xs[i]; adv++) ;
		addelement(outvalues,iny[adv-1]+(iny[adv]-iny[adv-1])*(xs[i]-inx[adv-1])/(inx[adv]-inx[adv-1]),nixes);
	}
	if(nixes != (unsigned int) nxs)
		cerr << "interpolation fail" << endl;
	return;
}

void splinterp (double *xs,int nxs,double *inx, double *iny, int nin, double **outvalues)
{
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *splointer = gsl_spline_alloc (gsl_interp_cspline, nin);
	gsl_spline_init(splointer, inx, iny, nin);
	for (int i = 0; i<nxs; i++)
		{
		(*outvalues)[i]=gsl_spline_eval(splointer,xs[i],acc);
		}
	gsl_spline_free(splointer);
	gsl_interp_accel_free(acc);
	return;
}
void linterp (double *xs,int nxs,double *inx, double *iny, int nin, double **outvalues)
{
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *splointer = gsl_spline_alloc (gsl_interp_linear, nin);
	gsl_spline_init(splointer, inx, iny, nin);
	for (int i = 0; i<nxs; i++)
		{
		(*outvalues)[i]=gsl_spline_eval(splointer,xs[i],acc);
		}
	gsl_spline_free(splointer);
	gsl_interp_accel_free(acc);
	return;
}
void akimterp (double *xs,int nxs,double *inx, double *iny, int nin, double **outvalues)
{
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *splointer = gsl_spline_alloc (gsl_interp_akima, nin);
	gsl_spline_init(splointer, inx, iny, nin);
	for (int i = 0; i<nxs; i++)
		{
		(*outvalues)[i]=gsl_spline_eval(splointer,xs[i],acc);
		}
	gsl_spline_free(splointer);
	gsl_interp_accel_free(acc);
	return;
}
bool addelement(double **myarray, double newelement, unsigned int &elements)
{
	double *newarray = new double[elements+1];
	if (!newarray) return false;
	//if (elements>0)
	for(unsigned int i=0; i < elements; i++)
		newarray[i]=(*myarray)[i];

	delete [] *myarray;

	*myarray=newarray;
	(*myarray)[elements]=newelement;
//	newarray=NULL;
	elements++;
	return true;
}


#endif /* GLOBALFUNCTIONS_H_ */
