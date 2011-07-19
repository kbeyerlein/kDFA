/*
 * globalfunctions.cpp
 *
 *  Created on: Jul 19, 2011
 *      Author: ken
 */
#include "globalfunctions.h"
double Samplemax (double size, double wv,  int pkPts)
/*
 * For a cubic particle, gives smallest peak size, well smallest
 */
{
	double ib = .83*wv/size;
	double sigma = sqrt(ib/(2*M_PI));
	return 6*sigma/pkPts;
}
/*
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
*/
void interp (double *inx, double *iny, int ninvalues, double *xs, double *outvalues, int nxs )
{
	for (int i = 0; i<nxs; i++)
	{
		int adv = 0;
		// Bound loop to prevent infinite loop for case when interpolation array does not have sufficient data.
		for ( ; inx[adv]<xs[i]&&adv<ninvalues; adv++) ;
		// Avoid bad address for iny[adv-1] when adv=0
		if (adv==0){
			outvalues[i]=iny[adv]+(iny[adv+1]-iny[adv])*(xs[i]-inx[adv])/(inx[adv+1]-inx[adv]);
		}
		// Most likely case
		else if (adv!=ninvalues-1){
			outvalues[i]=iny[adv-1]+(iny[adv]-iny[adv-1])*(xs[i]-inx[adv-1])/(inx[adv]-inx[adv-1]);
		}
		// Check if last point is really within bounds
		else if (adv==ninvalues-1 && inx[adv]>xs[i]){
			outvalues[i]=iny[adv-1]+(iny[adv]-iny[adv-1])*(xs[i]-inx[adv-1])/(inx[adv]-inx[adv-1]);
		}
		else{
			cerr << "interpolation fail" << endl;
			exit(0);
		}
	}

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
