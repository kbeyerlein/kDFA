/*
 * globalfunctions.h
 *
 *  Created on: Jul 19, 2011
 *      Author: fenrisulfr
 */

#ifndef GLOBALFUNCTIONS_H_
#define GLOBALFUNCTIONS_H_

#include "Includes.h"
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <fstream>
#include <iostream>
#include <cstring>

double Samplemax (double size, int pkPts, double wv);
//void interp (double *xs,int nxs,double *inx, double *iny, int ninvalues, double *outvalues);
void interp (double *, double *, int, double *, double *, int);
void splinterp (double *xs,int nxs,double *inx, double *iny, int nin, double *outvalues);
void linterp (double *xs,int nxs,double *inx, double *iny, int nin, double *outvalues);
bool addelement(double**myarray, double newelement, unsigned int &elements);




#endif /* GLOBALFUNCTIONS_H_ */
