#ifndef COMPARE_H
#define COMPARE_H

#include "Shape.h"
#include "Param.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

class Compare
{
public:
	Compare(string, string, Shape*, int);
	Compare(Shape *, int);
	~Compare(void);

	int nParam;
	Pattern **initCalcI, *obsI, *finalCalcI;
	Param **indParams;
	float rwp, rexp, chisq;
	double **gradI, **A, *g, *delta, nAtoms;
	gsl_matrix *covar;
	bool init;
	Shape *firstShape;
	string name, path;

	double *derivfSize_mu, *derivfSize_sig;

	int LevMarq(int);
	int InputObservedInt(string, string);
	void CalcRwp();
	int InitParams();
	int CalcPattern(bool, bool);
	int OutputPattern(string, string, int, bool);
	int CalcDerivW(Shape *);
	double GetDerivNumAtoms_mu(Shape *);
	double GetDerivNumAtoms_sig(Shape *);
	void CalcGrad();
	int CalcDeriv_B(double *);
	int CalcDeriv_t(double *);
	int CalcDeriv_a(double *, string);
	int CalcDeriv_a_D(double*, int, string);
	int CalcDeriv_mu(double *, string);
	int CalcDeriv_sigma(double *, string);
	int CalcDeriv_alpha(double *, string);
	int CalcDeriv_beta(double *, string);
	int CalcDeriv_deltaSize(double *, string);
	int CalcDeriv_f(double *, string);
	int CalcDeriv_kappa(double *, string);
	int CalcDeriv_wShape(double *, string);
	void CalcDeriv_scale(double *);
	void CalcDeriv_ChebScale(double *, int);
	int CalcDeriv_smplDisplace(double *);
	void CalcDeriv_gaussA(double *, int);
	void CalcDeriv_gaussMu(double *, int);
	void CalcDeriv_gaussSig(double *, int);
	void CalcDeriv_expAmp(double *);
	void CalcDeriv_expKappa(double *);
	void CalcDeriv_expDel(double *);
	void CopyShapeInfo(Shape *, Shape *);
	void CalcA();
	void Calcg();
	int CholeskyDecomp(double **, int);
	int CholeskyBackSub(double **, int, double *, double*);
	void CheckParamBounds();
	void AdjustShapeWeights();
	void NormalizeShapeWeights();
	double Cosc(double);
	int SymDiagDecomp(double **, int, double*, double*);
	int SymDiagBackSub(double**, double*, int, double*, double*);
	void CalcResidual();
	double CalcNumAtoms();
	void CleanGlobalParam(Param**);
	void CalcCovarMatrix(double);
	void SolveForDelta(double **, double *, double *, int);
};

#endif
