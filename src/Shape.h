#ifndef SHAPE_H
#define SHAPE_H

#include "Includes.h"
#include "Param.h"
#include "Distance.h"
#include "Position.h"
#include "Pattern.h"
#include "Global.h"




class Shape
{
public:
	Shape();
	~Shape(void);

	string name, type, posPath, distPath, pattPath, element;
	int maxShell, nRelaxShells, endSize, nPreCalcPatts;
	int nParams_a;
	Param **a, *wShape, *mu, *sigma, *f, *kappa, *n_R,*alpha, *beta, *deltaSize;
	Distance *firstDist, *surfRelaxDist, *derivSRD_mu, *derivSRD_sig;
	Position *firstPos;
	Shape *nextShape;
	Pattern **diffPatt, **preCalcPatts;

	double *wShell, *p, *dW_mu, *dW_sig, *D, *a_D, nAtoms;
	int *nA_D;
	string distribution;

	int LoadPositions();
	void LoadDistances();
	void FindPreCalcPatts();
	void CheckPreCalcPatts();
	void LoadPatterns();
	void LoadIntFromFile(Pattern *, double);
	int InitParamBounds();
	int SurfaceRelaxation(Distance *);
	int SurfaceRelaxation(double *, double *);
	void CopyPosition(Position *, Position *);
	void InitSRDistance(Distance *);
	double GetShapeShellD();
	double GetShapeShellL();
	double CalcTotVol();
	double GetUnitlessVolume(int);
	void CalcNumAtoms();
	void CalcSizeDistribution();
	void CalcShellWeights();
	void NormalizeSizeDistribution();
	int WriteIntFileHeader(string, string);
	void CalcEndSize();
	double GammLn(double);
	void GammaDist(double *, double *, double*);
	void LogNormalDist(double *, double *, double*);
	void DeltaDist(double *);
	void InitLatticeParam();
	void InitDiameterArray();
	int CalcPattern();
	void PolynomialLatticeParam();
	void ConstantLatticeParam();
	void CalcDistForShell(Distance *);
	int GetNumSRShells();
	void CheckSupportedShapes();
	void CheckShellThickness();
	int CalcDeltaShellSize();
	void CopyDistance(Distance *, Distance *);
	int GetNumSizes();
	int GetIndexOfShell(int);
};
#endif
