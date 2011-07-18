#ifndef PATTERN_H
#define PATTERN_H
#include "Includes.h"
#include "Param.h"
#include "Distance.h"

class Pattern
{
public:
	Pattern();
	Pattern(string);
	~Pattern(void);

	static int nInt, nChebPoly, nChebBG, nChebScale, nLambda, polQ, nGaussBG;
	static double min2Theta, max2Theta, *lambda, *wLambda, intStep, *theta2, *offset2Theta, **q, **offsetq, gonioRad;
	static double a1,b1,a2,b2,a3,b3,a4,b4,c;
	static Param *B, *scale, *sampleDisplace, *monoAngle2T;
	static Param **ChebBG, **gaussA, **gaussMu, **gaussSig, **ChebScale;
	static double **chebPolys, **gaussBG;
	static double ax,bx,cx,dx,ex;
	static Param *AExp, *kExp, *delExp;
	static bool synchThetaShift, multCos, divSin;
	static string geometry;
	static Param *t, *absCoef; //thickness (cm) and linear absorption coefficient (1/cm)

	string name, fileName;
	bool initBG;
	double *I, *stdDevSq, *residual, *background;
	int shellIndex;

	//int maxShell;

	void CalcPattern(Distance *, int , double *, double, Distance *, int);
	void CalcPattern(Distance *, int , double *, double *, int, Distance *, int);
	void CalcPattern(Distance *, int , double *, double *, int, Distance *, int, Pattern **);
	double Sinc(double);
	int OutputToFile(string, string);
	void ApplyDWFactor(int);
	void ApplyAtomicScatFactor(int);
	void ApplyPolarizationFactor();
	int ReadIntensityFile();
	void ReconfigPattBounds(int, double, double);
	void CalcChebyshevPolys();
	void CalcGaussBG();
	void ApplyChebyshevBG();
	void ApplyChebyshevScale();
	void ApplyGaussBG();
	void ApplySampleDisplacement();
	void ReconfigPattBoundsRaw();
	void ApplySynchThetaShift();
	void InitOffsetAngle();
	void ApplyExpBG();
	void InitBackground();
	void ApplyAbsorption();
	void ApplyWaveLenDependMultFactors(int);
	void ApplyGeomDependMultFactors();
	void ScalePattern(double);
	void AllocPattern();
	void ZeroPattern();
	void CheckGeometry();
	void MultCosTheta();
	void DivSinTheta();
	void LoadIntFromFile(string, int);
};
#endif
