#ifndef POSITION_H
#define POSITION_H
#include "Includes.h"
#include <xmmintrin.h>

class Position
{
public:
	Position();
	Position(char*);
	~Position(void);

	string fileName;
	int shell, nAtoms;
	float *x, *y, *z;
	double f, kappa;
	Position *nextPos;

	int ReadPosFile();
	void ApplyPlanarSurfaceRelax(int, double, double);
	void ApplyRadialSurfaceRelax(int, double);
	void ScalePos(double);
	void ApplyNormExpRadialSurfRelax(int, double, double, double);
	void ApplyNormOscExpRadialSurfRelax(int, double, double, double, double);

};
#endif
