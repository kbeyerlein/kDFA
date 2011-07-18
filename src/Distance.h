
#ifndef DISTANCE_H
#define DISTANCE_H

#include "Includes.h"
#include "Position.h"

#include <emmintrin.h>  // for SSE intrinsics and helper functions
#include <xmmintrin.h>

class Distance
{
public:
	Distance(char*);
	Distance();
	~Distance(void);

	string fileName;
	int shell, nDist;
	double *dist, *mult, f, kappa;
	Distance *nextDist;

	void ReadDistFile();
	int CalcDistances(int, Position*);
	void QuickSortImproved(double *x, double *y, int lb, int ub);
	void sortByInsertion(double *x, double *y, int lb, int ub);
	int partition(double *x, double *y, int lb, int ub);
	int MultiplicitiesCompute(double *d, double *m, int nD);
	void OutputDistToFile(string, string);
	void DebbyCalcSums(float*, float*, float*, int, float*, float*, float*, int, double*, double, int, bool, float);
	int CalcDistLatticeSym(int, Position*);
	void ConsolidateDistance();
	void InitDistance(int);
	void ScaleDistance(double);
};

#endif
