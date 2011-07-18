/*
 * Param.h
 *
 *  Created on: Jul 12, 2011
 *      Author: ken
 */

#ifndef PARAM_H_
#define PARAM_H_
#include "Includes.h"
class Param {
private:
	double val;
public:
	bool global, fixed;
	double min, max;
	string name, shapeName;

	Param();
	virtual ~Param();
	void InitParam(double, double, double, string, string, bool, bool);
	void CheckParam();
	double GetVal();
	void ChangeVal(double);
	void CopyParam(Param *);
	void CopyParam(Param *, bool);
};

#endif /* PARAM_H_ */
