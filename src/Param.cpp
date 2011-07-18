/*
 * Param.cpp
 *
 *  Created on: Jul 12, 2011
 *      Author: ken
 */

#include "Param.h"

Param::Param() {
	val=0;
	min=0;
	max=0;
	fixed=true;
	global=false;
	name="";
	shapeName="";
}

Param::~Param() {
	// TODO Auto-generated destructor stub
}
void Param::InitParam(double _v, double _min, double _max, string _name, string _sName, bool _fixed, bool _global)
{
	val=_v;
	min=_min;
	max=_max;
	name=_name;
	shapeName=_sName;
	fixed=_fixed;
	global=_global;
}
void Param ::CheckParam()
{
	if (min>max){
		cout<<"Error in input of parameter: "<<name<<", min>max."<<endl;
		double tempmin=min;
		min=getmin(min,max);
		max=getmax(tempmin,max);
		cout<<"Setting min: "<<min<<", setting max: "<<max<<endl;
	}
	if (min==max){
		cout<<"Error in input of parameter: "<<name<<", min=max."<<endl;
		exit(0);
	}
	ChangeVal(val);
}
void Param::ChangeVal(double newVal)
{
	if (newVal>max||newVal<min){
		cout<<"New Value of parameter: "<<name<<" is outside of bounds."<<endl;
	}
	val=getmin(newVal,max);
	val=getmax(newVal,min);
	if (val!=newVal){
		cout<<"Setting value to "<<val<<endl;
	}
}
void Param::CopyParam(Param *source)
{
	val=source->val;
	min=source->min;
	max=source->max;
	name=source->name;
	shapeName=source->shapeName;
	fixed=source->fixed;
	global=source->global;
}
void Param::CopyParam(Param *source, bool glob)
{
	val=source->val;
	min=source->min;
	max=source->max;
	name=source->name;
	shapeName=source->shapeName;
	fixed=source->fixed;
	global=glob;
}
double Param::GetVal()
{
	return val;
}
