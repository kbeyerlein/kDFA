/*
 * Shape_PreCalcPatts.cpp
 *
 *  Created on: Jul 26, 2011
 *      Author: ken
 *      Description: Functions working with pre-calculated intensity files
 */

#include "Shape.h"

void Shape::FindPreCalcPatts()
{
	//Allocate the precalculated intensity array
	nPreCalcPatts=GetNumSizes();
	try{
		preCalcPatts=new Pattern *[nPreCalcPatts];
		for (int i=0;i<nPreCalcPatts;i++){
			preCalcPatts[i]=0;
		}
	}
	catch(exception &e){
		cout<<"Error: exception in LoadPatterns()"<<e.what()<<endl;
		exit(0);
	}

	// Read the precalculated patterns and compare with the existing positions
	ifstream pattfiles;
	string temp=pattPath+"intensityFiles.txt";
	pattfiles.open(temp.c_str());
	if(pattfiles.is_open()){
		printf("Reading intensity files for shape %s ...\n", name.c_str());
		while(!pattfiles.eof()){
			pattfiles >> temp;
			if (!pattfiles.eof()){
				int intShell=0;
				ifstream calcedPatt;
				//Read the shellSize from the first line of the intensity file
				calcedPatt.open(temp.c_str());
				if (calcedPatt.is_open()){
					string firstLine;
					getline(calcedPatt, firstLine);
					stringstream test(firstLine);
					while (!test.eof()){
						string stuff;
						test>>stuff;
						if (stuff=="shell:"){
							test>>intShell;
						}
					}
					calcedPatt.close();
				}
				else{
					cout<<"Warning: Unable to open intensity file: "<<temp<<endl;
				}
				if (intShell!=0){
					int intIndex=GetIndexOfShell(intShell);
					if (intIndex>=0){
						//Allocate and Load Intensity
						try{
							preCalcPatts[intIndex]=new Pattern [Pattern::nLambda];
						}
						catch(exception &e){
							cout<<"Error: exception in LoadPatterns()"<<e.what()<<endl;
							exit(1);
						}
						for (int i=0;i<Pattern::nLambda;i++){
							preCalcPatts[intIndex][i].fileName=temp;
							preCalcPatts[intIndex][i].shellIndex=intIndex;

							//preCalcPatts[intIndex][i].LoadIntFromFile(temp, i);
						}
					}
					else{
						cout<<"Warning: Positions of shell "<<intShell<<" are not loaded..."<<endl;
						cout<<"Ignoring intensity file: "<<temp<<endl;
					}
				}
				else{
					cout<<"Warning: Unable to determine shell size of intensity from file: "<<temp<<endl;
					cout<<"Check file format."<<endl;
				}
			}
		}
		pattfiles.close();
	}
	else{
		cout<<"Warning: Unable to open pre-calculated pattern file: "<< temp<<endl;
	}
}
void Shape::LoadPatterns()
{
	for (int i=0;i<nPreCalcPatts;i++){
		if (preCalcPatts[i]!=0){
			//NOTE: Need to have the lattice parameter determined for a given size before loading patterns.
			double aD=a_D[preCalcPatts[i][0].shellIndex];
			LoadIntFromFile(preCalcPatts[i], aD);
		}
	}
}
void Shape::CheckPreCalcPatts()
{
	//Check list of intensities and calculate necessary intensities
	bool stop=false;
	for (int i=0;i<nPreCalcPatts&&!stop;i++){
		if (preCalcPatts[i]==0){
			cout<<"Note: Not all shells of shape "<<name<<" have valid pre-calculated intensity files."<<endl;
			cout<<"To create them and speed up calculation use the command: 'function preCalcInt'"<<endl;
			stop=true;
		}
	}
}
void Shape::LoadIntFromFile(Pattern *patts, double aD)
{
	string lattice, shape;
	int shell, nI, nAtoms;
	double minSa, maxSa, dsa;
	double *inSa, *inI;
	ifstream calcedPatt;
	//Read the important info from the first two lines of the intensity file
	// Input format:
	// 			(lattice) (shape) (shell #) (# atoms)
	//			(minS*a) (maxS*a) (ds*a) (# Int points)
	//			(I_1)
	//			(I_2)
	//			(...)
	calcedPatt.open(patts[0].fileName.c_str());
	if (calcedPatt.is_open()){
		string firstLine;
		getline(calcedPatt, firstLine);
		stringstream test(firstLine);
		while (!test.eof()){
			string stuff;
			test>>lattice>>shape>>shell>>nAtoms;
			test>>minSa>>maxSa>>dsa>>nI;
		}

		//Allocate and Input the intensity from file

		try{
			inSa=new double [nI];
			inI=new double [nI];
			for (int i=0;i<nI;i++){
				inSa[i]=minSa+i*dsa;
				test>>inI[i];
			}
		}
		catch(exception &e){
			cout<<"Error: Exception caught in Pattern::LoadIntFromFile"<<e.what()<<endl;
			exit(1);
		}

		calcedPatt.close();

	}
	else{
		cout<<"Warning: Unable to open intensity file: "<<patts[0].fileName<<endl;
		return;
	}

	//Allocate and interpolate each intensity file.
	for (int i=0;i<Pattern::nLambda;i++){
		double *neededSa;
		try{
			if (patts[i].I){
				delete [] patts[i].I;
			}
			patts[i].I=new double [Pattern::nInt];
			neededSa =new double [Pattern::nInt];
			for (int j=0;j<Pattern::nInt;j++){
				neededSa[j]=Pattern::offsetq[i][j]*aD/(2*PI);
			}
		}
		catch(exception &e){
			cout<<"Error: Exception caught in Pattern::LoadIntFromFile"<<e.what()<<endl;
			exit(1);
		}
		// call to interpolation routine
		// Interpolate(inQ, inI, neededSa,patts[i].I);
		Interp(inSa, inI, nI,neededSa, patts[i].I, Pattern::nInt);
		delete [] neededSa;
	}

	delete [] inSa;
	delete [] inI;
}
