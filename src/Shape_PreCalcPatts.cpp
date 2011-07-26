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
	for (int i=0;i<nPreCalcPatts;i++){
		if (preCalcPatts[i]==0){
			//Call routine to calculate the intensity
			// The x-axis must be in terms of q*a (dimensionless)
		}
	}
}
