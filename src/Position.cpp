#include "Position.h"

Position::Position(char *_fileName)
{
	fileName=_fileName;
	shell=0;
	nAtoms=0;
	x=0; y=0; z=0;
	f=0; kappa=0;
	nextPos=0;
}

Position::Position()
{
	shell=0;
	nAtoms=0;
	x=0; y=0; z=0;
	f=0; kappa=0;
	nextPos=0;
}

Position::~Position(void)
{
	if (x!=0)
		_mm_free (x);
	if (y!=0)
		_mm_free (y);
	if (z!=0)
		_mm_free (z);
}

int Position::ReadPosFile()
{
	string tempstring1, tempstring2;
	int counter=1,i, err=0, temp;
	ifstream infile;
	infile.open(fileName.c_str());
	if (infile.is_open()){
		//printf("file %s opened\n", fileName.c_str());
		while ((!infile.eof())&& err==0)
		{
			infile >> tempstring1>>tempstring2;
			if (tempstring2=="SHAPE")
				infile >> tempstring1;
			//Input number of atoms in simulation
			else if (tempstring2=="NUMBER"){
				infile >> tempstring1>>tempstring1>>nAtoms;
				//allocate mem for positions and other info

				x=(float*) _mm_malloc(nAtoms*sizeof(float), 16);
				y=(float*) _mm_malloc(nAtoms*sizeof(float), 16);
				z=(float*) _mm_malloc(nAtoms*sizeof(float), 16);
				/*
				x=(float*) malloc(nAtoms*sizeof(float));
				y=(float*) malloc(nAtoms*sizeof(float));
				z=(float*) malloc(nAtoms*sizeof(float));
				*/
			}
			//Input atom positions
			else if (tempstring2=="ATOMS"){
				for (i=0; i<nAtoms; i++){
					infile >>temp>>shell>>x[i]>> y[i] >>z[i];

				}
			}
			else {
				//If there is an error with reading the data from the file an error message will be displayed and
				//	the loop index that the error occured on.
				if (counter <5){
					printf("Error in file: %s\n", fileName.c_str());
					printf("Error with data input: cycle %s\n",tempstring2.c_str());
					infile.close();
					err= -1;
				}
			}
			counter++;
		}
		infile.close();
		err=0;
	}
	else {
		printf("Could not open file: %s \n", fileName.c_str());
		err=-1;
	}
	return err;
}

void Position::ApplyPlanarSurfaceRelax(int size, double f, double kappa)
{
	double delta=f*exp(-(size-shell)/kappa)/shell;
	double scale=delta+1;
	ScalePos(scale);
}
void Position::ApplyNormExpRadialSurfRelax(int size, double f, double k, double delShell)
{
	double maxR=size*delShell;
	// So f*delShell will give the displacement of the largest shell in units of shell thickness
	f*=delShell;
	// kappa is  relative to the particle size, so kappa*maxR gives real units
	// Could subtract zero to make sure the function goes to zero as r->0
	// However, this leads to a correlation between f and kappa for magnitude of the displacement
	double zero=f*exp(-size/k);
	for (int i=0;i<nAtoms;i++){
		double ratio=0;
		double r=sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
		if (r!=0){
			ratio=f*exp(-(maxR-r)/(k*maxR));
			//ratio-=zero;
			ratio/=r;
		}
		// in the next comment r> means the vector r
		// r'> = r> + delta>; delta> = r> delta/r; r'> = r>(1+delta/r)
		double scale=(1+ratio);
		x[i]*=scale;
		y[i]*=scale;
		z[i]*=scale;
	}
}
void Position::ApplyNormOscExpRadialSurfRelax(int size, double f, double k, double n, double delShell)
{
	double maxR=size*delShell;
	// So f*delShell will give the displacement of the largest shell in units of shell thickness
	f*=delShell;
	// kappa is  relative to the particle size, so kappa*maxR gives real units
	// Could subtract zero to make sure the function goes to zero as r->0
	// However, this leads to a correlation between f and kappa for magnitude of the displacement

	// Define constant which is used in the oscillating part.
	// This is then 2Pi*k where k=1/lam
	// so lam=n/R where n is the number of oscillations in the radial direction
	double p=2*PI*n/maxR;
	for (int i=0;i<nAtoms;i++){
		double ratio=0;
		double r=sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
		if (r!=0){
			ratio=f*sin(p*r)*exp(-(maxR-r)/(k*maxR));
			//ratio-=zero;
			ratio/=r;
		}
		// in the next comment r> means the vector r
		// r'> = r> + delta>; delta> = r> delta/r; r'> = r>(1+delta/r)
		double scale=(1+ratio);
		x[i]*=scale;
		y[i]*=scale;
		z[i]*=scale;
	}
}
void Position::ApplyRadialSurfaceRelax(int size, double shellRad)
{
	float delta, r, maxR=size*shellRad;

	for (int i=0;i<nAtoms;i++){
		r=sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
		if (r!=0)
			delta=f*shellRad*exp(-(maxR-r)/(kappa*shellRad))/r;
		else delta=0;
		x[i]=x[i]*(1+delta);
		y[i]=y[i]*(1+delta);
		z[i]=z[i]*(1+delta);
	}
}

void Position::ScalePos(double scale)
{
	for (int i=0;i<nAtoms;i++){
		x[i]*=scale;
		y[i]*=scale;
		z[i]*=scale;
	}
}
