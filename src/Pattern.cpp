#include "Pattern.h"

Pattern::Pattern()
{
	I=0;
	stdDevSq=0;
	residual=0;
	background=0;
	name="";
	initBG=false;
	fileName="";
	shellIndex=-1;
}
Pattern::Pattern(string _name)
{
	I=0;
	stdDevSq=0;
	residual=0;
	background=0;
	name=_name;
	initBG=false;
	fileName="";
	shellIndex=-1;
}

Pattern::~Pattern(void)
{
	delete [] stdDevSq;
	delete [] residual;
	delete [] I;
	delete [] background;
}

void Pattern::CalcPattern(Distance *curDist, int maxShell, double *wShell, double a, Distance *SRDist, int lambdaIndex)
{
	//Dynamic
	try{
		if (I==0)
			I= new double[nInt];
	}
	catch (exception& e){
		printf("\nException found in CalcPattern: %s\n", e.what());
		exit(1);
	}
	if (!curDist){
		cout<<"Error: Distance passed to Pattern->CalcPattern() is invalid\n";
		exit(0);
	}
	printf("Calculating pattern for %s\n", name.c_str());
	for (int i=0;i<nInt;i++)
		I[i]=0;
	//Intensity from the normal distances
	for (int i=0;(i<maxShell)&&(wShell[i]!=0);i++){
		for(int j=0;j<nInt;j++){
			for (int k=1;k<curDist->nDist;k++){
					I[j]+=curDist->mult[k]*wShell[i]*Sinc(offsetq[lambdaIndex][j]*a*curDist->dist[k]);
			}
			I[j]+=wShell[i]*curDist->mult[0];
		}
		curDist=curDist->nextDist;
	}
	// Add intensity from the surface relax distance list
	if (SRDist!=0){
		curDist=SRDist;
		for(int j=0;j<nInt;j++){
			for (int k=1;k<curDist->nDist;k++){
				//Distances already scaled by lattice param during calculation
				I[j]+=curDist->mult[k]*Sinc(offsetq[lambdaIndex][j]*curDist->dist[k]);
			}
			I[j]+=curDist->mult[0];
		}
	}
}

void Pattern::CalcPattern(Distance *curDist, int maxShell, double *p, double *a_D, int nSRShells, Distance *SRDist, int lambdaIndex)
{
	Distance *firstDist=0;
	if (!I){
		try{
			I= new double[nInt];
		}
		catch (exception& e){
			printf("\nException found in CalcPattern: %s\n", e.what());
			exit(1);
		}
	}
	if (!curDist){
		cout<<"Error: Distance passed to Pattern->CalcPattern() is invalid\n";
		exit(0);
	}
	else {
		firstDist=curDist;
	}
	printf("Calculating pattern for %s\n", name.c_str());
	for (int i=0;i<nInt;i++)
		I[i]=0;
	//Intensity from the normal distances
	for(int i=0; i<maxShell; i++){
		for (int j=0;j<i-nSRShells; j++){
			for (int k=0; k<nInt; k++){
				for (int l=1; l<curDist->nDist;l++){
					I[k]+=p[i]*curDist->mult[l]*Sinc(offsetq[lambdaIndex][k]*a_D[i]*curDist->dist[l]);
				}
				I[k]+=p[i]*curDist->mult[0];
			}
			curDist=curDist->nextDist;
		}
		curDist=firstDist;
	}
	// Add intensity from the surface relax distance list
	if (SRDist!=0){
		curDist=SRDist;
		for(int j=0;j<nInt;j++){
			for (int k=1;k<curDist->nDist;k++){
				//Distances already scaled by lattice parameter and weighted with size weight during SR calculation.
				I[j]+=curDist->mult[k]*Sinc(offsetq[lambdaIndex][j]*curDist->dist[k]);
			}
			I[j]+=curDist->mult[0];
		}
	}
}
void Pattern::CalcPattern(Distance *curDist, int maxShell, double *p, double *a_D, int nSRShells, Distance *SRDist, int lambdaIndex, Pattern **preCalcPatts)
{
	Distance *firstDist=0;
	if (!I){
		try{
			I= new double[nInt];
		}
		catch (exception& e){
			printf("\nException found in CalcPattern: %s\n", e.what());
			exit(1);
		}
	}
	if (!curDist){
		cout<<"Error: Distance passed to Pattern->CalcPattern() is invalid\n";
		exit(0);
	}
	else {
		firstDist=curDist;
	}
	printf("Calculating pattern for %s\n", name.c_str());
	for (int i=0;i<nInt;i++)
		I[i]=0;
	//Intensity from the precalculated Intensities or normal distances
	for(int i=0; i<maxShell; i++){
		if (preCalcPatts[i]==0){
			for (int j=0;j<i-nSRShells; j++){
				for (int k=0; k<nInt; k++){
					for (int l=1; l<curDist->nDist;l++){
						I[k]+=p[i]*curDist->mult[l]*Sinc(offsetq[lambdaIndex][k]*a_D[i]*curDist->dist[l]);
					}
					I[k]+=p[i]*curDist->mult[0];
				}
				curDist=curDist->nextDist;
			}
		}
		else{
			for (int k=0; k<nInt; k++){
				I[k]+=p[i]*preCalcPatts[i][lambdaIndex].I[k];
			}
		}
		curDist=firstDist;
	}
	// Add intensity from the surface relax distance list
	if (SRDist!=0){
		curDist=SRDist;
		for(int j=0;j<nInt;j++){
			for (int k=1;k<curDist->nDist;k++){
				//Distances already scaled by lattice parameter and weighted with size weight during SR calculation.
				I[j]+=curDist->mult[k]*Sinc(offsetq[lambdaIndex][j]*curDist->dist[k]);
			}
			I[j]+=curDist->mult[0];
		}
	}
}
double Pattern::Sinc(double x)
{
	if (x!=0){
		double sinc = sin(x)/x;
		return sinc;
	}
	else
		return 1.00;
}
int Pattern::OutputToFile(string _fPath, string _fName)
{
	int i,err=0;
	string filename=_fPath+_fName;
	fstream file(filename.c_str(), ios_base::out|ios_base::app);
	file.setf(ios::fixed, ios::floatfield);
	if (file.is_open()){
		printf("Intensity file exists.. Appending Intensity Data\n");
		file<<"\n2Theta Intensity ";
		if (residual!=0){
			file<<"Difference ";
		}

		if (background!=0){
			file<<"Background";
		}
		file<<endl;
		for (i=0;i<nInt;i++){
			file.precision(5);
			file <<theta2[i]<< " "<< I[i]<<" ";
			if (residual!=0){
				file<<residual[i]<<" ";
			}
			if (background!=0){
				file<<background[i];
			}
			file<<endl;
		}
		file.close();
		printf("\nIntensity written to file: %s\n",_fName.c_str());
	}
	else{
		printf("\nError: could not write to file: %s\n",filename.c_str());
		err=-1;
	}
	return err;
}
void Pattern::ApplyDWFactor(int lambdaIndex)
{
	double k;
	if (B->GetVal()!=0){
		for (int i=0;i<nInt;i++){
			k=offsetq[lambdaIndex][i]/(2*PI);
			I[i]*=exp(-B->GetVal()*k*k/2.0); // Is not correct for multiple wavelengths
		}

	}
	else printf("Not applying Debye Waller factor because B=0\n");
}
void Pattern::ApplyAtomicScatFactor(int lambdaIndex)
{
	double s, f;
	for (int i=0;i<nInt;i++){
		s=offsetq[lambdaIndex][i]/(4*PI);
		f=a1*exp(-b1*s*s)+a2*exp(-b2*s*s)+a3*exp(-b3*s*s)+a4*exp(-b4*s*s)+c;
		f=f*f;
		//TODO
		I[i]*=f;//Is not correct for multiple wavelengths
	}

	if (c==1) printf("Using a constant atomic scattering factor = 1\n");
}
void Pattern::ApplyPolarizationFactor()
{
	double cosMono, cos2Ang;;
	if (monoAngle2T!=0){
		cosMono=cos(monoAngle2T->GetVal()*PI/180);
	}
	else {
		cosMono=1.0;
	}
	for (int i=0; i<nInt;i++){
		cos2Ang=cos(offset2Theta[i]*PI/180);
		I[i]*=((1+polQ)+(1-polQ)*cosMono*cosMono*cos2Ang*cos2Ang)/(1+cosMono*cosMono);
	}
	//The correct form of Lorentz factor for powder diffraction pattern is 1/sin(theta)^2 (or 1/(sin(theta)*sin(theta_bragg)))
	//No Correction is needed for transmission (any difference is due to irradiated volume and is applied in the absorption term)
}

//Input Observed intensity file to be modeled
int Pattern::ReadIntensityFile()
{
	//TODO Clean-up ReadIntensityFile
	int fileType=-1, endHeader=0;
	int count=0, err=0, i;
	double minAngle=-1, maxAngle=0, temp, tempAngle;
	size_t found;
	ifstream infile;

	printf("Reading intensity file: %s\n", name.c_str());

	found=name.find(".debye");
	if (found!=string::npos) fileType=0;
	if (fileType==-1){
		found=name.find(".raw");
		if (found!=string::npos) fileType=1;
	}
	if (fileType==-1){
		found=name.find(".xye");
		if (found!=string::npos) fileType=2;
	}
	if(fileType==-1){
		found=name.find(".fit");
		if (found!=string::npos) fileType=3;
	}
	if (fileType==-1){
		found=name.find(".xy");
		if (found!=string::npos) fileType=4;
	}
	infile.open(name.c_str());
	if (infile.is_open()){
		if (fileType==0){
			string tempstring="";
			endHeader=infile.tellg();
			while (!infile.eof())
			{
				if (fileType==0&&endHeader==0){
					while(tempstring!="Intensity"&&!infile.eof()){
						infile>>tempstring;
					}
					endHeader=infile.tellg();
				}
				if (minAngle==-1){
					infile>>minAngle>>temp;
					maxAngle=minAngle;
					count++;
				}
				else{
					infile>>tempAngle>>temp;
					if (tempAngle>maxAngle){
						maxAngle=tempAngle;
						count++;
					}
				}
			}
			if (count>1){
		//Currently not accurate when trying to check if input intensity params = calc intensity params
				if (count!=nInt||minAngle!=min2Theta||(maxAngle!=max2Theta)){
					cout<<"Setting nInt: "<<count<<" min2theta: "<<minAngle<<" max2theta: "<<maxAngle<<endl;
					ReconfigPattBounds(count, minAngle, maxAngle);
				}
				try{
					if (I!=0) {
						delete [] I;
					}
					if (stdDevSq!=0){
						delete [] stdDevSq;
					}
					I= new double [Pattern::nInt];
					stdDevSq =new double [Pattern::nInt];

				}
				catch (exception& e){
					printf("\nException found in ReadIntensityFile: %s\n", e.what());
					return -1;
				}
				infile.clear();

				infile.seekg(endHeader, ios::beg);
				for (int i=0;i<Pattern::nInt;i++) {
					infile>>temp>>I[i];
					stdDevSq[i]=I[i];
				}
			}
			else{
				printf("\nError: Cannot find Intensity data\n");
				infile.close();
				return -1;
			}
		}
		else if(fileType==1){
			infile.ignore(256,'\n');
			/*Pattern::nLambda=1;
			if (Pattern::lambda){
				delete [] Pattern::lambda;
			}
			if (Pattern::wLambda){
				delete [] Pattern::wLambda;
			}
			Pattern::lambda=new double [Pattern::nLambda];
			Pattern::wLambda=new double [Pattern::nLambda];
			Pattern::wLambda[0]=1;
			infile>>Pattern::nInt>>Pattern::intStep>>Pattern::min2Theta>>Pattern::lambda[0]>>tempAngle;
			*/
			infile>>Pattern::nInt>>Pattern::intStep>>Pattern::min2Theta;
			infile.ignore(256,'\n');
			ReconfigPattBoundsRaw();
			offset2Theta=theta2;
			offsetq=q;

			try{
				if (I!=0) {
					delete [] I;
				}
				if (stdDevSq!=0){
					delete [] stdDevSq;
				}
				I= new double [Pattern::nInt];
				stdDevSq =new double [Pattern::nInt];

			}
			catch (exception& e){
				printf("\nException found in ReadIntensityFile: %s\n", e.what());
				return -1;
			}
			for(int i=0;i<Pattern::nInt;i++){
				infile>>I[i];
				stdDevSq[i]=I[i];
			}
		}
		else if (fileType==2){
			while ((!infile.eof())&& err==0)
			{
				if (minAngle==-1){
					infile>>minAngle>>temp>>temp;
					maxAngle=minAngle;
					count++;
				}
				else{
					infile>>tempAngle>>temp>>temp;
					if (tempAngle>maxAngle){
						maxAngle=tempAngle;
						count++;
					}
				}
			}
			infile.close();
			Pattern::nInt=count;
			Pattern::min2Theta=minAngle;
			Pattern::max2Theta=maxAngle;
			try{
				if (I!=0) {
					delete [] I;
				}
				if (stdDevSq!=0){
					delete [] stdDevSq;
				}
				stdDevSq=new double [Pattern::nInt];
				I= new double [Pattern::nInt];
				Pattern::offset2Theta=0;
				if (Pattern::theta2!=0){
					delete [] Pattern::theta2;
				}
				Pattern::theta2=new double [Pattern::nInt];
				Pattern::offset2Theta=Pattern::theta2;
				Pattern::offsetq=0;
				if (Pattern::q!=0){
					for (int k=0; k<Pattern::nLambda; k++){
						delete [] Pattern::q[k];
					}
				}
				for (int k=0;k<Pattern::nLambda; k++){
					Pattern::q[k]=new double [Pattern::nInt];
				}
				Pattern::offsetq=Pattern::q;
			}
			catch (exception& e){
				printf("\nException found in ReadIntensityFile: %s\n", e.what());
				infile.close();
				return -1;
			}
			infile.clear();
			infile.open(name.c_str());

			for (i=0;i<Pattern::nInt;i++) {
				infile>>Pattern::theta2[i]>>I[i]>>temp;
				stdDevSq[i]=temp*temp;
				for (int j=0;j<Pattern::nLambda;j++){
					Pattern::q[j][i]=4*PI*sin(Pattern::theta2[i]*PI/360.0f)/Pattern::lambda[j];
				}
			}
			Pattern::intStep=Pattern::theta2[1]-Pattern::theta2[0];

		}
		else if(fileType==3){
			string tempstring="";
			endHeader=infile.tellg();
			while (!infile.eof())
			{
				if (fileType==3&&endHeader==0){
					while(tempstring!="Background"&&!infile.eof()){
						infile>>tempstring;
					}
					endHeader=infile.tellg();
				}
				if (minAngle==-1){
					infile>>minAngle>>temp>>temp;
					maxAngle=minAngle;
					count++;
				}
				else{
					infile>>tempAngle>>temp>>temp;
					if (tempAngle>maxAngle){
						maxAngle=tempAngle;
						count++;
					}
				}
			}
			if (count>1){
		//Currently not accurate when trying to check if input intensity params = calc intensity params
				/*if (count!=nInt||minAngle!=min2Theta||(maxAngle!=max2Theta-intStep))
					err=ReconfigPattBounds(count, minAngle, maxAngle);
					*/
				try{
					if (I!=0) {
						delete [] I;
					}
					if (stdDevSq!=0){
						delete [] stdDevSq;
					}
					I= new double [Pattern::nInt];
					stdDevSq =new double [Pattern::nInt];

				}
				catch (exception& e){
					printf("\nException found in ReadIntensityFile: %s\n", e.what());
					return -1;
				}
				infile.clear();

				infile.seekg(endHeader, ios::beg);
				for (i=0;i<Pattern::nInt;i++) {
					infile>>temp>>I[i]>>temp;
					stdDevSq[i]=I[i];
				}
			}
			else{
				printf("\nError: Cannot find Intensity data\n");
				infile.close();
				return -1;
			}
		}
		else if (fileType==4){
			while ((!infile.eof())&& err==0)
			{
				if (minAngle==-1){
					infile>>minAngle>>temp;
					maxAngle=minAngle;
					count++;
				}
				else{
					infile>>tempAngle>>temp;
					if (tempAngle>maxAngle){
						maxAngle=tempAngle;
						count++;
					}
				}
			}
			infile.close();
			Pattern::nInt=count;
			Pattern::min2Theta=minAngle;
			Pattern::max2Theta=maxAngle;
			try{
				if (I!=0) {
					delete [] I;
				}
				if (stdDevSq!=0){
					delete [] stdDevSq;
				}
				stdDevSq=new double [Pattern::nInt];
				I= new double [Pattern::nInt];
				Pattern::offset2Theta=0;
				if (Pattern::theta2!=0){
					delete [] Pattern::theta2;
				}
				Pattern::theta2=new double [Pattern::nInt];
				Pattern::offset2Theta=Pattern::theta2;
				Pattern::offsetq=0;
				if (Pattern::q!=0){
					for (int k=0; k<Pattern::nLambda; k++){
						delete [] Pattern::q[k];
					}
				}
				for (int k=0;k<Pattern::nLambda; k++){
					Pattern::q[k]=new double [Pattern::nInt];
				}
				Pattern::offsetq=Pattern::q;
			}
			catch (exception& e){
				printf("\nException found in ReadIntensityFile: %s\n", e.what());
				infile.close();
				return -1;
			}
			infile.clear();
			infile.open(name.c_str());

			for (int i=0;i<Pattern::nInt;i++) {
				infile>>Pattern::theta2[i]>>I[i];
				stdDevSq[i]=I[i];
				for (int j=0;j<Pattern::nLambda;j++){
					Pattern::q[j][i]=4*PI*sin(Pattern::theta2[i]*PI/360.0f)/Pattern::lambda[j];
				}
			}
			Pattern::intStep=Pattern::theta2[1]-Pattern::theta2[0];

		}
		else{
			printf("\nError: Unknown observed intensity file type\n");
			infile.close();
			return -1;
		}
		infile.close();
	}
	else{
		printf("\nError: Could not open file %s\n",name.c_str());
		return -1;
	}
	printf("Intensity file successfully input\n");
	return 0;
}
void Pattern::ReconfigPattBounds(int _nInt, double _minTheta, double _maxTheta)
{
	Pattern::min2Theta=_minTheta;
	Pattern::max2Theta=_maxTheta;
	Pattern::nInt=_nInt;
	Pattern::intStep=(Pattern::max2Theta-Pattern::min2Theta)/(Pattern::nInt-1);

	try{
		if (Pattern::theta2)
			delete [] Pattern::theta2;
		if (Pattern::q){
			for(int i=0;i<nLambda;i++){
				if (q[i]){
					delete [] q[i];
				}
			}
			delete [] q;
		}
		Pattern::theta2= new double [Pattern::nInt];
		Pattern::q=new double*[Pattern::nLambda];
		for (int i=0;i<Pattern::nLambda;i++)
			Pattern::q[i]=new double [Pattern::nInt];
	}
	catch (exception& e){
		printf("\nException found in ReconfigPattBounds: %s\n", e.what());
		exit(1);
	}

	for(int i=0;i<Pattern::nInt;i++){
		Pattern::theta2[i]=Pattern::min2Theta+i*Pattern::intStep;
		for (int j=0;j<Pattern::nLambda; j++)
			Pattern::q[j][i]=4*PI*sin(Pattern::theta2[i]*PI/360.0f)/Pattern::lambda[j];
	}
	offset2Theta=theta2;
	offsetq=q;
}
void Pattern::CalcChebyshevPolys()
{
	try{
		chebPolys = new double* [nChebPoly];
		for (int i=0;i<nChebPoly; i++){
			chebPolys[i]= new double [nInt];
		}
	}
	catch (exception& e){
		printf("\nException found in CalcChebyshevPolys: %s\n", e.what());
		exit(1);
	}
	double shift=.5*(min2Theta+max2Theta);
	double scale=.5*(max2Theta-min2Theta);
	//Specify the first two and use the recursion relation to calculate higher orders with 2Theta in place of x
	for(int i=0;i<nInt;i++){
		chebPolys[0][i]=1;
		double x=(theta2[i]-shift)/scale;
		if (nChebPoly>1) chebPolys[1][i]=x;
		for (int j=2; j<nChebPoly; j++){
			chebPolys[j][i]=2*x*chebPolys[j-1][i]-chebPolys[j-2][i];
		}
	}
}
void Pattern::CalcGaussBG()
{
	try{
		if (nGaussBG!=0 && gaussBG==0){
			gaussBG= new double* [nGaussBG];
			for (int i=0;i<nGaussBG; i++){
				gaussBG[i]= new double [nInt];
			}
		}
	}
	catch (exception& e){
		printf("\nException found in CalcGaussBG: %s\n", e.what());
		exit(1);
	}
	for (int i=0; i<nInt; i++){
		for(int j=0;j<nGaussBG; j++){
			gaussBG[j][i]=gaussA[j]->GetVal()*exp(-(theta2[i]-gaussMu[j]->GetVal())*(theta2[i]-gaussMu[j]->GetVal())/(2*gaussSig[j]->GetVal()*gaussSig[j]->GetVal()));
		}
	}
}
void Pattern::ApplyChebyshevBG()
{
	bool negI=false;
	if (!initBG){
		InitBackground();
	}
	for (int i=0; i<nInt; i++){
		if (I[i]<0) negI=true;
		for (int j=0;j<nChebBG; j++){
			I[i]+=chebPolys[j][i]*ChebBG[j]->GetVal();
			background[i]+=chebPolys[j][i]*ChebBG[j]->GetVal();
		}
	}
	if (negI){
		printf("\nWarning: Negative Intensity in Pattern::ApplyChebyshevBG\n");
	}
}
void Pattern::ApplyChebyshevScale()
{
	double temp;
	if (chebPolys==0&&nChebScale>0) CalcChebyshevPolys();
	for (int i=0; i<nInt; i++){
		temp=scale->GetVal();
		for (int j=0;j<nChebScale; j++){
			temp+=chebPolys[j+1][i]*ChebScale[j]->GetVal();
		}
		I[i]*=temp;
	}
}
void Pattern::ApplyGaussBG()
{
	bool negI=false;
	if (!initBG){
		InitBackground();
	}
	for (int i=0; i<nInt; i++){
		if (I[i]<0) negI=true;
		for (int j=0;j<nGaussBG; j++){
			I[i]+=gaussBG[j][i];
			background[i]+=gaussBG[j][i];
		}
	}
	if (negI){
		printf("\nWarning: Negative Intensity in Pattern::ApplyGaussBG\n");
	}
}
void Pattern::ApplyExpBG()
{
	double tempExp;
	double del=0;
	bool negIb4=false, negIaftr=false;
	if (delExp){
		del=delExp->GetVal();
	}
	if (!initBG){
		InitBackground();
	}
	for (int i=0; i<nInt; i++){
		if (I[i]<0) negIb4=true;
		tempExp=AExp->GetVal()*exp(-(theta2[i]-del)/kExp->GetVal());
		I[i]+=tempExp;
		background[i]+=tempExp;
		if (I[i]<0) negIaftr=true;
	}
	if (negIb4){
		cout<<"\nWarning: Negative Intensity BEFORE application of Exp BG! \n";
	}
	if (negIaftr){
		cout<<"\nWarning: Negative Intensity AFTER application of Exp BG! \n";
	}

}

void Pattern::ApplySampleDisplacement()
{
	if (offset2Theta==theta2 && offsetq==q){
		InitOffsetAngle();
	}
	if (gonioRad!=0){
		synchThetaShift=false;
		for (int i=0; i<nInt;i++){
			// initialize offset2theta
			offset2Theta[i]=theta2[i];
			// Apply angular shift, is in units of degrees (without 180/PI it is in units of radians)
			offset2Theta[i]-=2*180*sampleDisplace->GetVal()*cos(theta2[i]*PI/360)/(gonioRad*PI);
			for (int j=0; j<nLambda; j++)
				offsetq[j][i]=4*PI*sin(offset2Theta[i]*PI/360.0)/lambda[j];
		}
	}
	else {
		printf("\nError: Goniometer Radius set to zero, cannot apply Sample Displacement correction.\n");
		exit(0);
	}
}
void Pattern::ApplySynchThetaShift()
{
	double tanTheta;
	if (offset2Theta==theta2 && offsetq==q){
		InitOffsetAngle();
	}
	synchThetaShift=true;
	for (int i=0; i<nInt;i++){
		//Not sure if shift is in correct direction, tested it both ways and (-) seemed to be the most reasonable
		tanTheta=tan(offset2Theta[i]*PI/360.0);
		offset2Theta[i]-=ax/tanTheta+bx+cx*tanTheta+dx*tanTheta*tanTheta+ex*tanTheta*tanTheta*tanTheta;
		for (int j=0; j<nLambda; j++)
			offsetq[j][i]=4*PI*sin(offset2Theta[i]*PI/360.0)/lambda[j];
	}
}
void Pattern::InitOffsetAngle()
{
	try{
		offset2Theta= new double [nInt];
		offsetq= new double* [nLambda];
		for (int i=0;i<nLambda;i++)
			offsetq[i]=new double [nInt];
	}
	catch (exception& e){
		printf("\nException found in ApplySynchThetaDisplacement: %s\n", e.what());
		exit(1);
	}
	for (int i=0;i<nInt;i++){
		offset2Theta[i]=theta2[i];
	}
}
void Pattern::ReconfigPattBoundsRaw()
{
	Pattern::max2Theta=min2Theta+(nInt-1)*intStep;

	try{
		if (Pattern::theta2!=0){
			delete [] Pattern::theta2;
		}

		if (Pattern::q!=0){
			delete [] Pattern::q;
		}
		Pattern::theta2= new double [Pattern::nInt];
		Pattern::q=new double* [Pattern::nLambda];
		for (int i=0;i<Pattern::nLambda;i++)
			Pattern::q[i]=new double [Pattern::nInt];
	}
	catch (exception& e){
		printf("\nException found in ReconfigPattBounds: %s\n", e.what());
		exit(1);
	}

	for(int i=0;i<Pattern::nInt;i++){
		Pattern::theta2[i]=Pattern::min2Theta+i*Pattern::intStep;
		for (int j=0;j<Pattern::nLambda;j++)
			Pattern::q[j][i]=4*PI*sin(Pattern::theta2[i]*PI/360.0f)/Pattern::lambda[j];
	}
}
void Pattern::InitBackground()
{
	if (background==0){
		try{
			background=new double [nInt];
		}
		catch (exception &e){
			cout<<"Exception in ApplyChebyshevBG: "<<e.what();
			exit(1);
		}
	}
	for (int i=0;i<nInt;i++){
		background[i]=0;
	}
	initBG=true;
}
void Pattern::ApplyAbsorption()
{
	//Absorption is multiplied by a scaling term (sc) so that as the values of the absorption coefficients,
	//are changed, the actually scale of the pattern is minimally influenced.
	//This decreases the dependance of the scale parameter on the absorption params.

	//double sc=exp(absCoef->GetVal()*t->GetVal());
	if (geometry=="transmission"){
		for (int i=0;i<nInt;i++){
			//I[i]*=sc*exp(-absCoef->GetVal()*t->GetVal()/cos(offset2Theta[i]*PI/360.0))/cos(offset2Theta[i]*PI/360);
			I[i]*=t->GetVal()*exp(-absCoef->GetVal()*t->GetVal()/cos(offset2Theta[i]*PI/360.0))/cos(offset2Theta[i]*PI/360);
		}
	}
	/*
	else if(geometry=="transmission_bathed"){
		for (int i=0;i<nInt;i++){
			//I[i]*=sc*exp(-absCoef->GetVal()*t->GetVal()/cos(offset2Theta[i]*PI/360.0));
			I[i]*=sc*exp(-absCoef->GetVal()*t->GetVal()/cos(offset2Theta[i]*PI/360.0));
		}
	}
	*/
	else if(geometry=="reflection_thinFilm"){
		for (int i=0;i<nInt; i++){
			//I[i]*=(1-exp(-absCoef->GetVal()*t->GetVal()/sin(offset2Theta[i]*PI/360.0)));
			I[i]*=(1-exp(-absCoef->GetVal()*t->GetVal()/sin(offset2Theta[i]*PI/360.0)))/absCoef->GetVal();

		}
	}
	else if(geometry=="reflection"){
		for (int i=0;i<nInt;i++){
			//I[i]*=(1-exp(-2*absCoef->GetVal()*t->GetVal()/sin(offset2Theta[i]*PI/360.0)));
			I[i]*=(1-exp(-2*absCoef->GetVal()*t->GetVal()/sin(offset2Theta[i]*PI/360.0)))/(2*absCoef->GetVal());
		}
	}
	else{
		cout<<"Geometry: "<<geometry<<" is not supported in Absorption Correction."<<endl;
		exit(0);
	}
}
void Pattern::ApplyWaveLenDependMultFactors(int j)
{
	ApplyAtomicScatFactor(j);
	ApplyDWFactor(j);
}
void Pattern::ApplyGeomDependMultFactors()
{
	ApplyPolarizationFactor();
	if (t&&absCoef){
		ApplyAbsorption();
	}
	if(multCos){
		MultCosTheta();
	}
	if (divSin){
		DivSinTheta();
	}
}
void Pattern::ScalePattern(double scale)
{
	for (int i=0;i<nInt;i++){
		I[i]*=scale;
	}
}
void Pattern::AllocPattern()
{
	try{
		if (I){
			delete [] I;
		}
		I=new double [nInt];
		if (background){
			delete [] background;
		}
		background =new double [nInt];
		if (residual){
			delete [] residual;
		}
		residual = new double [nInt];
	}
	catch(exception &e){
		cout<<"Exception in Pattern::AllocI : "<<e.what()<<endl;
		exit(1);
	}
	ZeroPattern();
}
void Pattern::ZeroPattern()
{
	for (int i=0;i<nInt;i++){
		I[i]=0;
		background[i]=0;
		residual[i]=0;
	}
}
void Pattern::CheckGeometry()
{
	if (Pattern::geometry=="reflection")
		;
	else if(Pattern::geometry=="reflection_thinFilm")
		;
	else if(Pattern::geometry=="transmission")
		;
	else if(Pattern::geometry=="transmission_bathed")
		;
	else{
		cout<<endl<<"Error: Specified unsupported measurement geometry: "<<Pattern::geometry<<endl;
	}
}
void Pattern::MultCosTheta()
{
	for (int i=0;i<nInt;i++){
		I[i]*=cos(offset2Theta[i]*PI/360.0);
	}
}
void Pattern::DivSinTheta()
{
	for (int i=0;i<nInt;i++){
		I[i]/=sin(offset2Theta[i]*PI/360.0);
	}
}

