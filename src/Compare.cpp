#include "Compare.h"

Compare::Compare(Shape *_firstShape, int _nParam)
{
	nParam=_nParam;
	firstShape=_firstShape;
	indParams=0;
	initCalcI=0;
	finalCalcI=0;
	obsI = 0;
	derivfSize_mu=0;
	derivfSize_sig=0;
	gradI=0;
	A=0;
	g=0;
	delta=0;
	name="";
	path="";
	nAtoms=1;
}
Compare::Compare(string _path, string _name, Shape *_firstShape, int _nParam)
{
	path=_path;
	name=_name;
	nParam=_nParam;
	firstShape=_firstShape;
	indParams=0;
	initCalcI=0;
	finalCalcI=0;
	obsI = 0;
	derivfSize_mu=0;
	derivfSize_sig=0;
	gradI=0;
	A=0;
	g=0;
	delta=0;
	nAtoms=1;
	covar=0;
}

Compare::~Compare(void)
{
	if (initCalcI){
		for (int i=0;i<Pattern::nLambda;i++){
			delete initCalcI[i];
		}
		delete [] initCalcI;
	}
	delete finalCalcI;
	delete obsI;
	delete [] indParams;
	for (int i=0;i<nParam; i++){
		if (gradI)
			delete [] gradI[i];
		if (A)
			delete [] A[i];
	}
	delete [] gradI;
	delete [] A;
	delete [] g;
	delete [] delta;
	delete [] derivfSize_mu;
	delete [] derivfSize_sig;
	if (covar){
		gsl_matrix_free(covar);
	}

//Dealloc all other shapes then firstShape to avoid deletion of global parameters
	Shape *curShape=firstShape->nextShape, *nextShape;
	while (curShape!=0){
		//Dealloc pos and dist info here so that temp shapes can be correctly deallocated
		Position *curPos=curShape->firstPos, *nextPos=0;
		while(curPos){
			nextPos=curPos->nextPos;
			delete curPos;
			curPos=nextPos;
		}
		Distance *curDist=curShape->firstDist, *nextDist=0;
		while(curDist){
			nextDist=curDist->nextDist;
			delete curDist;
			curDist=nextDist;
		}
		nextShape=curShape->nextShape;
		delete curShape;
		curShape=nextShape;
	}
	//Clean Global Parameters
	if (firstShape->a[0]->global){
		for (int i=0;i<firstShape->nParams_a;i++){
			delete firstShape->a[i];
		}
		delete [] firstShape->a;
		firstShape->a=0;
	}
	CleanGlobalParam(&firstShape->alpha);
	CleanGlobalParam(&firstShape->beta);
	CleanGlobalParam(&firstShape->mu);
	CleanGlobalParam(&firstShape->sigma);
	CleanGlobalParam(&firstShape->deltaSize);
	CleanGlobalParam(&firstShape->f);
	CleanGlobalParam(&firstShape->kappa);

	Position *curPos=firstShape->firstPos, *nextPos=0;
	while(curPos){
		nextPos=curPos->nextPos;
		delete curPos;
		curPos=nextPos;
	}
	Distance *curDist=firstShape->firstDist, *nextDist=0;
	while(curDist){
		nextDist=curDist->nextDist;
		delete curDist;
		curDist=nextDist;
	}
	delete firstShape;
}

//Levenberg-Marquardt routine
int Compare::LevMarq(int nIter)
{
	//Local Variables
	bool converged=false, e=false, outputFinalPatt=false;
	int err=0, nConvParam=0, nInitialParam=nParam, iter=0;
	float oldRwp=0, bestRexp=0, bestChiSq=0;
	double lambda=.1, gamma, tempsum, convlimit=.00001;
	double *saveParams=0;
	double **tempA;
	string outFileName=path+"levMarqOut."+name+".txt";
	ofstream outFile(outFileName.c_str(), ios::app);
	//Dynamic allocation
	try {
		saveParams= new double [nParam];
	}
	catch (exception& e){
		cout<<"\nException found in LevMarq alloc: "<< e.what()<<endl;
		return -1;
	}

	cout<<"\nStarting LevMarq routine with "<<nParam<<" free parameters.\n";
	outFile<<"Starting LevMarq routine with "<<nParam<<" free parameters.\n";

	//Making a copy of initial parameters
	for (int i=0;i<nParam;i++){
		saveParams[i]=indParams[i]->GetVal();
	}

	cout<<"Calculating Pattern from Initial Parameters\n";
	err=CalcPattern(false, true);

	if (!err){
		CalcRwp();

		cout<<"\nRwp equals: "<<rwp<<"\n";
		cout<<"Rexp equals: "<<rexp<<"\n";
		cout<<"ChiSq equals: "<<chisq<<"\n";
		outFile<<"\nRwp equals: "<<rwp<<"\n";
		outFile<<"Rexp equals: "<<rexp<<"\n";
		outFile<<"ChiSq equals: "<<chisq<<"\n";


		oldRwp=rwp;
		bestRexp=rexp;
		bestChiSq=chisq;
		for (iter=0; iter<nIter&&!err&&!converged; iter++){
			CalcResidual();
			if (!err) err=OutputPattern(path,name, iter, converged);

			cout<<"\n----------------------------\nStarting iteration "<<iter+1<<"\n";
			outFile<<"\n----------------------------\nStarting iteration "<<iter+1<<"\n";

			nParam=nInitialParam;
			e = false;
			lambda*=.1;			// Other suggested lambda step is lambda*=.31622777;
			//lambda*=.31622777;
			//Check if any shapes have been set to zero and then avoid calculating gradient for variables of those shapes
			/*curShape=firstShape;
			while (curShape!=0){
				if (curShape->wShape->val<.01){

					cout<<"\nWarning: Very low shape vol. weight, Consider refining excluding "<<curShape->name<<"\n";
					outFile<<"\nWarning: Very low shape vol. weight, Consider refining excluding "<<curShape->name<<"\n";

					for (int m=0; m<nParam;m++){
						if (indParams[m]->shapeName==curShape->name&&indParams[m]->name!=curShape->name+"_w"){
							Param *temp=indParams[nParam-1];
							double tempSave=saveParams[nParam-1];
							indParams[nParam-1]=indParams[m];
							saveParams[nParam-1]=saveParams[m];
							indParams[m]=temp;
							saveParams[m]=tempSave;
							nParam--;
							m--;
						}
					}
				}
				curShape=curShape->nextShape;
			}
			*/
			try{
				delta=new double [nParam];
				tempA=new double* [nParam];
				for (int m=0; m<nParam; m++)
					tempA[m]=new double [(m+1)];
			}
			catch (exception& e){
				cout<<"\nException found in LevMarq alloc: "<<e.what()<<"\n";
				return -1;
			}
			//Calc gradI vector wrt free parameters
			CalcGrad();
			//Calc gradI_i*gradI_j symmetric matrix
			CalcA();
			//Calc gradient of chisq wrt free parameters
			Calcg();
			//Save tempA that can be changed so to avoid recalculating the gradient
			for (int k=0;k<nParam;k++){
				for (int j=k; j<nParam;j++){
					tempA[j][k]=A[j][k];
				}
			}
			do{
				//Find delta to solve solution to (1+lam)A*delta=g where g is the gradient vector
				nConvParam=0;
				for (int j=0;j<nParam;j++) tempA[j][j]=(1+lambda)*A[j][j];
				// Decompose A into two triangular matrices (old algorithm)
				//err=SymDiagDecomp(tempA, nParam, g, delta);
				SolveForDelta(tempA, g, delta, nParam);

				// Unscale the delta which was found
				//for (j=0;j<nParam;j++) delta[j]/=sqrtDiagA[j];
				//gamma is the angle between the search vector and the gradient vector
				gamma=0;
				for (int j=0; j<nParam;j++) gamma+=delta[j]*g[j];
				tempsum=0;
				for (int j=0;j<nParam;j++) tempsum+=delta[j]*delta[j];
				if (tempsum!=0){
					gamma=gamma/sqrt(tempsum);
				}
				else{
					err=1;
					cout<< "Error in Lev-Marq: Undefined delta vector\n";
					exit(1);
				}
				tempsum=0;
				for (int j=0;j<nParam;j++) tempsum+=g[j]*g[j];
				if (tempsum!=0){
					gamma/=sqrt(tempsum);
				}
				else{
					err=1;
					cout<< "Error in Lev-Marq: Undefined g vector\n";
					exit(1);
				}
				gamma=acos(gamma);

				cout<<"\nLambda "<<lambda<<"\nGamma "<<gamma<<"\n\nDelta and New Params\n";
				outFile<<"\nLambda "<<lambda<<"\nGamma "<<gamma<<"\n\nDelta and New Params\n";

				for (int j=0; j<nParam&&!err; j++){
					//Calculating quantity used in the Marquardt (1967) paper to measure convergence
					if (fabs(delta[j])<convlimit*(fabs(gamma)+fabs(indParams[j]->GetVal()))) nConvParam++;
					//updating parameters
					indParams[j]->ChangeVal(saveParams[j]+delta[j]);
					cout <<"delta "<<indParams[j]->name<<": "<<delta[j]<<endl;
					cout<<indParams[j]->name<<": "<<indParams[j]->GetVal()<<endl;
					outFile <<"delta "<<indParams[j]->name<<": "<<delta[j]<<endl;
					outFile<<indParams[j]->name<<": "<<indParams[j]->GetVal()<<endl;
				}

				cout<<"Num Converged Parameters: "<<nConvParam<<"\n\n";
				outFile<<"Num Converged Parameters: "<<nConvParam<<"\n\n";
				NormalizeShapeWeights();

				//Calc new pattern and chisq with updated parameters
				cout<<"Calculating Pattern:\n";
				err=CalcPattern(false, true);
				if (!err){
					CalcRwp();

					cout<<"\nRwp equals: "<<rwp<<"\n";
					cout<<"Rexp equals: "<<rexp<<"\n";
					cout<<"ChiSq equals: "<<chisq<<"\n";
					outFile<<"\nRwp equals: "<<rwp<<"\n";
					outFile<<"Rexp equals: "<<rexp<<"\n";
					outFile<<"ChiSq equals: "<<chisq<<"\n";


					if (rwp<oldRwp){
						oldRwp=rwp;
						bestRexp=rexp;
						bestChiSq=chisq;
						for(int j=0;j<nParam;j++){
							saveParams[j]=indParams[j]->GetVal();
						}
						outputFinalPatt=true;
						e=true;
					}
					//Changing lambda on the diagonal and looping again
					else {
						lambda*=10; //Other suggested step of lambda*=3.1622777;
						//lambda*=3.1622777;
						outputFinalPatt=false;
						//if (!err) err=OutputPattern(path,name, iter, converged);
					}
					if (nConvParam==nParam){
						converged=true;
						cout<<"\nLevenberg-Marquardt Converged!\n";
						outFile<<"\nLevenberg-Marquardt Converged!\n";
					}
				}
			}while(!e&&!err&&!converged);
			for (int m=0; m<nParam; m++){
				if (tempA[m]) delete [] tempA[m];
			}
			if (tempA) delete [] tempA;
		}
		if (!err){
			CalcCovarMatrix(bestChiSq);
			cout<<"----------------------------"<<endl;
			cout<<"Best Rwp: "<<oldRwp<<"\n";
			cout<<"Rexp: "<<bestRexp<<"\n";
			cout<<"Best ChiSq: "<<bestChiSq<<"\n";
			outFile<<"----------------------------"<<endl;
			outFile<<"Best Rwp: "<<oldRwp<<"\n";
			outFile<<"Rexp: "<<bestRexp<<"\n";
			outFile<<"Best ChiSq: "<<bestChiSq<<"\n";

			cout<<"\nIndex Name  Val  Std.Dev.\n";
			for(int j=0;j<nParam;j++){
				double stdev=sqrt(gsl_matrix_get(covar,j, j));
				cout<<j+1<<"   "<<indParams[j]->name<<"   "<< saveParams[j]<<"   "<<stdev;
				if (stdev>abs(saveParams[j])||gsl_matrix_get(covar,j,j)<0){
					cout<<" "<<"possible correlation";
				}
				cout<<endl;
				outFile<<indParams[j]->name<<"   "<< saveParams[j]<<endl;
			}
			cout<<"\nCovariance by Index\n";
			cout<<"XXX ";
			for (int j=0;j<nParam-1;j++){
				cout<<j+1<<" ";
			}
			cout<<endl;
			for (int j=1;j<nParam;j++){
				cout<<j+1<<" ";
				for (int k=0;k<j;k++){
					cout<<gsl_matrix_get(covar,j, k)<<" ";
				}
				cout<<endl;
			}
			cout<<"----------------------------"<<endl;
			outFile<<"----------------------------"<<endl;
			if (outputFinalPatt){
				CalcResidual();
				if (!err) OutputPattern(path,name, iter, converged);
			}
		}
	}
	outFile.close();
	delete [] saveParams;
	return err;
}


void Compare::CalcRwp()
{
	float sumSq=0;
	rwp=0;
	for(int i=0; i<Pattern::nInt; i++){
		rwp+=(obsI->I[i]-finalCalcI->I[i])*(obsI->I[i]-finalCalcI->I[i])/obsI->stdDevSq[i];

	}
	if (rwp>0){
		for (int i=0;i<Pattern::nInt; i++){
			sumSq+=obsI->I[i]*obsI->I[i]/obsI->stdDevSq[i];
		}
		rwp=sqrt(rwp/sumSq);

		rexp=sqrt((Pattern::nInt-nParam)/sumSq);

		chisq=rwp/rexp;

	}
	else {
		printf("\nError: Undefined Rwp.. check parameter and observed intensity values\n");
		exit(0);
	}
}

// Minimizer keeps a personal array of the free parameters to be minimized,
// This function creates that array
int Compare::InitParams()
{
	//Local Variables
	int i=0;
	Shape *curShape=firstShape;
	//Dynamic Alloc
	try{
		indParams = new Param* [nParam];
	}
	catch (exception& e){
		printf("\nException found in InitParams alloc: %s\n", e.what());
		return -1;
	}
	// Goes through possible variables and if they are free then adds them to a central list called indParams
	if (!Pattern::monoAngle2T->fixed) {;} //No derivative function written yet for this function.
	if (!Pattern::scale->fixed) {indParams[i]=Pattern::scale; i++;}
	if (!Pattern::B->fixed){indParams[i]=Pattern::B; i++;}
	if (Pattern::t){
		if (!Pattern::t->fixed){indParams[i]=Pattern::t; i++;}
	}
	if (Pattern::ChebBG!=0){
		for (int j=0; j<Pattern::nChebBG;j++){
			if (!Pattern::ChebBG[j]->fixed){indParams[i]=Pattern::ChebBG[j];i++;}
		}
	}
	if (Pattern::ChebScale!=0){
		for (int j=0; j<Pattern::nChebScale;j++){
			if (!Pattern::ChebScale[j]->fixed){indParams[i]=Pattern::ChebScale[j];i++;}
		}
	}
	if (Pattern::nGaussBG!=0){
		for (int j=0; j<Pattern::nGaussBG; j++){
			if (!Pattern::gaussA[j]->fixed){indParams[i]=Pattern::gaussA[j]; i++;}
			if (!Pattern::gaussMu[j]->fixed){indParams[i]=Pattern::gaussMu[j]; i++;}
			if (!Pattern::gaussSig[j]->fixed){indParams[i]=Pattern::gaussSig[j]; i++;}
		}
	}
	if (Pattern::AExp!=0&&Pattern::kExp!=0){
		if (!Pattern::AExp->fixed){indParams[i]=Pattern::AExp; i++;}
		if (!Pattern::kExp->fixed){indParams[i]=Pattern::kExp; i++;}
	}
	if(Pattern::delExp){
		if(!Pattern::delExp->fixed){indParams[i]=Pattern::delExp; i++;}
	}
	if (Pattern::sampleDisplace!=0){
		if (!Pattern::sampleDisplace->fixed){indParams[i]=Pattern::sampleDisplace;i++;}
	}
	for (int j=0;j<curShape->nParams_a;j++){
		if (!curShape->a[j]->fixed){indParams[i]=curShape->a[j]; i++;}
	}
	if (curShape->distribution=="LogNorm"){
		if (!curShape->mu->fixed){indParams[i]=curShape->mu; i++;}
		if (!curShape->sigma->fixed) {indParams[i]=curShape->sigma; i++;}
	}
	if (curShape->distribution=="Gamma"){
		if (!curShape->alpha->fixed){indParams[i]=curShape->alpha; i++;}
		if (!curShape->beta->fixed) {indParams[i]=curShape->beta; i++;}
	}
	if (curShape->distribution=="Delta"){
		if (!curShape->deltaSize->fixed) {indParams[i]=curShape->deltaSize; i++;}
	}
	if (!curShape->f->fixed) {indParams[i]=curShape->f; i++;}
	if (!curShape->kappa->fixed) {indParams[i]=curShape->kappa; i++;}
	if (!curShape->wShape->fixed) {indParams[i]=curShape->wShape; i++;}
	curShape=curShape->nextShape;
	while (curShape!=0){

		for (int j=0;j<curShape->nParams_a;j++){
			if (!curShape->a[j]->fixed&&!curShape->a[j]->global){indParams[i]=curShape->a[j]; i++;}
		}
		if (curShape->distribution=="LogNorm"){
			if (!curShape->mu->global&&!curShape->mu->fixed) {indParams[i]=curShape->mu; i++;}
			if (!curShape->sigma->global&&!curShape->sigma->fixed) {indParams[i]=curShape->sigma; i++;}
		}
		if (curShape->distribution=="Gamma"){
			if (!curShape->alpha->global&&!curShape->alpha->fixed) {indParams[i]=curShape->alpha; i++;}
			if (!curShape->beta->global&&!curShape->beta->fixed) {indParams[i]=curShape->beta; i++;}
		}
		if (curShape->distribution=="Delta"){
			if (!curShape->deltaSize->global&&!curShape->deltaSize->fixed) {indParams[i]=curShape->deltaSize; i++;}
		}
		if (!curShape->f->global&&!curShape->f->fixed) {indParams[i]=curShape->f; i++;}
		if (!curShape->kappa->global&&!curShape->kappa->fixed) {indParams[i]=curShape->kappa; i++;}
		if (!curShape->wShape->fixed) {indParams[i]=curShape->wShape; i++;}
		curShape=curShape->nextShape;
	}
	return 0;

}

//Calculates pattern assuming only one atomic species... If needed code to account for more atomic species must be added
// In order to do this care in calculating distances between different atomic species and different atomic scattering factors must be implimented.
int Compare::CalcPattern(bool calcSizeDerivMult, bool volNorm)
{
	//Locals
	int err=0;
	Shape *curShape=firstShape;
	//Dynamic Alloc
	try{
		if (!initCalcI){
			initCalcI = new Pattern* [Pattern::nLambda];
			for (int i=0;i<Pattern::nLambda;i++){
				initCalcI[i]=new Pattern("PartofCalcPatt");
				initCalcI[i]->AllocPattern();
			}
		}
		if (!finalCalcI){
			finalCalcI=new Pattern("TotCalcPattern");
			finalCalcI->AllocPattern();
		}
	}
	catch (exception& e){
		printf("\nException found in CalcPattern alloc: %s\n", e.what());
		exit(1);
	}
//firstCalc is pattern before the background is added and finalCalc is with background
	//Zero patterns
	for (int j=0; j<Pattern::nLambda;j++){
		initCalcI[j]->ZeroPattern();
	}
	finalCalcI->ZeroPattern();

	// Calc offset of 2theta is necessary
	if (Pattern::sampleDisplace!=0){
		initCalcI[0]->ApplySampleDisplacement();
	}
	if ((Pattern::ax!=0||Pattern::bx!=0||Pattern::cx!=0||Pattern::dx!=0||Pattern::ex!=0)&&(!Pattern::synchThetaShift)){
		initCalcI[0]->ApplySynchThetaShift();
	}
	//Calc pattern for each shape and weight accordingly.
	while (curShape!=0 && !err){
		if (curShape->wShape->GetVal()!=0){
			if (!err) err=curShape->CalcPattern();
			if (!err){
				for (int j=0; j<Pattern::nLambda; j++){
					for (int i=0;i<Pattern::nInt;i++){
						initCalcI[j]->I[i]+=curShape->wShape->GetVal()*curShape->diffPatt[j]->I[i];
					}
				}
			}
		}
		curShape=curShape->nextShape;
	}
	if (!err)
	{
		nAtoms=CalcNumAtoms();
		for (int j=0;j<Pattern::nLambda;j++){
			initCalcI[j]->ApplyWaveLenDependMultFactors(j);
			initCalcI[j]->ScalePattern(1.0/nAtoms);
			initCalcI[j]->ApplyGeomDependMultFactors();
			for (int i=0;i<Pattern::nInt&&!err;i++){
				finalCalcI->I[i]+=Pattern::wLambda[j]*initCalcI[j]->I[i];
			}
		}
		//finalCalcI->ApplyGeomDependMultFactors();
		finalCalcI->ApplyChebyshevScale();

	//Apply Background
		finalCalcI->initBG=false;
		if (Pattern::nChebBG!=0){
			if (Pattern::chebPolys==0) finalCalcI->CalcChebyshevPolys();
			finalCalcI->ApplyChebyshevBG();
		}
		if (Pattern::nGaussBG!=0){
			finalCalcI->CalcGaussBG();
			finalCalcI->ApplyGaussBG();
		}
		if (Pattern::AExp!=0){
			finalCalcI->ApplyExpBG();
		}
	}
	return err;
}
int Compare::OutputPattern(string _intOutputPath, string _intOutputName, int iter, bool conv)
{
	int err=0;
	//Write to File
	Shape *curShape=firstShape;
	string filePath=_intOutputPath+"/";
	string fileName="intensity."+_intOutputName+".fit";
	string filename=filePath+fileName;
	fstream file(filename.c_str(), ios_base::out|ios_base::app);
	if(file.is_open()){
		if (conv){
			file<<"Converged Calculated Pattern\n";
		}
		file<<"Calculated Pattern for Iteration: "<<iter<<endl;
	}
	else {
		printf("Error: Could not output to file: %s\n", filename.c_str());
	}

	while(curShape!=0 &&!err){
		if (curShape->wShape->GetVal()!=0){
			err=curShape->WriteIntFileHeader(filePath, fileName);
		}
		curShape=curShape->nextShape;
	}
	if (!err)
		err=finalCalcI->OutputToFile(filePath, fileName);
	return err;
}
/*
//Still need to Determine the format of the obs (xy or debye)
int Compare::InputObservedInt(string _obsPattPath, string _obsPattName)
{
	//Locals
	int err=0;
	char path[_MAX_DIR];
	struct _finddata_t file;
	intptr_t hFile;
	string search=_obsPattName;
	//Mem alloc
	try {
		obsI=new Pattern(_obsPattName);
	}
	catch (exception& e){
		printf("Exception found in InputObservedInt alloc: %s", e.what());
		return -1;
	}

	if (!_chdir(_obsPattPath.c_str())){
		//Change dir to the input path
		_getdcwd (0, path, _MAX_DIR);
		printf("\nCurrent directory is: %s \n", path);
		//Search in dir for filename
		if((hFile=_findfirst(search.c_str(),&file))!=-1L){
			printf("Observed File Found: \n%s\n", file.name);
			err=obsI->ReadIntensityFile();
			if (!err)
				printf("File %s Input\n", file.name);
		}
		else{
			printf("/nError: Cannot find desired intensity file.. Check name and path\n");
			err=-1;
		}
	}
	else{
		printf("\nError: Cannot change current working directory... Desired directory may not exist\n");
		err=-1;
	}

	return err;

}
*/
//Currently is not accurate
int Compare::CalcDerivW(Shape *curShape)
{
	int err=0,i,j;
	double scalefactor = curShape->GetShapeShellD();
	if (derivfSize_mu==0) derivfSize_mu= new double [curShape->maxShell];//(double*) malloc(curShape->maxShell*sizeof(double));
	if (derivfSize_sig==0) derivfSize_sig= new double [curShape->maxShell];//(double*) malloc(curShape->maxShell*sizeof(double));
	if (derivfSize_mu!=0&&derivfSize_sig!=0){
		//For other shapes other than spheres Diameter is replaced by the maximum distance in the shape before any distrotions are applied.
		for (i=0;i<curShape->maxShell;i++){
			derivfSize_mu[i]=curShape->p[i]*(log((i+1)*scalefactor)-curShape->mu->GetVal())/(curShape->sigma->GetVal()*curShape->sigma->GetVal());
			derivfSize_sig[i]=curShape->p[i]*((log((i+1)*scalefactor)-curShape->mu->GetVal())*(log((i+1)*scalefactor)-curShape->mu->GetVal())/(curShape->sigma->GetVal()*curShape->sigma->GetVal()*curShape->sigma->GetVal())-1/curShape->sigma->GetVal());
		}
		//Combine the derivative of the size distribution by use of the quotient rule (The normalization factor of the number of atoms is also dependant on the size distribution)
		double dNdMu=GetDerivNumAtoms_mu(curShape), dNdSig=GetDerivNumAtoms_sig(curShape), totN=curShape->nAtoms;
		for (i=0;i<curShape->maxShell;i++){
			derivfSize_mu[i]=(totN*derivfSize_mu[i]-curShape->p[i]*dNdMu)/(totN*totN);
			derivfSize_sig[i]=(totN*derivfSize_sig[i]-curShape->p[i]*dNdSig)/(totN*totN);
		}
		//int endSize=int(ceil(exp(curShape->mu->GetVal()+3*curShape->sigma->GetVal())/scalefactor));

		if (curShape->dW_mu==0) curShape->dW_mu=new double [curShape->maxShell];//(double*) malloc(curShape->maxShell*sizeof(double));
		if (curShape->dW_sig==0) curShape->dW_sig=new double [curShape->maxShell];//(double*) malloc(curShape->maxShell*sizeof(double));
		if (curShape->dW_mu!=0 && curShape->dW_sig!=0){
			for (i=0;i<curShape->maxShell;i++){
				curShape->dW_mu[i]=0;
				curShape->dW_sig[i]=0;
				//Weights for shells without SR
				for (j=i+curShape->nRelaxShells;j<curShape->endSize;j++){
					curShape->dW_mu[i]+=derivfSize_mu[j];
					curShape->dW_sig[i]+=derivfSize_sig[j];
				}
			}
		}
		else {
			printf("Error in size distribution derivative memory allocation\n");
			err=-1;
		}
	}
	else {
		printf("Size Function Value vector memory allocation error\n");
		err=-1;
	}
	return err;
}

double Compare::GetDerivNumAtoms_mu(Shape *curShape)
{
	double nTot=0;
	Position *curPos=curShape->firstPos;
	int i,j;
	for (j=0;j<curShape->maxShell;j++){
		for (i=j;i<curShape->maxShell;i++){
			nTot+=curPos->nAtoms*derivfSize_mu[i];
		}
		curPos=curPos->nextPos;
	}
	return nTot;
}
double Compare::GetDerivNumAtoms_sig(Shape *curShape)
{
	double nTot=0;
	Position *curPos=curShape->firstPos;
	int i,j;
	for (j=0;j<curShape->maxShell;j++){
		for (i=j;i<curShape->maxShell;i++){
			nTot+=curPos->nAtoms*derivfSize_sig[i];
		}
		curPos=curPos->nextPos;
	}
	return nTot;
}

//Function to call the necesary routines to calculate the gradient of the intensity
void Compare::CalcGrad()
{
	//Locals
	int err=0;
	//Dynamic alloc
	try{
		if (gradI==0){
			gradI= new double* [nParam];
			for (int i=0;i<nParam;i++){
				gradI[i]=new double [Pattern::nInt];
			}
		}
	}
	catch (exception& e){
		printf("\nException found in CalcGrad alloc: %s\n", e.what());
		exit(1);
	}

	//Loop through independant parameters
	//need to refine so that it does not have to match up the the way they are input in the list
	if (indParams!=0){
		for (int i=0;i<nParam&&!err;i++){

			if (indParams[i]->name.find("Cheb_BG")!=string::npos){
				printf("\nCalculating deriv wrt Chebyshev background\n");
				char num[1];
				indParams[i]->name.copy(num,1,7);
				int k=int(num[0])-48;//subtracting 48 b/c '0'=48
				for (int j=0; j<Pattern::nInt;j++){
					gradI[i][j]=Pattern::chebPolys[k][j];
				}
			}
			else if (indParams[i]->name.find("Cheb_Scale")!=string::npos){
				char num[1];
				indParams[i]->name.copy(num,1,10);
				CalcDeriv_ChebScale(gradI[i],int(num[0])-48);
			}
			else if (indParams[i]->name.find("Gauss Amp")!=string::npos){
				printf("\nCalculating deriv wrt %s\n", indParams[i]->name.c_str());
				char num[1];
				indParams[i]->name.copy(num,1,20);
				CalcDeriv_gaussA(gradI[i],int(num[0])-48); //subtracting 48 b/c '0'=48
			}
			else if (indParams[i]->name.find("Gauss Mu")!=string::npos){
				printf("\nCalculating deriv wrt %s\n", indParams[i]->name.c_str());
				char num[1];
				indParams[i]->name.copy(num,1,19);
				CalcDeriv_gaussMu(gradI[i],int(num[0])-48); //subtracting 48 b/c '0'=48
			}
			else if (indParams[i]->name.find("Gauss Sigma")!=string::npos){
				printf("\nCalculating deriv wrt %s\n", indParams[i]->name.c_str());
				char num[1];
				indParams[i]->name.copy(num,1,22);
				CalcDeriv_gaussSig(gradI[i],int(num[0])-48); //subtracting 48 b/c '0'=48
			}
			else if (indParams[i]->name.find("ExpDecayBG Amp")!=string::npos){
				CalcDeriv_expAmp(gradI[i]);
			}
			else if (indParams[i]->name.find("ExpDecayBG Kappa")!=string::npos){
				CalcDeriv_expKappa(gradI[i]);
			}
			else if (indParams[i]->name.find("ExpDecayBG Del")!=string::npos){
				CalcDeriv_expDel(gradI[i]);
			}
			else if (indParams[i]->name.find("Scale")!=string::npos){
				CalcDeriv_scale(gradI[i]);
			}
			else if (indParams[i]->name.find("Sample Displacement")!=string::npos)
				err=CalcDeriv_smplDisplace(gradI[i]);
			else if (indParams[i]->name.find("AbsorptionThickness")!=string::npos)
				err=CalcDeriv_t(gradI[i]);
			else if (indParams[i]->name.find("a_Poly")!=string::npos){
				char num[1];
				indParams[i]->name.copy(num,1,6);
				err=CalcDeriv_a_D(gradI[i], int(num[0])-48, indParams[i]->shapeName);
			}
			else if (indParams[i]->name.find("mu")!=string::npos)
				err=CalcDeriv_mu(gradI[i], indParams[i]->shapeName);
			else if (indParams[i]->name.find("sigma")!=string::npos)
				err=CalcDeriv_sigma(gradI[i], indParams[i]->shapeName);
			else if (indParams[i]->name.find("alpha")!=string::npos)
				err=CalcDeriv_alpha(gradI[i], indParams[i]->shapeName);
			else if (indParams[i]->name.find("beta")!=string::npos)
				err=CalcDeriv_beta(gradI[i], indParams[i]->shapeName);
			else if (indParams[i]->name.find("deltaSize")!=string::npos)
				err=CalcDeriv_deltaSize(gradI[i], indParams[i]->shapeName);
			else if (indParams[i]->name.find("kappa")!=string::npos)
				err=CalcDeriv_kappa(gradI[i], indParams[i]->shapeName);
			else if (indParams[i]->name.find("w")!=string::npos)
				err=CalcDeriv_wShape(gradI[i], indParams[i]->shapeName);
			else if (indParams[i]->name.find("f")!=string::npos)
				err=CalcDeriv_f(gradI[i], indParams[i]->shapeName);
			else if (indParams[i]->name.find("a")!=string::npos)
				err=CalcDeriv_a(gradI[i], indParams[i]->shapeName);
			else if (indParams[i]->name.find("B")!=string::npos)
				err=CalcDeriv_B(gradI[i]);
			else{
				printf("Error in finding correct function to calculate derivative of parameter: %s", indParams[i]->name.c_str());
				exit(0);
			}
		}
	}
}
int Compare::CalcDeriv_B(double *outputGrad)
{
	Pattern derivPatt;
	derivPatt.AllocPattern();
	printf("\nCalculating gradient wrt B\n");
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
		for (int j=0;j<Pattern::nLambda;j++){
			derivPatt.I[i]+=Pattern::wLambda[j]*initCalcI[j]->I[i]*(-Pattern::offsetq[j][i]*Pattern::offsetq[j][i]/(4*PI*PI*2));
		}
	}
	//derivPatt.ApplyGeomDependMultFactors();
	derivPatt.ApplyChebyshevScale();
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=derivPatt.I[i];
	}
	return 0;
}
int Compare::CalcDeriv_t(double *outputGrad)
{
	//Form of derivative is to account for scaled absorption expression.
	//(see Pattern::ApplyAbsorption)
	Pattern derivPatt;
	derivPatt.AllocPattern();
	printf("\nCalculating gradient wrt t\n");
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
		for (int j=0;j<Pattern::nLambda;j++){
			derivPatt.I[i]+=Pattern::wLambda[j]*initCalcI[j]->I[i];
		}
	}
	derivPatt.ApplyPolarizationFactor();
	if(Pattern::multCos){
		derivPatt.MultCosTheta();
	}
	if (Pattern::divSin){
		derivPatt.DivSinTheta();
	}
	derivPatt.ApplyChebyshevScale();
	if (Pattern::geometry=="transmission"||Pattern::geometry=="transmission_bathed"){
		derivPatt.ApplyAbsorption();
	}

	for (int i=0;i<Pattern::nInt;i++){
		if (Pattern::geometry=="transmission"||Pattern::geometry=="transmission_bathed"){
			derivPatt.I[i]*=(1/Pattern::t->GetVal())-Pattern::absCoef->GetVal()/cos(Pattern::offset2Theta[i]*PI/360);
		}
		else if (Pattern::geometry=="reflection"){
			double abs=exp(-2*Pattern::absCoef->GetVal()*Pattern::t->GetVal()/sin(Pattern::offset2Theta[i]*PI/360.0));
			//derivPatt.I[i]*=2*Pattern::absCoef->GetVal()*abs/sin(Pattern::offset2Theta[i]*PI/360.0);
			derivPatt.I[i]*=-abs/sin(Pattern::offset2Theta[i]*PI/360.0);
		}
		else if (Pattern::geometry=="reflection_thinFilm"){
			double abs=exp(-Pattern::absCoef->GetVal()*Pattern::t->GetVal()/sin(Pattern::offset2Theta[i]*PI/360.0));
			//derivPatt.I[i]*=Pattern::absCoef->GetVal()*abs/sin(Pattern::offset2Theta[i]*PI/360.0);
			derivPatt.I[i]*=-abs/sin(Pattern::offset2Theta[i]*PI/360.0);
		}
		else {
			cout<<"Geometry :"<<Pattern::geometry<<" not supported in CalcDeriv_t.\n";
			exit(0);
		}
	}

	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=derivPatt.I[i];
	}
	return 0;
}
/*
int Compare::CalcDeriv_mu(double *outputGrad, string shapeName)
{
	int err=0;
	Shape *curShape=firstShape;
	double nA=0;
	double delta=0.001;
	Pattern *newPattern=new Pattern[Pattern::nLambda];
	try{
		for (int j=0;j<Pattern::nLambda; j++)
			newPattern[j].I=new double[Pattern::nInt];
	}
	catch(exception &e){
		cout<<"Exception in Compare::CalcDeriv_mu: "<<e.what();
		exit(1);
	}
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
		for (int j=0;j<Pattern::nLambda;j++){
			newPattern[j].I[i]=0;
		}
	}
	printf("\nCalculating gradient wrt mu\n");
	while (curShape!=0 && !err){
		if (curShape->name==shapeName||shapeName=="global"){
			if (curShape->wShape->GetVal()!=0){
				Shape tempShape;
				err=CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/dMu";
				if (!err){
					tempShape.mu->GetVal()*=1.001;
					delta=tempShape.mu->GetVal()-curShape->mu->GetVal();
					err=tempShape.CalcPattern();
				}
				nA+=tempShape.GetNumAtoms()*tempShape.wShape->GetVal();
				for(int i=0;i<Pattern::nInt;i++){
					for (int j=0;j<Pattern::nLambda;j++){
						newPattern[j].I[i]+=tempShape.wShape->GetVal()*tempShape.diffPatt[j]->I[i];
					}
				}
					/*
					for (int j=0;j<Pattern::nLambda&&!err;j++){
						for (i=0;i<Pattern::nInt&&!err;i++){
							tempShape.diffPatt[j]->I[i]-=curShape->diffPatt[j]->I[i];
						}
						tempShape.diffPatt[j]->ApplyMultFactors(j);
						err=tempShape.diffPatt[j]->ApplyChebyshevScale();
						for (i=0;i<Pattern::nInt&&!err;i++){
							outputGrad[i]+=Pattern::wLambda[j]*curShape->wShape->GetVal()*tempShape.diffPatt[j]->I[i]/(curShape->mu->GetVal()*.001);
						}
					}

					//}
				}*//*
			}
		}
		else{
			nA+=curShape->GetNumAtoms()*curShape->wShape->GetVal();
			for(int i=0;i<Pattern::nInt;i++){
				for (int j=0;j<Pattern::nLambda;j++){
					newPattern[j].I[i]+=curShape->wShape->GetVal()*curShape->diffPatt[j]->I[i];
				}
			}
		}
		curShape=curShape->nextShape;
	}
	for (int j=0;j<Pattern::nLambda&&!err;j++){
		newPattern[j].ApplyMultFactors(j);
		newPattern[j].ScalePattern(1.0/nA);
		newPattern[j].ApplyChebyshevScale();
		for (int i=0;i<Pattern::nInt&&!err;i++){
			outputGrad[i]+=Pattern::wLambda[j]*(newPattern[j].I[i]-initCalcI[j]->I[i])/(delta);
		}
	}
	delete [] newPattern;
	return err;
}
*/

int Compare::CalcDeriv_mu(double *outputGrad, string shapeName)
{
	int err=0;
	Shape *curShape=firstShape;
	double delNAtoms=0;
	double delta=0.001;
	Pattern *newPattern=new Pattern[Pattern::nLambda];
	for (int j=0;j<Pattern::nLambda; j++){
		newPattern[j].AllocPattern();
	}

	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
	}
	printf("\nCalculating gradient wrt mu\n");
	while (curShape!=0 && !err){
		if (curShape->name==shapeName||shapeName=="global"){
			if (curShape->wShape->GetVal()!=0){
				Shape tempShape;
				CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/dMu";
				if (!err){
					if (tempShape.mu->GetVal()!=0){
						delta=tempShape.mu->GetVal()*0.001;
					}
					tempShape.mu->ChangeVal(tempShape.mu->GetVal()+delta);
					err=tempShape.CalcPattern();
				}
				tempShape.CalcNumAtoms();
				delNAtoms+=tempShape.nAtoms*tempShape.wShape->GetVal();

				for(int i=0;i<Pattern::nInt;i++){
					for (int j=0;j<Pattern::nLambda;j++){
						newPattern[j].I[i]+=tempShape.wShape->GetVal()*tempShape.diffPatt[j]->I[i];
					}
				}

			}
		}
		else{
			for(int i=0;i<Pattern::nInt;i++){
				for (int j=0;j<Pattern::nLambda;j++){
					newPattern[j].I[i]+=curShape->wShape->GetVal()*curShape->diffPatt[j]->I[i];
				}
			}
			delNAtoms+=curShape->wShape->GetVal()*curShape->nAtoms;
		}


		curShape=curShape->nextShape;
	}
	for (int j=0;j<Pattern::nLambda&&!err;j++){
		newPattern[j].ScalePattern(1/delNAtoms);
		newPattern[j].ApplyWaveLenDependMultFactors(j);
		newPattern[j].ApplyGeomDependMultFactors();
		for (int i=0;i<Pattern::nInt;i++){
			newPattern[j].I[i]-=initCalcI[j]->I[i];
		}
		//newPattern[j].ApplyGeomDependMultFactors();
		newPattern[j].ApplyChebyshevScale();
		for (int i=0;i<Pattern::nInt&&!err;i++){
			outputGrad[i]+=Pattern::wLambda[j]*newPattern[j].I[i]/delta;
		}
	}
	delete [] newPattern;
	return err;
}
/*
int Compare::CalcDeriv_sigma(double *outputGrad, string shapeName)
{
	int err=0, i;
	Shape *curShape=firstShape;

	for (i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
	}
	printf("\nCalculating gradient wrt sigma\n");
	while (curShape!=0 && !err){
		if (curShape->wShape->GetVal()!=0){
			if (curShape->name==shapeName||shapeName=="global"){
				Shape tempShape;
				err=CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/dSigma";
				if (!err){
					tempShape.sigma->val*=1.001;
					err=tempShape.CalcPattern();
					for (int j=0;j<Pattern::nLambda&&!err;j++){
						for (i=0;i<Pattern::nInt&&!err;i++){
							tempShape.diffPatt[j]->I[i]-=curShape->diffPatt[j]->I[i];
						}
						tempShape.diffPatt[j]->ApplyMultFactors(j);
						tempShape.diffPatt[j]->ScalePattern(1.0/nAtoms);
						err=tempShape.diffPatt[j]->ApplyChebyshevScale();
						for (i=0;i<Pattern::nInt&&!err;i++){
							outputGrad[i]+=Pattern::wLambda[j]*curShape->wShape->val*tempShape.diffPatt[j]->I[i]/(curShape->sigma->val*.001);
						}
					}
					//}
				}
			}
		}
		curShape=curShape->nextShape;
	}
	return err;
}
*/
/*
int Compare::CalcDeriv_sigma(double *outputGrad, string shapeName)
{
	int err=0;
	Shape *curShape=firstShape;
	double nA=0;
	double delta=0.001;
	Pattern *newPattern=new Pattern[Pattern::nLambda];
	try{
		for (int j=0;j<Pattern::nLambda; j++)
			newPattern[j].I=new double[Pattern::nInt];
	}
	catch(exception &e){
		cout<<"Exception in Compare::CalcDeriv_sigma: "<<e.what();
		exit(1);
	}
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
		for (int j=0;j<Pattern::nLambda;j++){
			newPattern[j].I[i]=0;
		}
	}
	printf("\nCalculating gradient wrt sigma\n");
	while (curShape!=0 && !err){
		if (curShape->name==shapeName||shapeName=="global"){
			if (curShape->wShape->val!=0){
				Shape tempShape;
				err=CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/dSigma";
				if (!err){
					tempShape.sigma->val*=1.001;
					delta=tempShape.sigma->val-curShape->sigma->val;
					err=tempShape.CalcPattern();
				}
				nA+=tempShape.GetNumAtoms()*tempShape.wShape->val;
				for(int i=0;i<Pattern::nInt;i++){
					for (int j=0;j<Pattern::nLambda;j++){
						newPattern[j].I[i]+=tempShape.wShape->val*tempShape.diffPatt[j]->I[i];
					}
				}
			}
		}
		else{
			nA+=curShape->GetNumAtoms()*curShape->wShape->val;
			for(int i=0;i<Pattern::nInt;i++){
				for (int j=0;j<Pattern::nLambda;j++){
					newPattern[j].I[i]+=curShape->wShape->val*curShape->diffPatt[j]->I[i];
				}
			}
		}
		curShape=curShape->nextShape;
	}
	for (int j=0;j<Pattern::nLambda&&!err;j++){
		newPattern[j].ApplyMultFactors(j);
		newPattern[j].ScalePattern(1.0/nA);
		newPattern[j].ApplyChebyshevScale();
		for (int i=0;i<Pattern::nInt&&!err;i++){
			outputGrad[i]+=Pattern::wLambda[j]*(newPattern[j].I[i]-initCalcI[j]->I[i])/(delta);
		}
	}
	delete [] newPattern;
	return err;
}
*/
/*
int Compare::CalcDeriv_sigma(double *outputGrad, string shapeName)
{
	int err=0;
	Shape *curShape=firstShape;
	double delNAtoms=0;
	double delta=0.001;
	Pattern *newPattern=new Pattern[Pattern::nLambda];
	Pattern *oldPatt=new Pattern[Pattern::nLambda];
	try{
		for (int j=0;j<Pattern::nLambda; j++){
			newPattern[j].I=new double[Pattern::nInt];
			oldPatt[j].I=new double [Pattern::nInt];
		}

	}
	catch(exception &e){
		cout<<"Exception in Compare::CalcDeriv_sigma: "<<e.what();
		exit(1);
	}
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
		for (int j=0;j<Pattern::nLambda;j++){
			newPattern[j].I[i]=0;
			oldPatt[j].I[i]=0;
		}
	}
	printf("\nCalculating gradient wrt sigma\n");
	while (curShape!=0 && !err){
		if (curShape->name==shapeName||shapeName=="global"){
			if (curShape->wShape->val!=0){
				Shape tempShape;
				err=CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/dSigma";
				if (!err){
					tempShape.sigma->val*=1.01;
					delta=tempShape.sigma->val-curShape->sigma->val;
					err=tempShape.CalcPattern();
				}
				delNAtoms+=(tempShape.GetNumAtoms()-curShape->GetNumAtoms())*tempShape.wShape->val/delta;
				for(int i=0;i<Pattern::nInt;i++){
					for (int j=0;j<Pattern::nLambda;j++){
						newPattern[j].I[i]+=tempShape.wShape->val*(tempShape.diffPatt[j]->I[i]-curShape->diffPatt[j]->I[i])/delta;
					}
				}
			}
		}
		for(int i=0;i<Pattern::nInt;i++){
			for (int j=0;j<Pattern::nLambda;j++){
				oldPatt[j].I[i]+=curShape->wShape->val*curShape->diffPatt[j]->I[i];
			}
		}
		curShape=curShape->nextShape;
	}
	for (int j=0;j<Pattern::nLambda&&!err;j++){
		newPattern[j].ScalePattern(1/nAtoms);
		oldPatt[j].ScalePattern(delNAtoms/(nAtoms*nAtoms));
		for (int i=0;i<Pattern::nInt;i++){
			newPattern[j].I[i]-=oldPatt[j].I[i];
		}
		newPattern[j].ApplyMultFactors(j);
		newPattern[j].ApplyChebyshevScale();
		for (int i=0;i<Pattern::nInt&&!err;i++){
			outputGrad[i]+=Pattern::wLambda[j]*newPattern[j].I[i];
		}
	}
	delete [] newPattern;
	delete [] oldPatt;
	return err;
}
*/
int Compare::CalcDeriv_sigma(double *outputGrad, string shapeName)
{
	int err=0;
	Shape *curShape=firstShape;
	double delNAtoms=0;
	double delta=0.001;
	Pattern *newPattern=new Pattern[Pattern::nLambda];
	for (int j=0;j<Pattern::nLambda; j++){
		newPattern[j].AllocPattern();
	}

	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
	}
	printf("\nCalculating gradient wrt sigma\n");
	while (curShape!=0 && !err){
		if (curShape->name==shapeName||shapeName=="global"){
			if (curShape->wShape->GetVal()!=0){
				Shape tempShape;
				CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/dsigma";
				if (!err){
					if (tempShape.sigma->GetVal()!=0){
						delta=tempShape.sigma->GetVal()*.001;
					}
					tempShape.sigma->ChangeVal(tempShape.sigma->GetVal()+delta);
					err=tempShape.CalcPattern();
				}
				tempShape.CalcNumAtoms();
				delNAtoms+=tempShape.nAtoms*tempShape.wShape->GetVal();

				for(int i=0;i<Pattern::nInt;i++){
					for (int j=0;j<Pattern::nLambda;j++){
						newPattern[j].I[i]+=tempShape.wShape->GetVal()*tempShape.diffPatt[j]->I[i];
					}
				}

			}
		}
		else{
			for(int i=0;i<Pattern::nInt;i++){
				for (int j=0;j<Pattern::nLambda;j++){
					newPattern[j].I[i]+=curShape->wShape->GetVal()*curShape->diffPatt[j]->I[i];
				}
			}
			delNAtoms+=curShape->wShape->GetVal()*curShape->nAtoms;
		}


		curShape=curShape->nextShape;
	}
	for (int j=0;j<Pattern::nLambda&&!err;j++){
		newPattern[j].ScalePattern(1/delNAtoms);
		newPattern[j].ApplyWaveLenDependMultFactors(j);
		newPattern[j].ApplyGeomDependMultFactors();
		for (int i=0;i<Pattern::nInt;i++){
			newPattern[j].I[i]-=initCalcI[j]->I[i];
		}
		//newPattern[j].ApplyGeomDependMultFactors();
		newPattern[j].ApplyChebyshevScale();
		for (int i=0;i<Pattern::nInt&&!err;i++){
			outputGrad[i]+=Pattern::wLambda[j]*newPattern[j].I[i]/delta;
		}
	}
	delete [] newPattern;
	return err;
}
/*
int Compare::CalcDeriv_alpha(double *outputGrad, string shapeName)
{
	int err=0, i;
	Shape *curShape=firstShape;

	for (i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
	}
	printf("\nCalculating gradient wrt alpha\n");
	while (curShape!=0 && !err){
		if (curShape->name==shapeName||shapeName=="global"){
			if (curShape->wShape->val!=0){
				Shape tempShape;
				err=CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/dAlpha";
				if (!err){
					tempShape.alpha->val*=1.001;
					err=tempShape.CalcPattern();
					for (int j=0;j<Pattern::nLambda;j++){
						for (i=0;i<Pattern::nInt&&!err;i++){
							tempShape.diffPatt[j]->I[i]-=curShape->diffPatt[j]->I[i];
						}
						tempShape.diffPatt[j]->ApplyMultFactors(j);
						tempShape.diffPatt[j]->ScalePattern(1.0/nAtoms);
						err=tempShape.diffPatt[j]->ApplyChebyshevScale();
						for (i=0;i<Pattern::nInt&&!err;i++){
							outputGrad[i]+=Pattern::wLambda[j]*curShape->wShape->val*tempShape.diffPatt[j]->I[i]/(curShape->alpha->val*.001);
						}
					}
					//}
				}
			}
		}
		curShape=curShape->nextShape;
	}

	return err;
}
*/
/*
int Compare::CalcDeriv_alpha(double *outputGrad, string shapeName)
{
	int err=0;
	Shape *curShape=firstShape;
	double nA=0;
	double delta=0.001;
	Pattern *newPattern=new Pattern[Pattern::nLambda];
	try{
		for (int j=0;j<Pattern::nLambda; j++)
			newPattern[j].I=new double[Pattern::nInt];
	}
	catch(exception &e){
		cout<<"Exception in Compare::CalcDeriv_alpha: "<<e.what();
		exit(1);
	}
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
		for (int j=0;j<Pattern::nLambda;j++){
			newPattern[j].I[i]=0;
		}
	}
	printf("\nCalculating gradient wrt Alpha\n");
	while (curShape!=0 && !err){
		if (curShape->name==shapeName||shapeName=="global"){
			if (curShape->wShape->val!=0){
				Shape tempShape;
				err=CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/dAlpha";
				if (!err){
					tempShape.alpha->val*=1.001;
					delta=tempShape.alpha->val-curShape->alpha->val;
					err=tempShape.CalcPattern();
				}
				nA+=tempShape.GetNumAtoms()*tempShape.wShape->val;
				for(int i=0;i<Pattern::nInt;i++){
					for (int j=0;j<Pattern::nLambda;j++){
						newPattern[j].I[i]+=tempShape.wShape->val*tempShape.diffPatt[j]->I[i];
					}
				}
			}
		}
		else{
			nA+=curShape->GetNumAtoms()*curShape->wShape->val;
			for(int i=0;i<Pattern::nInt;i++){
				for (int j=0;j<Pattern::nLambda;j++){
					newPattern[j].I[i]+=curShape->wShape->val*curShape->diffPatt[j]->I[i];
				}
			}
		}
		curShape=curShape->nextShape;
	}
	for (int j=0;j<Pattern::nLambda&&!err;j++){
		newPattern[j].ApplyMultFactors(j);
		newPattern[j].ScalePattern(1.0/nA);
		newPattern[j].ApplyChebyshevScale();
		for (int i=0;i<Pattern::nInt&&!err;i++){
			outputGrad[i]+=Pattern::wLambda[j]*(newPattern[j].I[i]-initCalcI[j]->I[i])/(curShape->alpha->val*.001);
		}
	}
	delete [] newPattern;
	return err;
}
*/
/*
int Compare::CalcDeriv_alpha(double *outputGrad, string shapeName)
{
	int err=0;
	Shape *curShape=firstShape;
	double delNAtoms=0;
	double delta=0.001;
	Pattern *newPattern=new Pattern[Pattern::nLambda];
	Pattern *oldPatt=new Pattern[Pattern::nLambda];
	try{
		for (int j=0;j<Pattern::nLambda; j++){
			newPattern[j].I=new double[Pattern::nInt];
			oldPatt[j].I=new double [Pattern::nInt];
		}

	}
	catch(exception &e){
		cout<<"Exception in Compare::CalcDeriv_alpha: "<<e.what();
		exit(1);
	}
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
		for (int j=0;j<Pattern::nLambda;j++){
			newPattern[j].I[i]=0;
			oldPatt[j].I[i]=0;
		}
	}
	printf("\nCalculating gradient wrt alpha\n");
	while (curShape!=0 && !err){
		if (curShape->name==shapeName||shapeName=="global"){
			if (curShape->wShape->val!=0){
				Shape tempShape;
				err=CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/dalpha";
				if (!err){
					tempShape.alpha->val*=1.01;
					delta=tempShape.alpha->val-curShape->alpha->val;
					err=tempShape.CalcPattern();
				}
				tempShape.CalcNumAtoms();
				delNAtoms+=(tempShape.nAtoms-curShape->nAtoms)*tempShape.wShape->val/delta;
				for(int i=0;i<Pattern::nInt;i++){
					for (int j=0;j<Pattern::nLambda;j++){
						newPattern[j].I[i]+=tempShape.wShape->val*(tempShape.diffPatt[j]->I[i]-curShape->diffPatt[j]->I[i])/delta;
					}
				}
			}
		}
		for(int i=0;i<Pattern::nInt;i++){
			for (int j=0;j<Pattern::nLambda;j++){
				oldPatt[j].I[i]+=curShape->wShape->val*curShape->diffPatt[j]->I[i];
			}
		}
		curShape=curShape->nextShape;
	}
	for (int j=0;j<Pattern::nLambda&&!err;j++){
		newPattern[j].ScalePattern(1/nAtoms);
		oldPatt[j].ScalePattern(delNAtoms/(nAtoms*nAtoms));
		for (int i=0;i<Pattern::nInt;i++){
			newPattern[j].I[i]-=oldPatt[j].I[i];
		}
		newPattern[j].ApplyMultFactors(j);
		newPattern[j].ApplyChebyshevScale();
		for (int i=0;i<Pattern::nInt&&!err;i++){
			outputGrad[i]+=Pattern::wLambda[j]*newPattern[j].I[i];
		}
	}
	delete [] newPattern;
	delete [] oldPatt;
	return err;
}
*/
int Compare::CalcDeriv_alpha(double *outputGrad, string shapeName)
{
	int err=0;
	Shape *curShape=firstShape;
	double delNAtoms=0;
	double delta=0.001;
	Pattern *newPattern=new Pattern[Pattern::nLambda];
	for (int j=0;j<Pattern::nLambda; j++){
		newPattern[j].AllocPattern();
	}

	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
	}
	printf("\nCalculating gradient wrt alpha\n");
	while (curShape!=0 && !err){
		if (curShape->name==shapeName||shapeName=="global"){
			if (curShape->wShape->GetVal()!=0){
				Shape tempShape;
				CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/dalpha";
				if (!err){
					if (tempShape.alpha->GetVal()!=0){
						delta=tempShape.alpha->GetVal()*.001;
					}
					tempShape.alpha->ChangeVal(tempShape.alpha->GetVal()+delta);
					err=tempShape.CalcPattern();
				}
				tempShape.CalcNumAtoms();
				delNAtoms+=tempShape.nAtoms*tempShape.wShape->GetVal();

				for(int i=0;i<Pattern::nInt;i++){
					for (int j=0;j<Pattern::nLambda;j++){
						newPattern[j].I[i]+=tempShape.wShape->GetVal()*tempShape.diffPatt[j]->I[i];
					}
				}

			}
		}
		else{
			for(int i=0;i<Pattern::nInt;i++){
				for (int j=0;j<Pattern::nLambda;j++){
					newPattern[j].I[i]+=curShape->wShape->GetVal()*curShape->diffPatt[j]->I[i];
				}
			}
			delNAtoms+=curShape->wShape->GetVal()*curShape->nAtoms;
		}


		curShape=curShape->nextShape;
	}
	for (int j=0;j<Pattern::nLambda&&!err;j++){
		newPattern[j].ScalePattern(1/delNAtoms);
		newPattern[j].ApplyWaveLenDependMultFactors(j);
		newPattern[j].ApplyGeomDependMultFactors();
		for (int i=0;i<Pattern::nInt;i++){
			newPattern[j].I[i]-=initCalcI[j]->I[i];
		}
		//newPattern[j].ApplyGeomDependMultFactors();
		newPattern[j].ApplyChebyshevScale();
		for (int i=0;i<Pattern::nInt&&!err;i++){
			outputGrad[i]+=Pattern::wLambda[j]*newPattern[j].I[i]/delta;
		}
	}
	delete [] newPattern;
	return err;
}
/*
int Compare::CalcDeriv_beta(double *outputGrad, string shapeName)
{
	int err=0, i;
	Shape *curShape=firstShape;

	for (i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
	}
	printf("\nCalculating gradient wrt beta\n");
	while (curShape!=0 && !err){
		if (curShape->wShape->val!=0){
			if (curShape->name==shapeName||shapeName=="global"){
				Shape tempShape;
				err=CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/dBeta";
				if (!err){
					tempShape.beta->val*=1.001;
					err=tempShape.CalcPattern();
					for (int j=0;j<Pattern::nLambda;j++){
						for (i=0;i<Pattern::nInt&&!err;i++){
							tempShape.diffPatt[j]->I[i]-=curShape->diffPatt[j]->I[i];
						}
						tempShape.diffPatt[j]->ApplyMultFactors(j);
						tempShape.diffPatt[j]->ScalePattern(1.0/nAtoms);
						err=tempShape.diffPatt[j]->ApplyChebyshevScale();
						for (i=0;i<Pattern::nInt&&!err;i++){
							outputGrad[i]+=Pattern::wLambda[j]*curShape->wShape->val*tempShape.diffPatt[j]->I[i]/(curShape->beta->val*.001);
						}
					}
					//}
				}
			}
		}
		curShape=curShape->nextShape;
	}
	return err;
}
*/
/*
int Compare::CalcDeriv_beta(double *outputGrad, string shapeName)
{
	int err=0;
	Shape *curShape=firstShape;
	double nA=0;
	double delta =0.001;
	Pattern *newPattern=new Pattern[Pattern::nLambda];
	try{
		for (int j=0;j<Pattern::nLambda; j++)
			newPattern[j].I=new double[Pattern::nInt];
	}
	catch(exception &e){
		cout<<"Exception in Compare::CalcDeriv_beta: "<<e.what();
		exit(1);
	}
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
		for (int j=0;j<Pattern::nLambda;j++){
			newPattern[j].I[i]=0;
		}
	}
	printf("\nCalculating gradient wrt beta\n");
	while (curShape!=0 && !err){
		if (curShape->name==shapeName||shapeName=="global"){
			if (curShape->wShape->val!=0){
				Shape tempShape;
				err=CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/dBeta";
				if (!err){
					tempShape.beta->val*=1.001;
					delta=tempShape.beta->val-curShape->beta->val;
					err=tempShape.CalcPattern();
				}
				nA+=tempShape.GetNumAtoms()*tempShape.wShape->val;
				for(int i=0;i<Pattern::nInt;i++){
					for (int j=0;j<Pattern::nLambda;j++){
						newPattern[j].I[i]+=tempShape.wShape->val*tempShape.diffPatt[j]->I[i];
					}
				}
			}
		}
		else{
			nA+=curShape->GetNumAtoms()*curShape->wShape->val;
			for(int i=0;i<Pattern::nInt;i++){
				for (int j=0;j<Pattern::nLambda;j++){
					newPattern[j].I[i]+=curShape->wShape->val*curShape->diffPatt[j]->I[i];
				}
			}
		}
		curShape=curShape->nextShape;
	}
	for (int j=0;j<Pattern::nLambda&&!err;j++){
		newPattern[j].ApplyMultFactors(j);
		newPattern[j].ScalePattern(1.0/nA);
		newPattern[j].ApplyChebyshevScale();
		for (int i=0;i<Pattern::nInt&&!err;i++){
			outputGrad[i]+=Pattern::wLambda[j]*(newPattern[j].I[i]-initCalcI[j]->I[i])/(delta);
		}
	}
	delete [] newPattern;
	return err;
}
*/
/*
int Compare::CalcDeriv_beta(double *outputGrad, string shapeName)
{
	int err=0;
	Shape *curShape=firstShape;
	double delNAtoms=0;
	double delta=0.001;
	Pattern *newPattern=new Pattern[Pattern::nLambda];
	Pattern *oldPatt=new Pattern[Pattern::nLambda];
	try{
		for (int j=0;j<Pattern::nLambda; j++){
			newPattern[j].I=new double[Pattern::nInt];
			oldPatt[j].I=new double [Pattern::nInt];
		}

	}
	catch(exception &e){
		cout<<"Exception in Compare::CalcDeriv_beta: "<<e.what();
		exit(1);
	}
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
		for (int j=0;j<Pattern::nLambda;j++){
			newPattern[j].I[i]=0;
			oldPatt[j].I[i]=0;
		}
	}
	printf("\nCalculating gradient wrt beta\n");
	while (curShape!=0 && !err){
		if (curShape->name==shapeName||shapeName=="global"){
			if (curShape->wShape->val!=0){
				Shape tempShape;
				err=CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/dbeta";
				if (!err){
					tempShape.beta->val*=1.01;
					delta=tempShape.beta->val-curShape->beta->val;
					err=tempShape.CalcPattern();
				}
				tempShape.CalcNumAtoms();
				delNAtoms+=(tempShape.nAtoms-curShape->nAtoms)*tempShape.wShape->val/delta;
				for(int i=0;i<Pattern::nInt;i++){
					for (int j=0;j<Pattern::nLambda;j++){
						newPattern[j].I[i]+=tempShape.wShape->val*(tempShape.diffPatt[j]->I[i]-curShape->diffPatt[j]->I[i])/delta;
					}
				}
			}
		}
		for(int i=0;i<Pattern::nInt;i++){
			for (int j=0;j<Pattern::nLambda;j++){
				oldPatt[j].I[i]+=curShape->wShape->val*curShape->diffPatt[j]->I[i];
			}
		}
		curShape=curShape->nextShape;
	}
	for (int j=0;j<Pattern::nLambda&&!err;j++){
		newPattern[j].ScalePattern(1/nAtoms);
		oldPatt[j].ScalePattern(delNAtoms/(nAtoms*nAtoms));
		for (int i=0;i<Pattern::nInt;i++){
			newPattern[j].I[i]-=oldPatt[j].I[i];
		}
		newPattern[j].ApplyMultFactors(j);
		newPattern[j].ApplyChebyshevScale();
		for (int i=0;i<Pattern::nInt&&!err;i++){
			outputGrad[i]+=Pattern::wLambda[j]*newPattern[j].I[i];
		}
	}
	delete [] newPattern;
	delete [] oldPatt;
	return err;
}
*/
int Compare::CalcDeriv_beta(double *outputGrad, string shapeName)
{
	int err=0;
	Shape *curShape=firstShape;
	double delNAtoms=0;
	double delta=0.001;
	Pattern *newPattern=new Pattern[Pattern::nLambda];
	for (int j=0;j<Pattern::nLambda; j++){
		newPattern[j].AllocPattern();
	}

	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
	}
	printf("\nCalculating gradient wrt beta\n");
	while (curShape!=0 && !err){
		if (curShape->name==shapeName||shapeName=="global"){
			if (curShape->wShape->GetVal()!=0){
				Shape tempShape;
				CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/dbeta";
				if (!err){
					if (tempShape.beta->GetVal()!=0){
						delta=tempShape.beta->GetVal()*.001;
					}
					tempShape.beta->ChangeVal(tempShape.beta->GetVal()+delta);
					err=tempShape.CalcPattern();
				}
				tempShape.CalcNumAtoms();
				delNAtoms+=tempShape.nAtoms*tempShape.wShape->GetVal();

				for(int i=0;i<Pattern::nInt;i++){
					for (int j=0;j<Pattern::nLambda;j++){
						newPattern[j].I[i]+=tempShape.wShape->GetVal()*tempShape.diffPatt[j]->I[i];
					}
				}

			}
		}
		else{
			for(int i=0;i<Pattern::nInt;i++){
				for (int j=0;j<Pattern::nLambda;j++){
					newPattern[j].I[i]+=curShape->wShape->GetVal()*curShape->diffPatt[j]->I[i];
				}
			}
			delNAtoms+=curShape->wShape->GetVal()*curShape->nAtoms;
		}


		curShape=curShape->nextShape;
	}
	for (int j=0;j<Pattern::nLambda&&!err;j++){
		newPattern[j].ScalePattern(1/delNAtoms);
		newPattern[j].ApplyWaveLenDependMultFactors(j);
		newPattern[j].ApplyGeomDependMultFactors();
		for (int i=0;i<Pattern::nInt;i++){
			newPattern[j].I[i]-=initCalcI[j]->I[i];
		}
		//newPattern[j].ApplyGeomDependMultFactors();
		newPattern[j].ApplyChebyshevScale();
		for (int i=0;i<Pattern::nInt&&!err;i++){
			outputGrad[i]+=Pattern::wLambda[j]*newPattern[j].I[i]/delta;
		}
	}
	delete [] newPattern;
	return err;
}
int Compare::CalcDeriv_deltaSize(double *outputGrad, string shapeName)
{
	int err=0;
	Shape *curShape=firstShape;
	double tempNAtoms=0;
	double delta=0.01;
	Pattern *newPattern=new Pattern[Pattern::nLambda];
	for (int j=0;j<Pattern::nLambda; j++){
		newPattern[j].AllocPattern();
	}

	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
	}
	if (shapeName=="global"){
		//TODO If Possible implement global deltaSize gradient calc
		cout<<"Calculation of GLOBAL deltaSize gradient is not supported.\n";
		err=-1;
	}
	else{
		printf("\nCalculating gradient wrt deltaSize\n");
	}
	while (curShape!=0 && !err){
		if (curShape->name==shapeName){
			if (curShape->wShape->GetVal()!=0){
				Shape tempShape;
				CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/dSize";
				int shell=tempShape.CalcDeltaShellSize();
				if (shell==-1){
					cout<<"Cannot Calc Gradient for "<<tempShape.name<<" size is out of range.\n";
					err=-1;
				}
				else if (shell==tempShape.maxShell){
					//Decrease by one shell
					delta=-tempShape.GetShapeShellL()*tempShape.a[0]->GetVal();
				}
				else{
					//Only increase by one shell
					delta=tempShape.GetShapeShellL()*tempShape.a[0]->GetVal();
				}
				if (!err){
					tempShape.deltaSize->ChangeVal(tempShape.deltaSize->GetVal()+delta);
					err=tempShape.CalcPattern();
					tempShape.CalcNumAtoms();
					tempNAtoms+=tempShape.nAtoms*tempShape.wShape->GetVal();

					for(int i=0;i<Pattern::nInt;i++){
						for (int j=0;j<Pattern::nLambda;j++){
							newPattern[j].I[i]+=tempShape.wShape->GetVal()*tempShape.diffPatt[j]->I[i];
						}
					}
				}
			}
		}
		else{
			for(int i=0;i<Pattern::nInt;i++){
				for (int j=0;j<Pattern::nLambda;j++){
					newPattern[j].I[i]+=curShape->wShape->GetVal()*curShape->diffPatt[j]->I[i];
				}
			}
			tempNAtoms+=curShape->wShape->GetVal()*curShape->nAtoms;
		}


		curShape=curShape->nextShape;
	}
	for (int j=0;j<Pattern::nLambda&&!err;j++){
		newPattern[j].ScalePattern(1.0/tempNAtoms);
		newPattern[j].ApplyWaveLenDependMultFactors(j);
		newPattern[j].ApplyGeomDependMultFactors();
		for (int i=0;i<Pattern::nInt;i++){
			newPattern[j].I[i]-=initCalcI[j]->I[i];
		}
		//newPattern[j].ApplyGeomDependMultFactors();
		newPattern[j].ApplyChebyshevScale();
		for (int i=0;i<Pattern::nInt&&!err;i++){
			outputGrad[i]+=Pattern::wLambda[j]*newPattern[j].I[i]/delta;
		}
	}
	delete [] newPattern;
	return err;
}
//Could make a routine to calc derivative by considering the derivative of sinc
int Compare::CalcDeriv_a(double *outputGrad, string shapeName)
{
	int err=0;//, j, k;
	Shape *curShape=firstShape;
	Pattern **tempI;
	try{
		tempI=new Pattern* [Pattern::nLambda];
		for (int i=0;i<Pattern::nLambda;i++)
			tempI[i]=new Pattern ("dI/da");
	}
	catch (exception& e){
		printf("\nException found in CalcDeriv_a alloc: %s\n", e.what());
		return -1;
	}
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
	}
	printf("\nCalculating gradient wrt a\n");

	while (curShape!=0 && !err){
		if (curShape->wShape->GetVal()!=0){
			if (curShape->name==shapeName||shapeName=="global"){
				double del=.0001;
				if (curShape->a[0]->GetVal()!=0){
					del=curShape->a[0]->GetVal()*.001;
				}
				Shape tempShape;
				CopyShapeInfo(curShape, &tempShape);
				tempShape.a[0]->ChangeVal(curShape->a[0]->GetVal()+del);
				if (curShape->nPreCalcPatts>0){
					tempShape.InitLatticeParam();
					tempShape.LoadPatterns();
				}
				if (curShape->nRelaxShells>0){
					//Save calculation time by only rescaling the surface relax distances
					Distance tempDist;
					curShape->CopyDistance(&tempDist, curShape->surfRelaxDist);
					tempDist.ScaleDistance(tempShape.a[0]->GetVal()/curShape->a[0]->GetVal());
					for (int i=0;i<Pattern::nLambda&&!err;i++){
						if (tempShape.preCalcPatts==0){
							// Function to use if there are no pre-calculated patterns
							tempI[i]->CalcPattern(curShape->firstDist, curShape->maxShell, curShape->wShell, tempShape.a[0]->GetVal(), &tempDist, i);
						}
						else{
							// Function to use if there are pre-calculated patterns
							tempI[i]->CalcPattern(curShape->firstDist, curShape->maxShell, curShape->p, tempShape.a_D, curShape->nRelaxShells, &tempDist, i, tempShape.preCalcPatts);
						}
					}
				}
				else{
					//Functions if there are no surface relax distances calculated
					for (int i=0;i<Pattern::nLambda&&!err;i++){
						if (tempShape.preCalcPatts==0){
							tempI[i]->CalcPattern(curShape->firstDist, curShape->maxShell,  curShape->wShell, tempShape.a[0]->GetVal(), 0,i);
						}
						else{
							tempI[i]->CalcPattern(curShape->firstDist, curShape->maxShell, curShape->p, tempShape.a_D, curShape->nRelaxShells, 0, i, tempShape.preCalcPatts);
						}
					}

				}
				for (int j=0;j<Pattern::nLambda;j++){
					for (int i=0;i<Pattern::nInt&&!err;i++){
						tempI[j]->I[i]-=curShape->diffPatt[j]->I[i];
					}
					tempI[j]->ApplyWaveLenDependMultFactors(j);
					tempI[j]->ScalePattern(1.0/nAtoms);
					tempI[j]->ApplyGeomDependMultFactors();
					tempI[j]->ApplyChebyshevScale();
					for (int i=0;i<Pattern::nInt&&!err;i++){
						outputGrad[i]+=Pattern::wLambda[j]*curShape->wShape->GetVal()*tempI[j]->I[i]/(del);
					}
				}
			}
		}
		curShape=curShape->nextShape;
	}
	/*

	//Analytical Derivative (under construction... causes an error in the cholesky Decon if initial lambda val is not good)
		try{
			if (tempI.I==0)
				tempI.I= new double[tempI.nInt];
		}
		catch (exception& e){
			printf("\nException found in CalcDeriv_a: %s\n", e.what());
			return -1;
		}

		//printf("Calculating pattern for %s\n", name.c_str());
		for (i=0;i<tempI.nInt;i++)
			tempI.I[i]=0;
		while (curShape!=0 && !err){
			if(curShape->wShape->val!=0){
				curDist=curShape->firstDist;
				for (i=0;(i<curShape->maxShell)&&(curShape->wShell[i]!=0);i++){
					for(j=0;j<tempI.nInt;j++){
						for (k=0;k<curDist->nDist;k++)
							{

								tempI.I[j]+=curDist->mult[k]*curShape->wShell[i]*Pattern::offsetq[j]*curDist->dist[k]*Cosc(Pattern::offsetq[j]*curShape->a->val*curDist->dist[k]);
							}
							//tempI.I[j]+=curShape->wShell[i]*curDist->mult[0];
					}
					curDist=curDist->nextDist;
				}
				// Add intensity from the surface relax distance list
				if (curShape->surfRelaxDist!=0){
					curDist=curShape->surfRelaxDist;
					for(j=0;j<nInt;j++){
						for (k=0;k<curDist->nDist;k++)
							{
								//I[j]+=curDist->mult[k]*Sinc(offsetq[j]*a*curDist->dist[k]);
								tempI.I[j]+=curDist->mult[k]*Pattern::offsetq[j]*curDist->dist[k]*Cosc(Pattern::offsetq[j]*curShape->a->val*curDist->dist[k]);
							}
							//tempI.I[j]+=curDist->mult[0];
					}
				}
				for (i=0;i<Pattern::nInt&&!err;i++){
					outputGrad[i]=curShape->wShape->val*tempI.I[i];
				}
			}
			curShape=curShape->nextShape;
		}
		*/
	//Apply Debye Waller and Atomic Scattering factor
	/*
	if (!err){
		for (i=0;i<Pattern::nInt; i++){
			tempI.I[i]=1.0f;
		}
		tempI.ApplyAtomicScatFactor();
		tempI.ApplyDWFactor();
		for (i=0; i<Pattern::nInt; i++){
			outputGrad[i]*=Pattern::scale->val*tempI.I[i];
		}
	}
	*/
	for (int i=0;i<Pattern::nLambda;i++)
		delete tempI[i];
	delete [] tempI;
	return err;
}
int Compare::CalcDeriv_a_D(double *outputGrad, int n, string shapeName)
{
	int err=0;
	Shape *curShape=firstShape;

	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
	}
	cout<<"\nCalculating gradient wrt a_"<<n<<" Poly"<<endl;
	while (curShape!=0 && !err){
		if (curShape->name==shapeName||shapeName=="global"){
			if (curShape->wShape->GetVal()!=0){
				Shape tempShape;
				CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/da";
				double del=0.00001;
				if (!err){
					if (tempShape.a[n]->GetVal()!=0){
						del=tempShape.a[n]->GetVal()*.001;
					}
					tempShape.a[n]->ChangeVal(tempShape.a[n]->GetVal()+del);
					err=tempShape.CalcPattern();
				}
				for (int j=0;j<Pattern::nLambda&&!err;j++){
					for (int i=0;i<Pattern::nInt&&!err;i++){
						tempShape.diffPatt[j]->I[i]-=curShape->diffPatt[j]->I[i];
					}
					tempShape.diffPatt[j]->ApplyWaveLenDependMultFactors(j);
					tempShape.diffPatt[j]->ScalePattern(1.0/nAtoms);
					tempShape.diffPatt[j]->ApplyGeomDependMultFactors();
					tempShape.diffPatt[j]->ApplyChebyshevScale();
					for (int i=0;i<Pattern::nInt&&!err;i++){
						outputGrad[i]+=Pattern::wLambda[j]*curShape->wShape->GetVal()*tempShape.diffPatt[j]->I[i]/(del);
					}
				}
			}
		}
		curShape=curShape->nextShape;
	}

	return err;
}
/*
int Compare::CalcDeriv_f(double *outputGrad, string shapeName)
{
	int err=0, i;
	Shape *curShape=firstShape;

	for (i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
	}
	printf("\nCalculating gradient wrt f\n");
	while (curShape!=0 && !err){
		if (curShape->name==shapeName||shapeName=="global"){
			if (curShape->wShape->val!=0){
				if (curShape->kappa->val>0){
					Shape tempShape;
					err=CopyShapeInfo(curShape, &tempShape);
					tempShape.name="dI/df";
					if (!err){
						tempShape.f->val*=(1+5*SRPREC);
						err=tempShape.CalcPattern();

						for (int j=0;j<Pattern::nLambda;j++){
							for (i=0;i<Pattern::nInt&&!err;i++){
								tempShape.diffPatt[j]->I[i]-=curShape->diffPatt[j]->I[i];
							}
							tempShape.diffPatt[j]->ApplyMultFactors(j);
							tempShape.diffPatt[j]->ScalePattern(1.0/nAtoms);
							err=tempShape.diffPatt[j]->ApplyChebyshevScale();
							for (i=0;i<Pattern::nInt&&!err;i++){
								outputGrad[i]+=Pattern::wLambda[j]*curShape->wShape->val*tempShape.diffPatt[j]->I[i]/(curShape->f->val*5*SRPREC);
							}
						}
						//}
					}
				}
				else{
					printf("\nError: Cannot calculate derivative of f when kappa=0\n");
					err=-1;
				}
			}
		}
		curShape=curShape->nextShape;
	}
	return err;
}
*/
int Compare::CalcDeriv_f(double *outputGrad, string shapeName)
{
	int err=0;
	Shape *curShape=firstShape;
	Pattern *newPattern=new Pattern[Pattern::nLambda];
	try{
		for (int j=0;j<Pattern::nLambda; j++)
			newPattern[j].I=new double[Pattern::nInt];
	}
	catch(exception &e){
		cout<<"Exception in Compare::CalcDeriv_f: "<<e.what();
		exit(1);
	}
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
		for (int j=0;j<Pattern::nLambda;j++){
			newPattern[j].I[i]=0;
		}
	}
	printf("\nCalculating gradient wrt f\n");
	while (curShape!=0 && !err){
		if (curShape->name==shapeName||shapeName=="global"){
			if (curShape->wShape->GetVal()!=0){
				Shape tempShape;
				CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/df";
				double delta;
				if (curShape->f->GetVal()+2*SRPREC>curShape->f->max){
					delta=-2*SRPREC;
				}
				else{
					delta=2*SRPREC;
				}
				tempShape.f->ChangeVal(tempShape.f->GetVal()+delta);

				err=tempShape.CalcPattern();
				for(int i=0;i<Pattern::nInt;i++){
					for (int j=0;j<Pattern::nLambda;j++){
						newPattern[j].I[i]+=tempShape.wShape->GetVal()*(tempShape.diffPatt[j]->I[i]-curShape->diffPatt[j]->I[i])/delta;
					}
				}
			}
		}
		curShape=curShape->nextShape;
	}
	for (int j=0;j<Pattern::nLambda&&!err;j++){
		newPattern[j].ApplyWaveLenDependMultFactors(j);
		newPattern[j].ScalePattern(1.0/nAtoms);
		newPattern[j].ApplyGeomDependMultFactors();
		newPattern[j].ApplyChebyshevScale();
		for (int i=0;i<Pattern::nInt&&!err;i++){
			outputGrad[i]+=Pattern::wLambda[j]*newPattern[j].I[i];
		}
	}
	delete [] newPattern;
	return err;
}
/*
int Compare::CalcDeriv_kappa(double *outputGrad, string shapeName)
{
	int err=0, i;
	Shape *curShape=firstShape;

	for (i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
	}
	printf("\nCalculating gradient wrt kappa\n");
	while (curShape!=0 && !err){
		if (curShape->name==shapeName||shapeName=="global"){
			if (curShape->wShape->val!=0){
				Shape tempShape;
				err=CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/dKappa";
				if (!err){
					tempShape.kappa->val+=(5*SRPREC);
					err=tempShape.CalcPattern();
					for (int j=0;j<Pattern::nLambda&&!err;j++){
						for (i=0;i<Pattern::nInt&&!err;i++){
								tempShape.diffPatt[j]->I[i]-=curShape->diffPatt[j]->I[i];
						}
						tempShape.diffPatt[j]->ApplyMultFactors(j);
						tempShape.diffPatt[j]->ScalePattern(1.0/nAtoms);
						err=tempShape.diffPatt[j]->ApplyChebyshevScale();
						for (i=0;i<Pattern::nInt&&!err;i++){
							outputGrad[i]+=Pattern::wLambda[j]*curShape->wShape->val*tempShape.diffPatt[j]->I[i]/(5*SRPREC);
						}
					}

				}

			}
		}
		curShape=curShape->nextShape;
	}
	return err;
}
*/
int Compare::CalcDeriv_kappa(double *outputGrad, string shapeName)
{
	int err=0;
	Shape *curShape=firstShape;
	Pattern *newPattern=new Pattern[Pattern::nLambda];
	try{
		for (int j=0;j<Pattern::nLambda; j++)
			newPattern[j].I=new double[Pattern::nInt];
	}
	catch(exception &e){
		cout<<"Exception in Compare::CalcDeriv_kappa: "<<e.what();
		exit(1);
	}
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
		for (int j=0;j<Pattern::nLambda;j++){
			newPattern[j].I[i]=0;
		}
	}
	printf("\nCalculating gradient wrt kappa\n");
	while (curShape!=0 && !err){
		if (curShape->name==shapeName||shapeName=="global"){
			if (curShape->wShape->GetVal()!=0){
				Shape tempShape;
				CopyShapeInfo(curShape, &tempShape);
				tempShape.name="dI/dKappa";
				double delta;
				if (curShape->kappa->GetVal()+5*SRPREC>curShape->kappa->max){
					delta=-10*SRPREC;
				}
				else{
					delta=10*SRPREC;
				}
				tempShape.kappa->ChangeVal(tempShape.kappa->GetVal()+delta);
				err=tempShape.CalcPattern();
				for(int i=0;i<Pattern::nInt;i++){
					for (int j=0;j<Pattern::nLambda;j++){
						newPattern[j].I[i]+=tempShape.wShape->GetVal()*(tempShape.diffPatt[j]->I[i]-curShape->diffPatt[j]->I[i])/delta;
					}
				}
			}
		}
		curShape=curShape->nextShape;
	}
	for (int j=0;j<Pattern::nLambda&&!err;j++){
		newPattern[j].ApplyWaveLenDependMultFactors(j);
		newPattern[j].ScalePattern(1.0/nAtoms);
		newPattern[j].ApplyGeomDependMultFactors();
		newPattern[j].ApplyChebyshevScale();
		for (int i=0;i<Pattern::nInt&&!err;i++){
			outputGrad[i]+=Pattern::wLambda[j]*newPattern[j].I[i];
		}
	}
	delete [] newPattern;
	return err;
}
/*
int Compare::CalcDeriv_wShape(double *outputGrad, string shapeName)
{
	int err=0;
	int nShapes=0;
	float step=.01;
	Shape *curShape=firstShape;
	Pattern **tempI;
	try{
		tempI=new Pattern* [Pattern::nLambda];
		for (int i=0;i<Pattern::nLambda;i++){
			tempI[i]=new Pattern ("dI/dwShape");
			if (tempI[i]->I==0) tempI[i]->I= new double [Pattern::nInt];
		}
	}
	catch (exception& e){
		printf("\nException found in CalcDeriv_wShape alloc: %s\n", e.what());
		return -1;
	}
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
	}
	curShape=firstShape;
	while(curShape!=0){
		if (curShape->wShape->val!=0){
			nShapes++;
			curShape=curShape->nextShape;
		}
	}
	printf("\nCalculating gradient wrt wShape\n");
	for (int j=0;j<Pattern::nLambda;j++){
		curShape=firstShape;
		for (int i=0;i<Pattern::nInt; i++){
			tempI[j]->I[i]=0;
		}
		while (curShape!=0 && !err){
			if (curShape->wShape->val!=0){
				if (curShape->name==shapeName){
					for (int i=0;i<Pattern::nInt&&!err;i++){
						tempI[j]->I[i]+=(curShape->wShape->val+(nShapes-1)*step)*curShape->diffPatt[j]->I[i];
					}
				}
				else{
					for (int i=0;i<Pattern::nInt&&!err;i++){
						tempI[j]->I[i]+=(curShape->wShape->val-step)*curShape->diffPatt[j]->I[i];
					}
				}
				/* Analytical Derivative
				if (curShape==firstShape){
					for(i=0;i<Pattern::nInt&&!err;i++){
						tempI[j]->I[i]=-curShape->diffPatt[j]->I[i];
					}
				}
				if (curShape->shape==shapeName){
					for (i=0;i<Pattern::nInt&&!err;i++){
						tempI[j]->I[i]+=curShape->diffPatt[j]->I[i];
					}
				}
				*/
/*
			}
			curShape=curShape->nextShape;
		}
		/*tempI[j]->ApplyAtomicScatFactor(j);
		tempI[j]->ApplyDWFactor(j);
		tempI[j]->ApplyPolarizationFactor();*/
/*
		tempI[j]->ApplyMultFactors(j);
		//err=tempI[j]->ApplyChebyshevScale();
		for (int i=0;i<Pattern::nInt;i++)
			tempI[j]->I[i]=Pattern::wLambda[j]*(tempI[j]->I[i]-initCalcI[j]->I[i])/(step*(nShapes-1));
		err=tempI[j]->ApplyChebyshevScale();
		if (!err){
			for (int i=0;i<Pattern::nInt;i++){
				outputGrad[i]+=tempI[j]->I[i];
			}
		}
	}
	for (int i=0;i<Pattern::nLambda;i++)
		delete tempI[i];
	delete [] tempI;
	return err;
}
*/
/*
int Compare::CalcDeriv_wShape(double *outputGrad, string shapeName)
{
	int err=0;
	int nShapes=0;
	float step=.01;
	double nA=0;
	Shape *curShape=firstShape;
	Pattern **tempI;
	try{
		tempI=new Pattern* [Pattern::nLambda];
		for (int i=0;i<Pattern::nLambda;i++){
			tempI[i]=new Pattern ("dI/dwShape");
			if (tempI[i]->I==0) tempI[i]->I= new double [Pattern::nInt];
		}
	}
	catch (exception& e){
		printf("\nException found in CalcDeriv_wShape alloc: %s\n", e.what());
		return -1;
	}
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
	}
	curShape=firstShape;
	while(curShape!=0){
		if (curShape->wShape->val!=0){
			nShapes++;
			curShape=curShape->nextShape;
		}
	}
	printf("\nCalculating gradient wrt wShape\n");
	for (int j=0;j<Pattern::nLambda;j++){
		for (int i=0;i<Pattern::nInt; i++){
			tempI[j]->I[i]=0;
		}
	}
	curShape=firstShape;

	while (curShape!=0 && !err){
		if (curShape->wShape->val!=0){
			if (curShape->name==shapeName){
				for (int j=0;j<Pattern::nLambda;j++){
					for (int i=0;i<Pattern::nInt&&!err;i++){
						tempI[j]->I[i]+=(curShape->wShape->val+(nShapes-1)*step)*curShape->diffPatt[j]->I[i];
					}
				}
				nA+=(curShape->wShape->val+(nShapes-1)*step)*curShape->GetNumAtoms();
			}
			else{
				for (int j=0;j<Pattern::nLambda;j++){
					for (int i=0;i<Pattern::nInt&&!err;i++){
						tempI[j]->I[i]+=(curShape->wShape->val-step)*curShape->diffPatt[j]->I[i];
					}
				}
				nA+=(curShape->wShape->val-step)*curShape->GetNumAtoms();
			}
		}
		curShape=curShape->nextShape;
	}
	for (int j=0;j<Pattern::nLambda;j++){
		tempI[j]->ApplyMultFactors(j);
		tempI[j]->ScalePattern(1.0/nA);
		tempI[j]->ApplyChebyshevScale();
		for (int i=0;i<Pattern::nInt;i++)
			outputGrad[i]=Pattern::wLambda[j]*(tempI[j]->I[i]-initCalcI[j]->I[i])/(step*(nShapes-1));
	}
	for (int i=0;i<Pattern::nLambda;i++)
		delete tempI[i];
	delete [] tempI;
	return err;
}
*/
int Compare::CalcDeriv_wShape(double *outputGrad, string shapeName)
{
	int err=0;
	Shape *curShape=firstShape;
	double delta=0.001;
	double derivNAtoms=0;
	Pattern *derivPatt=new Pattern[Pattern::nLambda];
	for (int j=0;j<Pattern::nLambda; j++){
		derivPatt[j].AllocPattern();
	}
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
	}
	printf("\nCalculating gradient wrt wShape\n");
	double newW=0, otherWs=0;
	while (curShape!=0 && !err){
		if (curShape->name==shapeName){
			newW=curShape->wShape->GetVal()*1.001;
		}
		else{
			otherWs+=curShape->wShape->GetVal();
		}
		curShape=curShape->nextShape;
	}
	double sumNewWs=newW+otherWs;
	newW/=sumNewWs;
	curShape=firstShape;
	while (curShape!=0 && !err){
		if (curShape->name==shapeName){
			if (curShape->wShape->GetVal()!=0){
				delta=curShape->wShape->GetVal()-newW;
				for(int i=0;i<Pattern::nInt;i++){
					for (int j=0;j<Pattern::nLambda;j++){
						derivPatt[j].I[i]+=newW*curShape->diffPatt[j]->I[i];
					}
				}
				derivNAtoms+=newW*curShape->nAtoms;
			}
		}
		else{
			for(int i=0; i<Pattern::nInt;i++){
				for(int j=0;j<Pattern::nLambda; j++){
					derivPatt[j].I[i]+=curShape->wShape->GetVal()*curShape->diffPatt[j]->I[i]/sumNewWs;
				}
			}
			derivNAtoms+=curShape->wShape->GetVal()*curShape->nAtoms;
		}
		curShape=curShape->nextShape;
	}
	Pattern totCurPatt, totDerivPatt;
	totCurPatt.AllocPattern();
	totDerivPatt.AllocPattern();
	for (int j=0;j<Pattern::nLambda;j++){
		derivPatt[j].ApplyWaveLenDependMultFactors(j);
		derivPatt[j].ApplyGeomDependMultFactors();
	}
	for(int i=0; i<Pattern::nInt;i++){
		for(int j=0;j<Pattern::nLambda; j++){
			totCurPatt.I[i]+=Pattern::wLambda[j]*initCalcI[j]->I[i];
			totDerivPatt.I[i]+=Pattern::wLambda[j]*derivPatt[j].I[i];
		}
	}
	totDerivPatt.ScalePattern(1.0/derivNAtoms);
	//totDerivPatt.ApplyGeomDependMultFactors();
	//totCurPatt.ApplyGeomDependMultFactors();
	totCurPatt.ApplyChebyshevScale();
	totDerivPatt.ApplyChebyshevScale();

	for (int i=0;i<Pattern::nInt&&!err;i++){
		outputGrad[i]+=(totDerivPatt.I[i]-totCurPatt.I[i])/delta;
	}
	delete [] derivPatt;
	return err;
}
int Compare::CalcDeriv_smplDisplace(double *outputGrad)
{
	int err=0, i,j;//, j,k ;
	double saveSmplDisp;
	Shape *curShape=firstShape;
	Pattern tempI("dI/dsmplDisplace");
	for (i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
	}
	printf("\nCalculating gradient wrt sample displacement\n");
	if (Pattern::sampleDisplace->GetVal()==0){
		err=-1;
		printf("\nError: cannot calculate numerical derivative (divide by zero)\n");
	}
	else{

		//Numerical Derivative
		saveSmplDisp=tempI.sampleDisplace->GetVal();
		Pattern::sampleDisplace->ChangeVal(tempI.sampleDisplace->GetVal()+saveSmplDisp*.001);
		if (!err&&Pattern::sampleDisplace!=0){
			tempI.ApplySampleDisplacement();
			if ((Pattern::ax!=0||Pattern::bx!=0||Pattern::cx!=0||Pattern::dx!=0||Pattern::ex!=0)&&(!Pattern::synchThetaShift)){
				tempI.ApplySynchThetaShift();
			}
		}
		//Change back to original value
		Pattern::sampleDisplace->ChangeVal(saveSmplDisp);

		Shape tempShape;
		CopyShapeInfo(curShape, &tempShape);
		if (curShape->nPreCalcPatts>0){
			tempShape.InitLatticeParam();
			tempShape.LoadPatterns();
		}
		while (curShape!=0 && !err){
			if (curShape->wShape->GetVal()!=0){
				for (j=0;j<Pattern::nLambda&&!err;j++)
				{
					if (curShape->nPreCalcPatts==0){
						if (curShape->nParams_a>1){
							tempI.CalcPattern(curShape->firstDist, curShape->maxShell, curShape->p, curShape->a_D, curShape->nRelaxShells, curShape->surfRelaxDist, j);
						}
						else{
							tempI.CalcPattern(curShape->firstDist, curShape->maxShell, curShape->wShell, curShape->a[0]->GetVal(), curShape->surfRelaxDist, j);
						}
					}
					else{
						//Routine to use when precalculated patterns are loaded.
						tempI.CalcPattern(curShape->firstDist, curShape->endSize, curShape->p, curShape->a_D, curShape->nRelaxShells, curShape->surfRelaxDist, j, tempShape.preCalcPatts);
					}
					for (i=0;i<Pattern::nInt&&!err;i++){
						tempI.I[i]-=curShape->diffPatt[j]->I[i];
					}
					/*tempI.ApplyAtomicScatFactor(j);
					tempI.ApplyDWFactor(j);
					tempI.ApplyPolarizationFactor();
					*/
					tempI.ApplyWaveLenDependMultFactors(j);
					tempI.ApplyGeomDependMultFactors();
					tempI.ScalePattern(1.0/nAtoms);
					if (tempI.nChebScale){
						tempI.ApplyChebyshevScale();
					}
					//tempI.ApplyGeomDependMultFactors();
					for (i=0;i<Pattern::nInt&&!err;i++){
						outputGrad[i]+=Pattern::wLambda[j]*curShape->wShape->GetVal()*tempI.I[i]/(saveSmplDisp*.001);

					}
				}
			}
			curShape=curShape->nextShape;
		}
		/*
		//Analytical Derivative
		try{
			if (tempI.I==0)
				tempI.I= new double[tempI.nInt];
		}
		catch (exception& e){
			printf("\nException found in CalcDeriv_smplDisplace: %s\n", e.what());
			return -1;
		}

		//printf("Calculating pattern for %s\n", name.c_str());
		for (i=0;i<tempI.nInt;i++)
			tempI.I[i]=0;
		while (curShape!=0 && !err){
			if(curShape->wShape->val!=0){
				curDist=curShape->firstDist;
				for (i=0;(i<curShape->maxShell)&&(curShape->wShell[i]!=0);i++){
					for(j=0;j<tempI.nInt;j++){
						for (k=0;k<curDist->nDist;k++)
							{

								tempI.I[j]+=curDist->mult[k]*curShape->wShell[i]*curShape->a->val*curDist->dist[k]*Cosc(Pattern::offsetq[j]*curShape->a->val*curDist->dist[k]);
							}
							//tempI.I[j]+=curShape->wShell[i]*curDist->mult[0];
					}
					curDist=curDist->nextDist;
				}
				// Add intensity from the surface relax distance list
				if (curShape->surfRelaxDist!=0){
					curDist=curShape->surfRelaxDist;
					for(j=0;j<nInt;j++){
						for (k=0;k<curDist->nDist;k++)
							{
								//I[j]+=curDist->mult[k]*Sinc(offsetq[j]*a*curDist->dist[k]);
								tempI.I[j]+=curDist->mult[k]*curShape->a->val*curDist->dist[k]*Cosc(Pattern::offsetq[j]*curShape->a->val*curDist->dist[k]);
							}
							//tempI.I[j]+=curDist->mult[0];
					}
				}
				for (i=0;i<Pattern::nInt&&!err;i++){
					outputGrad[i]=curShape->wShape->val*tempI.I[i]*(-8*PI*PI*cos(Pattern::offset2Theta[i]*PI/360)*cos(Pattern::theta2[i]*PI/360))/(Pattern::lambda*360*Pattern::gonioRad);
				}
			}
			curShape=curShape->nextShape;
		}
		*/
		//Atomic scattering factor, Debye Waller factor
			/*
		if (!err) {
			for (i=0;i<Pattern::nInt; i++){
				tempI.I[i]=1.0f;
			}
			tempI.ApplyAtomicScatFactor();
			tempI.ApplyDWFactor();
			for (i=0; i<Pattern::nInt; i++){
				outputGrad[i]*=Pattern::scale->val*tempI.I[i];
			}
		}
		*/
	}
	return err;
}
void Compare::CalcDeriv_scale(double *outputGrad)
{
	Pattern derivPatt;
	derivPatt.AllocPattern();
	printf("\nCalculating deriv wrt scale\n");
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=0;
		for (int j=0; j<Pattern::nLambda;j++){
			derivPatt.I[i]+=Pattern::wLambda[j]*initCalcI[j]->I[i];
		}
	}
	//derivPatt.ApplyGeomDependMultFactors();
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=derivPatt.I[i];
	}
}
void Compare::CalcDeriv_ChebScale(double *outputGrad, int k)
{
	Pattern derivPatt;
	derivPatt.AllocPattern();
	cout<<endl<<"Calculating deriv wrt ChebScale "<<k<<endl;
	for (int i=0;i<Pattern::nInt;i++){
		for (int j=0; j<Pattern::nLambda; j++){
			derivPatt.I[i]+=Pattern::wLambda[j]*initCalcI[j]->I[i]*Pattern::chebPolys[k][i];
		}
	}
	//derivPatt.ApplyGeomDependMultFactors();
	for (int i=0;i<Pattern::nInt;i++){
		outputGrad[i]=derivPatt.I[i];
	}
}
void Compare::CalcDeriv_gaussA(double *outputGrad, int n)
{
	printf("Calculating deriv wrt Background Gaussian Amplitude\n");
	for (int i=0; i<Pattern::nInt;i++){
		outputGrad[i]=Pattern::gaussBG[n][i]/Pattern::gaussA[n]->GetVal();
	}
}
void Compare::CalcDeriv_gaussMu(double *outputGrad, int n)
{
	printf("Calculating deriv wrt Background Gaussian Mu\n");
	for (int i=0; i<Pattern::nInt;i++){
		outputGrad[i]=Pattern::gaussBG[n][i]*((Pattern::theta2[i]-Pattern::gaussMu[n]->GetVal())/(Pattern::gaussSig[n]->GetVal()*Pattern::gaussSig[n]->GetVal()));
	}
}
void Compare::CalcDeriv_gaussSig(double *outputGrad, int n)
{
	printf("Calculating deriv wrt Background Gaussian Sigma\n");
	for (int i=0; i<Pattern::nInt;i++){
		outputGrad[i]=Pattern::gaussBG[n][i]*((Pattern::theta2[i]-Pattern::gaussMu[n]->GetVal())*(Pattern::theta2[i]-Pattern::gaussMu[n]->GetVal())/(Pattern::gaussSig[n]->GetVal()*Pattern::gaussSig[n]->GetVal()*Pattern::gaussSig[n]->GetVal()));
	}
}
void Compare::CalcDeriv_expAmp(double *outputGrad)
{
	cout<<endl<<"Calculating deriv wrt Exp Decay Background Amplitude"<<endl;
	double del=0;
	if (Pattern::delExp){
		del=Pattern::delExp->GetVal();
	}
	for (int i=0; i<Pattern::nInt;i++){
		outputGrad[i]=exp(-(Pattern::theta2[i]-del)/Pattern::kExp->GetVal());
	}
}
void Compare::CalcDeriv_expKappa(double *outputGrad)
{
	cout<<endl<<"Calculating deriv wrt Exp Decay Background Kappa"<<endl;
	double del=0;
	if (Pattern::delExp){
		del=Pattern::delExp->GetVal();
	}
	for (int i=0; i<Pattern::nInt;i++){
		outputGrad[i]=((Pattern::theta2[i]-del)*Pattern::AExp->GetVal()*exp(-(Pattern::theta2[i]-del)/Pattern::kExp->GetVal()))/(Pattern::kExp->GetVal()*Pattern::kExp->GetVal());
	}
}
void Compare::CalcDeriv_expDel(double *outputGrad)
{
	cout<<endl<<"Calculating deriv wrt Exp Decay Background Shift"<<endl;
	for (int i=0; i<Pattern::nInt;i++){
		outputGrad[i]=Pattern::AExp->GetVal()*exp(-(Pattern::theta2[i]-Pattern::delExp->GetVal())/Pattern::kExp->GetVal())/Pattern::kExp->GetVal();
	}
}
void Compare::CopyShapeInfo(Shape *from, Shape *to)
{
	try{
		to->nParams_a=from->nParams_a;
		to->a=new Param* [from->nParams_a];
		for (int i=0; i<from->nParams_a;i++){
			to->a[i] = new Param();
		}
		if (from->distribution=="LogNorm"){
			to->mu=new Param();
			to->sigma=new Param();
		}
		else if (from->distribution=="Gamma"){
			to->alpha=new Param();
			to->beta=new Param();
		}
		else if (from->distribution=="Delta"){
			to->deltaSize=new Param();
		}
		to->wShape=new Param();
		to->f=new Param();
		to->kappa=new Param();
		if (from->n_R!=0){
			to->n_R=new Param();
		}
		if (from->nRelaxShells>0){
			to->surfRelaxDist=new Distance();
		}
		if (from->preCalcPatts!=0){
			to->nPreCalcPatts=from->nPreCalcPatts;
			to->preCalcPatts= new Pattern*[to->nPreCalcPatts];
		}
	}
	catch (exception& e){
		printf("\nException found in CopyShapeInfo: %s\n", e.what());
		exit(1);
	}
	for (int i=0;i<from->nParams_a;i++){
		to->a[i]->CopyParam(from->a[i],false);
	}
	to->distribution=from->distribution;
	to->firstDist=from->firstDist;
	to->firstPos=from->firstPos;
	to->maxShell=from->maxShell;
	to->type=from->type;
	to->name=from->name;
	//to->shape=from->shape;
	if (from->distribution=="LogNorm"){
		to->mu->CopyParam(from->mu,false);
		to->sigma->CopyParam(from->sigma,false);
	}
	else if (from->distribution=="Gamma"){
		to->alpha->CopyParam(from->alpha,false);
		to->beta->CopyParam(from->beta,false);
	}
	else if (from->distribution=="Delta"){
		to->deltaSize->CopyParam(from->deltaSize,false);
	}
	to->wShape->CopyParam(from->wShape,false);
	to->f->CopyParam(from->f,false);
	to->kappa->CopyParam(from->kappa,false);
	if (to->n_R!=0){
		to->n_R->CopyParam(from->n_R,false);
	}
	if (from->nRelaxShells>0){
		to->nRelaxShells=from->nRelaxShells;
	}
	if (to->nPreCalcPatts>0){
		for (int i=0;i<to->nPreCalcPatts;i++){
			if(from->preCalcPatts[i]!=0){
				to->preCalcPatts[i]=new Pattern [Pattern::nLambda];
				for (int j=0;j<Pattern::nLambda;j++){
					to->preCalcPatts[i][j].fileName=from->preCalcPatts[i][j].fileName;
					to->preCalcPatts[i][j].shellIndex=from->preCalcPatts[i][j].shellIndex;
					//TODO see if it is worthwhile to copy the loaded intensities.
				}
			}
		}
	}
}

void Compare::CalcA()
{
	try{
		if (A){
			for (int i=0;i<nParam;i++){
				delete [] A[i];
			}
			delete [] A;
		}
		A=new double* [nParam];
		for (int i=0;i<nParam;i++){
			A[i]= new double [(i+1)];
		}
	}
	catch (exception& e){
		printf("\nException found in Compare::CalcA: %s\n", e.what());
		exit(1);
	}
	for(int j=0;j<nParam;j++){
		for (int i=j;i<nParam;i++){
			A[i][j]=0;
			for (int k=0;k<Pattern::nInt;k++){
				A[i][j]+=gradI[i][k]*gradI[j][k]/obsI->stdDevSq[k];
			}
		}
	}
}
void Compare::Calcg()
{
	try{
		if (g){
			delete [] g;
		}
		g= new double [nParam];
	}
	catch (exception& e){
		printf("\nException found in Calcg: %s\n", e.what());
		exit(1);
	}

	for (int i=0;i<nParam;i++){
		g[i]=0;
		for (int j=0;j<Pattern::nInt;j++){
			g[i]+=((obsI->I[j]-finalCalcI->I[j])/obsI->stdDevSq[j])*gradI[i][j];
		}
	}
}
int Compare::CholeskyDecomp(double **m, int n)
{
	int i,j,k;
	double temp;
	double **L;
	try{
		L = new double* [n];
		for (i=0;i<n;i++)
			L[i]= new double [(i+1)];
	}
	catch (exception& e){
		printf("\nException found in CholeskyDecomp alloc: %s\n", e.what());
		return -1;
	}

	for(j=0;j<n;j++){
		for(i=j;i<n;i++){
			temp=0;
			if (i==j){
				for (k=0;k<i;k++)
					temp+=L[i][k]*L[i][k];
				temp=m[i][j]-temp;
				if (temp>=0){
					L[i][j]=sqrt(temp);
				}
				else{
					printf("\nError: neg value in Cholesky Decomposition\n");
					return -1;
				}
			}
			else{
				for (k=0;k<j;k++) temp+=L[i][k]*L[j][k];
				if (L[j][j]!=0){
					L[i][j]=(m[i][j]-temp)/L[j][j];
				}
				else{
					printf("\nError: zero value in Cholesky Decomposition\n");
					return -1;
				}
			}
		}
	}
	/*
	Beginning of more general diagonal decomp function
	for(j=0;j<n;j++){
		for(i=j;i<n;i++){
			temp=0;
			for (k=0;k<i;k++)
				temp+=L[i][k]*L[i][k]*D[k];
			D[i]=m[i][i]-temp;
		}
	}
	*/
	for (j=0;j<n;j++){
		for (i=j;i<n;i++){
			m[i][j]=L[i][j];
		}
	}
	for (i=0;i<n;i++)
		delete [] L[i];
	delete [] L;
	return 0;
}

int Compare::CholeskyBackSub(double **L, int n, double *b, double *x)
{
	int i,j;
	double temp;
	double *y;
	try{
		y=new double [n];
	}
	catch (exception& e){
		printf("\nException found in CholeskyBackSub alloc: %s\n", e.what());
		return -1;
	}
	if (x!=0){
		for (i=0;i<n;i++){
			temp=0;
			for (j=0;j<i;j++){
				temp+=L[i][j]*y[j];
			}
			y[i]=(b[i]-temp)/L[i][i];
		}
		for (i=(n-1); i>=0;i--){
			temp=0;
			for (j=(i+1);j<n;j++){
				temp+=L[j][i]*x[j];
			}
			x[i]=(y[i]-temp)/L[i][i];
		}
	}
	else{
		printf("\nError in Cholesky Back Sub memory allocation\n");
		return -1;
	}
	delete [] y;
	return 0;
}
void Compare::AdjustShapeWeights()
{
	size_t pos;
	double temp=0;
	Shape *curShape=firstShape;
	for (int i=0;i<nParam; i++)
	{
		pos=indParams[i]->name.find("_w");
		if (pos!=string::npos) temp+=this->delta[i];
	}
	curShape->wShape->ChangeVal(curShape->wShape->GetVal()-temp);
}
void Compare::NormalizeShapeWeights()
{
	float wTot=0, temp=0;
	Shape *curShape=firstShape;
	while (curShape!=0){
		if (curShape->wShape->GetVal()!=0){
			wTot+=curShape->wShape->GetVal();
			curShape=curShape->nextShape;
		}
	}
	bool renorm=false;
	//Normalize and output change
	curShape=firstShape;
	while(curShape!=0){
		temp=curShape->wShape->GetVal();
		if (curShape->wShape->GetVal()!=0){
			double temp=curShape->wShape->GetVal()/wTot;
			curShape->wShape->ChangeVal(temp);
			if (curShape->wShape->GetVal()!=temp){
				renorm=true;
			}
			cout<<"The weight of "<<curShape->name<<" has been normalized from: "<<temp<<" to: "<<curShape->wShape->GetVal()<<endl;
			curShape=curShape->nextShape;
		}
	}
	if (renorm){
		NormalizeShapeWeights();
	}

}
double Compare::Cosc(double x)
{
	if (x!=0)
		return (cos(x)*x-sin(x))/(x*x);
	else
		return 0;
}
int Compare::SymDiagDecomp(double **m, int n, double *b, double *x)
{
	int err=0;
	double temp;
	double *D;
	double **L;
	try{
		L = new double* [n];
		D = new double [n];
		for (int i=0;i<n;i++)
			L[i]= new double [(i+1)];
	}
	catch (exception& e){
		printf("\nException found in SymDiagDecomp alloc: %s\n", e.what());
		return -1;
	}

	for(int j=0;j<n;j++){
		temp=0;
		L[j][j]=1;
		for (int k=0;k<j;k++)
			temp+=L[j][k]*L[j][k]*D[k];
		D[j]=m[j][j]-temp;
		if (D[j]!=0){
			for(int i=j+1;i<n;i++){
				temp=0;
				for (int k=0;k<j;k++)
					temp+=L[i][k]*L[j][k]*D[k];
				L[i][j]=(m[i][j]-temp)/D[j];
			}
		}
		else {
			printf("\nError in SymDiagDecomp: Cannot divide by zero\n");
			cout<<"index number: "<<j<<endl;
			return -1;
		}
	}

	err=SymDiagBackSub(L,D,n,b,x);

	for (int i=0;i<n;i++)
		delete [] L[i];
	delete [] L;
	delete [] D;
	return err;
}
int Compare::SymDiagBackSub(double **L, double *D, int n, double *b, double *x)
{
	int i,j;
	double temp;
	double *y;
	try{
		y=new double[n];
	}
	catch (exception& e){
		printf("\nException found in SymDiagBackSub alloc: %s\n", e.what());
		return -1;
	}
	for (i=0;i<n;i++){
		temp=0;
		for (j=0;j<i;j++){
			temp+=L[i][j]*D[j]*y[j];
		}
		y[i]=(b[i]-temp)/D[i];
	}
	for (i=(n-1); i>=0;i--){
		temp=0;
		for (j=(i+1);j<n;j++){
			temp+=L[j][i]*x[j];
		}
		x[i]=(y[i]-temp);
	}
	delete [] y;
	return 0;
}
void Compare::CalcResidual()
{
	if (finalCalcI->residual==0){
		try{
			finalCalcI->residual=new double [Pattern::nInt];
		}
		catch (exception &e){
			cout<<"Exception in ComputeResidual: "<<e.what();
			exit(1);
		}
	}
	for (int i=0;i<Pattern::nInt;i++){
		finalCalcI->residual[i]=obsI->I[i]-finalCalcI->I[i];
	}
}
double Compare::CalcNumAtoms()
{
	double nA=0;
	Shape *shape=firstShape;
	while(shape!=0){
		shape->CalcNumAtoms();
		nA+=shape->nAtoms*shape->wShape->GetVal();
		shape=shape->nextShape;
	}
	return nA;
}
void Compare::CleanGlobalParam(Param **p)
{
	try{
		if (*p){
			if ((*p)->global){
				delete (*p);
				*p=0;
			}
		}
	}
	catch(exception &e){
		cout<<"Exception caught in Compare::CleanGlobalParam : "<<e.what();
		exit(1);
	}
}
void Compare::CalcCovarMatrix(double Chisq)
{
	//Use Calculated A matrix and Chisq to calculate covariance matrix
	// sig_ij = A-1_ij*ChiSq
	//See: Scott, H.G.(1983) J. Appl. Cryst. 16, 159-163
	gsl_matrix *A1=gsl_matrix_alloc(nParam, nParam);
	for (int i=0;i<nParam;i++){
		for(int j=0;j<i+1;j++){
			gsl_matrix_set(A1, i,j,A[i][j]);
			if(i!=j){
				gsl_matrix_set(A1,j,i,A[i][j]);
			}
		}
	}
	//LU decomposition of A
	gsl_permutation *p=gsl_permutation_alloc(nParam);
	int sig;
	gsl_linalg_LU_decomp(A1, p, &sig);

	try{
		if (covar){
			gsl_matrix_free(covar);
		}
		covar=gsl_matrix_alloc(nParam,nParam);
	}
	catch (exception &e){
		cout<<"Exception caught in Compare::CalcCovarMatrix: "<<e.what();
		exit(1);
	}
	//Invert A and multiply with chisq
	gsl_linalg_LU_invert(A1, p, covar);
	// Prince (1993) in "The Reitveld Method" by R.A. Young says that this is not always necessary
	// but that it tries to account for an overweighting of all data points
	gsl_matrix_scale(covar, Chisq);
	//Cleanup
	gsl_matrix_free(A1);
	gsl_permutation_free(p);
}
void Compare::SolveForDelta(double **A1, double *g1, double *del, int n)
{
	gsl_matrix *a=0;
	gsl_permutation *p=0;
	gsl_vector *grad=0;
	gsl_vector *d=0;

	try{
		a=gsl_matrix_alloc(n,n);
		p=gsl_permutation_alloc(n);
		grad=gsl_vector_alloc(n);
		d=gsl_vector_alloc(n);
	}
	catch(exception &e){
		cout<<"Exception caught in Compare::SolveForDelta: "<<e.what();
		exit(1);
	}
	for (int i=0;i<n;i++){
		for(int j=0;j<=i;j++){
			gsl_matrix_set(a, i,j, A1[i][j]);
			if(i!=j){
				gsl_matrix_set(a,j,i,A1[i][j]);
			}
		}
		gsl_vector_set(grad, i,g1[i]);
	}

	//LU decomposition of A
	int sig;
	gsl_linalg_LU_decomp(a, p, &sig);
	//Gsl solve using LU decomp
	gsl_linalg_LU_solve(a, p, grad, d);

	for(int i=0;i<n;i++){
		del[i]=gsl_vector_get(d,i);
	}
	//Cleanup
	gsl_matrix_free(a);
	gsl_permutation_free(p);
	gsl_vector_free(grad);
	gsl_vector_free(d);
}
