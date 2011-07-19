/*
	Program created to fit experimental diffraction patterns with calculated Debye patterns
	Ken Beyerlein
	10/4/2008
	Universita degli Studi di Trento  / Georgia Institute of Technology
*/


#include "Compare.h"
bool minimize=false, calcPatt=false, volNorm=true;
int Pattern::nInt=0, err=0, nIter=0, Pattern::nLambda=0;
Shape *firstShape=0;
double Pattern::min2Theta=0, Pattern::max2Theta=0, *Pattern::wLambda=0, *Pattern::lambda=0, Pattern::intStep=0, *Pattern::theta2=0, **Pattern::q=0, *Pattern::offset2Theta=0, **Pattern::offsetq=0;
double Pattern::a1=0.0f, Pattern::b1=0.0f, Pattern::a2=0.0f, Pattern::b2=0.0f, Pattern::a3=0.0f, Pattern::b3=0.0f, Pattern::a4=0.0f, Pattern::b4=0.0f, Pattern::c=1.0f;
Param *Pattern::B=0;
Param **Pattern::ChebBG=0, **Pattern::ChebScale=0;
int Pattern::nChebPoly=0, Pattern::nChebBG=0, Pattern::nChebScale=0;
double **Pattern::chebPolys=0, **Pattern::gaussBG=0;
Param *Pattern::scale=0;
Param *Pattern::sampleDisplace=0;
double Pattern::gonioRad=0;
int Pattern::polQ=0;
Param *Pattern::monoAngle2T=0;
int Pattern::nGaussBG=0;
Param **Pattern::gaussA=0, **Pattern::gaussMu=0, **Pattern::gaussSig=0;
double Pattern::ax=0, Pattern::bx=0, Pattern::cx=0, Pattern::dx=0, Pattern::ex=0;
bool Pattern::synchThetaShift=false;
Param *Pattern::AExp=0, *Pattern::kExp=0, *Pattern::delExp=0;
string Pattern::geometry="";
Param *Pattern::t=0, *Pattern::absCoef=0;
bool Pattern::multCos=false, Pattern::divSin=false;

string obsIntPath, obsIntName, intOutputPath, intOutputName;
/*
void CopyParams(Param * from, Param *to)
{
	to->name=from->name;
	to->shapeName=from->shapeName;
	to->val=from->val;
	to->min=from->min;
	to->max=from->max;
	to->fixed=from->fixed;
	to->global=from->global;
}
*/
void InputPositions()
{
	Shape *curShape=firstShape;
	while ((curShape !=0) && (err==0) ){
		err=curShape->LoadPositions();
		curShape=curShape->nextShape;
	}
}

void InputDistances()
{
	Shape *curShape=firstShape;
	while ((curShape !=0) && (err==0) ){
		curShape->LoadDistances();
		curShape=curShape->nextShape;
	}
}
void InputPreCalcPatterns()
{
	Shape *curShape=firstShape;
	while((curShape!=0)&& (err==0)){
		if (curShape->pattPath!=""){
			curShape->FindPreCalcPatts();
		}
		curShape=curShape->nextShape;

	}
}


void SetDiffPatternParams()
{
	//Initialize unspecified parameters with default values
	if (Pattern::min2Theta==0)
		Pattern::min2Theta=0;
	if (Pattern::max2Theta==0)
		Pattern::max2Theta=175;
	if (Pattern::nInt==0)
		Pattern::nInt=1000;
	if (Pattern::lambda==0){
		Pattern::lambda=new double[1];
		Pattern::lambda[0]=1.54056;
	}
	if (Pattern::intStep==0)
		Pattern::intStep=(Pattern::max2Theta-Pattern::min2Theta)/(Pattern::nInt-1);

	try{
		Pattern::theta2= new double [Pattern::nInt];
		Pattern::q=new double* [Pattern::nLambda];
		for (int i=0;i<Pattern::nLambda;i++)
			Pattern::q[i]=new double [Pattern::nInt];
	}
	catch (exception& e){
		printf("\nException found in SetDiffParams alloc: %s\n", e.what());
		err=-1;
		return;
	}
	//Normalize relative wavelength intensities
	double totLamWeight=0;
	for (int i=0; i< Pattern::nLambda; i++){
		totLamWeight+=Pattern::wLambda[i];
	}
	for (int i=0;i<Pattern::nLambda; i++){
		Pattern::wLambda[i]/=totLamWeight;
	}

	Pattern::offset2Theta=Pattern::theta2;
	Pattern::offsetq=Pattern::q;
	for(int i=0;i<Pattern::nInt;i++){
		Pattern::theta2[i]=Pattern::min2Theta+i*Pattern::intStep;
		for (int j=0; j<Pattern::nLambda;j++){
			Pattern::q[j][i]=4*PI*sin(Pattern::theta2[i]*PI/360.0f)/Pattern::lambda[j];
		}
	}
}
void DeletePattParams()
{
	delete Pattern::AExp;
	delete Pattern::kExp;
	delete Pattern::delExp;
	for (int i=0;i<Pattern::nGaussBG;i++){
		delete [] Pattern::gaussBG[i];
		delete Pattern::gaussA[i];
		delete Pattern::gaussMu[i];
		delete Pattern::gaussSig[i];
	}
	delete [] Pattern::gaussBG;
	delete [] Pattern::gaussA;
	delete [] Pattern::gaussMu;
	delete [] Pattern::gaussSig;
	if (Pattern::chebPolys){
		for (int i=0;i<Pattern::nChebPoly;i++){
			delete [] Pattern::chebPolys[i];
		}
		delete [] Pattern::chebPolys;
	}
	if (Pattern::ChebBG){
		for (int i=0;i<Pattern::nChebBG;i++){
			delete Pattern::ChebBG[i];
		}
		delete [] Pattern::ChebBG;
	}
	if (Pattern::ChebScale){
		for (int i=0;i<Pattern::nChebScale;i++){
			delete Pattern::ChebScale[i];
		}
		delete [] Pattern::ChebScale;
	}
	delete Pattern::B;
	delete Pattern::scale;
	delete Pattern::sampleDisplace;
	delete Pattern::monoAngle2T;
	delete [] Pattern::wLambda;
	delete [] Pattern::lambda;
	if (Pattern::offset2Theta!=Pattern::theta2){
		delete [] Pattern::offset2Theta;
	}
	delete [] Pattern::theta2;
	//delete [] Pattern::offset2Theta;
	if (Pattern::offsetq!=Pattern::q){
		for (int i=0;i<Pattern::nLambda; i++){
			delete [] Pattern::offsetq[i];
		}
		delete [] Pattern::offsetq;
	}
	if (Pattern::q){
		for (int i=0;i<Pattern::nLambda; i++){
			delete [] Pattern::q[i];
		}
		delete [] Pattern::q;
	}
	delete Pattern::t;
	delete Pattern::absCoef;

}
//Counts the number of free parameters used in minimization
int CountIndParams()
{
	int nParams=0;
	Shape *curShape=firstShape;
	if (!Pattern::B->fixed) nParams++;
	if (Pattern::t){
		if (!Pattern::t->fixed) nParams++;
	}

	if (curShape->a!=0){
		for (int i=0;i<curShape->nParams_a;i++){
			if (!curShape->a[i]->fixed) nParams++;
		}
	}
	if (curShape->distribution=="LogNorm"){
		if (!curShape->mu->fixed) nParams++;
		if (!curShape->sigma->fixed) nParams++;
	}
	if (curShape->distribution=="Gamma"){
		if (!curShape->alpha->fixed) nParams++;
		if (!curShape->beta->fixed) nParams++;
	}
	if (curShape->distribution=="Delta"){
		if (!curShape->deltaSize->fixed) nParams++;
	}
	if (!curShape->f->fixed) nParams++;
	if (!curShape->kappa->fixed) nParams++;
	if (!curShape->wShape->fixed) nParams++;
	curShape=curShape->nextShape;
	while (curShape!=0){
		if (!curShape->wShape->global&&!curShape->wShape->fixed) nParams++;
		if (curShape->a!=0){
			for (int i=0;i<curShape->nParams_a;i++){
				if (!curShape->a[i]->global&&!curShape->a[i]->fixed) nParams++;
			}
		}
		if (curShape->distribution=="LogNorm"){
			if (!curShape->mu->global&&!curShape->mu->fixed) nParams++;
			if (!curShape->sigma->global&&!curShape->sigma->fixed) nParams++;
		}
		if (curShape->distribution=="Gamma"){
			if (!curShape->alpha->global&&!curShape->alpha->fixed) nParams++;
			if (!curShape->beta->global&&!curShape->beta->fixed) nParams++;
		}
		if (curShape->distribution=="Delta"){
			if (!curShape->deltaSize->global&&!curShape->deltaSize->fixed) nParams++;
		}
		if (!curShape->f->global&&!curShape->f->fixed) nParams++;
		if (!curShape->kappa->global&&!curShape->kappa->fixed) nParams++;
		curShape=curShape->nextShape;
	}
	if (Pattern::ChebBG!=0){
		for (int i=0; i<Pattern::nChebBG; i++){
			if (!Pattern::ChebBG[i]->fixed) nParams++;
		}
	}
	if (Pattern::ChebScale!=0){
		for (int i=0; i<Pattern::nChebScale; i++){
			if (!Pattern::ChebScale[i]->fixed) nParams++;
		}
	}
	if (Pattern::nGaussBG!=0){
		for (int i=0;i<Pattern::nGaussBG; i++){
			if (!Pattern::gaussA[i]->fixed) nParams++;
			if (!Pattern::gaussMu[i]->fixed) nParams++;
			if (!Pattern::gaussSig[i]->fixed) nParams++;
		}
	}
	if (!Pattern::scale->fixed) nParams++;
	if (Pattern::sampleDisplace!=0){
		if (!Pattern::sampleDisplace->fixed) nParams++;
	}
	if (Pattern::AExp!=0&&Pattern::kExp!=0){
		if(!Pattern::AExp->fixed) nParams++;
		if(!Pattern::kExp->fixed) nParams++;
	}
	if (Pattern::delExp){
		if(!Pattern::delExp->fixed) nParams++;
	}
	return nParams;
}
void ReadInputFile(string inFilePath, string inFile)
{
	err=0;
	bool exitShape=false;
	string temp;
	int global_nParams_a=0;
	Param **global_a=0, *global_mu=0, *global_sigma=0, *global_f=0, *global_kappa=0;
	Shape *curShape=0, *prevShape=0;
	Pattern::scale=new Param;
	Pattern::scale->fixed=false; Pattern::scale->min=0; Pattern::scale->max=1000000; Pattern::scale->ChangeVal(1.0); Pattern::scale->name="Scale";
	//if (!_chdir(inFilePath.c_str())){
	temp=inFilePath+"/"+inFile;
	ifstream input;
	input.open(temp.c_str(), ios::in);
	if(input.is_open()){
		printf("Reading from input file: %s\n", inFile.c_str());
		while(!input.eof() && err==0){
			err=1;
			exitShape=false;
			input >> temp;

	//Input Information about particle shape
			if (temp=="Shape"){
				err=0;
				curShape = new Shape();
				if (firstShape==0) firstShape=curShape;
				input>>curShape->name;
				while(!input.eof() &&!err&&!exitShape){
					err=1;
					input>>temp;
				//Type of shape ie: sphere, cube ...
					if(temp=="type"){
						err=0;
						input>>curShape->type;
						curShape->CheckSupportedShapes();
					}
				//Lattice parameter (constant for all sizes)
					if (temp=="a"){
						err=0;
						curShape->nParams_a=1;
						try {
							curShape->a = new Param* [curShape->nParams_a];
							for (int i=0;i<curShape->nParams_a;i++){
								curShape->a[i]=new Param;
							}
						}
						catch (exception& e){
							printf("\nException found in parameter initialization: %s\n", e.what());
							err=-1;
						}
						if (!err&&curShape->a){
							for (int i=0;i<curShape->nParams_a;i++){
								string fix;
								double min, max, val;
								input>>val>>min>>max>>fix;
								curShape->a[i]->InitParam(val, min, max, temp, curShape->name, fix=="fixed", false);
								curShape->a[i]->CheckParam();
							}
						}
					}
				//Lattice parameter function of particle size
					if (temp=="a(D)"){
						err=0;
						string temp1;
						input>>temp1>>curShape->nParams_a;
						if (temp1=="Poly"){
							try {
								curShape->a = new Param* [curShape->nParams_a];
								for (int i=0;i<curShape->nParams_a;i++){
									curShape->a[i]=new Param;
								}
							}
							catch (exception& e){
								printf("\nException found in parameter initialization: %s\n", e.what());
								err=-1;
							}
							if (!err&&curShape->a){
								for (int i=0;i<curShape->nParams_a;i++){
									ostringstream temp2;
									temp2<<i;
									string fix;
									double min, max, val;
									input>>val>>min>>max>>fix;
									curShape->a[i]->InitParam(val, min, max, "a_"+temp1+temp2.str(), curShape->name, fix=="fixed", false);
									curShape->a[i]->CheckParam();
								}
							}
						}
						else{
							cout<<"Error a(D): "<<temp1<<"not supported.\n";
							err=-1;
						}
					}
				//Parameters for LogNormal Size Distribution
					if (temp=="mu"){
						err=0;
						curShape->mu = new Param;
						string fix;
						double min, max, val;
						input>>val>>min>>max>>fix;
						curShape->mu->InitParam(val, min, max, temp, curShape->name, fix=="fixed", false);
						curShape->mu->CheckParam();
					}
					if (temp=="sigma"){
						err=0;
						curShape->sigma = new Param;
						string fix;
						double min, max, val;
						input>>val>>min>>max>>fix;
						curShape->sigma->InitParam(val, min, max, temp, curShape->name, fix=="fixed", false);
						curShape->sigma->CheckParam();
					}
				//Parameters for Gamma Size Distribution
					if (temp=="alpha"){
						err=0;
						curShape->alpha = new Param;
						string fix;
						double min, max, val;
						input>>val>>min>>max>>fix;
						curShape->alpha->InitParam(val, min, max, temp, curShape->name, fix=="fixed", false);
						curShape->alpha->CheckParam();
					}
					if (temp=="beta"){
						err=0;
						curShape->beta = new Param;
						string fix;
						double min, max, val;
						input>>val>>min>>max>>fix;
						curShape->beta->InitParam(val, min, max, temp, curShape->name, fix=="fixed", false);
						curShape->beta->CheckParam();
					}
				//Parameters controlling exponential surface strain
					if (temp=="f"){
						err=0;
						curShape->f = new Param;
						string fix;
						double min, max, val;
						input>>val>>min>>max>>fix;
						curShape->f->InitParam(val, min, max, temp, curShape->name, fix=="fixed", false);
						curShape->f->CheckParam();
					}
					if (temp=="kappa"){
						err=0;
						curShape->kappa = new Param;
						string fix;
						double min, max, val;
						input>>val>>min>>max>>fix;
						curShape->kappa->InitParam(val, min, max, temp, curShape->name, fix=="fixed", false);
						curShape->kappa->CheckParam();
					}
					if (temp=="n_R"){
						err=0;
						curShape->n_R = new Param;
						string fix;
						double min, max, val;
						input>>val>>min>>max>>fix;
						curShape->n_R->InitParam(val, min, max, temp, curShape->name, fix=="fixed", false);
						curShape->n_R->CheckParam();
					}
				//Shape Fraction (Numerical)
					if (temp=="weight"){
						err=0;
						curShape->wShape = new Param;
						string fix;
						double min, max, val;
						input>>val>>min>>max>>fix;
						curShape->wShape->InitParam(val, min, max, curShape->name+"_w", curShape->name, fix=="fixed", false);
						curShape->wShape->CheckParam();
					}
				//Folder containing position files of each shell
					if(temp=="posPath"){
						err=0;
						input>>curShape->posPath;
					}
				//Folder containing distance files of each shell
					if(temp=="distPath"){
						err=0;
						input>>curShape->distPath;
					}
					if(temp=="pattPath"){
						err=0;
						input>>curShape->pattPath;
					}
				//Type of size distribution to use
					if(temp=="distribution"){
						err=0;
						input>>curShape->distribution;
					}
				//Particle size for Delta size distribution
					if (temp=="deltaSize"){
						err=0;
						curShape->deltaSize = new Param;
						string fix;
						double min, max, val;
						input>>val>>min>>max>>fix;
						curShape->deltaSize->InitParam(val, min, max, curShape->name+"_deltaSize", curShape->name, fix=="fixed", false);
						curShape->deltaSize->CheckParam();
					}
				// Marker to indicate end of input for a given shape
					if (temp==";"){
						err=0;
						exitShape=true;
					}
				// Comment line
					if (temp=="#"){
						err=0;
						input.ignore(256, '\n');
					}
				}
				if(curShape->distPath==""||curShape->posPath==""){
					err=-1;
					printf("Distance or Postion file path not specified for %s\n", curShape->name.c_str());
				}
				if (curShape->name==""){
					err=-1;
					cout<<"Shape name has not been specified.\n";
				}
				if (curShape->type==""){
					err=-1;
					cout<<"Shape type has not been specified.\n";
				}
			// If a global variable has already been declared then replace the variables with global
				if(global_mu){
					if (curShape==firstShape){
						if (!curShape->mu){
							curShape->mu=new Param;
						}
						curShape->mu->CopyParam(global_mu);
					}
					else{
						if (curShape->mu){
							cout<<"\nWarning overwriting parameter "<<curShape->mu->name<< " of shape "<<curShape->name<<" with global value.\n";
							delete curShape->mu;
						}
						curShape->mu=firstShape->mu;
					}
				}
				if(global_sigma){
					if (curShape==firstShape){
						if (!curShape->sigma){
							curShape->sigma=new Param;
						}
						curShape->sigma->CopyParam(global_sigma);
					}
					else{
						if (curShape->sigma){
							cout<<"\nWarning overwriting parameter "<<curShape->sigma->name<< " of shape "<<curShape->name<<" with global value.\n";
							delete curShape->sigma;
						}
						curShape->sigma=firstShape->sigma;
					}
				}
				if(global_f){
					if (curShape==firstShape){
						if (!curShape->f){
							curShape->f=new Param;
						}
						curShape->f->CopyParam(global_f);
					}
					else{
						if (curShape->f){
							cout<<"\nWarning overwriting parameter "<<curShape->f->name<< " of shape "<<curShape->name<<" with global value.\n";
							delete curShape->f;
						}
						curShape->f=firstShape->f;
					}
				}
				if(global_kappa){
					if (curShape==firstShape){
						if (!curShape->kappa){
							curShape->kappa=new Param;
						}
						curShape->kappa->CopyParam(global_kappa);
					}
					else{
						if (curShape->kappa){
							cout<<"\nWarning overwriting parameter "<<curShape->kappa->name<< " of shape "<<curShape->name<<" with global value.\n";
							delete curShape->kappa;
						}
						curShape->kappa=firstShape->kappa;
					}
				}
				if (global_a){
					try{
						if (curShape->a){
							cout<<"\nWarning overwriting parameter "<<curShape->a[0]->name<< " of shape "<<curShape->name<<" with global value.\n";
							for (int i=0;i<curShape->nParams_a;i++){
								delete curShape->a[i];
							}
							delete [] curShape->a;
						}
						curShape->nParams_a=global_nParams_a;
						if (curShape==firstShape){
							curShape->a=new Param*[curShape->nParams_a];
							for (int i=0;i<curShape->nParams_a;i++){
								curShape->a[i]= new Param;
								curShape->a[i]->CopyParam(global_a[i]);
							}
						}
						else{
							curShape->a=firstShape->a;
						}
					}
					catch(exception &e){
						cout<<"Exception in lattice parameter alloc: "<<e.what()<<endl;
						exit(1);
					}
				}

				if(curShape->a==0||curShape->wShape==0) {
					err=-1;
					printf("Must specify lattice parameter and weight for %s\n", curShape->name.c_str());
				}
				if(prevShape!=0) prevShape->nextShape=curShape;
				prevShape=curShape;
			}

		// Input of global variables (see above for function)
			if (temp=="a"){
				err=0;
				global_nParams_a=1;
				try {
					global_a = new Param* [global_nParams_a];
					for (int i=0;i<global_nParams_a;i++){
						global_a[i]=new Param;
					}
				}
				catch (exception& e){
					printf("\nException found in parameter initialization: %s\n", e.what());
					err=-1;
				}
				if (!err&&global_a){
					for (int i=0;i<global_nParams_a;i++){
						string fix;
						double min, max, val;
						input>>val>>min>>max>>fix;
						global_a[i]->InitParam(val, min, max, "a", "global", fix=="fixed", true);
						global_a[i]->CheckParam();
						//global_a[i]->name="a"; global_a[i]->shapeName="global"; double val; input>>val; input>>global_a[i]->min; input>>global_a[i]->max; global_a[i]->ChangeVal(val);
						//input>>temp; if(temp=="fixed") global_a[i]->fixed=true; else global_a[i]->fixed=false;
						//global_a[i]->global=true;
					}
				}
			}
			if (temp=="a(D)"){
				err=0;
				string temp1;
				input>>temp1>>global_nParams_a;
				if (temp1=="Poly"){
					try {
						global_a = new Param* [global_nParams_a];
						for (int i=0;i<global_nParams_a;i++){
							global_a[i]=new Param;
						}
					}
					catch (exception& e){
						printf("\nException found in parameter initialization: %s\n", e.what());
						err=-1;
					}
					if (!err&&global_a){
						for (int i=0;i<global_nParams_a;i++){
							ostringstream temp2;
							temp2<<i;

							string fix;
							double min, max, val;
							input>>val>>min>>max>>fix;
							global_a[i]->InitParam(val, min, max, "a"+temp1+temp2.str(), "global", fix=="fixed", true);
							global_a[i]->CheckParam();
							/*
							global_a[i]->name="a_"+temp1+temp2.str();
							global_a[i]->shapeName="global";
							input>>global_a[i]->val>>global_a[i]->min>>global_a[i]->max;
							input>>temp; if(temp=="fixed") global_a[i]->fixed=true; else global_a[i]->fixed=false;
							global_a[i]->global=true;
							*/
						}
					}
				}
				else{
					err=-1;
					cout<<"Error a(D): " <<temp1<<" not supported."<<endl;
				}
			}
			if (temp=="mu"){
				err=0;
				try {
					global_mu = new Param;
				}
				catch (exception& e){
					printf("\nException found in parameter initialization: %s\n", e.what());
					err=-1;
				}
				if (!err&&global_mu){
					string fix;
					double min, max, val;
					input>>val>>min>>max>>fix;
					global_mu->InitParam(val, min, max, temp, "global", fix=="fixed", true);
					global_mu->CheckParam();
					/*
					global_mu->name="mu"; global_mu->shapeName="global"; input>>global_mu->val; input>>global_mu->min; input>>global_mu->max;
					input>>temp; if(temp=="fixed") global_mu->fixed=true; else global_mu->fixed=false;
					global_mu->global=true;
					*/
				}
			}
			if (temp=="sigma"){
				err=0;
				try {
					global_sigma = new Param;
				}
				catch (exception& e){
					printf("\nException found in parameter initialization: %s\n", e.what());
					err=-1;
				}
				if (!err&&global_sigma){
					string fix;
					double min, max, val;
					input>>val>>min>>max>>fix;
					global_sigma->InitParam(val, min, max, temp, "global", fix=="fixed", true);
					global_sigma->CheckParam();
					/*
					global_sigma->name="sigma"; global_sigma->shapeName="global"; input>>global_sigma->val; input>>global_sigma->min; input>>global_sigma->max;
					input>>temp; if(temp=="fixed") global_sigma->fixed=true; else global_sigma->fixed=false;
					global_sigma->global=true;
					*/
				}
			}
			/*
			 // TODO Create support for Global Alpha or Beta
			if (temp=="alpha"){
				err=0;
				try {
					global_sigma = new Param;
				}
				catch (exception& e){
					printf("\nException found in parameter initialization: %s\n", e.what());
					err=-1;
				}
				if (!err&&global_sigma){
					global_sigma->name="sigma"; global_sigma->shapeName="global"; input>>global_sigma->val; input>>global_sigma->min; input>>global_sigma->max;
					input>>temp; if(temp=="fixed") global_sigma->fixed=true; else global_sigma->fixed=false;
					global_sigma->global=true;
				}
			}
			if (temp=="beta"){
				err=0;
				try {
					global_sigma = new Param;
				}
				catch (exception& e){
					printf("\nException found in parameter initialization: %s\n", e.what());
					err=-1;
				}
				if (!err&&global_sigma){
					global_sigma->name="sigma"; global_sigma->shapeName="global"; input>>global_sigma->val; input>>global_sigma->min; input>>global_sigma->max;
					input>>temp; if(temp=="fixed") global_sigma->fixed=true; else global_sigma->fixed=false;
					global_sigma->global=true;
				}
			}
			*/
			//TODO input of Delta size as a global parameter
			if (temp=="f"){
				err=0;
				try {
					global_f = new Param;
				}
				catch (exception& e){
					printf("\nException found in parameter initialization: %s\n", e.what());
					err=-1;
				}
				if (!err&&global_f){
					string fix;
					double min, max, val;
					input>>val>>min>>max>>fix;
					global_f->InitParam(val, min, max, temp, "global", fix=="fixed", true);
					global_f->CheckParam();
					/*
					global_f->name="f"; global_f->shapeName="global";
					input>>global_f->val;
					input>>global_f->min;
					input>>global_f->max;
					input>>temp; if(temp=="fixed") global_f->fixed=true; else global_f->fixed=false;
					global_f->global=true;
					*/
				}
			}
			if (temp=="kappa"){
				err=0;
				try {
					global_kappa = new Param;
				}
				catch (exception& e){
					printf("\nException found in parameter initialization: %s\n", e.what());
					err=-1;
				}
				if (!err&&global_kappa){
					string fix;
					double min, max, val;
					input>>val>>min>>max>>fix;
					global_kappa->InitParam(val, min, max, temp, "global", fix=="fixed", true);
					global_kappa->CheckParam();
					/*
					global_kappa->name="kappa"; global_kappa->shapeName="global";
					input>>global_kappa->val;
					input>>global_kappa->min;
					input>>global_kappa->max;
					input>>temp; if(temp=="fixed") global_kappa->fixed=true; else global_kappa->fixed=false;
					global_kappa->global=true;
					*/
				}
			}

		//Debye-Waller factor
			if (temp=="B"){
				err=0;
				try {
					Pattern::B=new Param;
				}
				catch (exception& e){
					printf("\nException found in parameter initialization: %s\n", e.what());
					err=-1;
				}
				if (!err&&Pattern::B){
					string fix;
					double min, max, val;
					input>>val>>min>>max>>fix;
					Pattern::B->InitParam(val, min, max, temp, "global", fix=="fixed", true);
					Pattern::B->CheckParam();
					/*Pattern::B->name="B"; input>>Pattern::B->val; input>>Pattern::B->min; input>>Pattern::B->max;
					input>>temp; if(temp=="fixed") Pattern::B->fixed=true; else Pattern::B->fixed=false;
					Pattern::B->global=true;
					*/
				}
			}
		//Atomic Scattering factor described by 4 gaussians and a constant
			if (temp=="atomicScatFact"){
				err=0;
				input>>Pattern::a1>>Pattern::b1>>Pattern::a2>>Pattern::b2>>Pattern::a3>>Pattern::b3>>Pattern::a4>>Pattern::b4>>Pattern::c;
			}
		//Input the tyepe of calculation you intend to carry out
			if (temp=="Function"){
				err=0;
				input>>temp;
				if (temp=="minimize")
					minimize=true;
				else if (temp=="calcPattern")
					calcPatt=true;
				else
					printf("Function must be 'minimize' or 'calcPatt'\n");
			}
		//Folder containing observed intensity
			if (temp=="obsIntPath"){
				err=0;
				input>>obsIntPath;
			}
		//Name of observed intensity file
			if (temp=="obsIntName"){
				err=0;
				input>>obsIntName;
			}
		//Folder to output intensity files
			if (temp=="intOutputPath"){
				err=0;
				input>>intOutputPath;
			}
		//Name to give output intensity files
			if (temp=="intOutputName"){
				err=0;
				input>>intOutputName;
			}
		//Number of iterations to perform
			if (temp=="nIter"){
				err=0;
				input>>nIter;
			}
		//Addition of Chebyshev Background
			if (temp=="Chebyshev_BG"){
				err=0;
				string fix;
				input >>Pattern::nChebBG>>fix;
				Pattern::nChebPoly=getmax(Pattern::nChebPoly, Pattern::nChebBG);
				try{
					Pattern::ChebBG=new Param* [Pattern::nChebBG];
					for (int i=0;i<Pattern::nChebBG; i++){
						Pattern::ChebBG[i]=new Param;
					}
				}
				catch (exception& e){
					printf("\nException found in parameter initialization: %s\n", e.what());
					err=-1;
				}

				for (int i=0; i<Pattern::nChebBG&&err==0; i++){
					double val;
					input>>val;
					ostringstream temp2;
					temp2<<i;
					Pattern::ChebBG[i]->InitParam(val, -1000000, 1000000, "Cheb_BG"+temp2.str(), "", fix=="fixed", true);
					Pattern::ChebBG[i]->CheckParam();
				}
			}
		//Multiplication of pattern by scale factor which is a function of the diffraction angle
			if (temp=="Chebyshev_Scale"){
				err=0;
				string fix;
				input >>Pattern::nChebScale>>fix;
				Pattern::nChebPoly=getmax(Pattern::nChebPoly, Pattern::nChebScale+1);
				try{
					Pattern::ChebScale=new Param* [Pattern::nChebScale];
					for (int i=0;i<Pattern::nChebScale; i++){
						Pattern::ChebScale[i]=new Param;
					}
				}
				catch (exception& e){
					printf("\nException found in parameter initialization: %s\n", e.what());
					err=-1;
				}
				for (int i=0; i<Pattern::nChebScale&&err==0; i++){
					double val;
					input>>val;
					ostringstream temp2;
					temp2<<i;
					Pattern::ChebScale[i]->InitParam(val, -1000000, 1000000, "Cheb_Scale"+temp2.str(), "", fix=="fixed", true);
					Pattern::ChebScale[i]->CheckParam();
				}
			}
		//Multiplication of pattern by constant scale factor
			if (temp=="scale"){
				err=0;
				string fix;
				double val;
				input>>val>>fix;
				Pattern::scale->InitParam(val, 0.00000000000000001, 1000000, temp, "", fix=="fixed", true);
				Pattern::scale->CheckParam();
			}
		//Correction to account for sample displacement from perfect focusing
			if (temp=="sampleDisplace"){
				err=0;
				Pattern::sampleDisplace=new Param;
				string fix;
				double val;
				input>> val>> Pattern::gonioRad>>fix;
				Pattern::sampleDisplace->InitParam(val,-Pattern::gonioRad, Pattern::gonioRad, "SampleDisplacement", "", fix=="fixed", true);
				Pattern::sampleDisplace->CheckParam();
				/*
				Pattern::sampleDisplace->name="Sample Displacement";
				input>>Pattern::sampleDisplace->val>>Pattern::gonioRad>>temp;
				if (temp=="fixed") Pattern::sampleDisplace->fixed=true; else Pattern::sampleDisplace->fixed=false;
				Pattern::sampleDisplace->min=-Pattern::gonioRad; Pattern::sampleDisplace->max=Pattern::gonioRad;
				*/
			}
		//Range of angles to consider in 2 Theta when calculating pattern
			if (temp=="2ThetaRange"){
				err=0;
				input>>Pattern::min2Theta>>Pattern::max2Theta;
			}
		//Number of intensity points to calculate
			if (temp=="nInt"){
				err=0;
				input>>Pattern::nInt;
			}
		//Wavelengths for which to calculate pattern along with relative ratio
			if (temp=="lambda"){
				err=0;
				input>>Pattern::nLambda;
				try{
					Pattern::lambda= new double [Pattern::nLambda];
					Pattern::wLambda= new double [Pattern::nLambda];
				}
				catch (exception& e){
					printf("\nException found in lambda initialization: %s\n", e.what());
					err=-1;
				}
				if (!err){
					for (int i=0;i<Pattern::nLambda;i++){
						input>>Pattern::lambda[i]>>Pattern::wLambda[i];
					}
				}
			}
		//Multplication of pattern by polarization effect
			if (temp=="monoPolarization"){
				err=0;
				double val;
				string fix;
				input>>Pattern::polQ>>val>>fix;
				Pattern::monoAngle2T=new Param;
				Pattern::monoAngle2T->InitParam(val,0, 180, "MonoAngle", "", fix=="fixed", true);
				Pattern::monoAngle2T->CheckParam();
				/*
				Pattern::monoAngle2T->name="Monochromator Angle 2T";
				input>>Pattern::monoAngle2T->val>>temp;
				if (temp=="fixed") Pattern::monoAngle2T->fixed=true; else Pattern::monoAngle2T->fixed=false;
				Pattern::monoAngle2T->min=0; Pattern::monoAngle2T->max=180;
				*/
			}
		//Addition of Gaussian peak in pattern
			if (temp=="Gaussian_BG"){
				err=0;
				input>>Pattern::nGaussBG;
				try{
					Pattern::gaussA=new Param* [Pattern::nGaussBG];
					Pattern::gaussMu=new Param* [Pattern::nGaussBG];
					Pattern::gaussSig=new Param* [Pattern::nGaussBG];
					for (int i=0;i<Pattern::nGaussBG; i++){
						Pattern::gaussA[i]=new Param;
						Pattern::gaussMu[i]=new Param;
						Pattern::gaussSig[i]=new Param;
					}

				}
				catch (exception& e){
					printf("\nException found in Background Gaussian initialization: %s\n", e.what());
					err=-1;
				}
				if (!err){
					for (int i=0;i<Pattern::nGaussBG; i++){
						double val, min, max;
						string fix;
						ostringstream temp2;
						temp2<<i;
						input>>val>> min>>max>>fix;
						Pattern::gaussA[i]->InitParam(val,min, max, "BG_Gauss_A"+temp2.str(), "", fix=="fixed", true);
						Pattern::gaussA[i]->CheckParam();

						input>>val>> min>>max>>fix;
						Pattern::gaussMu[i]->InitParam(val,min, max, "BG_Gauss_Mu"+temp2.str(), "", fix=="fixed", true);
						Pattern::gaussMu[i]->CheckParam();

						input>>val>> min>>max>>fix;
						Pattern::gaussSig[i]->InitParam(val,min, max, "BG_Gauss_Sigma"+temp2.str(), "", fix=="fixed", true);
						Pattern::gaussSig[i]->CheckParam();
					}
				}
			}
		// ??
			if (temp=="normalize"){
				err=0;
				input>>temp;
				if (temp=="false"){
					volNorm=false;
				}
			}
		//Apply tangent polynomial correction to peak shift
			//TODO consider refining tangent poly theta shift
			if (temp=="synchThetaShift"){
				err=0;
				input>>Pattern::ax>>Pattern::bx>>Pattern::cx>>Pattern::dx>>Pattern::ex;
			}
		//Addition of background function which is an exponential function
			if (temp=="ExpDecayBG"){
				err=0;
				try{
					Pattern::AExp=new Param;
					Pattern::kExp=new Param;
				}
				catch (exception &e){
					printf("\nException found in Exponential Background initialization: %s\n", e.what());
					err=-1;
				}
				for (int i=0;i<2&&!err;i++){
					input>>temp;
					if (temp=="A"){
						double val, min, max;
						string fix;
						input>>val>> min>>max>>fix;
						Pattern::AExp->InitParam(val,min, max, "BG_ExpDecay_A", "", fix=="fixed", true);
						Pattern::AExp->CheckParam();
					}
					else if (temp=="k"){
						double val, min, max;
						string fix;
						input>>val>> min>>max>>fix;
						Pattern::kExp->InitParam(val,min, max, "BG_ExpDecay_k", "", fix=="fixed", true);
						Pattern::kExp->CheckParam();
					}
					else {
						cout<<"Error in Background Exp Decay- Unrecognized parameter name: "<<temp<<endl;
						err=-1;
					}
				}
			}
		//Addition of Shifted exponential background function
			if (temp=="ShiftedExpDecayBG"){
				err=0;
				try{
					Pattern::AExp=new Param;
					Pattern::kExp=new Param;
					Pattern::delExp=new Param;
				}
				catch (exception &e){
					printf("\nException found in Exponential Background initialization: %s\n", e.what());
					err=-1;
				}
				for (int i=0;i<3&&!err;i++){
					input>>temp;
					if (temp=="A"){
						double val, min, max;
						string fix;
						input>>val>> min>>max>>fix;
						Pattern::AExp->InitParam(val,min, max, "BG_ExpDecay_A", "", fix=="fixed", true);
						Pattern::AExp->CheckParam();
					}
					else if (temp=="k"){
						double val, min, max;
						string fix;
						input>>val>> min>>max>>fix;
						Pattern::kExp->InitParam(val,min, max, "BG_ExpDecay_k", "", fix=="fixed", true);
						Pattern::kExp->CheckParam();
					}
					else if (temp=="del"){
						double val, min, max;
						string fix;
						input>>val>> min>>max>>fix;
						Pattern::delExp->InitParam(val,min, max, "BG_ExpDecay_del", "", fix=="fixed", true);
						Pattern::delExp->CheckParam();
					}
					else {
						cout<<"Error in Shifted Background Exp Decay- Unrecognized parameter name: "<<temp<<endl;
						err=-1;
					}
				}
			}
		//Type of geometry used to collect the diffraction pattern
			if (temp=="Geometry"){
				err=0;
				input>>Pattern::geometry;
				Pattern temp;
				temp.CheckGeometry();

			}
		//Multiplication of pattern by absorption effect
			if (temp=="Absorption"){
				err=0;
				try{
					Pattern::t=new Param;
					Pattern::absCoef =new Param;
				}
				catch (exception &e){
					printf("\nException found in Absorption Parameter initialization: %s\n", e.what());
					err=-1;
				}
				for (int i=0;i<2&&!err;i++){
					input>>temp;
					if (temp=="t"){
						double val, min, max;
						string fix;
						input>>val>> min>>max>>fix;
						Pattern::t->InitParam(val,min, max, "Abs_t", "", fix=="fixed", true);
						Pattern::t->CheckParam();
					}
					else if (temp=="mu_l"){
						double val;
						input>>val;
						Pattern::absCoef->InitParam(val,0, 100, "Abs_mu_l", "", true, true);
					}
					else {
						cout<<"Error in Absorption - Unrecognized parameter name: "<<temp<<endl;
						err=-1;
					}
				}
			}
			if (temp=="multCos"){
				err=0;
				input>>temp;
				if (temp=="true"){
					Pattern::multCos=true;
				}
			}
			if (temp=="#"){
				err=0;
				input.ignore(256, '\n');
			}
		}
		input.close();
	}
	else{
		printf("Cannot find input file %s\n", inFile.c_str());
		err=-1;
	}
	if (err==1) printf("Unrecognized input keyword %s\n", temp.c_str());
	//Cleanup global
	if (global_a){
		for (int i=0;i<global_nParams_a;i++){
			delete global_a[i];
		}
		delete [] global_a;
	}
	delete global_mu;
	delete global_sigma;
	delete global_f;
	delete global_kappa;
}

int main(int argc, char**argv)
{
	int i;
	string infile="", path="";

	for (i = 0; i < argc; i++)
	{
		if (string(argv[i])=="-f")
		{
			infile = argv[i+1];
			i++;
		}
		else if (string(argv[i])=="-p")
		{
			path=argv[i+1];
			i++;
			cout<<"Path to input file is set to: "<<path<<endl;
		}
	}
	if (path==""){
		//path="C:/Users/Ken/Documents/Input/Fitting_With_Debye";
		//path="/harbor/kbeyer/DebyeFunctionAnalysis_Input";
		path=".";
	}

	if (infile == "")
	{
		cout<<"No input file specified ... Exiting\n";
		return 0;
	}

//Reading from input file
	else{
		ReadInputFile(path, infile);
		SetDiffPatternParams();
	}
//Reading Position and Distance Files
	if (!err) InputPositions();
	if (!err) InputDistances();
	if (!err) InputPreCalcPatterns();
//Fit pattern
	if (!err && minimize){
		Compare minimizer(intOutputPath, intOutputName, firstShape, CountIndParams());
		minimizer.obsI=new Pattern(obsIntPath+obsIntName);
		err=minimizer.obsI->ReadIntensityFile();
		if (!err) err=minimizer.InitParams();
		if (!err) err=minimizer.LevMarq(nIter);
	}
//Only calculate pattern for desired parameters
	else if (!err && calcPatt){
		Compare pattern (firstShape, 0);
		err=pattern.CalcPattern(false, volNorm);
		if (!err) err=pattern.OutputPattern(intOutputPath, intOutputName,0,false);
	}
	DeletePattParams();

	return 0;

}

