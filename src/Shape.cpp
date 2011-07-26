#include "Shape.h"

//#define DEBUG

Shape::Shape()
{
	name="";
	type="";
	firstDist=0;
	surfRelaxDist=0;
	firstPos=0;
	nextShape=0;
	maxShell=0;
	diffPatt=0;
	/*
	try{
		diffPatt=new Pattern* [Pattern::nLambda];
		for (int i=0;i<Pattern::nLambda;i++){
			diffPatt[i]=new Pattern(shape);
		}
	}
	catch (exception& e){
		printf("\nException found in Shape diffPatt initialization: %s\n", e.what());
	}*/
	element="";
	a=0;
	nParams_a=0;
	a_D=0;
	wShape=0;
	mu=0;
	sigma=0;
	alpha=0;
	beta=0;
	f=0;
	kappa=0;
	n_R=0;
	wShell=0;
	nRelaxShells=0;
	derivSRD_mu=0;
	derivSRD_sig=0;
	distribution="";
	deltaSize=0;
	p=0;
	dW_mu=0;
	dW_sig=0;
	endSize=0;
	posPath="";
	distPath="";
	D=0;
	nA_D=0;
	nAtoms=0;
	preCalcPatts=0;
	nPreCalcPatts=0;
}
Shape::~Shape(void)
{
	delete surfRelaxDist;
	if (diffPatt){
		for (int i=0;i<Pattern::nLambda; i++){
			delete diffPatt[i];
		}
		delete [] diffPatt;
	}
	delete [] wShell;
	delete [] p;
	delete [] dW_mu;
	delete [] dW_sig;
	if (a){
		if (!(a[0]->global)){
			for (int i=0;i<nParams_a; i++){
				delete a[i];
			}
			delete [] a;
		}
	}
	delete wShape;
	//NOTE: when adding a global parameter deallocate it in ~Compare()
	if (mu){
		if (!mu->global){
			delete mu;
		}
	}
	if (sigma){
		if (!sigma->global){
			delete sigma;
		}
	}
	if (alpha){
		if (!alpha->global){
			delete alpha;
		}
	}
	if (beta){
		if (!beta->global){
			delete beta;
		}
	}
	if (f){
		if (!f->global){
			delete f;
		}
	}
	if (kappa){
		if (!kappa->global){
			delete kappa;
		}
	}
	if (n_R){
		if (!n_R->global){
			delete n_R;
		}
	}
	if (deltaSize){
		if (!deltaSize->global){
			delete deltaSize;
		}
	}

	delete [] D;
	delete [] a_D;
	delete [] nA_D;
	if (preCalcPatts){
		for (int i=0;i<nPreCalcPatts;i++){
			if (preCalcPatts[i]){
				//for (int j=0;j<Pattern::nLambda;j++){
				//	delete preCalcPatts[i][j];
				//}
				delete [] preCalcPatts[i];
			}
		}
		delete [] preCalcPatts;
	}
}
int Shape::LoadPositions()
{
	//TODO Change so that function kills program in the case of errors
	int err=0;
	Position *curPos, *prevPos=0;
	ifstream posfiles;
	string temp=posPath+"positionFiles.txt";
	posfiles.open(temp.c_str());
	if(posfiles.is_open()){
		printf("Reading position files for shape %s ...\n", name.c_str());
		int count=1;
		while(!posfiles.eof()){
			posfiles >> temp;
			if (!posfiles.eof()){
				curPos = new Position(const_cast<char *>(temp.c_str()));
				if (firstPos==0){
					firstPos=curPos;
					prevPos=curPos;
				}
				else{
					if (prevPos){
						prevPos->nextPos=curPos;
						prevPos=curPos;
					}
				}
				err=curPos->ReadPosFile();
				curPos->shell=count;
				maxShell=getmax(maxShell,curPos->shell);
				if (err!=0) {
					printf("Could not open position file: %s", temp.c_str());
					return -1;
				}
				count++;
			}
		}
		posfiles.close();
	}
	else{
		printf("Error: Cannot find list of position files %s\n", temp.c_str());
		return -1;
	}
	printf("MaxShell: %i\n",maxShell);
	//Check to max sure pos files of all shells are found and sort linked list
	/*
	//TODO Below is pointless because pos->shell is given counter above in loop
	int shell, place, missing=0;
	Position *orderPos;
	orderPos=firstPos;
	for (shell=1;shell<=maxShell;shell++){
		curPos=orderPos;
		prevPos=orderPos;
		if (orderPos==firstPos)
			place=1+missing;
		else
			place=orderPos->shell+missing;
		flag=false;
		while(curPos!=0&&!flag)
		{
			if (curPos->shell==shell){
				flag=true;
				if(place!=shell){
					prevPos->nextPos=curPos->nextPos;
					if (shell==1){
						curPos->nextPos=firstPos;
						firstPos=curPos;
					}
					else{
						curPos->nextPos=orderPos->nextPos;
						orderPos->nextPos=curPos;
					}
				}
				orderPos=curPos;
			}
			prevPos=curPos;
			curPos=curPos->nextPos;
			place++;
		}
		if (!flag){
			printf("Position file for shell %i is missing", shell);
			err=-1;
			missing++;
		}
	}
	*/
	return err;
}

void Shape::LoadDistances()
{
	int nDist=0;
	Distance *curDist, *prevDist=0;
	ifstream distfiles;
	string temp=distPath+"distanceFiles.txt";
	distfiles.open(temp.c_str());
	if(distfiles.is_open()){
		printf("Reading distance files for shape %s ...\n", name.c_str());
		while(!distfiles.eof()){
			distfiles >> temp;
			if (!distfiles.eof()){
				curDist = new Distance(const_cast<char *>(temp.c_str()));
				if (firstDist==0){
					firstDist=curDist;
				}
				else{
					if (prevDist){
						prevDist->nextDist=curDist;
					}
				}
				prevDist=curDist;
				curDist->ReadDistFile();
				nDist++;
			}
		}
		distfiles.close();
	}
	else{
		printf("Could not find list of distance files: %s\n", temp.c_str());
		exit(1);
	}
	curDist=firstDist;
	int maxDistShell=0;
	for (int i=0;i<nDist;i++){
		maxDistShell=getmax(curDist->shell,maxDistShell);
		curDist=curDist->nextDist;
	}
	if (maxDistShell>maxShell){
		cout<<"Error: maximum distance shell is larger than maximum position shell for shape: "<< name<<endl;
	}
	else{
		maxDistShell=getmax(maxDistShell,maxShell);
	}
	//Init ordered list
	Distance **orderedDistList=new Distance * [maxDistShell];
	for (int i=0; i<maxDistShell; i++){
		orderedDistList[i]=0;
	}
	//Fill in already existing distances
	curDist=firstDist;
	for (int i=0; i<nDist; i++){
		orderedDistList[curDist->shell-1]=curDist;
		curDist=curDist->nextDist;
	}
	//Find missing distances, calculate and output them
	string prefix="distances."+type+"Shell", suffix=".debye";
	for (int i=0;i<maxDistShell;i++){
		if (orderedDistList[i]==0){
			int shell=i+1;
			printf("Distance file for shell %i is missing\n", shell);
			ostringstream num;
			num<<shell;
			string filename=prefix+num.str()+suffix;
			orderedDistList[i]= new Distance();
			orderedDistList[i]->fileName=filename;
			orderedDistList[i]->shell=shell;
			CalcDistForShell(orderedDistList[i]);
#define OUTDIST
#ifdef OUTDIST
			orderedDistList[i]->OutputDistToFile(distPath, type);
			printf("Distance file calculated for shell %i\n",shell);
#endif
		}
	}
	//Rebuild Linked List
	firstDist=orderedDistList[0];
	curDist=firstDist;
	for (int i=1;i<maxDistShell;i++){
		curDist->nextDist=orderedDistList[i];
		curDist=curDist->nextDist;
	}
	curDist->nextDist=0;
	delete [] orderedDistList;
}

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
/*
int Shape::InitParamBounds()
{
	int err=0;
	if (wShape->val==-1) err=-1;
	wShape->name=name+"Weight";
	wShape->min=0;
	wShape->max=1;
	wShape->global=false;
	if (a[0]->val==-1) err=-1;
	if (a[0]->min==-1) a[0]->min=0;
	if (a[0]->max==-1) a[0]->max=10;
	// Want to rewrite this to calc min and max given the range of shell sizes
	if (distribution=="LogNorm"){
		if (mu->val==-1) err=-1;
		if (mu->min==-1) mu->min=-1000;
		if (mu->max==-1) mu->max=1000;
		if (sigma->val==-1) err=-1;
		if (sigma->min==-1) sigma->min=0;
		if (sigma->max==-1) sigma->max=1000;
	}
	// Not sure that kappa has no real bounds, need to look into it some more.
	if (f->val==-1) err=-1;
	if (f->min==-1) f->min=0;
	if (f->max==-1) f->max=1;
	if (kappa->val==-1) err=-1;
	if (kappa->min==-1) kappa->min=0;
	if (kappa->max==-1) kappa->max=1000;

	if (err==-1)
		printf("Parameter not initalized for shape %s\n",name.c_str());
	return err;
}
*/
/*
int Shape::SurfaceRelaxation(double *derivfSize_mu, double *derivfSize_sig)
{

		Currently not used and needs some work before use!!

	int i,j, k,err=0;
	printf("Calculating surface relaxation for %s\n", shape.c_str());
	//Initialize temp surface relax position files
	Position **surfRelaxPos=0;
	Distance *tempSRDist=0;
	surfRelaxPos=(Position**) malloc(nRelaxShells*sizeof(Position*));
	if (surfRelaxPos==0) err=-1;
	//for (i=0;i<nRelaxShells&&!err;i++){
		surfRelaxPos[i] = new Position();
		if (surfRelaxPos[i]!=0){
			surfRelaxPos[i]->f=f->val;
			surfRelaxPos[i]->kappa=kappa->val;
		}
		else err=-1;
	}
	//Initlaize surface relax distance files
	if (!err){
		tempSRDist = new Distance ();
		if (tempSRDist==0) err=-1;
	}
	if (!err){
		if (surfRelaxDist==0) surfRelaxDist= new Distance();
		if (surfRelaxDist==0) err=-1;
	}
	if (!err) InitSRDistance(*tempSRDist);
	if (!err) InitSRDistance(*surfRelaxDist);

	if (derivSRD_mu==0&&!err) derivSRD_mu=new Distance();
	if (derivSRD_sig==0&&!err) derivSRD_sig=new Distance();
	if (derivSRD_mu==0||derivSRD_sig==0) err=-1;
	if (!err) err=InitSRDistance(*derivSRD_mu);
	if (!err) err=InitSRDistance(*derivSRD_sig);
	int startSize=1;
	if (!err){
		if (distribution=="LogNorm"){
			//
			if (startSize<1||endSize>maxShell){
				printf("Needed range of Prob function is not possible with given distances\n");
				err=-1;
			}
		}
		else if (distribution=="Delta"){
			startSize=size;
			endSize=startSize;
		}
	}
	//calculate SR and distances for this shape
	Position **posPlaceHold=(Position **) malloc(nRelaxShells*sizeof(Position*));
	for (i=0;i<nRelaxShells&&!err;i++){
		if (posPlaceHold)
			posPlaceHold[i]=firstPos;
		else {printf("Error in memory allocation in Surface Relaxation function\n"); err=-1;}
		surfRelaxPos[i]= new Position();
	}
	for (i=startSize; i<=endSize&&!err; i++){

		if (distribution=="LogNorm"){
			for (j=0;j<nRelaxShells && j<i && !err;j++){
				//if (surfRelaxPos[j]->nAtoms==0) err=CopyPosition(surfRelaxPos[j], firstPos);
				//else err=CopyPosition(surfRelaxPos[j], surfRelaxPos[j]->nextPos);
				if (surfRelaxPos[j]!=0){
					err=CopyPosition(*surfRelaxPos[j],posPlaceHold[j]);
				}
				else err=-1;
				if (!err) surfRelaxPos[j]->ApplyPlanarSurfaceRelax(i);
				else{
					printf("Error in Surface Relaxation position memory allocation\n");
					err=-1;
				}
			}
		}
		else if (distribution=="Delta"){
			Position *tempPos=firstPos;
			for (j=0;j<size-nRelaxShells; j++) tempPos=tempPos->nextPos;
			for (j=nRelaxShells-1; j>=0&&!err; j--){
				err=CopyPosition(*surfRelaxPos[j], tempPos);
				if (!err){
					surfRelaxPos[j]->ApplyPlanarSurfaceRelax(i);
					tempPos=tempPos->nextPos;
				}
				else{
					printf("Error in Surface Relaxation position memory allocation\n");
					err=-1;
				}
			}
		}
		if (!err){
			Position *curPos;
			for (j=0;j<nRelaxShells && j<i;j++){
				curPos=firstPos;
				// for each relaxed shell calc distances between other relaxed shells
				for (k=j;k<nRelaxShells&&k<i;k++){
					//Need to check that this produces correct results
					tempSRDist->DebbyCalcSums(surfRelaxPos[j]->x,surfRelaxPos[j]->y,surfRelaxPos[j]->z, surfRelaxPos[j]->nAtoms, surfRelaxPos[k]->x, surfRelaxPos[k]->y, surfRelaxPos[k]->z, surfRelaxPos[k]->nAtoms, tempSRDist->mult, SRPREC, tempSRDist->nDist, j==k, 0);
				}
				// for each relaxed shell calc distances between core shells
				for (k=0;k<i-nRelaxShells;k++){
					tempSRDist->DebbyCalcSums(surfRelaxPos[j]->x,surfRelaxPos[j]->y,surfRelaxPos[j]->z, surfRelaxPos[j]->nAtoms, curPos->x, curPos->y, curPos->z, curPos->nAtoms, tempSRDist->mult, SRPREC, tempSRDist->nDist, false, 0);
					curPos=curPos->nextPos;
				}
			}
			// For each size weight the mult by the size distribution
			int startind=0;
			if (surfRelaxDist->dist[0]==0){
				surfRelaxDist->mult[0]+=p[i-1]*tempSRDist->mult[0];

				derivSRD_mu->mult[0]+=derivfSize_mu[i-1]*tempSRDist->mult[0];
				derivSRD_sig->mult[0]+=derivfSize_sig[i-1]*tempSRDist->mult[0];
				tempSRDist->mult[0]=0;
				startind++;
			}
			for (j=startind;j<tempSRDist->nDist;j++){
				surfRelaxDist->mult[j]+=2*p[i-1]*tempSRDist->mult[j];

				derivSRD_mu->mult[j]+=2*derivfSize_mu[i-1]*tempSRDist->mult[j];
				derivSRD_sig->mult[j]+=2*derivfSize_sig[i-1]*tempSRDist->mult[j];
				tempSRDist->mult[j]=0;

			}
			//memset(tempSRDist->mult, 0, tempSRDist->nDist*sizeof(double));
		}

		for (j=0; j<nRelaxShells&&j<i;j++){
			posPlaceHold[j]=posPlaceHold[j]->nextPos;
		}
	}
	// Consolidate the relaxation distances
	if (!err){
		surfRelaxDist->ConsolidateDistance();
		derivSRD_mu->ConsolidateDistance();
		derivSRD_sig->ConsolidateDistance();
	}
	return err;
}
*/
// Working on this function
//TODO Debug SurfaceRelaxation
int Shape::SurfaceRelaxation(Distance *totSRDist)
{
	int err=0, startSize=1, startind=0;
	Position *curPos;
	Position **surfRelaxPos, **posPlaceHold;
	Distance *tempSRDist;
#ifdef DEBUG
	cout<<"Entering SurfaceRelaxation\n";
#endif
	//Initialize temp surface relax position files
	try{
		surfRelaxPos=new Position* [nRelaxShells];
		for (int i=0;i<nRelaxShells;i++)
			surfRelaxPos[i]=new Position();
		posPlaceHold = new Position* [nRelaxShells];
		tempSRDist = new Distance();
	}
	catch (exception& e){
		printf("\nException found in SurfaceRelaxation alloc: %s\n", e.what());
		exit(1);
	}

	//Initialize surface relax distance files
	InitSRDistance(tempSRDist);
	InitSRDistance(totSRDist);

	//Change startSize
	if (distribution=="LogNorm"||distribution=="Gamma"){
		/*
		//startSize=floor(exp(mu->val-3*sig->val)/GetShapeScaleFactor());
		startSize=1;
		if (startSize<1||endSize>maxShell){
			printf("Needed range of Prob function is not possible with given distances\n");
			err=-1;
		}
		*/
	}
	else if (distribution=="Delta"){
		int shell= CalcDeltaShellSize();
		if (shell==-1){
			//Do not calculate SR in this case
			cout<<"Not Calculating SurfaceRelaxation, Delta size not supported.\n";
			startSize=1;
			endSize=startSize-1;
		}
		else{
			startSize=shell;
		}


	}
	else{
		cout<< "Size Distribution: "<<distribution<<" is not supported in Surface Relaxation.\n";
		err=-1;
	}
	//calculate SR and distances for this shape
	printf("Calculating surface relaxation for %s\n", name.c_str());

	//Initialize place holder
	for (int i=0;i<nRelaxShells;i++)
		posPlaceHold[i]=firstPos;

	//Copy positions and Apply surface relax
	for (int i=startSize; i<=endSize&&!err; i++){
		startind=0;
		if (distribution=="LogNorm"||distribution=="Gamma"){
			for (int j=0;j<nRelaxShells && j<i && !err;j++){
				CopyPosition(surfRelaxPos[j],posPlaceHold[j]);
				if (type=="sphere"||type=="sphere_Cerv"){
				//	surfRelaxPos[j]->ApplyPlanarSurfaceRelax(i, f->val, kappa->val);//surfRelaxPos[j]->ApplyRadialSurfaceRelax(i, GetShapeScaleFactor()/(2.0*a->val));
					if (n_R==0){
						surfRelaxPos[j]->ApplyNormExpRadialSurfRelax(i, f->GetVal(), kappa->GetVal(), GetShapeShellD()*10.0/2.0);
					}
					else{
						surfRelaxPos[j]->ApplyNormOscExpRadialSurfRelax(i, f->GetVal(), kappa->GetVal(), n_R->GetVal(), GetShapeShellD()*10.0/2.0);
					}
				}
				else
					surfRelaxPos[j]->ApplyPlanarSurfaceRelax(i, f->GetVal(), kappa->GetVal());
				surfRelaxPos[j]->ScalePos(a_D[i-1]);
			}
		}
		//Not Tested yet!!
		else if (distribution=="Delta"){
			Position *tempPos=firstPos;
			int shell= CalcDeltaShellSize();
			for (int j=0;j<shell-nRelaxShells; j++) tempPos=tempPos->nextPos;
			//for (int j=nRelaxShells-1; j>=0&&!err; j--){
			for (int j=0;j<nRelaxShells&&!err;j++){
				CopyPosition(surfRelaxPos[j], tempPos);
				if (type=="sphere"||type=="sphere_Cerv"){
					if (n_R==0){
						surfRelaxPos[j]->ApplyNormExpRadialSurfRelax(i, f->GetVal(), kappa->GetVal(), GetShapeShellD()*10.0/2.0);
					}
					else{
						surfRelaxPos[j]->ApplyNormOscExpRadialSurfRelax(i, f->GetVal(), kappa->GetVal(), n_R->GetVal(), GetShapeShellD()*10.0/2.0);
					}
				}
				else
					surfRelaxPos[j]->ApplyPlanarSurfaceRelax(i, f->GetVal(), kappa->GetVal());
				surfRelaxPos[j]->ScalePos(a_D[i-1]);
				tempPos=tempPos->nextPos;
			}
		}
		else{
			cout<<"\nERROR in SurfaceRelaxation()\n";
			err=-1;
		}
		//Calculate distances including surface relax shells
		if (!err){
			for (int j=0;j<nRelaxShells && j<i;j++){
				curPos=firstPos;
				// for each relaxed shell calc distances between other relaxed shells
				for (int k=j;k<nRelaxShells&&k<i;k++){
					tempSRDist->DebbyCalcSums(surfRelaxPos[j]->x,surfRelaxPos[j]->y,surfRelaxPos[j]->z, surfRelaxPos[j]->nAtoms, surfRelaxPos[k]->x, surfRelaxPos[k]->y, surfRelaxPos[k]->z, surfRelaxPos[k]->nAtoms, tempSRDist->mult, float(SRPREC), tempSRDist->nDist, j==k, 0);
				}
				// for each relaxed shell calc distances between core shells
				for (int k=0;k<i-nRelaxShells;k++){
					Position tempPos;
					CopyPosition(&tempPos, curPos);
					//Scale position so can use a(D) and surface relax at same time.
					tempPos.ScalePos(a_D[i-1]);
					tempSRDist->DebbyCalcSums(surfRelaxPos[j]->x,surfRelaxPos[j]->y,surfRelaxPos[j]->z, surfRelaxPos[j]->nAtoms, tempPos.x, tempPos.y, tempPos.z, tempPos.nAtoms, tempSRDist->mult, float(SRPREC), tempSRDist->nDist, false, 0);
					curPos=curPos->nextPos;
				}

			}
			// For each size weight the mult by the size distribution
			if (totSRDist->dist[0]==0){
				totSRDist->mult[0]+=p[i-1]*tempSRDist->mult[0];
				tempSRDist->mult[0]=0;
				startind=1;
			}
			for (int j=startind;j<tempSRDist->nDist;j++){
				totSRDist->mult[j]+=2*p[i-1]*tempSRDist->mult[j];
				tempSRDist->mult[j]=0;
			}

		}
		//Increase placeholder will not increase more holders than there are shells for the given size
		for (int j=0; j<nRelaxShells&& j<i;j++){
			posPlaceHold[j]=posPlaceHold[j]->nextPos;
		}
	}
	// Consolidate the relaxation distances
	if (!err){
		totSRDist->ConsolidateDistance();
	}

	delete tempSRDist;
	for (int i=0;i<nRelaxShells;i++)
		delete surfRelaxPos[i];
	delete [] surfRelaxPos;
	delete [] posPlaceHold;

	return err;
}

//Problem in release mode!!
void Shape::CopyPosition(Position* to, Position* from)
{
	to->nAtoms=from->nAtoms;
	to->shell= from->shell;
	to->nextPos=from->nextPos;

	if (to->x!=0)
		_mm_free (to->x);
	to->x= (float *) _mm_malloc (to->nAtoms*sizeof(float),16);

	if (to->y!=0)
		_mm_free (to->y);
	to->y= (float *) _mm_malloc (to->nAtoms*sizeof(float),16);

	if (to->z!=0)
		 _mm_free (to->z);
	to->z= (float *) _mm_malloc (to->nAtoms*sizeof(float),16);

	if (to->x==0||to->y==0||to->z==0){
		printf("\nException found in CopyPosition alloc: error in aligned memory alloc\n");
		exit(1);
	}

	for (int i=0;i<to->nAtoms;i++){
		to->x[i]=from->x[i];
		to->y[i]=from->y[i];
		to->z[i]=from->z[i];
	}
}
void Shape::CopyDistance(Distance * to, Distance *from)
{
	to->f=from->f;
	to->kappa=from->kappa;
	to->fileName=from->fileName;
	to->shell=from->shell;
	to->nDist=from->nDist;
	try{
		if (to->dist){
			delete [] to->dist;
		}
		if (to->mult){
			delete [] to->mult;
		}
		to->dist=new double [to->nDist];
		to->mult=new double [to->nDist];
	}
	catch(exception &e){
		cout<<"Error in CopyDistance: "<<e.what()<<endl;
		exit(-1);
	}
	for (int i=0;i<to->nDist;i++){
		to->dist[i]=from->dist[i];
		to->mult[i]=from->mult[i];
	}
}
void Shape::InitSRDistance(Distance *currSRDist)
{
	double maxDistSq, maxDist;
	double delD=GetShapeShellD()*10/2.0;
	Distance *currDist=firstDist;
	for (int i=0;i<endSize-1;i++) currDist=currDist->nextDist;
	if (currDist->shell==endSize){
		maxDist=currDist->dist[currDist->nDist-1];
		maxDist/=2.0;
		maxDist=2*sqrt(maxDist*maxDist+2*f->GetVal()*delD*maxDist+f->GetVal()*f->GetVal()*delD*delD);
		//Scale by lattice parameter so can use a_D with SR.
		maxDist*=a_D[endSize-1];
	}
	else {
		printf("Error in finding maxDist for InitSRDistance\n");
		exit(0);
	}
	maxDistSq=(maxDist+2)*(maxDist+2);
	currSRDist->InitDistance((int)(ceil(maxDistSq/SRPREC)));
}

// TODO make recalculation of diffraction pattern slightly more efficient by using preloaded intensities
// As long as lattice param has not changed.
int Shape::CalcPattern()
{
	int err=0;
	if (!diffPatt){
		try{
			diffPatt=new Pattern* [Pattern::nLambda];
			for (int i=0;i<Pattern::nLambda;i++){
				diffPatt[i]=new Pattern(name);
			}
		}
		catch (exception& e){
			printf("\nException found in Shape diffPatt initialization: %s\n", e.what());
			return -1;
		}
	}
	// Initialize lattice parameter for different size particles
	InitLatticeParam();
	// Initialize the size distribution
	CalcSizeDistribution();
	// Load any pre-calculated intensities from file
	if (nPreCalcPatts>0){
		LoadPatterns();
	}
	// Calculated then number of SurfaceRelaxation shells and the wieghts of each shell
	nRelaxShells=GetNumSRShells();
	CalcShellWeights();
	cout<<"maxShell: " <<maxShell<<" nRelaxShells: "<<nRelaxShells<<endl;
	// Calculate the SurfaceRelaxation distances
	if (nRelaxShells>0){
		if (surfRelaxDist==0){
			try{
				surfRelaxDist=new Distance();
			}
			catch (exception& e){
				printf("\nException found in Shape->CalcPattern: %s\n", e.what());
				exit(1);
			}
		}

		err=SurfaceRelaxation(surfRelaxDist);
	}
	else {
		if (surfRelaxDist!=0){
			delete surfRelaxDist;
			surfRelaxDist=0;
		}

	}
	if (!err){
		//Calc pattern for each lambda
		for (int j=0;j<Pattern::nLambda&&!err;j++){
			if(preCalcPatts==0){
				//Routines to use when preCalculated patterns are NOT loaded.
				if (nParams_a>1){
					//Calculates for each size separately, scaling by a(D)
					diffPatt[j]->CalcPattern(firstDist,endSize, p, a_D, nRelaxShells, surfRelaxDist, j);
				}
				else{
					//Faster when only using constant lattice parameter because using wShell
					diffPatt[j]->CalcPattern(firstDist, endSize, wShell, a[0]->GetVal(), surfRelaxDist, j);
				}
			}
			else{
				//Routine to use when precalculated patterns are loaded.
				diffPatt[j]->CalcPattern(firstDist,endSize, p, a_D, nRelaxShells, surfRelaxDist, j, preCalcPatts);
			}
		}
	}
	return err;
}

//"Diameters" of shapes
double Shape::GetShapeShellD()
{
	//Returns the maximum distance stepping for a given shape.
	//Currently unitless.
	//when multiplied by the lattice parameter in Ang,=>D(nm)
	if (type=="tetrahedron") return (2*sqrt(2.0)/10.0);
	else if (type=="octahedron") return (2/10.0);
	else if (type=="cuboctahedron") return  (sqrt(2.0)/10.0);
	else if (type=="cube") return (sqrt(3.0)/10.0);
	else if (type=="sphere") return (sqrt(2.0)/10.0);//NN Def
	else if (type=="sphere_Cerv") return (.7815926418);//Cerv Def
	else if (type=="sphereDefFaultAtCenter") return (sqrt(2.0)/10.0);
	else if (type=="sphereTwinAtCenter") return (sqrt(2.0)/10.0);
	else if (type=="cuboid1x1_1x1_21") return (sqrt(3.31)/10.0);
	else if (type=="cuboid1x1_2x1_44") return (sqrt(3.64)/10.0);
	else if (type=="icosahedron") return (sqrt(2.0)*(0.95105)*(1.05146)/10.0);
	else if (type=="decahedron") return (sqrt(2.0)*(0.8506508)/10.0);
	else if (type=="centroDecahedron") return (2*sqrt(2.0)*(0.8506508)/10.0);
	else return 0;
}
//"Edge Lengths" of shapes
double Shape::GetShapeShellL()
{
	//Returns the maximum distance stepping for a given shape.
	//Currently unitless.
	//when multiplied by the lattice parameter =>D(units of a)
	if (type=="tetrahedron") return (2*sqrt(2.0));
	else if (type=="octahedron") return ((sqrt(2.0)));
	else if (type=="cuboctahedron") return (1.0/(sqrt(2.0)));
	else if (type=="cube") return (1.0);
	else if (type=="sphere") return (sqrt(2.0));//NN Def
	else if (type=="sphereDefFaultAtCenter") return (sqrt(2.0));
	else if (type=="sphereTwinAtCenter") return (sqrt(2.0));
	else if (type=="cuboid1x1_1x1_21") return (1.0);
	else if (type=="cuboid1x1_2x1_44") return (1.0);
	else if (type=="icosahedron") return (1.05146/(sqrt(2.0)));
	else if (type=="decahedron") return (1.0/(sqrt(2.0)));
	else if (type=="centroDecahedron") return (sqrt(2.0));
	//else if (type=="sphere_Cerv") return (.7815926418);//Cerv Def
	else return 0;
}

double Shape::CalcTotVol()
{
	double vTot=0;
	int i;
	for (i=1;i<=maxShell; i++)
	{
		vTot+=p[i-1]*GetUnitlessVolume(i)*a_D[i-1]*a_D[i-1]*a_D[i-1];
	}
	return vTot;
}
double Shape::GetUnitlessVolume(int n)
{
	//Assume Lattice Parameter is unity
	//Actual volume is then (Unitless volume) * a^3 (units of a^3)
	double vol=0, sideLength=n*GetShapeShellL();
	if (type=="tetrahedron"){
		vol=sqrt(2.0)*sideLength*sideLength*sideLength/12.0;
	}
	else if (type =="octahedron"){
		vol=sqrt(2.0)*sideLength*sideLength*sideLength/3.0;
	}
	else if (type=="cuboctahedron"){
		vol=5*sqrt(2.0)*sideLength*sideLength*sideLength/3.0;
	}
	else if (type=="sphere"||type=="sphere_Cerv"||type=="sphereDefFaultAtCenter"||type=="sphereTwinAtCenter"){
		vol=PI*sideLength*sideLength*sideLength/6.0;
	}
	else if (type=="cube"){
		vol=sideLength*sideLength*sideLength;
	}
	else if (type=="cuboid1x1_1x1_21"){
		vol=sideLength*1.1*sideLength*1.21*sideLength;
	}
	else if (type=="cuboid1x1_2x1_44"){
		vol=sideLength*1.2*sideLength*1.44*sideLength;
	}
	else if (type=="icosahedron"){
		vol=5*(3+sqrt(5.0))*sideLength*sideLength*sideLength/12.0;
	}
	else if (type=="decahedron"){
		vol=(5+sqrt(5.0))*sideLength*sideLength*sideLength/12.0;
	}
	else if (type=="centroDecahedron"){
		vol=(5+sqrt(5.0))*sideLength*sideLength*sideLength/12.0;
	}
	return vol;
}
void Shape::CalcNumAtoms()
{
	Position *curPos=firstPos;
	if (!nA_D){
		try{
			nA_D=new int [maxShell];
		}
		catch(exception &e){
			cout<<"Exception found in Shape::GetNumAtoms: "<<e.what()<<endl;
		}
		for (int i=0;i<maxShell;i++){
			nA_D[i]=0;
		}
		int total=0;
		for (int i=0;i<maxShell;i++){
			total+=curPos->nAtoms;
			nA_D[i]=total;
			curPos=curPos->nextPos;
		}
	}
	nAtoms=0;
	for (int i=0;i<maxShell;i++){
		nAtoms+=nA_D[i]*p[i];
	}
	printf("Number of atoms: %f\n", nAtoms);
}


void Shape::CalcSizeDistribution()
{
	if (!a_D){
		cout<<"Error in CalcWeights: Lattice Parameter array is not initialized\n";
		exit(0);
	}
	//Use lattice parameter in nanometers
	double *a_nm;
	try{
		if (p){
			delete [] p;
		}
		p = new double [maxShell];
		a_nm = new double [maxShell];
	}
	catch (exception &e){
		cout<<"\nException found in CalcWeights: "<<e.what()<<endl;
		exit(1);
	}
	//Initialize arrays
	for (int i=0;i<maxShell;i++){
		p[i]=0;
		a_nm[i]=a_D[i]/10.0;
	}
	//Check for Diameter Array
	if (!D) {
		InitDiameterArray();
	}
	CheckShellThickness();
	//Calculate probabilty distribution function (p[i])
	if (distribution=="Delta") DeltaDist(&p[0]);
	else if(distribution=="LogNorm") LogNormalDist(&D[0], &p[0], &a_nm[0]);
	else if (distribution=="Gamma") GammaDist(&D[0], &p[0], &a_nm[0]);
	else {
		cout<<endl<<"Error: Unsupported Size distribution: "<<distribution<<endl;
		exit(0);
	}
	NormalizeSizeDistribution();
	// Get endsize  by looking at volume*weights
	CalcEndSize();
	//Check Number of Surf Relax Shells.
	cout<<"Max size cutoff implied: "<<endSize<<endl;
	cout<<"Calculating shell weights for "<< name<<endl;

	try{
		delete [] a_nm;
	}
	catch(exception &e){
		cout<<"Error in CalcSizeDistributions: "<<e.what()<<"\n";
		exit(1);
	}
}
void Shape::CalcShellWeights()
{
	try{
			if (wShell){
				delete [] wShell;
			}
			wShell = new double [maxShell];
	}
	catch (exception &e){
		cout<<"\nException found in CalcWeights: "<<e.what()<<endl;
		exit(1);
	}
	for (int i=0;i<maxShell;i++){
		wShell[i]=0;
		//Weights for shells without SR
		for (int j=i+nRelaxShells;j<endSize;j++){
			wShell[i]+=p[j];
		}
	}
}
void Shape::NormalizeSizeDistribution()
{
	double tot=0;
	for (int i=0;i<maxShell; i++){
		tot+=p[i];
	}
	for (int i=0;i<maxShell; i++){
		p[i]/=tot;
	}
}
int Shape::WriteIntFileHeader(string _fPath, string _fName)
{
	int err=0;
	string filename=_fPath+_fName;
	ifstream infile;
	infile.open(filename.c_str());
	if (infile.is_open()) printf("Intensity file %s exists.. Appending Header\n", _fName.c_str());
	else printf("Creating Intensity file %s and writing header\n", _fName.c_str());
	infile.close();
	fstream file(filename.c_str(), ios_base::out|ios_base::app);  //|ios_base::trunc);
	file.setf(ios::fixed, ios::floatfield);
	if (file.is_open()){
		if (element!="") file<<"ELEM	";
		file<<"SHAPE	A	DISTRIBUTION	";
		if (distribution=="Delta") file<<"SHELLS	";
		else if (distribution=="LogNorm") file<<"MU	  SIGMA	";
		if (f!=0&&kappa!=0){
			if (f->GetVal()!=0) file<<"SR_F	 SR_KAPPA	";
		}
		file<<"WEIGHT";
		file<<endl;
		if(element!="") file<<element<<"	";
		file<<type<<"		"<<a_D[maxShell-1]<<"		"<<distribution<<"		";
		if (distribution=="Delta") file<<deltaSize->GetVal()<<"	";
		else if (distribution=="LogNorm") file<<mu->GetVal()<<"	"<<sigma->GetVal()<<"		";
		if (f!=0&&kappa!=0){
			if (f->GetVal()!=0) file<<f->GetVal()<<"	"<<kappa->GetVal()<<"		";
		}
		file<<wShape->GetVal();
		file<<endl;
	}
	else{
		printf ("Error writing to file %s\n", filename.c_str());
		err=-1;
	}
	file.close();
	return err;
}
void Shape::CalcEndSize()
{
	bool stop=false;
	int n=0;
	double *volWeightSizeDist, maxVal=0;
	Position *curPos=firstPos;
	try{
		volWeightSizeDist=new double [maxShell];
	}
	catch (exception& e){
		printf("\nException found in CalcEndSize: %s\n", e.what());
		exit(1);
	}

	for (int i=0;i<maxShell;i++){
		n+=curPos->nAtoms;
		volWeightSizeDist[i]=p[i]*n;
		maxVal=getmax(maxVal, volWeightSizeDist[i]);
		curPos=curPos->nextPos;
	}
	endSize=maxShell;
	for (int i=(maxShell-1);i>=0&&!stop; i--)
	{
		if (volWeightSizeDist[i]<(.001*maxVal)){
			endSize--;
		}
		else{
			stop=true;
		}
	}

	delete [] volWeightSizeDist;
}

double Shape::GammLn(const double xx)
{
	// Returns the value ln(Gamma(xx)) for xx>0
	//Taken from Numerical Recipes 3rd Edition pg. 257
	int j;
	double x, tmp, y, ser;
	static const double cof[14]={57.1562356658629235, -59.5979603554754912,
			14.1360979747417471, -0.491913816097620199, .339946499848118887e-4,
			.465236289270485756e-4,-.983744753048795646e-4, .158088703224912494e-3,
			-.210264441724104883e-3, .217439618115212643e-3, -.164318106536763890e-3,
			.844182239838527433e-4, -.261908384015814087e-4, .368991826595316234e-5};
	if (xx<0) throw("bad arg in GammLn");
	y=x=xx;
	tmp = x+5.2421875000;
	tmp = (x+0.5)*log(tmp)-tmp;
	ser=0.999999999999997092;
	for (j=0;j<14;j++) ser+=cof[j]/++y;
	return tmp+log(2.5066282746310005*ser/x);
}
void Shape::GammaDist(double *d, double *P, double *scale)
{
	//Param *alpha=0, *beta=0;
	if (!alpha||!beta){
		cout<<"Uninitialized alpha and beta, cannot calculate Gamma distribution"<<endl;
		exit(0);
	}
	double delD=(d[1]-d[0]);
	if (delD==0){
		cout<<"Error: shell diameter array in GammaDist"<<endl;
		exit(0);
	}
	double c1=alpha->GetVal()*log(beta->GetVal()), c2=GammLn(alpha->GetVal());
	for(int i=0;i<maxShell; i++){
		P[i]=c1-c2+(alpha->GetVal()-1)*log(d[i]*scale[i])-(beta->GetVal()*d[i]*scale[i]);
		P[i]=exp(P[i]);
		if (i==0){
			delD=d[i]*scale[i];
		}
		else{
			delD=d[i]*scale[i]-d[i-1]*scale[i-1];
		}
		P[i]*=delD;
	}
}
void Shape::LogNormalDist(double *d, double *P, double *scale)
{
	if (!mu||!sigma){
		cout<<"Uninitialized mu and sigma, cannot calculate Lognormal distribution"<<endl;
		exit(0);
	}
	double delD=(d[1]-d[0]);
	if (delD==0){
		cout<<"Error: shell diameter array in LogNormalDist"<<endl;
		exit(0);
	}
	//For other shapes other than spheres Diameter is replaced by the characteristic size parameter (L) in the shape before any strain is applied.
	for (int i=0;i<maxShell;i++){
		P[i]=exp(-(log(d[i]*scale[i])-mu->GetVal())*(log(d[i]*scale[i])-mu->GetVal())/(2*sigma->GetVal()*sigma->GetVal()))/(d[i]*scale[i]*sigma->GetVal()*sqrt(2*PI));
		if (i==0){
			delD=d[i]*scale[i];
		}
		else{
			delD=d[i]*scale[i]-d[i-1]*scale[i-1];
		}
		P[i]*=delD;
	}
}
void Shape::DeltaDist(double *P)
{
	int shell=CalcDeltaShellSize();

	for (int i=0;i<maxShell;i++){
		P[i]=0.0;
	}
	if (shell==-1){
		cout<<"\nWARNING: Particle size "<<deltaSize->GetVal()<<" is not supported by Delta Distribution, setting all size probabilities to 0.\n";
	}
	else{
		//Adjust the particle size to a multiple of the shell thickness for clarity.
		deltaSize->ChangeVal(shell*a[0]->GetVal()*GetShapeShellL());
		cout<<"Delta Size Set to: "<<deltaSize->GetVal()<<endl;
		P[shell-1]=1.0;
	}
}
void Shape::InitLatticeParam()
{
	if (!a_D){
		try{
			a_D=new double [maxShell];
		}
		catch(exception &e){
			cout<<"Exception in InitLatticeParameter: "<<e.what()<<endl;
			exit(1);
		}
	}
	if (a[0]->name.find("Poly")!=string::npos){
		PolynomialLatticeParam();
	}
	else{
		ConstantLatticeParam();
	}
}
// Generates array of lattice parameter following polynomial a+bx+cx^2
// Diameter of sphere with equivalent volume is used a size parameter
// This is so same lattice parmeter array can be applied to multiple shapes
void Shape::PolynomialLatticeParam()
{
	for (int i=0;i<maxShell;i++){
		double d=exp(log(GetUnitlessVolume(i+1)*6/PI)/3.0);
		a_D[i]=a[0]->GetVal();
		double temp=d;
		for (int j=1; j<nParams_a;j++){
			a_D[i]+=a[j]->GetVal()*temp;
			temp*=d;
		}
	}
}

void Shape::ConstantLatticeParam()
{
	for (int i=0;i<maxShell;i++){
		a_D[i]=a[0]->GetVal();
	}
}
//
void Shape::InitDiameterArray()
{
	try{
		D=new double [maxShell];
	}
	catch(exception &e){
		cout<<"Exception found in InitDiameterArray: "<<e.what()<<endl;
		exit(0);
	}
	//double deltaD=GetShapeShellD();
	double deltaD=GetShapeShellL();
	if (!deltaD){
		cout<<"Unrecognized shape name in GetShapeShellD: "<<name<<endl;
		exit(0);
	}
	for (int i=0;i<maxShell;i++){
		D[i]=(i+1)*deltaD;
	}
}
void Shape::CalcDistForShell(Distance *curDist)
{
	float minx=0, miny=0, minz=0, maxx=0, maxy=0, maxz=0;
	Position *outsideShell=firstPos;
	for (int i=1;i<curDist->shell;i++){
		outsideShell=outsideShell->nextPos;
	}
	minx=outsideShell->x[0];
	maxx=minx;
	miny=outsideShell->y[0];
	maxy=miny;
	minz=outsideShell->z[0];
	maxz=minz;
	for (int i=1;i<outsideShell->nAtoms;i++){
		maxx=getmax(maxx,outsideShell->x[i]);
		minx=getmin(minx,outsideShell->x[i]);
		maxy=getmax(maxy,outsideShell->y[i]);
		miny=getmin(miny,outsideShell->y[i]);
		maxz=getmax(maxz,outsideShell->z[i]);
		minz=getmin(minz,outsideShell->z[i]);
	}
	if (type=="decahedron"){
		Position *tempShell=firstPos;
		for (int j=1;j<curDist->shell-1;j++){
			for (int i=1;i<tempShell->nAtoms;i++){
				maxx=getmax(maxx,outsideShell->x[i]);
				minx=getmin(minx,outsideShell->x[i]);
				maxy=getmax(maxy,outsideShell->y[i]);
				miny=getmin(miny,outsideShell->y[i]);
				maxz=getmax(maxz,outsideShell->z[i]);
				minz=getmin(minz,outsideShell->z[i]);
			}
			tempShell=tempShell->nextPos;
		}
	}
	curDist->nDist=(int) ceil(((maxx-minx)*(maxx-minx)+(maxy-miny)*(maxy-miny)+(maxz-minz)*(maxz-minz))/DISTPREC)+1;
	try{
		curDist->dist= new double [curDist->nDist];
		curDist->mult=new double [curDist->nDist];
	}
	catch (exception &e){
		cout<<"Error in distance array memory allocation.\n";
		exit(1);
	}
	for (int i=0;i<curDist->nDist;i++){
		curDist->dist[i]=sqrt(i*DISTPREC);
		curDist->mult[i]=0;
	}
	Position *curPos=firstPos;
	for (int i=1;i<curDist->shell;i++){
		curDist->DebbyCalcSums(curPos->x, curPos->y, curPos->z, curPos->nAtoms, outsideShell->x, outsideShell->y, outsideShell->z, outsideShell->nAtoms, curDist->mult, DISTPREC, curDist->nDist, false, 0);
		curPos=curPos->nextPos;
	}
	curDist->DebbyCalcSums(curPos->x, curPos->y, curPos->z, curPos->nAtoms, outsideShell->x, outsideShell->y, outsideShell->z, outsideShell->nAtoms, curDist->mult, DISTPREC, curDist->nDist, true, 0);
	curDist->ConsolidateDistance();
	for (int i=1; i<curDist->nDist;i++ ){
		curDist->mult[i]*=2.0;
	}
}
int Shape::GetNumSRShells()
{
	// Determine if there is surface relaxation
	// Only apply to a shell if the amplitude is as large as the distance precision
	if (f!=0&&kappa!=0){
		if (f->GetVal()!=0){
			//TODO Check that this SR function is correct..
			double shellT=GetShapeShellD()*10.0*a_D[maxShell-1]/2.0;
			if (shellT==0){
				cout<<"Warning: Shape not supported in GetNumSRShells routine, using lattice parameter as shell thickness...\n";
				shellT=a_D[maxShell-1];
			}
			//int nSR=(int)ceil(kappa->GetVal()*log(shellT*abs(f->GetVal())/SRPREC));
			//For now just calculate all shells
			int nSR=endSize;
			if (nSR>maxShell){
				cout<<"Warning: surface relaxation extends throughout particle and does not go to zero at center. \n";
				nSR=maxShell;
			}
			if (nRelaxShells>endSize){
				nRelaxShells=endSize;
				cout<<"Number of Surf. Relax Shells is larger than largest considered particle..\n";
				cout<<"Setting Number of Surf. Relax Shells to: "<<nRelaxShells<<endl;
			}
			else if (nSR<2){
				cout<<"Warning: number of surface relaxation shells <2 !!!!"<<endl;
				cout<<"Setting number of surface relaxation shells = 2\n";
				nSR=2;
			}
			else{
				cout<<"Number of Surface Relaxed Shells = "<<nSR<<endl;
			}

			return nSR;
			//return (int)ceil(-kappa->GetVal()*log(SRPREC));
		}
		else {
			cout<<"No Surface Relaxation applied (f=0)\n";
			return 0;
		}
	}
	else{
		cout<<"No Surface Relaxation applied (f and/or kappa are uninitialized).\n";
		return 0;
	}
}
void Shape::CheckSupportedShapes()
{
	if (type=="tetrahedron")
		;
	else if (type=="cuboctahedron")
		;
	else if (type=="octahedron")
		;
	else if (type=="sphere")
		;
	else if (type=="sphere_Cerv")
		;
	else if (type=="cube")
		;
	else if (type=="sphereDefFaultAtCenter")
		;
	else if (type=="sphereTwinAtCenter")
		;
	else if (type=="cuboid1x1_1x1_21")
		;
	else if (type=="cuboid1x1_2x1_44")
		;
	else if (type=="icosahedron")
		;
	else if (type=="decahedron")
		;
	else if (type=="centroDecahedron")
		;
	else {
		cout<<"Unsupported Shape Type: "<<type<<endl;
		exit(0);
	}
}
void Shape::CheckShellThickness()
{
	for (int i=1;i<maxShell;i++){
		if ((D[i]*a_D[i])<(D[i-1]*a_D[i-1])){
			cout<<"Error in Shape::CheckShellThickness(), lattice parameter results in negative shell thickness. "<<i<<endl;
			exit(0);
		}
	}
}
int Shape::CalcDeltaShellSize()
{
	//Assuming that lattice parameter is independent of size.
	double N=deltaSize->GetVal()/(a[0]->GetVal()*GetShapeShellL());
	//Consider +/- 0.5
	int shell=getmax((int)N, (int)(N+0.5));
	if (shell<1||shell>maxShell){
		return -1;
	}
	else{
		return shell;
	}
}
int Shape::GetNumSizes()
{
	int nSizes=0;
	Position *curPos=firstPos;
	while(curPos!=0){
		nSizes++;
		curPos=curPos->nextPos;
	}
	return nSizes;
}
int Shape::GetIndexOfShell(int shell)
{
	int index=0;
	Position *curPos=firstPos;
	while(curPos!=0){
		if (curPos->shell==shell){
			return index;
		}
		else{
			index++;
			curPos=curPos->nextPos;
		}
	}
	return -1;
}
void Shape::LoadIntFromFile(Pattern *patts, double aD)
{
	string lattice, shape;
	int shell, nI, nAtoms;
	double minSa, maxSa, dsa;
	double *inSa, *inI;
	ifstream calcedPatt;
	//Read the important info from the first line of the intensity file
	calcedPatt.open(patts[0].fileName.c_str());
	if (calcedPatt.is_open()){
		string firstLine;
		getline(calcedPatt, firstLine);
		stringstream test(firstLine);
		while (!test.eof()){
			string stuff;
			test>>stuff;
			if (stuff=="shell:"){
				test>>shell;
			}
			else if(stuff=="latt:"){
				test>>lattice;
			}
			else if (stuff=="shape:"){
				test>>shape;
			}
			else if (stuff=="nAtoms:"){
				test>>nAtoms;
			}
			else if (stuff=="minS*a:"){
				test>>minSa;
			}
			else if (stuff=="maxS*a;"){
				test>>maxSa;
			}
			else if (stuff=="d(s*a):"){
				test>>dsa;
			}
			else if (stuff=="nInt:"){
				test>>nI;
			}
			else{

			}
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
