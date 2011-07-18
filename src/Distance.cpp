#include "Distance.h"

Distance::Distance(char *_fileName)
{
	fileName=_fileName;
	shell=0;
	f=0;
	kappa=0;
	nDist=0;
	dist=0;
	mult=0;
	nextDist=0;
}

Distance::Distance()
{
	fileName="";
	shell=0;
	f=0;
	kappa=0;
	nDist=0;
	dist=0;
	mult=0;
	nextDist=0;
}

Distance::~Distance(void)
{
	delete [] dist;
	delete [] mult;
}

void Distance::ReadDistFile()
{
	string tempstring1;
	int counter=1;
	ifstream infile;
	infile.open(fileName.c_str());
	if (infile.is_open()){
		//printf("file %s opened\n", fileName.c_str());
		while ((!infile.eof()))
		{
			infile >> tempstring1;
			if (tempstring1=="SHAPE")
				infile >> tempstring1>>tempstring1>>shell;
			//Input surface relax parameters
			else if (tempstring1=="SURFACE"){
				infile >>tempstring1>>f>>kappa;
			}
			//Input number of distances in file
			else if (tempstring1=="NUMBER"){
				infile >> tempstring1>>tempstring1>>nDist;
				//allocate mem for distances and mult
				try{
					dist=new double [nDist];
					mult=new double [nDist];
				}
				catch (exception& e){
					printf("\nException found in ReadDistFile alloc: %s\n", e.what());
					exit(-1);
				}
			}
			//Input distances and multiplicites
			else if (tempstring1=="DISTANCES"){
				for (int i=0; i<nDist; i++){
					infile >>dist[i]>>mult[i];
				}
			}
			else {
				//If there is an error with reading the data from the file an error message will be displayed and
				//	the loop index that the error occured on.
				/*
				if (counter <5){
					printf("Error in reading file: %s\n", fileName.c_str());
					printf("Error with data input: cycle %s\n",tempstring1.c_str());
					infile.close();
					err= -1;
				}
				*/
				if (counter <5){
					cout<<"Error in reading file: "<<fileName<<endl;
					cout<<"Unrecognized keyword: " <<tempstring1<<endl;
					infile.close();
					exit (1);
				}

			}
			counter++;
		}
		infile.close();
	}
	else {
		printf("Could not open file: %s \n", fileName.c_str());
		exit(1);
	}
}


int Distance::CalcDistances(int _shell, Position *firstPos)
{
	int i, j, k, nTot=0, n1, n2, nD, place;
	Position *curPos=firstPos;
	float *x1, *y1, *z1, *x2, *y2, *z2;
	double *d, *m, *tempd, *tempm;
	//Calc max number of possible dist
	shell=_shell;
	for (i=1;i<shell;i++){
		nTot+=curPos->nAtoms;
		curPos=curPos->nextPos;
	}
	x1=curPos->x;
	y1=curPos->y;
	z1=curPos->z;
	n1=curPos->nAtoms;
	nD=(n1*(n1-1)/2)+n1*nTot+1;
	try {
		d=new double [nD];
		m=new double [nD];
	}
	catch (exception& e){
		printf("\nException found in CalcDistances alloc: %s\n", e.what());
		return -1;
	}
	//Calc distances within shell
	x2=curPos->x;
	y2=curPos->y;
	z2=curPos->z;
	d[0]=0;
	m[0]=double(n1);
	place=1;
	for (i=0;i<n1;i++){
		for (j=i+1;j<n1;j++){
			d[place]=(x1[i]-x2[j])*(x1[i]-x2[j])+(y1[i]-y2[j])*(y1[i]-y2[j])+(z1[i]-z2[j])*(z1[i]-z2[j]);
			m[place]=2;
			place++;
		}
	}
	//Calc distances btw shell and smaller shells
	curPos=firstPos;
	for (k=1;k<shell;k++){
		n2=curPos->nAtoms;
		x2=curPos->x;
		y2=curPos->y;
		z2=curPos->z;
		for (i=0;i<n1;i++){
			for (j=0;j<n2;j++){
				d[place]=(x1[i]-x2[j])*(x1[i]-x2[j])+(y1[i]-y2[j])*(y1[i]-y2[j])+(z1[i]-z2[j])*(z1[i]-z2[j]);
				m[place]=2;
				place++;
			}
		}
		curPos=curPos->nextPos;
	}
	//compress distances
	QuickSortImproved(d,m,0,nD-1);
	nDist=MultiplicitiesCompute(d,m,nD);
	try{
		tempd=new double [nDist];
		tempm=new double [nDist];
		for (i=0;i<nDist;i++){
			tempd[i]=d[i];
			tempm[i]=m[i];
		}
		delete [] d;
		delete [] m;
		dist=new double [nDist];
		mult=new double [nDist];
	}
	catch (exception& e){
		printf("\nException found in CalcDistances alloc: %s\n", e.what());
		return -1;
	}
	//Set the
	for (i=0;i<nDist;i++){
		dist[i]=sqrt(tempd[i]);
		mult[i]=tempm[i];
	}

	delete [] tempd;
	delete [] tempm;

	return 0;
}

void Distance:: QuickSortImproved(double *x, double *y, int lb, int ub)
{
	 while (lb < ub) {
        int m;

        /* quickly sort short lists */
        if (ub - lb <= 50) {
            this->sortByInsertion(x, y, lb, ub);
            return;
        }

        m = partition(x, y, lb, ub);

        /* eliminate tail recursion and */
        /* sort the smallest partition first */
        /* to minimize stack requirements    */
        if (m - lb <= ub - m) {
            QuickSortImproved(x, y,lb, m);
            lb = m + 1;
        } else {
            QuickSortImproved(x, y,m + 1, ub);
            ub = m;
        }
    }
}
void Distance :: sortByInsertion(double *x, double *y, int lb, int ub) {
    int i, j;

    for (i = lb + 1; i <= ub; i++) {
        double t = x[i];
		double u = y[i];

        /* shift down until insertion point found */
		for (j = i-1; j >= lb && (x[j]> t); j--){
            x[j+1] = x[j];
			y[j+1] = y[j];
		}

        /* insert */
        x[j+1] = t;
		y[j+1] = u;
    }
}

int Distance :: partition(double *x, double *y, int lb, int ub) {

    /* select a pivot */
    double pivot = x[(lb+ub)/2];

    /* work from both ends, swapping to keep   */
    /* values less than pivot to the left, and */
    /* values greater than pivot to the right  */
    int i = lb - 1;
    int j = ub + 1;
    while (1) {
        double t;
		double u;

        while ((x[--j]> pivot));
        while ((x[++i]< pivot));
        if (i >= j) break;

        /* swap x[i], x[j] */
        t = x[i];
		u = y[i];
        x[i] = x[j];
		y[i]=y[j];
        x[j] = t;
		y[j]=u;
    }

    return j;
}

int Distance::MultiplicitiesCompute(double *d, double *m, int nD)
{
	int i,j=0;
	double prevVal, currVal=0;

	prevVal=d[0];
	for (i=1; i<nD;i++){
		currVal=d[i];
		if (currVal<=(prevVal+PREC))
		{
			m[j]+=m[i];
		}
		else
		{
			j++;
			d[j]=currVal;
			m[j]=m[i];
			prevVal=currVal;
		}
	}
	return (j+1);
}
void Distance::OutputDistToFile(string path, string shape)
{
	string _file = path + fileName;
	fstream file(_file.c_str(), ios_base::out|ios_base::trunc);
	if (file.is_open()){
		file.setf(ios::fixed, ios::floatfield);
		file << "SHAPE		SHELL\n";
		file << shape<<"	"<<shell<<"\n";
		file << "SURFACE RELAXATION\n";
		file<<f<<"		"<<kappa<<"\n";
		file <<"NUMBER OF DISTANCES\n";
		file<<nDist<<"\n";
		file<< "DISTANCES" <<"\n";
		file.precision(5);
		for (int i=0;i<nDist;i++){
			file << dist[i] << " ";
			file << mult[i]<< "\n";
		}
		file.close();
	}
	else {
		printf("Error writing file... check path.");
		exit (1);
	}
}

void Distance::DebbyCalcSums(float* pArray1x,
			  float* pArray1y,
			  float* pArray1z,
			  int p1Size,
			  float* pArray2x,
			  float* pArray2y,
			  float* pArray2z,
			  int p2Size,
			  double* pResult,                   // [out] result array
			  double stepd,						// [in] step size in d
			  int nSize,						// [in] size of result array
			  bool sameArray,                   // [in] if pArray1==pArray2
			  float startd)

{
#ifndef NOUSESTEP
	float istep=1.0f/float(stepd);
	__m128  pvals = _mm_load1_ps(&istep);
	__m128  d0= _mm_load1_ps(&startd);
#endif

#ifdef DBLSPD
	int DSPD = 2;
#else
	int DSPD = 1;
#endif

	int DIV = 4*DSPD;
	int p2optsize = p2Size/DIV;
	int p2optrem  = p2Size%DIV;
	int maxBegVect = p1Size-p2optrem;
	int flag=0;



//#pragma omp parallel for //schedule(static, 2)
	for(int i=0;i<p1Size;i++)
	{
		__m128* p2x = (__m128*) pArray2x;
		__m128* p2y = (__m128*) pArray2y;
		__m128* p2z = (__m128*) pArray2z;
		__m128 p1x = _mm_load1_ps(pArray1x++);
		__m128 p1y = _mm_load1_ps(pArray1y++);
		__m128 p1z = _mm_load1_ps(pArray1z++);

		int ibeg;

		if(sameArray){
#ifdef JUSTENOUGH
			int irem = i%DIV;
			ibeg = i/DIV;

			p2x+=(DSPD*ibeg);
			p2y+=(DSPD*ibeg);
			p2z+=(DSPD*ibeg);

			// LEFTOVERS: BEGINNING OF VECTOR (if i>0 and not %4)
			// The second control is to impliment an upper bound of when to use the beginnings (do not need in the bottom corner of the distance matrix)
			if (irem>0&&(i<maxBegVect))
			{
				ibeg++;

			// get 4 deltacoords

				__m128 dx,dy,dz, d2;
				__m128i jumps;
				if (irem<4){
				dx = _mm_sub_ps(p1x, *(p2x++));
				dy = _mm_sub_ps(p1y, *(p2y++));
				dz = _mm_sub_ps(p1z, *(p2z++));
				// square

				dx = _mm_mul_ps(dx, dx);
				dy = _mm_mul_ps(dy, dy);
				dz = _mm_mul_ps(dz, dz);

				// sum dx+dy+dz

				d2 = _mm_add_ps(dx, dy);
				d2 = _mm_add_ps(d2, dz);

				//subtract the offset of the distances

				d2 = _mm_sub_ps(d2,d0);

				}
				else{
					p2x++;
					p2y++;
					p2z++;
				}


	#ifdef DBLSPD
				dx = _mm_sub_ps(p1x, *(p2x++));
				dy = _mm_sub_ps(p1y, *(p2y++));
				dz = _mm_sub_ps(p1z, *(p2z++));
				dx = _mm_mul_ps(dx, dx);
				dy = _mm_mul_ps(dy, dy);
				dz = _mm_mul_ps(dz, dz);
	#endif

				// square root

	#ifdef USESQRT
				d2 = _mm_sqrt_ps(d2);
	#endif

				// get entry points in the PDF table
				if (irem<4){
	#ifndef NOUSESTEP
					d2 = _mm_mul_ps(d2,pvals);
	#endif
					jumps = _mm_cvtps_epi32(d2);
				}

	#ifdef DBLSPD
				__m128 d22 = _mm_add_ps(dx, dy);
				d22 = _mm_add_ps(d22, dz);
				d22 = _mm_sub_ps(d22,d0);
	#ifndef NOUSESTEP
				d22 = _mm_mul_ps(d22,pvals);
	#endif
				__m128i jumps2 = _mm_cvtps_epi32(d22);
	#endif
				// update PDF
				for (int i=0;i<4;i++){
					if (i>=irem)
						pResult[_mm_cvtsi128_si32(jumps)]++;
					jumps=_mm_srli_si128(jumps,4);
				}


	#ifdef DBLSPD
				for (int i=0;i<4;i++){
					if (i>=(irem-4))
						pResult[_mm_cvtsi128_si32(jumps2)]++;
					jumps2=_mm_srli_si128(jumps2,4);
				}
	#endif
			}
#endif
		}
		else ibeg=0;


// MAIN CORE: 4-tuples

//#pragma omp parallel for schedule(static, 2)

		for(int j=ibeg;j<p2optsize;j++)
		{
			// get 4 deltacoords

			__m128 dx = _mm_sub_ps(p1x, *(p2x++));
			__m128 dy = _mm_sub_ps(p1y, *(p2y++));
			__m128 dz = _mm_sub_ps(p1z, *(p2z++));

			// square

			dx = _mm_mul_ps(dx, dx);
			dy = _mm_mul_ps(dy, dy);
			dz = _mm_mul_ps(dz, dz);

			// sum dx+dy+dz

			__m128 d2 = _mm_add_ps(dx, dy);
			d2 = _mm_add_ps(d2, dz);
			//subtract the offset of the distances

			d2 = _mm_sub_ps(d2,d0);

#ifdef DBLSPD
			dx = _mm_sub_ps(p1x, *(p2x++));
			dy = _mm_sub_ps(p1y, *(p2y++));
			dz = _mm_sub_ps(p1z, *(p2z++));
			dx = _mm_mul_ps(dx, dx);
			dy = _mm_mul_ps(dy, dy);
			dz = _mm_mul_ps(dz, dz);
#endif

			// square root

#ifdef USESQRT
			d2 = _mm_sqrt_ps(d2);
#endif

			// get entry points in the PDF table
#ifndef NOUSESTEP
			d2 = _mm_mul_ps(d2,pvals);
#endif
			__m128i jumps = _mm_cvtps_epi32(d2);

#ifdef DBLSPD
			__m128 d22 = _mm_add_ps(dx, dy);
			d22 = _mm_add_ps(d22, dz);
			//subtract the offset of the distances
			d22 = _mm_sub_ps(d22,d0);
#ifndef NOUSESTEP
			d22 = _mm_mul_ps(d22,pvals);
#endif
			__m128i jumps2 = _mm_cvtps_epi32(d22);

#endif

			// update PDF

			for (int i=0;i<4;i++){
			//	int temp=_mm_cvtsi128_si32(jumps);
				pResult[_mm_cvtsi128_si32(jumps)]++;
				jumps=_mm_srli_si128(jumps,4);
			}

#ifdef DBLSPD
			for (int i=0;i<4;i++){
				pResult[_mm_cvtsi128_si32(jumps2)]++;
				jumps2=_mm_srli_si128(jumps2,4);
			}
#endif
		}

#ifdef JUSTENOUGH
// LEFTOVERS: END OF VECTOR
		if (p2optrem>0)
		{
			//flag is a control to correctly shrink the Leftover vector when in the bottom corner of the distance matrix Only necessary when p1Array=p2Array
			if (i>maxBegVect && sameArray){
				flag++;
			}

			// get 4 deltacoords
#ifndef DBLSPD
			__m128 dx = _mm_sub_ps(p1x, *(p2x));
			__m128 dy = _mm_sub_ps(p1y, *(p2y));
			__m128 dz = _mm_sub_ps(p1z, *(p2z));
#else
			__m128 dx = _mm_sub_ps(p1x, *(p2x++));
			__m128 dy = _mm_sub_ps(p1y, *(p2y++));
			__m128 dz = _mm_sub_ps(p1z, *(p2z++));
#endif

			// square
			dx = _mm_mul_ps(dx, dx);
			dy = _mm_mul_ps(dy, dy);
			dz = _mm_mul_ps(dz, dz);

			// sum dx+dy+dz
			__m128 d2 = _mm_add_ps(dx, dy);
			d2 = _mm_add_ps(d2, dz);
			//subtract the offset of the distances
			d2 = _mm_sub_ps(d2,d0);

#ifdef DBLSPD
			if(p2optrem>4){
				dx = _mm_sub_ps(p1x, *(p2x));
				dy = _mm_sub_ps(p1y, *(p2y));
				dz = _mm_sub_ps(p1z, *(p2z));
				dx = _mm_mul_ps(dx, dx);
				dy = _mm_mul_ps(dy, dy);
				dz = _mm_mul_ps(dz, dz);
			}
#endif
			// square root
#ifdef USESQRT
			d2 = _mm_sqrt_ps(d2);
#endif
			// get entry points in the PDF table
#ifndef NOUSESTEP
			d2 = _mm_mul_ps(d2,pvals);
#endif
			__m128i jumps = _mm_cvtps_epi32(d2);
#ifdef DBLSPD
			__m128 d22;
			__m128i jumps2;
			if (p2optrem>4){
				d22= _mm_add_ps(dx, dy);
				d22 = _mm_add_ps(d22, dz);
				//subtract the offset of the distances
				d22 = _mm_sub_ps(d22,d0);
#ifndef NOUSESTEP
				d22 = _mm_mul_ps(d22,pvals);
#endif
				jumps2 = _mm_cvtps_epi32(d22);
			}
#endif

			// update PDF

			//Start at flag and only increase the PDF for the flag<=k<=p2optrem
			for (int i=0;i<4;i++){
				if (i<p2optrem && i>=flag)
					pResult[_mm_cvtsi128_si32(jumps)]++;
				jumps=_mm_srli_si128(jumps,4);
			}
#ifdef DBLSPD
			for (int i=0;i<4;i++){
				if (i<p2optrem-4 && i>=flag-4)
					pResult[_mm_cvtsi128_si32(jumps2)]++;
				jumps2=_mm_srli_si128(jumps2,4);
			}

#endif
		}
#endif

	}

}

int Distance::CalcDistLatticeSym(int _shell, Position *firstPos)
{
	int i, j, k, n1, n2, place, nbins;
	Position *curPos=firstPos;
	float *x1, *y1, *z1, *x2, *y2, *z2;
	float minx=0, miny=0, minz=0, maxx=0, maxy=0, maxz=0;
	double *d, *m;
	//Calc max number of possible dist
	shell=_shell;
	for (i=1;i<shell;i++){
		curPos=curPos->nextPos;
	}
	x1=curPos->x;
	y1=curPos->y;
	z1=curPos->z;
	n1=curPos->nAtoms;
	for (i=0;i<n1;i++){
		maxx=getmax(maxx,x1[i]);
		minx=getmin(minx,x1[i]);
		maxy=getmax(maxy,y1[i]);
		miny=getmin(miny,y1[i]);
		maxz=getmax(maxz,z1[i]);
		minz=getmin(minz,z1[i]);
	}
	nbins=int(ceil(((maxx-minx)*(maxx-minx)+(maxy-miny)*(maxy-miny)+(maxz-minz)*(maxz-minz))*2)+1);

	try {
		d=new double [nbins];
		m=new double [nbins];
	}
	catch (exception& e){
		printf("\nException found in CalcDistances alloc: %s\n", e.what());
		return -1;
	}
	for(i=0;i<nbins;i++){
		d[i]=i/2.0;
		m[i]=0;
	}
	//Calc distances within shell
	x2=curPos->x;
	y2=curPos->y;
	z2=curPos->z;
	m[0]+=n1;
	for (i=0;i<n1;i++){
		for (j=i+1;j<n1;j++){
			place=int(2*((x1[i]-x2[j])*(x1[i]-x2[j])+(y1[i]-y2[j])*(y1[i]-y2[j])+(z1[i]-z2[j])*(z1[i]-z2[j])));
			if(place<nbins)
				m[place]+=2.0;
			else{
				printf("\nError: Calculated distance larger than allowed in histogram array\n");
				return -1;
			}
		}
	}
	//Calc distances btw shell and smaller shells
	curPos=firstPos;
	for (k=1;k<shell;k++){
		n2=curPos->nAtoms;
		x2=curPos->x;
		y2=curPos->y;
		z2=curPos->z;
		for (i=0;i<n1;i++){
			for (j=0;j<n2;j++){
				place=int(2*((x1[i]-x2[j])*(x1[i]-x2[j])+(y1[i]-y2[j])*(y1[i]-y2[j])+(z1[i]-z2[j])*(z1[i]-z2[j])));
				m[place]+=2.0;
			}
		}
		curPos=curPos->nextPos;
	}
	//compress distances
	for (i=0;i<nbins;i++){
		if (m[i]!=0)
			nDist++;
	}

	try{
		dist=new double [nDist];
		mult=new double [nDist];
	}
	catch (exception& e){
		printf("\nException found in CalcDistances alloc: %s\n", e.what());
		return -1;
	}
	//Set the class variables
	place=0;
	for (i=0;i<nbins;i++){
		if (m[i]!=0){
			dist[place]=sqrt(d[i]);
			mult[place]=m[i];
			place++;
		}
	}

	delete [] d;
	delete [] m;

	return 0;
}
void Distance::ConsolidateDistance()
{
	int nNZDist=0;
	double *d, *m;

	for (int i=0;i<nDist;i++){
		if (mult[i]!=0) nNZDist++;
	}

	try{
		d= new double [nNZDist];
		m= new double [nNZDist];
	}
	catch (exception& e){
		printf("\nException found in ConsolidateDistance alloc: %s\n", e.what());
		exit(-1);
	}

	nNZDist=0;
	for (int i=0;i<nDist;i++){
		if (mult[i]!=0) {
			d[nNZDist]=dist[i];
			m[nNZDist]=mult[i];
			nNZDist++;
		}
	}
	try{
		delete [] dist;
		delete [] mult;
		dist=new double [nNZDist];
		mult=new double [nNZDist];
	}
	catch (exception& e){
		printf("\nException found in ConsolidateDistance alloc: %s\n", e.what());
		exit(-1);
	}
	nDist=nNZDist;
	for (int i=0;i<nDist;i++){
		dist[i]=d[i];
		mult[i]=m[i];
	}
	delete [] d;
	delete [] m;
}
void Distance::InitDistance(int nD)
{
	nDist=nD;
	try{
		if (dist!=0)
			delete [] dist;
		dist= new double [nDist];
		if (mult!=0)
			delete [] mult;
		mult= new double [nDist];
	}
	catch (exception& e){
		printf("\nException found in InitDistance alloc: %s", e.what());
		exit(1);
	}
	//Do not need to recalculate sqrts for distances that already exist
	for (int i=0;i<nDist;i++){
		mult[i]=0;
		dist[i]=sqrt(i*SRPREC);
	}
}
void Distance::ScaleDistance(double scale)
{
	for (int i=0;i<nDist;i++){
		dist[i]*=scale;
	}
}
