//=============================================================================
// filename: random.cpp
// version: 1.01
// date: 13.8.2014
// author: Robert Hess
// description: a set of random number generators
//=============================================================================

#include "random.h"
#include <iostream>
using namespace std;

#define IM1 2147483563l
#define IM2 2147483399l
#define AM (1.0f/IM1)
#define IMM1 2147483562l   // (IM1-1)
#define IA1 40014l
#define IA2 40692l
#define IQ1 53668l
#define IQ2 52774l
#define IR1 12211
#define IR2 3791
#define NDIV 67108862l     // (1+IMM1/RAND_NTAB)
// #define EPS 1.2e-7
#define RNMX .99999988f     // (1.0-EPS)


//-------------------------------------------------------------------------
cRandU::cRandU()
//-------------------------------------------------------------------------
/** This is the default constructor for the random number generator.
It initialises the generator with seed-value 0, see method sRand(). */
//-------------------------------------------------------------------------
{
	sRand(0);
}


//-------------------------------------------------------------------------
cRandU::cRandU(
	long seed)	///< [in] value to initialize the random number generator
//-------------------------------------------------------------------------
/** This constructor initialises the random number generator with the given
seet-value, see method sRand(). */
//-------------------------------------------------------------------------
{
	sRand(seed);
}


//-------------------------------------------------------------------------
void cRandU::sRand(
	long seed)	///< [in] value to initialize the random number generator
//-------------------------------------------------------------------------
/** Initialises the random number generator to a given point. The generated
random numbers are in fact a long chain of deterministic numbers, i.e. each
time the generator is started under same conditions it will generate exactly
the same numbers. This method sRand() may be used to initialize the generator
to different points in the random number chain. The parameter seed is not
the position in the chain. E.g. seed values 0 and 1 do not represent positions
one and two in the random number chain. Hence, each value for the parameter
can be thought of a completely new position in the random number chain.

Commonly there are three ways to prepare the generator with this method:
-# Seed-value zero: This initialises the generator in order that it always
generates pecisely the same results.
-# User input seed-value: The user can perform experiments with different
seed values to perform statistics.
-# Seed-value by actual time: Foir the user this will give result in
non-predictable random numbers. */
//-------------------------------------------------------------------------
{
	int    j;
	long   k;

	idum = seed+1;
	idum2 = idum;
	for(j=RAND_NTAB+7; j>=0; j--) {
		k = idum/IQ1;
		idum = IA1*(idum-k*IQ1)-k*IR1;
		if(idum<0) idum += IM1;
		if(j<RAND_NTAB) iv[j] = idum;
	}
	iy = iv[0];
}


//-------------------------------------------------------------------------
float cRandU::rand(void)
//-------------------------------------------------------------------------
/** This method generaties a random number in the range of (including) zero
to (excluding) one. The generator may be initialised with method sRand(). */
//-------------------------------------------------------------------------
// A 'good' random number generator!
{
	static int j;
	static long k;
	static float temp;

	k = idum/IQ1;						// compute idum=(IA1*idum)%IM1
	idum = IA1*(idum-k*IQ1)-k*IR1;		//     by Schrage's method (no overflow)
	if(idum<0) idum += IM1;

	k = idum2/IQ2;						// compute idum2=(IA2*idum2)%IM2
	idum2 = IA2*(idum2-k*IQ2)-k*IR2;	//     by Schrage's method (no overflow)
	if(idum2<0) idum2 += IM2;

	j = (int)(iy/NDIV);					// use previous random number as index
	iy = iv[j]-idum2;					// the difference if idum and idum2 to
										//   generate new output
	iv[j] = idum;                        // store idum in array
	if(iy<1) iy += IMM1;				// get iy in the range from 1 to IM1-1
	if((temp=AM*iy)>RNMX) return RNMX;	// users do not expect endpoint value
	else return temp;
}


//-------------------------------------------------------------------------
bool cRandU::save(void)
//-------------------------------------------------------------------------
/** Saves the actual state of the random number generator in the file
'RandU.sav' at the default directory. Together with restore() it is possible
to interrupt long calculations and to continue them at precisely the same
state. */
//-------------------------------------------------------------------------
{
	bool error=false;
	ofstream fout;

	// open file
	fout.open("RandU.sav");
	if(!fout.is_open()) {
		error = true;
	} else {

		// write content
		error = (fout << *this).fail()!=0;

		// close file
		fout.close();
	}

	return error;
}


//-------------------------------------------------------------------------
bool cRandU::restore(void)
//-------------------------------------------------------------------------
/** Restores the actual state of the random number generator from the file
'RandU.sav' at the default directory. Together with save() it is possible
to interrupt long calculations and to continue them at precisely the same
state. */
//-------------------------------------------------------------------------
{
	ifstream finp;
	bool  error=false;

	// open file
	finp.open("RandU.sav");
	if(!finp.is_open()) {
		error = true;
	} else {

		// read content
		error = (finp >> *this).good()==0;

		// close file
		finp.close();
	}

	return error;
}

//-------------------------------------------------------------------------
ofstream &operator<<(ofstream &fout, cRandU &randU)
//-------------------------------------------------------------------------
/** Saves the actual state of the random number generator in a given stream.
Together with operator>>() it is possible to interrupt long calculations
and to continue them at precisely the same state. */
//-------------------------------------------------------------------------
{
	int   i;

	if(fout.good()) fout << randU.idum2 << "\n";
	if(fout.good()) fout << randU.iy << "\n";
	for(i=0; i<RAND_NTAB; i++)
		if(fout.good()) fout << randU.iv[i] << "\n";
	if(fout.good()) fout << randU.idum << "\n";

	return fout;
}

//-------------------------------------------------------------------------
ifstream &operator>>(ifstream &finp, cRandU &randU)
//-------------------------------------------------------------------------
/** Restores the actual state of the random number generator from in a
given stream. Together with operator<<() it is possible to interrupt long
calculations and to continue them at precisely the same state. */
//-------------------------------------------------------------------------
{
	int   i;

	if(finp.good()) finp >> randU.idum2;
	if(finp.good()) finp >> randU.iy;
	for(i=0; i<RAND_NTAB; i++)
		if(finp.good()) finp >> randU.iv[i];
	if(finp.good()) finp >> randU.idum;

	return finp;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
cRandContAR::cRandContAR(cRandU *pRandU)
//-------------------------------------------------------------------------
{
	this->pRandU = pRandU;
	reset();
}


//-------------------------------------------------------------------------
cRandContAR::~cRandContAR()
//-------------------------------------------------------------------------
{
	reset();
}


//-------------------------------------------------------------------------
void cRandContAR::reset()
//-------------------------------------------------------------------------
{
	error = false;
	prepared = false;
}


//-------------------------------------------------------------------------
bool cRandContAR::setDistrib(const std::vector<double> distribution[2])
//-------------------------------------------------------------------------
{
	// assumption: no error
	error = false;

	// copy points
	try {
		dist[0] = distribution[0];
		dist[1] = distribution[1];
	} catch(...) {
		error = true;
	}

	// reset prepared status
	prepared = false;

	return error;
}


//-------------------------------------------------------------------------
bool cRandContAR::prepare()
//-------------------------------------------------------------------------
{
	int i;				// local index variable
	double value;		// actual value
	int n;				// number of points

	// check points
	n = dist[0].size();
	error = error || n<2;		// not less than 2 points
	for(i=1; !error&&i<n; i++)	// ordered x-values
		error = error || dist[0][i-1]>dist[0][i];
	for(i=0; !error&&i<n; i++)	// no y-value below zero
		error = error || dist[1][i]<0;

	// get memory for inverted curve
	if(!error) {
		try {
			inv[0].resize(n);
			inv[1].resize(n);
		} catch(...) {
			error = true;
		}
	}

	// integrate covering curve
	if(!error) {
		inv[0][0] = 0;
		inv[1][0] = dist[0][0];
		for(i=1; !error&&i<n; i++) {
			value = dist[1][i-1]>dist[1][i]?dist[1][i-1]:dist[1][i];
			value *= dist[0][i]-dist[0][i-1];
			inv[0][i] = inv[0][i-1]+value;
			inv[1][i] = dist[0][i];
		}
		error = inv[0][n-1]<=0;
	}
	if(!error) {
		for(i=1; i<n; i++) inv[0][i] /= inv[0][n-1];
		invX = &inv[0][0];
		invY = &inv[1][0];
		maxInvIndex = n-1;
		distY = &dist[1][0];
		step0 = 1<<(int)(log(maxInvIndex+0.5)/log(2.0));
	}

	// set prepared status
	if(!error) prepared = true;

	return error;
}


//-------------------------------------------------------------------------
double cRandContAR::rand(void)
//-------------------------------------------------------------------------
{
	int i;			// index in inversed curve
	int step;		// step to find index in inversed curve
	double value=0;	// random number to be evaluated
	double tmp;
	bool finished;

	// if not prepared return
	if(!prepared) return 0;

	// assumption: not finished
	finished=false;

	// repeat until valid number found
	do {
		// get uniform random numberr
		value = pRandU->rand();

		// find index in inverted array
		i = 0;
		step = step0;
		while(step) {
			if((i|step)<maxInvIndex && invX[i|step]<value) i |= step;
			step >>= 1;
		}
//		i = maxInvIndex-1;
//		step = maxInvIndex;
//		while(step>1) {
//			step = (step+1)/2;
//			if(step>i) step = i;
//			if(invX[i-step+1]>value) i -= step;
//		}

		// if within step, finished
//		if(invX[i]==invX[i+1]) {
//			value = (invY[i]+invY[i+1])/2;
//			finished = true;
//		} else {
			// fraction between two points
			value = (value-invX[i])/(invX[i+1]-invX[i]);
			// correct density in distribution
			tmp = distY[i]+value*(distY[i+1]-distY[i]);
			// fraction underneeth covering curve
			tmp /= distY[i]>distY[i+1]?distY[i]:distY[i+1];
			// if a random number is underneeth correct value, finish
			if(tmp>=pRandU->rand()) {
				value = invY[i]+(invY[i+1]-invY[i])*value;
				finished = true;
			}
//		}
	} while(!finished);

	return value;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
cRandCont::cRandCont()
//-------------------------------------------------------------------------
{
	this->pRandU = NULL;
	reset();
}


//-------------------------------------------------------------------------
cRandCont::cRandCont(cRandU *pRandU)
//-------------------------------------------------------------------------
{
	this->pRandU = pRandU;
	reset();
}


//-------------------------------------------------------------------------
cRandCont::cRandCont(cRandCont &orig)
//-------------------------------------------------------------------------
{
	copy(orig);
}


//-------------------------------------------------------------------------
cRandCont::~cRandCont()
//-------------------------------------------------------------------------
{
	reset();
}


//-------------------------------------------------------------------------
void cRandCont::reset()
//-------------------------------------------------------------------------
{
	error = false;
	prepared = false;
}


//-------------------------------------------------------------------------
void cRandCont::copy(cRandCont &orig)
//-------------------------------------------------------------------------
{
	error = orig.error;
	pRandU = orig.pRandU;
	try {
		dist[0] = orig.dist[0];
		dist[1] = orig.dist[1];
		inv[0] = orig.inv[0];
		inv[1] = orig.inv[1];
	} catch(...) {
		error = true;
	}
	invX = (error||inv[0].size()<1)?NULL:&inv[0][0];
	invY = (error||inv[1].size()<1)?NULL:&inv[1][0];
	maxInvIndex = orig.maxInvIndex;
	prepared = orig.prepared && !error;
	step0 = orig.step0;
}

//-------------------------------------------------------------------------
cRandCont &cRandCont::operator=(cRandCont &orig)
//-------------------------------------------------------------------------
{
	copy(orig);

	return *this;
}

//-------------------------------------------------------------------------
bool cRandCont::setDistrib(const std::vector<double> distribution[2])
//-------------------------------------------------------------------------
{
	// assumption: no error
	error = false;

	// copy points
	try {
		dist[0] = distribution[0];
		dist[1] = distribution[1];
	} catch(...) {
		error = true;
	}

	// reset prepared status
	prepared = false;

	return error;
}


//-------------------------------------------------------------------------
bool cRandCont::prepare()
//-------------------------------------------------------------------------
{
	int i;				// local index variable
	double value;		// actual value
	int n;				// number of points

	// check points
	n = dist[0].size();
	error = error || n<2;		// not less than 2 points
	for(i=1; !error&&i<n; i++)	// ordered x-values
		error = error || dist[0][i-1]>dist[0][i];
	for(i=0; !error&&i<n; i++)	// no y-value below zero
		error = error || dist[1][i]<0;

	// get memory for inverted curve
	if(!error) {
		try {
			inv[0].resize(n);
			inv[1].resize(n);
		} catch(...) {
			error = true;
		}
	}

	// integrate curve
	if(!error) {
		inv[0][0] = 0;
		inv[1][0] = dist[0][0];
		for(i=1; !error&&i<n; i++) {
			value = (dist[1][i-1]+dist[1][i])*(dist[0][i]-dist[0][i-1]);
			inv[0][i] = inv[0][i-1]+value;
			inv[1][i] = dist[0][i];
		}
		error = inv[0][n-1]<=0;
	}
	if(!error) {
		for(i=1; i<n; i++) inv[0][i] /= inv[0][n-1];
		invX = &inv[0][0];
		invY = &inv[1][0];
		maxInvIndex = n-1;
		step0 = 1<<(int)(log(maxInvIndex+0.5)/log(2.0));
	}

	// set prepared status
	if(!error) prepared = true;

	return error;
}


//-------------------------------------------------------------------------
bool cRandCont::prepareRelAcc(
	double acc/*=0.01*/)		// relative accuracy, i.e. 0.01 for 1%
//-------------------------------------------------------------------------
{
	unsigned i, j;		// local index variables
	unsigned j0;		// offset index for actual interval
	double v1, v2;		// actual values
	unsigned n;			// number of points
	unsigned steps;		// number of steps within an x-interval

	// check points
	error = error || dist[0].size()<2;		// not less than 2 points
	error = error || dist[0].size()!=dist[1].size();
	for(i=1; !error&&i<dist[0].size(); i++)	// ordered x-values
		error = error || dist[0][i-1]>dist[0][i];
	for(i=0; !error&&i<dist[0].size(); i++)	// no y-value below zero
		error = error || dist[1][i]<=0;

	// check accuracy
	error = error || acc<=0;		// no negative accuracies

	// evaluate number of elements
	n = 1;
	for(i=1; !error && i<dist[0].size(); i++) {
		steps = (unsigned)(fabs(0.5/acc*log(dist[1][i]/dist[1][i-1]))+0.5);
		if(steps<1) steps = 1;
//cout << "step: " << step << endl;
		n += steps;
	}

	// get memory for inverted curve
	if(!error) {
		try {
			inv[0].resize(n);
			inv[1].resize(n);
		} catch(...) {
			error = true;
		}
	}

	// fill vector
	if(!error) {
		inv[0][0] = dist[1][0];
		inv[1][0] = dist[0][0];
		j0 = 0;
		for(i=1; i<dist[0].size(); i++) {
			steps = (unsigned)(fabs(0.5/acc*log(dist[1][i]/dist[1][i-1]))+0.5);
			if(steps<1) steps = 1;
			for(j=1; j<=steps; j++) {
				inv[0][j0+j] = dist[1][i-1]*exp(j*log(dist[1][i]/dist[1][i-1])/steps);
				inv[1][j0+j] = dist[0][i-1]+(dist[0][i]-dist[0][i-1])
					*(inv[0][j0+j]-dist[1][i-1])/(dist[1][i]-dist[1][i-1]);
			}
			j0 += steps;
		}
	}

	// integrate and normalize curve
	if(!error) {
		v2 = inv[0][0];
		inv[0][0] = 0;
		for(i=1; !error&&i<inv[0].size(); i++) {
			v1 = v2;
			v2 = inv[0][i];
			inv[0][i] = inv[0][i-1]+(v1+v2)*(inv[1][i]-inv[1][i-1]);
		}
		for(i=1; !error&&i<inv[0].size(); i++)
			inv[0][i] /= inv[0][inv[0].size()-1];
		invX = &inv[0][0];
		invY = &inv[1][0];
		maxInvIndex = n-1;
		step0 = 1<<(int)(log(maxInvIndex+0.5)/log(2.0));
	}

	// set prepared status
	if(!error) prepared = true;

	return error;
}


//-------------------------------------------------------------------------
double cRandCont::rand(void)
//-------------------------------------------------------------------------
{
	static int i;			// index in inversed curve
	static int step;		// step to find index in inversed curve
	static double value;	// random value

	if(!prepared) return 0;

	// get random number
	value = pRandU->rand();

	// find index in inversed curve
	i = 0;
	step = step0;
	while(step) {
		if((i|step)<maxInvIndex && invX[i|step]<value) i |= step;
		step >>= 1;
	}

	// interpolate
	value = (value-invX[i])/(invX[i+1]-invX[i]);
	value = invY[i]+(invY[i+1]-invY[i])*value;

	return value;
}


//-------------------------------------------------------------------------
double cRandCont::rand(float randU)
//-------------------------------------------------------------------------
{
	static int i;			// index in inversed curve
	static int step;		// step to find index in inversed curve
	static double value;	// random value

	if(!prepared) return 0;

	// get random number
	value = randU;

	// find index in inversed curve
	i = 0;
	step = step0;
	while(step) {
		if((i|step)<maxInvIndex && invX[i|step]<value) i |= step;
		step >>= 1;
	}

	// interpolate
	value = (value-invX[i])/(invX[i+1]-invX[i]);
	value = invY[i]+(invY[i+1]-invY[i])*value;

	return value;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
cRandDisc::cRandDisc(cRandU *pRandU)
//-------------------------------------------------------------------------
{
	this->pRandU = pRandU;
	reset();
}


//-------------------------------------------------------------------------
cRandDisc::~cRandDisc()
//-------------------------------------------------------------------------
{
	reset();
}


//-------------------------------------------------------------------------
void cRandDisc::reset()
//-------------------------------------------------------------------------
{
	dist[0].clear();
	dist[1].clear();
	error = false;
	prepared = false;
}


//-------------------------------------------------------------------------
bool cRandDisc::setDistrib(const std::vector<double> distribution[2])
//-------------------------------------------------------------------------
{
	// assumption: no error
	error = false;

	// copy points
	try {
		dist[0] = distribution[0];
		dist[1] = distribution[1];
	} catch(...) {
		error = true;
	}

	// reset prepared status
	prepared = false;

	return error;
}


//-------------------------------------------------------------------------
bool cRandDisc::add(double x, double y)
//-------------------------------------------------------------------------
{
	// add point to distribution
	try {
		dist[0].push_back(x);
		dist[1].push_back(y);
	} catch(...) {
		error = true;
	}
	
	// return error flag
	return error;
}


//-------------------------------------------------------------------------
bool cRandDisc::prepare()
//-------------------------------------------------------------------------
{
	unsigned i;			// local index variable
	unsigned n=0;		// number of points

	// check points
	if(!error) n = dist[0].size();
	error = error || n<2;
	for(i=1; !error&&i<n; i++)	// ordered x-values
		error = error || dist[0][i-1]>dist[0][i];
	for(i=0; !error&&i<n; i++)	// no y-value below zero
		error = error || dist[1][i]<0;
	// remove points with zero value
	for(i=0; i<dist[0].size(); i++) {
		if(dist[1][i]<=0) {
			dist[0].erase(dist[0].begin()+i);
			dist[1].erase(dist[1].begin()+i);
		}
	}
	// combine points of same position
	for(i=1; i<dist[0].size(); i++) {
		if(dist[0][i-1]==dist[0][i]) {
			dist[1][i-1] += dist[1][i];
			dist[0].erase(dist[0].begin()+i);
			dist[1].erase(dist[1].begin()+i);
		}
	}
	if(!error) n = dist[0].size();

	// get memory for inverted curve
	if(!error) {
		try {
			inv[0].resize(n);
			inv[1].resize(n);
		} catch(...) {
			error = true;
		}
	}

	// integrate curve
	if(!error) {
		inv[0][0] = dist[1][0];
		inv[1][0] = dist[0][0];
		for(i=1; !error&&i<n; i++) {
			inv[0][i] = inv[0][i-1]+dist[1][i];
			inv[1][i] = dist[0][i];
		}
		error = inv[0][n-1] <= 0;
	}
	if(!error) {
		for(i=0; i<n; i++) inv[0][i] /= inv[0][n-1];
		invX = &inv[0][0];
		invY = &inv[1][0];
		maxInvIndex = n-1;
	}

	// set prepared status
	if(!error) prepared = true;

	// return error flag
	return error;
}

//-------------------------------------------------------------------------
double cRandDisc::rand(void)
//-------------------------------------------------------------------------
{
	static int i;			// index in inversed curve
	static int step;		// step to find index in inversed curve
	static double value;

	if(!prepared) return 0;

	step = i = maxInvIndex;
	value = pRandU->rand();
	while(step>1) {
		step = (step+1)/2;
		if(step>i) step = i;
		if(invX[i-step]>value) i -= step;
	}

	return invY[i];
}
/*
 N i0 step
 2  1  1
 3  2  2 1
 4  3  2 1
 5  4  3 2 1
 6  5  3 2 1
 7  6  4 2 1
 8  7  4 2 1
 9  8  5 3 2 1
10  9  5 3 2 1
11 10  6 3 2 1
12 11  6 3 2 1
13 12  7 4 2 1
14 13  7 4 2 1
15 14  8 4 2 1
16 15  8 4 2 1
17 16  9 5 3 2 1
18 17  9 5 3 2 1
19 18 10 5 3 2 1
20 19 10 5 3 2 1
*/
