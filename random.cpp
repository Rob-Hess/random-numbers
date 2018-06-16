//=============================================================================
// filename: random.cpp
// version: 1.06
// date: June 16th 2018
// author: Robert Hess, HAW-Hamburg, Germany
// description: a set of random number generators
//=============================================================================

#include "random.h"
#include <iostream>
#include <cmath>
#include <exception>
using namespace std;

// error messages
#define RAND_BAD_PROB     runtime_error("Invalid probability data.")
#define RAND_NO_RAND_BASE runtime_error("No random generator base found.")


//=============================================================================
//=============================================================================
//=============================================================================


//-----------------------------------------------------------------------------
cRandKiss::cRandKiss()
//-----------------------------------------------------------------------------
/** This is the default constructor for the random number generator.
It initialises the generator with seed-value 0, see method srand(). */
//-----------------------------------------------------------------------------
{
	// set default seed value
	srand(0);
}


//-----------------------------------------------------------------------------
cRandKiss::cRandKiss(
	const cRandKiss &source)	///< [in] generator to be copied
//-----------------------------------------------------------------------------
/** Copy constructor for the random number generator.
It takes an exact copy of the generator without calling srand(). */
//-----------------------------------------------------------------------------
{
	// copy all attributes of this class
	x32 = source.x32;
	y32 = source.y32;
	z32 = source.z32;
	c32 = source.c32;
}


//-----------------------------------------------------------------------------
cRandKiss::cRandKiss(
	uint32_t seed)	///< [in] value to initialise the random number generator
//-----------------------------------------------------------------------------
/** This constructor initialises the random number generator with the given
seed-value, see method srand(). */
//-----------------------------------------------------------------------------
{
	// set seed value for this generator
	srand(seed);
}


//-----------------------------------------------------------------------------
void cRandKiss::srand(
	uint32_t seed)	///< [in] value to initialise the random number generator
//-----------------------------------------------------------------------------
/** Initialises the random number generator to a given point. The generated
random numbers are in fact a long chain of deterministic numbers, i.e. each
time the generator is started under equal conditions it will generate exactly
the same numbers. This method srand() may be used to initialise the generator
to different points in the random number chain. The parameter seed is not
the position in the chain. E.g. seed values 0 and 1 do not represent positions
one and two in the random number chain. Hence, each value for the parameter
can be thought of a completely new position in the random number chain.

Commonly there are three ways to prepare the generator with this method:
-# Zero seed-value: This initialises the generator in order that it always
generates precisely the same results.
-# User input seed-value: The user can perform experiments with different
seed values to perform statistics.
-# Seed-value by actual time: From the users perspective this will result in
non-predictable random numbers. */
//-----------------------------------------------------------------------------
{
	int i;	// local index variable

	// set seed
	x32 = seed;
	y32 = seed + 1;
	if (!y32) y32 = 1;	// for xorshift it must not be zero
	z32 = seed + 2;
	c32 = seed + 3;

	// warm up generator
	for (i = 0; i < 10; i++) rand();
}


//-----------------------------------------------------------------------------
uint32_t cRandKiss::rand(void)
//-----------------------------------------------------------------------------
/** This method generates a uniform 32 bit unsigned random number. The
generator may be initialised with method srand(). */
//-----------------------------------------------------------------------------
{
	static uint64_t t;	// temporary value

	// linear congruential generator
	x32 = 69069 * x32 + 12345;

	// xorshift generator
	y32 ^= y32 << 13;
	y32 ^= y32 >> 17;
	y32 ^= y32 << 5;

	// multiply-with-carry generator
	t = 698769069ULL * z32 + c32;
	c32 = t >> 32;
	z32 = (uint32_t)t;

	// return sum of three generators
	return x32 + y32 + z32;
}


//-----------------------------------------------------------------------------
double cRandKiss::randU(void)
//-----------------------------------------------------------------------------
/** This method generates a uniform random number on the interval [0,1).
The generator may be initialised with method srand(). */
//-----------------------------------------------------------------------------
{
	// evaluate and return uniform random value on interval [0,1)
	return rand() / 4294967296.0;
}


//-----------------------------------------------------------------------------
void cRandKiss::getState(
	std::vector<uint8_t> &state)	///< [out] state of random number generator
//-----------------------------------------------------------------------------
/** Stores the four state variables of type uint32_t in the vector 'state'
with 16 elements. */
//-----------------------------------------------------------------------------
{
	int i, j;				// local index variable
	uint32_t *data[] = {	// pointer array on state attributes
		&x32,
		&y32,
		&z32,
		&c32
	};

	state.resize(4 * sizeof(uint32_t));

	// loop over all attributes
	for (i = 0; i < 4; i++) {
		// loop over the four bytes of uint32_t
		for (j = 0; j < 4; j++)
			// copy actual byte
			state[4 * i + j] = (uint8_t)(*data[i] >> 8 * j);
	}
}

//-----------------------------------------------------------------------------
void cRandKiss::setState(
	const std::vector<uint8_t> &state)	///< [out] state of random number generator
//-----------------------------------------------------------------------------
/** Restores the four state variables of type uint32_t from the vector 'state'
with at least 16 elements. */
//-----------------------------------------------------------------------------
{
	int i, j;				// local index variable
	uint32_t *data[] = {	// pointer array on state attributes
		&x32,
		&y32,
		&z32,
		&c32
	};

	// loop over all attributes
	for (i = 0; state.size() >= 16 && i < 4; i++) {
		// loop over four bytes of uint32_t
		for (j = 0; j < 4; j++)
			// copy actual byte
			((uint8_t*)(data[i]))[j] = state[4 * i + j];
	}
}


//=============================================================================
//=============================================================================
//=============================================================================


//-----------------------------------------------------------------------------
cRandDistribution::cRandDistribution()
//-----------------------------------------------------------------------------
{
	// reset reference to base random generator
	pRand = NULL;

	// bring class into default state
	reset();
}


//-----------------------------------------------------------------------------
cRandDistribution::cRandDistribution(
	const cRandDistribution &source)	///< [in] distribution to be copied
//-----------------------------------------------------------------------------
{
	// copy all attributes of class
	copy(source);
}


//-----------------------------------------------------------------------------
cRandDistribution::cRandDistribution(
	cRandBase &rand)	///< [in] reference to base random generator
//-----------------------------------------------------------------------------
{
	// keep reference to base class
	pRand = &rand;

	// bring class into default state
	reset();
}


//-----------------------------------------------------------------------------
void cRandDistribution::setRandomBase(
	cRandBase &rand)	///< [in] reference to base random generator
//-----------------------------------------------------------------------------
{
	// keep reference to base class
	pRand = &rand;
}


//-----------------------------------------------------------------------------
void  cRandDistribution::srand(
	uint32_t seed)	///< [in] random seed value
//-----------------------------------------------------------------------------
{
	// set seed of base class
	if (pRand) pRand->srand(seed);
}


//-----------------------------------------------------------------------------
void cRandDistribution::reset()
//-----------------------------------------------------------------------------
{
	// bring class into default state
	prepared = false;
	cumulative.clear();
	dist[0].clear();
	dist[1].clear();
}


//-----------------------------------------------------------------------------
cRandDistribution &cRandDistribution::operator=(
	cRandDistribution &source)	///< distribution to be copied
//-----------------------------------------------------------------------------
{
	// copy all attributes of class
	copy(source);

	// return this class
	return *this;
}


//-----------------------------------------------------------------------------
void cRandDistribution::addDistPoint(
	double x,		///< [in] value of distribution point
	double prob)	///< [in] probability of distribution point
//-----------------------------------------------------------------------------
{
	// reset prepared status
	prepared = false;

	// store point
	dist[0].push_back(x);
	dist[1].push_back(prob);
}


//-----------------------------------------------------------------------------
void cRandDistribution::setDistribution(
	const std::vector<double> &x,	///< [in] values for random generator
	const std::vector<double> &p)	///< [in] probabilities for random generator
//-----------------------------------------------------------------------------
{
	// reset prepared status
	prepared = false;

	// copy points
	dist[0] = x;
	dist[1] = p;
}


//-----------------------------------------------------------------------------
void cRandDistribution::setDistribution(
	const std::vector<double> distribution[2])	///< [in] probability function
//-----------------------------------------------------------------------------
{
	// reset prepared status
	prepared = false;

	// copy points
	dist[0] = distribution[0];
	dist[1] = distribution[1];
}


//-----------------------------------------------------------------------------
void cRandDistribution::prepare(
	bool discrete)	///< [in] flag for discrete (true) or continuous (false) preparation
//-----------------------------------------------------------------------------
{
	unsigned n;			// number of points
	unsigned i;			// local index variable
	vector<double> cum;	// cumulative distribution function

	// check random generator base
	if (!pRand) throw RAND_NO_RAND_BASE;

	// check points
	n = dist[0].size();
	// equal number of points
	if (n != dist[1].size()) throw RAND_BAD_PROB;
	// not less than 2 points
	if (n < 2) throw RAND_BAD_PROB;
	// ordered x-values
	for (i = 1; !discrete && i < n; i++)
		if (dist[0][i - 1] > dist[0][i]) throw RAND_BAD_PROB;
	// no probability below zero
	for (i = 0; i < n; i++)
		if (dist[1][i] < 0) throw RAND_BAD_PROB;

	// get new memory
	cumulative.resize(dist[0].size());
	cum.resize(dist[0].size());

	// integrate curve
	if (discrete) {
		cum[0] = dist[1][0];
		for (i = 1; i<n; i++) {
			cum[i] = cum[i - 1] + dist[1][i];
		}
	}
	else {
		cum[0] = 0;
		for (i = 1; i < n; i++) {
			cum[i] = cum[i-1] + (dist[1][i-1]+dist[1][i])/2 * (dist[0][i]-dist[0][i-1]);
		}
	}
	if (cum[n - 1] <= 0) throw RAND_BAD_PROB;

	// normalise probabilities
	for (i = 0; i<n; i++) cum[i] /= cum[n - 1];

	// keep 2^32 multiple
	for (i = 0; i<n; i++) cumulative[i] = (uint32_t)(4294967295.9 * cum[i]);
	mapInp = &cumulative[0];
	mapOut = &dist[0][0];

	// maximum index and number highest bit in maximum index
	maxInvIndex = n - 1;
	for (step0 = 1 << 31; !(step0 & maxInvIndex); step0 >>= 1);
}


//-----------------------------------------------------------------------------
void cRandDistribution::copy(
	const cRandDistribution &source)	///< [in] distribution to be copied
//-----------------------------------------------------------------------------
{
	// mark as not prepared
	prepared = false;

	// copy all attributes of class
	cumulative = source.cumulative;
	pRand = source.pRand;
	dist[0] = source.dist[0];
	dist[1] = source.dist[1];
	maxInvIndex = source.maxInvIndex;
	step0 = source.step0;
	mapInp = source.mapInp;
	mapOut = source.mapOut;
	prepared = source.prepared;
}


//=============================================================================
//=============================================================================
//=============================================================================


//-----------------------------------------------------------------------------
cRandDisc::cRandDisc()
//-----------------------------------------------------------------------------
{
	// bring class into default state
	reset();
}


//-----------------------------------------------------------------------------
cRandDisc::cRandDisc(
	const cRandDisc &source)	///< [in] generator to be copied
//-----------------------------------------------------------------------------
{
	// copy all attributes of class
	copy(source);
}


//-----------------------------------------------------------------------------
cRandDisc::cRandDisc(
	cRandBase &rand)	///< [in] reference to base random generator
//-----------------------------------------------------------------------------
{
	// keep reference to base random generator
	pRand = &rand;

	// bring class into default state
	reset();
}


//-----------------------------------------------------------------------------
void cRandDisc::reset()
//-----------------------------------------------------------------------------
{
	// bring base class into default state
	cRandDistribution::reset();
}


//-----------------------------------------------------------------------------
cRandDisc &cRandDisc::operator=(
	cRandDisc &source)	///< [in] generator to be copied
//-----------------------------------------------------------------------------
{
	// copy all attributes of base class
	copy(source);

	// return class
	return *this;
}


//-----------------------------------------------------------------------------
void cRandDisc::prepare()
//-----------------------------------------------------------------------------
{
	// prepare base class
	cRandDistribution::prepare(true);

	// mark as prepared
	prepared = true;
}


//-----------------------------------------------------------------------------
double cRandDisc::rand(void)
//-----------------------------------------------------------------------------
{
	unsigned i;		// index in inverse curve
	unsigned step;	// step to find index in inverse curve
	uint32_t value;	// current uniform random value

	// if not prepared return zero
	if (!prepared) return 0.0;

	// find index in inverse curve
	value = pRand->rand();
	i = 0;
	step = step0;
	while (step) {
		if ((i | step)<=maxInvIndex && cRandDistribution::mapInp[i | step]<value) i |= step;
		step >>= 1;
	}
	if (cRandDistribution::mapInp[0] < value) i++;

	// return discrete random value
	return cRandDistribution::mapOut[i];
}


//=============================================================================
//=============================================================================
//=============================================================================


//-----------------------------------------------------------------------------
cRandDiscParam::cRandDiscParam()
//-----------------------------------------------------------------------------
{
	// reset reference to base random generator
	pRand = NULL;

	// bring class into default state
	reset();
}


//-----------------------------------------------------------------------------
cRandDiscParam::cRandDiscParam(
	const cRandDiscParam &source)	///< [in] generator to be copied
//-----------------------------------------------------------------------------
{
	// copy all attributes of class
	copy(source);
}


//-----------------------------------------------------------------------------
cRandDiscParam::cRandDiscParam(
	cRandBase &rand)	///< [in] reference to base random generator
//-----------------------------------------------------------------------------
{
	// keep reference to base random generator
	pRand = &rand;

	// bring class into default state
	reset();
}


//-----------------------------------------------------------------------------
void cRandDiscParam::reset()
//-----------------------------------------------------------------------------
{
	// bring class into default state
	prepared = false;
	value.clear();
	prob.clear();
	parameter.clear();
	cum.clear();
}


//-----------------------------------------------------------------------------
void cRandDiscParam::setRandomBase(
	cRandBase &rand)	///< [in] reference to base random generator
//-----------------------------------------------------------------------------
{
	// keep reference to base random generator
	pRand = &rand;
}


//-----------------------------------------------------------------------------
void cRandDiscParam::setRandomValues(
	const std::vector<double> &x)	///< [in] list of all values for generator
//-----------------------------------------------------------------------------
{
	// mark as not prepared
	prepared = false;

	// store random values
	value = x;
}


//-----------------------------------------------------------------------------
void cRandDiscParam::addProbabilities(
	const std::vector<double> &p,	///< [in] list of probabilities for given parameter
	double param)					///< [in] associated parameter to probabilities
//-----------------------------------------------------------------------------
{
	// mark as not prepared
	prepared = false;

	// store parameter and probabilities
	parameter.push_back(param);
	prob.push_back(p);
}


//-----------------------------------------------------------------------------
void cRandDiscParam::prepare()
//-----------------------------------------------------------------------------
{
	unsigned iParam;	// index for parameter
	unsigned iProb;		// index for probability

	// assumption: bad preparation
	prepared = false;

	// check random generator base
	if (!pRand) throw RAND_NO_RAND_BASE;

	// check number of parameters
	if (parameter.size() < 2) throw RAND_BAD_PROB;
	if (parameter.size() != prob.size()) throw RAND_BAD_PROB;

	// check order of parameters
	for (iParam = 1; iParam < parameter.size(); iParam++)
		if (parameter[iParam] < parameter[iParam - 1])
			throw RAND_BAD_PROB;

	// check number of probabilities
	if (value.size() < 2) throw RAND_BAD_PROB;
	for (iParam = 0; iParam < parameter.size(); iParam++)
		if (prob[iParam].size() != value.size())
			throw RAND_BAD_PROB;

	// integrate and normalise probabilities
	cum.resize(prob.size());
	for (iParam = 0; iParam < parameter.size(); iParam++) {
		cum[iParam].resize(prob[iParam].size());
		cum[iParam][0] = 0;
		for (iProb = 1; iProb < prob[iParam].size(); iProb++)
			cum[iParam][iProb] = cum[iParam][iProb - 1] + prob[iParam][iProb - 1];
		for (iProb = 0; iProb < prob[iParam].size(); iProb++)
			cum[iParam][iProb] /= cum[iParam].back() + prob[iParam].back();
	}

	// some further settings
	maxIndexParam = parameter.size() - 1;
	for (stepParam0 = 1 << 31; (maxIndexParam&stepParam0) == 0; stepParam0 >>= 1);
	maxIndexProb = value.size() - 1;
	for (stepProb0 = 1 << 31; (maxIndexProb&stepProb0) == 0; stepProb0 >>= 1);

	// set prepared
	prepared = true;
}


//-----------------------------------------------------------------------------
void  cRandDiscParam::srand(
	uint32_t seed)	///< [in] seed value for base random generator
//-----------------------------------------------------------------------------
{
	// set seed value for base random generator
	if (pRand) pRand->srand(seed);
}


//-----------------------------------------------------------------------------
double cRandDiscParam::rand(
	double param)	///< [in] parameter to generate random number for
//-----------------------------------------------------------------------------
{
	// evaluate and return random value
	return rand(pRand->randU(), param);
}


//-----------------------------------------------------------------------------
double cRandDiscParam::rand(
	double randU,	///< [in] uniform distributed random number on interval [0,1)
	double param)	///< [in] parameter to generate random number for
//-----------------------------------------------------------------------------
{
	unsigned iParam;	// index for parameter
	unsigned iProb;		// index for probability
	double tmp;			// temporary value
	unsigned step;		// step size to find index
	double frac;		// fraction between two parameters

	// check if prepared
	if (!prepared) return 0.0;

	// if parameter below range
	if (param <= parameter[0]) {
		iParam = 0;
		param = parameter[0];
		frac = 0;
	}

	// if parameter above range
	else if (param >= parameter.back()) {
		iParam = parameter.size() - 2;
		param = parameter.back();
		frac = 1;
	}
	// find index in between
	else {
		iParam = 0;
		for (step = stepParam0; step; step >>= 1)
			if ((iParam | step) < maxIndexParam && param >= parameter[iParam | step])
				iParam |= step;
		frac = (param - parameter[iParam]) / (parameter[iParam + 1] - parameter[iParam]);
	}

	// find index in probabilities
	iProb = 0;
	for (step = stepProb0; step; step >>= 1) {
		if ((iProb | step) <= maxIndexProb) {
			iProb ^= step;
			tmp = cum[iParam][iProb] + frac * (cum[iParam + 1][iProb] - cum[iParam][iProb]);
			if (randU < tmp) iProb ^= step;
		}
	}

	// return random value
	return value[iProb];
}


//-----------------------------------------------------------------------------
void cRandDiscParam::copy(
	const cRandDiscParam &source)	///< [in] generator to be copied
//-----------------------------------------------------------------------------
{
	// assumption: not prepared
	prepared = false;

	// copy all items
	pRand = source.pRand;
	value = source.value;
	prob = source.prob;
	parameter = source.parameter;
	cum = source.cum;
	maxIndexParam = source.maxIndexParam;
	stepParam0 = source.stepParam0;
	maxIndexProb = source.maxIndexProb;
	stepProb0 = source.stepProb0;
	prepared = source.prepared;
}


//=============================================================================
//=============================================================================
//=============================================================================


//-----------------------------------------------------------------------------
cRandCont::cRandCont()
//-----------------------------------------------------------------------------
{
	// bring class into default state
	reset();
}


//-----------------------------------------------------------------------------
cRandCont::cRandCont(
	const cRandCont &source)	///< [in] generator to be copied
//-----------------------------------------------------------------------------
{
	// copy all attributes of class
	copy(source);
}


//-----------------------------------------------------------------------------
cRandCont::cRandCont(
	cRandBase &rand)	///< [in] reference to base random generator
//-----------------------------------------------------------------------------
{
	// keep reference to base random generator
	pRand = &rand;

	// bring class into default state
	reset();
}


//-----------------------------------------------------------------------------
cRandCont &cRandCont::operator=(
	const cRandCont &source)	///< [in] generator to be copied
//-----------------------------------------------------------------------------
{
	// copy all attributes of class
	copy(source);

	// return this class
	return *this;
}


//-----------------------------------------------------------------------------
void cRandCont::reset()
//-----------------------------------------------------------------------------
{
	// bring base class into default state
	cRandDistribution::reset();

	// bring this class into default state
	a1.clear();
	a2.clear();
	a3.clear();
	a4.clear();
	linear.clear();
}


//-----------------------------------------------------------------------------
void cRandCont::prepare()
//-----------------------------------------------------------------------------
{
	std::vector<double> tmp;	// temporary array to prepare data
	unsigned i;					// local index variable
	double m, b;				// coefficients for an equation of a line
	unsigned size;				// number of points in distribution

	// assumption: bad preparation
	prepared = false;

	// prepare distribution
	cRandDistribution::prepare(false);

	// get memory
	size = dist[0].size();
	tmp.resize(size);
	a1.resize(size - 1);
	a2.resize(size - 1);
	a3.resize(size - 1);
	a4.resize(size - 1);
	linear.resize(size);

	// integrate
	tmp[0] = 0;
	for (i = 1; i < size; i++)
		tmp[i] = tmp[i - 1] + (dist[1][i - 1] + dist[1][i]) / 2 * (dist[0][i] - dist[0][i - 1]);

	// prepare vector for linear interpolation
	for (i = 0; i < size - 1; i++) {
		m = (dist[1][i + 1] - dist[1][i]) / (dist[0][i + 1] - dist[0][i]) / tmp[size - 1];
		linear[i] = abs(m) < 1e-8;
		if (linear[i]) {
			a1[i] = pow(2, -32) / dist[1][i] * tmp[size - 1];
			a2[i] = dist[0][i] - mapInp[i] * pow(2, -32) / dist[1][i] * tmp[size - 1];
		}
		else {
			b = dist[1][i] / tmp[size - 1] - m * dist[0][i];
			a1[i] = -b / m;
			a2[i] = 1 / m;
			a3[i] = 2 * m * pow(2, -32);
			a4[i] = pow(dist[1][i] / tmp[size - 1], 2) - 2 * m * mapInp[i] * pow(2, -32);
		}
	}

	// mark generator as prepared
	prepared = true;
}


//-----------------------------------------------------------------------------
double cRandCont::rand()
//-----------------------------------------------------------------------------
{
	// if no data prepared return zero
	if (!prepared) return 0.0;

	// return random number
	return rand(pRand->rand());
}


//-----------------------------------------------------------------------------
double cRandCont::rand(
	const uint32_t &input)	///< [in] uniform distributed random number on interval [0,2^32)
//-----------------------------------------------------------------------------
{
	static int mask;		// mask to set and reset bits
	register unsigned k;	// index of value to be returned
	double value;			// random value to be evaluated

	// if no data prepared return zero
	if (!prepared) return 0.0;

	// find k
	k = 0;
	for (mask = step0; mask; mask >>= 1) {
		k ^= mask;	// set bit
		if (k >= maxInvIndex || input < mapInp[k])
			k ^= mask;	// reset bit
	}

	// apply interpolation
	if (linear[k]) value = a1[k] * input + a2[k];
	else value = a1[k] + a2[k] * sqrt(a3[k] * input + a4[k]);

	// return desired value
	return value;
}


//-----------------------------------------------------------------------------
void cRandCont::copy(
	const cRandCont &source)	///< [in] generator to be copied
//-----------------------------------------------------------------------------
{
	// reset generator
	pRand = NULL;
	reset();

	// copy pointer to base generator
	pRand = source.pRand;

	// copy distribution
	dist[0] = source.dist[0];
	dist[1] = source.dist[1];

	// copy interpolation coefficients
	a1 = source.a1;
	a2 = source.a2;
	a3 = source.a3;
	a4 = source.a4;
	linear = source.linear;
}


//=============================================================================
//=============================================================================
//=============================================================================


//-----------------------------------------------------------------------------
cRandContParam::cRandContParam()
//-----------------------------------------------------------------------------
{
	// reset reference to base random generator
	pRand = NULL;

	// bring class into default state
	reset();
}


//-----------------------------------------------------------------------------
cRandContParam::cRandContParam(
	const cRandContParam &source)	///< [in] generator to be copied
//-----------------------------------------------------------------------------
{
	// copy all attributes of class
	copy(source);
}


//-----------------------------------------------------------------------------
cRandContParam::cRandContParam(
	cRandBase *pRandBase)	///< [in] reference to base random generator
//-----------------------------------------------------------------------------
{
	// store reference to base random generator
	pRand = pRandBase;

	// bring class into default state
	reset();
}


//-----------------------------------------------------------------------------
void cRandContParam::reset()
//-----------------------------------------------------------------------------
{
	// bring class into default state
	prepared = false;
	parameter.clear();
	generator.clear();
}


//-----------------------------------------------------------------------------
void cRandContParam::setRandomBase(
	cRandBase &rand)	///< [in] reference to base random generator
//-----------------------------------------------------------------------------
{
	// store reference of base random generator
	pRand = &rand;
}


//-----------------------------------------------------------------------------
void cRandContParam::setDistribution(
	const std::vector<double> &x,	///< [in] list of random values
	const std::vector<double> &p,	///< [in] list of probabilities for given values
	double param)					///< [in] parameter for given probability density function
//-----------------------------------------------------------------------------
{
	// reset prepared state
	prepared = false;

	// store parameter and probability distribution
	parameter.push_back(param);
	generator.push_back(cRandCont());
	generator.back().setDistribution(x, p);
}


//-----------------------------------------------------------------------------
void cRandContParam::setDistribution(
	const std::vector<double> distribution[2],	///< [in] probability density function
	double param)								///< [in] parameter for given probability density function
//-----------------------------------------------------------------------------
{
	// reset prepared state
	prepared = false;

	// store parameter and probability distribution
	parameter.push_back(param);
	generator.push_back(cRandCont());
	generator.back().setDistribution(distribution);
}


//-----------------------------------------------------------------------------
void cRandContParam::add(
	double x,		///< [in] value for random distribution
	double p,		///< [in] probability for random distribution
	double param)	///< [in] parameter for given distribution
//-----------------------------------------------------------------------------
{
	unsigned iParam;	// index for parameter

	// reset prepared state
	prepared = false;

	// find index for parameter
	for (iParam = 0; iParam < parameter.size() && parameter[iParam] != param; iParam++);

	// if new parameter extend arrays
	if (iParam >= parameter.size()) {
		parameter.push_back(param);
		generator.push_back(cRandCont());
	}

	// store new point
	generator[iParam].addDistPoint(x, p);
}


//-----------------------------------------------------------------------------
void cRandContParam::prepare()
//-----------------------------------------------------------------------------
{
	unsigned iParam;	// index for parameter

	// reset prepared state
	prepared = false;

	// check random generator base
	if (!pRand) throw RAND_NO_RAND_BASE;

	// check number of parameters
	if (parameter.size() < 2) throw RAND_BAD_PROB;
	if (parameter.size() != generator.size()) throw RAND_BAD_PROB;

	// prepare all generators
	for (iParam = 0; iParam < parameter.size(); iParam++) {
		generator[iParam].setRandomBase(*pRand);
		generator[iParam].prepare();
	}

	// maximum index and initial step size to find parameter index
	maxInd = parameter.size() - 2;
	for (step0 = 1 << 31; step0 && !(maxInd&step0); step0 >>= 1);

	// set prepared
	prepared = true;
}


//-----------------------------------------------------------------------------
void cRandContParam::srand(
	uint32_t seed)	///< [in] seed value base random generator
//-----------------------------------------------------------------------------
{
	// set seed value for base random generator
	if (pRand) pRand->srand(seed);
}


//-----------------------------------------------------------------------------
double cRandContParam::rand(
	double param)	///< [in] parameter to generate random value for
//-----------------------------------------------------------------------------
{
	unsigned iParam = 0;	// index for parameter
	uint32_t u;				// uniform uint32_t random number
	unsigned step;			// actual step size to find parameter index
	double r1, r2;			// points for linear interpolation

	// check if prepared
	if (!prepared) return 0.0;

	// get uniform random number
	u = pRand->rand();

	// if parameter below range
	if (param <= parameter[0]) return generator[0].rand(u);

	// if parameter above range
	if (param >= parameter.back()) return generator.back().rand(u);

	// find index of parameter
	for (step = step0; step; step >>= 1)
		if ((iParam | step) <= maxInd && param>parameter[iParam | step]) iParam |= step;

	// evaluate random number
	r1 = generator[iParam].rand(u);
	r2 = generator[iParam + 1].rand(u);

	// perform interpolation and return result
	return r1 + (r2 - r1) * (param - parameter[iParam])
		/ (parameter[iParam + 1] - parameter[iParam]);
}


//-----------------------------------------------------------------------------
void cRandContParam::copy(
	const cRandContParam &source)	///< [in] generator to be copied
//-----------------------------------------------------------------------------
{
	// copy all attributes of class
	prepared = false;
	pRand = source.pRand;
	generator = source.generator;
	parameter = source.parameter;
	prepared = source.prepared;
}


//=============================================================================
//=============================================================================
//=============================================================================


//-----------------------------------------------------------------------------
cRandExponential::cRandExponential()
//-----------------------------------------------------------------------------
{
	// reset reference to base random generator
	pRand = NULL;

	// bring class into default state
	reset();
}


//-----------------------------------------------------------------------------
cRandExponential::cRandExponential(
	const cRandExponential &source)	///< [in] generator to be copied
//-----------------------------------------------------------------------------
{
	// copy all attributes of class
	copy(source);
}


//-----------------------------------------------------------------------------
cRandExponential::cRandExponential(
	cRandBase &rand)	///< [in] reference to base random generator
//-----------------------------------------------------------------------------
{
	// keep reference to base random generator
	pRand = &rand;

	// bring class into default state
	reset();
}


//-----------------------------------------------------------------------------
void cRandExponential::reset()
//-----------------------------------------------------------------------------
{
	// bring class into default state
	nmean = -1;
}


//-----------------------------------------------------------------------------
cRandExponential &cRandExponential::operator=(
	cRandExponential &source)	///< [in] generator to be copied
//-----------------------------------------------------------------------------
{
	// copy all attributes of class
	copy(source);

	return *this;
}


//-----------------------------------------------------------------------------
void cRandExponential::setRandomBase(
	cRandBase &rand)	///< [in] reference to base random generator
//-----------------------------------------------------------------------------
{
	// keep reference of base random generator
	pRand = &rand;
}


//-----------------------------------------------------------------------------
void cRandExponential::setExpectation(
	double expectation)	///< [in] expectation of exponential distribution
//-----------------------------------------------------------------------------
{
	// store negative expectation
	nmean = -expectation;
}


//-----------------------------------------------------------------------------
double cRandExponential::rand()
//-----------------------------------------------------------------------------
{
	// evaluate and return random value
	return pRand ? nmean*log(1 - pRand->randU()) : 0.0;
}


//-----------------------------------------------------------------------------
void cRandExponential::copy(
	const cRandExponential &source)	///< [in] generator to be copied
//-----------------------------------------------------------------------------
{
	// copy attributes of this class
	pRand = source.pRand;
	nmean = source.nmean;
}


//=============================================================================
//=============================================================================
//=============================================================================


//-----------------------------------------------------------------------------
cRandNormal::cRandNormal()
//-----------------------------------------------------------------------------
{
	// bring class into default state
	reset();
}

//-----------------------------------------------------------------------------
cRandNormal::cRandNormal(
	const cRandNormal &source)	///< [in] generator to be copied
//-----------------------------------------------------------------------------
{
	// copy all items of current class
	copy(source);
}


//-----------------------------------------------------------------------------
cRandNormal::cRandNormal(
	cRandBase &rand)	///< [in] reference to base random generator
//-----------------------------------------------------------------------------
{
	// bring class into default state
	reset();

	// keep base generator in base class
	setRandomBase(rand);
}


//-----------------------------------------------------------------------------
void cRandNormal::reset()
//-----------------------------------------------------------------------------
{
	// bring class into default state
	cRandCont::reset();
	intervals = 100;
	mean = 0;
	stdev = 1.0;
	nStdev = 10;
}


//-----------------------------------------------------------------------------
cRandNormal &cRandNormal::operator=(
	cRandNormal &source)	///< [in] generator to be copied
//-----------------------------------------------------------------------------
{
	// copy class
	copy(source);

	// return current class
	return *this;
}


//-----------------------------------------------------------------------------
void cRandNormal::setRandomBase(
	cRandBase &rand)	///< [in] reference to base random generator
//-----------------------------------------------------------------------------
{
	// keep reference of base generator in base class
	cRandDistribution::setRandomBase(rand);
}


//-----------------------------------------------------------------------------
void cRandNormal::setExpectation(
	double expectation)	///< [in] expectation of normal random distribution
//-----------------------------------------------------------------------------
{
	// reset prepared status
	prepared = false;

	// store expectation
	mean = expectation;
}


//-----------------------------------------------------------------------------
void cRandNormal::setStandardDeviation(
	double sd)	///< [in] standard deviation of normal random distribution
//-----------------------------------------------------------------------------
{
	// reset prepared status
	prepared = false;

	// store standard deviation
	stdev = abs(sd);
}


//-----------------------------------------------------------------------------
void cRandNormal::setNoOfIntervals(
	unsigned n)	///< [in] number of linear equi-distant steps for normal distribution
//-----------------------------------------------------------------------------
{
	// reset prepared status
	prepared = false;

	// store number of intervals
	intervals = n;
}


//-----------------------------------------------------------------------------
void cRandNormal::setWidth(
	double nStdev)	///< [in] width of normal distribution in standard deviations
//-----------------------------------------------------------------------------
{
	// reset prepared status
	prepared = false;

	// store number of standard deviations
	this->nStdev = abs(nStdev);
}


//-----------------------------------------------------------------------------
void cRandNormal::prepare()
//-----------------------------------------------------------------------------
{
	unsigned i;				// local index variable
	double value;			// actual value

	// reset prepared status
	prepared = false;

	// reset previous distribution
	cRandCont::reset();

	// set distribution
	for (i = 0; i <= intervals; i++) {
		value = i * nStdev / intervals - nStdev / 2;
		addDistPoint(stdev * value + mean, exp(-0.5 * value * value));
	}

	// prepare generator
	if(stdev>0) cRandCont::prepare();
	else prepared = true;
}


//-----------------------------------------------------------------------------
double cRandNormal::rand()
//-----------------------------------------------------------------------------
{
	// return random value if prepared and zero otherwise
	return prepared ? (stdev>0?cRandCont::rand():mean) : 0;
}


//-----------------------------------------------------------------------------
void cRandNormal::copy(
	const cRandNormal &source)	///< [in] generator to be copied
//-----------------------------------------------------------------------------
{
	// copy all elements of this class
	intervals = source.intervals;
	mean = source.mean;
	stdev = source.stdev;
	nStdev = source.nStdev;

	// copy all elements of base class
	cRandCont::copy(source);
}
