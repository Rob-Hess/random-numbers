//=============================================================================
// filename: random.h
// version: 1.06
// date: June 16th 2018
// author: Robert Hess, HAW-Hamburg, Germany
// description: a set of random number generators
//=============================================================================

// compile header only once
#pragma once
#ifndef LIBRARY_RANDOM_BY_ROBERT_HESS_INCUDED
#define LIBRARY_RANDOM_BY_ROBERT_HESS_INCUDED

/** \file random.h
\brief Header-file for a set of random number generators.

See description of class cRandU, cRandCont, cRandContAR, and cRandDisc for
more details. */

#include <vector>
#include <fstream>


//=============================================================================
/// \brief Base for random number generator with uniform distribution on uint32_t.

/** Base for random number generator for numbers from zero to 2^32 - 1 of type
uint32_t, i.e. values from 0 to 4,294,967,295. This is the base class for
elementary generators like cRandKISS. */
//=============================================================================
class cRandBase {
public:
	/// \brief Default empty destructor
	virtual ~cRandBase() {};
	/// \brief Initialises the generator to an arbitrary point.
	virtual void  srand(uint32_t seed) = 0;
	/// \brief Generates a uniform 32 bit unsigned random number.
	virtual uint32_t rand(void) = 0;
	/// \brief Generates a uniform random number on the interval [0,1).
	virtual double randU(void) = 0;
	/// \brief Get current state of the generator
	virtual void getState(std::vector<uint8_t> &state) = 0;
	/// \brief Set the state of the generator
	virtual void setState(const std::vector<uint8_t> &state) = 0;
};


//=============================================================================
/// \brief Random number generator with uniform distribution on uint32_t.

/** Random number generator for numbers from zero to 2^32 - 1 of type uint32_t,
i.e. values from 0 to 4,294,967,295. This class serves as a base for some other
random generators.

A seed may be placed by the method srand(). At start the seed zero is placed.
See srand() for more explanation on preparing the generator. */
//=============================================================================
class cRandKiss : public cRandBase
{
public:
	/// \brief Default constructor.
	cRandKiss();
	/// \brief Copy constructor.
	cRandKiss(const cRandKiss &source);
	/// \brief Initialises the generator with the given seed-value.
	cRandKiss(uint32_t seed);
	/// \brief Initialises the generator to an arbitrary point.
	void srand(uint32_t seed);
	/// \brief Generates a uniform 32 bit unsigned random number.
	uint32_t rand(void);
	/// \brief Generates a uniform random number on the interval [0,1).
	double randU(void);
	/// \brief Get current state of the generator
	void getState(std::vector<uint8_t> &state);
	/// \brief Set the state of the generator
	void setState(const std::vector<uint8_t> &state);
private:
	uint32_t x32;
	uint32_t y32;
	uint32_t z32;
	uint32_t c32;
};


//=============================================================================
/// \brief Base for random variables containing a probability function.
//=============================================================================
class cRandDistribution {
private:
	std::vector<uint32_t> cumulative;		///< cumulative distribution function
protected:
	cRandBase *pRand;						///< pointer to uniform random number generator
	std::vector<double> dist[2];			///< target distribution of generator
	unsigned maxInvIndex;					///< highest index of integrated inverse distribution
	unsigned step0;							///< fist step to search index
	bool prepared;							///< indicator if generator is prepared to be handled by derived class
	std::vector<uint32_t>::pointer mapInp;	///< pointer to x-values of integrated inverse distribution
	std::vector<double>::pointer mapOut;	///< pointer to x-values of integrated inverse distribution

public:
	/// \brief Standard Constructor.
	cRandDistribution();
	/// \brief Copy constructor.
	cRandDistribution(const cRandDistribution &source);
	/// \brief Constructor with base generator.
	cRandDistribution(cRandBase &rand);
	/// \brief sets the random generator base
	void setRandomBase(cRandBase &rand);
	/// \brief Initialises the generator to an arbitrary point.
	void srand(uint32_t seed);
	/// \brief Resets the present distribution.
	void reset();
	/// \brief copies an object of this class
	cRandDistribution &operator=(cRandDistribution &orig);
	/// \brief Add a point to the probability distribution
	void addDistPoint(double x, double prob);
	/// \brief Sets the desired distribution.
	void setDistribution(const std::vector<double> &x, const std::vector<double> &p);
	/// \brief Sets the desired distribution.
	void setDistribution(const std::vector<double> distribution[2]);
protected:
	/// \brief Prepares the mapping function
	void prepare(bool discrete);
	/// \brief Copies an object of this class
	void copy(const cRandDistribution &source);

#ifdef TEST_RANDOM_CLASSES
	friend void testRandDistribution();
#endif
};


//=============================================================================
/// \brief Random number generator for arbitrary discrete probability distribution.

/// Random number generator which generates only numbers and probabilities
/// defined before.
//=============================================================================
class cRandDisc : public cRandDistribution {
public:
	/// \brief Standard Constructor.
	cRandDisc();
	/// \brief Copy constructor.
	cRandDisc(const cRandDisc &source);
	/// \brief Constructor with base generator.
	cRandDisc(cRandBase &rand);
	/// \brief Resets the present distribution.
	void reset();
	/// \brief Returns the current state of the generator
	bool isPrepared() { return prepared; };
	/// \brief copies an object of this class
	cRandDisc &operator=(cRandDisc &source);
	/// \brief Prepares the random number generator.
	void prepare();
	/// \brief Generates a random number as part of the desired distribution.
	double rand();
};


//=============================================================================
/// \brief Discrete random number generator with parameter.

/** Random number generator for arbitrary discrete probability distributions.
The generator gets several distributions with parameters. When generating a
random number a parameter is required to perform linear interpolation between
distributions.

The possible output values are defined by setRandomValues(). The probabilities
for these values are defined by addProbabilities(). This random number
generator is based on a uniform random number generator, e.g. cRandKiss
which must be passed by reference. */
//=============================================================================
class cRandDiscParam
{
private:
	cRandBase *pRand;						///< pointer to uniform random number generator
	std::vector<double> value;				///< list of possible random values
	std::vector<std::vector<double>> prob;	///< probabilities for different parameters
	std::vector<double> parameter;			///< list of parameters
	bool prepared;							///< indicator if generator is prepared
	std::vector<std::vector<double>> cum;	///< cumulative distributions
	unsigned maxIndexParam;					///< maximum index for parameter
	unsigned stepParam0;					///< initial step size to find index for parameter
	unsigned maxIndexProb;					///< maximum index for probability
	unsigned stepProb0;						///< initial step size to find index for probability
public:
	/// \brief Default constructor.
	cRandDiscParam();
	/// \brief Copy constructor.
	cRandDiscParam(const cRandDiscParam &source);
	/// \brief Constructor with base random generator.
	cRandDiscParam(cRandBase &rand);
	/// \brief Resets the present distribution.
	void reset();
	/// \brief Returns the current state of the generator
	bool isPrepared() { return prepared; };
	/// \brief Sets the uniform uint32_t base random generator
	void setRandomBase(cRandBase &rand);
	/// \brief Sets the possible random values
	void setRandomValues(const std::vector<double> &x);
	/// \brief Sets a new set of probabilities
	void addProbabilities(const std::vector<double> &p, double param);
	/// \brief Prepares the random number generator.
	void prepare();
	/// \brief Initialises the generator to an arbitrary point.
	void  srand(uint32_t seed);
	/// \brief Generates a random number as part of the desired distribution.
	double rand(double param);
	/// \brief Generates a random number as part of the desired distribution.
	double rand(double randU, double param);
private:
	/// \brief Copies an object of this class
	void copy(const cRandDiscParam &source);
};


//=============================================================================
/// \brief Random number generator for arbitrary continuous probability distribution.

/// Random number generator which generates numbers with the given probability density function.
//=============================================================================
class cRandCont : public cRandDistribution {
private:
	std::vector<double> a1;		///< coefficient 1 for interpolation
	std::vector<double> a2;		///< coefficient 2 for interpolation
	std::vector<double> a3;		///< coefficient 3 for interpolation
	std::vector<double> a4;		///< coefficient 4 for interpolation
	std::vector<bool> linear;	///< flag for linear/polynomial interpolation
public:
	/// \brief Default constructor.
	cRandCont();
	/// \brief Copy constructor.
	cRandCont(const cRandCont &source);
	/// \brief Constructor with base random generator.
	cRandCont(cRandBase &rand);
	/// \brief copy operator.
	cRandCont &operator=(const cRandCont &source);
	/// \brief Resets the present distribution.
	void reset();
	/// \brief Returns the current state of the generator
	bool isPrepared() { return prepared; };
	/// \brief Prepares the random number generator.
	void prepare();
	/// \brief Generates a random number as part of the desired distribution.
	double rand();
	/// \brief Generates a random number with given uniform input.
	double rand(const uint32_t &input);

protected:
	/// \brief Copies an object of this class
	void copy(const cRandCont &source);
};


//=============================================================================
/// \brief Continuous random number generator with parameter.

/** Random number generator for arbitrary continuous probability distributions.
The generator gets several distributions with parameter. When generating a
random number a parameter is required to perform linear interpolation between
the distributions. This generator is based on a vector of the continuous
generators cRandCont.

The probability distribution is defined with setDistribution(). This random
number generator is based on a uniform random number generator, e.g. cRandKiss
which must be passed by reference. */
//=============================================================================
class cRandContParam
{
private:
	cRandBase *pRand;					///< pointer to uniform random number generator
	std::vector<cRandCont> generator;	///< a vector of generators for each parameter
	std::vector<double> parameter;		///< list of parameter values
	bool prepared;						///< indicator if generator is prepared
	unsigned maxInd;					///< maximum index of parameter
	unsigned step0;						///< initial step size to find index of parameter
public:
	/// \brief Default constructor.
	cRandContParam();
	/// \brief Copy constructor.
	cRandContParam(const cRandContParam &source);
	/// \brief Constructor with base random generator.
	cRandContParam(cRandBase *pRandBase);
	/// \brief Resets the present distribution.
	void reset();
	/// \brief Returns the current state of the generator
	bool isPrepared() { return prepared; };
	/// \brief Sets the uniform uint32_t base random generator
	void setRandomBase(cRandBase &rand);
	/// \brief Sets a desired distribution.
	void setDistribution(const std::vector<double> &x, const std::vector<double> &p, double param);
	/// \brief Sets a desired distribution.
	void setDistribution(const std::vector<double> distribution[2], double param);
	/// \brief Sets a point for the desired distribution
	void add(double x, double p, double param);
	/// \brief Prepares the random number generator.
	void prepare();
	/// \brief Initialises the generator to an arbitrary point.
	void  srand(uint32_t seed);
	/// \brief Generates a random number as part of the desired distribution.
	double rand(double param);
private:
	/// \brief Copies an object of this class
	void copy(const cRandContParam &source);
};


//=============================================================================
/// \brief Continuous random number generator with normal distribution.

/** Random number generator with exponential distribution using the natural
logarithm. Expectation is adjustable with default value 1. A large number of
random numbers will give the expected exponential distribution. */
//=============================================================================
class cRandExponential
{
private:
	cRandBase *pRand;		///< pointer to uniform random number generator
	double nmean;			///< expectation of random numbers with negative sign
public:
	/// \brief Standard constructor.
	cRandExponential();
	/// \brief Copy constructor.
	cRandExponential(const cRandExponential &source);
	/// \brief Constructor with uniform generator.
	cRandExponential(cRandBase &rand);
	/// \brief Resets the present distribution.
	void reset();
	/// \brief copies an object of this class
	cRandExponential &operator=(cRandExponential &source);
	/// \brief Sets the random generator base
	void setRandomBase(cRandBase &rand);
	/// \brief Sets the mean value of the normal distribution
	void setExpectation(double expectation);
	//	/// \brief Prepares the random number generator.
	//	void prepare();
	/// \brief Generates a random number as part of the desired distribution.
	double rand();
private:
	/// \brief Copies an object of this class
	void copy(const cRandExponential &source);

#ifdef TEST_RANDOM_CLASSES
	friend void testRandExponential();
#endif
};


//=============================================================================
/// \brief Continuous random number generator with normal distribution.

/** Random number generator with normal distribution using the accept or
reject method. Expectation, standard deviation, width of distribution and
number of intervals are adjustable with default values 0, 1, 10 and 100,
respectively. A large number of random numbers will give a normal
distribution. */
//=============================================================================
class cRandNormal : private cRandCont
{
private:
	unsigned intervals;		///< number of intervals for interpolation
	double mean;			///< expectation of random numbers
	double stdev;			///< standard deviation of normal distribution
	double nStdev;			///< total width of distribution in stdev
public:
	/// \brief Standard constructor.
	cRandNormal();
	/// \brief Copy constructor.
	cRandNormal(const cRandNormal &source);
	/// \brief Constructor with uniform generator.
	cRandNormal(cRandBase &rand);
	/// \brief Resets the present distribution.
	void reset();
	/// \brief Returns the current state of the generator
	bool isPrepared() { return prepared; };
	/// \brief copies an object of this class
	cRandNormal &operator=(cRandNormal &source);
	/// \brief Sets the random generator base
	void setRandomBase(cRandBase &rand);
	/// \brief Sets the mean value of the normal distribution
	void setExpectation(double expectation);
	/// \brief Sets the standard deviation of the normal distribution
	void setStandardDeviation(double sd);
	/// \brief Sets the number of intervals to approximate normal distribution
	void setNoOfIntervals(unsigned n);
	/// \brief Sets the maximum width of normal distribution in multiples of standard deviation
	void setWidth(double nStdev);
	/// \brief Prepares the random number generator.
	void prepare();
	/// \brief Generates a random number as part of the desired distribution.
	double rand();
private:
	/// \brief Copies an object of this class
	void copy(const cRandNormal &source);

#ifdef TEST_RANDOM_CLASSES
	friend void testRandNormal();
#endif
};


#endif // #ifndef LIBRARY_RANDOM_BY_ROBERT_HESS_INCUDED
