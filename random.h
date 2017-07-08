//=============================================================================
// filename: random.h
// version: 1.01
// date: 13.8.2014
// author: Robert Hess
// description: a set of random number generators
//=============================================================================

#ifndef LIBRARY_RANDOM_BY_ROBERT_HESS_INCUDED
#define LIBRARY_RANDOM_BY_ROBERT_HESS_INCUDED

/** \file random.h
\brief Header-file for a set of random number generators.

See description of class cRandU, cRandCont, cRandContAR, and cRandDisc for
more details. */

#include <vector>
#include <fstream>

/// \brief number of random value for shufle mechanism
#define RAND_NTAB 32
//===========================================================================
/// \brief Random number generator with uniform distribution between zero an one.

/** Random number generator for numbers from (including) zero to (excluding)
one of type float, i.e. values from 0 to 0.999... .

<b>Note:</b> An application should contain only one object of this class. This
avoids unwanted correlations between random numbers.

After construction this generator may be used right away. However, see sRand()
for more explanation on preparing the generator. */
//===========================================================================
class cRandU {
public:
	/// \brief Default constructor.
	cRandU();
	/// \brief Initialises the generator with the given seed-value.
	cRandU(long seed);
	/// \brief Generates a random number from the interval [0,1).
    float rand(void);
	/// \brief Initialises the generator to an abitrary point.
    void  sRand(long seed);
	/// \brief Stores the current state of the generator in a file.
    bool  save(void);
	/// \brief Restores the current state of the generator from a file.
    bool  restore(void);
private:
    long idum2;
    long iy;
    long iv[RAND_NTAB];
    long idum;
	friend std::ofstream &operator<<(std::ofstream &fout, cRandU &randU);
	friend std::ifstream &operator>>(std::ifstream &finp, cRandU &randU);
};
/// \brief Stores the actual state of the generator in the given output stream.
std::ofstream &operator<<(std::ofstream &fout, cRandU &randU);
/// \brief Restores the actual state of the generator from the given input stream.
std::ifstream &operator>>(std::ifstream &finp, cRandU &randU);

//===========================================================================
/// \brief Continuous random number generator (accept/reject method).

/** Random number generator for abitrary continuous probability distributions
implementing the accept or reject method. The probability distribution is
defined with setDistribution(). A large number of random numbers will give
an distribution as given with setDistribution() with a linear interpolation
between the given points. This random number generator is based on the
uniform random number generator cRandU which must be passed by reference
during construction. */
//===========================================================================
class cRandContAR {
public:
	/// \brief Constructor.
	cRandContAR(cRandU *pRandU);
	/// \brief Destructor.
	~cRandContAR();
	/// \brief Resets the present distribution.
	void reset();
	/// \brief Sets the desired distribution.
	bool setDistrib(const std::vector<double> distribution[2]);
	/// \brief Prepares the random number generator.
	bool prepare();
	/// \brief Generates a random number as part of the desired distribution.
	double rand(void);
private:
	cRandU *pRandU;							///< pointer to uniform random number generator
	std::vector<double> dist[2];			///< target distribution of generator
	std::vector<double>::pointer distY;		///< pointer to y-values of distribution
	std::vector<double> inv[2];				///< target distr. integrated and inversed
	std::vector<double>::pointer invX;		///< pointer to x-values of integrated inversed distribution
	std::vector<double>::pointer invY;		///< pointer to y-values of integrated inversed distribution
	int maxInvIndex;						///< highest index of integrated inversed distribution
	bool error;								///< error indicator
	bool prepared;							///< indicator if generator is prepared
	unsigned step0;							///< fist step to search index
};


//===========================================================================
/// \brief Continuous random number generator (step distribution).

/// Random number generator for abitrary continuous probability distributions
/// implementing a stepwise rectangular distribution.
//===========================================================================
class cRandCont {
public:
	/// \brief Standard Constructor.
	cRandCont();
	/// \brief Constructor with uniform generator.
	cRandCont(cRandU *pRandU);
	/// \brief Copy constructor.
	cRandCont(cRandCont &orig);
	/// \brief Destructor.
	~cRandCont();
	/// \brief Resets the present distribution.
	void reset();
	/// \brief copies an object of this class
	cRandCont &operator=(cRandCont &orig);
	/// \brief Sets the desired distribution.
	bool setDistrib(const std::vector<double> distribution[2]);
	/// \brief Prepares the random number generator.
	bool prepare();
	/// \brief Prepares the random number generator with given accuracy.
	bool prepareRelAcc(double acc=0.01);
	/// \brief Generates a random number as part of the desired distribution.
	double rand(void);
	/// \brief Converts a uniform random value into a desired random variable.
	double rand(float randU);
private:
	cRandU *pRandU;							///< pointer to uniform random number generator
	std::vector<double> dist[2];			///< target distribution of generator
	std::vector<double> inv[2];				///< target distr. integrated and inversed
	std::vector<double>::pointer invX;		///< pointer to x-values of integrated inversed distribution
	std::vector<double>::pointer invY;		///< pointer to y-values of integrated inversed distribution
	int maxInvIndex;						///< highest index of integrated inversed distribution
	bool error;								///< error indicator
	bool prepared;							///< indicator if generator is prepared
	unsigned step0;							///< fist step to search index

	/// \brief Copies an object of this class
	void copy(cRandCont &orig);
};


//===========================================================================
/// \brief Random number generator for abitrary discrete probability distributions.

/// Random number generator which generates only numbers and probabilities
/// defined before.
//===========================================================================
class cRandDisc {
public:
	/// \brief Constructor.
	cRandDisc(cRandU *pRandU);
	/// \brief Destructor.
	~cRandDisc();
	/// \brief Resets the present distribution.
	void reset();
	/// \brief Sets the desired distribution.
	bool setDistrib(const std::vector<double> distribution[2]);
	/// \brief Adds a point to the distribution
	bool add(double x, double y);
	/// \brief Prepares the random number generator.
	bool prepare();
	/// \brief Generates a random number as part of the desired distribution.
	double rand(void);
private:
	cRandU *pRandU;							///< pointer to uniform random number generator
	std::vector<double> dist[2];			///< target distribution of generator
	std::vector<double> inv[2];				///< target distr. integrated and inversed
	std::vector<double>::pointer invX;		///< pointer to x-values of integrated inversed distribution
	std::vector<double>::pointer invY;		///< pointer to y-values of integrated inversed distribution
	int maxInvIndex;						///< highest index of integrated inversed distribution
	bool error;								///< error indicator
	bool prepared;							///< indicator if generator is prepared
};

#endif // #ifndef LIBRARY_RANDOM_BY_ROBERT_HESS_INCUDED
