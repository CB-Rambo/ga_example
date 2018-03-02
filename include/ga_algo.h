/** 
 * @file
 * 
 * @authors Krishna Vedala (Implementation in C)
 */

#ifndef __GAALGO__
#define __GAALGO__

#include <stdbool.h>
#include <stdint.h>

#define		MAX_ITER	1000	/**< maximum number of iterations to perform */
#define		EPSILON 	1e-10	/**< minimum value for objective function */
#define		NUM_GENES 	20		/**< number of genes for this cost-function example */
#define		RHO			0.6		/**< crossover probability */
#define		MUT_RATE	0.1		/**< gene mutation rate in percentage */
#define		DEBUG		0		/**< macro to print debug messages or not */


/**
	@brief Chromosome structure. Each variable identifies a gene
	of the chromosome.
*/
struct s_chromosome
{
	bool *genes;		/**< pointer to an array of genes */
	uint16_t num_genes;	/**< number of genes in the chromosome */
	double objective;	/**< result of objective value for the chromosome */
	double fit;			/**< fitness value of the chromosome */
	double prob;		/**< probability of the chromosome to be selected*/
};
typedef struct s_chromosome s_chromosome;

typedef s_chromosome* chromosome;

/**
	@brief Cost function that needs to be minimized by the algorithm.
	The function uses the genes from the chromosome and returns the
	value of the cost function for that chromosome.
	@f[ f(x)=\left|(a+2b+3c+4d)-30\right|@f]
	In this cost function, the variables @f$a@f$, @f$b@f$, @f$c@f$ and
	@f$d@f$ are the four genes that form one chromosome.
	@param[in] c 	chromosome to compute cost-function for
*/
void cost_function(chromosome c);

/**
	@brief Function to initialize a chromosome with an array of
	genes of random values.
	@param[in,out] 	c 	s_chromosome structure to initialize
	@param[in]			num	number of genes in the chromosome
*/
void chromosome_init(chromosome c, unsigned int num);

/**
	@brief Function to de-initialize a chromosome and de-allocate
	@param[in,out] 	c 	s_chromosome structure to de-allocate
*/
void chromosome_deinit(chromosome c);

/**
	@brief Function to swap genetic data between two chromosome sets
	@param[in,out] 	c1 	array 1 of s_chromosome structure
	@param[in,out] 	c2 	array 2 of s_chromosome structure
	@param[in] 	N 	number of chromosomes to swap
*/
void swap_chromosomes(chromosome c1, chromosome c2, int N);


/**
	@brief Function to copy genetic data from one chromosome to another
	@param[out] c1 destination s_chromosome structure
	@param[in] 	c2 source s_chromosome structure
*/
void copy_chromosome(chromosome c1, chromosome c2);


float decode_f(bool*);
void encode_f(float, bool*);
double decode_d(bool*);
void encode_d(double, bool*);
int decode_i(bool*);
void encode_i(int, bool*);

void cost_function(chromosome c);

#endif
