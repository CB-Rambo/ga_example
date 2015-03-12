/**
	@file ga_test.c Genetic algorithm example
	@brief Example of Genetic algorithm implementation.

	@authors Denny Hermawanto (Source: http://arxiv.org/ftp/arxiv/papers/1308/1308.4675.pdf)
	@authors Krishna Vedala (Implementation in C)
*/
/**
	@mainpage ga_example.c

	This is a simple implementation of a genetic algorithm used to solve
	an objective function.
*/

/*
#if OPENMP_FOUND == FALSE
#	error "Please re-run cmake from the parent directory"
#endif
*/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>


#define		MAX_ITER	1000	/**< maximum number of iterations to perform */
#define		EPSILON 	1e-6	/**< minimum value for objective function */
#define		NUM_GENES 2			/**< number of genes for this cost-function example */
#define		RHO				0.8		/**< crossover probability */
#define		MUT_RATE	0.1		/**< gene mutation rate in percentage */
#define		DEBUG			0			/**< macro to print debug messages or not */


/**
	@brief Chromosome structure. Each variable identifies a gene
	of the chromosome.
*/
typedef struct s_chromosome
{
	int *genes;							/**< pointer to an array of genes */
	unsigned int num_genes;	/**< number of genes in the chromosome */
	double objective;				/**< result of objective value for the chromosome */
	double fit;							/**< fitness value of the chromosome */
	double prob;						/**< probability of the chromosome to be selected*/
} s_chromosome;


/**
	@brief Cost function that needs to be minimized by the algorithm.
	The function uses the genes from the chromosome and returns the
	value of the cost function for that chromosome.
	@f[ f(x)=\left|(a+2b+3c+4d)-30\right|@f]
	In this cost function, the variables @f$a@f$, @f$b@f$, @f$c@f$ and
	@f$d@f$ are the four genes that form one chromosome.
	@param[in] c 	chromosome to compute cost-function for
*/
static void cost_function(s_chromosome *c)
{
	int i;
	double r = 0.F, temp = 1.F;

	for(i = 0; i < NUM_GENES; i++)
	{
		r += c->genes[i] * temp;
		temp /= (double)10.F;
	}
	r = r * r - 3.F;

	c->objective = (r > 0.0F) ? r : -r;		// absolute value
	c->fit = 1.F / (1.F + c->objective);	// fitness value
}


static inline int init_gene()
{
	return (rand() % 90);
}


/**
	@brief Function to initialize a chromosome with an array of
	genes of random values.
	@param[in,out] 	c 	s_chromosome structure to initialize
	@param[in]			num	number of genes in the chromosome
*/
static void chromosome_init(s_chromosome *c, unsigned int num)
{
	int i; 
	
	c->genes = NULL;
	c->num_genes = 0;
	c->genes = malloc(num * sizeof(int));	// allocate memory for genes
	if(!c->genes)	// if unable to allocate
	{
		printf("\n\tError initializing memory.");
		return;
	}
	for (i = 0; i < num; i++)
	{
		c->genes[i] = init_gene();		    // initialize i^th gene
	}
	c->num_genes = num;
	c->objective = 1.F;
}


/**
	@brief Function to de-initialize a chromosome and de-allocate
	@param[in,out] 	c 	s_chromosome structure to de-allocate
*/
static inline void chromosome_deinit(s_chromosome *c)
{
	if(c->genes)
		free(c->genes);		// deallocate memory for genes
}


/**
	@brief Function to print genes from a chromosome
	@param[in] 	c 	s_chromosome structure to print from
*/
static inline void print_genes(s_chromosome *c)
{
	int i;
	for (i = 0; i < c->num_genes; i++)
	{
		printf("%d ", c->genes[i]);
	}
}


/**
	@brief Function to swap genetic data between two chromosome sets
	@param[in,out] 	c1 	array 1 of s_chromosome structure
	@param[in,out] 	c2 	array 2 of s_chromosome structure
	@param[in] 	N 	number of chromosomes to swap
*/
static void swap_chromosomes(s_chromosome *c1, s_chromosome *c2, int N)
{
	int g, n, i;
	for (n = 0; n < N; n++)
	{
		for (i = 0; i < c1[n].num_genes; i++)
		{
			g = c1[n].genes[i];
			c1[n].genes[i] = c2[n].genes[i];
			c2[n].genes[i] = g;
		}
		c1[n].objective = c2[n].objective;
		c1[n].fit = c2[n].fit;
	}
}


/**
	@brief Function to copy genetic data from one chromosome to another
	@param[out] c1 destination s_chromosome structure
	@param[in] 	c2 source s_chromosome structure
*/
static void copy_chromosome(s_chromosome *c1, s_chromosome *c2)
{
	int i;
	for (i = 0; i < c1->num_genes; i++)
	{
		c1->genes[i] = c2->genes[i];
	}
	c1->objective = c2->objective;
	c1->fit = c2->fit;
}


/**
	@brief Generate a floating point random number in the interval @f$[0,B)@f$.
	@param[in] 	B	upper bound of output random number
	@return 	return a random number in the interval @f$[0,B)@f$
*/
static inline float rand2(float B)
{
	B *= 1000;
	return (float)(rand() % (int)B) / 1000.F;
}


/**
	Main Function
*/
int main(int argc, char *argv[])
{
	srand(time(NULL));		// initialize random number generator

	unsigned int pop_num, i;

  if (argc == 1)
  {
	  printf("Enter the chromosome population: ");
	  scanf ("%u", &pop_num);
  } else {
		pop_num = atoi(argv[1]);
  }

#if	DEBUG == 1
	printf("\nInitializing and chromosome population...");
#endif
	s_chromosome chromos1[pop_num];
	s_chromosome chromos2[pop_num];	// temporary storage
	for (i = 0; i < pop_num; i++)
	{
		chromosome_init(chromos1 + i, NUM_GENES);
		chromosome_init(chromos2 + i, NUM_GENES);
	}
#if	DEBUG == 1
	printf("Done!\n\n");
#endif

	unsigned long iter = 1;
	double start_time, end_time;

	start_time = omp_get_wtime();
	while ((iter < MAX_ITER) && (chromos1[0].objective > EPSILON))
	{
#if	DEBUG == 1
		printf("\n\t\t\tIter# %lu\n", iter);
		for (i = 0; i < pop_num; i++)
		{
			printf("\tChromosome# %d:\t", i);
			print_genes(chromos1 + i);
			printf("\n");
		}

		printf("Evaluating...");
#endif

		double F = 0.F;
#pragma omp parallel shared (F, chromos1) private(i)
		{
		#pragma omp for schedule(dynamic) nowait
		for (i = 0; i < pop_num; i++)
		{
			cost_function(chromos1 + i);	// compute cost function for each chromosome
			F += chromos1[i].fit;
		}
		}
#if	DEBUG == 1
		printf("Done:\t\t%f\n", chromos1[0].fit);
#endif
		for(i = 0; i < pop_num; i++)
			if (chromos1[i].objective <= EPSILON)
			{
				copy_chromosome(chromos1, chromos1 + i);
				goto DONE_GA;
			}
#if	DEBUG == 1
		printf("Selection probabilities...");
#endif
		double C = 0.F;
#pragma omp parallel shared (C, chromos1) private(i)
		{
		#pragma omp for schedule(dynamic) nowait
		for (i = 0; i < pop_num; i++)		// Chromosome selection probabilities
		{
			chromos1[i].prob = chromos1[i].fit / F;
			C += chromos1[i].prob;
			chromos1[i].prob = C;
		}
		}
#if	DEBUG == 1
		printf("Done\n");
		printf("Selecting...");
#endif
#pragma omp parallel shared (chromos1, chromos2) private(i)
		{
		#pragma omp for schedule(dynamic) nowait
		for (i = 0; i < pop_num; i++)
		{
			float R = rand2(C);		// generate a random number between 0 and C
			int ii;
			for (ii = 0; ii < pop_num; ii++)
				if (R <= chromos1[ii].prob)
					break;
			copy_chromosome(chromos2 + i, chromos1 + ii);
		}
		}
		swap_chromosomes(chromos1, chromos2, pop_num);
#if	DEBUG == 1
		printf("Done\n");
		printf("Applying crossover...");
#endif
		unsigned int R[pop_num];
		for (i = 0; i < pop_num; i++)
		{
			float r = rand2(1.F);		// generate a random number between 0 and 1.0
			R[i] = (r < RHO) ? 0 : 1;
		}
		unsigned int posi = rand() % NUM_GENES;		// crossing over position
#pragma omp parallel shared (R, chromos1) private(i)
		{
		#pragma omp for schedule(dynamic)
		for (i = 0; i < pop_num; i++)
		{
			if (R[i])
			{
				int ii = i, j;
				do
				{
					ii++;
					if (ii >= pop_num)	ii = 0;
				} while(!R[ii] && (ii != i));
				for (j = posi; j < NUM_GENES; j++)
        {
          int temp = chromos1[i].genes[j];
					chromos1[i].genes[j] = chromos1[ii].genes[j];
          chromos1[ii].genes[j] = temp;
        }
			}
		}
		}
#if	DEBUG == 1
		printf("Done\n");
		printf("Mutating...");
#endif
		unsigned long L = pop_num * NUM_GENES;
		for (i = 0; i < floor(L * MUT_RATE); i++)
		{
			unsigned long m = rand2(L);
			unsigned long q, r;
			q = floor(m / NUM_GENES);
			r = m % NUM_GENES;
			chromos1[q].genes[r] = init_gene();
		}
#if	DEBUG == 1
		printf("Done\n");
#endif
		iter++;
	}

DONE_GA:
	
	end_time = omp_get_wtime();
	if (chromos1[0].objective < EPSILON)
	{
		printf("\nFound solution in %lu iterations.\n", iter);
		printf("Solution chromosome (Fitness: %g%%):\n", chromos1[0].fit * 100.F);
		print_genes(chromos1);
		printf("\n");
	} else
		printf("\nCould not converge to solution after %lu iterations!\n", iter);
	printf("Actual elapsed time = %.6g s\n", end_time - start_time);

#if	DEBUG == 1
	printf("\nDe-initializing chromosome population...");
#endif
	for (i = 0; i < pop_num; i++)
	{
		chromosome_deinit(chromos1 + i);
		chromosome_deinit(chromos2 + i);
	}
#if	DEBUG == 1
	printf("Done!\n");
#endif
	return 0;
}
