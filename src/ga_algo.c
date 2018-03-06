/**
 * @file
 * 
 * @authors Krishna Vedala (Implementation in C)
 */

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "macrologger.h"
#include "ga_algo.h"



bool init_gene(void)
{
	if (rand() % 2 == 0)
		return true;
	else
		return false;
}

/**
 * @brief Function to initialize a chromosome with an array of
 * genes of random values.
 * @param[in,out] 	c 	s_chromosome structure to initialize
 * @param[in]			num	number of genes in the chromosome
 */
void chromosome_init(chromosome c, unsigned int num)
{
	unsigned int i;

 	// Log and print something with level 0.
    LOG_DEBUG("DEBUG message");
    // Log and print something with level 1.
    LOG_INFO("INFO message");
    // Log and print something with level 2.
    LOG_ERROR("ERROR message");

	c->genes = NULL;
	c->num_genes = 0;
	c->genes = (bool*) malloc(num * sizeof(int));	// allocate memory for genes
	if(!c->genes)	// if unable to allocate
	{
		perror("\n\tError initializing memory.");
		exit(-1);
	}

	#pragma omp parallel for
	for (i = 0; i < num; i++)
	{
		c->genes[i] = init_gene();		    // initialize i^th gene
	}
	c->num_genes = num;
	c->objective = 1.F;
}

/**
 * @brief Function to de-initialize a chromosome and de-allocate
 * @param[in,out] 	c 	s_chromosome structure to de-allocate
 */
void chromosome_deinit(chromosome c)
{
	if(c->genes)
		free(c->genes);		// deallocate memory for genes
}


/**
 * @brief Function to swap genetic data between two chromosome sets
 * @param[in,out] 	c1 	array 1 of s_chromosome structure
 * @param[in,out] 	c2 	array 2 of s_chromosome structure
 * @param[in] 	N 	number of chromosomes to swap
 */
void swap_chromosomes(chromosome c1, chromosome c2, int N)
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
 * @brief Function to copy genetic data from one chromosome to another
 * @param[out] c1 destination s_chromosome structure
 * @param[in] 	c2 source s_chromosome structure
 */
void copy_chromosome(chromosome c1, chromosome c2)
{
	int i;
	for (i = 0; i < c1->num_genes; i++)
	{
		c1->genes[i] = c2->genes[i];
	}
	c1->objective = c2->objective;
	c1->fit = c2->fit;
}


double decode_d(bool* gene)
{
	return (double)(*gene);
}

void encode_d(double val, bool* gene)
{
	memcpy((char*) gene, (char*) &val, sizeof(val));
}
