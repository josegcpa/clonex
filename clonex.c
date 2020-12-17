/*
clonex.c

Version : 0.1.23
Author  : Niko Beerenwinkel
    Moritz Gerstung

(c) Niko Beerenwinkel, Moritz Gerstung 2009
 */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <getopt.h>
#include <omp.h>
#include "merge.c"

//#define genes 1000000 // dimension (number of loci)
#define MAX_GENOTYPES 200000  // max. 10^6
#define MAX_K 1000  // max. no. of mutations per genotype


typedef int (*compfunc) (const void*, const void*);  // type for comparison function in qsort


struct Genotype
{
	int mutation[MAX_K+1];
	int k;  // number of mutations
	int count;
};


struct Genotype genotype[MAX_GENOTYPES];
double geno_prob[MAX_GENOTYPES];
unsigned int geno_count[MAX_GENOTYPES];
int N_g;  // current number of genotypes (make sure that N_g <= MAX_GENOTYPES)
double fitness[MAX_K+1];

int k_abs_freq[MAX_K+1];
double k_rel_freq[MAX_K+1];

gsl_rng *RNG;  // random number generator


void summary_to_R(FILE *DB)
{
	int i, k;
    fprintf(DB, "structure(list(n=c(");
	for (i=0; i<N_g; i++)
	{
		fprintf(DB, "%d", genotype[i].count);
		if(i<N_g-1)
			fprintf(DB, ", ");
	}
	fprintf(DB, "), genotype=c(");
	for (i=0; i<N_g; i++)
	{
		fprintf(DB, "\"");
		for (k=0; k<genotype[i].k; k++)
		{
			fprintf(DB, "%d", genotype[i].mutation[k]);
			if(k<genotype[i].k-1)
				fprintf(DB, ":");
		}
		fprintf(DB, "\"");
		if(i<N_g-1)
			fprintf(DB, ", ");
	}
	fprintf(DB, ")), row.names=c(NA, %i), class=\"data.frame\")", N_g);
}

void summary_to_tsv(FILE *DB, int gen)
{
	int i, k;
    //fprintf(DB, "n\tgenotype\n");
	for (i=0; i<N_g; i++)
	{
		//fprintf(DB, "%d\t%d\t", gen, genotype[i].count);
		if(genotype[i].k==0)
			fprintf(DB, "%d\t%d\t%d\t%d\n", gen, genotype[i].count, i+1, 0);
		for (k=0; k<genotype[i].k; k++)
		{
			fprintf(DB, "%d\t%d\t%d\t%d\n", gen, genotype[i].count, i+1, genotype[i].mutation[k]);
		}
		//fprintf(DB, "\n");
	}
}

void fitness_distribution_to_tsv(int d,float *fitness_effects,FILE *file_handle) {
    int i;
    for (i=1; i<=d; i++){
        fprintf(file_handle,"%d,%f\n",i,fitness_effects[i-1]);
    }
}

int count_cmp (const struct Genotype *g, const struct Genotype *h)
{
	return (h->count - g->count);  // decreasing
}


int no_of_mut_cmp (const struct Genotype *g, const struct Genotype *h)
{
	return (h->k - g->k);  // decreasing
}


int pos_zero_cmp (const struct Genotype *g, const struct Genotype *h)
{
	if ((h->count > 0) && (g->count == 0))
		return 1;
	if ((h->count == 0) && (g->count > 0))
		return -1;

	return 0;
}


void remove_zeros(int total_sort)
{

	if (!total_sort) //normally (0) only sort out zeros
		qsort(genotype, N_g, sizeof(struct Genotype), (compfunc) pos_zero_cmp);
	else
		qsort(genotype, N_g, sizeof(struct Genotype), (compfunc) count_cmp);

	while (genotype[N_g-1].count == 0)
		N_g--;

}

void remove_zeros_fast()
{
	int i;
		for(i=0; i<N_g;i++){
			if(genotype[i].count==0){
				while(genotype[N_g-1].count==0) N_g--;
				if(N_g == i) break;
				genotype[i] = genotype[N_g-1];
				N_g--;
			}
		}
}


int int_cmp (const int *a, const int *b)
{
	return (*a - *b);  // increasing
}


int mutation_cmp (const struct Genotype *g, const struct Genotype *h)
{
	int j;

	if (h->k != g->k)
		return (h->k - g->k);  // decreasing
	else
	{
		for (j=0; j<g->k; j++)
		{
			if (g->mutation[j] < h->mutation[j])
				return -1;
			if (g->mutation[j] > h->mutation[j])
				return 1;
		}
		return 0;
	}


}

void count_N()
{
	int i,j;
	j=0;
	for (i=0; i<N_g; i++)
		{
			j+=genotype[i].count;
		}
	printf("N=%i\n",j);

}

void remove_duplicates()
{
	int i;

	if(N_g > 10000){
	int split = (int) floor((double)N_g/2);

	#pragma omp parallel num_threads(2) // pre-sort on two cores
		{
	#pragma omp single
			{
	#pragma omp task
				{qsort(genotype, split, sizeof(struct Genotype), (compfunc) mutation_cmp);}  // sort population
	#pragma omp task
				{qsort(genotype + split, N_g - split, sizeof(struct Genotype), (compfunc) mutation_cmp);}  // sort population}
			}
		}
	}

	mergesort(genotype, N_g, sizeof(struct Genotype), (compfunc) mutation_cmp);  // sort population

	for (i=0; i<N_g-1; i++)
	{
		if (mutation_cmp(&genotype[i], &genotype[i+1]) == 0)  // collect counts
		{
			genotype[i+1].count += genotype[i].count;
			genotype[i].count = 0;
		}
	}

	remove_zeros_fast();
}

void insert_mutation(struct Genotype *g, int mutation){ // this makes use of the fact that g->mutation is sorted
	int i=0;
	int pos = g->k;
	for(i=0;i<pos;i++){
		if(g->mutation[i] == mutation) // duplicate; do nothing
			return;
		if(g->mutation[i] > mutation)
			pos=i; // find position to insert, break
	}
	(g->k)++; // increment count
	for(i=g->k -1; i>pos; i--){
		g->mutation[i] = g->mutation[i-1]; // move one up
	}
	g->mutation[pos]=mutation; // insert
}

int isNumeric(const char *str) 
{
    while(*str != '\0')
    {
        if(*str < '0' || *str > '9' )
            return 0;
        str++;
    }
    return 1;
}


int simulate(FILE* DB, int N_init, int N_max, double r, int gen_max, double u, double v, double s, double s1, int run, int genes, int d0, int d1, int verbose, int G, int X0, int X1, double Xr, int Xm, int gd, float *fitness_effects)
{

	int gen, i, j, c, N; //k;
	int mutation, index_mut; //, mutation_number, free_slots;
	int mut[MAX_K + 1], mut_k;

	double p, N_exp_growth;
	//double v = u;
	//double k_mean, p_more, k_plus, k_minus, prob_add, homo, div;
	//int k_min, k_max, k_median;
	//int k_obs;

	//double a = exp( ( log(N_max) - log(N_init) ) / gen_max );  // exponential growth rate
	/* ensures growth in gen_max generations from N_init to N_fin */

	//a = exp( log(2) / 60 );  // doubling time is 60 generations
	//a = exp( log(2) / 30 );  // doubling time is 30 generations

	//printf("a = %f\n", a);
	printf("doubling time = %f generations\n", log(2) / log(r));
  if (X0 < gen_max) {
      printf("clones will decrease between generation %d and %d\n", X0, X1);
  }

	for (j=0; j<MAX_K; j++)
		fitness[j] = pow(1.0 + s, (double) j);


	// initialize:
	N = N_init;
	N_exp_growth = (double) N;

	for (j=0; j<MAX_K; j++)
		genotype[0].mutation[j] = 0;
	genotype[0].k = 0;
	genotype[0].count = N;
	N_g = 1;

	// grow population:
	gen = 0;
	do
	{
		gen++;

		/* population growth */

		N_exp_growth *= r;
		N = (int) (fmin(N_max, N_exp_growth) + 0.5);
    if (gen > X0 && gen < X1)
    {
        N = N * (1.0 - Xr);
        N = fmax(N,1);
        if (Xm == 0)
        {
            N_max = N * (1.0 - Xr);
            N_max = fmax(N_max,1);
        }
    }
		//printf("%i\n",N_max);

		if (N > 2000000000)
		{
			printf("Clone has grown too large!\n");
			exit(1);
		}


		/* selection */
		double fit = 1.0;

		// compute probabilities:
    if (gd == 0) {
		  for (i=0; i<N_g; i++)
		  {
			  for (j=0;j<genotype[i].k;j++)
			  {
				  if (genotype[i].mutation[j] > 0 && genotype[i].mutation[j] <= genes - d0)
				  {
					  if (genotype[i].mutation[j] <= genes - d0 - d1 )
						  fit *= 1 + s;
					  else
						  fit *= 1 + s1;
				  }
				  //      printf("%f\t%i\n", fit, genotype[i].k);
			  }
			  //geno_prob[i] = fitness[genotype[i].k] * genotype[i].count;  // no need to normalize for gsl function
			  geno_prob[i] = fit * genotype[i].count;  // no need to normalize for gsl function
			  fit = 1.0; // reset fitness to 1.0 for next genotype.
		  }
      //we're here!
    } else {
      int mut;
		  for (i=0; i<N_g; i++)
		  {
			  for (j=0;j<genotype[i].k;j++)
			  {
          mut = genotype[i].mutation[j];
				  if (mut > 0 && mut <= genes - d0) {
            fit *= 1 + fitness_effects[mut];
          }
			  }
			  geno_prob[i] = fit * genotype[i].count;  // no need to normalize for gsl function
			  fit = 1.0; // reset fitness to 1.0 for next genotype.
		  }
    }

		for (i=0; i<N_g; i++)
			geno_count[i] = genotype[i].count;

		gsl_ran_multinomial(RNG, N_g, N, geno_prob, geno_count);

		for (i=0; i<N_g; i++)
			genotype[i].count = geno_count[i];

		remove_zeros_fast();  // because low-frequency mutants are likely not to be sampled
		// and we need to consider less genotypes for mutation


		/* mutation */

		double mu;
		int l,m;
		double a, b;

		for(m=0; m<2; m++){
			if(m == 0){ // driver
				mu = u;
				l = genes-d0;
				a = 1.0;
				b = (double)genes - (double)d0 + 1.0;

			}else{ // passenger mutation
				mu = v;
				l = d0;
				a = (double)genes - (double)d0 + 1.0;
				b = (double)genes + 1.0;

			}

			int N_mut_cells;

			//#pragma omp parallel for private(p, N_mut_cells,i,c,index_mut,mutation, mut, mut_k,j,)
			for (i=0; i<N_g; i++)
			{
				p = 1.0 - gsl_ran_binomial_pdf(0, mu, l);  // prob. of >= 1 mutation
				N_mut_cells = gsl_ran_binomial(RNG, p, genotype[i].count);  // number of mutated cells

				if (N_g + N_mut_cells > MAX_GENOTYPES)
				{
					printf("Too many genotypes: out of memory\n");
					printf("MAX_GENOTYPES = %d\n", MAX_GENOTYPES);
					summary_to_tsv(DB, gen);
					exit(1);
				}

				//printf("%i\n", N_g);
				for (c=0; c<N_mut_cells; c++)
				{
					// make a copy of genotype i for mutation, decrease count:
					(genotype[i].count)--;
					if (genotype[i].count == 0)
						index_mut = i;  // overwrite
					else
					{
						index_mut = N_g;
						genotype[index_mut] = genotype[i];  // copy
						N_g++;
					}
					genotype[index_mut].count = 1;


					if (genotype[index_mut].k - 2 > MAX_K)
					{
						printf("Warning: More than %d mutations generated!\n", MAX_K);
						genotype[index_mut].k = MAX_K - 2;  // reset
					}

					// add first mutation:
					mutation = (int) floor(gsl_ran_flat(RNG, a, b));
					if(mutation > genes) printf("Error!");


					//printf("%i\n", mutation);

					//genotype[index_mut].mutation[genotype[index_mut].k] = mutation;
					//(genotype[index_mut].k)++;

					insert_mutation(genotype + index_mut, mutation);


					// maybe add a second mutation:
					p = gsl_ran_binomial_pdf(0, mu, l);

					// prob of zero additional mutations
					if (gsl_ran_flat(RNG, 0.0, 1.0) > p)  // i.e., additional mutations, assume exactly 1
					{
						mutation = (int) floor(gsl_ran_flat(RNG, a, b));
						//printf("%i\n", mutation);
					    //genotype[index_mut].mutation[genotype[index_mut].k] = mutation;
						//(genotype[index_mut].k)++;
						insert_mutation(genotype + index_mut, mutation);
					}

					// remove duplicate mutations (mutation[] is just a list!):

					if (genotype[index_mut].k > 1)
					{
						qsort(genotype[index_mut].mutation, genotype[index_mut].k, sizeof(int), (compfunc) int_cmp);
						mut[0] = genotype[index_mut].mutation[0];
						mut_k = 1;
						for (j=1; j<genotype[index_mut].k; j++)
						{
							if (genotype[index_mut].mutation[j] != genotype[index_mut].mutation[j-1])
							{
								mut[mut_k] = genotype[index_mut].mutation[j];
								mut_k++;
							}
						}
						// overwrite:
						for (j=0; j<mut_k; j++)
							genotype[index_mut].mutation[j] = mut[j];
						genotype[index_mut].k = mut_k;

					}



				}

			}
		}


		//p_more = detect_mutations(&k_min, &k_max, &k_mean, &k_obs, 0.5);

		remove_zeros_fast();

		//fprintf(DB, "%d\t%d\t%d\t%f\t%d\t%d\t%f\n", gen, N, k_min, k_mean, k_max, k_obs, p_more);

		//count_N();

		if (gen % 10 == 0)
		{
			printf(".");
			fflush(stdout);
			//if(gen < gen_max)
			//fprintf(DB, ", ");
		}
		if ((gen % G ==0) | (gen == gen_max)){
			remove_duplicates();
			summary_to_tsv(DB, gen);
			printf("\b|");
			fflush(stdout);
		}
		if (gen == gen_max)
			printf("\n");
		//		{
		//			remove_duplicates();
		//			remove_zeros(1);
		//			summary(DB);

		//for (j=0; j<d; j++)
		//printf("%d\n", sum_obs[j]);

		//p_more = detect_mutations(&k_min, &k_max, &k_mean, &k_median, &k_obs, &k_plus, &k_minus, &prob_add, &homo, &div, 0.5);
		//fprintf(DB, "%d\t%d\t%d\t%d\t%g\t%g\t%d\t%g\t%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\n", run, gen, N_init, N, u, s, k_min, k_mean, k_median, k_max, k_obs, p_more, k_plus, k_minus, prob_add, homo, div);


		//}
		if(N_g == 0){
			printf("Error in g=%i\n",gen);
			exit(1);
		}


	}
	while(gen < gen_max);



	return 0;
}

int main(int argc, char **argv)
{

	// defaults:
	int    N = 1000000000;  // population size
	int N_init = 1;
	//  int    n_d = 100;  // drivers
	//  int    n_p = 100;  // passengers
	double u = 1e-7;  // mutation rate
	double v = -1.0; //mutation rate (passengers)
	double s = 1e-2;  // Selective advantage
	double s1 = 1.5*s; // Selective advantage of super drivers
	int    g = 1800;  // number of Generations
  int    X0 = 0; // generation at which population starts decreasing
  int    X1 = 0; // generation at which population stops decreasing
  double Xr = 1e-5; // rate at which population decreases
  int    Xm = 0; // whether cells are allowed to grow after decrease
	int    R = 1;  // number of simulations Runs
	int    d = 1000; // number of driver mutations
	int    d0 = 0; //number of passengers
	int    d1 = 0; //number of super drivers
	int    G = -1; //when to output
	double a = 2; // initial growth rate
  int    gd = 0; // should s follow be sample from a gamma dist.
	char *filestem;
	unsigned int seed = (unsigned) time(NULL);  // r, random seed
	int verbose = 0;

	int error_flag = 0;
	int f_flag = 0;

	int c = 0;
	while((c = getopt(argc, argv, "N:n:a:u:v:s:t:g:R:f:r:p:d:X:Y:Z:L:D:o:G:wh")) != EOF )
	{
		switch(c)
		{
		case 'N':
			if (atoi(optarg) > 0)
				N = atoi(optarg);
			else
				error_flag++;
			break;

		case 'n':
			if (atoi(optarg) > 0)
				N_init = atoi(optarg);
			else
				error_flag++;
			break;

		case 'a':
			if (atoi(optarg) > 0)
				a = atof(optarg);
			else
				error_flag++;
			break;

		case 'd':
			if (atoi(optarg) > 0)
				d = atoi(optarg);
			else
				error_flag++;
			break;

		case 'p':
			if (atoi(optarg) > 0)
				d0 = atoi(optarg);
			else
				error_flag++;
			break;

		case 'o':
			if (isNumeric(optarg) > 0)
				d1 = atoi(optarg);
			else
				error_flag++;
			break;

		case 'u':
			if (atof(optarg) > 0)
				u = atof(optarg);
			else
				error_flag++;
			break;

		case 'v':
			if (atof(optarg) > 0)
				v = atof(optarg);
			else
				error_flag++;
			break;

		case 's':
			if (atof(optarg) >= 0)
				s = atof(optarg);
			else
				error_flag++;
			break;

		case 't':
			if (atof(optarg) >= 0)
				s1 = atof(optarg);
			else
				error_flag++;
			break;

		case 'g':
			if (atoi(optarg) > 0)
				g = atoi(optarg);
			else
				error_flag++;
			break;

		case 'X':
			if (isNumeric(optarg) > 0)
				X0 = atoi(optarg);
			else
				error_flag++;
			break;

		case 'Y':
			if (isNumeric(optarg) > 0)
				X1 = atoi(optarg);
			else
				error_flag++;
			break;

    case 'Z':
			if (atof(optarg) > 0)
				Xr = atof(optarg);
			else
				error_flag++;
			break;

    case 'L':
			if (isNumeric(optarg) > 0)
				Xm = atoi(optarg);
			else
				error_flag++;
			break;

    case 'D':
      if (isNumeric(optarg) > 0)
        gd = atoi(optarg);
      else
        error_flag++;
      break;

		case 'R':
			if (atoi(optarg) >= 0)
				R = atoi(optarg);
			else
				error_flag++;
			break;

		case 'f':
			filestem = optarg;
			f_flag++;
			break;

		case 'r':
			seed = atoi(optarg);
			break;

		case 'G':
			if (atoi(optarg) >= 0)
				G = atoi(optarg);
			else
				error_flag++;
			break;

		case 'w':
			verbose = 1;
			break;

		case 'h':
			printf("usage: clonex [-N:n:u:v:s:t:g:R:f:r:p:d:o:wh]\n");
			printf("  N - Maximal population size (default = %d)\n", N);
			printf("  n - Initial population size (default = %d)\n", N_init);
			printf("  a - Initial growth rate (default = %g)\n", a);
			printf("  u - Mutation rate (default = %g)\n", u);
			printf("  v - Mutation rate passengers (default = u)\n");
			printf("  s - Selective advantage (default = %g)\n", s);
			printf("  t - Selective advantage of other drivers (default = %g)\n", s1);
			printf("  g - Number of generations (default = %d)\n", g);
			printf("  X - generation at which population starts decreasing (default = %d)\n", X0);
			printf("  Y - generation at which population stops decreasing (default = %d)\n", X1);
			printf("  Z - rate at which population decreases (default = %g)\n", Xr);
			printf("  L - Population increases after decreasing (default = %d (no))\n", Xm);
			printf("  D - Should s be sampled from a gamma distribution (default = %d (no))\n", gd);
			printf("  d - Number of drivers (default = %d)\n", d);
			printf("  p - Number of passengers (default = %d)\n", d0);
			printf("  o - Number of other drivers (default = %d)\n", d1);
			printf("  R - Replicates (default = %d)\n", R);
			printf("  r - Random seed (default = time)\n");
			printf("  f - File directory (Required! Make sure that the directory exists!)\n");
			printf("  G - Output every G generations (default = g)\n");
			printf("  h - This help\n");
			exit(0);

		default :
			exit(1);
		}
	}

	if(G == -1)
		G = g;

  if (X0 == 0)
    X0 = g + 1;

  if (X1 == 0)
    X1 = g + 1;

	if(v == -1.0)
		v=u;

	int genes = d + d0 + d1;

	if (error_flag || (! f_flag))
	{
		fprintf(stderr, "Error!\n");
		exit(1);
	}

	// random number generator:
	RNG = gsl_rng_alloc (gsl_rng_taus);  // global variable
	gsl_rng_set(RNG, seed);  // seed rng
	if (gsl_rng_max(RNG) - gsl_rng_min(RNG) < N)
	{
		printf("Population size N = %d too large for random number generator!\n", N);
		printf("RNG range = [%lu, %lu]\n", gsl_rng_min(RNG), gsl_rng_max(RNG));
		exit(1);
	}


	char summary_filename[255]; //, filename[255];
  int i;
	sprintf(summary_filename, "%s/sim.par", filestem);
	FILE *DB;
	if ((DB = fopen(summary_filename, "w")) == NULL)
	{
		fprintf(stderr, "Cannot open output file -- %s\n", summary_filename);
		exit(1);
	}
	
  for (i=0; i<argc; i++)
		fprintf(DB, "%s ", argv[i]);
	fprintf(DB, "\n");
	fclose(DB);

  float fitness_effects[d];
  if (gd == 1) {
    int i;
    double a;
    double b;
    double sample;
    a = 1;
    b = 1;
    for (i=1; i<=d; i++) {
      sample = gsl_ran_gamma(RNG,a,b) * s;
      fitness_effects[i-1] = sample;
    }
    FILE *file_handle;
    char fe_filename[255];

    sprintf(fe_filename, "%s/fitness_effects", filestem);
    if ((file_handle = fopen(fe_filename, "w")) == NULL) {
      fprintf(stderr, "Cannot open output file -- %s\n", fe_filename);
      exit(1);
    }

    sprintf(fe_filename, "%s/fitness_effects", filestem);
    file_handle = fopen(fe_filename,"w");
    fitness_distribution_to_tsv(d,fitness_effects,file_handle);
    fclose(file_handle);
  }

	int r;
	//#pragma omp parallel for private(DB, genotype, geno_prob, geno_count, N_g, fitness, sum_obs, k_abs_freq, k_rel_freq)
	for (r=0; r<R; r++)
	{
    printf("Output format is Gen, NMutClones, cloneID, mutID\n");
		printf("Sample %i/%i\n", r+1,R);
		sprintf(summary_filename, "%s/r%03d.csv", filestem, r+1);
		if ((DB = fopen(summary_filename, "w")) == NULL)
		{
			fprintf(stderr, "Cannot open output file -- %s\n", summary_filename);
			exit(1);
		}
		printf("%d\n", genes);
		simulate(DB, N_init, N, a, g, u, v, s, s1, r+1, genes, d0, d1, verbose, G, X0, X1, Xr, Xm,
             gd,fitness_effects);
		fflush(DB);
		fclose(DB);

	}



	return 0;
}



