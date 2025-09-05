/*********************************************************************
  RIsearch2   --   RNA-RNA interaction search

  Copyright (c) 2016-2025 by the contributors (see AUTHORS file)

  This file is part of RIsearch.

  RIsearch is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  RIsearch is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with RIsearch, see file COPYING.
  If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <sys/param.h>
#include <omp.h>
#include <errno.h>
#include <pcre.h>
#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <zlib.h>

#include "main.h"
#include "sa.h"
#include "dsm.h"
#include "seed.h"
#include "search.h"
#include "fasta.h"
#include "lists.h"
#include "weights.h"

#define min(a,b) \
    ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
      _a < _b ? _a : _b; })
#define max(a,b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
       _a > _b ? _a : _b; })

#define MATVERSION 1

Qt_type qt_type = RRNA;
char qt_mat_defaults[4][20] = {
	"t04.v3",
	"slh04.v3",
	"su95.v3",
	"su95.v3"
};

extern int print_debug;
int comp[6] = { 0, 4, 3, 2, 1, 5 };

int pair[6][6] = {
  /*  a, g, c, u */
  {0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 1, 0},		/* a */
  {0, 0, 0, 1, 1, 0},		/* g */
  {0, 0, 1, 0, 0, 0},		/* c */
  {0, 1, 1, 0, 0, 0},		/* u */
  {0, 0, 0, 0, 0, 0}
};

int pair_noGU[6][6] = {
  /*  a, g, c, u */
  {0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 1, 0},		/* a */
  {0, 0, 0, 1, 0, 0},		/* g */
  {0, 0, 1, 0, 0, 0},		/* c */
  {0, 1, 0, 0, 0, 0},		/* u */
  {0, 0, 0, 0, 0, 0}
};

const dsm_t *S;
const char *completed, *fasta_raw, *suffix, *output, *query;
const char *matname = NULL, *matname2 = NULL;
const char *matpath = MATPATH;
double min_energy = -20.0f;
float tempK[3] = {
	TEMPERATURE1, /* temperature to scale to */
	TEMPERATURE1, /* temperature valid for matname */
	TEMPERATURE2  /* temperature valid for matname2 */
};
char **name;
int threads = 0, seed_flag = 0, seed_orig[3] = { 6, 0, 0 }, extPen = 0;
int bands_flag = 0, bands = 0;

int show_alignment = 0;
int all_vs_all=1;
char *three_prime_opt , *five_prime_opt;
pcre *five_prime_match = NULL;
pcre *three_prime_match = NULL;
int noGUseed = 0;
double min_seed_energy_per_length = 0.0;
int seed_mismatch[3] = { 0, 0, 0};

int seed_mismatch_flag = 0, seed_threshold_flag = 0;
int max_ext_len = 20;
saidx64_t *idx_lengths;
saidx64_t num_idxs;
saidx64_t *sum_l;

int weight_stacks = 0;
char *weights_arr_name;
float *weights;


/* usage information */
void usage(const char *progname, const char *message)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "=============================== RIsearch2 v%s ===============================\n", PACKAGE_VERSION);
	fprintf(stderr, "================ Energy based RNA-RNA interaction predictions ================\n\n");
	fprintf(stderr, "Usage:         %s [ARGUMENTS]\n\n", progname);
	fprintf(stderr, "  -h,         --help\n");
	fprintf(stderr, "     show this message\n");
	fprintf(stderr, "--------------------------- SUFFIX ARRAY CREATION ----------------------------\n");
	fprintf(stderr, "  -c <FILE>,  --create=FILE (.fa or .fa.gz)\n");
	fprintf(stderr, "     create suffix array for target sequence(s) together with\n");
	fprintf(stderr, "     their reverse complements, FASTA format, use '-' for stdin\n");
	fprintf(stderr, "  -o <FILE>,  --output=FILE\n");
	fprintf(stderr, "     save created suffix array to given index file path \n");
	fprintf(stderr, "--------------------------- INTERACTION PREDICTION ---------------------------\n");
	fprintf(stderr, "  -q <FILE>, --query=FILE (.fa or .fa.gz)\n");
	fprintf(stderr, "     FASTA file for query sequence(s), use '-' for stdin\n");
	fprintf(stderr, "  -i <FILE>, --index=FILE\n");
	fprintf(stderr, "     pregenerated suffix array file for target sequence(s)\n");
	fprintf(stderr, "-------------------------------- SEED OPTIONS --------------------------------\n");
	fprintf(stderr, "  -s n:m/l,   --seed=n:m/l\n");
	fprintf(stderr, "     set seed length (-s l = length only; -s n:m = full interval;\n");
	fprintf(stderr, "     -s n:m/l = length in interval; default -s 6)\n");
	fprintf(stderr, "  --noGUseed     consider G-U wobble pairs as mismatch within the seed\n");
	fprintf(stderr, "     (only for locating seeds, energy model is not affected)\n");
	fprintf(stderr, "  -m c        --mismatch=c    Mismatched seeds with max num of mismatches (c)\n");
	fprintf(stderr, "  -m c:p      --mismatch=c:p    min consecutive matches at seed start/end (p)\n");
	fprintf(stderr, "  -m c:ps:pe  --mismatch=c:ps:pe or ps matches at the start & pe at the end\n");
	fprintf(stderr, "     These seeds will not overlap with perfect complementary seeds\n");
	fprintf(stderr, "     (default -m 0:0 means no mismatch                          )\n");
	fprintf(stderr, "     (when c>0, please set p>0 to avoid overlaps, p=c at default)\n");
	fprintf(stderr, "  -x <float>, --seed_energy=F\n");
	fprintf(stderr, "     set energy per length threshold to filters seeds (default=0)\n");
	fprintf(stderr, "------------------------------- BANDED RIsearch ------------------------------\n");
	fprintf(stderr, "  -b band,     --band=band\n");
	fprintf(stderr, "     Integer, size of the bands limiting the search.\n");
	fprintf(stderr, "     The min size is 1; use the seed option to avoid any bulge.\n");
	fprintf(stderr, "------------------------------- MODEL OPTIONS --------------------------------\n");
	fprintf(stderr, "  -d <int>, --penalty=dP\n");
	fprintf(stderr, "     per-nucleotide extension penalty given in dacal/mol\n");
	fprintf(stderr, "     (recommended: 30, default: 0)\n");
	fprintf(stderr, "  -l <int>,   --extension=L \n");
	fprintf(stderr, "     max extension length(L) on the seed (do DP for max this\n");
	fprintf(stderr, "     length up- and downstream of seed) (default L=20)\n");
	fprintf(stderr, "-------------------------- ENERGY MATRIX OPTIONS -----------------------------\n");
	fprintf(stderr, "  -N <str>, --qt_type\n");
	fprintf(stderr, "	     Query/Target type, one of RNA/RNA(default), DNA/DNA, RNA/DNA or \n");
	fprintf(stderr, "            DNA/RNA, makes the matrix default to t04.v3, slh04.v3 and su95.v3 or\n");
	fprintf(stderr, "            su95.v3, respectively.; \n");
	fprintf(stderr, "  -z mat, --matrix=mat\n");
	fprintf(stderr, "  -y mat2, --matrix2=mat2\n");
	fprintf(stderr, "     Only needed if you design your own energy matrices\n");
	fprintf(stderr, "  -K T1[,T2,T3], --temperature=T0[,T1,T2]\n");
	fprintf(stderr, "     Temperatures in Kelvin for scaling of energies. T0 is\n");
	fprintf(stderr, "     the temperature to be scaled to. T1 and T2 are only\n");
	fprintf(stderr, "     needed if you design your own energy parameters.\n");
	fprintf(stderr, "  -M PATH,   Path to directory holding the energy matrices\n");
	fprintf(stderr, "------------------------------- OUTPUT AND FILTERING -------------------------\n");
	fprintf(stderr, "  -e <float>, --energy=dG\n");
	fprintf(stderr, "     set deltaG energy threshold in kcal/mol to filter predictions\n");
	fprintf(stderr, "     (default=-20 for RIsearch2)\n");
	fprintf(stderr, "  -1, --one_vs_one\n");
	fprintf(stderr, "     only print results where query and target name matches\n");
	fprintf(stderr, "  -p1,-p     one line per hit, incl. interaction string, still header for each pair (query / target)\n");
	fprintf(stderr, "  -p2, --report_alignment=2   \n");
	fprintf(stderr, "     report in a simple format with CIGAR-like \n");
	fprintf(stderr, "     string for interaction structure\n");
	fprintf(stderr, "  -p3, --report_alignment=3   \n");
	fprintf(stderr, "     report in a simple format together with \n");
	fprintf(stderr, "     binding site (3'->5'), flanking 5'end (3'->5') and \n");
	fprintf(stderr, "     flanking 3'end (5'->3') sequences of the target\n");
	fprintf(stderr, "     (for post-processing of CRISPR off-target predictions)\n");
	fprintf(stderr, "  -p4, --report_alignment=3   \n");
	fprintf(stderr, "     report in the format target, start, strand, energy\n");
	fprintf(stderr, "--------------------------- CRISPR related options ---------------------------\n");
	fprintf(stderr, "  -3 <reg>,   --three_prime_match=PC\n");
	fprintf(stderr, "     report only predictions with a matching 3' PAM forward\n");
	fprintf(stderr, "     complement for cas9 this corresponds to PC=^(.cc|.uc|.cu)\n");
	fprintf(stderr, "     matching NGG/NAG/NGA 3' PAMs\n");
	fprintf(stderr, "     report mode 3/4 required\n");
	fprintf(stderr, "  -5 <reg>,   --five_prime_match=PC\n");
	fprintf(stderr, "     report only predictions with a matching 5' PAM reverse\n");
	fprintf(stderr, "     complement for cas12a this corresponds to PC=^([^a]aaa)\n");
	fprintf(stderr, "     matching TTTV 5'PAMs\n");
	fprintf(stderr, "     report mode 3/4 required\n");
	fprintf(stderr, "  -w arr,     --weights=arr\n");
	fprintf(stderr, "     CRISPR_gRNApPAM to weight gRNA-target interactions by\n");
	fprintf(stderr, "     CRISPR/Cas9 impact.\n");
	fprintf(stderr, "--------------------------------- MISC ------------------------------------\n");
	fprintf(stderr, "  -C dir,        Move results to dir (must exists) when they are completed.\n");
	fprintf(stderr, "  -t <int>,   --threads=N\n");
	fprintf(stderr, "     set maximum number of threads to use (default=1)\n");
	fprintf(stderr, "  --verbose      verbose output\n");
	fprintf(stderr, "\n");
	if (message) {
		fprintf(stderr, message);
	}
	exit(1);
}

void
str_rev (char *s)
{
  int s_len = strlen (s);
  int i;
  char temp;

  for (i = 0; i < s_len / 2; i++)
    {
      temp = s[i];
      s[i] = s[s_len - i - 1];
      s[s_len - i - 1] = temp;
    }
}

/* parse options */
void options(int argc, char *argv[])
{
	int c, option_index;
	char *ptr;

	/* list of all allowed options */
	static struct option long_options[] = {
		// Shared risearch1/2
		{"query", required_argument, 0, 'q'},
		{"one_vs_one", no_argument, 0, '1'},
		{"penalty", required_argument, 0, 'd'},
		{"matrix", required_argument, 0, 'z'},
		{"matrix2", required_argument, 0, 'y'},
		{"temperature", required_argument, 0, 'K'},
		{"energy", required_argument, 0, 'e'},
		{"weights", required_argument, 0, 'w'},
		{"qt_type", required_argument, 0, 'N'},
		{"verbose", no_argument, 0, 'v'},
		// risearch2 only
		{"three_prime_match", required_argument, 0, '3'},
		{"five_prime_match", required_argument, 0, '5'},
		{"bands", required_argument, 0, 'b'},
		{"create", required_argument, 0, 'c'},
		{"help", no_argument, 0, 'h'},
		{"index", required_argument, 0, 'i'},
		{"extension", required_argument, 0, 'l'},
		{"mismatch", required_argument, 0, 'm'},
		{"matpath", required_argument, 0, 'M'},
		{"output", required_argument, 0, 'o'},
		{"report_alignment", optional_argument, 0, 'p'},
		{"seed", required_argument, 0, 's'},
		{"threads", required_argument, 0, 't'},
		{"noGUseed", no_argument, 0, 'U'},
		{"seed_energy", required_argument, 0, 'x'},
		{0, 0, 0, 0}
	};

	/* parse all options */
	while ((c = getopt_long(argc, argv, "13:5:b:c:C:d:e:hi:K:l:m:M:N:o:p:q:s:t:Uvw:x:y:z:", long_options, &option_index)) != -1) {
		switch (c) {
		case 0:
			break;
		case '1':
			all_vs_all = 0;
			break;
		case 'v':
			verbose = 1;
			break;
		case 'q':
			query = optarg;
			debug("opt: query=%s\n", query);
			break;
		case 'c':
			fasta_raw = optarg;
			debug("opt: fasta_raw=%s\n", fasta_raw);
			break;
		case 'd':
			extPen = atoi(optarg) * SCALEFACTOR / 100;
			debug("opt: extPen=%d\n", extPen);
			break;
		case 'M':
			matpath = strdup(optarg);
			debug("opt: matpath=%s\n", matpath);
			break;
		//case 'm':
		//	fprintf(stderr, "Parameter -m is deprecated, please use -z.\n");
		//	matname = optarg;
		//	break;
		case '3':
			{
				const char *error;
				int error_offset;
				three_prime_opt = optarg;
				three_prime_match = pcre_compile(optarg, 0, &error, &error_offset, NULL);
				if (three_prime_match == 0) {
					usage(argv[0], "Could not compile anti-PAM regex\n");
				}
				debug("opt: compiled_prime_match_from=%s\n", three_prime_opt);
			}
			break;
		case '5':
			{
				const char *error;
				int error_offset;
				five_prime_opt = optarg;
				five_prime_match = pcre_compile(optarg, 0, &error, &error_offset, NULL);
				if (five_prime_match == 0) {
					usage(argv[0], "Could not compile anti-PAM regex\n");
				}
				debug("opt: compiled_five_prime_match_from=%s\n", five_prime_opt);
			}
			break;
		case 'C':
			completed = optarg;
			debug("opt: completed=%s\n", completed);
			break;
		case 'b':
			// bands
			bands_flag = 1;
			bands = atoi(optarg);
			if (bands < 1) {
				fprintf(stderr, "Invalid bands parameter. Please use a positive integer > 0.\n");
				abort();
			}
			debug("opt: bands param =(%d)\n", bands);
			break;
		case 'e':
			min_energy = strtof(optarg, &ptr);
			if (*ptr) {
				// Did not read until \0, came across something that could not be converted to float
				usage(argv[0], "Non-numeric value passed with option -e");
			}
			debug("opt: min_energy=%f\n", min_energy);
			break;
		case 'h':
			usage(argv[0], NULL);
			break;
		case 'i':
			suffix = optarg;
			debug("opt: suffix=%s\n", suffix);
			break;
		case 'y':
			matname2 = optarg;
			debug("opt: matrix2=%s\n", matname2);
			break;
		case 'z':
			matname = optarg;
			debug("opt: matrix=%s\n", matname);
			break;
		case 'N':
			if (strcmp(optarg, "RNA/RNA") == 0) {
				qt_type = RRNA;
			} else if (strcmp(optarg, "DNA/DNA") == 0) {
				qt_type = DDNA;
			} else if (strcmp(optarg, "RNA/DNA") == 0) {
				qt_type = RDNA;
			} else if (strcmp(optarg, "DNA/RNA") == 0) {
				qt_type = DRNA;
			} else {
				fprintf(stderr, "Parameter -N should be one of RNA/RNA, DNA/DNA, RNA/DNA and DNA/RNA.\n");
				exit(1);
			}
			break;
		case 'K':
			{
				char *next;
				char *tmp = optarg;
				int i = 0;
				debug("opt: temperatures=%s\n", optarg);
				next = strtok(tmp, ",");
				while (next && i < 3) {
					debug("opt: next=%s\n", next);
					tempK[i] = atof(next);
					next = strtok(NULL, ",");
					i++;
				}
				debug("opt: temperatures=(%f, %f, %f)\n", tempK[0], tempK[1], tempK[2]);
			}
			break;
		case 'w':
			weight_stacks = 1;
			weights_arr_name = optarg;
			break;
		case 'l':
			max_ext_len = MAX(0, atoi(optarg));
			debug("opt: max_ext_len=%d\n", max_ext_len);
			break;
		case 'm':
			/* "-m c:s:e" c:allowed number of mismatches. s/e is the number of consecutive matches required at the beginning/end of the seed. */
			seed_mismatch_flag = 1;
			seed_mismatch[0] = atoi(optarg);
			optarg = strchr(optarg, ':');
			if (!optarg) {
				seed_mismatch[2] = seed_mismatch[0];
				seed_mismatch[1] = seed_mismatch[0];
				debug("opt: seed_mismatch=(%d,%d,%d)\n", seed_mismatch[0], seed_mismatch[1], seed_mismatch[2]);
				break;
			}
			optarg++;
			seed_mismatch[1] = atoi(optarg);
			optarg = strchr(optarg, ':');
			if (!optarg) {
				seed_mismatch[2] = seed_mismatch[1];
				debug("opt: seed_mismatch=(%d,%d,%d)\n", seed_mismatch[0], seed_mismatch[1], seed_mismatch[2]);
				break;
			}
			optarg++;
			seed_mismatch[2] = atoi(optarg);
			debug("opt: seed_mismatch=(%d,%d,%d)\n", seed_mismatch[0], seed_mismatch[1], seed_mismatch[2]);
			break;
		case 'o':
			output = optarg;
			debug("opt: output=%s\n", output);
			break;
		case 'p':
			debug("opt: show_alignment=%s\n", optarg);
			debug("opt: show_alignment=%d\n", atoi(optarg));
			if (optarg != NULL)
				show_alignment = atoi(optarg);
			else
				show_alignment = 1;
			debug("opt: show_alignment=%d\n", show_alignment);
			break;
		case 's':
			/* -s l (length), -s m:n (interval), -s m:n/l (interval+length) */
			seed_orig[0] = atoi(optarg);
			optarg = strchr(optarg, ':');
			if (optarg) {
				optarg++;
				seed_orig[1] = atoi(optarg);
				seed_flag = 1;
				optarg = strchr(optarg, '/');
				if (optarg) {
					optarg++;
					seed_orig[2] = atoi(optarg);
				}
			}
			debug("opt: seed_orig=(%d,%d,%d)\n", seed_orig[0], seed_orig[1], seed_orig[2]);
			break;
		case 't':
			//threads = MAX(1, atoi(optarg));
			if ((threads = MAX(1, atoi(optarg))) > 0) {
				omp_set_dynamic(0);
				omp_set_num_threads(threads);
			}
			debug("opt: threads=%d\n", threads);
			break;
		case 'U':
			noGUseed = 1;
			break;
		case 'x':
			seed_threshold_flag = 1;
			min_seed_energy_per_length = atof(optarg);
			debug("opt: min_seed_energy_per_length=%f\n", min_seed_energy_per_length);
			break;
		case '?':
			fprintf(stderr, "Parameter -%c: unknown option or missing argument.\n", optopt);
			break;
		default:
			abort();
		}

	}
	if (!matname) {
		matname = qt_mat_defaults[qt_type];
		debug("default: matrix=%s\n", matname);
	}
}


void getWeights(char *weigths_arr_name)
{
    if (!strcmp(weigths_arr_name, "CRISPR_gRNApPAM")){
        weights = (float *) malloc(size_wC23_3p_5p*sizeof(float));
        for(int i = 0; i < size_wC23_3p_5p; ++i){
            weights[i] = wC23_3p_5p[i];
        }
    }
    else{
        fprintf(stderr, "Undefined name for array of weights\n");
		exit(EXIT_FAILURE);
    }
}

/** @brief Reads queries from fasta file and sets it up
 *  
 *  @param ffp The fasta (opened) from where to read
 *  @param max_to_read The maximum number of queries that are read
 *  @param n_read The address to hold the number of sequences read (#entries in returned)
 *  
 *  @return queries or 0 if none were read
 */
query_t *read_queries(fasta_t *ffp, int max_to_read, int *n_read)
{
	query_t *queries = calloc(max_to_read, sizeof(query_t));
	if (queries == NULL) {
		fprintf(stderr, "Failed allocating queries\n");
		exit(EXIT_FAILURE);
	}
	int n_seq = 0;

	while(fasta_read(ffp, &queries[n_seq].seq, &queries[n_seq].name, &queries[n_seq].length)) {
		query_t *query = &queries[n_seq];
		memset(query->seed, 0, 3 * sizeof(int));
		//fprintf(stderr, "sequence: %s, %ld\n" ,query->seq, query->length );
		seeds_setup(query->seed, &query->length);
		if (query->seed[2] > query->length) {
			//printf("#%s\n", query->name);
			continue;
		}
		sa_create_partial_reverse(query->seq, query->seed[0]-1, query->seed[1]-1, query->seed[2], &query->qsa, &query->length);
		n_seq += 1;
		*n_read = n_seq;
		if (n_seq == max_to_read) {
			return queries;
		}
	}

	if (n_seq == 0) {
		free(queries);
		return 0;
	}

	return queries;
}

typedef struct {
	int tid;
	saidx64_t *sa, n;
	query_t *queries;
	int n_queries;
} thread_aux_t;

static void *match_worker(void *data)
{
	int i;
	thread_aux_t *d = (thread_aux_t *)data;
	sa_interval_list_t *intervals;
	int tid;

#pragma omp parallel private(tid, i)
	{
		int nthreads = omp_get_num_threads();
		char tmpname[100];
		char tmpname2[120];

		for (i = 0; i < d->n_queries; i++) {
			tid = omp_get_thread_num();

			//fprintf(stderr, "i: %d tid: %d\n", i, tid);
			if (i % nthreads == tid) {
				query_t *query = &d->queries[i];

				{
					debug("processing query %s \n", query->name);
					// should probably get a basename as parameter and combine to: path/foo_{query->name} or so 
					// possibly also warn if file is present before open "w" and overwrite what was there...
					assert(strlen(query->name) < 80);
					snprintf(tmpname, 100, "risearch_%s.out.gz", query->name);
					query->out = gzopen(tmpname, "wb");
					if (query->out == NULL) {
						printf("gzopen file %s failed, errno = %d\n", tmpname, errno);
					}
					//if (i % 10 == 0)
					//fprintf(stderr, "starting: %d tid: %d len: %d\n", i, tid, query->length);
					else if (query->length >= query->seed[2]) {
						sa_parallel_match_neg(query->qsa, 0, query->length, d->sa, 0, d->n, 0, query->seed[2], 0, 0, &query->intervals);
						if (!weight_stacks){
						    sa_evaluate_interval(query->intervals, query, d->sa, '-');
						}
						else{
						    sa_evaluate_interval_weighted(query->intervals, query, d->sa, '-');
						}
						gzclose(query->out);
						if (completed) {
							snprintf(tmpname2, 120, "completed/risearch_%s.out.gz", query->name);
							rename(tmpname, tmpname2);
						}
					}
					//fprintf(stderr, "ending: %d tid: %d\n", i, tid);
				}

				//fprintf(stderr, "**i %d tid: %d nthreads: %d\n", i, tid, nthreads);
				//fprintf(stderr, "**current entry: %d name: %s\n", i, query->name);

				free(query->qsa);
				free(query->seq);
				free(query->name);
				intervals = query->intervals;

#pragma omp critical
				{
					while(intervals) {
						void *p = (void *)intervals;
						intervals = intervals->next;
						free(p);
					}
				}
			}
		}
	}
	return 0;
}

/** @brief entry point....
 * 
 *  parses options
 *  gets matrix
 *  given -c and -o : creates suffix from fasta
 *  given -i and -q : reads suffix, prepares queries, starts threads...
 */
int main(int argc, char *argv[])
{
	saidx64_t *sa, n;
	fasta_t *ffp;
	int i;
	dsm_t dsm;

	options(argc, argv);

	debug("initializing\n");

	if ((show_alignment != 3 && show_alignment != 4) && (three_prime_match || five_prime_match)) {
		usage(argv[0], "Pam matching only works with mode 3/4\n");
		return 0;
	}

	if(fasta_raw && output) {
		debug("creating suffix array from \"%s\"\n", fasta_raw);
		sa_fasta_to_file(fasta_raw, output);
	} else if(suffix && query) {
		getMat(qt_type, MATVERSION, matname, matname2, matpath, tempK, &dsm[0][0][0][0]);
		if (weight_stacks == 1){
		    getWeights(weights_arr_name);
		}
		S = (const dsm_t*) &dsm;
		debug("reading suffix array from \"%s\"\n", suffix);
		sa_read(suffix, &sa, &n, &idx_lengths, &sum_l, &num_idxs, &name);
		// n is twice the number of nucleotides in input, because of reverse complement
		// k is twice the number of sequences in input, because of reverse complement
		debug("N=%li, K=%li\n", n, num_idxs);
		//sa_length_to_index(&l, k);

		ffp = fasta_open(query);
		if (ffp == NULL) {
			fprintf(stderr, "Query file %s is not readable\n", query);
			exit(EXIT_FAILURE);
		}
		query_t *queries;
		int n_queries = 0;
		while ((queries = read_queries(ffp, 0x40000, &n_queries)) != 0) {
			thread_aux_t data;
			data.n = n;
			data.sa = sa;
			data.queries = queries;
			data.n_queries = n_queries;

			match_worker(&data);

			free(queries);
		}

		fasta_close(ffp);
		free(sa);
		free(idx_lengths);
		free(sum_l);
		for(i = 0; i < num_idxs; ++i)
			free(name[i]);
		free(name);
	} else {
		usage(argv[0], "Missing options -q, -i, -o or -c");
	}

	debug("exiting\n");

	return 0;
}

