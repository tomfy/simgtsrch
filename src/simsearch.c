// C version of 'k-mer' search for pairs of similar genotype sets.
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h> // needed for getopt
#include <assert.h>

#include "gtset.h"
#include "various.h"

//***********************************************************************************************
// **************  typedefs  ********************************************************************

typedef struct{ // one of these for each chunk
  long capacity;
  long size;
  Vlong** a; // array of Vlong*; indices are patterns, values are Vlong*s containing indices of accessions with that pattern
} Pattern_ids;

typedef struct{
  long capacity;
  long size; // number of chunks
  long chunk_size; // number of chunks used.
  long n_patterns; // 3^chunk_size
  Pattern_ids** a; // array of Pattern_ids, one element for each chunk
} Chunk_pattern_ids;

typedef struct{
  long query_index;
  long match_index;
  double usable_chunks; // estimate or actual count
  long n_matching_chunks;
  double est_agmr;
  double agmr;
  double hgmr; 
} Mci; // 'Mci = Matching chunk info'

typedef struct{
  long capacity;
  long size;
  Mci** a;
} Vmci;

int  do_checks_flag = 0;


// *********************** function declarations ************************************************

char* ipat_to_strpat(long len, long ipat); // unused
long strpat_to_ipat(long len, char* strpat); // unused
double agmr(Accession* gts1, Accession* gts2, double* hgmr);

// *****  Mci  ********
Mci* construct_mci(long qidx, long midx, double n_usable_chunks, long n_matching_chunks,
		   double est_agmr, double agmr, double hgmr);
// *****  Vmci  *********************************************************************************
Vmci* construct_vmci(long init_size);
void add_mci_to_vmci(Vmci* the_vmci, Mci* the_mci);
void free_vmci(Vmci* the_vmci);

// *****  Pattern_ids; indices are patterns; elements are Vlong* of accids having that pattern.
Pattern_ids* construct_pattern_ids(long n_patterns);
void free_pattern_ids(Pattern_ids*);

// *****  Chunk_pattern_ids; indices are chunk numbers; elements are Pattern_ids*
Chunk_pattern_ids* construct_chunk_pattern_ids(long n_chunks, long chunk_size);
// Vlong** get_all_match_counts(long n_accessions, Chunk_pattern_ids* the_cpi);
void populate_chunk_pattern_ids_from_vaccession(Vaccession* the_accessions, Chunk_pattern_ids* the_cpi);
void print_chunk_pattern_ids(Chunk_pattern_ids* the_cpi, FILE* ostream);
void free_chunk_pattern_ids(Chunk_pattern_ids* the_cpi);

// *****  Gts and Chunk_pattern_ids  ***********
Vlong* find_chunk_match_counts(Accession* the_gts, Chunk_pattern_ids* the_cpi);
Vmci** find_matches(long n_ref_accessions, Vaccession* the_accessions, Chunk_pattern_ids* the_cpi, double max_est_agmr);
// Vmci** find_matches_alt(Vgts* the_accessions, Vlong** match_counts, long n_chunks, long chunk_size);
long print_results(Vaccession* the_accessions, Vmci** query_vmcis, FILE* ostream);

// *************************  end of declarations  **********************************************


// **********************************************************************************************
// ****************************** main **********************************************************
// **********************************************************************************************

int
main(int argc, char *argv[])
{
  double start0 = hi_res_time();
  long n_chunks = 100000; // default number of chunks (large number -> use all markers)
  long chunk_size = 8; // default number of genotype in each chunk 
  long max_number_of_accessions = -1;
  //  double fraction_to_analyze = 1;
  unsigned rand_seed = (unsigned)time(0);
  double max_md_factor = 100.0; // multiplies n_markers/chunk_size to give max number of md gts allowed in gt set of an accession.
  long max_md_gts; // accessions with > this number of missing data gts are omitted from analysis.
  // long min_usable_chunks = 5;
  double max_est_agmr = 0.2;
  //  long genotype_file_type = GENOTYPES; // DOSAGES;
  double delta = 0.05;
 
  double max_marker_missing_data_fraction = 0.2;

  char* rparam_buf;
  size_t rparam_len;
  FILE* rparam_stream = open_memstream(&rparam_buf, &rparam_len);

  
  fprintf(stderr, "# command line:  ");
  for(int i=0; i<argc; i++){
    fprintf(stderr, "%s  ", argv[i]);
  }fprintf(stderr, "\n");
  
  // ***** process command line *****
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <file> options\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  char* input_filename = NULL;
  FILE *in_stream = NULL;
   char* reference_set_filename = NULL;
   FILE *ref_in_stream = NULL;
  char* output_filename = NULL;
  FILE* out_stream = stdout;
    
  int c;
  while((c = getopt(argc, argv, "i:o:p:n:k:e:s:x:r:")) != -1){
    // i: input file name (required).
    // r: reference set file name.
    // o: output file name. Default: output goes to stdout.
    // p not implemented   // p: keep each accession with probability p. (For test purposes) Default = 1
    // n: number of chunks to use. Default: use each marker ~once.
    // k: chunk size (number of markers per chunk). Default: 8
    // e: max estimated agmr. Default: 0.2 (Calculate agmr only if quick est. is < this value.)
    // s: random number seed. Default: get seed from clock.
    // x: max missing data factor (default: 1 -> max_md_gts = n_markers/chunk_size )
     
    switch(c){
    case 'i':
      input_filename = optarg;
      in_stream = fopen(input_filename, "r");
      if(in_stream == NULL){
	fprintf(stderr, "Failed to open %s for reading.\n", input_filename);
	exit(EXIT_FAILURE);
      }
      fclose(in_stream);
      break;
       case 'r':
      reference_set_filename = optarg;
      ref_in_stream = fopen(reference_set_filename, "r");
      if(ref_in_stream == NULL){
	fprintf(stderr, "Failed to open %s for reading.\n", reference_set_filename);
	exit(EXIT_FAILURE);
      }
      fclose(ref_in_stream);
      break;
    case 'o':
      output_filename = optarg;
      if(output_filename != NULL){
	out_stream = fopen(output_filename, "w");
	if(out_stream == NULL){
	  fprintf(stderr, "Failed to open %s for writing.\n", output_filename);
	  exit(EXIT_FAILURE);
	}
      }
      break;
    /* case 'p': // keep each accession with probability p (for testing with random smaller data set) */
    /*   fprintf(stderr, "%s\n", optarg); */
    /*   fraction_to_analyze = (double)atof(optarg); */
    /*   if(fraction_to_analyze <= 0  || fraction_to_analyze > 1.0){ */
    /* 	fprintf(stderr, "# fraction_to_analyze specified as %6.3f. Out of valid range; exiting.\n", fraction_to_analyze); */
    /* 	exit(EXIT_FAILURE); */
    /*   } */
    /*   break; */
    case 'n': 
      n_chunks = (long)atoi(optarg);
      if(n_chunks <= 0){
	fprintf(stderr, "option n (n_chunks) requires an integer argument > 0\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'k':
      chunk_size = (long)atoi(optarg);
      if(chunk_size <= 0){
	fprintf(stderr, "option k (chunk_size) requires a integer argument > 0\n");
	exit(EXIT_FAILURE);
      }
      break; 
    case 'e': 
      max_est_agmr = (double)atof(optarg);
      if(max_est_agmr <= 0  || max_est_agmr > 1.0){
	fprintf(stderr, "option e (max_est_agmr) requires a numerical argument 0<x<=1 \n");
	exit(EXIT_FAILURE);
      }
      break;
    case 's': // random number generator seed
      rand_seed = (unsigned)atoi(optarg);
      if(rand_seed <= 0){
	fprintf(stderr, "option s (rng seed) requires a integer argument > 0\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'x': 
      max_md_factor = (double)atof(optarg);
      if(max_md_factor <= 0){
	fprintf(stderr, "option x (max_md_factor) requires an real argument > 0\n");
	exit(EXIT_FAILURE);
      }
      break;
    case '?':
      fprintf(stderr, "? case in command line processing switch.\n");
      if ((optopt == 'i') || (optopt == 'n') || (optopt == 'k') ||
	  (optopt == 'o') || (optopt == 'p') || (optopt == 'e') ||
	  (optopt == 's') || (optopt == 'h') || (optopt == 'z') )
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (isprint (optopt))
	fprintf (stderr, "Unknown option `-%c'.\n", optopt);
      else
	fprintf (stderr, "Unknown option character: %d\n", optopt);
      exit(EXIT_FAILURE);
    default:
      fprintf(stderr, "default case (abort)\n");
      exit(EXIT_FAILURE);
    } // end of switch block
  } // end of loop over c.l. arguments

  if(optind < argc){
    fprintf(stderr, "Non-option arguments. Bye.\n");
    exit(EXIT_FAILURE);
  }  if(input_filename == NULL){
    perror("must specify input filename: -i <filename>");
    exit(EXIT_FAILURE);
  }
  
  srand(rand_seed);

  fprintf(rparam_stream, "# input file: %s  output to: %s\n", input_filename, (output_filename == NULL)? "stdout" : output_filename);
  fprintf(stderr, "# input file: %s \n", input_filename);	  
 
  // *****  done processing command line  *****

  long n_ref_accessions = 0;
  long n_accessions = 0;
  long n_markers = 0;

  GenotypesSet* the_gtsset = construct_empty_genotypesset(delta, max_marker_missing_data_fraction);
  Vaccession* the_accessions; // = construct_vaccession(INIT_VACC_CAPACITY);
  if(reference_set_filename != NULL){ // load the reference set, if one was specified.
      add_accessions_to_genotypesset_from_file(reference_set_filename, the_gtsset);
      n_ref_accessions = the_gtsset->accessions->size;
      the_gtsset->n_ref_accessions = n_ref_accessions;
    }
  add_accessions_to_genotypesset_from_file(input_filename, the_gtsset); // load the new set of accessions
    clean_genotypesset(the_gtsset); 
    the_accessions = the_gtsset->accessions;
   
  n_markers = the_gtsset->n_markers;
   if(n_chunks*chunk_size > n_markers){
    n_chunks = n_markers/chunk_size;
  }

  fprintf(rparam_stream, "# n_chunks: %ld  chunk_size: %ld  max_est_agmr: %5.3lf rng seed: %u\n", n_chunks, chunk_size, max_est_agmr, rand_seed);
  fclose(rparam_stream);
  fprintf(stderr, "%s", rparam_buf);
  fprintf(out_stream, "%s", rparam_buf);

  // *****  done reading and storing input  **********
  
  double start = hi_res_time();
  fprintf(stderr, "# n_markers: %ld\n", n_markers);
  Vlong* marker_indices = construct_vlong_whole_numbers(n_markers);
  shuffle_vlong(marker_indices);
  set_vaccession_chunk_patterns(the_accessions, marker_indices, n_chunks, chunk_size);
  Chunk_pattern_ids* the_cpi = construct_chunk_pattern_ids(n_chunks, chunk_size);
   populate_chunk_pattern_ids_from_vaccession(the_accessions, the_cpi);
   
  fprintf(stderr, "# time to construct chunk_pattern_ids structure: %12.6f\n", hi_res_time() - start);
  
  start = hi_res_time(); 
  Vmci** query_vmcis = find_matches(n_ref_accessions, the_accessions, the_cpi, max_est_agmr);
  long true_agmr_count = print_results(the_accessions, query_vmcis, out_stream);
  fprintf(stderr, "# time to find candidate matches and %ld true agmrs: %9.3f\n", true_agmr_count, hi_res_time() - start);
  
  // *****  clean up  *****
  for(long i=0; i< the_accessions->size; i++){
    free_vmci(query_vmcis[i]);
  }
  free(query_vmcis); 
  fclose(out_stream);
  free_vlong(marker_indices);
  free_chunk_pattern_ids(the_cpi);
  free_vaccession(the_accessions);
  fprintf(stderr, "# total simsearch run time: %9.3f\n", hi_res_time() - start0);
  exit(EXIT_SUCCESS);
}

// **********************************************************************************************
// **********************  end of main  *********************************************************
// **********************************************************************************************


// *******************  function definitions  ***************************************************
// *****  Vaccession  ***********************************************************


void populate_chunk_pattern_ids_from_vaccession(Vaccession* the_accessions, Chunk_pattern_ids* the_cpi){
  long n_patterns = the_cpi->n_patterns;
 
  for(long i_gts=0; i_gts<the_accessions->size; i_gts++){
    Accession* the_gts = the_accessions->a[i_gts];
    if(DO_ASSERT) assert(i_gts == the_gts->index);
    Vlong* the_chunk_patterns = the_gts->chunk_patterns; // the gt patterns (longs) occurring in each chunk of this gts 
    long mdcount = 0;
    for(long i=0; i<the_chunk_patterns->size; i++){
      if(the_chunk_patterns->a[i] == n_patterns){ mdcount++; }
    }
    for(long i_chunk=0; i_chunk<the_chunk_patterns->size; i_chunk++){
      long the_pat = the_chunk_patterns->a[i_chunk];
      if(DO_ASSERT) assert(the_pat >= 0);

      if(the_cpi->a[i_chunk]->a[the_pat] == NULL){
	the_cpi->a[i_chunk]->a[the_pat] = construct_vlong(1);
      }
      Vlong* the_accidxs = the_cpi->a[i_chunk]->a[the_pat];	
      add_long_to_vlong(the_accidxs, the_gts->index);
    }
  }

  long total_mdchunk_count = 0;
  fprintf(stderr, "# the_cpi->size: %ld\n", the_cpi->size);
  for(long i=0; i<the_cpi->size; i++){
    long chunk_md_count = (the_cpi->a[i]->a[n_patterns] == NULL)?
      0 :  // there are no accessions having missing data for this chunk
      the_cpi->a[i]->a[n_patterns]->size;
    total_mdchunk_count += chunk_md_count;
  } 
}

// *****  Pattern_ids; indices are patterns; elements are Vlong* of accidxs having that pattern.

Pattern_ids* construct_pattern_ids(long n_patterns){
  Pattern_ids* pat_ids = (Pattern_ids*)malloc(1*sizeof(Pattern_ids));
  pat_ids->capacity = n_patterns+1; // 0..n_patterns-1 are the indices of the n_patterns (=3^k) good patterns, and index n_pattern is for the missing data case.
  pat_ids->size = n_patterns+1;
  pat_ids->a = (Vlong**)malloc((n_patterns+1)*sizeof(Vlong*));
  for(long ipat=0; ipat< pat_ids->size; ipat++){
    pat_ids->a[ipat] = NULL; // construct_vlong(1); // waste of memory? set to NULL until needed?
  }
  return pat_ids;
}

void free_pattern_ids(Pattern_ids* pat_ids){
    if(pat_ids == NULL) return;
  for(long i=0; i<pat_ids->size; i++){
    if(pat_ids->a[i] != NULL) free_vlong(pat_ids->a[i]);
  }
  free(pat_ids->a);
  free(pat_ids);
}

// *****  Chunk_pattern_ids; indices are chunk numbers; elements are Pattern_ids*

Chunk_pattern_ids* construct_chunk_pattern_ids(long n_chunks, long chunk_size){ // needed size is known at construct time, so one param for both cap and size
  Chunk_pattern_ids* chunk_pat_ids = (Chunk_pattern_ids*)malloc(1*sizeof(Chunk_pattern_ids));
  chunk_pat_ids->capacity = n_chunks;
  chunk_pat_ids->size = n_chunks;
  chunk_pat_ids->chunk_size = chunk_size;
  long n_patterns = int_power(3, chunk_size);
  chunk_pat_ids->n_patterns = n_patterns;
  chunk_pat_ids->a = (Pattern_ids**)malloc(n_chunks*sizeof(Pattern_ids*));
  for(int i=0; i< chunk_pat_ids->size; i++){
    chunk_pat_ids->a[i] = construct_pattern_ids(n_patterns);
  }
  return chunk_pat_ids;
}

void free_chunk_pattern_ids(Chunk_pattern_ids* the_cpi){
    if(the_cpi == NULL) return;
  for(long i=0; i< the_cpi->size; i++){
    free_pattern_ids(the_cpi->a[i]);
  }
  free(the_cpi->a);
  free(the_cpi);
}

void print_chunk_pattern_ids(Chunk_pattern_ids* the_cpi, FILE* ostream){
  for(long i_chunk=0; i_chunk<the_cpi->size; i_chunk++){
    fprintf(ostream, "i_chunk: %ld\n", i_chunk);
    Pattern_ids* the_pi = the_cpi->a[i_chunk];
    for(long i_pat=0; i_pat<the_pi->size; i_pat++){
      fprintf(ostream, "  i_pat: %ld \n", i_pat);
      Vlong* the_idxs = the_pi->a[i_pat];
      if(the_idxs->size > 0){
	fprintf(ostream, "     matches: ");
	for(long ii=0; ii<the_idxs->size; ii++){
	  fprintf(ostream, "%ld  ", the_idxs->a[ii]);
	}
	fprintf(ostream, "\n");
      }
    }
  }
}

// *****  Accession and Chunk_pattern_ids  ***********

Vlong* find_chunk_match_counts(Accession* the_gts, Chunk_pattern_ids* the_cpi){ //, long n_accessions){ //, Vlong** accidx_hmatchcounts){
  // fprintf(stderr, "top of find_chunk_match_counts\n");
  long n_patterns = the_cpi->n_patterns;
  Vlong* chunk_pats = the_gts->chunk_patterns;
  Vlong* accidx_matchcounts = construct_vlong_zeroes(the_gts->index); //   n_accessions);
 
  for(long i_chunk=0; i_chunk < chunk_pats->size; i_chunk++){
    long the_pat = chunk_pats->a[i_chunk];  
    Vlong* chunk_match_idxs = the_cpi->a[i_chunk]->a[the_pat]; // array of indices of the matches to this chunk & pat
    // (patterns 0..n_patterns-1 are good, n_patterns=3^chunk_size is the pattern for missing data )
    if(the_pat == n_patterns){ // missing data in this chunk 
    }else{ // the_pat = 0..n_patterns-1 (good data)   
      // just get the counts for matches with index < index of the_gts
      // since those are the only ones used in find_matches. (slightly faster)
      for(long i=0; i<chunk_match_idxs->size; i++){	  
	long accidx = chunk_match_idxs->a[i]; // index of one of the accessions matching on this chunk
	if(accidx >= the_gts->index) break;
	accidx_matchcounts->a[accidx]++; // accidx here is < the_gts->index
      }
    }
  }
  // fprintf(stderr, "bottom of find_chunk_match_counts\n");
  return accidx_matchcounts; 
}

// ***** Mci  *****

Mci* construct_mci(long qidx, long midx, double usable_chunks, long n_matching_chunks,
		   // double est_matching_chunk_fraction, double matching_chunk_fraction){
		   double est_agmr, double agmr, double hgmr){
  Mci* the_mci = (Mci*)calloc(1,sizeof(Mci));
  the_mci->query_index = qidx;
  the_mci->match_index = midx;
  the_mci->usable_chunks = usable_chunks;
  the_mci->n_matching_chunks = n_matching_chunks;
  the_mci->est_agmr = est_agmr;
  the_mci->agmr = agmr;
  the_mci->hgmr = hgmr;
  
  return the_mci;
}



// *****  Vmci  *********************************************************************************

Vmci* construct_vmci(long init_cap){
  Vmci* the_vmci = (Vmci*)malloc(1*sizeof(Vmci));
  the_vmci->capacity = init_cap;
  the_vmci->size = 0;
  the_vmci->a = (Mci**)malloc(init_cap*sizeof(Mci*));
  return the_vmci;
}

void add_mci_to_vmci(Vmci* the_vmci, Mci* mci){
  long cap = the_vmci->capacity;
  long n = the_vmci->size;
  // if necessary, resize w realloc
  if(n == cap){
    cap *= 2;
    the_vmci->a = (Mci**)realloc(the_vmci->a, cap*sizeof(Mci*));
    the_vmci->capacity = cap;
  }
  the_vmci->a[n] = mci;
  the_vmci->size++;
}

void free_vmci(Vmci* the_vmci){
    if(the_vmci == NULL) return;
  for(long i=0; i< the_vmci->size; i++){
    free(the_vmci->a[i]);
  }
  free(the_vmci->a);
  free(the_vmci);
}

// *********************************************
// *********************************************

double agmr(Accession* gtset1, Accession* gtset2, double* hgmr){
  char* gts1 = gtset1->genotypes->a;
  char* gts2 = gtset2->genotypes->a;
  long usable_pair_count = 0; // = agmr_denom
  long mismatches = 0; // = agmr_numerator
  long hgmr_denom = 0;
  long hgmr_numerator = 0;
  // fprintf(stderr, "strlen gts1, gts2: %ld %ld \n", strlen(gts1), strlen(gts2));
  for(long i=0; ;i++){
    char a1 = gts1[i];
    if(a1 == '\0') break; // end of 
    char a2 = gts2[i];
    if(DO_ASSERT) assert(a2 != '\0');
    if(a1 != MISSING_DATA_CHAR){
      if(a2 != MISSING_DATA_CHAR){
	usable_pair_count++;
	if(a1 != a2) mismatches++;
	if(a1 != '1' && a2 != '1'){
	  hgmr_denom++;
	  if(a1 != a2) hgmr_numerator++;
	}
      }
    }
  }
  *hgmr = (hgmr_denom > 0)? (double)hgmr_numerator/hgmr_denom : -1;
  return (usable_pair_count > 0)? (double)mismatches/(double)usable_pair_count : -1;
}


Vmci** find_matches(long n_ref_accessions, Vaccession* the_accessions, Chunk_pattern_ids* the_cpi,
		    //long min_usable_chunks,
		    double max_est_agmr)
{
  clock_t start = clock();
  clock_t fcmc_ticks = 0;
  
  long n_markers = the_accessions->a[0]->genotypes->length;
  long n_chunks = the_cpi->size;
  long chunk_size = the_cpi->chunk_size;
  
  long true_agmr_count = 0;
  double min_matching_chunk_fraction = pow(1.0 - max_est_agmr, chunk_size);
  Vmci** query_vmcis = (Vmci**)malloc(the_accessions->size * sizeof(Vmci*)); //
  for(long i = 0; i<the_accessions->size; i++){
    query_vmcis[i] = construct_vmci(4);
  }
  for(long i_query=n_ref_accessions; i_query< the_accessions->size; i_query++){
  
    Accession* q_gts = the_accessions->a[i_query];
    long q_md_chunk_count = q_gts->md_chunk_count;
    //clock_t ticks_before_fcmc = clock();
    Vlong* chunk_match_counts = find_chunk_match_counts(q_gts, the_cpi); //, n_ref_accessions);
    //fcmc_ticks += clock() - ticks_before_fcmc;
    
    for (long i_match = 0; i_match < i_query; i_match++){
      long matching_chunk_count = chunk_match_counts->a[i_match];
      long match_md_chunk_count = the_accessions->a[i_match]->md_chunk_count;
      // xxx
	
      double usable_chunk_count = (double)((n_chunks-q_md_chunk_count)*(n_chunks-match_md_chunk_count))/(double)n_chunks; // estimate
      
      if( //( usable_chunk_count >= 0*min_usable_chunks ) &&
	 (matching_chunk_count > min_matching_chunk_fraction*usable_chunk_count) ){
	double matching_chunk_fraction = (double)matching_chunk_count/usable_chunk_count; // fraction matching chunks
	double est_agmr = 1.0 - pow(matching_chunk_fraction, 1.0/chunk_size);
	double true_hgmr;
	double true_agmr = agmr(q_gts, the_accessions->a[i_match], &true_hgmr);
	true_agmr_count++;
	add_mci_to_vmci(query_vmcis[i_query],
			construct_mci(i_query, i_match, usable_chunk_count, matching_chunk_count, est_agmr, true_agmr, true_hgmr));	
	if(i_match >= n_ref_accessions) add_mci_to_vmci(query_vmcis[i_match],
			construct_mci(i_match, i_query, usable_chunk_count, matching_chunk_count, est_agmr, true_agmr, true_hgmr));
      }
    } // end loop over potential matches to query
    free_vlong(chunk_match_counts);
  } // end loop over queries.
  //clock_t find_matches_ticks = clock() - start; 
  /* fprintf(stderr, "# time in: find_chunk_match_count: %8.3lf; rest of find_matches: %8.3lf; find_matches total: %8.3lf\n", */
  /* 	  clock_ticks_to_seconds(fcmc_ticks), clock_ticks_to_seconds(find_matches_ticks - fcmc_ticks), clock_ticks_to_seconds(find_matches_ticks)); */
  return query_vmcis;
}

long print_results(Vaccession* the_accessions, Vmci** query_vmcis, FILE* ostream){
  long true_agmr_count = 0;
  for(long i_q=0; i_q<the_accessions->size; i_q++){
    Vmci* the_vmci = query_vmcis[i_q];
    for(long i_m=0; i_m < the_vmci->size; i_m++){
      Mci* the_mci = the_vmci->a[i_m];
      //  long match_idx = the_mci->match_index; //(the_mci->query_index == i_q)? the_mci->match_index : the_mci->query_index;
      Accession* q_gts = the_accessions->a[i_q];
      Accession* m_gts = the_accessions->a[the_mci->match_index];
      fprintf(ostream, "%5ld %30s %30s  %5.2f  %4ld  %7.4f  %7.4f    %7.4f  ", //  %7.4f\n",
	      i_q,  the_accessions->a[i_q]->id->a,  the_accessions->a[the_mci->match_index]->id->a,
	      the_mci->usable_chunks,  the_mci->n_matching_chunks,
	      the_mci->est_agmr,  the_mci->agmr,
	      the_mci->hgmr);
      //   fprintf(ostream, "%5ld %5ld %5ld %5ld", q_gts->missing_data_count, q_gts->md_chunk_count, m_gts->missing_data_count, m_gts->md_chunk_count );
      fprintf(ostream, "\n");
      true_agmr_count++;
    }
  }
  return true_agmr_count;
}

// given a vlong of matches to the gts so far in chunk (first j gts, say),
// take the next gt in query (i012) and add matching
// gts to get vlong of matches to first (j+1) gts.
// i.e. multiply each by 3, and add 0,1, or 2 as appropriate
// matching here <->  0 matched by 0, 2 by 2, 1 by anything (0,1,2) 
Vlong* matching_ipats(Vlong* ipats, long i012){
  if(i012 == 1){
    Vlong* new_vlong = construct_vlong(3*ipats->size);
    for(long i=0; i<ipats->size; i++){
      add_long_to_vlong(new_vlong, 3*ipats->a[i] + 0);
      add_long_to_vlong(new_vlong, 3*ipats->a[i] + 1);
      add_long_to_vlong(new_vlong, 3*ipats->a[i] + 2);
    }
    free_vlong(ipats);
    return new_vlong;
  }else{
    if(DO_ASSERT) assert(i012 == 0  || i012 == 2);
    for(long i=0; i<ipats->size; i++){
      ipats->a[i] = 3*ipats->a[i] + i012;
    }
    return ipats;
  }
}

char* ipat_to_strpat(long k, long ipat){ // note: this generates a string which looks like the base-3 representation of ipat, EXCEPT in reverse order
  // i.e. (for k=4) ipat=1 -> pat = '1000'
  char* pat = (char*)malloc((k+1)*sizeof(char));
  if(ipat >= 0){
    for(long i=0; i<k; i++){
      pat[i] = 48 + ipat % 3;
      ipat /= 3;
    }
  }else{
    for(long i=0; i<k; i++){
      pat[i] = 'X';
    }
  }
  pat[k] = '\0';
  // printf("ipat %ld   strpat: %s \n", ipat, pat);
  return pat;
}

long strpat_to_ipat(long len, char* strpat){ // the inverse of ipat_to_strpat
  long pat = 0;
  long f = 1;
  for(long j=0; j < len; j++){
    char a = strpat[j] - 48;
    if((a>=0) && (a<=2)){
      pat += f*a;
      f*=3;
    }else{
      pat = -1;
      break;
    }
  }
  return pat;
}

// *****  end of function definitions  *****
