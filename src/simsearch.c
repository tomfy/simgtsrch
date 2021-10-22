// C version of 'k-mer' search for pairs of similar genotype sets.
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h> // needed for getopt
#include <assert.h>
#include "vect.h"

#define DO_ASSERT 1

long do_dbl_md_chunk_counts = 0;
long min_usable_chunks = 25;
double max_est_agmr = 0.2;
// long match_count_increments_count = 0;

long do_H = 0; 
double min_hmatch_fraction = 0.2;
Vlong** ipat_hmatches; // ipat_hmatches->a[i] is a Vlong* of patterns which match i in hgmr sense.
long min_homozyg_count = 3;


//***********************************************************************************************
// **************  typedefs  ********************************************************************

typedef struct{ // genotype set
  char* id;
  long index;
  long n_markers;
  char* gtset;
  Vlong* chunk_patterns;
  long md_chunk_count; // number of chunks with missing data (in at least one gt)

  Vlong* chunk_homozyg_counts; // array of chunks with suff. many homozygous gts.
  long h_usable_chunk_count;
  long h_usable_homozygs_count;
  // long missing_count; // number of gts with missing data.
} Gts;

typedef struct{
  long capacity; // allocated size
  long size; // number of elements
  Gts** a; // array of Gts* (genotype set)
} Vgts;

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
  Vlong** hmatches; // array of Vlong*; indices are patterns, values are Vlong*s containing patterns which count as a match for hgmr purposes.
} Chunk_pattern_ids;

typedef struct{
  long query_index;
  long match_index;
  double usable_chunks; // estimate or actual count
  long n_matching_chunks;
  double est_agmr;
  double agmr;

  long n_hmatching_chunks;
  double n_h_usable_chunks;
  double hgmr; 
} Mci; // 'Mci = Matching chunk info'

typedef struct{
  long capacity;
  long size;
  Mci** a;
} Vmci; 


// *********************** function declarations ************************************************

long int_power(long base, long power);
char* ipat_to_strpat(long len, long ipat); // unused
long strpat_to_ipat(long len, char* strpat); // unused
double agmr(Gts* gts1, Gts* gts2, double* hgmr);
double hi_res_time(void);
double clock_ticks_to_seconds(clock_t nticks);
// double agmr_quick(Gts* gtset1, Gts* gtset2, double* hgmr);

// *****  Gts  **********************************************************************************
Gts* construct_gts(long index, char* id, char* gtset);
long set_gts_chunk_patterns(Gts* the_gts, Vlong* m_indices, long n_chunks, long k);
char* print_gts(Gts* the_gts, FILE* ostream);
void free_gts(Gts* the_gts);

// *****  Vgts  *********************************************************************************
Vgts* construct_vgts(long min_size);
void add_gts_to_vgts(Vgts* the_vgts, Gts* gts);
void set_vgts_chunk_patterns(Vgts* the_vgts, Vlong* m_indices, long n_chunks, long k);
void populate_chunk_pattern_ids_from_vgts(Vgts* the_vgts, Chunk_pattern_ids* the_cpi);
void print_vgts(Vgts* the_vgts, FILE* ostream);
void check_gts_indices(Vgts* the_vgts);
void free_vgts(Vgts* the_vgts);

// *****  Mci  ********
Mci* construct_mci(long qidx, long midx, double n_usable_chunks, long n_matching_chunks,
		   double est_agmr, double agmr, double hgmr);
// *****  Vmci  *********************************************************************************
Vmci* construct_vmci(long init_size);
void add_mci_to_vmci(Vmci* the_vmci, Mci* the_mci);
void free_vmci(Vmci* the_vmci);

// *****  Pattern_ids; indices are patterns; elements are Vlong* of accids having that pattern.
Pattern_ids* construct_pattern_ids(long chunk_size);
void free_pattern_ids(Pattern_ids*);

// *****  Chunk_pattern_ids; indices are chunk numbers; elements are Pattern_ids*
Chunk_pattern_ids* construct_chunk_pattern_ids(long n_chunks, long chunk_size);
void get_all_match_counts(Chunk_pattern_ids* the_cpi, Vlong** match_counts);
void print_chunk_pattern_ids(Chunk_pattern_ids* the_cpi, FILE* ostream);
void free_chunk_pattern_ids(Chunk_pattern_ids* the_cpi);

// *****  Gts and Chunk_pattern_ids  ***********
Vlong* find_chunk_match_counts(Gts* the_gts, Chunk_pattern_ids* the_cpi, long n_accessions, Vlong* accidx_doublemissingdatacounts); //, Vlong** accidx_hmatchcounts);
Vmci** find_matches(Vgts* the_vgts, Chunk_pattern_ids* the_cpi);
Vmci** find_matches_alt(Vgts* the_vgts, Vlong** match_counts, long n_chunks, long chunk_size);
long print_results(Vgts* the_vgts, Vmci** query_vmcis, FILE* ostream);

// hgmr stuff
Vlong* find_chunk_hmatch_counts(Gts* the_gts, Chunk_pattern_ids* the_cpi, long n_accessions, Vlong** accidx_hgtmatchcounts);
long find_matches_ah(Vgts* the_vgts, Chunk_pattern_ids* the_cpi, FILE* ostream);

Vlong* get_hmatches(long chunk_size, long ipat, long min_homozyg_count);
Vlong* matching_ipats(Vlong* ipats, long i012);
Vlong** get_all_hmatches(long chunk_size, long min_homozyg_count);
long hmatches01(long chunk_size, long pat1, long pat2, long* hnumer);
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
  long use_alt = 0; // if true, use find_matches_ah
  double fraction_to_analyze = 1;
  unsigned rand_seed = (unsigned)time(0);

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
  char* output_filename = NULL;
  FILE* out_stream = stdout;
    
  int c;
  while((c = getopt(argc, argv, "i:o:p:n:k:e:s:amHh:z:")) != -1){
    // i: input file name (required).
    // o: output file name. Default: output goes to stdout.
    // p: keep each accession with probability p. (For test purposes) Default = 1
    // n: number of chunks to use. Default: use each marker ~once.
    // k: chunk size (number of markers per chunk). Default: 8
    // e: max estimated agmr. Default: 0.2 (Calculate agmr only if quick est. is < this value.)
    // s: random number seed. Default: get seed from clock.
    // a: switch to invoke alternative (slightly slower) method.
    // m: switch to calculate and use exact number of usable chunks for each accession pair (slow).
    // H: switch to do quick approx. hgmr calculation (doesn't work well).
    // h: min hmatch fraction (ignored unless -H)
    // z: min homozyg count (ignored unless -H)
    switch(c){
    case 'i':
      input_filename = optarg;
      in_stream = fopen(input_filename, "r");
      if(in_stream == NULL){
	fprintf(stderr, "Failed to open %s for reading.\n", input_filename);
	exit(EXIT_FAILURE);
      }
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
    case 'p': // keep each accession with probability p (for testing with random smaller data set)
      fprintf(stderr, "%s\n", optarg);
      fraction_to_analyze = (double)atof(optarg);
      if(fraction_to_analyze <= 0  || fraction_to_analyze > 1.0){
	fprintf(stderr, "# fraction_to_analyze specified as %6.3f. Out of valid range; exiting.\n", fraction_to_analyze);
	exit(EXIT_FAILURE);
      }
      break;
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
    case 'a': // 1 to use find_matches_ah (to 'quickly' find low hgmr pairs
      use_alt = 1;
      break;
    case 'm': // 1 to get exact number of usable chunks
      do_dbl_md_chunk_counts = 1;
      break;
    case 'H': // 1 to use find_matches_ah (to 'quickly' find low hgmr pairs
      do_H = 1;
      break;
  
    case 'h': 
      min_hmatch_fraction = (double)atof(optarg);
      if(min_hmatch_fraction < 0  || min_hmatch_fraction >= 1.0){
	fprintf(stderr, "option h requires a numerical argument 0<=x<1 \n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'z': 
      min_homozyg_count = (long)atoi(optarg);
      if(min_homozyg_count <= 0){
	fprintf(stderr, "option z requires an integer argument > 0\n");
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
  }

  if(input_filename == NULL){
    perror("must specify input filename: -i <filename>");
    exit(EXIT_FAILURE);
  }
  srand(rand_seed);

  fprintf(rparam_stream, "# input file: %s  output to: %s  fraction of accessions to analyze: %5.3lf  get exact md chunk counts: %ld\n",
	  input_filename, (output_filename == NULL)? "stdout" : output_filename, fraction_to_analyze, do_dbl_md_chunk_counts);	  
 
 
  // *****  done processing command line  *****

  double start;
  start = hi_res_time();
  
  // ***** *****  read in genotype matrix file  ***** *****
 
  char *line = NULL;
  size_t len = 0;
  ssize_t nread;
 
  if (in_stream == NULL) {
    perror("fopen");
    exit(EXIT_FAILURE);
  }

  // *****  read first line: 'MARKER followed by marker ids  *****
  nread = getline(&line, &len, in_stream); // read first line. should have 'MARKER' and marker ids. 
  char mrkr[64];
  sscanf(line, "%s", mrkr);
  if(strcmp(mrkr, "MARKER") !=0){
    fprintf(stderr, "Start of first line: %s ; should be MARKER. Bye.\n", mrkr); // should start with 'MARKER'
    exit(EXIT_FAILURE);
  } 
  
  // ***** read rest of file and store genotype sets and ids in Gts objs.  *****
  int max_accession_id_length = 200;
  int init_accessions_capacity = 100; //
  Vgts* the_vgts = construct_vgts(init_accessions_capacity);
  
  long n_markers;
  long gtsets_count = 0;
  if((nread = getline(&line, &len, in_stream)) != -1){ // process first line with genotypes 
    char* id = (char*)malloc(max_accession_id_length*sizeof(char));
    char* gtstr = (char*)malloc(nread*sizeof(char)); 
    sscanf(line, "%s %s", id, gtstr);
    //    add_string_to_vstr(accession_ids, id);
    Gts* the_gts = construct_gts(gtsets_count, id, gtstr);
    n_markers = strlen(gtstr);
    add_gts_to_vgts(the_vgts, the_gts);
    gtsets_count++;
  }else{
    exit(EXIT_FAILURE);
  }

  while ((nread = getline(&line, &len, in_stream)) != -1) { // read the rest of the lines with genotypes
    if(fraction_to_analyze >= 1  || (double)rand()/RAND_MAX < fraction_to_analyze){
      char* id = (char*)malloc(max_accession_id_length*sizeof(char));
      char* gtstr = (char*)malloc(nread*sizeof(char)); 
      sscanf(line, "%s %s", id, gtstr);
      //   add_string_to_vstr(accession_ids, id);
      if(strlen(gtstr) != n_markers) exit(EXIT_FAILURE); // check number of genotypes is same.
      Gts* the_gts = construct_gts(gtsets_count, id, gtstr);
      add_gts_to_vgts(the_vgts, the_gts);
      if(DO_ASSERT) check_gts_indices(the_vgts);
      gtsets_count++;
      if(max_number_of_accessions > 0  &&  gtsets_count >= max_number_of_accessions) break;
    }
  }
  free(line);
  fclose(in_stream);
  
  fprintf(stderr, "# done reading genotypes data. %ld accessions, %ld markers.  Time to read input: %8.2f\n", the_vgts->size, n_markers, hi_res_time() - start);

  // *****  done reading in genotype sets for all accessions  **********

  start = hi_res_time();
  Vlong* marker_indices = construct_vlong_whole_numbers(n_markers);
  shuffle_vlong(marker_indices); 
  if(n_chunks*chunk_size > marker_indices->size){
    n_chunks = marker_indices->size/chunk_size;
    fprintf(stderr, "# Number of chunks will be reduced to %ld\n", n_chunks);
  }
  set_vgts_chunk_patterns(the_vgts, marker_indices, n_chunks, chunk_size);
   fprintf(rparam_stream, "# n_chunks: %ld  chunk_size: %ld  max_est_agmr: %5.3lf rng seed: %u\n", n_chunks, chunk_size, max_est_agmr, rand_seed);
	  if(do_H == 1) fprintf(rparam_stream, "# min_homozyg_count: %ld  min_hmatch_fractions: %5.3lf \n", min_homozyg_count, min_hmatch_fraction);
  fclose(rparam_stream);
  fprintf(stderr, "%s", rparam_buf);
  fprintf(out_stream, "%s", rparam_buf);
  
  //  fprintf(stderr, "after set_vgts_...\n");
  Chunk_pattern_ids* the_cpi = construct_chunk_pattern_ids(n_chunks, chunk_size);
  //  fprintf(stderr, "after construct_chunk_pattern_ids\n");
  populate_chunk_pattern_ids_from_vgts(the_vgts, the_cpi);
  fprintf(stderr, "# time to construct chunk_pattern_ids data structure: %12.6f\n", hi_res_time() - start);
  
  //Vlong** the_match_count = get_all_match_counts(the_vgts->size, the_cpi);

  start = hi_res_time();
  long true_agmr_count;
  fprintf(stderr, "# use_alt: %ld \n", use_alt);
  if(!do_H){
    Vmci** query_vmcis;
    if(!use_alt){
      query_vmcis = find_matches(the_vgts, the_cpi);
      true_agmr_count = print_results(the_vgts, query_vmcis, out_stream);
    }else{
      Vlong** match_counts = (Vlong**)malloc(the_vgts->size*sizeof(Vlong*));
      // match_counts[idx_a]->a[idx_b-idx_a]  will be the number of matching chunks between idx_a and idx_b
      for(long i_acc = 0; i_acc < the_vgts->size; i_acc++){
	match_counts[i_acc] = construct_vlong_zeroes(the_vgts->size - i_acc);
      }
    
      get_all_match_counts(the_cpi, match_counts);
      //    fprintf(stderr, "done with get_all_match_counts.\n");
      query_vmcis = find_matches_alt(the_vgts, match_counts, n_chunks, chunk_size);
      true_agmr_count = print_results(the_vgts, query_vmcis, out_stream);
    }
    for(long i=0; i< the_vgts->size; i++){
      free_vmci(query_vmcis[i]);
    }
    free(query_vmcis);    
  }else{
    ipat_hmatches = get_all_hmatches(chunk_size, min_homozyg_count);
    fprintf(stderr, "# time to get_all_hmatches %9.3lf\n", hi_res_time() - start);
    true_agmr_count = find_matches_ah(the_vgts, the_cpi, out_stream);
    fprintf(stderr, "time for find_matches_ah: %lf \n", hi_res_time() - start);
  }
  //  fprintf(stderr, "match_count_increments_count: %ld \n", match_count_increments_count);
  fprintf(stderr, "# time to find candidate matches and %ld true agmrs: %9.3f\n", true_agmr_count, hi_res_time() - start);
  fclose(out_stream);
  free_vlong(marker_indices);
  free_chunk_pattern_ids(the_cpi);
  free_vgts(the_vgts);
  fprintf(stderr, "# total simsearch run time: %9.3f\n", hi_res_time() - start0);
  exit(EXIT_SUCCESS);
}
  
// **********************  end of main  *********************************************************


// *******************  function definitions  ***************************************************
// *****  Gts  ****************************************************

Gts* construct_gts(long index, char* id, char* gtset){
  Gts* the_gts = (Gts*)malloc(1*sizeof(Gts));
  the_gts->index = index;
  the_gts->id = id;
  the_gts->gtset = gtset;
  the_gts->n_markers = strlen(gtset);
  the_gts->chunk_patterns = NULL;
  the_gts->md_chunk_count = 0;
  // long missing_count = 0;
  for(long i=0; ; i++){
    char a = the_gts->gtset[i];
    if(a == '\0'){ the_gts->n_markers = i;  break; }
  }
  // the_gts->missing_count = missing_count;
  // printf("id, missing count: %s  %ld \n", id, missing_count);
  return the_gts;
}

// for one accession's set of genotypes, loop over chunks and find the gt patterns. Store in the_gts->chunk_patterns
long set_gts_chunk_patterns(Gts* the_gts, Vlong* m_indices, long n_chunks, long k){
  long gts_mdchunk_count = 0;
  long n_patterns = int_power(3, k); // 3^k, the number of valid patterns, also there is a 'pattern' for missing data, making 3^k + 1 in all
  the_gts->chunk_homozyg_counts = construct_vlong_zeroes(n_chunks);
  Vlong* chunk_pats = construct_vlong(n_chunks); // (Vlong*)malloc(n_chunks*sizeof(Vlong));

  long n_h_usable_chunks = 0;
  long n_usable_homozygs = 0; // the number of homozygous gts in usable chunks.
  for(long i_chunk=0; i_chunk < n_chunks; i_chunk++){
    long n_homozygs = 0;
    long i_chunkstart = k*i_chunk;
    long i_pat = 0;
    long f = 1;

    // loop over characters in the chunk and construct a corresponding long index, in range [0..3^k] (3^k is the index for a chunk with any missing data)
    for(long j=0; j < k; j++){ 
      long idx = m_indices->a[i_chunkstart + j]; // 
      char a = the_gts->gtset[idx];
      long l = (long)a - 48;
      if((l>=0) && (l<=2)){ // this char is ok (0, 1, or 2, not missing data)
	i_pat += f*l;
	f *=3;
	if(l == 0  ||  l == 2){
	  n_homozygs++;
	}
      }else{ // missing data in (at least) one of the chunks
	i_pat = n_patterns;
	gts_mdchunk_count++;
	n_homozygs = -1;
	break;
      }
    } // end of loop over the k chars in a chunk.
    add_long_to_vlong(chunk_pats, i_pat);
    //   fprintf(stderr, "acc idx: %ld   i_chunk: %ld ipat: %ld  n_homozygs: %ld  n_usable_chunks: %ld \n", the_gts->index, i_chunk, i_pat, n_homozygs, n_h_usable_chunks);
    if(n_homozygs >= min_homozyg_count  &&  i_pat < n_patterns) {
      n_h_usable_chunks++; // exclude missing data chunks!
      n_usable_homozygs += n_homozygs;
    }
    the_gts->chunk_homozyg_counts->a[i_chunk] = n_homozygs;
  } // loop over chunks.
  the_gts->chunk_patterns = chunk_pats;
  /* for(long jj=0; jj< the_gts->chunk_patterns->size; jj++){ */
  /*   fprintf(stderr, "jj:  %ld    %ld  %ld \n", jj, chunk_pats->a[jj], the_gts->chunk_patterns->a[jj]); */
  /* } */
  the_gts->h_usable_chunk_count = n_h_usable_chunks;
  the_gts->h_usable_homozygs_count = n_usable_homozygs;
  the_gts->md_chunk_count = gts_mdchunk_count;
  //  fprintf(stderr, "returning from set_gts_chunk_patterns.\n");
  return gts_mdchunk_count;
}

char* print_gts(Gts* the_gts, FILE* ostream){
  fprintf(ostream, "Gts index: %ld  gtset: %s\n", the_gts->index, the_gts->gtset);
}

void free_gts(Gts* the_gts){
  free(the_gts->id);
  free(the_gts->gtset);
  free_vlong(the_gts->chunk_patterns);
  free_vlong(the_gts->chunk_homozyg_counts);
  free(the_gts);
}

// *****  Vgts  ***********************************************************

Vgts* construct_vgts(long init_cap){
  Vgts* the_vgts = (Vgts*)malloc(1*sizeof(Vgts));
  the_vgts->capacity = init_cap;
  the_vgts->size = 0;
  the_vgts->a = (Gts**)malloc(init_cap*sizeof(Gts*));
  //  printf("returning from construct_vgts\n");
  return the_vgts;
}

void add_gts_to_vgts(Vgts* the_vgts, Gts* gts){
  long cap = the_vgts->capacity;
  long n = the_vgts->size;
  if(n == cap){ // if necessary, resize w realloc
    cap *= 2;
    the_vgts->a = (Gts**)realloc(the_vgts->a, cap*sizeof(Gts*));
    the_vgts->capacity = cap;
  }
  the_vgts->a[n] = gts;
  the_vgts->size++;
}

// loop over the gts's , setting the chunk_patterns of each
void set_vgts_chunk_patterns(Vgts* the_vgts, Vlong* m_indices, long n_chunks, long k){
  long total_mdchunk_count = 0;
  for(long i=0; i < the_vgts->size; i++){
    long mdchcount = set_gts_chunk_patterns(the_vgts->a[i], m_indices, n_chunks, k);
    //  fprintf(stderr, "i: %ld  n_h_usable_chunks: %ld \n", i, the_vgts->a[i]->h_usable_chunk_count);
    //  getchar();
    total_mdchunk_count += mdchcount;
  }
}

void populate_chunk_pattern_ids_from_vgts(Vgts* the_vgts, Chunk_pattern_ids* the_cpi){
  long n_patterns = the_cpi->n_patterns;
  for(long i_gts=0; i_gts<the_vgts->size; i_gts++){
    Gts* the_gts = the_vgts->a[i_gts];
    Vlong* the_chunk_patterns = the_gts->chunk_patterns; // the gt patterns (longs) occurring in each chunk of this gts 
    long mdcount = 0;
    for(long i=0; i<the_chunk_patterns->size; i++){
      if(the_chunk_patterns->a[i] == n_patterns){ mdcount++; }
    }
    //  fprintf(stderr, "i_gts: %ld\n", i_gts);
    
    for(long i_chunk=0; i_chunk<the_chunk_patterns->size; i_chunk++){
      long the_pat = the_chunk_patterns->a[i_chunk];
      //  fprintf(stderr, "i_chunk: %ld  the_pat: %ld \n", i_chunk, the_pat);
      if(the_pat >= 0){
	Vlong* the_accidxs = the_cpi->a[i_chunk]->a[the_pat];
	if(i_gts != the_gts->index){
	  fprintf(stderr, "In populate_chunk_pattern_ids_from_vgts. indexing problem.\n"); exit(EXIT_FAILURE);
	}
	add_long_to_vlong(the_accidxs, the_gts->index);
      }else{
	fprintf(stderr, "negative pat: %ld \n", the_pat);
	exit(EXIT_FAILURE);
      }
    }
  }

  long total_mdchunk_count = 0;
  for(long i=0; i<the_cpi->size; i++){
    long chunk_md_count = the_cpi->a[i]->a[n_patterns]->size;
    total_mdchunk_count += chunk_md_count;
  } 
}

void print_vgts(Vgts* the_vgts, FILE* ostream){
  for(int i=0; i<the_vgts->size; i++){
    print_gts(the_vgts->a[i], ostream);
  }
}

void check_gts_indices(Vgts* the_vgts){
  for(long i=0; i<the_vgts->size; i++){
    Gts* a_gts = the_vgts->a[i];
    if(a_gts->index != i){
      fprintf(stderr, "In check_gts_indices. i: %ld  index: %ld\n", i, a_gts->index);
      exit(EXIT_FAILURE);
    }
  }
}

void free_vgts(Vgts* the_vgts){
  for(long i=0; i< the_vgts->size; i++){
    free_gts(the_vgts->a[i]);
  }
  free(the_vgts->a);
  free(the_vgts);
}

// *****  Pattern_ids; indices are patterns; elements are Vlong* of accids having that pattern.


Pattern_ids* construct_pattern_ids(long chunk_size){ // needed size is known at construct time, so one param for both cap and size
  long n_patterns = int_power(3, chunk_size);
  Pattern_ids* pat_ids = (Pattern_ids*)malloc(1*sizeof(Pattern_ids));
  pat_ids->capacity = n_patterns+1; // 0..n_patterns-1 are the indices of the n_patterns (=3^k) good patterns, and index n_pattern is for the missing data case.
  pat_ids->size = n_patterns+1;
  pat_ids->a = (Vlong**)malloc((n_patterns+1)*sizeof(Vlong*));
  for(long ipat=0; ipat< pat_ids->size; ipat++){
    pat_ids->a[ipat] = construct_vlong(2); // waste of memory? set to NULL until needed?
    
    //   char* spat = ipat_to_strpat(chunk_size, i);
    //  pat_ids->hmatches = get_hmatches(chunk_size, ipat, 2);
    //    long ipat = strpat_to_ipat(chunk_size, spat);
    // fprintf(stderr, "i, spat, ipat: %ld %s %ld\n", i, spat, ipat);
  }
  return pat_ids;
}

void free_pattern_ids(Pattern_ids* pat_ids){
  for(long i=0; i<pat_ids->size; i++){
    free_vlong(pat_ids->a[i]);
  }
  free(pat_ids->a);
  free(pat_ids);
}

// *****  Chunk_pattern_ids; indices are chunk numbers; elements are Pattern_ids*

Chunk_pattern_ids* construct_chunk_pattern_ids(long n_chunks, long chunk_size){ // needed size known at construct time, so one param for both cap and size
  Chunk_pattern_ids* chunk_pat_ids = (Chunk_pattern_ids*)malloc(1*sizeof(Chunk_pattern_ids));
  chunk_pat_ids->capacity = n_chunks;
  chunk_pat_ids->size = n_chunks;
  chunk_pat_ids->chunk_size = chunk_size;
  long n_patterns = int_power(3, chunk_size);
  chunk_pat_ids->n_patterns = n_patterns;
  chunk_pat_ids->a = (Pattern_ids**)malloc(n_chunks*sizeof(Pattern_ids*));
  for(int i=0; i< chunk_pat_ids->size; i++){
    chunk_pat_ids->a[i] = construct_pattern_ids(chunk_size);
  }
  return chunk_pat_ids;
}

void free_chunk_pattern_ids(Chunk_pattern_ids* the_cpi){
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

// *****  Gts and Chunk_pattern_ids  ***********


Vlong* find_chunk_match_counts(Gts* the_gts, Chunk_pattern_ids* the_cpi, long n_accessions, Vlong* accidx_dbl_md_counts){ //, Vlong** accidx_hmatchcounts){
  long n_patterns = the_cpi->n_patterns;
  Vlong* chunk_pats = the_gts->chunk_patterns;
  Vlong* accidx_matchcounts = construct_vlong_zeroes(n_accessions);
 
  if(do_dbl_md_chunk_counts) accidx_dbl_md_counts = construct_vlong_zeroes(n_accessions);
  for(long i_chunk=0; i_chunk < chunk_pats->size; i_chunk++){
    long the_pat = chunk_pats->a[i_chunk];  
    Vlong* chunk_match_idxs = the_cpi->a[i_chunk]->a[the_pat]; // array of indices of the matches to this chunk & pat
    // (patterns 0..n_patterns-1 are good, n_patterns=3^chunk_size is the pattern for missing data )
    if(the_pat == n_patterns){ // missing data in this chunk 
      if(accidx_dbl_md_counts != NULL){ // control whether to do the double missing data chunks.
	for(long i=0; i<chunk_match_idxs->size; i++){
	  long accidx = chunk_match_idxs->a[i]; // index of one of the accessions matching on this chunk
	  accidx_dbl_md_counts->a[accidx]++;
	}
      }
    }else{ // the_pat = 0..n_patterns-1 (good data)
      for(long i=0; i<chunk_match_idxs->size; i++){
	long accidx = chunk_match_idxs->a[i]; // index of one of the accessions matching on this chunk
	accidx_matchcounts->a[accidx]++; 
      }
    }
  }
  return accidx_matchcounts; 
}

void get_all_match_counts(Chunk_pattern_ids* the_cpi, Vlong** match_counts){
  double start = hi_res_time();
  long n_chunks = the_cpi->size;
  long n_patterns = the_cpi->n_patterns; 
  
  for(long i_chunk = 0; i_chunk < n_chunks; i_chunk++){
    for(long i_pat = 0; i_pat < n_patterns; i_pat++){
      Vlong* chunk_match_idxs = the_cpi->a[i_chunk]->a[i_pat];
      for(long i_1 = 0; i_1 < chunk_match_idxs->size; i_1++){
	long idx_1 = chunk_match_idxs->a[i_1];
	for(long i_2 = i_1+1; i_2 < chunk_match_idxs->size; i_2++){
	  long idx_2 = chunk_match_idxs->a[i_2];
	  //  assert(idx_2 > idx_1);
	  match_counts[idx_1]->a[idx_2 - idx_1]++;
	  // match_count_increments_count++;
	}
      }
    }
  }
  fprintf(stderr, "time in get_all_match_counts: %lf \n", hi_res_time() - start);
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
  for(long i=0; i< the_vmci->size; i++){
    free(the_vmci->a[i]);
  }
  free(the_vmci->a);
  free(the_vmci);
}

// *********************************************
// *********************************************



double agmr(Gts* gtset1, Gts* gtset2, double* hgmr){
  char* gts1 = gtset1->gtset;
  char* gts2 = gtset2->gtset;
  long usable_pair_count = 0; // = agmr_denom
  long mismatches = 0; // = agmr_numerator
  long hgmr_denom = 0;
  long hgmr_numerator = 0;
 
  for(long i=0; ;i++){
    char a1 = gts1[i];
    if(a1 == '\0') break; // end of 
    char a2 = gts2[i];
    if(DO_ASSERT) assert(a2 != '\0');
    //   if((a1 == '0') || (a1 == '1') || (a1 == '2')){
    if(a1 != '3'){
      //   if((a2 == '0') || (a2 == '1') || (a2 == '2')){
      if(a2 != '3'){
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


Vmci** find_matches(Vgts* the_vgts, Chunk_pattern_ids* the_cpi)
{
  clock_t start = clock();
  clock_t fcmc_ticks = 0;
  
  long n_markers = the_vgts->a[0]->n_markers;
  long n_chunks = the_cpi->size;
  long chunk_size = the_cpi->chunk_size;
  
  long true_agmr_count = 0;
  double min_matching_chunk_fraction = pow(1.0 - max_est_agmr, chunk_size);
  Vmci** query_vmcis = (Vmci**)malloc(the_vgts->size * sizeof(Vmci*)); //
  for(long i = 0; i<the_vgts->size; i++){
    query_vmcis[i] = construct_vmci(4);
  }
  for(long i_query=0; i_query< the_vgts->size; i_query++){
  
    Gts* q_gts = the_vgts->a[i_query];
    long q_md_chunk_count = q_gts->md_chunk_count;
    Vlong* accidx_dblmdcounts = NULL; // construct_vlong_zeroes(the_vgts->size);
    clock_t ticks_before_fcmc = clock();
    Vlong* chunk_match_counts = find_chunk_match_counts(q_gts, the_cpi, the_vgts->size, accidx_dblmdcounts);
    fcmc_ticks += clock() - ticks_before_fcmc;
    
    for(long i_match = i_query+1; i_match<the_vgts->size; i_match++){
      long matching_chunk_count = chunk_match_counts->a[i_match];
      long match_md_chunk_count = the_vgts->a[i_match]->md_chunk_count;    
      double usable_chunk_count = (accidx_dblmdcounts == NULL)?
	(double)((n_chunks-q_md_chunk_count)*(n_chunks-match_md_chunk_count))/(double)n_chunks : // estimate
	(double)(n_chunks - (q_md_chunk_count + match_md_chunk_count - accidx_dblmdcounts->a[i_match])); // exact
      
      if( ( usable_chunk_count >= min_usable_chunks )
	  &&
	  (matching_chunk_count > min_matching_chunk_fraction*usable_chunk_count) ){
	double matching_chunk_fraction = (double)matching_chunk_count/usable_chunk_count; // fraction matching chunks
	double est_agmr = 1.0 - pow(matching_chunk_fraction, 1.0/chunk_size);
	double true_hgmr;
	double true_agmr = agmr(q_gts, the_vgts->a[i_match], &true_hgmr);
	true_agmr_count++;
	add_mci_to_vmci(query_vmcis[i_query],
			construct_mci(i_query, i_match, usable_chunk_count, matching_chunk_count, est_agmr, true_agmr, true_hgmr));	
	add_mci_to_vmci(query_vmcis[i_match],
			construct_mci(i_match, i_query, usable_chunk_count, matching_chunk_count, est_agmr, true_agmr, true_hgmr));
      }
    } // end loop over potential matches to query
    if(accidx_dblmdcounts != NULL) free_vlong(accidx_dblmdcounts);
    free_vlong(chunk_match_counts);
  } // end loop over queries.
  clock_t find_matches_ticks = clock() - start; 
  fprintf(stderr, "# time in: find_chunk_match_count: %8.3lf; rest of find_matches: %8.3lf; find_matches total: %8.3lf\n",
	  clock_ticks_to_seconds(fcmc_ticks), clock_ticks_to_seconds(find_matches_ticks - fcmc_ticks), clock_ticks_to_seconds(find_matches_ticks));
  return query_vmcis;
}

long print_results(Vgts* the_vgts, Vmci** query_vmcis, FILE* ostream){
  long true_agmr_count = 0;
  for(long i_q=0; i_q<the_vgts->size; i_q++){
    Vmci* the_vmci = query_vmcis[i_q];
    for(long i_m=0; i_m < the_vmci->size; i_m++){
      Mci* the_mci = the_vmci->a[i_m];
      //  long match_idx = the_mci->match_index; //(the_mci->query_index == i_q)? the_mci->match_index : the_mci->query_index;
      fprintf(ostream, "%5ld %30s %30s  %5.2f  %4ld  %7.4f  %7.4f    %7.4lf %3ld  %7.4f\n", //  %7.4f\n",
	      i_q,  the_vgts->a[i_q]->id,  the_vgts->a[the_mci->match_index]->id,
	      the_mci->usable_chunks,  the_mci->n_matching_chunks,
	      the_mci->est_agmr,  the_mci->agmr,
	      the_mci->n_h_usable_chunks,  the_mci->n_hmatching_chunks,
	      the_mci->hgmr
	      );
      true_agmr_count++;
    }
  }
  return true_agmr_count;
}

Vmci** find_matches_alt(Vgts* the_vgts, Vlong** match_counts, long n_chunks, long chunk_size){
  double start = hi_res_time();
  double min_matching_chunk_fraction = pow(1.0 - max_est_agmr, chunk_size);

  Vmci** query_vmcis = (Vmci**)malloc(the_vgts->size * sizeof(Vmci*)); //
  for(long i = 0; i<the_vgts->size; i++){
    query_vmcis[i] = construct_vmci(4);
  }
  for(long i_q = 0; i_q < the_vgts->size; i_q++){
    Vlong* mc = match_counts[i_q];
    Gts* q_gts = the_vgts->a[i_q];
    long q_md_chunk_count = q_gts->md_chunk_count;
    for(long i_2 = 0; i_2 < mc->size; i_2++){
      long i_m = i_q + i_2;
      //  fprintf(stderr, "i_q, i_2, i_m: %ld %ld %ld \n", i_q, i_2, i_m);
      long matching_chunk_count = mc->a[i_2];
      
      //   fprintf(stderr, "   %ld %ld \n", mc->size, matching_chunk_count);
      Gts* m_gts = the_vgts->a[i_m];
      long m_md_chunk_count = m_gts->md_chunk_count;
      //  fprintf(stderr, "md_chunk_counts: %ld %ld \n", q_md_chunk_count, m_md_chunk_count);
      
      double usable_chunk_count = (double)((n_chunks-q_md_chunk_count)*(n_chunks-m_md_chunk_count))/(double)n_chunks;

      if(usable_chunk_count >= min_usable_chunks  &&  matching_chunk_count > min_matching_chunk_fraction*usable_chunk_count){
	//	fprintf(stderr, "will calc true agmr\n");
	double matching_chunk_fraction = (double)matching_chunk_count/usable_chunk_count; // fraction matching chunks
	//	Gts* match_gts = the_vgts->a[i_m];
	//	printf("%lf  %ld\n", matching_chunk_fraction, chunk_size);
	double est_agmr = 1.0 - pow(matching_chunk_fraction, 1.0/chunk_size);
	double true_hgmr;
	double true_agmr = agmr(q_gts, m_gts, &true_hgmr);
	//	true_agmr_count++;
	Mci* the_mci = construct_mci(i_q, i_m,  usable_chunk_count, matching_chunk_count, est_agmr, true_agmr, true_hgmr);
	add_mci_to_vmci(query_vmcis[i_q], the_mci);
	Mci* the_other_mci = construct_mci(i_m, i_q, usable_chunk_count, matching_chunk_count, est_agmr, true_agmr, true_hgmr);
	add_mci_to_vmci(query_vmcis[i_m], the_other_mci);
      }
      //   fprintf(stderr, "%ld %ld \n\n", i_q, i_m);
    }
  }
  fprintf(stderr, "time in find_matches_alt: %lf \n", hi_res_time() - start);
  return query_vmcis;
}


long int_power(long base, long power){ // calculate base^power using integer math.
  long result = 1;
  for(int i=0; i<power; i++){
    result *= base;
  }
  return result;
}

double hi_res_time(void){
  return (double)clock()/(double)CLOCKS_PER_SEC;
}

double clock_ticks_to_seconds(clock_t nticks){
  return (double)nticks/(double)CLOCKS_PER_SEC;
}

// **************************************************************************************
// ***** functions for finding, (or attempting to find), candidates for small hgmr  *****

Vlong* find_chunk_hmatch_counts(Gts* the_gts, Chunk_pattern_ids* the_cpi, long n_accessions, Vlong** accidx_hgtmatchcounts){
  long n_patterns = the_cpi->n_patterns;
  Vlong* chunk_pats = the_gts->chunk_patterns;
  Vlong* accidx_hmatchcounts = construct_vlong_zeroes(n_accessions);
 
  for(long i_chunk=0; i_chunk < chunk_pats->size; i_chunk++){
    long the_pat = chunk_pats->a[i_chunk];
    
    /*   */
    /*   long homozyg_count = 0; */
    /*    for(long i=0; i<g_chunk_size; i++){ */
    /*   long ipat012 = the_pat % 3; */
    /*   if(ipat012 != 1) homozyg_count++; */
    /*   the_pat /= 3; */
    /* } */
    ///   fprintf(stderr, "%ld  %ld  %s   ", i_chunk, the_pat, ipat_to_strpat(g_chunk_size, the_pat));
    long homozyg_count = the_gts->chunk_homozyg_counts->a[i_chunk];
    if(homozyg_count < 0) continue;
    //    fprintf(stderr, "%ld   \n", homozyg_count);
    /* Vlong* chunk_match_idxs = the_cpi->a[i_chunk]->a[the_pat]; // array of indices of the matches to this chunk & pat */
    /* // (patterns 0..n_patterns-1 are good, n_patterns=3^chunk_size is the pattern for missing data ) */
    /* if(the_pat == n_patterns){ // missing data in this chunk  */
    /*   if(accidx_dbl_md_counts == NULL) continue; // control whether to do the double missing data chunks. */
    /*   for(long i=0; i<chunk_match_idxs->size; i++){ */
    /* 	long accidx = chunk_match_idxs->a[i]; // index of one of the accessions matching on this chunk */
    /* 	accidx_dbl_md_counts->a[accidx]++; */
    /*   } */
    /* }else{ // the_pat = 0..n_patterns-1 (good data) */
    /*   for(long i=0; i<chunk_match_idxs->size; i++){ */
    /* 	long accidx = chunk_match_idxs->a[i]; // index of one of the accessions matching on this chunk */
    /* 	accidx_matchcounts->a[accidx]++; match_count_increments_count++; */
    /*   } */
    /* } */

    // ******************** hgmr chunk-match counting *********************
    // hgmr match counts (i.e. only 02, 20 cause non-match
    //  fprintf(stderr, "the_pat: %ld \n", the_pat, ipat_to_strpat();
    Vlong* hmatchpats = ipat_hmatches[the_pat];
    if(hmatchpats != NULL){
      for(long j=0; j<hmatchpats->size; j++){ // loop over pats that match the_pat (in hmgr sense).
	//	fprintf(stderr, "j: %ld  \n ", j);
	long jpat = hmatchpats->a[j];
	//	fprintf(stderr, "jpat: %ld \n", jpat);
	Vlong* h_chunk_match_idxs = the_cpi->a[i_chunk]->a[jpat]; // indices of accessions with pattern jpat for chunk i_chunk.
	//	fprintf(stderr, "XXX %ld %ld  %ld \n", i_chunk, jpat, h_chunk_match_idxs->size);
	for(long i=0; i<h_chunk_match_idxs->size; i++){
	  long accidx = h_chunk_match_idxs->a[i];
	  accidx_hmatchcounts->a[accidx]++;
	  (*accidx_hgtmatchcounts)->a[accidx] += homozyg_count;
	}
      }
    }
    // ******************** end of hgmr chunk-match counting ********************
  }
  return accidx_hmatchcounts; 
}


long find_matches_ah(Vgts* the_vgts, Chunk_pattern_ids* the_cpi, FILE* ostream) // slower
{
  long n_markers = the_vgts->a[0]->n_markers;
  long n_chunks = the_cpi->size;
  long chunk_size = the_cpi->chunk_size;
  long true_agmr_count = 0;
  double min_matching_chunk_fraction = pow(1.0 - max_est_agmr, chunk_size);
  
  for(long i_query=0; i_query< the_vgts->size; i_query++){
    Gts* q_gts = the_vgts->a[i_query];
    //   long q_md_gt_count = q_gts->missing_count;
    long q_md_chunk_count = q_gts->md_chunk_count;

    Vlong* accidx_dblmdcounts = NULL; // construct_vlong_zeroes(the_vgts->size);
      
    Vlong* chunk_match_counts = find_chunk_match_counts(q_gts, the_cpi, the_vgts->size, accidx_dblmdcounts); //, &chunk_hmatch_counts);
    Vlong* chunk_hgtmatch_counts = construct_vlong_zeroes(the_vgts->size);
    Vlong* chunk_hmatch_counts = find_chunk_hmatch_counts(q_gts, the_cpi, the_vgts->size, &chunk_hgtmatch_counts);
    double q_h_usable_chunk_count = q_gts->h_usable_chunk_count;
    long q_h_usable_homozygs_count = q_gts->h_usable_homozygs_count;

    for(long i_match=0; i_match<the_vgts->size; i_match++){
      if(i_match == i_query) continue;
      Gts* match_gts = the_vgts->a[i_match];
      long matching_chunk_count = chunk_match_counts->a[i_match];
      long match_md_chunk_count = match_gts->md_chunk_count;

      double usable_chunk_count;
      if(do_dbl_md_chunk_counts){
	// using correct number of usable pairs
	long dbl_md_count = accidx_dblmdcounts->a[i_match]; // correct double missing data chunk count (if done in find_...
	usable_chunk_count = (double)(n_chunks - (q_md_chunk_count + match_md_chunk_count - dbl_md_count));
      }else{ // using estimated number of usable chunks
	//	long est_md_gt_count = (q_md_gt_count >= match_md_gt_count)? q_md_gt_count : match_md_gt_count;
	usable_chunk_count = 
	  (double)((n_chunks-q_md_chunk_count)*(n_chunks-match_md_chunk_count))/(double)n_chunks;
      }

      double m_h_usable_chunk_count = match_gts->h_usable_chunk_count;
      long h_chunk_match_count = chunk_hmatch_counts->a[i_match];
      
      double h_usable_chunk_count = sqrt((double)q_h_usable_chunk_count * (double)m_h_usable_chunk_count); // using geom. avg to make it symmetric
      double est_h_usable_homozygs =  (double)q_h_usable_homozygs_count*( (double)(n_chunks - match_md_chunk_count)/(double)n_chunks );
      
      
      if(usable_chunk_count >= min_usable_chunks  &&  matching_chunk_count > min_matching_chunk_fraction*usable_chunk_count){
	double matching_chunk_fraction = (double)matching_chunk_count/usable_chunk_count; // fraction matching chunks
	double est_agmr = 1.0 - pow(matching_chunk_fraction, 1.0/chunk_size);
	double hgmr;
	double true_agmr = agmr(q_gts, match_gts, &hgmr);
	//	double m_nomd_fraction = (double)(n_chunks - match_md_chunk_count)/(double)n_chunks;
	if(true_agmr < 0.25  ||  hgmr < 0.03){
	  fprintf(ostream, "%5ld %30s %30s    %5.2f  %3ld  %6.4f  %6.4f   %6.2f %6.2f %6.2f %ld %ld %6.4f  %ld %6.4f\n",
		  i_query, the_vgts->a[i_query]->id, the_vgts->a[i_match]->id, // cols 1-3
		  usable_chunk_count, matching_chunk_count, est_agmr, true_agmr, // 4-7
		  h_usable_chunk_count, q_h_usable_chunk_count, m_h_usable_chunk_count, h_chunk_match_count, chunk_hgtmatch_counts->a[i_match], hgmr, // 8-13
		  q_h_usable_homozygs_count, est_h_usable_homozygs // 14-15
		  );
	}
	true_agmr_count++;
      }
    }
    if(accidx_dblmdcounts != NULL) free_vlong(accidx_dblmdcounts);
    free_vlong(chunk_match_counts);
  }
  return true_agmr_count;
} /**/

Vlong** get_all_hmatches(long chunk_size, long min_homozyg_count){
  long n_patterns = int_power(3, chunk_size);
  Vlong** ipat_hmatches = (Vlong**)malloc((n_patterns+1)*sizeof(Vlong*));
  for(long i=0; i<n_patterns; i++){
    //   fprintf(stderr, "Ipat: %ld\n", i);
    
    ipat_hmatches[i] = get_hmatches(chunk_size, i, min_homozyg_count);
    //  if(ipat_hmatches[i] != NULL)
    /* fprintf(stderr, "ipat: %ld  n hmatches: %ld \n", */
    /* 	    i, (ipat_hmatches[i] == NULL)? 0 : ipat_hmatches[i]->size); */
  }
  ipat_hmatches[n_patterns] = NULL;
  return ipat_hmatches;
}

//   for patterns with n homozygous gts >= min_homozyg_count
//   get 'matching' other patterns (i.e. homozyg gts in ipat_in must be matched 
Vlong* get_hmatches(long chunk_size, long ipat_in, long min_homozyg_count){ // get the patterns matching ipat_in
  long n_patterns = int_power(3, chunk_size);
  Vlong* hmatches = NULL; // this will hold the matching patterns, if ipat_in has sufficient homozygs
  Vlong* h0matches = NULL;
  Vlong* h1matches = NULL;
  if(ipat_in < n_patterns){ // skip if == n_patterns (chunks with missing data)
    long* pat_digits = (long*)malloc(chunk_size*sizeof(long));
    long homozyg_count = chunk_size;
    long ipat = ipat_in;
    for(long i=0; i<chunk_size; i++){
      long ipat012 = ipat % 3;
      if(ipat012 == 1) homozyg_count--;
      pat_digits[chunk_size-1-i] = ipat012;
      ipat /= 3;
    }  
    if(homozyg_count >= min_homozyg_count){
      hmatches = construct_vlong_zeroes(1);
      for(long i=0; i<chunk_size; i++){
	hmatches = matching_ipats(hmatches, pat_digits[i]);
      }
      long n_h0matches = 0;
      long n_h1matches = 0;
      //  fprintf(stderr, "pattern:  %5ld %s\n", ipat_in, ipat_to_strpat(chunk_size, ipat_in));
      /* for(long i=0; i<hmatches->size; i++){    */
      /*   fprintf(stderr, "  hmatch: %5ld %s\n", hmatches->a[i], ipat_to_strpat(chunk_size, hmatches->a[i])); */
      /* } */
      h0matches = construct_vlong(16);//
      for(long j=0; j<n_patterns; j++){
	
	long hnumer = 0;
	long hdenom = hmatches01(chunk_size, ipat_in, j, &hnumer);

	if(hdenom >= min_homozyg_count){
	  if( hnumer == 0){
	    //   if(h0matches == NULL) h0matches = construct_vlong(4);

	    /*q	fprintf(stderr, "%s\n", ipat_to_strpat(chunk_size, ipat_in));
	      fprintf(stderr, "%s  ", ipat_to_strpat(chunk_size, j));
	      fprintf(stderr, "%ld %ld \n\n", hnumer, hdenom); /**/
	    add_long_to_vlong(h0matches, j);
	  }else if( hnumer == 1){
	    if(h1matches == NULL) h1matches = construct_vlong(4);
	    add_long_to_vlong(h1matches, j);
	  }
	  
	}
      }
	   
      /*   fprintf(stderr, "ipat: %ld  n matches h0: %ld  h1: %ld \n", ipat_in,
	   (h0matches == NULL)? 0 : h0matches->size, (h1matches == NULL)? 0 : h1matches->size);
	   /**/
    }
    free(pat_digits);
  }
  return hmatches;
  // return h0matches;
}

long hmatches01(long chunk_size, long pat1, long pat2, long* hnumer){
  long hdenom = 0;
  for(long i=1; i<=chunk_size; i++){
    long i1 = pat1 % 3;
    pat1 /= 3;
    long i2 = pat2 % 3;
    pat2 /= 3;
    if(i1 == 0  || i1 == 2){
      if(i2 == 0  ||  i2 == 2){
	hdenom++;
	if(i1 != i2){
	  (*hnumer)++;
	  if(*hnumer > 0) return -1;
	}
      }
    }
    if(i - hdenom > chunk_size - min_homozyg_count){
      return -1; // too many pairs with a '1'
    }
  }
  return hdenom;
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
    assert(i012 == 0  || i012 == 2);
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

// *****  unused  *****
/* double agmr_quick(Gts* gtset1, Gts* gtset2, double* hgmr){ */
/*   // quit after enough markers to be confident whether */
/*   // larger or smaller than given thresholds */
/*   char* gts1 = gtset1->gtset; */
/*   char* gts2 = gtset2->gtset; */
/*   long usable_pair_count = 0; */
/*   long mismatches = 0; */
/*   long hgmr_denom = 0; */
/*   long hgmr_numerator = 0; */
/*   long i=0; */
/*   for(; ;i++){ */
/*     char a1 = gts1[i]; */
/*     if(a1 == '\0') break; */
/*     char a2 = gts2[i]; */
/*     if(a2 == '\0') break; */
/*     //   if((a1 == '0') || (a1 == '1') || (a1 == '2')){ */
/*     if(a1 != '3'){ */
/*       //   if((a2 == '0') || (a2 == '1') || (a2 == '2')){ */
/*       if(a2 != '3'){ */
/* 	usable_pair_count++; */
/* 	if(a1 != a2) mismatches++; */
/* 	if(a1 != '1' && a2 != '1'){ */
/* 	  hgmr_denom++; */
/* 	  if(a1 != a2) hgmr_numerator++; */
/* 	} */
/*       } */
/*     } */
/*     if(hgmr_denom > 0  && (i % 256 == 255)){ */
/*       int hgmr_done = 1; */
     
/*       /\*   double hgmr_t = (double)hgmr_numerator/(double)hgmr_denom; */
/* 	   double sigma_h = sqrt((double)(hgmr_numerator * (hgmr_denom - hgmr_numerator))*pow((double)hgmr_denom, -3.0)); */
/* 	   if( hgmr_t - 2*sigma_h > 0.02) hgmr_done = 1; /\**\/ */

/*       int agmr_done = 0; */
/*       double agmr_t = (double)mismatches/(double)usable_pair_count; */
/*       double sigma_a = sqrt((double)(mismatches * (usable_pair_count - mismatches))*pow((double)usable_pair_count, -3.0)); */
/*       if( */
/* 	 // agmr_t + 3*sigma_a < 0.15  || */
/* 	 agmr_t - 2*sigma_a > max_est_agmr) agmr_done = 1; */
/*       if(hgmr_done && agmr_done) { */
/* 	//	fprintf(stderr, "i: %ld dones: %2d %2d agmr: %7.3lf +- %lf  hgmr %7.3lf +- %lf \n", i, agmr_done, hgmr_done, agmr_t, sigma_a, hgmr_t, sigma_h);  */
/* 	break; */
/*       } */
/*     } */
/*   } */
/*   fprintf(stderr, "i: %ld \n", i); */
/*   *hgmr = (hgmr_denom > 0)? (double)hgmr_numerator/hgmr_denom : -1; */
/*   return (usable_pair_count > 0)? (double)mismatches/(double)usable_pair_count : -1; */
/* } */
// *****  end of function definitions  *****
