// C version of 'k-mer' search for pairs of similar genotype sets.
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>
#include "vect.h"


long do_dbl_md_chunk_counts = 0;
long min_usable_chunks = 40;
double max_est_agmr = 0.2;

//***********************************************************************************************
// **************  typedefs  ********************************************************************

typedef struct{ // genotype set
  char* id;
  long index;
  long n_markers;
  char* gtset;
  Vlong* chunk_patterns;
  long md_chunk_count; // number of chunks with missing data (in at least one gt)
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
} Chunk_pattern_ids;

typedef struct{
  long query_index;
  long match_index;
  double usable_chunks; // estimate or actual count
  long n_matching_chunks;
  double est_agmr;
  double agmr;
} Mci; // 'Mci = Matching chunk info'

typedef struct{
  long capacity;
  long size;
  Mci** a;
} Vmci; 


// *********************** function declarations ************************************************

long int_power(long base, long power);
// char* ipat_to_strpat(long len, long ipat); // unused
// long strpat_to_ipat(long len, char* strpat); // unused
double agmr(Gts* gts1, Gts* gts2);
double hi_res_time(void);


// *****  Gts  **********************************************************************************
Gts* construct_gts(long index, char* id, char* gtset);
long set_gts_chunk_patterns(Gts* the_gts, Vlong* m_indices, long n_chunks, long k);
char* print_gts(Gts* the_gts);
void free_gts(Gts* the_gts);

// *****  Vgts  *********************************************************************************
Vgts* construct_vgts(long min_size);
void add_gts_to_vgts(Vgts* the_vgts, Gts* gts);
void set_vgts_chunk_patterns(Vgts* the_vgts, Vlong* m_indices, long n_chunks, long k);
void populate_chunk_pattern_ids_from_vgts(Vgts* the_vgts, Chunk_pattern_ids* the_cpi);
void print_vgts(Vgts* the_vgts);
void check_gts_indices(Vgts* the_vgts);
void free_vgts(Vgts* the_vgts);

// *****  Mci  ********
Mci* construct_mci(long qidx, long midx, double n_usable_chunks, long n_matching_chunks,
		   double est_agmr, double agmr);
// *****  Vmci  *********************************************************************************
Vmci* construct_vmci(long init_size);
void add_mci_to_vmci(Vmci* the_vmci, Mci* the_mci);
void free_vmci(Vmci* the_vmci);

// *****  Pattern_ids; indices are patterns; elements are Vlong* of accids having that pattern.
Pattern_ids* construct_pattern_ids(long n_patterns);
void free_pattern_ids(Pattern_ids*);

// *****  Chunk_pattern_ids; indices are chunk numbers; elements are Pattern_ids*
Chunk_pattern_ids* construct_chunk_pattern_ids(long n_chunks, long chunk_size);
void print_chunk_pattern_ids(Chunk_pattern_ids* the_cpi);
void free_chunk_pattern_ids(Chunk_pattern_ids* the_cpi);

// *****  Gts and Chunk_pattern_ids  ***********
Vlong* find_chunk_match_counts(Gts* the_gts, Chunk_pattern_ids* the_cpi, long n_accessions, Vlong* accidx_doublemissingdatacounts);

long find_matches_alt(Vgts* the_vgts, Chunk_pattern_ids* the_cpi, FILE* ostream);
Vmci** find_matches(Vgts* the_vgts, Chunk_pattern_ids* the_cpi);
long print_results(Vgts* the_vgts, Vmci** query_vmcis, FILE* ostream);

// *************************  end of declarations  **********************************************


// **********************************************************************************************
// ****************************** main **********************************************************
// **********************************************************************************************

int
main(int argc, char *argv[])
{
  long n_chunks = 100000; // default number of chunks (large number -> use all markers)
  long chunk_size = 8; // default number of genotype in each chunk 
  long max_number_of_accessions = -1;
  long use_alt = 0; // if true, use find_matches_alt
  
  // ***** process command line *****
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <file> options\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  char* input_filename = NULL;
  FILE *stream = NULL;
    
  int c;
  while((c = getopt(argc, argv, "i:n:k:ma")) != -1){
    switch(c){
    case 'i':
      input_filename = optarg;
      stream = fopen(input_filename, "r");
      if(stream == NULL){
	printf("Failed to open %s for reading.\n", input_filename);
	exit(EXIT_FAILURE);
      }
      break;
    case 'n': 
      n_chunks = (long)atoi(optarg);
      if(n_chunks == 0){
	printf("option n requires a numerical argument > 0\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'k':
      chunk_size = (long)atoi(optarg);
      if(chunk_size == 0){
	printf("option k requires a numerical argument > 0\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'm': // 1 to get exact number of usable chunks
      do_dbl_md_chunk_counts = 1;
      break;
    case 'a': // 1 to use find_matches_alt
      use_alt = 1;
      break;
    case '?':
      printf("? case in command line processing switch.\n");
      if ((optopt == 'i') || (optopt == 'n') || (optopt == 'k'))
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (isprint (optopt))
	fprintf (stderr, "Unknown option `-%c'.\n", optopt);
      else
	fprintf (stderr,
		 "Unknown option character: %d\n",
		 optopt);
      exit(EXIT_FAILURE);
    default:
      printf("default case (abort)\n");
      abort ();
    } // end of switch block
  } // end of loop over c.l. arguments
  if(optind < argc){
    printf("Non-option arguments. Bye.\n");
  }

  if(input_filename == NULL){
    perror("must specify input filename: -i <filename>");
    exit(EXIT_FAILURE);
  }

  printf("# input file: %s  n_chunks: %ld  chunk size: %ld  get exact md chunk counts?: %ld \n",
	 input_filename, n_chunks, chunk_size, do_dbl_md_chunk_counts);
  // *****  done processing command line  *****

  double start;
  start = hi_res_time();
  
  // ***** *****  read in genotype matrix file  ***** *****
 
  char *line = NULL;
  size_t len = 0;
  ssize_t nread;
 
  if (stream == NULL) {
    perror("fopen");
    exit(EXIT_FAILURE);
  }

  // *****  read first line: 'MARKER followed by marker ids  *****
  nread = getline(&line, &len, stream); // read first line. should have 'MARKER' and marker ids. 
  char mrkr[64];
  sscanf(line, "%s", mrkr);
  if(strcmp(mrkr, "MARKER") !=0){
    printf("Start of first line: %s ; should be MARKER. Bye.\n", mrkr); // should start with 'MARKER'
    exit(EXIT_FAILURE);
  } 
  
  // ***** read rest of file and store genotype sets and ids in Gts objs.  *****
  int max_accession_id_length = 200;
  int init_accessions_capacity = 100; //
  Vgts* the_vgts = construct_vgts(init_accessions_capacity);
  
  long n_markers;
  long gtsets_count = 0;
  if((nread = getline(&line, &len, stream)) != -1){ // process first line with genotypes 
    char* id = (char*)malloc(100*sizeof(char));
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

  while ((nread = getline(&line, &len, stream)) != -1) {
    char* id = (char*)malloc(max_accession_id_length*sizeof(char));
    char* gtstr = (char*)malloc(nread*sizeof(char)); 
    sscanf(line, "%s %s", id, gtstr);
    //   add_string_to_vstr(accession_ids, id);
    if(strlen(gtstr) != n_markers) exit(EXIT_FAILURE); // check number of genotypes is same.
    Gts* the_gts = construct_gts(gtsets_count, id, gtstr);
    add_gts_to_vgts(the_vgts, the_gts);
    check_gts_indices(the_vgts);
    gtsets_count++;
    if(max_number_of_accessions > 0  &&  gtsets_count >= max_number_of_accessions) break;
  }
  free(line);
  fclose(stream);
  
  printf("# done reading genotypes data. %ld acessions, %ld markers.  Time to read input: %8.2f\n", the_vgts->size, n_markers, hi_res_time() - start);

  // *****  done reading in genotype sets for all accessions  **********

  start = hi_res_time();
  Vlong* marker_indices = construct_vlong_whole_numbers(n_markers);
  shuffle_vlong(marker_indices); 

  if(n_chunks*chunk_size > marker_indices->size){
    n_chunks = marker_indices->size/chunk_size;
    printf("Reducing number of chunks to %ld\n", n_chunks);
  }
  set_vgts_chunk_patterns(the_vgts, marker_indices, n_chunks, chunk_size); 

  Chunk_pattern_ids* the_cpi = construct_chunk_pattern_ids(n_chunks, chunk_size);
  populate_chunk_pattern_ids_from_vgts(the_vgts, the_cpi);
  printf("# time to construct chunk_pattern_ids data structure: %12.6f\n", hi_res_time() - start);

  start = hi_res_time();
  long true_agmr_count;
  if(!use_alt){  
    Vmci** query_vmcis = find_matches(the_vgts, the_cpi);
    true_agmr_count = print_results(the_vgts, query_vmcis, stderr);
      for(long i=0; i< the_vgts->size; i++){
    free_vmci(query_vmcis[i]);
  }
      free(query_vmcis);
  }else{
    true_agmr_count = find_matches_alt(the_vgts, the_cpi, stderr);
  }
  printf("time to find candidate matches and %ld true agmrs: %9.3f\n", true_agmr_count, hi_res_time() - start);
  free_vlong(marker_indices);
  free_chunk_pattern_ids(the_cpi);
  free_vgts(the_vgts);
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
  the_gts->chunk_patterns = NULL;
  the_gts->md_chunk_count = 0;
  // long missing_count = 0;
  for(long i=0; ; i++){
    char a = the_gts->gtset[i];
    if(a == '\0'){ the_gts->n_markers = i;  break; }
    /* if((a == '0') || (a == '1') || (a == '2')){ */
    /*   // do nothing */
    /* }else{ */
    /*   missing_count++; */
    /* } */
  }
  // the_gts->missing_count = missing_count;
  // printf("id, missing count: %s  %ld \n", id, missing_count);
  return the_gts;
}

// for one accession's set of genotypes, loop over chunks and find the gt patterns. Store in the_gts->chunk_patterns
long set_gts_chunk_patterns(Gts* the_gts, Vlong* m_indices, long n_chunks, long k){
  long gts_mdchunk_count = 0;
  long n_patterns = int_power(3, k); // 3^k, the number of valid patterns, also there is a 'pattern' for missing data, making 3^k + 1 in all
  Vlong* chunk_pats = construct_vlong(n_chunks); // (Vlong*)malloc(n_chunks*sizeof(Vlong));
  for(long i_chunk=0; i_chunk < n_chunks; i_chunk++){
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
      }else{ // missing data in (at least) one of the chunks
	i_pat = n_patterns;
	gts_mdchunk_count++;
	break;
      }
    } // end of loop over the k chars in a chunk.
    add_long_to_vlong(chunk_pats, i_pat);
  } // loop over chunks.
  the_gts->chunk_patterns = chunk_pats;
  the_gts->md_chunk_count = gts_mdchunk_count;
  return gts_mdchunk_count;
}

char* print_gts(Gts* the_gts){
  printf("Gts index: %ld  gtset: %s\n", the_gts->index, the_gts->gtset);
}

void free_gts(Gts* the_gts){
  free(the_gts->id);
  free(the_gts->gtset);
  free_vlong(the_gts->chunk_patterns);
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
    
    for(long i_chunk=0; i_chunk<the_chunk_patterns->size; i_chunk++){
      long the_pat = the_chunk_patterns->a[i_chunk];
      if(the_pat >= 0){
	Vlong* the_accidxs = the_cpi->a[i_chunk]->a[the_pat];
	if(i_gts != the_gts->index){
	  printf("In populate_chunk_pattern_ids_from_vgts. indexing problem.\n"); exit(EXIT_FAILURE);
	}
	add_long_to_vlong(the_accidxs, the_gts->index);
      }else{
	printf("negative pat: %ld \n", the_pat);
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

void print_vgts(Vgts* the_vgts){
  for(int i=0; i<the_vgts->size; i++){
    print_gts(the_vgts->a[i]);
  }
}

void check_gts_indices(Vgts* the_vgts){
  for(long i=0; i<the_vgts->size; i++){
    Gts* a_gts = the_vgts->a[i];
    if(a_gts->index != i){
      printf("in check_gts_indices. i: %ld  index: %ld\n", i, a_gts->index);
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

Pattern_ids* construct_pattern_ids(long n_patterns){ // needed size known at construct time, so one param for both cap and size
  Pattern_ids* pat_ids = (Pattern_ids*)malloc(1*sizeof(Pattern_ids));
  pat_ids->capacity = n_patterns+1; // 0..n_patterns-1 are the indices of the n_patterns (=3^k) good patterns, and index n_pattern is for the missing data case.
  pat_ids->size = n_patterns+1;
  pat_ids->a = (Vlong**)malloc((n_patterns+1)*sizeof(Vlong*));
  for(int i=0; i< pat_ids->size; i++){
    pat_ids->a[i] = construct_vlong(8); // waste of memory? set to NULL until needed?
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
    chunk_pat_ids->a[i] = construct_pattern_ids(n_patterns);
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

void print_chunk_pattern_ids(Chunk_pattern_ids* the_cpi){
  for(long i_chunk=0; i_chunk<the_cpi->size; i_chunk++){
    printf("i_chunk: %ld\n", i_chunk);
    Pattern_ids* the_pi = the_cpi->a[i_chunk];
    for(long i_pat=0; i_pat<the_pi->size; i_pat++){
      printf("  i_pat: %ld \n", i_pat);
      Vlong* the_idxs = the_pi->a[i_pat];
      if(the_idxs->size > 0){
	printf("     matches: ");
	for(long ii=0; ii<the_idxs->size; ii++){
	  printf("%ld  ", the_idxs->a[ii]);
	}
	printf("\n");
      }
    }
  }
}

// *****  Gts and Chunk_pattern_ids  ***********
 
Vlong* find_chunk_match_counts(Gts* the_gts, Chunk_pattern_ids* the_cpi, long n_accessions, Vlong* accidx_dbl_md_counts){
  long n_patterns = the_cpi->n_patterns;
  Vlong* chunk_pats = the_gts->chunk_patterns;
  Vlong* accidx_matchcounts = construct_vlong_zeroes(n_accessions);
  if(do_dbl_md_chunk_counts) accidx_dbl_md_counts = construct_vlong_zeroes(n_accessions);
  for(long i_chunk=0; i_chunk < chunk_pats->size; i_chunk++){
    long the_pat = chunk_pats->a[i_chunk];  
    Vlong* chunk_match_idxs = the_cpi->a[i_chunk]->a[the_pat]; // array of indices of the matches to this chunk & pat
    // (patterns 0..n_patterns-1 are good, n_patterns=3^chunk_size is the pattern for missing data )
    if(the_pat == n_patterns){ // missing data in this chunk 
      if(accidx_dbl_md_counts == NULL) continue; // control whether to do the double missing data chunks.
      for(long i=0; i<chunk_match_idxs->size; i++){
	long accidx = chunk_match_idxs->a[i]; // index of one of the accessions matching on this chunk
	accidx_dbl_md_counts->a[accidx]++;
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

// ***** Mci  *****

Mci* construct_mci(long qidx, long midx, double usable_chunks, long n_matching_chunks,
		   // double est_matching_chunk_fraction, double matching_chunk_fraction){
		   double est_agmr, double agmr){
  Mci* the_mci = (Mci*)malloc(1*sizeof(Mci));
  the_mci->query_index = qidx;
  the_mci->match_index = midx;
  the_mci->usable_chunks = usable_chunks;
  the_mci->n_matching_chunks = n_matching_chunks;
  the_mci->est_agmr = est_agmr;
  the_mci->agmr = agmr;
  
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

long int_power(long base, long power){ // calculate base^power using integer math.
  long result = 1;
  for(int i=0; i<power; i++){
    result *= base;
  }
  return result;
}

double agmr(Gts* gtset1, Gts* gtset2){
  char* gts1 = gtset1->gtset;
  char* gts2 = gtset2->gtset;
  long usable_pair_count = 0;
  long mismatches = 0;
  for(long i=0; ;i++){
    char a1 = gts1[i];
    if(a1 == '\0') break;
    char a2 = gts2[i];
    if(a2 == '\0') break;
    if((a1 == '0') || (a1 == '1') || (a1 == '2')){
      if((a2 == '0') || (a2 == '1') || (a2 == '2')){
	usable_pair_count++;
	if(a1 != a2) mismatches++;
      }
    }
  }
  return (usable_pair_count > 0)? (double)mismatches/(double)usable_pair_count : -1;
}

double hi_res_time(void){
  return (double)clock()/(double)CLOCKS_PER_SEC;
}

Vmci** find_matches(Vgts* the_vgts, Chunk_pattern_ids* the_cpi)
{
  double start = hi_res_time();
  long n_markers = the_vgts->a[0]->n_markers;
  long n_chunks = the_cpi->size;
  long chunk_size = the_cpi->chunk_size;
    long true_agmr_count = 0;
  double min_matching_chunk_fraction = pow(1.0 - max_est_agmr, chunk_size);
  
  Vmci** query_vmcis = (Vmci**)malloc(the_vgts->size * sizeof(Vmci*)); //
  for(long i = 0; i<the_vgts->size; i++){
    query_vmcis[i] = construct_vmci(8);
  }
  for(long i_query=0; i_query< the_vgts->size; i_query++){
    Gts* q_gts = the_vgts->a[i_query];
    //   long q_md_gt_count = q_gts->missing_count;
    long q_md_chunk_count = q_gts->md_chunk_count;

    Vlong* accidx_dblmdcounts = NULL; // construct_vlong_zeroes(the_vgts->size);
    Vlong* chunk_match_counts = find_chunk_match_counts(q_gts, the_cpi, the_vgts->size, accidx_dblmdcounts);
    for(long i_match = i_query+1; i_match<the_vgts->size; i_match++){

      Gts* match_gts = the_vgts->a[i_match];
      long matching_chunk_count = chunk_match_counts->a[i_match];
      //  long match_md_gt_count = match_gts->missing_count;
      long match_md_chunk_count = match_gts->md_chunk_count;


      
      double usable_chunk_count;
      if(accidx_dblmdcounts != NULL){
	// using correct number of usable pairs
	long dbl_md_count = accidx_dblmdcounts->a[i_match]; // correct double missing data chunk count (if done in find_...
	usable_chunk_count = (double)(n_chunks - (q_md_chunk_count + match_md_chunk_count - dbl_md_count));
      }else{ // using estimated number of usable chunks
	//	long est_md_gt_count = (q_md_gt_count >= match_md_gt_count)? q_md_gt_count : match_md_gt_count;
	usable_chunk_count = 
	  (double)((n_chunks-q_md_chunk_count)*(n_chunks-match_md_chunk_count))/(double)n_chunks;
      }

      if(usable_chunk_count >= min_usable_chunks  &&  matching_chunk_count > min_matching_chunk_fraction*usable_chunk_count){
	double matching_chunk_fraction = (double)matching_chunk_count/usable_chunk_count; // fraction matching chunks
	Gts* match_gts = the_vgts->a[i_match];
	//	printf("%lf  %ld\n", matching_chunk_fraction, chunk_size);
	double est_agmr = 1.0 - pow(matching_chunk_fraction, 1.0/chunk_size);
	double true_agmr = agmr(q_gts, match_gts);
	true_agmr_count++;
	Mci* the_mci = construct_mci(i_query, i_match,  usable_chunk_count, matching_chunk_count, est_agmr, true_agmr);
	add_mci_to_vmci(query_vmcis[i_query], the_mci);
	Mci* the_other_mci = construct_mci(i_match, i_query, usable_chunk_count, matching_chunk_count, est_agmr, true_agmr);
	add_mci_to_vmci(query_vmcis[i_match], the_other_mci);
      }
    }
    if(accidx_dblmdcounts != NULL) free_vlong(accidx_dblmdcounts);
    free_vlong(chunk_match_counts);
  }
  return query_vmcis;
}

long print_results(Vgts* the_vgts, Vmci** query_vmcis, FILE* ostream){
  long true_agmr_count = 0;
  for(long i_q=0; i_q<the_vgts->size; i_q++){
    Vmci* the_vmci = query_vmcis[i_q];
    for(long i_m=0; i_m < the_vmci->size; i_m++){
      Mci* the_mci = the_vmci->a[i_m];
      //  long match_idx = the_mci->match_index; //(the_mci->query_index == i_q)? the_mci->match_index : the_mci->query_index;
      fprintf(stderr, "%5ld %30s %30s  %5.2f  %4ld  %7.4f  %7.4f\n",
	      i_q, the_vgts->a[i_q]->id, the_vgts->a[the_mci->match_index]->id,
	      the_mci->usable_chunks, the_mci->n_matching_chunks,
	      the_mci->est_agmr, the_mci->agmr
	      );
      true_agmr_count++;
    }
  }
  return true_agmr_count;
}

long find_matches_alt(Vgts* the_vgts, Chunk_pattern_ids* the_cpi, FILE* ostream) // slower
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
    Vlong* chunk_match_counts = find_chunk_match_counts(q_gts, the_cpi, the_vgts->size, accidx_dblmdcounts);
    for(long i_match=0; i_match<the_vgts->size; i_match++){
      if(i_match == i_query) continue;
      Gts* match_gts = the_vgts->a[i_match];
      long matching_chunk_count = chunk_match_counts->a[i_match];
      //   long match_md_gt_count = match_gts->missing_count;
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

      
      if(usable_chunk_count >= min_usable_chunks  &&  matching_chunk_count > min_matching_chunk_fraction*usable_chunk_count){
	double matching_chunk_fraction = (double)matching_chunk_count/usable_chunk_count; // fraction matching chunks
	double est_agmr = 1.0 - pow(matching_chunk_fraction, 1.0/chunk_size);


	
	fprintf(stderr, "%5ld %30s %30s  %5.2f  %4ld  %7.4f  %7.4f\n",
		i_query, the_vgts->a[i_query]->id, the_vgts->a[i_match]->id,
		usable_chunk_count, matching_chunk_count,
		est_agmr, agmr(q_gts, match_gts)
		);
	true_agmr_count++;
      }
    }
    if(accidx_dblmdcounts != NULL) free_vlong(accidx_dblmdcounts);
    free_vlong(chunk_match_counts);
  }
  return true_agmr_count;
} /**/

// ***** unused functions *****

char* ipat_to_strpat(long k, long ipat){
  char* pat = (char*)malloc(4*sizeof(char));
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
  // printf("ipat %ld   strpat: %s \n", ipat, pat);
  return pat;
}

long strpat_to_ipat(long len, char* strpat){
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
