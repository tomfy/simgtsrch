// C version of 'k-mer' search for pairs of similar genotype sets.
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h> // needed for getopt
#include <assert.h>
// #include "vect.h"
#include "gtset.h"
#define UNKNOWN -1
#define DOSAGES 0
#define GENOTYPES 1
//#define UNKNOWN_FILE_TYPE 2
#define DO_ASSERT 1

//***********************************************************************************************
// **************  typedefs  ********************************************************************

typedef struct{ // genotype set
  char* id;
  long index;
  long n_markers;
  char* genotypes; // gtset; // should perhaps change name to avoid confusion with gtset.h
  Vlong* chunk_patterns;
  long md_chunk_count; // number of chunks with missing data (in at least one gt)
  long missing_data_count; // number of gts with missing data.
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
  double hgmr; 
} Mci; // 'Mci = Matching chunk info'

typedef struct{
  long capacity;
  long size;
  Mci** a;
} Vmci;

int  do_checks_flag = 0;


// *********************** function declarations ************************************************
long determine_file_format(char* filename);
long int_power(long base, long power);
char* ipat_to_strpat(long len, long ipat); // unused
long strpat_to_ipat(long len, char* strpat); // unused
double agmr(Gts* gts1, Gts* gts2, double* hgmr);
double hi_res_time(void);
double clock_ticks_to_seconds(clock_t nticks);

// *****  Gts  **********************************************************************************
Gts* construct_gts(char* id, char* gtset);
long set_gts_chunk_patterns(Gts* the_gts, Vlong* m_indices, long n_chunks, long k);
char* print_gts(Gts* the_gts, FILE* ostream);
void free_gts(Gts* the_gts);

// *****  Vgts  *********************************************************************************
Vgts* construct_vgts(long min_size);
Vgts* construct_vgts_from_genotypesset(GenotypesSet* gtset, long max_md_gts);
Vgts* construct_vgts_from_genotypes_file(char* filename, long chunk_size, double max_md_factor);
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
Pattern_ids* construct_pattern_ids(long n_patterns);
void free_pattern_ids(Pattern_ids*);

// *****  Chunk_pattern_ids; indices are chunk numbers; elements are Pattern_ids*
Chunk_pattern_ids* construct_chunk_pattern_ids(long n_chunks, long chunk_size);
// Vlong** get_all_match_counts(long n_accessions, Chunk_pattern_ids* the_cpi);
void print_chunk_pattern_ids(Chunk_pattern_ids* the_cpi, FILE* ostream);
void free_chunk_pattern_ids(Chunk_pattern_ids* the_cpi);

// *****  Gts and Chunk_pattern_ids  ***********
Vlong* find_chunk_match_counts(Gts* the_gts, Chunk_pattern_ids* the_cpi, long n_accessions); //, Vlong** accidx_hmatchcounts);
Vmci** find_matches(Vgts* the_vgts, Chunk_pattern_ids* the_cpi, long min_usable_chunks, double max_est_agmr);
// Vmci** find_matches_alt(Vgts* the_vgts, Vlong** match_counts, long n_chunks, long chunk_size);
long print_results(Vgts* the_vgts, Vmci** query_vmcis, FILE* ostream);

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
  double fraction_to_analyze = 1;
  unsigned rand_seed = (unsigned)time(0);
  double max_md_factor = 1.0; // multiplies n_markers/chunk_size
  long max_md_gts;
  long min_usable_chunks = 5;
  double max_est_agmr = 0.2;
  //  long genotype_file_type = GENOTYPES; // DOSAGES;
  double delta = 0.05;
 
  double max_marker_missing_data = 0.15;

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
  while((c = getopt(argc, argv, "i:o:p:n:k:e:s:x:")) != -1){
    // i: input file name (required).
    // o: output file name. Default: output goes to stdout.
    // p: keep each accession with probability p. (For test purposes) Default = 1
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
  }

  if(input_filename == NULL){
    perror("must specify input filename: -i <filename>");
    exit(EXIT_FAILURE);
  }
  srand(rand_seed);

  fprintf(rparam_stream, "# input file: %s  output to: %s  fraction of accessions to analyze: %5.3lf \n",
	  input_filename, (output_filename == NULL)? "stdout" : output_filename, fraction_to_analyze);	  
 
  // *****  done processing command line  *****


  
  double start;
  start = hi_res_time();

  long n_markers;
  Vgts* the_vgts;

  long genotype_file_type = determine_file_format(input_filename);
  if(genotype_file_type == UNKNOWN){
    fprintf(stderr, "# Input file has unknown format; exiting.\n");
    exit(EXIT_FAILURE);
  }else if(genotype_file_type == DOSAGES){
    fprintf(stderr, "# Input file type is dosages.\n");
     GenotypesSet* the_raw_genotypes_set = read_dosages_file_and_store(in_stream, delta);
    if(DBUG && do_checks_flag) check_genotypesset(the_raw_genotypes_set, max_marker_missing_data); 
    fprintf(stderr, "# Done reading in dosage data. %ld accessions and %ld markers. Time: %10.4lf sec.\n",
	    the_raw_genotypes_set->n_accessions, the_raw_genotypes_set->n_markers, hi_res_time() - start);
 
    // *****  clean genotypes set, i.e. remove markers with high missing data  ****
     GenotypesSet* the_genotypes_set = construct_cleaned_genotypesset(the_raw_genotypes_set, max_marker_missing_data);
  
    n_markers = the_genotypes_set->n_markers;
    fprintf(stderr, "# after construct_cleaned_genotypesset. n markers kept: %ld\n", n_markers);
    max_md_gts = (long)(max_md_factor*(double)n_markers/(double)chunk_size);
    the_vgts = construct_vgts_from_genotypesset(the_genotypes_set, max_md_gts);
  }else{
    fprintf(stderr, "# Input file type is genotypes.\n");
    // ***** *****  read in genotype matrix file  ***** *****
      the_vgts = construct_vgts_from_genotypes_file(input_filename, chunk_size, max_md_factor);
        n_markers = the_vgts->a[0]->n_markers;
  } 

  if(DO_ASSERT) check_gts_indices(the_vgts);
  fclose(in_stream);

  fprintf(rparam_stream, "# n_chunks: %ld  chunk_size: %ld  max_est_agmr: %5.3lf rng seed: %u\n", n_chunks, chunk_size, max_est_agmr, rand_seed);
  fclose(rparam_stream);
  fprintf(stderr, "%s", rparam_buf);
  fprintf(out_stream, "%s", rparam_buf);

  // *****  done reading and storing input  **********
  
  start = hi_res_time();
  Vlong* marker_indices = construct_vlong_whole_numbers(n_markers);
  shuffle_vlong(marker_indices);
  if(n_chunks*chunk_size > marker_indices->size){
    n_chunks = marker_indices->size/chunk_size;
    fprintf(stderr, "# Number of chunks will be reduced to %ld\n", n_chunks);
  }
  set_vgts_chunk_patterns(the_vgts, marker_indices, n_chunks, chunk_size);

  Chunk_pattern_ids* the_cpi = construct_chunk_pattern_ids(n_chunks, chunk_size);
  populate_chunk_pattern_ids_from_vgts(the_vgts, the_cpi);
  fprintf(stderr, "# time to construct chunk_pattern_ids structure: %12.6f\n", hi_res_time() - start);

  
  start = hi_res_time(); 
  Vmci** query_vmcis = find_matches(the_vgts, the_cpi, min_usable_chunks, max_est_agmr);
  long true_agmr_count = print_results(the_vgts, query_vmcis, out_stream);
  fprintf(stderr, "# time to find candidate matches and %ld true agmrs: %9.3f\n", true_agmr_count, hi_res_time() - start);
  
  // *****  clean up  *****
  for(long i=0; i< the_vgts->size; i++){
    free_vmci(query_vmcis[i]);
  }
  free(query_vmcis); 
  fclose(out_stream);
  free_vlong(marker_indices);
  free_chunk_pattern_ids(the_cpi);
  free_vgts(the_vgts);
  fprintf(stderr, "# total simsearch run time: %9.3f\n", hi_res_time() - start0);
  exit(EXIT_SUCCESS);
}

// **********************************************************************************************
// **********************  end of main  *********************************************************
// **********************************************************************************************


// *******************  function definitions  ***************************************************

Vgts* construct_vgts_from_genotypes_file(char* filename, long chunk_size,  double max_md_factor){
  char *line = NULL;
  size_t len = 0;
  ssize_t nread;

  FILE* in_stream = fopen(filename, "r");
  if (in_stream == NULL) {
    perror("fopen");
    exit(EXIT_FAILURE);
  }

  // *****  read up through first non-comment line, which should be: "MARKER" followed by marker ids  *****
  char a_string[64];
    while((nread = getline(&line, &len, in_stream)) != -1){ 
    sscanf(line, "%63s", a_string); // read first 63 non-whitespace chars (enough because just looking for "MARKER"
    if(a_string[0] == '#') { fprintf(stderr, "comment line\n"); continue; } // skip comments
    if(strcmp(a_string, "MARKER") == 0){
      break;
    }else{
      fprintf(stderr, "Start of first non-comment line: %s ; should be MARKER. Bye.\n", a_string);
      exit(EXIT_FAILURE);
    }
    }
  
  // ***** read rest of file and store genotype sets and ids in Gts objs.  *****
  int init_vgts_capacity = 1000; //
  Vgts* the_vgts = construct_vgts(init_vgts_capacity);
  Vgts* the_bad_vgts = construct_vgts(init_vgts_capacity); // for accessions with too much missing data
  
  long n_markers, max_md_gts;
  long n_genotype_lines_read = 0;
  while((nread = getline(&line, &len, in_stream)) != -1){ // process first line with genotypes 
    char* id = (char*)malloc(nread*sizeof(char));
    char* gtstr = (char*)malloc(nread*sizeof(char));
    long n_strs_read = sscanf(line, "%s %s %63s", id, gtstr, a_string);
    if(id[0] == '#') { free(id); free(gtstr); continue; }
    if(n_strs_read != 2){
      fprintf(stderr, "Read %ld strings on line; should be exactly 2.\n", n_strs_read);
      exit(EXIT_FAILURE);
    };

    n_genotype_lines_read++;
    if(n_genotype_lines_read == 1){
    n_markers = strlen(gtstr);
    max_md_gts = (long)(max_md_factor*(double)n_markers/(double)chunk_size);
    }else{
      assert(strlen(gtstr) == n_markers);
    }
    Gts* the_gts = construct_gts(id, gtstr); // uses memory allocated above for id, gtstr - don't free until freeing the_vgts, etc.
    if(the_gts->missing_data_count <= max_md_gts){
      add_gts_to_vgts(the_vgts, the_gts);
    }else{
      add_gts_to_vgts(the_bad_vgts, the_gts);
    }
  }
  if(n_genotype_lines_read == 0){
    fprintf(stderr, "Input file %s has no data.\n", filename);
    exit(EXIT_FAILURE);
  }

  free(line);
  fprintf(stderr, "# max missing data markers: %ld; %ld good accessions, %ld bad accessions \n", max_md_gts, the_vgts->size, the_bad_vgts->size);
  //  fprintf(stderr, "# done reading genotypes data. %ld accessions, %ld markers.  Time to read input: %8.2f\n", the_vgts->size, n_markers, hi_res_time() - start);
  free_vgts(the_bad_vgts); //fprintf(stderr, "bad vgts freed.\n");
  return the_vgts;
}


long determine_file_format(char* filename){
  // determine whether file is dosages, genotypes, or other
  FILE* in_stream = fopen(filename, "r");

  if(in_stream == NULL){
    fprintf(stderr, "Failed to open %s for reading.\n", filename);
    exit(EXIT_FAILURE);
  }
  char *line = NULL;
  size_t len = 0;
  ssize_t nread;
  long return_value = UNKNOWN;
  
  // read any comment lines and first non-comment line, which should start with 'MARKER'
  char string1[64];
  while((nread = getline(&line, &len, in_stream)) != -1){   
    sscanf(line, "%63s", string1);
    if(string1[0] == '#') continue; // skip comments
    if(strcmp(string1, "MARKER") == 0){
      break;
    }else{
      return return_value;
      fprintf(stderr, "Start of first line: %s ; should be MARKER. Bye.\n", string1); // should start with 'MARKER'
      exit(EXIT_FAILURE);
    }
  }
  nread = getline(&line, &len, in_stream);
  char* id = (char*)malloc(nread*sizeof(char));
  char* gtstr = (char*)malloc(nread*sizeof(char));
  char a_string[64];
  long n_strs_read = (long)sscanf(line, "%s %s %63s", id, gtstr, a_string);
  // fprintf(stderr, "n_strs_read: %ld\n", n_strs_read);
  if(n_strs_read == 2){
    return_value = GENOTYPES;
  }else if(n_strs_read == 3){ // at least 3 fields present; assume DOSAGES
    return_value = DOSAGES;
  }
  free(id); free(gtstr);
  return return_value;
}
   

// *****  Gts  ****************************************************
Gts* construct_gts(char* id, char* genotypes){ // does not alloc memory for id, genotypes and copy them, just uses already allocated memory.
  Gts* the_gts = (Gts*)malloc(1*sizeof(Gts));
  the_gts->id = id;
  the_gts->genotypes = genotypes;
  the_gts->n_markers = strlen(genotypes);
  the_gts->chunk_patterns = NULL;
  the_gts->md_chunk_count = 0;
  long md_count = 0;
  for(long i=0; ; i++){
    char a = the_gts->genotypes[i];
    if(a == '\0'){ the_gts->n_markers = i;  break; }
    if(a == '3'){
      md_count++;
    }else if(!(a == '0' || a == '1' || a == '2')){
      fprintf(stderr, "# disallowed genotype: %c found. Only 0,1,2,3 allowed.\n", a);
    }
  }
  the_gts->missing_data_count = md_count;
  // printf("id, missing count: %s  %ld \n", id, md_count);
  return the_gts;
}

// for one accession's set of genotypes, loop over chunks and find the gt patterns. Store in the_gts->chunk_patterns
long set_gts_chunk_patterns(Gts* the_gts, Vlong* m_indices, long n_chunks, long k){
  long gts_mdchunk_count = 0;
  long n_patterns = int_power(3, k); // 3^k, the number of valid patterns, also there is a 'pattern' for missing data, making 3^k + 1 in all
  //  the_gts->chunk_homozyg_counts = construct_vlong_zeroes(n_chunks);
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
      char a = the_gts->genotypes[idx];
      long l = (long)a - 48;
      if((l>=0) && (l<=2)){ // this char is ok (0, 1, or 2, not missing data)
	i_pat += f*l;
	f *=3;
	if(l == 0  ||  l == 2){
	  n_homozygs++;
	}
      }else{ // missing data in (at least) one of the markers in the chunk
	i_pat = n_patterns;
	gts_mdchunk_count++;
	n_homozygs = -1;
	break;
      }
    } // end of loop over the k chars in a chunk.
    add_long_to_vlong(chunk_pats, i_pat);
    
  } // loop over chunks.
  the_gts->chunk_patterns = chunk_pats;
  the_gts->md_chunk_count = gts_mdchunk_count;
  return gts_mdchunk_count;
}

char* print_gts(Gts* the_gts, FILE* ostream){
  // fprintf(ostream, "Gts index: %ld  gtset: %s\n", the_gts->index, the_gts->genotypes);
  // fprintf(ostream, "%s  %s\n", the_gts->id, the_gts->genotypes);
  fprintf(ostream, "XXXXXX %s %ld  %ld %ld %ld\n",
	  the_gts->id, the_gts->index, the_gts->n_markers, strlen(the_gts->genotypes), the_gts->missing_data_count);
}

void free_gts(Gts* the_gts){
    if(the_gts == NULL) return;
  free(the_gts->id);
  free(the_gts->genotypes);
  free_vlong(the_gts->chunk_patterns);
  free(the_gts);
}

// *****  Vgts  ***********************************************************

Vgts* construct_vgts(long init_cap){
  Vgts* the_vgts = (Vgts*)malloc(1*sizeof(Vgts));
  the_vgts->capacity = init_cap;
  the_vgts->size = 0;
  the_vgts->a = (Gts**)malloc(init_cap*sizeof(Gts*));
  return the_vgts;
}

Vgts* construct_vgts_from_genotypesset(GenotypesSet* gtset, long max_md_gts){
  Vgts* the_vgts = construct_vgts(gtset->n_accessions);
  long bad_accession_count = 0;
  for(long i = 0; i<gtset->accessions->size; i++){
    Accession* acc = gtset->accessions->a[i];
    Gts* a_gts = (Gts*)malloc(sizeof(Gts));
    a_gts->id = (char*)malloc((acc->id->length + 1)+sizeof(char)); // +1 is for the terminating null char
    strcpy(a_gts->id, acc->id->a);
    a_gts->index = acc->index;
    a_gts->n_markers = gtset->n_markers;
    a_gts->genotypes = (char*)malloc((acc->genotypes->length + 1)+sizeof(char)); // +1 is for the terminating null char
    strcpy(a_gts->genotypes, acc->genotypes->a);
    a_gts->missing_data_count = acc->missing_data_count;
    //  print_gts(a_gts, stderr);
    if(a_gts->missing_data_count <= max_md_gts){ 
      add_gts_to_vgts(the_vgts, a_gts);
    }else{
      bad_accession_count++;
    }
  }
  fprintf(stderr, "# good/bad accessions counts: %ld %ld\n", the_vgts->size, bad_accession_count);
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
  gts->index = n;
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
  if(the_vgts == NULL) return;
  for(long i=0; i< the_vgts->size; i++){
    free_gts(the_vgts->a[i]);
  }
  free(the_vgts->a);
  free(the_vgts);
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

// *****  Gts and Chunk_pattern_ids  ***********

Vlong* find_chunk_match_counts(Gts* the_gts, Chunk_pattern_ids* the_cpi, long n_accessions){ //, Vlong** accidx_hmatchcounts){
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

double agmr(Gts* gtset1, Gts* gtset2, double* hgmr){
  char* gts1 = gtset1->genotypes;
  char* gts2 = gtset2->genotypes;
  long usable_pair_count = 0; // = agmr_denom
  long mismatches = 0; // = agmr_numerator
  long hgmr_denom = 0;
  long hgmr_numerator = 0;
 
  for(long i=0; ;i++){
    char a1 = gts1[i];
    if(a1 == '\0') break; // end of 
    char a2 = gts2[i];
    if(DO_ASSERT) assert(a2 != '\0');
    if(a1 != '3'){
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


Vmci** find_matches(Vgts* the_vgts, Chunk_pattern_ids* the_cpi, long min_usable_chunks, double max_est_agmr)
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
    clock_t ticks_before_fcmc = clock();
    Vlong* chunk_match_counts = find_chunk_match_counts(q_gts, the_cpi, the_vgts->size);
    fcmc_ticks += clock() - ticks_before_fcmc;
    
    for (long i_match = 0; i_match < i_query; i_match++){
      long matching_chunk_count = chunk_match_counts->a[i_match];
      long match_md_chunk_count = the_vgts->a[i_match]->md_chunk_count;    
      double usable_chunk_count = (double)((n_chunks-q_md_chunk_count)*(n_chunks-match_md_chunk_count))/(double)n_chunks; // estimate
      
      if( //( usable_chunk_count >= 0*min_usable_chunks ) &&
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
      Gts* q_gts = the_vgts->a[i_q];
      Gts* m_gts = the_vgts->a[the_mci->match_index];
      fprintf(ostream, "%5ld %30s %30s  %5.2f  %4ld  %7.4f  %7.4f    %7.4f  ", //  %7.4f\n",
	      i_q,  the_vgts->a[i_q]->id,  the_vgts->a[the_mci->match_index]->id,
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

/* Vlong* find_chunk_hmatch_counts(Gts* the_gts, Chunk_pattern_ids* the_cpi, long n_accessions, Vlong** accidx_hgtmatchcounts){ */
/*   long n_patterns = the_cpi->n_patterns; */
/*   Vlong* chunk_pats = the_gts->chunk_patterns; */
/*   Vlong* accidx_hmatchcounts = construct_vlong_zeroes(n_accessions); */
 
/*   for(long i_chunk=0; i_chunk < chunk_pats->size; i_chunk++){ */
/*     long the_pat = chunk_pats->a[i_chunk]; */
    
/*     /\*   *\/ */
/*     /\*   long homozyg_count = 0; *\/ */
/*     /\*    for(long i=0; i<g_chunk_size; i++){ *\/ */
/*     /\*   long ipat012 = the_pat % 3; *\/ */
/*     /\*   if(ipat012 != 1) homozyg_count++; *\/ */
/*     /\*   the_pat /= 3; *\/ */
/*     /\* } *\/ */
/*     ///   fprintf(stderr, "%ld  %ld  %s   ", i_chunk, the_pat, ipat_to_strpat(g_chunk_size, the_pat)); */
/*     long homozyg_count = the_gts->chunk_homozyg_counts->a[i_chunk]; */
/*     if(homozyg_count < 0) continue; */
/*     //    fprintf(stderr, "%ld   \n", homozyg_count); */
/*     /\* Vlong* chunk_match_idxs = the_cpi->a[i_chunk]->a[the_pat]; // array of indices of the matches to this chunk & pat *\/ */
/*     /\* // (patterns 0..n_patterns-1 are good, n_patterns=3^chunk_size is the pattern for missing data ) *\/ */
/*     /\* if(the_pat == n_patterns){ // missing data in this chunk  *\/ */
/*     /\*   if(accidx_dbl_md_counts == NULL) continue; // control whether to do the double missing data chunks. *\/ */
/*     /\*   for(long i=0; i<chunk_match_idxs->size; i++){ *\/ */
/*     /\* 	long accidx = chunk_match_idxs->a[i]; // index of one of the accessions matching on this chunk *\/ */
/*     /\* 	accidx_dbl_md_counts->a[accidx]++; *\/ */
/*     /\*   } *\/ */
/*     /\* }else{ // the_pat = 0..n_patterns-1 (good data) *\/ */
/*     /\*   for(long i=0; i<chunk_match_idxs->size; i++){ *\/ */
/*     /\* 	long accidx = chunk_match_idxs->a[i]; // index of one of the accessions matching on this chunk *\/ */
/*     /\* 	accidx_matchcounts->a[accidx]++; match_count_increments_count++; *\/ */
/*     /\*   } *\/ */
/*     /\* } *\/ */

/*     // ******************** hgmr chunk-match counting ********************* */
/*     // hgmr match counts (i.e. only 02, 20 cause non-match */
/*     //  fprintf(stderr, "the_pat: %ld \n", the_pat, ipat_to_strpat(); */
/*     Vlong* hmatchpats = ipat_hmatches[the_pat]; */
/*     if(hmatchpats != NULL){ */
/*       for(long j=0; j<hmatchpats->size; j++){ // loop over pats that match the_pat (in hmgr sense). */
/* 	//	fprintf(stderr, "j: %ld  \n ", j); */
/* 	long jpat = hmatchpats->a[j]; */
/* 	//	fprintf(stderr, "jpat: %ld \n", jpat); */
/* 	Vlong* h_chunk_match_idxs = the_cpi->a[i_chunk]->a[jpat]; // indices of accessions with pattern jpat for chunk i_chunk. */
/* 	//	fprintf(stderr, "XXX %ld %ld  %ld \n", i_chunk, jpat, h_chunk_match_idxs->size); */
/* 	for(long i=0; i<h_chunk_match_idxs->size; i++){ */
/* 	  long accidx = h_chunk_match_idxs->a[i]; */
/* 	  accidx_hmatchcounts->a[accidx]++; */
/* 	  (*accidx_hgtmatchcounts)->a[accidx] += homozyg_count; */
/* 	} */
/*       } */
/*     } */
/*     // ******************** end of hgmr chunk-match counting ******************** */
/*   } */
/*   return accidx_hmatchcounts;  */
/* } */


/* long find_matches_ah(Vgts* the_vgts, Chunk_pattern_ids* the_cpi, FILE* ostream) // slower */
/* { */
/*   long n_markers = the_vgts->a[0]->n_markers; */
/*   long n_chunks = the_cpi->size; */
/*   long chunk_size = the_cpi->chunk_size; */
/*   long true_agmr_count = 0; */
/*   double min_matching_chunk_fraction = pow(1.0 - max_est_agmr, chunk_size); */
  
/*   for(long i_query=0; i_query< the_vgts->size; i_query++){ */
/*     Gts* q_gts = the_vgts->a[i_query]; */
/*     //   long q_md_gt_count = q_gts->missing_count; */
/*     long q_md_chunk_count = q_gts->md_chunk_count; */
      
/*     Vlong* chunk_match_counts = find_chunk_match_counts(q_gts, the_cpi, the_vgts->size); //, &chunk_hmatch_counts); */
/*     Vlong* chunk_hgtmatch_counts = construct_vlong_zeroes(the_vgts->size); */
/*     Vlong* chunk_hmatch_counts = find_chunk_hmatch_counts(q_gts, the_cpi, the_vgts->size, &chunk_hgtmatch_counts); */
/*     double q_h_usable_chunk_count = q_gts->h_usable_chunk_count; */
/*     long q_h_usable_homozygs_count = q_gts->h_usable_homozygs_count; */

/*     for(long i_match=0; i_match<the_vgts->size; i_match++){ */
/*       if(i_match == i_query) continue; */
/*       Gts* match_gts = the_vgts->a[i_match]; */
/*       long matching_chunk_count = chunk_match_counts->a[i_match]; */
/*       long match_md_chunk_count = match_gts->md_chunk_count; */

/*       double usable_chunk_count =  */
/* 	  (double)((n_chunks-q_md_chunk_count)*(n_chunks-match_md_chunk_count))/(double)n_chunks; */

/*       double m_h_usable_chunk_count = match_gts->h_usable_chunk_count; */
/*       long h_chunk_match_count = chunk_hmatch_counts->a[i_match]; */
      
/*       double h_usable_chunk_count = sqrt((double)q_h_usable_chunk_count * (double)m_h_usable_chunk_count); // using geom. avg to make it symmetric */
/*       double est_h_usable_homozygs =  (double)q_h_usable_homozygs_count*( (double)(n_chunks - match_md_chunk_count)/(double)n_chunks ); */
      
      
/*       if(usable_chunk_count >= min_usable_chunks  &&  matching_chunk_count > min_matching_chunk_fraction*usable_chunk_count){ */
/* 	double matching_chunk_fraction = (double)matching_chunk_count/usable_chunk_count; // fraction matching chunks */
/* 	double est_agmr = 1.0 - pow(matching_chunk_fraction, 1.0/chunk_size); */
/* 	double hgmr; */
/* 	double true_agmr = agmr(q_gts, match_gts, &hgmr); */
/* 	if(true_agmr < 0.25  ||  hgmr < 0.03){ */
/* 	  fprintf(ostream, "%5ld %30s %30s    %5.2f  %3ld  %6.4f  %6.4f   %6.2f %6.2f %6.2f %ld %ld %6.4f  %ld %6.4f\n", */
/* 		  i_query, the_vgts->a[i_query]->id, the_vgts->a[i_match]->id, // cols 1-3 */
/* 		  usable_chunk_count, matching_chunk_count, est_agmr, true_agmr, // 4-7 */
/* 		  h_usable_chunk_count, q_h_usable_chunk_count, m_h_usable_chunk_count, h_chunk_match_count, chunk_hgtmatch_counts->a[i_match], hgmr, // 8-13 */
/* 		  q_h_usable_homozygs_count, est_h_usable_homozygs // 14-15 */
/* 		  ); */
/* 	} */
/* 	true_agmr_count++; */
/*       } */
/*     } */
/*     free_vlong(chunk_match_counts); */
/*   } */
/*   return true_agmr_count; */
/* } /\**\/ */

/* Vlong** get_all_hmatches(long chunk_size, long min_homozyg_count){ */
/*   long n_patterns = int_power(3, chunk_size); */
/*   Vlong** ipat_hmatches = (Vlong**)malloc((n_patterns+1)*sizeof(Vlong*)); */
/*   for(long i=0; i<n_patterns; i++){ */
/*     //   fprintf(stderr, "Ipat: %ld\n", i); */
    
/*     ipat_hmatches[i] = get_hmatches(chunk_size, i, min_homozyg_count); */
/*     //  if(ipat_hmatches[i] != NULL) */
/*     /\* fprintf(stderr, "ipat: %ld  n hmatches: %ld \n", *\/ */
/*     /\* 	    i, (ipat_hmatches[i] == NULL)? 0 : ipat_hmatches[i]->size); *\/ */
/*   } */
/*   ipat_hmatches[n_patterns] = NULL; */
/*   return ipat_hmatches; */
/* } */

/* //   for patterns with n homozygous gts >= min_homozyg_count */
/* //   get 'matching' other patterns (i.e. homozyg gts in ipat_in must be matched  */
/* Vlong* get_hmatches(long chunk_size, long ipat_in, long min_homozyg_count){ // get the patterns matching ipat_in */
/*   long n_patterns = int_power(3, chunk_size); */
/*   Vlong* hmatches = NULL; // this will hold the matching patterns, if ipat_in has sufficient homozygs */
/*   Vlong* h0matches = NULL; */
/*   Vlong* h1matches = NULL; */
/*   if(ipat_in < n_patterns){ // skip if == n_patterns (chunks with missing data) */
/*     long* pat_digits = (long*)malloc(chunk_size*sizeof(long)); */
/*     long homozyg_count = chunk_size; */
/*     long ipat = ipat_in; */
/*     for(long i=0; i<chunk_size; i++){ */
/*       long ipat012 = ipat % 3; */
/*       if(ipat012 == 1) homozyg_count--; */
/*       pat_digits[chunk_size-1-i] = ipat012; */
/*       ipat /= 3; */
/*     }   */
/*     if(homozyg_count >= min_homozyg_count){ */
/*       hmatches = construct_vlong_zeroes(1); */
/*       for(long i=0; i<chunk_size; i++){ */
/* 	hmatches = matching_ipats(hmatches, pat_digits[i]); */
/*       } */
/*       long n_h0matches = 0; */
/*       long n_h1matches = 0; */
/*       //  fprintf(stderr, "pattern:  %5ld %s\n", ipat_in, ipat_to_strpat(chunk_size, ipat_in)); */
/*       /\* for(long i=0; i<hmatches->size; i++){    *\/ */
/*       /\*   fprintf(stderr, "  hmatch: %5ld %s\n", hmatches->a[i], ipat_to_strpat(chunk_size, hmatches->a[i])); *\/ */
/*       /\* } *\/ */
/*       h0matches = construct_vlong(16);// */
/*       for(long j=0; j<n_patterns; j++){ */
	
/* 	long hnumer = 0; */
/* 	long hdenom = hmatches01(chunk_size, ipat_in, j, &hnumer); */

/* 	if(hdenom >= min_homozyg_count){ */
/* 	  if( hnumer == 0){ */
/* 	    //   if(h0matches == NULL) h0matches = construct_vlong(4); */

/* 	    /\*q	fprintf(stderr, "%s\n", ipat_to_strpat(chunk_size, ipat_in)); */
/* 	      fprintf(stderr, "%s  ", ipat_to_strpat(chunk_size, j)); */
/* 	      fprintf(stderr, "%ld %ld \n\n", hnumer, hdenom); /\**\/ */
/* 	    add_long_to_vlong(h0matches, j); */
/* 	  }else if( hnumer == 1){ */
/* 	    if(h1matches == NULL) h1matches = construct_vlong(4); */
/* 	    add_long_to_vlong(h1matches, j); */
/* 	  } */
	  
/* 	} */
/*       } */
	   
/*       /\*   fprintf(stderr, "ipat: %ld  n matches h0: %ld  h1: %ld \n", ipat_in, */
/* 	   (h0matches == NULL)? 0 : h0matches->size, (h1matches == NULL)? 0 : h1matches->size); */
/* 	   /\**\/ */
/*     } */
/*     free(pat_digits); */
/*   } */
/*   return hmatches; */
/*   // return h0matches; */
/* } */

/* long hmatches01(long chunk_size, long pat1, long pat2, long* hnumer){ */
/*   long hdenom = 0; */
/*   for(long i=1; i<=chunk_size; i++){ */
/*     long i1 = pat1 % 3; */
/*     pat1 /= 3; */
/*     long i2 = pat2 % 3; */
/*     pat2 /= 3; */
/*     if(i1 == 0  || i1 == 2){ */
/*       if(i2 == 0  ||  i2 == 2){ */
/* 	hdenom++; */
/* 	if(i1 != i2){ */
/* 	  (*hnumer)++; */
/* 	  if(*hnumer > 0) return -1; */
/* 	} */
/*       } */
/*     } */
/*     if(i - hdenom > chunk_size - min_homozyg_count){ */
/*       return -1; // too many pairs with a '1' */
/*     } */
/*   } */
/*   return hdenom; */
/* } */

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

// *****  unused  *****

/* Pattern_ids* construct_pattern_ids_old(long n_patterns){ */
/*   Pattern_ids* pat_ids = (Pattern_ids*)malloc(1*sizeof(Pattern_ids)); */
/*   pat_ids->capacity = n_patterns+1; // 0..n_patterns-1 are the indices of the n_patterns (=3^k) good patterns, and index n_pattern is for the missing data case. */
/*   pat_ids->size = n_patterns+1; */
/*   pat_ids->a = (Vlong**)malloc((n_patterns+1)*sizeof(Vlong*)); */
/*   for(long ipat=0; ipat< pat_ids->size; ipat++){ */
/*     pat_ids->a[ipat] = construct_vlong(1); // waste of memory? set to NULL until needed? */
/*   } */
/*   return pat_ids; */
/* } */

/* Chunk_pattern_ids* construct_chunk_pattern_ids_old(long n_chunks, long chunk_size){ // needed size is known at construct time, so one param for both cap and size */
/*   Chunk_pattern_ids* chunk_pat_ids = (Chunk_pattern_ids*)malloc(1*sizeof(Chunk_pattern_ids)); */
/*   chunk_pat_ids->capacity = n_chunks; */
/*   chunk_pat_ids->size = n_chunks; */
/*   chunk_pat_ids->chunk_size = chunk_size; */
/*   long n_patterns = int_power(3, chunk_size); */
/*   chunk_pat_ids->n_patterns = n_patterns; */
/*   chunk_pat_ids->a = (Pattern_ids**)malloc(n_chunks*sizeof(Pattern_ids*)); */
/*   for(int i=0; i< chunk_pat_ids->size; i++){ // for each chunk */
/*     chunk_pat_ids->a[i] = construct_pattern_ids_old(n_patterns); // construct a pattern_ids */
/*   } */
/*   return chunk_pat_ids; */
/* } */

/* void populate_chunk_pattern_ids_from_vgts_old(Vgts* the_vgts, Chunk_pattern_ids* the_cpi){ */
/*   long n_patterns = the_cpi->n_patterns; */
/*   for(long i_gts=0; i_gts<the_vgts->size; i_gts++){ */
/*     Gts* the_gts = the_vgts->a[i_gts]; */
/*     Vlong* the_chunk_patterns = the_gts->chunk_patterns; // the gt patterns (longs) occurring in each chunk of this gts  */
/*     long mdcount = 0; */
/*     for(long i=0; i<the_chunk_patterns->size; i++){ */
/*       if(the_chunk_patterns->a[i] == n_patterns){ mdcount++; } */
/*     } */
/*     //  fprintf(stderr, "i_gts: %ld\n", i_gts); */
    
/*     for(long i_chunk=0; i_chunk<the_chunk_patterns->size; i_chunk++){ */
/*       long the_pat = the_chunk_patterns->a[i_chunk]; */
/*       //  fprintf(stderr, "i_chunk: %ld  the_pat: %ld \n", i_chunk, the_pat); */
/*       //  if(the_pat >= 0){ */
/*       assert(the_pat >= 0); */
/*       Vlong* the_accidxs = the_cpi->a[i_chunk]->a[the_pat]; */
/*       assert(i_gts == the_gts->index); */
/*       /\* if(i_gts != the_gts->index){ *\/ */
/*       /\*   fprintf(stderr, "In populate_chunk_pattern_ids_from_vgts. indexing problem.\n"); exit(EXIT_FAILURE); *\/ */
/*       /\* } *\/ */
/*       add_long_to_vlong(the_accidxs, the_gts->index); */
/*       /\* }else{ *\/ */
/*       /\* 	fprintf(stderr, "negative pat: %ld \n", the_pat); *\/ */
/*       /\* 	exit(EXIT_FAILURE); *\/ */
/*       /\* } *\/ */
/*     } */
/*   } */

/*   long total_mdchunk_count = 0; */
/*   for(long i=0; i<the_cpi->size; i++){ */
/*     long chunk_md_count = the_cpi->a[i]->a[n_patterns]->size; */
/*     total_mdchunk_count += chunk_md_count; */
/*   }  */
/* } */


// alt method (slower, more memory)

/* Vlong** get_all_match_counts(long n_accessions, Chunk_pattern_ids* the_cpi){ */
/*   long n_chunks = the_cpi->size; */
/*   long n_patterns = the_cpi->n_patterns; */
/*   Vlong** match_counts = (Vlong**)malloc(n_accessions*sizeof(Vlong*)); */
/*   // match_counts[idx_a]->a[idx_b-idx_a]  will be the number of matching chunks between idx_a and idx_b */
/*   for(long i_acc = 0; i_acc < n_accessions; i_acc++){ */
/*     match_counts[i_acc] = construct_vlong_zeroes(n_accessions - i_acc); */
/*   } */
  
/*   for(long i_chunk = 0; i_chunk < n_chunks; i_chunk++){ */
/*     for(long i_pat = 0; i_pat < n_patterns; i_pat++){ */
/*       Vlong* chunk_match_idxs = the_cpi->a[i_chunk]->a[i_pat]; */
/*       if(chunk_match_idxs != NULL){ */
/* 	for(long i_1 = 0; i_1 < chunk_match_idxs->size; i_1++){ */
/* 	  long idx_1 = chunk_match_idxs->a[i_1]; */
/* 	  for(long i_2 = i_1+1; i_2 < chunk_match_idxs->size; i_2++){ */
/* 	    long idx_2 = chunk_match_idxs->a[i_2]; */
/* 	    match_counts[idx_1]->a[idx_2 - idx_1]++; */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   return match_counts; */
/* } */

/* Vmci** find_matches_alt(Vgts* the_vgts, Vlong** match_counts, long n_chunks, long chunk_size){ // slower */
/*   double start = hi_res_time(); */
/*   double min_matching_chunk_fraction = pow(1.0 - max_est_agmr, chunk_size); */

/*   Vmci** query_vmcis = (Vmci**)malloc(the_vgts->size * sizeof(Vmci*)); // */
/*   for(long i = 0; i<the_vgts->size; i++){ */
/*     query_vmcis[i] = construct_vmci(4); */
/*   } */
/*   for(long i_q = 0; i_q < the_vgts->size; i_q++){ */
/*     Vlong* mc = match_counts[i_q]; */
/*     Gts* q_gts = the_vgts->a[i_q]; */
/*     long q_md_chunk_count = q_gts->md_chunk_count; */
/*     for(long i_2 = 0; i_2 < mc->size; i_2++){ */
/*       long i_m = i_q + i_2; */
/*       //  fprintf(stderr, "i_q, i_2, i_m: %ld %ld %ld \n", i_q, i_2, i_m); */
/*       long matching_chunk_count = mc->a[i_2]; */
      
/*       //   fprintf(stderr, "   %ld %ld \n", mc->size, matching_chunk_count); */
/*       Gts* m_gts = the_vgts->a[i_m]; */
/*       long m_md_chunk_count = m_gts->md_chunk_count; */
/*       //  fprintf(stderr, "md_chunk_counts: %ld %ld \n", q_md_chunk_count, m_md_chunk_count); */
      
/*       double usable_chunk_count = (double)((n_chunks-q_md_chunk_count)*(n_chunks-m_md_chunk_count))/(double)n_chunks; */

/*       if(usable_chunk_count >= min_usable_chunks  &&  matching_chunk_count > min_matching_chunk_fraction*usable_chunk_count){ */
/* 	double matching_chunk_fraction = (double)matching_chunk_count/usable_chunk_count; // fraction matching chunks */
/* 	double est_agmr = 1.0 - pow(matching_chunk_fraction, 1.0/chunk_size); */
/* 	double true_hgmr; */
/* 	double true_agmr = agmr(q_gts, m_gts, &true_hgmr); */
/* 	//	true_agmr_count++; */
/* 	Mci* the_mci = construct_mci(i_q, i_m,  usable_chunk_count, matching_chunk_count, est_agmr, true_agmr, true_hgmr); */
/* 	add_mci_to_vmci(query_vmcis[i_q], the_mci); */
/* 	Mci* the_other_mci = construct_mci(i_m, i_q, usable_chunk_count, matching_chunk_count, est_agmr, true_agmr, true_hgmr); */
/* 	add_mci_to_vmci(query_vmcis[i_m], the_other_mci); */
/*       } */
/*     } */
/*   } */
/*   fprintf(stderr, "time in find_matches_alt: %lf \n", hi_res_time() - start); */
/*   return query_vmcis; */
/* } */


// *****  end of function definitions  *****
