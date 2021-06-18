// C version of k-mer search for pairs of similar genotype sets.
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>

long do_dbl_md_chunk_counts = 0;

//***********************************************************************************************
// **************  typedefs  ********************************************************************
typedef struct{
  long capacity; // allocated size
  long size; // number of elements
  long* a; // array
} Vlong;

typedef struct{
  long capacity; // allocated size
  long size; // number of elements
  char** a; // array of strings
} Vstr;

typedef struct{ // genotype set
  long index;
  char* gtset;
  Vlong* chunk_patterns;
  long md_chunk_count; // number of chunks with missing data (in at least one gt)
  long missing_count; // number of gts with missing data.
} Gts;

typedef struct{
  long capacity; // allocated size
  long size; // number of elements
  Gts** a; // array of Gts* (genotype set)
} Vgts;

typedef struct{
  long capacity;
  long size;
  Vlong** a; // array of Vlong*
} Pattern_ids;

typedef struct{
  long capacity;
  long size;
  long n_patterns;
  Pattern_ids** a; // array Pattern_ids
} Chunk_pattern_ids;

typedef struct{
  long query_index;
  long match_index;
  long n_matching_chunks;
  //  long est_usable_chunk_count;
  //  long usable_chunk_count;
  double est_matching_chunk_fraction;
  double matching_chunk_fraction;
} Mci; // 'Mci = Matching chunk info'

typedef struct{
  long capacity;
  long size;
  Mci** a;
} Vmci; 


// *********************** function declarations ************************************************

long int_power(long base, long power);
char* ipat_to_strpat(long len, long ipat);
long strpat_to_ipat(long len, char* strpat);
double agmr(Gts* gts1, Gts* gts2);
double hi_res_time(void);

// ***** Vlong **********************************************************************************
Vlong* construct_vlong(long cap); // capacity = cap, size = 0
Vlong* construct_vlong_whole_numbers(long size); // initialize to 0,1,2,3,...size-1
void add_long_to_vlong(Vlong* the_vlong, long x);
void shuffle_vlong(Vlong* the_vlong);
void free_vlong(Vlong* the_vlong);

// *****  Vstr  *********************************************************************************
Vstr* construct_vstr(long min_size);
void add_string_to_vstr(Vstr* the_vstr, char* str);

// *****  Gts  **********************************************************************************
Gts* construct_gts(long index, char* gtset);
long set_gts_chunk_patterns(Gts* the_gts, Vlong* m_indices, long n_chunks, long k);
char* print_gts(Gts* the_gts);

// *****  Vgts  *********************************************************************************
Vgts* construct_vgts(long min_size);
void add_gts_to_vgts(Vgts* the_vgts, Gts* gts);
void set_vgts_chunk_patterns(Vgts* the_vgts, Vlong* m_indices, long n_chunks, long k);
void populate_chunk_pattern_ids_from_vgts(Vgts* the_vgts, Chunk_pattern_ids* the_cpi);
void print_vgts(Vgts* the_vgts);
void check_gts_indices(Vgts* the_vgts);

// *****  Mci  ********
Mci* construct_mci(long qidx, long midx, long n_matching_chunks,
		   // long est_n_usable_chunks, long n_usable_chunks);
		   double est_matching_chunk_fraction, double matching_chunk_fraction);
// *****  Vmci  *********************************************************************************
Vmci* construct_vmci(long init_size);
void add_mci_to_vmci(Vmci* the_vmci, Mci* the_mci); 

// *****  Pattern_ids; indices are patterns; elements are Vlong* of accids having that pattern.
Pattern_ids* construct_pattern_ids(long n_patterns);

// *****  Chunk_pattern_ids; indices are chunk numbers; elements are Pattern_ids*
Chunk_pattern_ids* construct_chunk_pattern_ids(long n_chunks, long n_patterns);
void print_chunk_pattern_ids(Chunk_pattern_ids* the_cpi);

// *****  Gts and Chunk_pattern_ids  ***********
Vlong* find_chunk_match_counts(Gts* the_gts, Chunk_pattern_ids* the_cpi, long n_accessions, Vlong* accidx_doublemissingdatacounts);

// *************************  end of declarations  **********************************************


// **********************************************************************************************
// ****************************** main **********************************************************
// **********************************************************************************************

int
main(int argc, char *argv[])
{
  long n_chunks = 1000; // number of chunks (default = 1000)
  long chunk_size = 5; // number of genotype in each chunk (default = 5)
  long max_number_of_accessions = 1000000;

  
  // ***** process command line *****
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <file> options\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  printf("argv[1]: [%s] \n", argv[1]);
  char* input_filename = NULL;
  FILE *stream = NULL;
  //  stream = fopen(argv[1], "r");
    
  int c;
  while((c = getopt(argc, argv, "i:n:k:m")) != -1){
    printf("c: %c %s %d\n", c, optarg, optind);
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
  printf("optind: %d argc: %d\n", optind, argc);
  if(optind < argc){
    printf("Non-option arguments. Bye.\n");
  }
  /* for (int index = optind; index < argc; index++){ */
  /*   printf ("Non-option argument %s\n", argv[index]); */
  /* } */
  //return 0;

  if(input_filename == NULL){
    perror("must specify input filename: -i <filename>");
    exit(EXIT_FAILURE);
  }

  printf("input file: %s  n_chunks: %ld  chunk size: %ld  get exact md chunk counts?: %ld \n",
	 input_filename, n_chunks, chunk_size, do_dbl_md_chunk_counts);
  //exit(1);
  // *****  done processing command line  *****

  double start;
  start = hi_res_time();
  
  // *****  read in genotype matrix file  *****
 
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
  int init_vgts_size = 100; // 
  Vgts* the_vgts = construct_vgts(init_vgts_size);

  long n_markers;
  long gtsets_count = 0;
  Vstr* accession_ids = construct_vstr(100);

  if((nread = getline(&line, &len, stream)) != -1){ // process first line with genotypes 
    char* id = (char*)malloc(100*sizeof(char));
    char* gtstr = (char*)malloc(nread*sizeof(char)); 
    sscanf(line, "%s %s", id, gtstr);
    add_string_to_vstr(accession_ids, id);
    Gts* the_gts = construct_gts(gtsets_count, gtstr);
    n_markers = strlen(gtstr); 
    add_gts_to_vgts(the_vgts, the_gts);
    gtsets_count++;
  }else{
    exit(EXIT_FAILURE);
  }

  while ((nread = getline(&line, &len, stream)) != -1) {
    char* id = (char*)malloc(100*sizeof(char));
    char* gtstr = (char*)malloc(nread*sizeof(char)); 
    sscanf(line, "%s %s", id, gtstr);
    add_string_to_vstr(accession_ids, id);
    if(strlen(gtstr) != n_markers) exit(EXIT_FAILURE); // check number of genotypes is same.
    Gts* the_gts = construct_gts(gtsets_count, gtstr);
    add_gts_to_vgts(the_vgts, the_gts);
    check_gts_indices(the_vgts);
    gtsets_count++;
    if(gtsets_count >= max_number_of_accessions) break;
  }
  printf("# done reading in genotypes data.\n");
  printf("# time to read input: %12.6f\n", hi_res_time() - start);
  printf("# number of genotype sets stored: %ld  %ld\n", gtsets_count, the_vgts->size);
  printf("# n_markers: %ld\n", n_markers);

  // *****  done reading in genotype sets for all accessions  **********

  start = hi_res_time();
  Vlong* marker_indices = construct_vlong_whole_numbers(n_markers);
  shuffle_vlong(marker_indices); 

  if(n_chunks*chunk_size > marker_indices->size){
    n_chunks = marker_indices->size/chunk_size;
    printf("Reducing number of chunks to %ld\n", n_chunks);
  }
  set_vgts_chunk_patterns(the_vgts, marker_indices, n_chunks, chunk_size); 
  printf("after set_vgts_chunk_patterns\n");
  printf("time to shuffle indices and set the genotype set chunk patterns: %12.6f\n",  hi_res_time() - start);
  // 

  start = hi_res_time();
  long n_patterns = int_power(3, chunk_size); // get 3^chunk_size using integer math.
  printf("n_patterns %ld\n", n_patterns);
  Chunk_pattern_ids* the_cpi = construct_chunk_pattern_ids(n_chunks, n_patterns);
  populate_chunk_pattern_ids_from_vgts(the_vgts, the_cpi);
  printf("after construct and populate_chunk_pattern_ids_from_vgts\n");
  printf("time to construct chunk_pattern_ids data structure: %12.6f\n", hi_res_time() - start);

  start = hi_res_time();
  if(0){
    long true_agmr_count = 0;
    for(long i_query=0; i_query< the_vgts->size; i_query++){
      Gts* q_gts = the_vgts->a[i_query];
      long q_md_gt_count = q_gts->missing_count;
      long q_md_chunk_count = q_gts->md_chunk_count;

      Vlong* accidx_dblmdcounts = construct_vlong(the_vgts->size);
      Vlong* chunk_match_counts = find_chunk_match_counts(q_gts, the_cpi, the_vgts->size, accidx_dblmdcounts);
      for(long i_match=i_query; i_match<the_vgts->size; i_match++){
	Gts* match_gts = the_vgts->a[i_match];
	long matching_chunk_count = chunk_match_counts->a[i_match];
	long match_md_gt_count = match_gts->missing_count;
	long match_md_chunk_count = match_gts->md_chunk_count;

	// using estimated number of usable chunks
	long est_md_gt_count = (q_md_gt_count >= match_md_gt_count)? q_md_gt_count : match_md_gt_count;
	double est_usable_chunk_count =
	  n_chunks*pow((1.0 - (double)est_md_gt_count/(double)n_markers), chunk_size);

	// using correct number of usable pairs (if was calculated)
	long usable_chunk_count = est_usable_chunk_count;
	if(do_dbl_md_chunk_counts){
	  long dbl_md_count = accidx_dblmdcounts->a[i_match]; // correct double missing data chunk count (if done in find_...
	  usable_chunk_count = n_chunks - (q_md_chunk_count + match_md_chunk_count - dbl_md_count); // actual number of chunks with no missing data in either query or match
	}
      
	if(usable_chunk_count > 50  &&  matching_chunk_count > 0.2*usable_chunk_count){
	  double est_matching_chunk_fraction = (double)matching_chunk_count/(double)est_usable_chunk_count;
	  double matching_chunk_fraction = (double)matching_chunk_count/(double)usable_chunk_count; // fraction matching chunks
	  double est_agmr_1 = 1.0 - pow(est_matching_chunk_fraction, 1.0/chunk_size);
	  double est_agmr_2 = 1.0 - pow(matching_chunk_fraction, 1.0/chunk_size);	  
	  fprintf(stderr, "%ld %ld %ld  %9.6f %9.6f %9.6f %9.6f %9.6f\n",
		  i_query, i_match, matching_chunk_count, 
		  est_matching_chunk_fraction, est_agmr_1,
		  matching_chunk_fraction, est_agmr_2,
		  agmr(q_gts, match_gts));
	  true_agmr_count++;
	}
      }
      // free_vlong(chunk_match_counts);
    }
    printf("time to find candidate matches and %ld true agmrs: %12.6f\n", true_agmr_count, hi_res_time() - start);
  }else{ // new way - get all good candidates for matches to all queries, and store them. Then 
    Vmci** query_vmcis = (Vmci**)malloc(the_vgts->size * sizeof(Vmci*)); // 
    
    for(long i_query=0; i_query< the_vgts->size; i_query++){
      Gts* q_gts = the_vgts->a[i_query];
      long q_md_gt_count = q_gts->missing_count;
      long q_md_chunk_count = q_gts->md_chunk_count;

      Vlong* accidx_dblmdcounts = construct_vlong(the_vgts->size);
      Vlong* chunk_match_counts = find_chunk_match_counts(q_gts, the_cpi, the_vgts->size, accidx_dblmdcounts);
      Vmci* the_vmci = construct_vmci(20);
      for(long i_match=i_query; i_match<the_vgts->size; i_match++){

	Gts* match_gts = the_vgts->a[i_match];
	long matching_chunk_count = chunk_match_counts->a[i_match];
	long match_md_gt_count = match_gts->missing_count;
	long match_md_chunk_count = match_gts->md_chunk_count;

	// using estimated number of usable chunks
	long est_md_gt_count = (q_md_gt_count >= match_md_gt_count)? q_md_gt_count : match_md_gt_count;
	double est_usable_chunk_count =
	  n_chunks*pow((1.0 - (double)est_md_gt_count/(double)n_markers), chunk_size);

	//   double est_agmr_1 = 1.0 - pow(est_matching_chunk_fraction, 1.0/chunk_size);
      
	// using correct number of usable pairs (if it was calculated)
	long usable_chunk_count = est_usable_chunk_count;
	if(do_dbl_md_chunk_counts){
	  long dbl_md_count = accidx_dblmdcounts->a[i_match]; // correct double missing data chunk count (if done in find_...
	  usable_chunk_count = n_chunks - (q_md_chunk_count + match_md_chunk_count - dbl_md_count); // actual number of chunks with no missing data in either query or match
	
	}
	//	if(matching_chunk_count > 50) printf("matching chunk count: %ld  usable_chunk_count %ld \n", matching_chunk_count, usable_chunk_count);
	if(usable_chunk_count > 50  &&  matching_chunk_count > 0.2*usable_chunk_count){
	  double est_matching_chunk_fraction = (double)matching_chunk_count/(double)est_usable_chunk_count;
	  double matching_chunk_fraction = (double)matching_chunk_count/(double)usable_chunk_count; // fraction matching chunks
	  Mci* the_mci = construct_mci(i_query, i_match, matching_chunk_count, est_matching_chunk_fraction, matching_chunk_fraction);
	  //est_usable_chunk_count, usable_chunk_count);
	  add_mci_to_vmci(the_vmci, the_mci);
	}
	//	printf("query: %ld   n good candidates:  %ld \n", i_query, the_vmci->size);
	query_vmcis[i_query] = the_vmci;
      }
      // free_vlong(accidx_dblmdcounts);
      // free_vlong(chunk_match_counts);
    }
    printf("time to find candidate matches and store good candidates: %12.6f\n", hi_res_time() - start);

    start = hi_res_time();
    long true_agmr_count = 0;
    for(long i_q=0; i_q<the_vgts->size; i_q++){
      Vmci* the_vmci = query_vmcis[i_q];
      Gts* q_gts = the_vgts->a[i_q];
      //  printf("i_q: %ld  the_vmci->size:  %ld \n", i_q, the_vmci->size);
      for(long i_m=0; i_m < the_vmci->size; i_m++){
	Mci* the_mci = the_vmci->a[i_m];
	
	Gts* match_gts = the_vgts->a[the_mci->match_index];
	double est_agmr_1 = 1.0 - pow(the_mci->est_matching_chunk_fraction, 1.0/chunk_size);
	double est_agmr_2 = 1.0 - pow(the_mci->matching_chunk_fraction, 1.0/chunk_size);
	double true_agmr = agmr(q_gts, match_gts);
	true_agmr_count++;
	fprintf(stderr, "%ld %ld  %ld  %9.6f %9.6f %9.6f %9.6f %9.6f\n", 
		i_q, i_m,
		the_mci->n_matching_chunks,
		the_mci->est_matching_chunk_fraction, est_agmr_1,
		the_mci->matching_chunk_fraction, est_agmr_2,
		true_agmr
		);
      }
    }
    printf("time to find %ld true agmrs: %12.6f\n", true_agmr_count, hi_res_time() - start);
  }
  printf("done getting match candidates for all queries, and true agmrs for good candidates.\n");
  // getchar();
  
  fclose(stream);
  exit(EXIT_SUCCESS);
}
  
// **********************  end of main  *********************************************************



// *******************  function definitions  ***************************************************

// *****  Vlong  *****

Vlong* construct_vlong(long cap){
  Vlong* the_vlong = (Vlong*)malloc(1*sizeof(Vlong));
  the_vlong->capacity = cap;
  the_vlong->size = 0;
  the_vlong->a = (long*)malloc(cap*sizeof(char*));
  return the_vlong;
}

Vlong* construct_vlong_whole_numbers(long size){ // initialize to 0,1,2,3,...size-1
  Vlong* the_vlong = (Vlong*)malloc(1*sizeof(Vlong));
  the_vlong->capacity = size;
  the_vlong->size = size;
  the_vlong->a = (long*)malloc(size*sizeof(char*));
  for(int i=0; i<size; i++){
    the_vlong->a[i] = i;
  }
  return the_vlong;
}

void add_long_to_vlong(Vlong* the_vlong, long x){
  long cap = the_vlong->capacity;
  long n = the_vlong->size;
  //  printf("cap n: %ld %ld\n", cap, n);
  // if necessary, resize w realloc
  if(n == cap){
    cap *= 2;
    the_vlong->a = (long*)realloc(the_vlong->a, cap*sizeof(long));
    //    printf("realloc in add_long_to_vlong. new cap: %ld\n", cap);
    the_vlong->capacity = cap;
  }
  //  printf("after realloc. cap: %ld \n", cap);
  the_vlong->a[n] = x;
  //  printf("after assignment to a[n]\n");
  the_vlong->size++;
  // printf("about to return from add_long_to_vlong. size %ld\n", the_vlong->size);
}

void shuffle_vlong(Vlong* the_vlong){
  long n = the_vlong->size;
  for(long i=0; i<n-1; i++){ // shuffle
    int j = i+1 + (long)((n-1-i)*(double)rand()/RAND_MAX); // get an integer in the range [i+1,n-1]
    long tmp = the_vlong->a[j];
    the_vlong->a[j] = the_vlong->a[i];
    the_vlong->a[i] = tmp;
  }
}

void free_vlong(Vlong* the_vlong){
  free(the_vlong->a);
  free(the_vlong);
}

// *****  Vstr  ***************************************************

Vstr* construct_vstr(long min_size){
  Vstr* the_vstr = (Vstr*)malloc(1*sizeof(Vstr));
  the_vstr->capacity = min_size;
  the_vstr->size = 0;
  the_vstr->a = (char**)malloc(min_size*sizeof(char*));
  printf("returning from construct_vstr\n");
  return the_vstr;
}

void add_string_to_vstr(Vstr* the_vstr, char* str){
  long cap = the_vstr->capacity;
  long n = the_vstr->size;
  // if necessary, resize w realloc
  if(n == cap){
    cap *= 2;
    the_vstr->a = (char**)realloc(the_vstr->a, cap*sizeof(char*));
    the_vstr->capacity = cap;
  }
  the_vstr->a[n] = str;
  the_vstr->size++;
}


// *****  Gts  ****************************************************

Gts* construct_gts(long index, char* gtset){
  Gts* the_gts = (Gts*)malloc(1*sizeof(Gts));
  the_gts->index = index;
  the_gts->gtset = gtset;
  the_gts->chunk_patterns = NULL;
  the_gts->md_chunk_count = 0;
  long missing_count = 0;
  for(long i=0; ; i++){
    char a = the_gts->gtset[i];
    if(a == '\0') break;
    if((a == '0') || (a == '1') || (a == '2')){
      // do nothing
    }else{
      missing_count++;
    }
  }
  the_gts->missing_count = missing_count;
  //  printf("In construct_gts. missing_count: %ld \n", the_gts->missing_count);
  return the_gts;
}

long set_gts_chunk_patterns(Gts* the_gts, Vlong* m_indices, long n_chunks, long k){
  //  printf("top of set_gts_chunk_patterns. n_chunks: %ld  k: %ld  gts index: %ld\n", n_chunks, k, the_gts->index);
  long gts_mdchunk_count = 0;
  long n_patterns = int_power(3, k); // 3^k, the number of valid patterns, also there is a 'pattern' for missing data, making 3^k + 1 in all
  Vlong* chunk_pats = construct_vlong(n_chunks); // (Vlong*)malloc(n_chunks*sizeof(Vlong));
  //  printf("after construct_vlong chunk_pats; n_chunks: %ld\n", n_chunks);
  //  printf("m_indices->a[0]: %ld \n", m_indices->a[0]);
  // printf("gtsidx: %ld  gtset: %s\n", the_gts->index, the_gts->gtset);
  for(long i_chunk=0; i_chunk < n_chunks; i_chunk++){
    long i_chunkstart = k*i_chunk;
    long i_pat = 0;
    long f = 1;

    // loop over characters in the chunk and construct an corresponding long index, in range [0..3^k] (3^k is the index for a chunk with any missing data)
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
    } // loop over the k chars in a chunk.
    //  printf("the_gts->index: %ld   sss: [%s]   ", the_gts->index, sss);
    //  printf("[%s] %ld  [%s]\n", sss, pat, ipat_to_strpat(k, pat));
    //   printf("gts index: %ld  chunk: %ld  pat: %ld \n", the_gts->index, i_chunk, pat);
    add_long_to_vlong(chunk_pats, i_pat);
    // printf("after add_long_to_vlong\n");
  } // loop over chunks.
  the_gts->chunk_patterns = chunk_pats;
  the_gts->md_chunk_count = gts_mdchunk_count;
  return gts_mdchunk_count;
}

char* print_gts(Gts* the_gts){
  printf("Gts index: %ld  gtset: %s\n", the_gts->index, the_gts->gtset);
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
  // if necessary, resize w realloc
  if(n == cap){
    cap *= 2;
    the_vgts->a = (Gts**)realloc(the_vgts->a, cap*sizeof(Gts*));
    the_vgts->capacity = cap;
  }
  the_vgts->a[n] = gts;
  the_vgts->size++;
}

void set_vgts_chunk_patterns(Vgts* the_vgts, Vlong* m_indices, long n_chunks, long k){
  long total_mdchunk_count = 0;
  // printf("In set_vgts_chunk_patterns. before gts %ld  %ld\n", i, the_vgts->a[i]->index);
  for(long i=0; i < the_vgts->size; i++){
    long mdchcount = set_gts_chunk_patterns(the_vgts->a[i], m_indices, n_chunks, k);
    total_mdchunk_count += mdchcount;
    /* for(long iiii=0; iiii<the_vgts->size; iiii++){ */
    /*   Gts* a_gts = the_vgts->a[iiii]; */
    /*   printf("after set_vgts... iiii: %ld  a_gts->index: %ld\n", iiii, a_gts->index); */
    /* } */
    // printf("In set_vgts_chunk_patterns. after gts %ld  %ld\n", i, the_vgts->a[i]->index);
    //    printf("gts idx: %ld  mdchunk_count: this gts: %ld  cumesofar %ld \n", i, mdchcount, total_mdchunk_count);
  }
}

void populate_chunk_pattern_ids_from_vgts(Vgts* the_vgts, Chunk_pattern_ids* the_cpi){
  long n_patterns = the_cpi->n_patterns;
  for(long i_gts=0; i_gts<the_vgts->size; i_gts++){
    Gts* the_gts = the_vgts->a[i_gts];
    Vlong* the_chunk_patterns = the_gts->chunk_patterns; // the gt patterns (longs) occurring in each chunk of this gts 
    //     printf("i_gts %ld  chunk_patterns size: %ld\n", i_gts, the_chunk_patterns->size);
    long mdcount = 0;
    for(long i=0; i<the_chunk_patterns->size; i++){
      if(the_chunk_patterns->a[i] == n_patterns){ mdcount++; }
    }
    //  printf("accidx: %ld  mdcount %ld \n", the_gts->index, mdcount);
    
    for(long i_chunk=0; i_chunk<the_chunk_patterns->size; i_chunk++){
      long the_pat = the_chunk_patterns->a[i_chunk];
      /* if(the_pat == n_patterns){ */
      /* 	printf("the_pat: %ld  i_chunk %ld  accidx: %ld\n", the_pat, i_chunk, i_gts);  */
      /* } */
      if(the_pat >= 0){
	Vlong* the_accidxs = the_cpi->a[i_chunk]->a[the_pat];
	/* printf("i_gts: %ld i_chunk: %ld  the_pat: %ld      the_accidxs cap: %ld \n", */
	/*        i_gts, i_chunk, the_pat, the_accidxs->capacity); */
	//     printf("the_gts->index: %ld\n", the_gts->index);
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
    //   printf("Anumber of missingdata accs for chunk %ld   is  %ld; totalsofar: %ld\n", i, chunk_md_count, total_mdchunk_count);
  }
  printf("bottom of populate. total_mdchunk_count: %ld\n", total_mdchunk_count);
  /* for(long i=0; i<the_cpi->size; i++){ */
  /*   printf("number of missingdata accs for chunk %ld   is  %ld\n", i, the_cpi->a[i]->a[n_patterns]->size); */
  /* } */
  // exit(1);
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

// *****  Chunk_pattern_ids; indices are chunk numbers; elements are Pattern_ids*

Chunk_pattern_ids* construct_chunk_pattern_ids(long n_chunks, long n_patterns){ // needed size known at construct time, so one param for both cap and size
  Chunk_pattern_ids* chunk_pat_ids = (Chunk_pattern_ids*)malloc(1*sizeof(Chunk_pattern_ids));
  chunk_pat_ids->capacity = n_chunks;
  chunk_pat_ids->size = n_chunks;
  chunk_pat_ids->n_patterns = n_patterns;
  chunk_pat_ids->a = (Pattern_ids**)malloc(n_chunks*sizeof(Pattern_ids*));
  for(int i=0; i< chunk_pat_ids->size; i++){
    chunk_pat_ids->a[i] = construct_pattern_ids(n_patterns);
  }
  return chunk_pat_ids;
}

void print_chunk_pattern_ids(Chunk_pattern_ids* the_cpi){
  for(long i_chunk=0; i_chunk<the_cpi->size; i_chunk++){
    printf("i_chunk: %ld\n", i_chunk);
    Pattern_ids* the_pi = the_cpi->a[i_chunk];
    for(long i_pat=0; i_pat<the_pi->size; i_pat++){
      printf("  i_pat: %ld \n", i_pat); //  strpat: %s\n", i_pat, ipat_to_strpat(3, i_pat));
      Vlong* the_idxs = the_pi->a[i_pat];
      //   printf("cpi->size: %ld  i_chunk: %ld   pi->size: %ld  nmatches:  %ld  ipat: %ld strpat: %s\n", the_cpi->size, i, the_pi->size, the_idxs->size, p, ipat_to_strpat(3,p));
      if(the_idxs->size > 0){
	//	printf("i_chunk: %ld   pat: %ld \n", i, p);
	printf("     matches: ");
	for(long ii=0; ii<the_idxs->size; ii++){
	  printf("%ld  ", the_idxs->a[ii]);
	  //	  printf("  imatch: %ld  idx: %ld;", ii, the_idxs->a[ii]);
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
  Vlong* accidx_matchcounts = construct_vlong(n_accessions);
  // Vlong*
 
  // printf("n chunks: %ld\n", chunk_pats->size);
  for(long i_chunk=0; i_chunk < chunk_pats->size; i_chunk++){
    long the_pat = chunk_pats->a[i_chunk];
    //   printf("i_chunk: %ld  the_pat: %ld\n", i_chunk, the_pat);
  
    Vlong* chunk_match_idxs = the_cpi->a[i_chunk]->a[the_pat]; // array of indices of the matches to this chunk & pat
    if(the_pat == n_patterns){ // missing data in this chunk
      /* printf("accidx: %ld  i_chunk: %ld  the_pat: %ld n_patterns %ld\n", the_gts->index, i_chunk, the_pat, n_patterns); */
      /* char achar = getchar(); */
      //  printf("Xnumber of missing data accs, this chunk: %ld \n", chunk_match_idxs->size);
      if(! do_dbl_md_chunk_counts) continue; // control whether to do the double missing data chunks.
      for(long i=0; i<chunk_match_idxs->size; i++){
	long accidx = chunk_match_idxs->a[i]; // index of one of the accessions matching on this chunk
	//	printf("accidx %ld    i: %ld  accidx: %ld \n", accidx, i, accidx);
	accidx_dbl_md_counts->a[accidx]++;
      }
    }else{ // the_pat = 0..n_patterns-1 (good data)
      for(long i=0; i<chunk_match_idxs->size; i++){
	long accidx = chunk_match_idxs->a[i]; // index of one of the accessions matching on this chunk
	accidx_matchcounts->a[accidx]++;
      }
    }
  }
  /* printf("query: %ld\n", the_gts->index); */
  /* printf("%ld %ld\n", accidx_dbl_md_counts->size, accidx_matchcounts->size); */
  /* for(long i=0; i< accidx_dbl_md_counts->size; i++){ */
  /*   printf("  ABC  %ld   %ld \n", i, accidx_matchcounts->a[i]); */
  /*   printf("    DEF  %ld %ld \n", i, accidx_dbl_md_counts->a[i]); */
  /* } */
  return accidx_matchcounts; 
}

// ***** Mci  *****

//Mci* construct_mci(long qidx, long midx, long n_matching_chunks, long est_n_usable_chunks, long n_usable_chunks){
Mci* construct_mci(long qidx, long midx, long n_matching_chunks, double est_matching_chunk_fraction, double matching_chunk_fraction){
  Mci* the_mci = (Mci*)malloc(1*sizeof(Mci));
  the_mci->query_index = qidx;
  the_mci->match_index = midx;
  the_mci->n_matching_chunks = n_matching_chunks;
  the_mci->est_matching_chunk_fraction = est_matching_chunk_fraction;
  the_mci->matching_chunk_fraction = matching_chunk_fraction;

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

// *********************************************

long int_power(long base, long power){ // calculate base^power using integer math.
  long result = 1;
  for(int i=0; i<power; i++){
    result *= base;
  }
  return result;
}

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
    //  printf("i_chunkstart: %ld  j: %ld m_indices->size: %ld\n", i_chunkstart, j, m_indices->size);
    //long idx = m_indices->a[i_chunkstart + j];
    //   printf("i_chunk: %ld ; j: %ld ; idx: %ld  size of gtset:  %ld\n", i_chunk, j, idx, (long)strlen(the_gts->gtset));
      
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
  

// ***** end of function definitions  *****


// *****  unused  *****

/*
  long* shuffled(int n);

  long* shuffled(int n){ // return an n-element array with number 0 through n-1 in random order.
  
  long* a = (long*)malloc(n*sizeof(long)); // allocate
  for(int i=0; i<n; i++){ // initialize to 0..n-1
  a[i] = i;
  }
  for(long i=0; i<n-1; i++){ // shuffle
  int j = i+1 + (long)((n-1-i)*(double)rand()/RAND_MAX); // Use a better rng here!
  long tmp = a[j];
  a[j] = a[i];
  a[i] = tmp;
  }
  return a;
  }
*/
