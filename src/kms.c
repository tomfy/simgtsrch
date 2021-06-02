// C version of k-mer search for pairs of similar genotype sets.
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <regex.h>

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
  Pattern_ids** a; // array Pattern_ids
} Chunk_pattern_ids;

// *********************** function declarations ************************************************

long int_power(long base, long power);
char* ipat_to_strpat(long len, long ipat);
long strpat_to_ipat(long len, char* strpat);

// ***** Vlong **********************************************************************************
Vlong* construct_vlong(long cap); // capacity = cap, size = 0
Vlong* construct_vlong_whole_numbers(long size); // initialize to 0,1,2,3,...size-1
void add_long_to_vlong(Vlong* the_vlong, long x);
void shuffle_vlong(Vlong* the_vlong);

// *****  Vstr  *********************************************************************************
Vstr* construct_vstr(long min_size);
void add_string_to_vstr(Vstr* the_vstr, char* str);

// *****  Gts  **********************************************************************************
Gts* construct_gts(long index, char* gtset);
void set_gts_chunk_patterns(Gts* the_gts, Vlong* m_indices, long n_chunks, long k);
char* print_gts(Gts* the_gts);

// *****  Vgts  *********************************************************************************
Vgts* construct_vgts(long min_size);
void add_gts_to_vgts(Vgts* the_vgts, Gts* gts);
void set_vgts_chunk_patterns(Vgts* the_vgts, Vlong* m_indices, long n_chunks, long k);
void populate_chunk_pattern_ids_from_vgts(Vgts* the_vgts, Chunk_pattern_ids* the_cpi);
void print_vgts(Vgts* the_vgts);
void check_gts_indices(Vgts* the_vgts);

// *****  Pattern_ids; indices are patterns; elements are Vlong* of accids having that pattern.
Pattern_ids* construct_pattern_ids(long n_patterns);

// *****  Chunk_pattern_ids; indices are chunk numbers; elements are Pattern_ids*
Chunk_pattern_ids* construct_chunk_pattern_ids(long n_chunks, long n_patterns);
void print_chunk_pattern_ids(Chunk_pattern_ids* the_cpi);

// *****  Gts and Chunk_pattern_ids  ***********
Vlong* find_kmer_match_counts(Gts* the_gts, Chunk_pattern_ids* the_cpi, long n_accessions);

// *************************  end of declarations  **********************************************


// **********************************************************************************************
// ****************************** main **********************************************************
// **********************************************************************************************

int
main(int argc, char *argv[])
{

  // read in genotype matrix file
  FILE *stream;
  char *line = NULL;
  size_t len = 0;
  ssize_t nread;
  long max_number_of_accessions = 1000000;

  long n_chunks = 1000;
  long kmer_size = 5;

  /* for(long i=0; i<27; i++){ */
  /*   char* strpat = ipat_to_strpat(kmer_size, i); */
  /*   long ipat = strpat_to_ipat(kmer_size, strpat); */
  /*   char* strpat2 = ipat_to_strpat(kmer_size, ipat); */
  /*   printf("%ld  %s  %ld  %s\n", i, strpat, ipat, strpat2); */
  /* } */
  /* exit(1); */
  
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <file>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  stream = fopen(argv[1], "r");
  if (stream == NULL) {
    perror("fopen");
    exit(EXIT_FAILURE);
  }

  // *****  read first line 'MARKER followed by marker ids  *****
  nread = getline(&line, &len, stream); // read first line. should have 'MARKER' and marker ids. 
  char mrkr[64];
  sscanf(line, "%s", mrkr);
  if(strcmp(mrkr, "MARKER") !=0){
    printf("Start of first line: %s ; should be MARKER. Bye.\n", mrkr); // should start with 'MARKER'
    exit(EXIT_FAILURE);
  } 
  
  // ***** read rest of file and store genotype sets and ids in Gts objs.  *****
  int init_vgts_size = 50; // 
  Vgts* the_vgts = construct_vgts(init_vgts_size);

  long n_markers;
  long gtsets_count = 0;
  if((nread = getline(&line, &len, stream)) == -1) exit(EXIT_FAILURE); // process first line with genotypes
  char* id = (char*)malloc(100*sizeof(char));
  char* gtstr = (char*)malloc(nread*sizeof(char)); 
  sscanf(line, "%s %s", id, gtstr);
  //  printf("gtsets_count: %ld  id %s;  gtstr: %s\n", gtsets_count, id, gtstr);
  Gts* the_gts = construct_gts(gtsets_count, gtstr);
  n_markers = strlen(gtstr); 
  add_gts_to_vgts(the_vgts, the_gts);
  gtsets_count++;
  
  while ((nread = getline(&line, &len, stream)) != -1) {
    char* id = (char*)malloc(100*sizeof(char));
    char* gtstr = (char*)malloc(nread*sizeof(char)); 
    sscanf(line, "%s %s", id, gtstr);
    //  printf("gtsets_count: %ld  id %s;  gtstr: %s\n", gtsets_count, id, gtstr);
    if(strlen(gtstr) != n_markers) exit(EXIT_FAILURE); // check number of genotypes is same.
    Gts* the_gts = construct_gts(gtsets_count, gtstr);
    //   printf("gts index: %ld\n", the_gts->index); 
    add_gts_to_vgts(the_vgts, the_gts);
    check_gts_indices(the_vgts);
    gtsets_count++;
    if(gtsets_count >= max_number_of_accessions) break;
  }
  //  print_vgts(the_vgts);
  printf("number of genotype sets stored: %ld\n", gtsets_count);
  printf("the_vgts size: %ld\n", the_vgts->size);
  printf("n_markers: %ld\n", n_markers);
  // construct_vlong_whole_numbers(
  /* for(long iiii=0; iiii<the_vgts->size; iiii++){ */
  /*   Gts* a_gts = the_vgts->a[iiii]; */
  /*   printf("after set_vgts... iiii: %ld  a_gts->index: %ld\n", iiii, a_gts->index); */
  /* } */
  //  exit(EXIT_FAILURE);
  Vlong* marker_indices = construct_vlong_whole_numbers(n_markers);
  shuffle_vlong(marker_indices); 

  // void set_vgts_chunk_patterns(Vgts* the_vgts, Vlong* m_indices, long n_chunks, long k){
  if(n_chunks*kmer_size > marker_indices->size){
    n_chunks = marker_indices->size/kmer_size;
    printf("Reducing number of chunks to %ld\n", n_chunks);
  }
  set_vgts_chunk_patterns(the_vgts, marker_indices, n_chunks, kmer_size); 

  /* for(long iiii=0; iiii<the_vgts->size; iiii++){ */
  /*   Gts* a_gts = the_vgts->a[iiii]; */
  /*   printf("after set_vgts... iiii: %ld  a_gts->index: %ld\n", iiii, a_gts->index); */
  /* } */

  printf("after set_vgts_chunk_patterns\n");
  
  long n_patterns = int_power(3, kmer_size); // get 3^kmer_size using integer math.
  printf("n_patterns %ld\n", n_patterns);
  Chunk_pattern_ids* the_cpi = construct_chunk_pattern_ids(n_chunks, n_patterns);
  printf("after construct_chunk_pattern_ids\n");
  
  populate_chunk_pattern_ids_from_vgts(the_vgts, the_cpi);

  //  print_chunk_pattern_ids(the_cpi);
  printf("after populate_...\n");
 

  printf("the_vgts->size: %ld\n", the_vgts->size);
  for(long i_query=0; i_query< the_vgts->size; i_query++){
    Gts* q_gts = the_vgts->a[i_query];
    // Vlong* find_kmer_match_counts(Gts* the_gts, Chunk_pattern_ids* the_cpi, long n_accessions){
    Vlong* kmer_match_counts = find_kmer_match_counts(q_gts, the_cpi, the_vgts->size);
    //  printf("after find_kmer_...\n");
      printf("query accession number: %ld\n", i_query);
    for(long i=i_query; i<the_vgts->size; i++){  
      if(kmer_match_counts->a[i] > 0){
	if(kmer_match_counts->a[i] > n_chunks/5){
	  printf("   matchidx: %ld  match counts: %ld \n", i, kmer_match_counts->a[i]);
	}
      }
    }
  }
  
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
  // printf("in construct_gts. index: %ld  gtset: %s \n", the_gts->index, the_gts->gtset);
  return the_gts;
}

void set_gts_chunk_patterns(Gts* the_gts, Vlong* m_indices, long n_chunks, long k){
  //  printf("top of set_gts_chunk_patterns. n_chunks: %ld  k: %ld  gts index: %ld\n", n_chunks, k, the_gts->index);
 
  Vlong* chunk_pats = construct_vlong(n_chunks); // (Vlong*)malloc(n_chunks*sizeof(Vlong));
  //  printf("after construct_vlong chunk_pats; n_chunks: %ld\n", n_chunks);
  //  printf("m_indices->a[0]: %ld \n", m_indices->a[0]);
  // printf("gtsidx: %ld  gtset: %s\n", the_gts->index, the_gts->gtset);
  for(long i_chunk=0; i_chunk < n_chunks; i_chunk++){
    long i_chunkstart = k*i_chunk;
    long pat = 0;
    long f = 1;
    char* sss = (char*)calloc(k+1,sizeof(char));
    sss[k] = '\0';
    for(long j=0; j < k; j++){
      //  printf("i_chunkstart: %ld  j: %ld m_indices->size: %ld\n", i_chunkstart, j, m_indices->size);
      long idx = m_indices->a[i_chunkstart + j]; // 
      //   printf("i_chunk: %ld ; j: %ld ; idx: %ld  size of gtset:  %ld\n", i_chunk, j, idx, (long)strlen(the_gts->gtset));
      
      char a = the_gts->gtset[idx];
      sss[j] = a;
      long l = (long)a - 48;
     
      if((l>=0) && (l<=2)){
	pat += f*l;
	f*=3;
      }else{
	pat = -1;
	sss = "XXX";
	break;
      }
    }
    //  printf("the_gts->index: %ld   sss: [%s]   ", the_gts->index, sss);
    //  printf("[%s] %ld  [%s]\n", sss, pat, ipat_to_strpat(k, pat));
    //   printf("gts index: %ld  chunk: %ld  pat: %ld \n", the_gts->index, i_chunk, pat);
    add_long_to_vlong(chunk_pats, pat);
    // printf("after add_long_to_vlong\n");
  }
  the_gts->chunk_patterns = chunk_pats;
  // printf("the_gts->index: %ld\n", the_gts->index);
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
 
    // printf("In set_vgts_chunk_patterns. before gts %ld  %ld\n", i, the_vgts->a[i]->index);
  for(long i=0; i < the_vgts->size; i++){
    set_gts_chunk_patterns(the_vgts->a[i], m_indices, n_chunks, k);
    /* for(long iiii=0; iiii<the_vgts->size; iiii++){ */
    /*   Gts* a_gts = the_vgts->a[iiii]; */
    /*   printf("after set_vgts... iiii: %ld  a_gts->index: %ld\n", iiii, a_gts->index); */
    /* } */
    // printf("In set_vgts_chunk_patterns. after gts %ld  %ld\n", i, the_vgts->a[i]->index);
  }
}

void populate_chunk_pattern_ids_from_vgts(Vgts* the_vgts, Chunk_pattern_ids* the_cpi){
  for(long i_gts=0; i_gts<the_vgts->size; i_gts++){
    Gts* the_gts = the_vgts->a[i_gts];
    Vlong* the_chunk_patterns = the_gts->chunk_patterns;
    //     printf("i_gts %ld  chunk_patterns size: %ld\n", i_gts, the_chunk_patterns->size);
    for(long i_chunk=0; i_chunk<the_chunk_patterns->size; i_chunk++){
      long the_pat = the_chunk_patterns->a[i_chunk];
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
	// printf("negative pat: %ld \n", the_pat);
      }
    }
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


// *****  Pattern_ids; indices are patterns; elements are Vlong* of accids having that pattern.

Pattern_ids* construct_pattern_ids(long n_patterns){ // needed size known at construct time, so one param for both cap and size
  Pattern_ids* pat_ids = (Pattern_ids*)malloc(1*sizeof(Pattern_ids));
  pat_ids->capacity = n_patterns;
  pat_ids->size = n_patterns;
  pat_ids->a = (Vlong**)malloc(n_patterns*sizeof(Vlong*));
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
 
Vlong* find_kmer_match_counts(Gts* the_gts, Chunk_pattern_ids* the_cpi, long n_accessions){
  Vlong* chunk_pats = the_gts->chunk_patterns;
  Vlong* accidx_matchcounts = construct_vlong(n_accessions);
  // printf("n chunks: %ld\n", chunk_pats->size);
  for(long i_chunk=0; i_chunk < chunk_pats->size; i_chunk++){
    long the_pat = chunk_pats->a[i_chunk];
    //   printf("i_chunk: %ld  the_pat: %ld\n", i_chunk, the_pat);
    if(the_pat >= 0){
      Vlong* kmer_match_idxs = the_cpi->a[i_chunk]->a[the_pat];
      for(long i=0; i<kmer_match_idxs->size; i++){
	long accidx = kmer_match_idxs->a[i]; // index of one of the accessions matching on this chunk
	accidx_matchcounts->a[accidx]++;
      }
    }
  }
  return accidx_matchcounts; 
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

