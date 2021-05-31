// C version of k-mer search for pairs of similar genotype sets.
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <regex.h>

long* shuffled(int n);

//************************************************************

typedef struct{
  long capacity; // allocated size
  long size; // number of elements
  long* a; // array
} Vlong;

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
  // if necessary, resize w realloc
  if(n == cap){
    cap *= 2;
    the_vlong->a = (long*)realloc(the_vlong->a, cap*sizeof(char*));
  }
  the_vlong->a[n] = x;
  the_vlong->size++;
}

// ************************************************************

typedef struct{
  long capacity; // allocated size
  long size; // number of elements
  char** a; // array of strings
} Vstr;

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
  }
  the_vstr->a[n] = str;
  the_vstr->size++;
}

//****************************************************************

typedef struct{ // genotype set
  long index;
  char* gtset;
  long* chunk_patterns;
} Gts;

Gts* construct_gts(long index, char* gtset){
  Gts* the_gts = (Gts*)malloc(1*sizeof(Gts));
  the_gts->index = index;
  the_gts->gtset = gtset;
  the_gts->chunk_patterns = NULL;
  return the_gts;
}

void gts_chunk_patterns(Gts* the_gts, Vlong* m_indices, long n_chunks, long k){
  long* chunk_pats = (long*)malloc(n_chunks*sizeof(long));
  long pat = 0;
  long f = 1;
  for(int i_chunk=0; i_chunk < n_chunks; i_chunk++){
    for(int j=0; j < k; j++){
      int jdx = (i_chunk*k+j) % m_indices->size;
      int idx = m_indices->a[jdx];
      char a = the_gts->gtset[idx] - 48;
      if((a>=0) && (a<=2)){
	pat += f*a;
	f*=3;
      }else{
	pat = -1;
	break;
      }
    }
    chunk_pats[i_chunk] = pat;
  }
  the_gts->chunk_patterns = chunk_pats;
}

char* print_gts(Gts* the_gts){
  printf("Gts index: %ld  gtset: %s\n", the_gts->index, the_gts->gtset);
}

// ****************************************************************

typedef struct{
  long capacity; // allocated size
  long size; // number of elements
  Gts** a; // array of Gts* (genotype set)
} Vgts;

Vgts* construct_vgts(long min_size){
  Vgts* the_vgts = (Vgts*)malloc(1*sizeof(Vgts));
  the_vgts->capacity = min_size;
  the_vgts->size = 0;
  the_vgts->a = (Gts**)malloc(min_size*sizeof(Gts*));
  printf("returning from construct_vgts\n");
  return the_vgts;
}

void add_gts_to_vgts(Vgts* the_vgts, Gts* gts){
  long cap = the_vgts->capacity;
  long n = the_vgts->size;
  // if necessary, resize w realloc
  if(n == cap){
    cap *= 2;
    the_vgts->a = (Gts**)realloc(the_vgts->a, cap*sizeof(Gts*));
  }
  the_vgts->a[n] = gts;
  the_vgts->size++;
}

void print_vgts(Vgts* the_vgts){
  for(int i=0; i<the_vgts->size; i++){
    print_gts(the_vgts->a[i]);
  }
}


// **************** main ***********************

int
main(int argc, char *argv[])
{

  // read in genotype matrix file
  FILE *stream;
  char *line = NULL;
  size_t len = 0;
  ssize_t nread;
  long max_number_of_accessions = 1000000;

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
  nread = getline(&line, &len, stream); // read first line. 
  char mrkr[64];
  sscanf(line, "%s", mrkr);
  if(strcmp(mrkr, "MARKER") !=0){
    printf("Start of first line: %s ; should be MARKER. Bye.\n", mrkr); // should start with 'MARKER'
    exit(EXIT_FAILURE);
  }
  
  // ***** read rest of file and store genotype sets and ids in Gts objs.  *****
  int init_vgts_size = 50; // 
  Vgts* the_vgts = construct_vgts(init_vgts_size); 

  int gtsets_count = 0;
  while ((nread = getline(&line, &len, stream)) != -1) {
    char* id = (char*)malloc(100*sizeof(char));
    char* gts = (char*)malloc(nread*sizeof(char)); 
    sscanf(line, "%s %s", id, gts);
    Gts* the_gts = construct_gts(gtsets_count, gts);

    add_gts_to_vgts(the_vgts, the_gts);
    
    gtsets_count++;
    if(gtsets_count >= max_number_of_accessions) break;
  }
  printf("number of genotype sets stored: %d\n", gtsets_count);
  printf("the_vstr size: %ld\n", the_vgts->size);
  // construct_vlong_whole_numbers(
 
  print_vgts(the_vgts);
  

  printf("%ld \n", (long)RAND_MAX);
  
  fclose(stream);
  exit(EXIT_SUCCESS);
}


//char* get_nonwhitespace_string(FILE* stream){
//  char* s = (char*)malloc(50*sizeof(char));
  
  
// *********************************************

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


