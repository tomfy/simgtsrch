// *****  Implementation of functions declared in vect.h  *****
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "vect.h"


// *****  Vlong  **********************************************

Vlong* construct_vlong(long cap){
  Vlong* the_vlong = (Vlong*)malloc(1*sizeof(Vlong));
  the_vlong->capacity = cap;
  the_vlong->size = 0;
  the_vlong->a = (long*)malloc(cap*sizeof(char*));
  return the_vlong;
}

Vlong* construct_vlong_from_array(long size, long* array){ // construct from (dynamically allocated) array of known size.
  Vlong* the_vlong = (Vlong*)malloc(1*sizeof(Vlong));
  the_vlong->capacity = size;
  the_vlong->size = size;
  the_vlong->a = array; // array must have been allocated dynamically.
  return the_vlong;
}

Vlong* construct_vlong_zeroes(long size){
  long* zeroes = (long*)calloc(size, sizeof(long));
  Vlong* the_vlong = (Vlong*)malloc(sizeof(Vlong));
  the_vlong->capacity = size;
  the_vlong->size = size;
  the_vlong->a = zeroes;
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

Vstr* construct_vstr(long cap){
  Vstr* the_vstr = (Vstr*)malloc(1*sizeof(Vstr));
  the_vstr->capacity = cap;
  the_vstr->size = 0;
  the_vstr->a = (char**)malloc(cap*sizeof(char*));
  // printf("returning from construct_vstr\n");
  return the_vstr;
}

Vstr* construct_vstr_copy(Vstr* the_vstr){
  Vstr* the_copy = (Vstr*)malloc(1*sizeof(Vstr));
  the_copy->capacity = the_vstr->capacity;
  the_copy->size = the_vstr->size;
  the_copy->a = (char**)malloc(the_copy->capacity * sizeof(char*));
  for(long i=0; i<the_vstr->size; i++){
    char* str = the_vstr->a[i];
    long len = strlen(str);
    char* str_copy = (char*)malloc((len+1)*sizeof(char));
    str_copy = strcpy(str_copy, str);
    the_copy->a[i] = str_copy;
  }
  return the_copy;
}

void add_string_to_vstr(Vstr* the_vstr, char* str){ // 
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

char* ith_str_from_vstr(Vstr* the_vstr, long i){ // returns a pointer to the ith string
  char* s = (i >= 0)?
    the_vstr->a[i] :
    the_vstr->a[the_vstr->size + i];
  return s;
}

char* copy_ith_str_from_vstr(Vstr* the_vstr, long i){ // returns a copy of the ith string 
  char* s = (i >= 0)?
    the_vstr->a[i] :
    the_vstr->a[the_vstr->size + i];
  return strcpy((char*)malloc((strlen(s)+1)*sizeof(char)), s);
}

void print_vstr(FILE* fh, Vstr* the_vstr){
  for(long i=0; i<the_vstr->size; i++){
    fprintf(fh, "%s ", the_vstr->a[i]);
  }
  fprintf(fh, "\n");
}

void free_vstr(Vstr* the_vstr){
  for(long i=0; i<the_vstr->size; i++){ 
    free(the_vstr->a[i]);
  }
  free(the_vstr->a);
  free(the_vstr);
}

// *****  IndexId *****
IndexId* construct_indexid(long idx, char* id){
  IndexId* the_idxid = (IndexId*)malloc(sizeof(IndexId));
  the_idxid->index = idx;
  the_idxid->id = strcpy((char*)malloc((strlen(id)+1) * sizeof(char)), id); 
  return the_idxid;
}

void free_indexid(IndexId* the_idxid){
  free(the_idxid->id);
  free(the_idxid);
}

// *****  Vidxid  *****
/* Vidxid* construct_vidxid_from_array(long size, char** strs){ */
/*   Vidxid* the_vidxid = (Vidxid*)malloc(sizeof(Vidxid)); */
/*   the_vidxid->capacity = size; */
/*   the_vidxid->size = size; */
/*   the_vidxid->a = (IndexId**)malloc(size*sizeof(IndexId*)); */
/*   for(long i=0; i<size; i++){ */
/*     IndexId* the_idxid = (IndexId*)malloc(sizeof(IndexId)); */
/*     the_idxid->index = i; */
/*     the_idxid->id = strcpy((char*)malloc((strlen(strs[i])+1)*sizeof(char)), strs[i]); */
/*     the_vidxid->a[i] = the_idxid; */
/*   } */
/*   return the_vidxid;   */
/* } */

Vidxid* construct_vidxid_from_vstr(Vstr* ids){
  Vidxid* the_vidxid = (Vidxid*)malloc(sizeof(Vidxid));
  the_vidxid->capacity = ids->size;
  the_vidxid->size = ids->size;
  the_vidxid->a = (IndexId**)malloc(ids->size*sizeof(IndexId*));
  for(long i=0; i<ids->size; i++){
    IndexId* the_idxid = (IndexId*)malloc(sizeof(IndexId));
    the_idxid->index = i;
    the_idxid->id = strcpy((char*)malloc((strlen(ids->a[i])+1)*sizeof(char)), ids->a[i]);
    the_vidxid->a[i] = the_idxid;
  }
  return the_vidxid;  
}

Vidxid* construct_sorted_vidxid_from_vstr(Vstr* ids){
  Vidxid* the_vidxid = construct_vidxid_from_vstr(ids);
  sort_vidxid_by_id(the_vidxid);
   if(DBUG) assert(check_idxid_map(the_vidxid, ids) == 1);
  return the_vidxid;
}

int strcmpx(const void* v1, const void* v2){
  const IndexId** s1 = (const IndexId**)v1;
  const IndexId** s2 = (const IndexId**)v2;
  return strcmp((*s1)->id, (*s2)->id);
}

void sort_vidxid_by_id(Vidxid* the_vidxid){ // sort in place
  qsort(the_vidxid->a, the_vidxid->size, sizeof(IndexId*), strcmpx);
}

long index_of_id_in_vidxid(Vidxid* the_vidxid, char* id){
  // find the index corresponding to id
  // by binary search
  long lb = 0;
  long ub = the_vidxid->size - 1;
  long idx_guess = (lb + ub)/2;
  char* id_guess = the_vidxid->a[idx_guess]->id;
  // fprintf(stderr, "A: lb: %ld  idx_guess: %ld  ub: %ld  %s\n", lb, idx_guess, ub, id);
   if(lb > ub) getchar();
  int icmp;
  while((icmp = strcmp(id_guess, id)) != 0){
    if(ub - lb > 1){
    if(icmp < 0){
      lb = idx_guess;
    }else if(icmp > 0){
      ub = idx_guess;
    }else{
      exit(EXIT_FAILURE);
    }
    idx_guess = (lb + ub)/2;
    id_guess = the_vidxid->a[idx_guess]->id;
    //   fprintf(stderr, "B: lb: %ld  idx_guess: %ld  ub: %ld   [%s]  [%s]\n", lb, idx_guess, ub, id, id_guess);
    }else if(ub - lb == 1){
      idx_guess = lb;
      id_guess = the_vidxid->a[idx_guess]->id;
      if(strcmp(id_guess, id) == 0) {
	//	fprintf(stderr, "C: lb: %ld  idx_guess: %ld  ub: %ld   [%s]  [%s]\n", lb, idx_guess, ub, id, id_guess);
	//	fprintf(stderr, "in index_of_... index found: %ld\n", the_vidxid->a[idx_guess]->index);
	return the_vidxid->a[idx_guess]->index;
      }

      idx_guess = ub;
      id_guess = the_vidxid->a[idx_guess]->id;
      if(strcmp(the_vidxid->a[ub]->id, id) == 0){
	//	fprintf(stderr, "D: lb: %ld  idx_guess: %ld  ub: %ld   [%s]  [%s]\n", lb, idx_guess, ub, id, id_guess);
	//	fprintf(stderr, "in index_of_... index found: %ld\n", the_vidxid->a[idx_guess]->index);
	return the_vidxid->a[idx_guess]->index;
      }
      //          fprintf(stderr, "E: lb: %ld  idx_guess: %ld  ub: %ld   %s\n", lb, idx_guess, ub, id);
      return -1;
    }else{
      //    fprintf(stderr, "F: lb: %ld  idx_guess: %ld  ub: %ld   %s\n", lb, idx_guess, ub, id);
      exit(EXIT_FAILURE);
    }
  }
  long result = the_vidxid->a[idx_guess]->index;
  //  fprintf(stderr, "in index_of_ ... index found: %ld \n", result);
  return result;
}

long check_idxid_map(Vidxid* vidxid, Vstr* accession_ids){
  for(long i=0; i<accession_ids->size; i++){
    char* id = accession_ids->a[i];
    long idx = index_of_id_in_vidxid(vidxid, id);
    if(idx != i) return 0;
  }
  return 1;
}

void print_vidxid(FILE* fh, Vidxid* the_vidxid){
  for(long i=0; i<the_vidxid->size; i++){
    IndexId* the_idxid = the_vidxid->a[i];
    fprintf(fh, "%ld  %s\n", the_idxid->index, the_idxid->id);
  }
}

void free_vidxid(Vidxid* the_vidxid){
  for(long i=0; i<the_vidxid->size; i++){
    free_indexid(the_vidxid->a[i]);
  }
  free(the_vidxid->a);
  free(the_vidxid);
}

