// *****  Implementation of functions declared in vect.h  *****
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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

void free_vstr(Vstr* the_vstr){
  for(long i=0; i<the_vstr->size; i++){ 
    free(the_vstr->a[i]);
  }
  free(the_vstr->a);
  free(the_vstr);
}

