// *********  typedefs  ********

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

// *****  Function declarations  ************************************************************
// ***** Vlong ******************************************************************************
Vlong* construct_vlong(long cap); // set capacity = cap, size = 0
Vlong* construct_vlong_zeroes(long size);
Vlong* construct_vlong_from_array(long size, long* array); // initialize with array of known size
Vlong* construct_vlong_whole_numbers(long size); // initialize to 0,1,2,3,...size-1
void add_long_to_vlong(Vlong* the_vlong, long x); // push, realloc if necessary
void shuffle_vlong(Vlong* the_vlong); // randomize order of array elements
void free_vlong(Vlong* the_vlong); // free memory


// *****  Vstr  *****************************************************************************
Vstr* construct_vstr(long cap); // set capacity = cap, size = 0
Vstr* construct_vstr_copy(Vstr* the_vstr);
void add_string_to_vstr(Vstr* the_vstr, char* str); // push, realloc if necessary
void free_vstr(Vstr* the_vstr); // free memory

// *****  Vchar  (string as vector of chars)  ***********************************************
// Vstr* construct_vstr_from_string(char* str); // initialize with null-terminated string
