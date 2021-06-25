#include <stdlib.h>
#include <stdio.h>
#include "gtset.h"
// #include "vect.h"
#include "pedigree.h"

#define PEDIGREE_FIELDS 7 // number of whitespace separated fields in pedigree file, with ids in last 3 fields.

// extern int do_checks_flag; // option -c sets this to 1 to do some checks.


// *****  Ahr  *****
void print_Ahr(FILE* fh, Ahr the_ahr){
  fprintf(fh, "%ld %8.5lf  %ld %8.5lf  %ld %8.5lf",
	  the_ahr.a.d, (the_ahr.a.d > 0)? (double)the_ahr.a.n/(double)the_ahr.a.d : -1,
	  the_ahr.h.d, (the_ahr.h.d > 0)? (double)the_ahr.h.n/(double)the_ahr.h.d : -1,
	  the_ahr.r.d, (the_ahr.r.d > 0)? (double)the_ahr.r.n/(double)the_ahr.r.d : -1
	  );
}

// *****  Pedigree  *****
Pedigree* construct_pedigree(IndexId* acc_idxid, IndexId* fempar_idxid, IndexId* malpar_idxid){
  Pedigree* the_pedigree = (Pedigree*)malloc(sizeof(Pedigree));
  the_pedigree->Fparent = fempar_idxid;
  the_pedigree->Mparent = malpar_idxid;
  the_pedigree->Accession = acc_idxid;
  return the_pedigree;
}

void agmr_hgmr_r(char* gts1, char* gts2, Ahr* the_ahr){
  char c1, c2;
  long n_00 = 0;
  long n_01 = 0;
  long n_02 = 0;
  long n_10 = 0;
  long n_11 = 0;
  long n_12 = 0;
  long n_20 = 0;
  long n_21 = 0;
  long n_22 = 0;
  long n_denom = 0;
  long i=0;
  while((c1 = gts1[i]) != '\0'){
    if(c1 != '3'){ // not missing data
      c2 = gts2[i];
      if(c2 != '3'){ // not missing data
	if(c1 == '0'){
	  if(c2 == '0'){
	    n_00++;
	  }else if(c2 == '1'){
	    n_01++;
	  }else if(c2 == '2'){
	    n_02++;
	  }
	}else if(c1 == '1'){
	  if(c2 == '0'){
	    n_10++;
	  }else if(c2 == '1'){
	    n_11++;
	  }else if(c2 == '2'){
	    n_12++;
	  }
	}else if(c1 == '2'){
	  if(c2 == '0'){
	    n_20++;
	  }else if(c2 == '1'){
	    n_21++;
	  }else if(c2 == '2'){
	    n_22++;
	  }
	}
      }
    }
    i++;
  }
  //  printf("%ld %ld %ld %ld\n", n_00, n_11, n_22, n_10);
  long a_numer = n_01 + n_02 + n_10 + n_12 + n_20 + n_21;
  long a_denom = a_numer + n_00 + n_11 + n_22;
  double agmr = (a_denom > 0)? (double)a_numer/(double)a_denom : -1; // agmr small: genotype sets are similar

  long h_numer = n_02 + n_20;
  long h_denom = h_numer + n_00 + n_22;
  double hgmr = (h_denom > 0)? (double)h_numer/(double)h_denom : -1; // hgmr small: gts1 likely to be parent of gts2 or vice versa

  long r_numer = n_01 + n_21;
  long r_denom = r_numer + n_00 + n_22;
  double r = (r_denom)? (double)r_numer/(double)r_denom : -1; // r small: if gts1 is parent of gts2, other parent likely also be gts1

  the_ahr->a.n = a_numer;
  the_ahr->a.d = a_denom;
  the_ahr->h.n = h_numer;
  the_ahr->h.d = h_denom;
  the_ahr->r.n = r_numer;
  the_ahr->r.d = r_denom;
}


double hgmr(char* gts1, char* gts2){
  char c1, c2;
  long n_numer = 0;
  long n_denom = 0;
  long i=0;
  while((c1 = gts1[i]) != '\0'){
    if((c1 == '0') || (c1 == '2')){
      c2 = gts2[i];
      if((c2 == '0') || (c2 == '2')){
	n_denom++;
	if(c1 != c2) n_numer++;
      }
    }
    i++;
  }
  return (n_denom > 0)? (double)n_numer/(double)n_denom : -1;  
}

void calculate_pedigree_test_info(Pedigree* the_pedigree, GenotypesSet* the_gtsset){

 char* fempar_gts = the_gtsset->genotype_sets->a[the_pedigree->Fparent->index];
  char* malpar_gts = the_gtsset->genotype_sets->a[the_pedigree->Mparent->index];
  char* acc_gts = the_gtsset->genotype_sets->a[the_pedigree->Accession->index];
    
  Ahr fp_ahr, mp_ahr, fm_ahr;
  agmr_hgmr_r(fempar_gts, acc_gts, &(the_pedigree->fp_ahr));
  agmr_hgmr_r(malpar_gts, acc_gts, &(the_pedigree->mp_ahr));
  agmr_hgmr_r(fempar_gts, malpar_gts, &(the_pedigree->fm_ahr));
}

void print_pedigree_test_info(FILE* fh, Pedigree* the_pedigree, GenotypesSet* the_gtsset){
  fprintf(fh, "%s %ld %s %s  ",
	 the_gtsset->accession_ids->a[the_pedigree->Accession->index],
	 the_gtsset->accession_missing_data_counts->a[the_pedigree->Accession->index],
	 the_gtsset->accession_ids->a[the_pedigree->Fparent->index],
	 the_gtsset->accession_ids->a[the_pedigree->Mparent->index]
	 );
	  
  print_Ahr(fh, the_pedigree->fp_ahr); fprintf(fh, "  ");
  print_Ahr(fh, the_pedigree->mp_ahr); fprintf(fh, "  ");
  print_Ahr(fh, the_pedigree->fm_ahr); fprintf(fh, "\n");
}

void free_pedigree(Pedigree* the_pedigree){
  free_indexid(the_pedigree->Fparent);
  free_indexid(the_pedigree->Mparent);
  free_indexid(the_pedigree->Accession);
  free(the_pedigree);
}

// *****  Vpedigree  *****

Vpedigree* read_the_pedigrees_file_and_store(FILE* p_stream, Vidxid* the_vidxid){
  Vpedigree* pedigrees = construct_vpedigree(1000);

  char* line = NULL;
  size_t len = 0;
  ssize_t nread;

  char* saveptr = NULL;
  if((nread = getline(&line, &len, p_stream)) != -1){
    char* token = strtok_r(line, "\t \n\r", &saveptr);
    // printf("0line (pointer): %p \n", line);
    if((token == NULL)  || (strcmp(token, "Accession") != 0)){
      exit(EXIT_FAILURE);
    }
  }
  // fprintf(stderr, "after reading first (header) line of pedigrees file\n");
  long i_pedigree = 0;
  while((nread = getline(&line, &len, p_stream)) != -1){
    Vstr* fields = construct_vstr(PEDIGREE_FIELDS);
    char* token = strtok_r(line, "\t \n\r", &saveptr);
    // fprintf(stderr, "token: %s\n", token);
    add_string_to_vstr(fields, strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token)); // store copy of accession id
    while(1){
      token = strtok_r(NULL, "\t \n\r", &saveptr);
      if(token == NULL) break;
      add_string_to_vstr(fields, strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token)); // store copy of accession id
    }
    // construct a Pedigree struct from last 3 fields  
    char* acc_id = ith_str_from_vstr(fields, -3); // ith_str ... copies the string, i.e. allocates more memory
    char* fempar_id = ith_str_from_vstr(fields, -2);
    char* malpar_id = ith_str_from_vstr(fields, -1);
    // fprintf(stderr, "[%s] [%s] [%s]\n", acc_id, fempar_id, malpar_id);
   
  
    long acc_idx, fempar_idx, malpar_idx;
    //  fprintf(stderr, "acc_idx: %ld\n", index_of_id_in_vidxid(the_vidxid, acc_id));
    if(
       (strcmp(acc_id, "NA") != 0) &&
       (strcmp(fempar_id, "NA") != 0) &&
       (strcmp(malpar_id, "NA") != 0) &&
       ((acc_idx = index_of_id_in_vidxid(the_vidxid, acc_id)) != -1) &&
       ((fempar_idx = index_of_id_in_vidxid(the_vidxid, fempar_id)) != -1) &&
       ((malpar_idx = index_of_id_in_vidxid(the_vidxid, malpar_id)) != -1)
       ){
      IndexId* acc_idxid = construct_indexid(acc_idx, acc_id);
      IndexId* fempar_idxid = construct_indexid(fempar_idx, fempar_id);
      IndexId* malpar_idxid = construct_indexid(malpar_idx, malpar_id);
      //  fprintf(stderr, "pedigree has 3 valid ids: %s %s %s\n", acc_id, fempar_id, malpar_id);
      add_pedigree_to_vpedigree(pedigrees, construct_pedigree(acc_idxid, fempar_idxid, malpar_idxid));
    }
    free_vstr(fields);
  } // done reading all lines
  free(line); // only needs to be freed once.
  fprintf(stderr, "size of Vpedigree pedigrees: %ld \n", pedigrees->size);
  return pedigrees;
}

Vpedigree* construct_vpedigree(long cap){
  Vpedigree* the_vped = (Vpedigree*)malloc(sizeof(Vpedigree));
  the_vped->capacity = cap;
  the_vped->size = 0;
  the_vped->a = (Pedigree**)malloc(the_vped->capacity*sizeof(Pedigree*));
}
  
void add_pedigree_to_vpedigree(Vpedigree* the_vped, Pedigree* the_ped){
  long cap = the_vped->capacity;
  long n = the_vped->size;
  if(n == cap){
    cap *= 2;
    the_vped->a = (Pedigree**)realloc(the_vped->a, cap*sizeof(Pedigree*));
    the_vped->capacity = cap;
    //  fprintf(stderr, "in add_pedigree_to... capacity increased to %ld. pointer to new memory block: %p\n", cap, the_vped->a);
  }
  //  printf("after realloc. cap: %ld \n", cap);
  the_vped->a[n] = the_ped;
  //  printf("after assignment to a[n]\n");
  the_vped->size++;
  // printf("about to return from add_pedigree_to_vpedigree. size %ld\n", the_vped->size);
}

void free_vpedigree(Vpedigree* the_vped){
  //fprintf(stderr, "in free_vpedigree. size: %ld\n", the_vped->size);
  for(long i=0; i<the_vped->size; i++){
    //fprintf(stderr, "freeing pedigree %ld \n", i);
    free_pedigree(the_vped->a[i]);
  }
  free(the_vped->a);
  free(the_vped);
}
