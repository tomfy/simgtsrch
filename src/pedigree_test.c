// C version of program to test pedigrees using genotype data.
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>
// #include "vect.h"
#include "gtset.h"

#define PEDIGREE_FIELDS 7

GenotypesSet* read_genotypes_file_and_store(FILE* g_stream, double delta, double max_missing_data_fraction);
Vpedigree* read_the_pedigrees_file_and_store(FILE* p_stream, Vidxid* the_vidxid); 
void print_pedigree_stats(Pedigree* the_pedigree, GenotypesSet* the_gtsset);
// **********************************************************************************************
// ***********************************  main  ***************************************************
// **********************************************************************************************

int
main(int argc, char *argv[])
{
  double delta = 0.05; // default; control this with -d command line option.
  double max_marker_missing_data = 0.2; // default; control this with -x command line option.
  long max_number_of_accessions = 1000000;

  
  // ***** process command line *****
  if (argc < 2) {
    fprintf(stderr, "Usage:  %s -g <genotypes_file>  -p <pedigree_file>  options -d -x \n", argv[0]);
    exit(EXIT_FAILURE);
  }

  char* genotypes_filename = NULL;
  FILE *g_stream = NULL;
  char* pedigrees_filename = NULL;
  FILE *p_stream = NULL;

  // g: genotypes filename, p pedigree filename,  d (delta for rounding), x max fraction of missing data for markers.
  int c;
  while((c = getopt(argc, argv, "g:p:d:x:")) != -1){
    // fprintf(stderr, "c: %c %s %d\n", c, optarg, optind);
    switch(c){
    case 'g':
      genotypes_filename = optarg;
      g_stream = fopen(genotypes_filename, "r");
      if(g_stream == NULL){
	fprintf(stderr, "Failed to open %s for reading.\n", genotypes_filename);
	exit(EXIT_FAILURE);
      }
      break;
    case 'p':
      pedigrees_filename = optarg;
      p_stream = fopen(pedigrees_filename, "r");
      if(p_stream == NULL){
	fprintf(stderr, "Failed to open %s for reading.\n", pedigrees_filename);
	exit(EXIT_FAILURE);
      }
      break;
    case 'd':
      if(optarg == 0){
	perror("option d requires a numerical argument > 0\n");
	exit(EXIT_FAILURE);
      }else{
	delta = atof(optarg);
      }
      break;
    case 'x':
      if(optarg == 0){
	perror("option x requires a numerical argument > 0\n");
	exit(EXIT_FAILURE);
      }else{
	max_marker_missing_data = atof(optarg);
      }
      break;
    case '?':
      printf("? case in command line processing switch.\n");
      if ((optopt == 'i') || (optopt == 'n') || (optopt == 'k'))
	fprintf(stderr, "Option -%c requires an argument.\n", optopt);
      else if (isprint (optopt))
	fprintf(stderr, "Unknown option `-%c'.\n", optopt);
      else
	fprintf(stderr, "Unknown option character: %d\n", optopt);
      exit(EXIT_FAILURE);
    default:
      perror("default case (abort)\n");
      abort ();
    } // end of switch block
  } // end of loop over c.l. arguments
  // printf("optind: %d argc: %d\n", optind, argc);
  if(optind < argc){
    perror("Non-option arguments. Bye.\n");
    exit(EXIT_FAILURE);
  }
  if(genotypes_filename == NULL){
    perror("must specify genotype filename: -i <filename>");
    exit(EXIT_FAILURE);
  }
  if(pedigrees_filename == NULL){
    perror("must specify pedigrees filename: -i <filename>");
    exit(EXIT_FAILURE);
  }

  fprintf(stderr, "# genotypes file: %s  pedigree file: %s  delta: %5.3lf  max marker missing data: %5.3lf \n",
	  genotypes_filename, pedigrees_filename, delta, max_marker_missing_data);

  // *****  done processing command line  *****

  // ***************  read the genotypes file  *******************************
  GenotypesSet* the_genotypes_set = read_genotypes_file_and_store(g_stream, delta, max_marker_missing_data);
  fclose(g_stream);
  fprintf(stderr, "after read_genotypes_ ...\n");
  // print_genotypesset(the_genotypes_set);
  

  long n_accessions = the_genotypes_set->n_accessions;
  long n_markers_all = the_genotypes_set->n_markers;

  Vidxid* the_vidxid = construct_vidxid_from_vstr(the_genotypes_set->accession_ids);
  sort_vidxid_by_id(the_vidxid);
  //  print_vidxid(the_vidxid);

  for(long i=0; i<n_accessions; i++){
    // fprintf(stderr, "i:  %ld \n", i);
    char* id = the_genotypes_set->accession_ids->a[i];
    // fprintf(stderr, "%ld  %s\n", i, id);
    long idx = index_of_id_in_vidxid(the_vidxid, id);
    //  char* id2 = 
    if(idx != i) fprintf(stderr, "%ld  %s  %ld  %s\n", i, id, idx, the_genotypes_set->accession_ids->a[idx]);
  }

  // ***************  read the pedigrees file  ***************************
 
  Vpedigree* pedigrees = read_the_pedigrees_file_and_store(p_stream, the_vidxid);
  fclose(p_stream);

  // ***************  Done reading input files  ******************************

  //print_genotypesset(the_genotypes_set);
  print_genotypesset_summary_info(the_genotypes_set);
  if(DBUG) check_genotypesset(the_genotypes_set, max_marker_missing_data);

  // *****  clean genotypes set, i.e. remove markers with high missing data  ****
  GenotypesSet* the_cleaned_genotypes_set = construct_cleaned_genotypesset(the_genotypes_set, max_marker_missing_data);
  long n_markers_good = the_cleaned_genotypes_set->n_markers;
  free_genotypesset(the_genotypes_set);
  print_genotypesset_summary_info(the_cleaned_genotypes_set);
  // print_genotypesset(the_cleaned_genotypes_set);
  if(DBUG) check_genotypesset(the_cleaned_genotypes_set, max_marker_missing_data);

  //
  fprintf(stderr, "sizeof(Pedigree*): %ld \n", (long) sizeof(Pedigree*));

  for(long i=0; i<pedigrees->size; i++){
    print_pedigree_stats(pedigrees->a[i], the_cleaned_genotypes_set);
  }

  // ********************  cleanup  **************************q
 
  free_genotypesset(the_cleaned_genotypes_set);
  free_vidxid(the_vidxid);
  free_vpedigree(pedigrees);
  // getchar();
}
// **********************************************************
// ********************  end of main  ***********************
// **********************************************************




// **********************************************************
// ******************  functions  ***************************
// **********************************************************

GenotypesSet* read_genotypes_file_and_store(FILE* g_stream, double delta, double max_missing_data_fraction){
  char* line = NULL;
  size_t len = 0;
  ssize_t nread;

  long markerid_count = 0;
  Vstr* marker_ids = construct_vstr(1000);
  char* saveptr = NULL;
  if((nread = getline(&line, &len, g_stream)) != -1){
    char* token = strtok_r(line, "\t \n\r", &saveptr);
    // printf("0line (pointer): %p \n", line);
    if((token == NULL)  || (strcmp(token, "MARKER") != 0)){
      exit(EXIT_FAILURE);
    }
    // printf("%s  ", token);
    while(1){
      token = strtok_r(NULL, "\t \n\r", &saveptr);
      //  printf("line (pointer): %p \n", line);
      if(token == NULL) break;
      char* mrkr_id = (char*)malloc((strlen(token)+1)*sizeof(char));
      add_string_to_vstr(marker_ids, strcpy(mrkr_id, token)); // store
      //  printf("%s ", token);
      markerid_count++;
    }
    // printf("\n");
  }
  //  free(line); line = NULL;
  fprintf(stderr, "# number of marker ids counted: %ld \n", markerid_count);

  long accession_count = 0;
  Vstr* accession_ids = construct_vstr(1000);
  Vstr* genotypes_strings = construct_vstr(1000);
  long* missing_data_counts = (long*)calloc(markerid_count, sizeof(long));
  Vlong* the_md_vlong = construct_vlong_from_array(markerid_count, missing_data_counts);
  while((nread = getline(&line, &len, g_stream)) != -1){
    accession_count++;
    char* token = strtok_r(line, "\t \n\r", &saveptr);
    //printf("accession number: %ld id: %s\n", accession_count, token);
    char* acc_id = (char*)malloc((strlen(token)+1)*sizeof(char));
    add_string_to_vstr(accession_ids, strcpy(acc_id, token)); // store copy of accession id
    //  printf("%s  ", token);
    long marker_count = 0;
    char* genotypes = (char*)malloc((markerid_count+1) * sizeof(char));
    genotypes[markerid_count] = '\0'; // terminate with null.
    while(1){
      token = strtok_r(NULL, "\t \n\r", &saveptr);
      if(token == NULL) break;
      double dosage = atof(token);

      if(dosage <= delta){
	genotypes[marker_count] = '0';
      }else if((dosage >= 1-delta) && (dosage <= 1+delta)){
	genotypes[marker_count] = '1';
      }else if(dosage >= 2-delta){
	genotypes[marker_count] = '2';
      }else{ // missing data
	genotypes[marker_count] = '3';
	missing_data_counts[marker_count]++;
      }
      marker_count++;
    } // done reading dosages for all markers, and rounding to 0, 1, 2, 3
    add_string_to_vstr(genotypes_strings, genotypes);
    if(marker_count != markerid_count) exit(EXIT_FAILURE);
    // printf("%s\n", genotypes);
  } // done reading all lines
  
  free(line); // only needs to be freed once.
  GenotypesSet* the_genotypes_set = construct_genotypesset(accession_ids, marker_ids, genotypes_strings, the_md_vlong);
  return the_genotypes_set;
}

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
    //  print_vstr(fields); printf("\n");
    
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

void print_pedigree_stats(Pedigree* the_pedigree, GenotypesSet* the_gtsset){

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


