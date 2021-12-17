#include <math.h>
#include "gtset.h"

extern int do_checks_flag; // option -c sets this to 1 to do some checks.

// *****  Accession  implementation *****
Accession* construct_accession(char* id, long idx, char* genotypes, long accession_md_count){
  Accession* the_accession = (Accession*)malloc(1*sizeof(Accession));
  the_accession->id = construct_vchar_from_str(id);
  the_accession->index = idx;
  the_accession->genotypes = construct_vchar_from_str(genotypes);
  the_accession->missing_data_count = accession_md_count;
  return the_accession;
}
void set_accession_missing_data_count(Accession* the_accession, long missing_data_count){
  the_accession->missing_data_count = missing_data_count;
}
void free_accession(Accession* the_accession){
  free_vchar(the_accession->id);
  free_vchar(the_accession->genotypes);
  free(the_accession);
}
void free_accession_innards(Accession* the_accession){ // doesn't free the_accession itself
  // use to free each element of an array of Accession (not an array of Accession*)
  free_vchar(the_accession->id);
  free_vchar(the_accession->genotypes);
}

// *****  Vaccession implementation  *****
Vaccession* construct_vaccession(long cap){ // construct Vaccession with capacity of cap, but empty.
  Vaccession* the_vacc = (Vaccession*)malloc(sizeof(Vaccession));
  the_vacc->capacity = cap;
  the_vacc->size = 0;
  the_vacc->a = (Accession**)malloc(cap*sizeof(Accession*));
}

void add_accession_to_vaccession(Vaccession* the_vacc, Accession* the_acc){
  long cap = the_vacc->capacity;
  long n = the_vacc->size;
  if(n == cap){   // if necessary, resize w realloc
    cap *= 2;
    the_vacc->a = (Accession**)realloc(the_vacc->a, cap*sizeof(Accession*));
    the_vacc->capacity = cap;
  }
  the_vacc->a[n] = the_acc;
  the_vacc->size++;
}

void free_vaccession(Vaccession* the_vacc){
  for(long i=0; i<the_vacc->size; i++){
    free_accession(the_vacc->a[i]);
  }
  free(the_vacc->a);
  free(the_vacc);
}

// *****  GenotypesSet implementation *****
GenotypesSet* read_dosages_file_and_store(FILE* g_stream, double delta){
  char* line = NULL;
  size_t len = 0;
  ssize_t nread;

  long accessions_capacity = 2000;
  Vaccession* the_accessions = construct_vaccession(accessions_capacity); // (Vaccession*)malloc(accessions_capacity*sizeof(Vaccession)); 

  long markerid_count = 0;
  Vstr* marker_ids = construct_vstr(1000);
  char* saveptr = NULL;
  if((nread = getline(&line, &len, g_stream)) != -1){
    char* token = strtok_r(line, "\t \n\r", &saveptr);
    if((token == NULL)  || (strcmp(token, "MARKER") != 0)){
      fprintf(stderr, "token: %s (should be MARKER)\n", token);
      exit(EXIT_FAILURE);
    }
    while(1){
      token = strtok_r(NULL, "\t \n\r", &saveptr);
      if(token == NULL) break;
      char* mrkr_id = (char*)malloc((strlen(token)+1)*sizeof(char));
      add_string_to_vstr(marker_ids, strcpy(mrkr_id, token)); // store
      markerid_count++;
    }
  }

  long accession_count = 0;
  // long* missing_data_counts = (long*)calloc(markerid_count, sizeof(long));
  Vlong* the_md_vlong = construct_vlong_zeroes(markerid_count); // construct_vlong_from_array(markerid_count, missing_data_counts);
  while((nread = getline(&line, &len, g_stream)) != -1){
  
    char* token = strtok_r(line, "\t \n\r", &saveptr);
    char* acc_id = strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token);
    long marker_count = 0;
    long accession_missing_data_count = 0;
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
	the_md_vlong->a[marker_count]++;
	accession_missing_data_count++;
      }
      marker_count++;
    } // done reading dosages for all markers, and rounding to 0, 1, 2, 3
    if(marker_count != markerid_count) exit(EXIT_FAILURE);
    /* if(accession_count == accessions_capacity){ // accession_count is the number of accessions already stored in the_accessions */
    /*   accessions_capacity *= 2; */
    /*   //    fprintf(stderr, "increasing capacity of accessions array to %ld \n", accessions_capacity); */
    /*   the_accessions = (Accession*)realloc(the_accessions, accessions_capacity*sizeof(Accession)); */
    /* } */
    /* the_accessions[accession_count].id = construct_vchar_from_str(acc_id);  free(acc_id); // or cut out the middleman (acc_id)? */
    /* the_accessions[accession_count].genotypes = construct_vchar_from_str(genotypes); free(genotypes); */
    /* the_accessions[accession_count].missing_data_count = accession_missing_data_count; */

    Accession* the_accession = construct_accession(acc_id, accession_count, genotypes, accession_missing_data_count);
    free(acc_id); // or cut out the middleman (acc_id)?
    free(genotypes);
    //  the_accession->missing_data_count = accession_missing_data_count;
    add_accession_to_vaccession(the_accessions, the_accession);
  
    accession_count++;
  } // done reading all lines

  free(line); // only needs to be freed once.
  GenotypesSet* the_genotypes_set = //construct_genotypesset(accession_ids, marker_ids, genotypes_strings, the_md_vlong);
    construct_genotypesset(the_accessions, marker_ids, the_md_vlong, delta, 1);
 
  return the_genotypes_set;
}

/* */
GenotypesSet* read_genotypes_file_and_store(FILE* g_stream){ // 
  char* line = NULL;
  size_t len = 0;
  ssize_t nread;

  long accessions_capacity = 2000;
  Vaccession* the_accessions = construct_vaccession(accessions_capacity); //(Vaccession*)malloc(accessions_capacity*sizeof(ccession));

  double delta;
  double max_marker_missing_data_fraction;
  char* saveptr = NULL;
  if((nread = getline(&line, &len, g_stream)) != -1){
    char* token = strtok_r(line, "\t \n\r", &saveptr);
    if((token == NULL)  || (strcmp(token, "#") != 0)){
      exit(EXIT_FAILURE);
    }
    long field_count = 1;
    while(1){
      token = strtok_r(NULL, "\t \n\r", &saveptr);
      if(token == NULL) break;
      if(field_count == 2){
	delta = atof(token);
      }else if(field_count == 4){
	max_marker_missing_data_fraction = atof(token);
      }
    }
  }

  long markerid_count = 0;
  Vstr* marker_ids = construct_vstr(1000);
  saveptr = NULL;
  if((nread = getline(&line, &len, g_stream)) != -1){
    char* token = strtok_r(line, "\t \n\r", &saveptr);
    if((token == NULL)  || (strcmp(token, "MARKER") != 0)){
      exit(EXIT_FAILURE);
    }
    while(1){
      token = strtok_r(NULL, "\t \n\r", &saveptr);
      if(token == NULL) break;
      char* mrkr_id = (char*)malloc((strlen(token)+1)*sizeof(char));
      add_string_to_vstr(marker_ids, strcpy(mrkr_id, token)); // store
      markerid_count++;
    }
  }

  long accession_count = 0;
  long* marker_missing_data_counts = (long*)calloc(markerid_count, sizeof(long));
  Vlong* the_md_vlong = construct_vlong_zeroes(markerid_count); // construct_vlong_from_array(markerid_count, marker_missing_data_counts);
  saveptr = NULL;
  while((nread = getline(&line, &len, g_stream)) != -1){
  
    char* token = strtok_r(line, "\t \n\r", &saveptr);
    char* acc_id = strcpy((char*)malloc((strlen(token)+1)*sizeof(char)), token);
    long token_count = 1;
    long accession_missing_data_count = 0;
    char* genotypes = (char*)malloc((markerid_count+1) * sizeof(char));
    genotypes[markerid_count] = '\0'; // terminate with null.
    while(1){
      token = strtok_r(NULL, "\t \n\r", &saveptr);
      if(token == NULL) break;
      token_count++;
      genotypes = strcpy(genotypes, token);
    }
    assert(token_count == 2); // because line should be: accession_id 00102011000201020....\n  just 2 whitespace-separated tokens.
    assert(strlen(genotypes) == markerid_count);
    for(long igt = 0; igt < markerid_count; igt++){ // count missing data     
      if(genotypes[igt] == '3'){
	marker_missing_data_counts[igt]++;
	accession_missing_data_count++;
      }
    } // done reading genotypes for all markers
  
    Accession* the_accession = construct_accession( acc_id, accession_count, genotypes, accession_missing_data_count);
    free(acc_id); // or cut out the middleman (acc_id)?
    free(genotypes);
    //   the_accession->missing_data_count = accession_missing_data_count;
    add_accession_to_vaccession(the_accessions, the_accession);
    accession_count++;
  } // done reading all lines
  
  free(line); // only needs to be freed once.
  GenotypesSet* the_genotypes_set = //construct_genotypesset(accession_ids, marker_ids, genotypes_strings, the_md_vlong);
    construct_genotypesset(the_accessions, marker_ids, the_md_vlong, delta, max_marker_missing_data_fraction);
  return the_genotypes_set;
}
/* */


GenotypesSet* construct_genotypesset(Vaccession* accessions, Vstr* marker_ids, Vlong* md_counts, double delta, double max_marker_md_fraction){
 
  GenotypesSet* the_gtsset = (GenotypesSet*)malloc(1*sizeof(GenotypesSet));
  the_gtsset->delta = delta;
  the_gtsset->max_marker_missing_data_fraction = max_marker_md_fraction;
  the_gtsset->n_accessions = accessions->size;
  the_gtsset->n_markers = marker_ids->size;
  the_gtsset->accessions = accessions;
  the_gtsset->marker_ids = marker_ids;
  the_gtsset->marker_missing_data_counts = md_counts;
  return the_gtsset;
}

void check_genotypesset(GenotypesSet* gtss, double max_marker_md_fraction){
  assert(gtss->marker_ids->size == gtss->n_markers);
  long* md_counts = (long*)calloc(gtss->n_markers, sizeof(long)); 
  for(long i=0; i<gtss->n_accessions; i++){
    Accession* an_acc = gtss->accessions->a[i];
    assert(strlen(an_acc->genotypes->a) == gtss->n_markers);
    for(long j=0; j<gtss->n_markers; j++){
      if(an_acc->genotypes->a[j] == '3')md_counts[j]++;
    }
  }
  for(long j=0; j<gtss->n_markers; j++){
    assert(md_counts[j] == gtss->marker_missing_data_counts->a[j]);
  }
  free(md_counts);
  fprintf(stderr, "# Successfully completed check_genotypesset\n");
}

GenotypesSet* construct_cleaned_genotypesset(const GenotypesSet* the_gtsset, double max_marker_md_fraction){
  Vlong* md_counts = the_gtsset->marker_missing_data_counts;
  long n_accs = the_gtsset->n_accessions;
  long* mdcount_histogram = (long*)calloc(n_accs+1, sizeof(long));
  if(max_marker_md_fraction < 0){
    double factor = -1*max_marker_md_fraction;
    for(long i=0; i< md_counts->size; i++){
      mdcount_histogram[md_counts->a[i]]++;
    }
    long mrkrs_so_far = 0;
    long median_md_count;
    for(long i=0; i<=n_accs; i++){
      mrkrs_so_far += mdcount_histogram[i];
      if(mrkrs_so_far > 0.5*md_counts->size){
	median_md_count = i;
	break;
      }
    }
    max_marker_md_fraction = factor*(double)median_md_count/(double)n_accs;
  }
  // identify the markers to keep:
  long n_markers_to_keep = 0;
  Vlong* md_ok = construct_vlong_zeroes(md_counts->size);
  Vstr* cleaned_marker_ids = construct_vstr(1000);
  Vlong* cleaned_md_counts = construct_vlong(1000);
  long mdsum_all = 0;
  long mdsum_kept = 0;
  for(long i=0; i<md_counts->size; i++){
    mdsum_all += md_counts->a[i];
    if(md_counts->a[i] <= max_marker_md_fraction*the_gtsset->n_accessions){
      md_ok->a[i] = 1;
      n_markers_to_keep++;
      mdsum_kept += md_counts->a[i];
      long marker_id_length = strlen(the_gtsset->marker_ids->a[i]);
      char* marker_id_to_keep = strcpy((char*)malloc((marker_id_length+1)*sizeof(char)), the_gtsset->marker_ids->a[i]);
      add_string_to_vstr(cleaned_marker_ids, marker_id_to_keep);
      add_long_to_vlong(cleaned_md_counts, md_counts->a[i]);
    }
  }
  double raw_md_fraction = (double)mdsum_all/(double)(md_counts->size*n_accs);
  double cleaned_md_fraction = (double)mdsum_kept/(double)(n_markers_to_keep*n_accs);
  fprintf(stderr, "# Removing markers with fraction of missing data greater than: %7.3lf\n", max_marker_md_fraction);
  fprintf(stderr, "# Raw data has %ld markers, with missing data fraction of %7.3lf\n", md_counts->size, raw_md_fraction);
  fprintf(stderr, "# Cleaned data has %ld markers, with missing data fraction of: %7.3lf\n", n_markers_to_keep, cleaned_md_fraction);

  Vaccession* the_accessions = construct_vaccession(the_gtsset->n_accessions); //(Accession*)malloc(the_gtsset->n_accessions*sizeof(Accession)); 
  for(long i=0; i<the_gtsset->n_accessions; i++){ // loop over accessions
    char* raw_gts = the_gtsset->accessions->a[i]->genotypes->a; // the string with all the genotypes for accession i
    char* cleaned_gts = (char*)malloc((n_markers_to_keep+1)*sizeof(char));
    long k=0; // k: index of kept markers
    long acc_md_count = 0;
    for(long j=0; j<the_gtsset->n_markers; j++){ // j: index of original markers
      if(md_ok->a[j] == 1){
	cleaned_gts[k] = raw_gts[j];
	k++;
	if(raw_gts[j] == '3') acc_md_count++;
      }
    }
    cleaned_gts[k] = '\0'; // terminate with null.
    if(DBUG && do_checks_flag) assert(k == n_markers_to_keep);
    
    Accession* the_accession = construct_accession( the_gtsset->accessions->a[i]->id->a, i, cleaned_gts, acc_md_count ); //(Accession*)malloc(sizeof(Accession));
    free(cleaned_gts);
    //   the_accession->missing_data_count = acc_md_count;
    add_accession_to_vaccession(the_accessions, the_accession);
    //   fprintf(stderr, "xx: %ld %s\n", the_accession->index, the_accession->id->a);
  }
  free_vlong(md_ok);

  GenotypesSet* cleaned_gtsset = construct_genotypesset(
							the_accessions,
							cleaned_marker_ids, cleaned_md_counts,
							the_gtsset->delta, max_marker_md_fraction
							);
  fprintf(stderr, "# About to return cleaned_gtsset \n");
  return cleaned_gtsset;
}

void print_genotypesset(FILE* fh, GenotypesSet* the_gtsset){
  fprintf(fh, "# delta: %8.4lf  max_marker_missing_data_fraction: %8.4lf \n", the_gtsset->delta, the_gtsset->max_marker_missing_data_fraction);
  fprintf(fh, "MARKER  ");
  for(long i=0; i<the_gtsset->n_markers; i++){
    fprintf(fh, "%s ", the_gtsset->marker_ids->a[i]);
  }fprintf(fh, "\n");
  for(long i=0; i<the_gtsset->n_accessions; i++){
    fprintf(fh, "%s  %s\n", the_gtsset->accessions->a[i]->id->a, the_gtsset->accessions->a[i]->genotypes->a);
  }
}

void print_genotypesset_summary_info(FILE* fh, GenotypesSet* the_gtsset){
  fprintf(fh, "# n_accessions: %ld\n", the_gtsset->n_accessions);
  fprintf(fh, "# delta: %6.4lf\n", the_gtsset->delta);
  fprintf(fh, "# max marker missing data fraction: %6.4f\n", the_gtsset->max_marker_missing_data_fraction);
  fprintf(fh, "# n_markers: %ld\n", the_gtsset->n_markers);
}

void free_genotypesset(GenotypesSet* the_gtsset){
  free_vstr(the_gtsset->marker_ids);
  free_vlong(the_gtsset->marker_missing_data_counts);
  //  fprintf(stderr, "n_accessions: %ld \n", the_gtsset->n_accessions);
  /* for(long i=0; i<the_gtsset->n_accessions; i++){ */
  /*   free_accession_innards(the_gtsset->accessions+i); */
  /* } */
  free_vaccession(the_gtsset->accessions);
  free(the_gtsset);
}

Vidxid* construct_vidxid(const GenotypesSet* the_gtsset){
  Vidxid* the_vidxid = (Vidxid*)malloc(sizeof(Vidxid));
  the_vidxid->capacity = the_gtsset->n_accessions;
  the_vidxid->size = the_vidxid->capacity;
  the_vidxid->a = (IndexId**)malloc(the_vidxid->size*sizeof(IndexId*));
  for(long i=0; i<the_vidxid->size; i++){
    IndexId* the_idxid = (IndexId*)malloc(sizeof(IndexId));
    the_idxid->index = i;
    char* the_id =  the_gtsset->accessions->a[i]->id->a;
    the_idxid->id = strcpy((char*)malloc((strlen(the_id)+1)*sizeof(char)), the_id);
    the_vidxid->a[i] = the_idxid;
  }
  return the_vidxid;  
}

Vidxid* construct_sorted_vidxid(const GenotypesSet* the_gtsset){
  Vidxid* the_vidxid = construct_vidxid(the_gtsset);
  sort_vidxid_by_id(the_vidxid);
  for(long i=0; i<the_vidxid->size; i++){
    IndexId* x = the_vidxid->a[i];
    // fprintf(stderr, "%ld %ld %s\n", i, x->index, x->id);
  }
  if(DBUG) assert(check_idxid_map(the_vidxid, the_gtsset) == 1);
  return the_vidxid;
}

