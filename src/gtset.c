#include "gtset.h"


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

GenotypesSet* construct_genotypesset(Vstr* acc_ids, Vstr* marker_ids, Vstr* gsets, Vlong* md_counts){
  if(gsets->size != acc_ids->size){ fprintf(stderr, "Inconsistency in construct_genotypesset\n"); exit(EXIT_FAILURE); }
  GenotypesSet* the_gtsset = (GenotypesSet*)malloc(1*sizeof(GenotypesSet));
  the_gtsset->n_accessions = acc_ids->size; 
  the_gtsset->n_markers = marker_ids->size;
  the_gtsset->accession_ids = acc_ids;
  the_gtsset->accession_missing_data_counts = construct_vlong_zeroes(the_gtsset->n_accessions);
  the_gtsset->marker_ids = marker_ids;
  the_gtsset->genotype_sets = gsets;
  the_gtsset->marker_missing_data_counts = md_counts;
  return the_gtsset;
}

void check_genotypesset(GenotypesSet* gtss, double max_marker_md_fraction){
  assert(gtss->accession_ids->size == gtss->n_accessions);
  assert(gtss->genotype_sets->size == gtss->n_accessions);
  assert(gtss->marker_ids->size == gtss->n_markers);
  long* md_counts = (long*)calloc(gtss->n_markers, sizeof(long)); 
  for(long i=0; i<gtss->genotype_sets->size; i++){
    assert(strlen(gtss->genotype_sets->a[i]) == gtss->n_markers);
    for(long j=0; j<gtss->n_markers; j++){
      if(gtss->genotype_sets->a[i][j] == '3')md_counts[j]++;
    }
  }
  for(long j=0; j<gtss->n_markers; j++){
    assert(md_counts[j] == gtss->marker_missing_data_counts->a[j]);
  }
  free(md_counts);
  set_accession_missing_data_counts(gtss);
  fprintf(stderr, "Successfully completed check_genotypesset\n");
}

GenotypesSet* construct_cleaned_genotypesset(GenotypesSet* the_gtsset, double max_md_fraction){
 
  Vstr* gsets = the_gtsset->genotype_sets;
  Vlong* md_counts = the_gtsset->marker_missing_data_counts;

  // identify the markers to keep:
  long n_markers_to_keep = 0;
  Vlong* md_ok = construct_vlong_zeroes(md_counts->size);
  Vstr* cleaned_marker_ids = construct_vstr(1000);
  Vlong* cleaned_md_counts = construct_vlong(1000);
  for(long i=0; i<md_counts->size; i++){
    if(md_counts->a[i] <= max_md_fraction*the_gtsset->n_accessions){
      md_ok->a[i] = 1;
      n_markers_to_keep++;
      long marker_id_length = strlen(the_gtsset->marker_ids->a[i]);
      char* marker_id_to_keep = strcpy((char*)malloc((marker_id_length+1)*sizeof(char)), the_gtsset->marker_ids->a[i]);
      add_string_to_vstr(cleaned_marker_ids, marker_id_to_keep);
      add_long_to_vlong(cleaned_md_counts, md_counts->a[i]);
    }
  }
  // fprintf(stderr, "after 1st loop\n");
  //GenotypesSet* the_cleaned_gtsset = (GenotypesSet*)malloc(sizeof(GenotypesSet));
  Vstr* copy_of_accids = construct_vstr_copy(the_gtsset->accession_ids);
  Vstr* cleaned_gsets = construct_vstr(the_gtsset->n_accessions);
  for(long i=0; i<gsets->size; i++){ // loop over accessions
    char* raw_gts = gsets->a[i]; // the string with all the genotypes for accession i
    char* cleaned_gts = (char*)malloc((n_markers_to_keep+1)*sizeof(char));
    // char* cleaned_marker_ids = (char*)malloc((n_markers_to_keep+1)*sizeof(char));
    long k=0; // k: index of kept markers
    for(long j=0; j<the_gtsset->n_markers; j++){ // j: index of original markers
      if(md_ok->a[j] == 1){
	cleaned_gts[k] = raw_gts[j];
	//	long len = strlen(the_gtsset->marker_ids[j]);
	// cleaned_marker_ids[k] = strcpy((char*)malloc((len+1)*sizeof(char), the_gtsset->marker_ids[j]);
	k++;
      }
    }
    cleaned_gts[k] = '\0'; // terminate with null.
    add_string_to_vstr(cleaned_gsets, cleaned_gts);
    if(DBUG) assert(k == n_markers_to_keep);
  }
  free_vlong(md_ok);
   set_accession_missing_data_counts(the_gtsset);
  //  GenotypesSet* construct_genotypesset(Vstr* acc_ids, Vstr* marker_ids, Vstr* gsets, Vlong* md_counts){
  GenotypesSet* cleaned_gtsset = construct_genotypesset(copy_of_accids, cleaned_marker_ids, cleaned_gsets, cleaned_md_counts);
  return cleaned_gtsset;
}

void set_accession_missing_data_counts(GenotypesSet* the_gtsset){
  for(long i=0; i<the_gtsset->accession_ids->size; i++){
    long md_count = 0;
   char* the_gtsstr = the_gtsset->genotype_sets->a[i];
    for(long j=0; the_gtsstr[j] != '\0'; j++){
      if(the_gtsstr[j] == '3'){
	md_count++;
      }
    }
    the_gtsset->accession_missing_data_counts->a[i] = md_count;
  }
}

void print_genotypesset(FILE* fh, GenotypesSet* the_gtsset){
  printf("MARKER  ");
  for(long i=0; i<the_gtsset->n_markers; i++){
    fprintf(fh, "%s ", the_gtsset->marker_ids->a[i]);
  }fprintf(fh, "\n");
  for(long i=0; i<the_gtsset->n_accessions; i++){
    fprintf(fh, "%s  %s\n", the_gtsset->accession_ids->a[i], the_gtsset->genotype_sets->a[i]);
  }
}

void print_genotypesset_summary_info(FILE* fh, GenotypesSet* the_gtsset){
  fprintf(fh, "# n_accessions: %ld\n", the_gtsset->n_accessions);
  fprintf(fh, "# n_markers: %ld\n", the_gtsset->n_markers);
}

void free_genotypesset(GenotypesSet* the_gtsset){
  free_vstr(the_gtsset->accession_ids);
  free_vstr(the_gtsset->marker_ids);
  free_vstr(the_gtsset->genotype_sets);
  free_vlong(the_gtsset->marker_missing_data_counts);
  free(the_gtsset);
}
