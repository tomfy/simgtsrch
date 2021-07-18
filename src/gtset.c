#include "gtset.h"

extern int do_checks_flag; // option -c sets this to 1 to do some checks.

// *****  Accession  implementation *****
Accession* construct_accession(char* id, long idx, char* genotypes){
  Accession* the_accession = (Accession*)malloc(1*sizeof(Accession));
  the_accession->id = construct_vchar_from_str(id);
  the_accession->index = idx;
  the_accession->genotypes = construct_vchar_from_str(genotypes);
  return the_accession;
}
void set_accession_missing_data_count(Accession* the_accession, long missing_data_count){
  the_accession->missing_data_count = missing_data_count;
}
void free_accession(Accession* the_accession){
  free_accession_innards(the_accession);
  free(the_accession);
}
void free_accession_innards(Accession* the_accession){ // doesn't free the_accession itself
  // use to free each element of an array of Accession (not an array of Accession*)
  free_vchar(the_accession->id);
  free_vchar(the_accession->genotypes);
}


// *****  GenotypesSet implementation *****
GenotypesSet* read_genotypes_file_and_store(FILE* g_stream, double delta, double max_missing_data_fraction){
  char* line = NULL;
  size_t len = 0;
  ssize_t nread;

  long accessions_capacity = 10000;
  Accession* the_accessions = (Accession*)malloc(accessions_capacity*sizeof(Accession)); 

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
  
    char* token = strtok_r(line, "\t \n\r", &saveptr);
    //printf("accession number: %ld id: %s\n", accession_count, token);
    char* acc_id = (char*)malloc((strlen(token)+1)*sizeof(char));
    add_string_to_vstr(accession_ids, strcpy(acc_id, token)); // store copy of accession id
    //  printf("%s  ", token);
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
	missing_data_counts[marker_count]++;
	accession_missing_data_count++;
      }
      marker_count++;
    } // done reading dosages for all markers, and rounding to 0, 1, 2, 3
    add_string_to_vstr(genotypes_strings, genotypes);
    if(marker_count != markerid_count) exit(EXIT_FAILURE);
    // printf("%s\n", genotypes);
    the_accessions[accession_count].id = construct_vchar_from_str(acc_id);
    the_accessions[accession_count].genotypes = construct_vchar_from_str(genotypes);
    the_accessions[accession_count].missing_data_count = accession_missing_data_count;
    accession_count++;
  } // done reading all lines
  
  free(line); // only needs to be freed once.
  GenotypesSet* the_genotypes_set = //construct_genotypesset(accession_ids, marker_ids, genotypes_strings, the_md_vlong);
    construct_genotypesset(accession_count, the_accessions, marker_ids, the_md_vlong);

  // the_genotypes_set->accessions = the_accessions;
  the_genotypes_set->capacity = accessions_capacity;
  return the_genotypes_set;
}

/* void check_gtsset(GenotypesSet* gtsset){ */
/*   for(long i=0; i<gtsset->n_accessions; i++){ */
/*     //  fprintf(stderr, "id(old): %s \n", gtsset->accession_ids->a[i]); */
/*     //  fprintf(stderr, "id(new): %s \n", gtsset->accessions[i].id->a); */
/*     assert(strcmp(gtsset->accessions[i].id->a, gtsset->accession_ids->a[i]) == 0); */
/*     assert(strcmp(gtsset->accessions[i].genotypes->a, gtsset->genotype_sets->a[i]) == 0); */
/*   }   */
/* } */

/* GenotypesSet* construct_genotypesset(Vstr* acc_ids, Vstr* marker_ids, Vstr* gsets, Vlong* md_counts){ */
/*   if(gsets->size != acc_ids->size){ fprintf(stderr, "Inconsistency in construct_genotypesset\n"); exit(EXIT_FAILURE); } */
/*   GenotypesSet* the_gtsset = (GenotypesSet*)malloc(1*sizeof(GenotypesSet)); */
/*   the_gtsset->n_accessions = acc_ids->size;  */
/*   the_gtsset->n_markers = marker_ids->size; */
/*   the_gtsset->accession_ids = acc_ids; */
/*   the_gtsset->accession_missing_data_counts = construct_vlong_zeroes(the_gtsset->n_accessions); */
/*   the_gtsset->marker_ids = marker_ids; */
/*   the_gtsset->genotype_sets = gsets; */
/*   the_gtsset->marker_missing_data_counts = md_counts; */
 
/*   return the_gtsset; */
/* } */

GenotypesSet* construct_genotypesset(long n_accessions, Accession* accessions, Vstr* marker_ids, Vlong* md_counts){
  GenotypesSet* the_gtsset = (GenotypesSet*)malloc(1*sizeof(GenotypesSet));
  the_gtsset->n_accessions = n_accessions;
  the_gtsset->n_markers = marker_ids->size;
  the_gtsset->accessions = accessions;
  //  the_gtsset->accession_missing_data_counts = construct_vlong_zeroes(the_gtsset->n_accessions);
  the_gtsset->marker_ids = marker_ids;
  //  the_gtsset->genotype_sets = gsets;
  the_gtsset->marker_missing_data_counts = md_counts;
 
  return the_gtsset;
}

void check_genotypesset(GenotypesSet* gtss, double max_marker_md_fraction){
  /* assert(gtss->accession_ids->size == gtss->n_accessions); */
  /* assert(gtss->genotype_sets->size == gtss->n_accessions); */
  assert(gtss->marker_ids->size == gtss->n_markers);
  long* md_counts = (long*)calloc(gtss->n_markers, sizeof(long)); 
  for(long i=0; i<gtss->n_accessions; i++){
    Accession* an_acc = gtss->accessions + i;
    assert(strlen(an_acc->genotypes->a) == gtss->n_markers);
    for(long j=0; j<gtss->n_markers; j++){
      if(an_acc->genotypes->a[j] == '3')md_counts[j]++;
    }
  }
  for(long j=0; j<gtss->n_markers; j++){
    assert(md_counts[j] == gtss->marker_missing_data_counts->a[j]);
  }
  free(md_counts);
  // set_accession_missing_data_counts(gtss);
  fprintf(stderr, "Successfully completed check_genotypesset\n");
}

GenotypesSet* construct_cleaned_genotypesset(GenotypesSet* the_gtsset, double max_md_fraction){
 
  //  Vstr* gsets = the_gtsset->genotype_sets;
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

  Accession* the_accessions = (Accession*)malloc(the_gtsset->n_accessions*sizeof(Accession)); 

  //  Vstr* copy_of_accids = construct_vstr_copy(the_gtsset->accession_ids);
  Vstr* cleaned_gsets = construct_vstr(the_gtsset->n_accessions);
  for(long i=0; i<the_gtsset->n_accessions; i++){ // loop over accessions
    char* raw_gts = the_gtsset->accessions[i].genotypes->a; // the string with all the genotypes for accession i
    char* cleaned_gts = (char*)malloc((n_markers_to_keep+1)*sizeof(char));
    //char* cleaned_gts2 = (char*)malloc((n_markers_to_keep+1)*sizeof(char));
    // char* cleaned_marker_ids = (char*)malloc((n_markers_to_keep+1)*sizeof(char));
    long k=0; // k: index of kept markers
    long acc_md_count = 0;
    for(long j=0; j<the_gtsset->n_markers; j++){ // j: index of original markers
      if(md_ok->a[j] == 1){
	cleaned_gts[k] = raw_gts[j];
	//cleaned_gts2[k] = raw_gts[j];
	//	long len = strlen(the_gtsset->marker_ids[j]);
	// cleaned_marker_ids[k] = strcpy((char*)malloc((len+1)*sizeof(char), the_gtsset->marker_ids[j]);
	k++;
	if(raw_gts[j] == '3') acc_md_count++;
      }
    }
    cleaned_gts[k] = '\0'; // terminate with null.
    // cleaned_gts2[k] = '\0';
    //  add_string_to_vstr(cleaned_gsets, cleaned_gts);
    // free(the_gtsset->accessions[i].genotypes->a);
   
    if(DBUG && do_checks_flag) assert(k == n_markers_to_keep);

    the_accessions[i].id = construct_vchar_from_str(the_gtsset->accessions[i].id->a);
    the_accessions[i].index = i; 
    the_accessions[i].genotypes = construct_vchar_from_str(cleaned_gts);
    the_accessions[i].missing_data_count = acc_md_count;
    //free(cleaned_gts2);
  }
  free_vlong(md_ok);
  //  GenotypesSet* construct_genotypesset(Vstr* acc_ids, Vstr* marker_ids, Vstr* gsets, Vlong* md_counts){
  GenotypesSet* cleaned_gtsset = //construct_genotypesset(copy_of_accids, cleaned_marker_ids, cleaned_gsets, cleaned_md_counts);
    construct_genotypesset(the_gtsset->n_accessions, the_accessions, cleaned_marker_ids, cleaned_md_counts);
  // set_accession_missing_data_counts(cleaned_gtsset);

  // cleaned_gtsset->accessions = the_accessions;
  cleaned_gtsset->capacity = the_gtsset->n_accessions;

  // the_gtsset->accessions[i].genotypes->a = cleaned_gts2;
  fprintf(stderr, "about to return cleaned_gtsset \n");
  return cleaned_gtsset;
}

/* void set_accession_missing_data_counts(GenotypesSet* the_gtsset){ */
/*   // fprintf(stderr, "In set_accession_missing_data_counts. \n"); */
/*   for(long i=0; i<the_gtsset->accession_ids->size; i++){ */
/*     long md_count = 0; */
/*     char* the_gtsstr = the_gtsset->genotype_sets->a[i]; */
/*     for(long j=0; the_gtsstr[j] != '\0'; j++){ */
/*       if(the_gtsstr[j] == '3'){ */
/* 	md_count++; */
/*       } */
/*     } */
/*     //  fprintf(stderr, "accid %s acc index: %ld, md_count: %ld\n", the_gtsset->accession_ids->a[i], i, md_count); */
/*     the_gtsset->accession_missing_data_counts->a[i] = md_count; */
/*   } */
/* } */

/* void print_genotypesset(FILE* fh, GenotypesSet* the_gtsset){ */
/*   printf("MARKER  "); */
/*   for(long i=0; i<the_gtsset->n_markers; i++){ */
/*     fprintf(fh, "%s ", the_gtsset->marker_ids->a[i]); */
/*   }fprintf(fh, "\n"); */
/*   for(long i=0; i<the_gtsset->n_accessions; i++){ */
/*     fprintf(fh, "%s  %s\n", the_gtsset->accession_ids->a[i], the_gtsset->genotype_sets->a[i]); */
/*   } */
/* } */

void print_genotypesset_summary_info(FILE* fh, GenotypesSet* the_gtsset){
  fprintf(fh, "# n_accessions: %ld\n", the_gtsset->n_accessions);
  fprintf(fh, "# n_markers: %ld\n", the_gtsset->n_markers);
}

void free_genotypesset(GenotypesSet* the_gtsset){
  //  free_vstr(the_gtsset->accession_ids);
  free_vstr(the_gtsset->marker_ids);
  // free_vstr(the_gtsset->genotype_sets);
  //  free_vlong(the_gtsset->accession_missing_data_counts);
  free_vlong(the_gtsset->marker_missing_data_counts);
  fprintf(stderr, "n_accessions: %ld \n", the_gtsset->n_accessions);
  for(long i=0; i<the_gtsset->n_accessions; i++){
    free_accession_innards(the_gtsset->accessions+i);
  }
  free(the_gtsset->accessions);
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
    char* the_id =  the_gtsset->accessions[i].id->a;
    the_idxid->id = strcpy((char*)malloc((strlen(the_id)+1)*sizeof(char)), the_id);
    the_vidxid->a[i] = the_idxid;
  }
  return the_vidxid;  
}

Vidxid* construct_sorted_vidxid(const GenotypesSet* the_gtsset){
  Vidxid* the_vidxid = construct_vidxid(the_gtsset);
  sort_vidxid_by_id(the_vidxid);
   if(DBUG) assert(check_idxid_map(the_vidxid, the_gtsset) == 1);
  return the_vidxid;
}

long check_idxid_map(Vidxid* vidxid, const GenotypesSet* the_gtsset){
  for(long i=0; i<the_gtsset->n_accessions; i++){
    char* id = the_gtsset->accessions[i].id->a;
    long idx = index_of_id_in_vidxid(vidxid, id);
    if(idx != i) return 0;
  }
  return 1;
}
