#include "gtset.h"

/* typedef struct{ */
/*   long n_accessions; // */
/*   long n_markers; // */
/*   Vstr* accession_ids; // array of accession_ids */
/*   Vstr* genotype_sets; // */
/*   Vstr* marker_ids; // array of marker_ids */
/*   Vlong* marker_missing_data_counts; // */
/* }GenotypesSet; */



GenotypesSet* construct_genotypesset(Vstr* acc_ids, Vstr* marker_ids, Vstr* gsets, Vlong* md_counts){
  if(gsets->size != acc_ids->size){ fprintf(stderr, "Inconsistency in construct_genotypesset\n"); exit(EXIT_FAILURE); }
  GenotypesSet* the_gtsset = (GenotypesSet*)malloc(1*sizeof(GenotypesSet));
  the_gtsset->n_accessions = acc_ids->size; 
  the_gtsset->n_markers = marker_ids->size;
  the_gtsset->accession_ids = acc_ids;
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
  //  GenotypesSet* construct_genotypesset(Vstr* acc_ids, Vstr* marker_ids, Vstr* gsets, Vlong* md_counts){
  GenotypesSet* cleaned_gtsset = construct_genotypesset(copy_of_accids, cleaned_marker_ids, cleaned_gsets, cleaned_md_counts);
  return cleaned_gtsset;
}


void print_genotypesset(GenotypesSet* the_gtsset){
  printf("MARKER  ");
  for(long i=0; i<the_gtsset->n_markers; i++){
    printf("%s ", the_gtsset->marker_ids->a[i]);
  }printf("\n");
  for(long i=0; i<the_gtsset->n_accessions; i++){
    printf("%s  %s\n", the_gtsset->accession_ids->a[i], the_gtsset->genotype_sets->a[i]);
  }
}

void print_genotypesset_summary_info(GenotypesSet* the_gtsset){
  fprintf(stderr, "# n_accessions: %ld\n", the_gtsset->n_accessions);
  fprintf(stderr, "# n_markers: %ld\n", the_gtsset->n_markers);
}

void free_genotypesset(GenotypesSet* the_gtsset){
  free_vstr(the_gtsset->accession_ids);
  free_vstr(the_gtsset->marker_ids);
  free_vstr(the_gtsset->genotype_sets);
  free_vlong(the_gtsset->marker_missing_data_counts);
  free(the_gtsset);
}
