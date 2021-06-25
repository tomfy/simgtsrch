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
#include "pedigree.h"

int do_checks_flag = 0; // option -c sets this to 1 to do some checks.

double hi_res_time(void);

// **********************************************************************************************
// ***********************************  main  ***************************************************
// **********************************************************************************************

int
main(int argc, char *argv[])
{
  double delta = 0.05; // default; control this with -d command line option.
  double max_marker_missing_data = 0.2; // default; control this with -x command line option.
  char* output_filename = "pedigree_check_info";

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
  while((c = getopt(argc, argv, "cg:p:d:x:o:")) != -1){
    // fprintf(stderr, "c: %c %s %d\n", c, optarg, optind);
    switch(c){
    case 'c':
      do_checks_flag = 1;
      break;
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
    case 'o':
      output_filename = optarg;
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
      if ((optopt == 'g') || (optopt == 'p') || (optopt == 'd') || (optopt == 'x') || (optopt == 'o'))
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
  
  FILE *o_stream = NULL;
  o_stream = fopen(output_filename, "w");
  if(o_stream == NULL){
    fprintf(stderr, "Failed to open %s for writing.\n", output_filename);
    exit(EXIT_FAILURE);
  }
      
  fprintf(stderr, "# genotypes file: %s  pedigree file: %s  delta: %5.3lf  max marker missing data: %5.3lf  output file: %s\n",
	  genotypes_filename, pedigrees_filename, delta, max_marker_missing_data, output_filename);

  // *****  done processing command line  *****

  // ***************  read the genotypes file  *******************************
  double t_start = hi_res_time();
  GenotypesSet* the_genotypes_set = read_genotypes_file_and_store(g_stream, delta, max_marker_missing_data);
  fclose(g_stream);

  print_genotypesset_summary_info(stderr, the_genotypes_set);
  if(DBUG && do_checks_flag) check_genotypesset(the_genotypes_set, max_marker_missing_data);
  
  long n_accessions = the_genotypes_set->n_accessions;
  long n_markers_all = the_genotypes_set->n_markers;
  
  fprintf(stderr, "Done reading in genotype data. %ld accession and %ld markers. Time: %lf sec.\n", n_accessions, n_markers_all, hi_res_time() - t_start);

  t_start = hi_res_time();
  Vidxid* the_vidxid = construct_sorted_vidxid_from_vstr(the_genotypes_set->accession_ids);
  fprintf(stderr, "Time to set up id index map: %lf \n", hi_res_time() - t_start);
 

  // ***************  read the pedigrees file  ***************************

  t_start = hi_res_time();
  Vpedigree* pedigrees = read_the_pedigrees_file_and_store(p_stream, the_vidxid);
  fclose(p_stream);
  fprintf(stderr, "Done reading pedigree file. Stored %ld pedigrees. Time: %lf sec.\n",
	  pedigrees->size, hi_res_time() - t_start);

  // ***************  Done reading input files  ******************************
  

  // *****  clean genotypes set, i.e. remove markers with high missing data  ****
  t_start = hi_res_time();
  GenotypesSet* the_cleaned_genotypes_set = construct_cleaned_genotypesset(the_genotypes_set, max_marker_missing_data);
  long n_markers_good = the_cleaned_genotypes_set->n_markers;
  free_genotypesset(the_genotypes_set);
  print_genotypesset_summary_info(stderr, the_cleaned_genotypes_set);
  // print_genotypesset(the_cleaned_genotypes_set);
 
  if(DBUG && do_checks_flag) check_genotypesset(the_cleaned_genotypes_set, max_marker_missing_data);
  fprintf(stderr, "Done cleaning marker set. Keeping %ld markers. Time: %lf sec.\n",
	  n_markers_good, hi_res_time() - t_start);

  t_start = hi_res_time();
  for(long i=0; i<pedigrees->size; i++){
    calculate_pedigree_test_info(pedigrees->a[i], the_cleaned_genotypes_set);
    print_pedigree_test_info(o_stream, pedigrees->a[i], the_cleaned_genotypes_set);
  }
  
  fprintf(stderr, "Done checking %ld pedigrees. Time: %lf sec.\n",
	  pedigrees->size, hi_res_time() - t_start);

  

  // ********************  cleanup  **************************
  t_start = hi_res_time();
  fclose(o_stream);
  free_genotypesset(the_cleaned_genotypes_set);
  free_vidxid(the_vidxid);
  free_vpedigree(pedigrees);
  fprintf(stderr, "Done with cleanup. Time %lf sec.\n", hi_res_time() - t_start);
  // getchar();
}
// **********************************************************
// ********************  end of main  ***********************
// **********************************************************




// **********************************************************
// ******************  functions  ***************************
// **********************************************************

double hi_res_time(void){
  return (double)clock()/(double)CLOCKS_PER_SEC;
}


