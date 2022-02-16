// C version of program to test pedigrees using genotype data.
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>
#include "gtset.h"
#include "pedigree.h"
#define UNKNOWN -1
#define DOSAGES 0
#define GENOTYPES 1

int do_checks_flag = 0; // option -c sets this to 1 to do some checks.

double hi_res_time(void);

// **********************************************************************************************
// ***********************************  main  ***************************************************
// **********************************************************************************************


int
main(int argc, char *argv[])
{
  double t_begin_main = hi_res_time();
  
  int do_alternative_pedigrees = 0; // 0: none, 1: only when given pedigree is 'bad', 2: all
  double delta = 0.05; // default; control this with -d command line option.
  double max_marker_missing_data = 0.2; // default; control this with -x command line option.
  char* pedigree_test_output_filename = "pedigree_test_info";
  char* genotypes_matrix_output_filename = "genotype_matrix_out";
  double max_self_agmr12 = 1; // need to specify if doing alternative pedigrees 
  double max_ok_hgmr = 1; // accept everything as ok
  double max_self_r = 1; // need to specify if doing alternative pedigrees 
  double max_ok_d = 1; // accept everything as ok
  // ***** process command line *****
  if (argc < 2) {
    fprintf(stderr, "Usage:  %s -g <genotypes_file>  -p <pedigree_file> \n", argv[0]);
    exit(EXIT_FAILURE);
  }
  //  char* dosages_filename = NULL;
  char* genotypes_filename = NULL;
  FILE *g_stream = NULL;
  char* pedigrees_filename = NULL;
  FILE *p_stream = NULL;

  // A: pedigree alternatives 0, 1, or 2
  // c: do checks,
  // d: dosages filename, g: genotypes filename (must have either d or g)
  // p: pedigree filename,
  // w: (width for rounding is +-w),
  // x: max fraction of missing data for markers,
  // o: output filename.
  // a: max 'self' agmr, h: max ok hgmr,  r: max 'self' r, D: max ok d;
  int c;
  int genotype_file_type = UNKNOWN;
  while((c = getopt(argc, argv, "A:cd:g:p:w:x:o:a:h:r:D:")) != -1){
    switch(c){
    case 'A':
      if(optarg == 0){
	perror("option A requires an integer argument; 0, 1, or >=2\n");
	exit(EXIT_FAILURE);
      }else{
	do_alternative_pedigrees = atoi(optarg);
	if(do_alternative_pedigrees < 0) exit(EXIT_FAILURE);
	if(do_alternative_pedigrees > 2) do_alternative_pedigrees = 2;
      }
      break;
    case 'c':
      do_checks_flag = 1;
      break;
    case 'd':
      genotypes_filename = optarg;
      g_stream = fopen(genotypes_filename, "r");
      if(g_stream == NULL){
	fprintf(stderr, "Failed to open %s for reading.\n", genotypes_filename);
	exit(EXIT_FAILURE);
      }
      genotype_file_type = DOSAGES;
      break;
    case 'g':
      genotypes_filename = optarg;
      g_stream = fopen(genotypes_filename, "r");
      if(g_stream == NULL){
	fprintf(stderr, "Failed to open %s for reading.\n", genotypes_filename);
	exit(EXIT_FAILURE);

      }
      genotype_file_type = GENOTYPES;
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
      pedigree_test_output_filename = optarg;
      break;
    case 'w':
      if(optarg == 0){
	perror("option w requires a numerical argument > 0\n");
	exit(EXIT_FAILURE);
      }else{
	delta = atof(optarg);
	if(delta < 0) exit(EXIT_FAILURE);
      }
      break;
    case 'x':
      if(optarg == 0){
	perror("option x requires a numerical argument > 0\n");
	exit(EXIT_FAILURE);
      }else{
	max_marker_missing_data = atof(optarg);
	if (max_marker_missing_data < 0){
	  printf("# max missing data fraction in markers will be set automatically.\n");
	}}
      break;
    case 'a': // agmr12 < this -> looks like self
      if(optarg == 0){
	perror("option a requires a numerical argument > 0\n");
	exit(EXIT_FAILURE);
      }else{
	max_self_agmr12 = atof(optarg);
	if (max_self_agmr12 < 0) exit(EXIT_FAILURE);
      }
    case 'h': // hgmr > this -> poor parent-offspring candidate
      if(optarg == 0){
	perror("option h requires a numerical argument > 0\n");
	exit(EXIT_FAILURE);
      }else{
	max_ok_hgmr = atof(optarg);
	if (max_ok_hgmr < 0) exit(EXIT_FAILURE);
      }
      break;
    case 'r': // r > this -> probably biparental
      if(optarg == 0){
	perror("option x requires a numerical argument > 0\n");
	exit(EXIT_FAILURE);
      }else{
	max_self_r = atof(optarg);
	if (max_self_r < 0) exit(EXIT_FAILURE);
      }
      break;
    case 'D': // d > this argument means this a is poor candidate triple of parents and offspring
      if(optarg == 0){
	perror("option x requires a numerical argument > 0\n");
	exit(EXIT_FAILURE);
      }else{
	max_ok_d = atof(optarg);
	if (max_ok_d < 0) exit(EXIT_FAILURE);
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

  // exit(EXIT_FAILURE);
  
  FILE *o_stream = NULL;
  o_stream = fopen(pedigree_test_output_filename, "w");
  if(o_stream == NULL){
    fprintf(stderr, "Failed to open %s for writing.\n", pedigree_test_output_filename);
    exit(EXIT_FAILURE);
  }
  // fprintf(stderr, "# genotypes file type: %d\n", genotype_file_type);
  char* geno_file_type = (genotype_file_type == DOSAGES)? "dosages" : "genotypes";
  fprintf(stderr, "# genotypes file type/name: %s / %s.  delta: %5.3lf  max marker missing data: %5.3lf\n",
	  geno_file_type, genotypes_filename, delta, max_marker_missing_data);
  fprintf(stderr, "# pedigrees filename: %s, output filename: %s \n", pedigrees_filename, pedigree_test_output_filename);

  
  // *****  done processing command line  *****

  // ***************  read the genotypes file  *******************************
  double t_start = hi_res_time();
  GenotypesSet* the_genotypes_set;
  if(genotype_file_type == DOSAGES){
    //  fprintf(stderr, "# dosages branch \n");
    GenotypesSet* the_raw_genotypes_set = read_dosages_file_and_store(g_stream, delta);
  
    fclose(g_stream);
    // print_genotypesset_summary_info(stderr, the_raw_genotypes_set); 
   
    if(DBUG && do_checks_flag) check_genotypesset(the_raw_genotypes_set, max_marker_missing_data);
  
    fprintf(stderr, "# Done reading in genotype data. %ld accessions and %ld markers. Time: %10.4lf sec.\n",
	    the_raw_genotypes_set->n_accessions, the_raw_genotypes_set->n_markers, hi_res_time() - t_start);
 
    // *****  clean genotypes set, i.e. remove markers with high missing data  ****
    t_start = hi_res_time();
    the_genotypes_set = construct_cleaned_genotypesset(the_raw_genotypes_set, max_marker_missing_data);

    
    //  check_gtsset(the_cleaned_genotypes_set);
    long n_markers_good = the_genotypes_set->n_markers;
    free_genotypesset(the_raw_genotypes_set); // free the raw genotypes set; use the cleaned one.
    print_genotypesset_summary_info(stderr, the_genotypes_set);

    //  fprintf(stderr, "# Removing markers with excessive fraction (> %5.3lf) missing data.\n", max_marker_missing_data);
    if(DBUG && do_checks_flag) check_genotypesset(the_genotypes_set, max_marker_missing_data);
    //   fprintf(stderr, "# Done cleaning marker set. Keeping %ld markers. Time: %10.4lf sec.\n", n_markers_good, hi_res_time() - t_start);
    FILE* og_stream = fopen(genotypes_matrix_output_filename, "w");
    print_genotypesset(og_stream, the_genotypes_set);
    fclose(og_stream);
  }else if(genotype_file_type == GENOTYPES){
    the_genotypes_set = read_genotypes_file_and_store(g_stream);
  }else{
    fprintf(stderr, "genotype_file_type: %d (should be 0 or 1)\n", genotype_file_type);
    exit(EXIT_FAILURE); // no input file specified.
  }
 
  
  t_start = hi_res_time();
  Vidxid* the_vidxid = construct_sorted_vidxid(the_genotypes_set);
  //  fprintf(stderr, "# Time to set up id index map: %10.4lf \n", hi_res_time() - t_start);
  // ***************  read the pedigrees file  ***************************
  t_start = hi_res_time();
  const Vpedigree* pedigrees = read_the_pedigrees_file_and_store(p_stream, the_vidxid, the_genotypes_set);
  fclose(p_stream);
  fprintf(stderr, "# Done reading pedigree file. Stored %ld pedigrees.\n", pedigrees->size);

  // ***************  Done reading input files  ******************************

  
  const Vlong* parent_idxs = accessions_with_offspring(pedigrees, the_genotypes_set->n_accessions);
  printf("# According to pedigree file there are %ld accession with offspring.\n", parent_idxs->size);

  t_start = hi_res_time();
  /* long* nd0 = (long*)calloc(the_genotypes_set->n_markers, sizeof(long)); */
  /* long* nd1 = (long*)calloc(the_genotypes_set->n_markers, sizeof(long)); */
  /* long* nd2 = (long*)calloc(the_genotypes_set->n_markers, sizeof(long)); */
  /* long n_markers; */
  /* long n_ok_pedigrees = 0; */
  for(long i=0; i<pedigrees->size; i++){
    if(i % 1000  == 0){
      fprintf(stderr, "# Done testing %ld pedigrees.\n", i);
    }
    Pedigree_stats* the_pedigree_stats = calculate_pedigree_stats(pedigrees->a[i]); //, nd0, nd1, nd2); //, the_cleaned_genotypes_set);
    // assert(strcmp(pedigrees->a[i]->Accession->id, pedigrees->a[i]->A->id->a) == 0);
    fprintf(o_stream, "%20s %5ld %20s %20s  ",
	    pedigrees->a[i]->A->id->a, pedigrees->a[i]->A->missing_data_count,
	    pedigrees->a[i]->F->id->a, pedigrees->a[i]->M->id->a);
    print_pedigree_stats(o_stream, the_pedigree_stats);
    /* if(0 && get_d(the_pedigree_stats) < 0.03){ */
    /* n_markers = marker_d_counts(pedigrees->a[i], nd0, nd1, nd2); */
    /* n_ok_pedigrees++; */
    /* } */
    if(do_alternative_pedigrees > 0){
      if((do_alternative_pedigrees == 2) || (pedigree_ok(the_pedigree_stats, max_self_agmr12, max_ok_hgmr, max_self_r, max_ok_d) == 0)){
	//   fprintf(stderr, "About to get alt. pedigrees. max_ok_hgmr: %8.4lf,  max_ok_d: %8.4lf \n", max_ok_hgmr, max_ok_d);
	Vpedigree* alt_pedigrees = pedigree_alternatives(pedigrees->a[i], the_genotypes_set, parent_idxs, max_ok_hgmr, max_ok_d);
	print_pedigree_alternatives(o_stream, alt_pedigrees);
      }
    }
    fprintf(o_stream, "\n");
    free(the_pedigree_stats);
    // if(i == 500) break; // for testing purposes only.
  }
  /* for(long i=0; i< the_genotypes_set->n_markers; i++){ */
  /*   double denomx = (double)(nd0[i]+nd1[i]+nd2[i]); */
  /*   printf("%ld %ld %ld  %8.5lf %8.5lf %8.5lf %8.5lf \n", n_ok_pedigrees, n_ok_pedigrees - (nd0[i]+nd1[i]+nd2[i]), */
  /* 	   i, denomx, (double)nd0[i]/denomx, (double)nd1[i]/denomx, (double)nd2[i]/denomx); */
  /* } */

  
  fprintf(stderr, "# Done testing %ld pedigrees.\n", pedigrees->size);

  // ********************  cleanup  **************************
  t_start = hi_res_time();
  fclose(o_stream);
  free_genotypesset(the_genotypes_set);
  //  fprintf(stderr, "after free_genotypesset\n");
  free_vidxid(the_vidxid);
  free_vpedigree(pedigrees);
  //  fprintf(stderr, "after free_vpedigree\n");
  free_vlong(parent_idxs);
  // fprintf(stderr, "# Done with cleanup\n");
  // getchar();
  fprintf(stderr, "# Total time: %10.4lf sec.\n", hi_res_time() - t_begin_main);
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


