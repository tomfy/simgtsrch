#include <stdlib.h>
#include <stdio.h>
#include "gtset.h"
// #include "vect.h"
#include "pedigree.h"

#define PEDIGREE_FIELDS 7 // number of whitespace separated fields in pedigree file, with ids in last 3 fields.

// extern int do_checks_flag; // option -c sets this to 1 to do some checks.


// *****  Ahr  *****
void print_Ahr(FILE* fh, Ahr the_ahr){
  fprintf(fh, "%4ld %7.5lf  %4ld %7.5lf  %4ld %7.5lf",
	  the_ahr.a.d, (the_ahr.a.d > 0)? (double)the_ahr.a.n/(double)the_ahr.a.d : -1,
	  the_ahr.h.d, (the_ahr.h.d > 0)? (double)the_ahr.h.n/(double)the_ahr.h.d : -1,
	  the_ahr.r.d, (the_ahr.r.d > 0)? (double)the_ahr.r.n/(double)the_ahr.r.d : -1
	  );
}

// *****  Ahr  *****
void print_Ahrx(FILE* fh, Ahr the_ahr){
  fprintf(fh, "%6.4lf %6.4lf %6.4lf",
	  (the_ahr.a.d > 0)? (double)the_ahr.a.n/(double)the_ahr.a.d : -1,
	  (the_ahr.h.d > 0)? (double)the_ahr.h.n/(double)the_ahr.h.d : -1,
	  (the_ahr.r.d > 0)? (double)the_ahr.r.n/(double)the_ahr.r.d : -1
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
  return (n_denom > 0)? (double)n_numer/(double)n_denom : -1.0;  
}

void triple_counts_x(char* gts1, char* gts2, char* proggts){
  char c1, c2, c3;
  long n_00_0 = 0, n_00_1 = 0, n_00_2 = 0;
  long n_01_0 = 0, n_01_1 = 0, n_01_2 = 0;
  long n_02_0 = 0, n_02_1 = 0, n_02_2 = 0;
  long n_10_0 = 0, n_10_1 = 0, n_10_2 = 0;
  long n_11_0 = 0, n_11_1 = 0, n_11_2 = 0;
  long n_12_0 = 0, n_12_1 = 0, n_12_2 = 0;
  long n_20_0 = 0, n_20_1 = 0, n_20_2 = 0;
  long n_21_0 = 0, n_21_1 = 0, n_21_2 = 0;
  long n_22_0 = 0, n_22_1 = 0, n_22_2 = 0;
 
  long i=0;
  while((c1 = gts1[i]) != '\0'){
    if(c1 != '3'){ // not missing data
      c2 = gts2[i];
      if(c2 != '3'){ // not missing data
	c3 = proggts[i];
	if(c3 != '3'){
	  if(c1 == '0'){ // 0xy
	    if(c2 == '0'){
	      if(c3 == '0'){
		n_00_0++;
	      }else if(c3 == '1'){
		n_00_1++;
	      }else if(c3 == '2'){
		n_00_2++;
	      }
	    }else if(c2 == '1'){ 
	      if(c3 == '0'){ // 010
		n_01_0++;
	      }else if(c3 == '1'){
		n_01_1++;
	      }else if(c3 == '2'){
		n_01_2++;
	      }
	    }else if(c2 == '2'){
	      if(c3 == '0'){
		n_02_0++;
	      }else if(c3 == '1'){
		n_02_1++;
	      }else if(c3 == '2'){
		n_02_2++;
	      }
	    }
	  }else if(c1 == '1'){ // 1xy
	    if(c2 == '0'){
	      if(c3 == '0'){
		n_10_0++;
	      }else if(c3 == '1'){
		n_10_1++;
	      }else if(c3 == '2'){
		n_10_2++;
	      }
	    }else if(c2 == '1'){
	      if(c3 == '0'){
		n_11_0++;
	      }else if(c3 == '1'){
		n_11_1++;
	      }else if(c3 == '2'){
		n_11_2++;
	      }
	    }else if(c2 == '2'){
	      if(c3 == '0'){
		n_12_0++;
	      }else if(c3 == '1'){
		n_12_1++;
	      }else if(c3 == '2'){
		n_12_2++;
	      }
	    }
	  }else if(c1 == '2'){ // 2xy
	    if(c2 == '0'){
	      if(c3 == '0'){
		n_20_0++;
	      }else if(c3 == '1'){
		n_20_1++;
	      }else if(c3 == '2'){
		n_20_2++;
	      }
	    }else if(c2 == '1'){
	      if(c3 == '0'){
		n_21_0++;
	      }else if(c3 == '1'){
		n_21_1++;
	      }else if(c3 == '2'){
		n_21_2++;
	      }
	    }else if(c2 == '2'){
	      if(c3 == '0'){
		n_22_0++;
	      }else if(c3 == '1'){
		n_22_1++;
	      }else if(c3 == '2'){
		n_22_2++;
	      }
	    }
	  }
	}
      }
    }
    i++;
  }
  long n_0 = n_00_0 + n_01_0 + n_01_1 + n_02_1 + n_10_0 + n_10_1 +
    n_12_1 + n_12_2 + n_20_1 + n_21_1 + n_21_2 + n_22_2; // 12
  long n_1 = n_00_1 + n_01_2 + n_02_0 + n_02_2 + n_10_2 + n_12_0 +
    n_20_0 + n_20_2 + n_21_0 + n_22_1;  // 10
  long n_2 = n_00_2 + n_22_0;
  // not using n_11_x
  double d1 = (n_0 > 0)? (double)n_1/(double)(n_0 + n_1) : 2;
  double d2 = (n_0 > 0)? (double)n_2/(double)(n_0 + n_2) : 2;
  long hgmr1_numer = n_00_2 + n_01_2 + n_02_2 + n_20_0 + n_21_0 + n_22_0;
  long hgmr2_numer = n_00_2 + n_10_2 + n_20_2 + n_02_0 + n_12_0 + n_22_0;
  long hgmr1_denom = n_00_0 + n_01_0 + n_02_0 + n_20_2 + n_21_2 + n_22_2 + hgmr1_numer;
  long hgmr2_denom = n_00_0 + n_10_0 + n_20_0 + n_02_2 + n_12_2 + n_22_2 + hgmr2_numer;
  double hgmr1 = (hgmr1_denom > 0)? (double)hgmr1_numer/(double)hgmr1_denom : 2;
  double hgmr2 = (hgmr2_denom > 0)? (double)hgmr2_numer/(double)hgmr2_denom : 2;

  fprintf(stderr, "%4ld %4ld %4ld %7.5lf %7.5lf %7.5lf %7.5lf  ", n_0, n_1, n_2, d1, d2, hgmr1, hgmr2);
}

void triple_counts(char* id1, char* id2, char* gts1, char* gts2, char* proggts, double max_hgmr, double max_d1){
  char c1, c2, c3;
  long n_00_0 = 0, n_00_1 = 0, n_00_2 = 0;
  long n_01_0 = 0, n_01_1 = 0, n_01_2 = 0;
  long n_02_0 = 0, n_02_1 = 0, n_02_2 = 0;
  long n_03_0 = 0, n_03_1 = 0, n_03_2 = 0;
  
  long n_10_0 = 0, n_10_1 = 0, n_10_2 = 0;
  long n_11_0 = 0, n_11_1 = 0, n_11_2 = 0;
  long n_12_0 = 0, n_12_1 = 0, n_12_2 = 0;
  long n_13_0 = 0, n_13_1 = 0, n_13_2 = 0;
  
  long n_20_0 = 0, n_20_1 = 0, n_20_2 = 0;
  long n_21_0 = 0, n_21_1 = 0, n_21_2 = 0;
  long n_22_0 = 0, n_22_1 = 0, n_22_2 = 0;
  long n_23_0 = 0, n_23_1 = 0, n_23_2 = 0;

  long n_30_0 = 0, n_30_1 = 0, n_30_2 = 0;
  long n_31_0 = 0, n_31_1 = 0, n_31_2 = 0;
  long n_32_0 = 0, n_32_1 = 0, n_32_2 = 0;
  
  long i=0;
  while((c3 = proggts[i]) != '\0'){
    if(c3 != '3'){
      c1 = gts1[i];
      c2 = gts2[i];
      if(c1 == '0'){ // 0xy
	if(c2 == '0'){ // 00y
	  if(c3 == '0'){ 
	    n_00_0++;
	  }else if(c3 == '1'){
	    n_00_1++;
	  }else if(c3 == '2'){
	    n_00_2++;
	  }
	}else if(c2 == '1'){ // 01y
	  if(c3 == '0'){ 
	    n_01_0++;
	  }else if(c3 == '1'){
	    n_01_1++;
	  }else if(c3 == '2'){
	    n_01_2++;
	  }
	}else if(c2 == '2'){ // 02y
	  if(c3 == '0'){
	    n_02_0++;
	  }else if(c3 == '1'){
	    n_02_1++;
	  }else if(c3 == '2'){
	    n_02_2++;
	  }
	}
	else if(c2 == '3'){ // 03y
	  if(c3 == '0'){
	    n_03_0++;
	  }else if(c3 == '1'){
	    n_03_1++;
	  }else if(c3 == '2'){
	    n_03_2++;
	  }
	}
      }else if(c1 == '1'){ // 1xy
	if(c2 == '0'){ // 10y
	  if(c3 == '0'){
	    n_10_0++;
	  }else if(c3 == '1'){
	    n_10_1++;
	  }else if(c3 == '2'){
	    n_10_2++;
	  }
	}else if(c2 == '1'){ // 11y
	  if(c3 == '0'){
	    n_11_0++;
	  }else if(c3 == '1'){
	    n_11_1++;
	  }else if(c3 == '2'){
	    n_11_2++;
	  }
	}else if(c2 == '2'){ // 12y
	  if(c3 == '0'){
	    n_12_0++;
	  }else if(c3 == '1'){
	    n_12_1++;
	  }else if(c3 == '2'){
	    n_12_2++;
	  }
	}else if(c2 == '3'){ // 13y
	  if(c3 == '0'){
	    n_13_0++;
	  }else if(c3 == '1'){
	    n_13_1++;
	  }else if(c3 == '2'){
	    n_13_2++;
	  }
	}
      }else if(c1 == '2'){ // 2xy
	if(c2 == '0'){ // 20y
	  if(c3 == '0'){
	    n_20_0++;
	  }else if(c3 == '1'){
	    n_20_1++;
	  }else if(c3 == '2'){
	    n_20_2++;
	  }
	}else if(c2 == '1'){ // 21y
	  if(c3 == '0'){
	    n_21_0++;
	  }else if(c3 == '1'){
	    n_21_1++;
	  }else if(c3 == '2'){
	    n_21_2++;
	  }
	}else if(c2 == '2'){ // 22y
	  if(c3 == '0'){
	    n_22_0++;
	  }else if(c3 == '1'){
	    n_22_1++;
	  }else if(c3 == '2'){
	    n_22_2++;
	  }
	}else if(c2 == '3'){ // 23y
	  if(c3 == '0'){
	    n_23_0++;
	  }else if(c3 == '1'){
	    n_23_1++;
	  }else if(c3 == '2'){
	    n_23_2++;
	  }
	}
      }else if(c1 == '3'){ // 3xy
	if(c2 == '0'){ // 30y
	  if(c3 == '0'){
	    n_30_0++;
	  }else if(c3 == '1'){
	    n_30_1++;
	  }else if(c3 == '2'){
	    n_30_2++;
	  }
	}else if(c2 == '1'){ // 31y
	  if(c3 == '0'){
	    n_31_0++;
	  }else if(c3 == '1'){
	    n_31_1++;
	  }else if(c3 == '2'){
	    n_31_2++;
	  }
	}else if(c2 == '2'){ // 32y
	  if(c3 == '0'){
	    n_32_0++;
	  }else if(c3 == '1'){
	    n_32_1++;
	  }else if(c3 == '2'){
	    n_32_2++;
	  }
	}else{ // 33y
	  // do nothing
	}
      }
    }
    i++;
  }
  long n_0 =
    n_00_0 + n_01_0 + n_01_1 + n_02_1 +
    n_10_0 + n_10_1 + n_11_0 + n_11_1 + n_11_2 + n_12_1 + n_12_2 +
    n_20_1 + n_21_1 + n_21_2 + n_22_2; // can happen in no-error case
  long n_1 =
    n_00_1 + n_01_2 + n_02_0 + n_02_2 +
    n_10_2 + n_12_0 +
    n_20_0 + n_20_2 + n_21_0 + n_22_1;  // can happen if just one error of 0<->1 of 1<->2 type.
  long n_2 = n_00_2 + n_22_0; // can happen if one 0<->2 error, or two errors of 0<->1 or 1<->2 type.
  double d1 = (n_0 > 0)? (double)n_1/(double)(n_0 + n_1) : 2;
  double d2 = (n_0 > 0)? (double)n_2/(double)(n_0 + n_2) : 2;
  
  long n_0x_0 = n_00_0 + n_01_0 + n_02_0 + n_03_0;
  long n_0x_1 = n_00_1 + n_01_1 + n_02_1 + n_03_1;
  long n_0x_2 = n_00_2 + n_01_2 + n_02_2 + n_03_2;
  
  long n_1x_0 = n_10_0 + n_11_0 + n_12_0 + n_13_0;
  long n_1x_1 = n_10_1 + n_11_1 + n_12_1 + n_13_1;
  long n_1x_2 = n_10_2 + n_11_2 + n_12_2 + n_13_2;
  
  long n_2x_0 = n_20_0 + n_21_0 + n_22_0 + n_23_0;
  long n_2x_1 = n_20_1 + n_21_1 + n_22_1 + n_23_1;
  long n_2x_2 = n_20_2 + n_21_2 + n_22_2 + n_23_2;

  long n_x0_0 = n_00_0 + n_10_0 + n_20_0 + n_30_0;
  long n_x0_1 = n_00_1 + n_10_1 + n_20_1 + n_30_1;
  long n_x0_2 = n_00_2 + n_10_2 + n_20_2 + n_30_2;
  
  long n_x1_0 = n_01_0 + n_11_0 + n_21_0 + n_31_0;
  long n_x1_1 = n_01_1 + n_11_1 + n_21_1 + n_31_1;
  long n_x1_2 = n_01_2 + n_11_2 + n_21_2 + n_31_2;
  
  long n_x2_0 = n_02_0 + n_12_0 + n_22_0 + n_32_0;
  long n_x2_1 = n_02_1 + n_12_1 + n_22_1 + n_32_1;
  long n_x2_2 = n_02_2 + n_12_2 + n_22_2 + n_32_2;

  long hgmr1_numer = n_0x_2 + n_2x_0;
  long hgmr2_numer = n_x0_2 + n_x2_0;
  long hgmr1_denom = n_0x_0 + n_2x_2 + hgmr1_numer;
  long hgmr2_denom = n_x0_0 + n_x2_2 + hgmr2_numer;
  double hgmr1 = (hgmr1_denom > 0)? (double)hgmr1_numer/(double)hgmr1_denom : 2;
  double hgmr2 = (hgmr2_denom > 0)? (double)hgmr2_numer/(double)hgmr2_denom : 2;

  long r0x1_numer = n_0x_1 + n_2x_1;
  long r0x1_denom = r0x1_numer + n_0x_0 + n_2x_2;
  long rx01_numer = n_x0_1 + n_x2_1;
  long rx01_denom = rx01_numer + n_x0_0 + n_x2_2;
  double r0x1 = (r0x1_denom > 0)? (double)r0x1_numer/(double)r0x1_denom : 2;
  double rx01 = (rx01_denom > 0)? (double)rx01_numer/(double)rx01_denom : 2;
 
  /* long agmr1_numer = hgmr1_numer + n_0x_1 + n_2x_1 + n_1x_0 + n_1x_2; */
  /* long agmr2_numer = hgmr2_numer + n_x0_1 + n_x2_1 + n_x1_0 + n_x1_2; */
  /* long agmr1_denom = agmr1_numer + n_0x_0 + n_1x_1 + n_2x_2; */
  /* long agmr2_denom = agmr2_numer + n_x0_0 + n_x1_1 + n_x2_2; */
  /* double agmr1 = (agmr1_denom > 0)? (double)agmr1_numer/(double)agmr1_denom : 2; */
  /* double agmr2 = (agmr2_denom > 0)? (double)agmr2_numer/(double)agmr2_denom : 2; */
  /* fprintf(stderr, "%6.5lf %6.5lf %6.5lf  ", agmr1,  hgmr1, r0x1); */
  /* fprintf(stderr, "%6.5lf %6.5lf %6.5lf  ", agmr2,  hgmr2, rx01); */
  /* fprintf(stderr, "%6.5lf %6.5lf  ", d1, d2); */
  if(hgmr1 <= max_hgmr  &&  hgmr2 <= max_hgmr  &&  d1 <= max_d1){
    fprintf(stderr, "  %-30s %-30s  ", id1, id2);
    fprintf(stderr, "%6.5lf %6.5lf  ", hgmr1, r0x1);
    fprintf(stderr, "%6.5lf %6.5lf  ", hgmr2, rx01);
    fprintf(stderr, "%6.5lf %6.5lf  ", d1, d2);
    fprintf(stderr, "\n");
  }
}

Pedigree_stats* triple_counts_z(char* gts1, char* gts2, char* proggts){ // Pedigree* the_pedigree, GenotypesSet* the_gtsset){

  char c1, c2, c3;
  long n_00_0 = 0, n_00_1 = 0, n_00_2 = 0;
  long n_01_0 = 0, n_01_1 = 0, n_01_2 = 0;
  long n_02_0 = 0, n_02_1 = 0, n_02_2 = 0;
  long n_03_0 = 0, n_03_1 = 0, n_03_2 = 0;
  
  long n_10_0 = 0, n_10_1 = 0, n_10_2 = 0;
  long n_11_0 = 0, n_11_1 = 0, n_11_2 = 0;
  long n_12_0 = 0, n_12_1 = 0, n_12_2 = 0;
  long n_13_0 = 0, n_13_1 = 0, n_13_2 = 0;
  
  long n_20_0 = 0, n_20_1 = 0, n_20_2 = 0;
  long n_21_0 = 0, n_21_1 = 0, n_21_2 = 0;
  long n_22_0 = 0, n_22_1 = 0, n_22_2 = 0;
  long n_23_0 = 0, n_23_1 = 0, n_23_2 = 0;

  long n_30_0 = 0, n_30_1 = 0, n_30_2 = 0;
  long n_31_0 = 0, n_31_1 = 0, n_31_2 = 0;
  long n_32_0 = 0, n_32_1 = 0, n_32_2 = 0;
  
  long i=0;
  while((c3 = proggts[i]) != '\0'){
    if(c3 != '3'){
      c1 = gts1[i];
      c2 = gts2[i];
      if(c1 == '0'){ // 0xy
	if(c2 == '0'){ // 00y
	  if(c3 == '0'){ 
	    n_00_0++;
	  }else if(c3 == '1'){
	    n_00_1++;
	  }else if(c3 == '2'){
	    n_00_2++;
	  }
	}else if(c2 == '1'){ // 01y
	  if(c3 == '0'){ 
	    n_01_0++;
	  }else if(c3 == '1'){
	    n_01_1++;
	  }else if(c3 == '2'){
	    n_01_2++;
	  }
	}else if(c2 == '2'){ // 02y
	  if(c3 == '0'){
	    n_02_0++;
	  }else if(c3 == '1'){
	    n_02_1++;
	  }else if(c3 == '2'){
	    n_02_2++;
	  }
	}
	else if(c2 == '3'){ // 03y
	  if(c3 == '0'){
	    n_03_0++;
	  }else if(c3 == '1'){
	    n_03_1++;
	  }else if(c3 == '2'){
	    n_03_2++;
	  }
	}
      }else if(c1 == '1'){ // 1xy
	if(c2 == '0'){ // 10y
	  if(c3 == '0'){
	    n_10_0++;
	  }else if(c3 == '1'){
	    n_10_1++;
	  }else if(c3 == '2'){
	    n_10_2++;
	  }
	}else if(c2 == '1'){ // 11y
	  if(c3 == '0'){
	    n_11_0++;
	  }else if(c3 == '1'){
	    n_11_1++;
	  }else if(c3 == '2'){
	    n_11_2++;
	  }
	}else if(c2 == '2'){ // 12y
	  if(c3 == '0'){
	    n_12_0++;
	  }else if(c3 == '1'){
	    n_12_1++;
	  }else if(c3 == '2'){
	    n_12_2++;
	  }
	}else if(c2 == '3'){ // 13y
	  if(c3 == '0'){
	    n_13_0++;
	  }else if(c3 == '1'){
	    n_13_1++;
	  }else if(c3 == '2'){
	    n_13_2++;
	  }
	}
      }else if(c1 == '2'){ // 2xy
	if(c2 == '0'){ // 20y
	  if(c3 == '0'){
	    n_20_0++;
	  }else if(c3 == '1'){
	    n_20_1++;
	  }else if(c3 == '2'){
	    n_20_2++;
	  }
	}else if(c2 == '1'){ // 21y
	  if(c3 == '0'){
	    n_21_0++;
	  }else if(c3 == '1'){
	    n_21_1++;
	  }else if(c3 == '2'){
	    n_21_2++;
	  }
	}else if(c2 == '2'){ // 22y
	  if(c3 == '0'){
	    n_22_0++;
	  }else if(c3 == '1'){
	    n_22_1++;
	  }else if(c3 == '2'){
	    n_22_2++;
	  }
	}else if(c2 == '3'){ // 23y
	  if(c3 == '0'){
	    n_23_0++;
	  }else if(c3 == '1'){
	    n_23_1++;
	  }else if(c3 == '2'){
	    n_23_2++;
	  }
	}
      }else if(c1 == '3'){ // 3xy
	if(c2 == '0'){ // 30y
	  if(c3 == '0'){
	    n_30_0++;
	  }else if(c3 == '1'){
	    n_30_1++;
	  }else if(c3 == '2'){
	    n_30_2++;
	  }
	}else if(c2 == '1'){ // 31y
	  if(c3 == '0'){
	    n_31_0++;
	  }else if(c3 == '1'){
	    n_31_1++;
	  }else if(c3 == '2'){
	    n_31_2++;
	  }
	}else if(c2 == '2'){ // 32y
	  if(c3 == '0'){
	    n_32_0++;
	  }else if(c3 == '1'){
	    n_32_1++;
	  }else if(c3 == '2'){
	    n_32_2++;
	  }
	}else{ // 33y
	  // do nothing
	}
      }
    }
    i++;
  }
  long n_0 =
    n_00_0 + n_01_0 + n_01_1 + n_02_1 +
    n_10_0 + n_10_1 + n_11_0 + n_11_1 + n_11_2 + n_12_1 + n_12_2 +
    n_20_1 + n_21_1 + n_21_2 + n_22_2; // can happen in no-error case
  long n_1 =
    n_00_1 + n_01_2 + n_02_0 + n_02_2 +
    n_10_2 + n_12_0 +
    n_20_0 + n_20_2 + n_21_0 + n_22_1;  // can happen if just one error of 0<->1 of 1<->2 type.
  long n_2 = n_00_2 + n_22_0; // can happen if one 0<->2 error, or two errors of 0<->1 or 1<->2 type.
  double d1 = (n_0 > 0)? (double)n_1/(double)(n_0 + n_1) : 2;
  double d2 = (n_0 > 0)? (double)n_2/(double)(n_0 + n_2) : 2;
  
  long n_0x_0 = n_00_0 + n_01_0 + n_02_0 + n_03_0;
  long n_0x_1 = n_00_1 + n_01_1 + n_02_1 + n_03_1;
  long n_0x_2 = n_00_2 + n_01_2 + n_02_2 + n_03_2;
  
  long n_1x_0 = n_10_0 + n_11_0 + n_12_0 + n_13_0;
  long n_1x_1 = n_10_1 + n_11_1 + n_12_1 + n_13_1;
  long n_1x_2 = n_10_2 + n_11_2 + n_12_2 + n_13_2;
  
  long n_2x_0 = n_20_0 + n_21_0 + n_22_0 + n_23_0;
  long n_2x_1 = n_20_1 + n_21_1 + n_22_1 + n_23_1;
  long n_2x_2 = n_20_2 + n_21_2 + n_22_2 + n_23_2;

  long n_x0_0 = n_00_0 + n_10_0 + n_20_0 + n_30_0;
  long n_x0_1 = n_00_1 + n_10_1 + n_20_1 + n_30_1;
  long n_x0_2 = n_00_2 + n_10_2 + n_20_2 + n_30_2;
  
  long n_x1_0 = n_01_0 + n_11_0 + n_21_0 + n_31_0;
  long n_x1_1 = n_01_1 + n_11_1 + n_21_1 + n_31_1;
  long n_x1_2 = n_01_2 + n_11_2 + n_21_2 + n_31_2;
  
  long n_x2_0 = n_02_0 + n_12_0 + n_22_0 + n_32_0;
  long n_x2_1 = n_02_1 + n_12_1 + n_22_1 + n_32_1;
  long n_x2_2 = n_02_2 + n_12_2 + n_22_2 + n_32_2;

  long hgmr1_numer = n_0x_2 + n_2x_0;
  long hgmr2_numer = n_x0_2 + n_x2_0;
  long hgmr1_denom = n_0x_0 + n_2x_2 + hgmr1_numer;
  long hgmr2_denom = n_x0_0 + n_x2_2 + hgmr2_numer;
  // double hgmr1 = (hgmr1_denom > 0)? (double)hgmr1_numer/(double)hgmr1_denom : 2;
  // double hgmr2 = (hgmr2_denom > 0)? (double)hgmr2_numer/(double)hgmr2_denom : 2;

  long r0x1_numer = n_0x_1 + n_2x_1;
  long r0x1_denom = r0x1_numer + n_0x_0 + n_2x_2;
  long rx01_numer = n_x0_1 + n_x2_1;
  long rx01_denom = rx01_numer + n_x0_0 + n_x2_2;
  // double r0x1 = (r0x1_denom > 0)? (double)r0x1_numer/(double)r0x1_denom : 2;
  // double rx01 = (rx01_denom > 0)? (double)rx01_numer/(double)rx01_denom : 2;
 
  /* long agmr1_numer = hgmr1_numer + n_0x_1 + n_2x_1 + n_1x_0 + n_1x_2; */
  /* long agmr2_numer = hgmr2_numer + n_x0_1 + n_x2_1 + n_x1_0 + n_x1_2; */
  /* long agmr1_denom = agmr1_numer + n_0x_0 + n_1x_1 + n_2x_2; */
  /* long agmr2_denom = agmr2_numer + n_x0_0 + n_x1_1 + n_x2_2; */
  /* double agmr1 = (agmr1_denom > 0)? (double)agmr1_numer/(double)agmr1_denom : 2; */
  /* double agmr2 = (agmr2_denom > 0)? (double)agmr2_numer/(double)agmr2_denom : 2; */
  /* fprintf(stderr, "%6.5lf %6.5lf %6.5lf  ", agmr1,  hgmr1, r0x1); */
  /* fprintf(stderr, "%6.5lf %6.5lf %6.5lf  ", agmr2,  hgmr2, rx01); */
  /* fprintf(stderr, "%6.5lf %6.5lf  ", d1, d2); */

  Pedigree_stats* pedigree_stats = (Pedigree_stats*)malloc(sizeof(Pedigree_stats));
  ND hgmr1_nd = {hgmr1_numer, hgmr1_denom};
  pedigree_stats->par1_hgmr = hgmr1_nd;
  ND r1_nd = {r0x1_numer, r0x1_denom};
  pedigree_stats->par1_r = r1_nd;
  ND hgmr2_nd = {hgmr2_numer, hgmr2_denom};
  pedigree_stats->par2_hgmr = hgmr2_nd; 
  ND r2_nd = {rx01_numer, rx01_denom};
  pedigree_stats->par2_r = r2_nd;
  ND d1_nd = {n_1, n_0 + n_1};
  pedigree_stats->d1 = d1_nd;
  ND d2_nd = {n_2, n_0 + n_2};
  pedigree_stats->d2 = d2_nd;
  
  /* if(hgmr1 <= max_hgmr  &&  hgmr2 <= max_hgmr  &&  d1 <= max_d1){ */
  /*   fprintf(stderr, "  %-30s %-30s  ", id1, id2); */
  /*   fprintf(stderr, "%6.5lf %6.5lf  ", hgmr1, r0x1); */
  /*   fprintf(stderr, "%6.5lf %6.5lf  ", hgmr2, rx01); */
  /*   fprintf(stderr, "%6.5lf %6.5lf  ", d1, d2); */
  /*   fprintf(stderr, "\n"); */
  /* } */

  return pedigree_stats;
}


/* void calculate_pedigree_test_info(Pedigree* the_pedigree, GenotypesSet* the_gtsset){ */

/*   char* fempar_gts = the_gtsset->genotype_sets->a[the_pedigree->Fparent->index]; */
/*   char* malpar_gts = the_gtsset->genotype_sets->a[the_pedigree->Mparent->index]; */
/*   char* acc_gts = the_gtsset->genotype_sets->a[the_pedigree->Accession->index]; */
    
/*   Ahr fp_ahr, mp_ahr, fm_ahr; */
/*   agmr_hgmr_r(fempar_gts, acc_gts, &(the_pedigree->fp_ahr)); */
/*   agmr_hgmr_r(malpar_gts, acc_gts, &(the_pedigree->mp_ahr)); */
/*   agmr_hgmr_r(fempar_gts, malpar_gts, &(the_pedigree->fm_ahr)); */
/* } */

Pedigree_stats* calculate_pedigree_stats(Pedigree* the_pedigree, GenotypesSet* the_gtsset){

  char* fempar_gts = the_gtsset->genotype_sets->a[the_pedigree->Fparent->index];
  char* malpar_gts = the_gtsset->genotype_sets->a[the_pedigree->Mparent->index];
  char* acc_gts = the_gtsset->genotype_sets->a[the_pedigree->Accession->index];

  Pedigree_stats* the_pedigree_stats = triple_counts_z(fempar_gts, malpar_gts, acc_gts);
  return the_pedigree_stats;
}


Vlong* accessions_with_offspring(Vpedigree* the_vped, long n_accessions){
  Vlong* offspring_counts = construct_vlong_zeroes(n_accessions);
  for(long i=0; i<the_vped->size; i++){
    Pedigree* the_ped = the_vped->a[i];
    long fp_index = the_ped->Fparent->index;
    offspring_counts->a[fp_index]++;
    long mp_index = the_ped->Mparent->index;
    offspring_counts->a[mp_index]++; 
  }
  Vlong* accids_with_offspring = construct_vlong(100);
  for(long i=0; i<offspring_counts->size; i++){
    if(offspring_counts->a[i] > 0){
      add_long_to_vlong(accids_with_offspring, i);
    }
  }
  free_vlong(offspring_counts);
  return accids_with_offspring;
}

void print_pedigree_test_info(FILE* fh, Pedigree* the_pedigree, GenotypesSet* the_gtsset, Vlong* parent_idxs){

  char* acc_id = the_pedigree->Accession->id;
  long acc_idx = the_pedigree->Accession->index;
  long fparent_idx = the_pedigree->Fparent->index;
  long mparent_idx = the_pedigree->Mparent->index;
  Vlong* best_parent_candidates = construct_vlong(10);
  add_long_to_vlong(best_parent_candidates, fparent_idx);
  add_long_to_vlong(best_parent_candidates, mparent_idx);

  fprintf(fh, "%-20s %4ld %-20s %-20s  ",
	  the_gtsset->accession_ids->a[acc_idx],
	  the_gtsset->accession_missing_data_counts->a[the_pedigree->Accession->index],
	  the_gtsset->accession_ids->a[fparent_idx],
	  the_gtsset->accession_ids->a[mparent_idx]
	  );

  print_Ahrx(fh, the_pedigree->fp_ahr); fprintf(fh, "   ");
  print_Ahrx(fh, the_pedigree->mp_ahr); fprintf(fh, "   ");
  //  print_Ahr(fh, the_pedigree->fm_ahr); fprintf(fh, "   ");
}

void print_pedigree_stats(FILE* fh, Pedigree_stats* the_pedigree_stats){
  fprintf(fh, "%6.5lf ", (the_pedigree_stats->par1_hgmr.d > 0)? (double)the_pedigree_stats->par1_hgmr.n/(double)the_pedigree_stats->par1_hgmr.d : 2);
  fprintf(fh, "%6.5lf ", (the_pedigree_stats->par1_r.d > 0)? (double)the_pedigree_stats->par1_r.n/(double)the_pedigree_stats->par1_r.d : 2);
  fprintf(fh, "%6.5lf ", (the_pedigree_stats->par2_hgmr.d > 0)? (double)the_pedigree_stats->par2_hgmr.n/(double)the_pedigree_stats->par2_hgmr.d : 2);
  fprintf(fh, "%6.5lf ", (the_pedigree_stats->par2_r.d > 0)? (double)the_pedigree_stats->par2_r.n/(double)the_pedigree_stats->par2_r.d : 2);
  fprintf(fh, "%6.5lf ", (the_pedigree_stats->d1.d > 0)? (double)the_pedigree_stats->d1.n/(double)the_pedigree_stats->d1.d : 2);
  fprintf(fh, "%6.5lf ", (the_pedigree_stats->d2.d > 0)? (double)the_pedigree_stats->d2.n/(double)the_pedigree_stats->d2.d : 2);
}

double get_hgmr1(Pedigree_stats* p){
  return (p->par1_hgmr.d > 0)? (double)p->par1_hgmr.n/(double)p->par1_hgmr.d : 2;
} 
double get_r1(Pedigree_stats* p){
   return (p->par1_r.d > 0)? (double)p->par1_r.n/(double)p->par1_r.d : 2;
}
double get_hgmr2(Pedigree_stats* p){
  return (p->par2_hgmr.d > 0)? (double)p->par2_hgmr.n/(double)p->par2_hgmr.d : 2;
} 
double get_r2(Pedigree_stats* p){
   return (p->par2_r.d > 0)? (double)p->par2_r.n/(double)p->par2_r.d : 2;
} 
double get_d1(Pedigree_stats* p){
 return (p->d1.d > 0)? (double)p->d1.n/(double)p->d1.d : 2;
}
double get_d2(Pedigree_stats* p){
 return (p->d2.d > 0)? (double)p->d2.n/(double)p->d2.d : 2;
}

void print_pedigree_alternatives(FILE* fh, Pedigree* the_pedigree, GenotypesSet* the_gtsset, Vlong* parent_idxs){
    long n_parents = parent_idxs->size;

     char* acc_id = the_pedigree->Accession->id;
  long acc_idx = the_pedigree->Accession->index;
  long fparent_idx = the_pedigree->Fparent->index;
  long mparent_idx = the_pedigree->Mparent->index;
  Vlong* best_parent_candidates = construct_vlong(10);
  add_long_to_vlong(best_parent_candidates, fparent_idx);
  if(mparent_idx != fparent_idx) add_long_to_vlong(best_parent_candidates, mparent_idx);
  Idxhgmr* the_idxhgmrs = (Idxhgmr*)malloc(n_parents*sizeof(Idxhgmr));
  char* acc_gts = the_gtsset->genotype_sets->a[the_pedigree->Accession->index];
  for(long i=0; i<parent_idxs->size; i++){
    long idx = parent_idxs->a[i];
    char* pgts = the_gtsset->genotype_sets->a[idx];
    //   agmr_hgmr_r(gts1, gts2, &the_ahr);
    //  double the_hgmr = hgmr(gts1, gts2); //(the_ahr.h.d > 0)? (double)the_ahr.h.n/(double)the_ahr.h.d : -1.0;
    the_idxhgmrs[i].idx = idx;
    the_idxhgmrs[i].hgmr = hgmr(acc_gts, pgts); //the_hgmr; //the_idxhgmr = {idx, the_hgmr);
  }
  //  fprintf(stderr, "just before sort idxhgmr\n");
  sort_idxhgmr_by_hgmr(n_parents, the_idxhgmrs);
  //  fprintf(stderr, "just after sort idxhgmr\n");
  long n_printed = 0;
  long hgmr_ge0_count = 0;
  long fpar_rank = -1;
  long mpar_rank = -1;
  //char alternative_candidates[4000];
  //  Vstr* alt_candidates = construct_vstr(100);
  for(long i=0; i<n_parents; i++){
    long the_idx = the_idxhgmrs[i].idx;
    if(the_idx != acc_idx){
      double the_hgmr = the_idxhgmrs[i].hgmr;
      if(the_hgmr >= 0.06) break; 
      if(the_hgmr >= 0){
	hgmr_ge0_count++;
	long is_a_parent = 0;
	if(the_idx == fparent_idx){
	  fpar_rank = hgmr_ge0_count;
	  is_a_parent = 1;
	}
	if(the_idx == mparent_idx){
	  mpar_rank = hgmr_ge0_count;
	  is_a_parent = 1;
	}
	if(is_a_parent == 0){
	  char* alt_cand = (char*)malloc(100*sizeof(char));
	  /* sprintf(alt_cand, "%s %7.5lf  ", the_gtsset->accession_ids->a[the_idx], the_hgmr); */
	  /* add_string_to_vstr(alt_candidates, alt_cand); */
	  add_long_to_vlong(best_parent_candidates, the_idx);
	  n_printed++;
	}
      }
    }
  }
  //  fprintf(fh, "%ld %ld %ld   ", n_printed, long_min(fpar_rank, mpar_rank), long_max(fpar_rank, mpar_rank));
  long ub = long_min(best_parent_candidates->size, 8);

  char* fparent_gts = the_gtsset->genotype_sets->a[fparent_idx];
  char* mparent_gts = the_gtsset->genotype_sets->a[mparent_idx];
  //  fprintf(stderr, "# %-30s \n", acc_id);
  /* fprintf(stderr, "   %-30s %-30s  %7.5lf %7.5lf \n", */
  /* 	  the_gtsset->accession_ids->a[fparent_idx], the_gtsset->accession_ids->a[mparent_idx], */
  /* 	  hgmr(fparent_gts, acc_gts), hgmr(mparent_gts, acc_gts)); */
  //  fprintf(stderr, "   %-30s %-30s  ", the_gtsset->accession_ids->a[fparent_idx], the_gtsset->accession_ids->a[mparent_idx]);
  /* triple_counts(the_gtsset->accession_ids->a[fparent_idx], the_gtsset->accession_ids->a[mparent_idx], */
  /* 		  fparent_gts, mparent_gts, acc_gts, 1, 1); // fprintf(stderr, "\n"); */
  for(long i=0; i<ub; i++){
    long idx1 = best_parent_candidates->a[i];
    char* id1 = the_gtsset->accession_ids->a[idx1];
    char* gts1 = the_gtsset->genotype_sets->a[idx1];
    for(long j=i; j<ub; j++){
      long idx2 = best_parent_candidates->a[j];
      char* id2 = the_gtsset->accession_ids->a[idx2];
      char* gts2 = the_gtsset->genotype_sets->a[idx2];
      //  fprintf(stderr, "   %-30s %-30s  ", id1, id2);
      //  triple_counts(id1, id2, gts1, gts2, acc_gts, 0.05, 0.015);
      Pedigree_stats* alt_pedigree_stats = triple_counts_z(gts1, gts2, acc_gts);
    
      if(get_hgmr1(alt_pedigree_stats) <= 0.05  && get_hgmr2(alt_pedigree_stats) <= 0.05  &&  get_d1(alt_pedigree_stats) <= 0.015){
	fprintf(fh, "%s %s ", id1, id2);
	print_pedigree_stats(fh, alt_pedigree_stats);
      }
      //  fprintf(stderr, "\n");
    }
  }// fprintf(fh, "\n");
  
  // fprintf(fh, "%s\n", alternative_candidates);
  // print_vstr(fh, alt_candidates);
  // free_vstr(alt_candidates);
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
  // fprintf(stderr "after reading first (header) line of pedigrees file\n");
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
  fprintf(stderr, "# size of Vpedigree pedigrees: %ld \n", pedigrees->size);
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

// *****  sorting an array of Idxhgmr  *****
int cmpidxhgmr(const void* v1, const void* v2){
  const Idxhgmr* s1 = (const Idxhgmr*)v1;
  const Idxhgmr* s2 = (const Idxhgmr*)v2;
  return (s1->hgmr < s2->hgmr)? -1 : (s1->hgmr > s2->hgmr)? 1 : 0;
}

void sort_idxhgmr_by_hgmr(long size, Idxhgmr* array){ // sort in place
  qsort(array, size, sizeof(Idxhgmr), cmpidxhgmr);
}

long long_min(long a, long b){
  return (a <= b)? a : b;
}

long long_max(long a, long b){
  return (a >= b)? a : b;
}
