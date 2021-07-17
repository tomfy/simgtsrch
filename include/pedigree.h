// *****  typedefs for Pedigree, Vpedigree *****
typedef struct{
  long n; // numerator
  long d; // denominator
}ND; // numerator and denominator

typedef struct{
  ND a;
  ND h;
  ND r;
}Ahr;

typedef struct{
  ND par1_hgmr;
  ND par1_r;
  ND par2_hgmr;
  ND par2_r;
  ND d1;
  ND d2;
}Pedigree_stats; // 

typedef struct{
  /* IndexId* Fparent; */
  /* IndexId* Mparent; */
  /* IndexId* Accession; */
  Accession* F;
  Accession* M;
  Accession* A;
  /* Ahr fp_ahr; // female parent - progeny */
  /* Ahr mp_ahr; // male parent - progeny */
  /* Ahr fm_ahr; // female - male parents */
  /* ND d1;  */
  /* ND d2; */
}Pedigree;

typedef struct{
  long capacity;
  long size;
  Pedigree** a;
}Vpedigree;

typedef struct{
  long idx;
  double hgmr;
}Idxhgmr;

// *****  function declarations  *****

// *****  Ahr  ******
void print_Ahr(FILE* fh, Ahr the_ahr);

// *****  Pedigree  *****
Pedigree* construct_pedigree(Accession* Acc, Accession* Fparent, Accession* Mparent);
Pedigree* construct_pedigree_from_idxids(IndexId* acc_idxid, IndexId* fempar_idxid, IndexId* malpar_idxid);
void agmr_hgmr_r(char* gts1, char* gts2, Ahr* the_nd3);
double hgmr(char* gts1, char* gts2);
// void calculate_pedigree_test_info(Pedigree* the_pedigree); // , GenotypesSet* the_gtsset);
Pedigree_stats* calculate_pedigree_stats(Pedigree* the_pedigree); // , GenotypesSet* the_gtsset);
// void print_pedigree_test_info(FILE* fh, Pedigree* the_pedigree, GenotypesSet* the_gtsset, Vlong* parents_idxs);
void print_pedigree_alternatives(FILE* fh, Pedigree* the_pedigree, const GenotypesSet* const the_gtsset, Vlong* parent_idxs);
void free_pedigree(Pedigree* the_pedigree);

// *****  Vpedigree  *****
Vpedigree* read_the_pedigrees_file_and_store(FILE* p_stream, Vidxid* the_vidxid, GenotypesSet* the_gtsset); 
Vpedigree* construct_vpedigree(long cap);
Vlong* accessions_with_offspring(Vpedigree* the_Vped, long n_accessions);
void add_pedigree_to_vpedigree(Vpedigree* the_vped, const Pedigree* const the_ped);
void free_vpedigree(Vpedigree* the_vped);

// *****  array of Idxhgmr  *****
int cmpidxhgmr(const void* v1, const void* v2);
void sort_idxhgmr_by_hgmr(long size, Idxhgmr* array);

// *****  miscellaneous  *****
long long_min(long a, long b);
long long_max(long a, long b);

void triple_counts_x(char* gts1, char* gts2, char* proggts);
void triple_counts(char* id1, char* id2, char* gts1, char* gts2, char* proggts, double max_hgmr, double max_d1);
Pedigree_stats* triple_counts_z(char* gts1, char* gts2, char* proggts);
// Pedigree_stats* calculate_pedigree_stats(Pedigree* the_pedigree, GenotypesSet* the_gtsset);
void print_pedigree_stats(FILE* fh, Pedigree_stats* the_pedigree_stats);
double get_hgmr1(Pedigree_stats* p);
double get_r1(Pedigree_stats* p);
double get_hgmr2(Pedigree_stats* p);
double get_r2(Pedigree_stats* p);
double get_d1(Pedigree_stats* p);
double get_d2(Pedigree_stats* p);
