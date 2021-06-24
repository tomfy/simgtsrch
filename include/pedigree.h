// *****  typedefs for Pedigree, Vpedigree *****
typedef struct{
  long n;
  long d;
}ND; // numerator and denominator

typedef struct{
  ND a;
  ND h;
  ND r;
}Ahr;   

typedef struct{
  IndexId* Fparent;
  IndexId* Mparent;
  IndexId* Accession;
  Ahr fp_ahr; // female parent - progeny
  Ahr mp_ahr; // male parent - progeny
  Ahr fm_ahr; // female - male parents
}Pedigree;

typedef struct{
  long capacity;
  long size;
  Pedigree** a;
}Vpedigree;

// *****  function declarations  *****

// *****  Ahr  ******
void print_Ahr(FILE* fh, Ahr the_ahr);

// *****  Pedigree  *****
Pedigree* construct_pedigree(IndexId* acc_idxid, IndexId* fempar_idxid, IndexId* malpar_idxid);
void agmr_hgmr_r(char* gts1, char* gts2, Ahr* the_nd3);
double hgmr(char* gts1, char* gts2);
void calculate_pedigree_test_info(Pedigree* the_pedigree, GenotypesSet* the_gtsset);
void print_pedigree_test_info(FILE* fh, Pedigree* the_pedigree, GenotypesSet* the_gtsset);
void free_pedigree(Pedigree* the_pedigree);

// *****  Vpedigree  *****
Vpedigree* read_the_pedigrees_file_and_store(FILE* p_stream, Vidxid* the_vidxid); 
Vpedigree* construct_vpedigree(long cap);
void add_pedigree_to_vpedigree(Vpedigree* the_vped, Pedigree* the_ped);
void free_vpedigree(Vpedigree* the_vped);
