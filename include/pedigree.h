// *****  typedefs for Pedigree, Vpedigree *****

typedef struct{
  IndexId* Fparent;
  IndexId* Mparent;
  IndexId* Accession;
}Pedigree;

typedef struct{
  long capacity;
  long size;
  Pedigree** a;
}Vpedigree;

// *****  functions  *****

// *****  Pedigree  *****
Pedigree* construct_pedigree(IndexId* acc_idxid, IndexId* fempar_idxid, IndexId* malpar_idxid);
void agmr_hgmr_r(char* gts1, char* gts2, ND3* the_nd3);
double hgmr(char* gts1, char* gts2);
void print_pedigree_stats(Pedigree* the_pedigree, GenotypesSet* the_gtsset);
void free_pedigree(Pedigree* the_pedigree);

// *****  Vpedigree  *****
Vpedigree* read_the_pedigrees_file_and_store(FILE* p_stream, Vidxid* the_vidxid); 
Vpedigree* construct_vpedigree(long cap);
void add_pedigree_to_vpedigree(Vpedigree* the_vped, Pedigree* the_ped);
void free_vpedigree(Vpedigree* the_vped);
