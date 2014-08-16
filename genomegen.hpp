#ifndef __GENOME_GEN__
#define __GENOME_GEN__

#include <string>
// #include <vector>


namespace genome_gen{

using namespace std;

void parse_args(int argc, char** argv);

char** create_genome_container();
void cleanup_genome_container();

void write_ref_genome(char** genome);
void write_private_genome(char** genome);
void write_reads(char** genome);

void generate_copies(char** genome);
void generate_inversions(char** genome);
void generate_insertions(char** genome);
void generate_deletions(char** genome);
void generate_snps(char** genome);
void generate_alus(char** genome);
void generate_strs(char** genome);

char random_snp(char base);
double rand_num();

};

#endif