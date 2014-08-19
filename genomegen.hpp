#ifndef __GENOME_GEN__
#define __GENOME_GEN__

#include <string>
#include <vector>


namespace genome_gen{

using namespace std;

void parse_args(int argc, char** argv);

void write_ref_genome(vector<vector<char>>& genome);
void write_private_genome(vector<vector<char>>& genome);
void write_reads(char** genome);

void generate_insertions(vector<vector<char>>& genome);
void generate_deletions(vector<vector<char>>& genome);
void generate_snps(vector<vector<char>>& genome);
// void generate_copies(char** genome);
// void generate_inversions(char** genome);
// void generate_alus(char** genome);
// void generate_strs(char** genome);

string get_slice(int chromosome, 
	unsigned long start, unsigned long end, vector<vector<char>>& genome);
char random_snp(char base);
double random();

};

#endif