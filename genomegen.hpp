#ifndef __GENOME_GEN__
#define __GENOME_GEN__

#include <string>
#include <vector>


namespace genome_gen{

using namespace std;

void parse_args(int argc, char** argv);
void print_args();

void write_ref_genome(vector<vector<char>>& genome);
void write_private_genome(vector<vector<char>>& genome);
void write_reads(char** genome);

void generate_mutations(vector<vector<char>>& genome);
void generate_insertions(vector<char>& chromosome);
void generate_deletions(vector<char>& chromosome);
void generate_snps(vector<char>& chromosome);
// void generate_copies(vector<char>& chromosome);
// void generate_inversions(vector<char>& chromosome);
// void generate_alus(vector<char>& chromosome);
// void generate_strs(vector<char>& chromosome);

string get_slice(int chromosome, 
	unsigned long start, unsigned long end, const vector<char>& genome);
char random_snp(char base);
double random();

};

#endif