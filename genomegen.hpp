#ifndef __GENOME_GEN__
#define __GENOME_GEN__

#include <string>
#include <vector>


namespace genome_gen{

using namespace std;

struct cmd_args{
	string genome_id;
	int num_chroms;
	unsigned long chrom_size;
	char scale;
};

cmd_args parse_args(int argc, char** argv);
vector<char> create_genome_container(int num_chroms, unsigned long chrom_size);
void write_ref_genome(
	string id, int num_chroms, unsigned long chrom_size, vector<char>& genome);
void write_private_genome(string id, vector<char>& genome);
void write_reads(string id, vector<char>& genome);
void generate_copies(vector<char>& genome);
void generate_inversions(vector<char>& genome);
void generate_insertions(vector<char>& genome);
void generate_deletions(vector<char>& genome);
void generate_snps(vector<char>& genome);
void generate_alus(vector<char>& genome);
void generate_strs(vector<char>& genome);
char random_snp(char base);
double rand_num();

};

#endif