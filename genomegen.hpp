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
string generate_ref_genome(string id, int num_chroms, unsigned long chrom_size);
void generate_copies(vector<char>* genome);
void generate_inversions(vector<char>* genome);
void generate_insertions(vector<char>* genome);
void generate_deletions(vector<char>* genome);
void generate_snps(vector<char>* genome);
void generate_alus(vector<char>* genome);
void generate_strs(vector<char>* genome);

};

#endif