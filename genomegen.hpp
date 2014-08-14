#ifndef __GENOME_GEN__
#define __GENOME_GEN__

#include <string>

using namespace std;

struct cmd_args{
	string genome_id;
	int num_chroms;
	unsigned long chrom_size;
	char scale;
};

cmd_args parse_args(int argc, char** argv);
string generate_ref_genome(string id, int num_chroms, unsigned long chrom_size);

#endif