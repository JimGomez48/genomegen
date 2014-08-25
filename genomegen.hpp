#pragma once
#ifndef __GENOME_GEN__
#define __GENOME_GEN__

#include <string>
#include <vector>


namespace genome_gen {

using namespace std;

void parse_args(int argc, char** argv);
void print_args();

void write_ref_genome(vector<vector<char>>& genome);
void write_private_genome(vector<vector<char>>& genome);
void write_reads(char** genome);

void generate_mutations(vector<vector<char>>& genome);
void generate_deletes(vector<vector<char>>& genome);
void generate_inserts(vector<vector<char>>& genome);
void generate_snps(vector<vector<char>>& genome);
// void generate_copies(vector<vector<char>>& genome);
// void generate_inversions(vector<vector<char>>& genome);
// void generate_alus(vector<vector<char>>& genome);
// void generate_strs(vector<vector<char>>& genome);

string get_garbage_read();
string get_slice(
	int chromosome, 
	unsigned long start, 
	unsigned long end, 
	const vector<vector<char>>& genome
);

string random_alleles(unsigned int length);
char random_allele();
char random_snp(char base);
double random();

};

#endif