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
void write_donor_genome(vector<vector<char>>& genome);
void write_reads(char** genome);

void generate_mutations(vector<vector<char>>& genome);
void generate_indels(vector<vector<char>>& genome);
void generate_snps(vector<vector<char>>& genome);

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
bool overlaps(const unsigned long start1, const unsigned int length1,
	const unsigned long start2, const unsigned int length2);

};

#endif