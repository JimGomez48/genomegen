#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <tuple>
// #include <utility>
// #include <random>

#include "genomegen.hpp"
#include "tclap/CmdLine.h"

using namespace std;
using namespace genome_gen;

namespace genome_gen{

static const string REF_PRE = "ref_";
static const string PRIV_PRE = "private_";
static const string READS_PRE = "reads_";
static const string ANS_PRE = "ans_";
static const string SUFFIX = ".txt";
static const unsigned int LINE_WIDTH = 80;

struct cmd_args{
	bool quiet;
	string genome_id;
	unsigned int num_chroms;
	unsigned long chrom_size;
	unsigned long genome_size;
	char scale;
	double coverage;
	unsigned int read_length;
	unsigned int read_gap;
	double read_error_rate;
	double read_garbage_rate;
	double snp_rate;
	double indel_rate;
} args;

string ref_file_name;
string priv_file_name;
string reads_file_name;
string ans_file_name;


void parse_args(int argc, char** argv){
	try{
		// Define Command line parser and description
		TCLAP::CmdLine cmd(
			"DESCRIPTION\nThis program generates a reference genome, a mutated "
			"donor genome, and a set of reads from the donor genome. The "
			"generated genomic data is written to the following set of files: "
			"\n1) reference genome \'ref_<genome_id>.txt\'\n2) mutated donor "
			"genome \'private_<genome_id>.txt\'\n3) paired-end reads "
            "\'reads_<genome_id>.txt\' from donor genome\n4) mutation answer "
			"key \'ans_<genome_id>.txt\'", 
			' '
		);
		// define args
		TCLAP::SwitchArg quiet(
			"q",
			"quiet",
			"If set, no ouput will be printed to stdout", 
			false
		);
		TCLAP::UnlabeledValueArg<string> genome_id(
			"genome-id", 
			"Genome-ID: The name or ID of this genome. (e.g. 'genome1')", 
			true, 
			"genome1",
			"string"
		);
		TCLAP::UnlabeledValueArg<unsigned int> num_chroms(
			"num-chromosomes",
			"Num-Chromosomes: The number of chromosomes to generate for the "
			"genome",
			true,
			1,
			"int"
		);
		TCLAP::UnlabeledValueArg<double> chrom_size(
			"chromosome-size",
			"Chromosome-Size: The size of each chromosome, by default in "
			"thousands of bp. (Change scale with -s option)",
			true,
			1,
			"float"
		);
		vector<char> allowed;
		allowed.push_back('k');
		allowed.push_back('m');
		allowed.push_back('b');
		TCLAP::ValuesConstraint<char> allowed_vals(allowed);
		TCLAP::ValueArg<char> scale(
			"s", "scale",
			"The amount to scale chromosome-size by. 'k': 1-thousand, 'm': "
			"1-million, 'b': 1-billion. Default: 'k' (1-thousand).",
			false,
			'k',
			// "k, m, b"
			&allowed_vals
		);
		TCLAP::ValueArg<double> coverage(
			"c", "coverage",
			"The amount of genome coverage to reflect in the reads. "
			"Default: 30.0 (30x).",
			false,
			30.0,
			"float"
		);
		TCLAP::ValueArg<unsigned int> read_length(
			"", "read-length",
			"The length of each read in bp. Default: 50",
			false,
			50,
			"int"
		);
		TCLAP::ValueArg<unsigned int> read_gap(
			"", "read-gap",
			"The average gap between paired end reads in bp. Actual gaps will "
			" vary by [-10, 10] bp. Default: 100",
			false,
			100,
			"int"
		);
		TCLAP::ValueArg<double> read_error_rate(
			"", "read-error",
			"The rate at which single base errors will be introduced into any "
			"given read. Default: 0.00",
			false,
			0.0,
			"float"
		);
		TCLAP::ValueArg<double> read_garbage_rate(
			"", "read-garbage",
			"The rate at which any given read is entirely erroneous (i.e. "
			"a garbage read). Default: 0.0",
			false,
			0.0,
			"float"
		);
		TCLAP::ValueArg<double> snp_rate(
			"", "snp",
			"The rate at which SNPs will be introduced into the mutated genome."
			" Default: 0.003",
			false,
			0.003,
			"float"
		);
		TCLAP::ValueArg<double> indel_rate(
			"", "indel",
			"The rate at which indels will be introduced into the mutated "
			"genome. Default: 0.001",
			false,
			0.001,
			"float"
		);

		// bind args to parser and parse
		cmd.add(quiet);
		cmd.add(genome_id);
		cmd.add(num_chroms);
		cmd.add(chrom_size);
		cmd.add(scale);
		cmd.add(coverage);
		cmd.add(read_length);
		cmd.add(read_gap);
		cmd.add(read_error_rate);
		cmd.add(read_garbage_rate);
		cmd.add(snp_rate);
		cmd.add(indel_rate);
		cmd.parse(argc, argv);
		if (coverage.getValue() < 0)
			throw TCLAP::ArgParseException("value cannot be < 0", "coverage");
		if (read_error_rate.getValue() < 0 || read_error_rate.getValue() > 1)
			throw TCLAP::ArgParseException("out of valid range: [0, 1]", "read-error");
		if (read_garbage_rate.getValue() < 0 || read_garbage_rate.getValue() > 1)
			throw TCLAP::ArgParseException("out of valid range: [0, 1]", "read-garbage");
		if (snp_rate.getValue() < 0 || snp_rate.getValue() > 1)
			throw TCLAP::ArgParseException("out of valid range: [0, 1]", "snp");
		if (indel_rate.getValue() < 0 || indel_rate.getValue() > 1)
			throw TCLAP::ArgParseException("out of valid range: [0, 1]", "indel");

    	// determine scale factor
    	double factor;
    	switch(scale.getValue()){
    		case 'k': factor = 1000; break;
    		case 'm': factor = 1000000; break;
    		case 'b': factor = 1000000000; break;
    	}

    	// save arg values in global struct
		args.quiet 				= quiet.getValue();
		args.genome_id 			= genome_id.getValue();
		args.num_chroms 		= num_chroms.getValue();
		args.chrom_size 		= factor * chrom_size.getValue();
		args.scale 				= scale.getValue();
		args.coverage 			= coverage.getValue();
		args.read_length 		= read_length.getValue();
		args.read_gap 			= read_gap.getValue();
		args.read_error_rate 	= read_error_rate.getValue();
		args.read_garbage_rate 	= read_garbage_rate.getValue();
    	args.genome_size 		= args.chrom_size * args.num_chroms;
    	args.snp_rate 			= snp_rate.getValue();
    	args.indel_rate 		= indel_rate.getValue();

    	// set file names
    	ref_file_name = REF_PRE + args.genome_id + SUFFIX;
    	priv_file_name = PRIV_PRE + args.genome_id + SUFFIX;
    	reads_file_name = READS_PRE + args.genome_id + SUFFIX;
    	ans_file_name = ANS_PRE + args.genome_id + SUFFIX;
    	
    	if (!args.quiet)
    		print_args();
	}
	catch(TCLAP::ArgException &e){
		if (!args.quiet)
			cout<<"PARSE ERROR:\n\t"<<e.what()<<endl;
		exit(1);
	}
}

void print_args(){
	cout<<"================================="<<endl;
	cout<<"genome-ID:\t"<<args.genome_id<<endl;
	cout<<"num-chroms:\t"<<args.num_chroms<<endl;
	switch(args.scale){
		case 'k': 
			cout<<"chrom-size:\t"<<args.chrom_size / 1000.0<<"-thousand bp"
			<<endl; break;
		case 'm': 
			cout<<"chrom-size:\t"<<args.chrom_size / 1000000.0<<"-million bp"
			<<endl; break;
		case 'b': 
			cout<<"chrom-size:\t"<<args.chrom_size / 1000000000.0<<"-billion bp"
			<<endl; break;
	}
	switch(args.scale){
		case 'k': 
			cout<<"genome-size:\t"<<args.genome_size / 1000.0<<"-thousand bp"
			<<endl; break;
		case 'm': 
			cout<<"genome-size:\t"<<args.genome_size / 1000000.0<<"-million bp"
			<<endl; break;
		case 'b': 
			cout<<"genome-size:\t"<<args.genome_size / 1000000000.0<<
			"-billion bp"<<endl; break;
	}
	cout<<"coverage:\t"<<args.coverage<<'x'<<endl;
	cout<<"read-length:\t"<<args.read_length<<" bp"<<endl;
	cout<<"read-gap:\t"<<args.read_gap<<" bp"<<endl;
	cout<<"read-error:\t"<<args.read_error_rate<<endl;
	cout<<"garbage-rate:\t"<<args.read_garbage_rate<<endl;
	cout<<"snp-rate:\t"<<args.snp_rate<<endl;
	cout<<"indel-rate:\t"<<args.indel_rate<<endl;
	cout<<"================================="<<endl;
}

void write_ref_genome(vector<vector<char>>& genome){
	ofstream outfile((char*)ref_file_name.c_str());
	char base;
	if (outfile.is_open()){
		if (!args.quiet)
			cout<<"Generating reference genome..."<<endl;
		outfile<<">" + args.genome_id;
    	for (unsigned int i = 0; i < args.num_chroms; i++){
    		string chrom_num = static_cast<ostringstream*>( 
    			&(ostringstream() << i) )->str();
    		outfile<<"\n>chromosome_" + chrom_num;
    		for (unsigned long j = 0; j < args.chrom_size; j++){
    			if (j % LINE_WIDTH == 0)
    				outfile<<'\n';
    			switch(rand() % 4){
    				case 0: base = 'A'; break;
    				case 1: base = 'C'; break;
    				case 2: base = 'G'; break;
    				case 3: base = 'T'; break;
    			}
    			genome[i][j] = base;
    			outfile<<base;
    		}
    	}
    	outfile<<"\n";
    	outfile.close();
  	}
}

void write_private_genome(vector<vector<char>>& genome){
	ofstream outfile((char*)priv_file_name.c_str());
	if (outfile.is_open()){
		if (!args.quiet)
			cout<<"Writing private genome..."<<endl;
		outfile<<">" + args.genome_id;
		for (unsigned int i = 0; i < args.num_chroms; i++){
			string chrom_num = static_cast<ostringstream*>( 
    			&(ostringstream() << i) )->str();
			outfile<<"\n>chromosome_" + chrom_num;
			for (unsigned long j = 0; j < genome[i].size(); j++){
				if (j % LINE_WIDTH == 0)
					outfile<<'\n';
				outfile<<genome[i][j];
			}
		}
		outfile<<"\n";
		outfile.close();
	}
}

void write_reads(vector<vector<char>>& genome){
	if (!args.quiet)
		cout<<"Generating reads..."<<endl;

	ofstream outfile((char*)reads_file_name.c_str());
	if (outfile.is_open()){
		outfile<<'>'<<args.genome_id<<'\n';
		// generate reads from genome
		unsigned long num_reads = 
			((args.coverage * args.genome_size) / args.read_length);
		string read1, read2;
		for (unsigned long i = 0; i < num_reads / 2; i++){
			unsigned int chrom_num = rand() % genome.size();
			unsigned int gap = 90 + (rand() % 20);
			unsigned long index1 = (unsigned long)
				(random()*(genome[chrom_num].size() - gap - 2*args.read_length));
			unsigned long index2 = index1 + args.read_length + gap;
			if (random() <= args.read_garbage_rate){
				read1 = get_garbage_read();	
				// read1 = "garbage";	
			}
			else{
				read1 = get_slice(
					chrom_num, index1, 
					index1 + args.read_length, 
					genome[chrom_num]
				);
				// introduce errors
				for (unsigned int i = 0; i < read1.size(); i++){
					if (random() < args.read_error_rate){
						// read1[i] = random_snp(read[i]);
						read1[i] = '~';
					}
				}
			}
			if (random() <= args.read_garbage_rate){
				read2 = get_garbage_read();	
				// read2 = "garbage";
			}
			else{
				read2 = get_slice(
					chrom_num, index2,
					index2 + args.read_length,
					genome[chrom_num]
				);
				// introduce errors
				for (unsigned int i = 0; i < read2.size(); i++){
					if (random() < args.read_error_rate){
						read2[i] = random_snp(read2[i]);
						// read2[i] = '~';
					}
				}
			}
			outfile<<read1<<','<<read2<<'\n';
		}
	}
	outfile.close();
	remove(".temp");
}

void generate_mutations(vector<vector<char>>& genome){
	ofstream outfile((char*)ans_file_name.c_str());
	if (outfile.is_open()){
		outfile<<'>'<<args.genome_id<<'\n';
	}
	outfile.close();
	generate_indels(genome);
	generate_snps(genome);
}

/*void generate_indels(vector<vector<char>>& genome){
	ofstream outfile((char*)ans_file_name.c_str(), ios::app);
	if (outfile.is_open()){
		unsigned long num_indels = args.genome_size * args.indel_rate;
		if (!args.quiet)
			cout<<"Generating Inserts..."<<endl;
		outfile<<">INSERT\n";
		for (unsigned int chrom = 0; chrom < genome.size(); chrom++){
			if (!args.quiet)
				cout<<"\tChromosome "<<chrom<<endl;
			unsigned long num_inserts = num_indels / 2;
				// (genome[chrom].size() * args.indel_rate) / 2;
			vector<unsigned long> ins_indeces(num_inserts);
			for (unsigned long j = 0; j < num_inserts; j++){
				unsigned long index = random() * (genome[chrom].size() - 1);
				ins_indeces[j] = index;
			}
			sort(ins_indeces.begin(), ins_indeces.end());
			for (unsigned long j = 0; j < ins_indeces.size(); j++){
				outfile<<chrom<<','<<genome[chrom][ins_indeces[j]]<<',';
				genome[chrom].insert(
					genome[chrom].begin() + ins_indeces[j], random_allele());
				// genome[chrom][ins_indeces[j]] = 'i';
				outfile<<genome[chrom][ins_indeces[j]]<<','<<ins_indeces[j]<<'\n';
			}
		}
		if (!args.quiet)
			cout<<"Generating Deletes..."<<endl;
		outfile<<">DELETE\n";
		for (unsigned int chrom = 0; chrom < genome.size(); chrom++){
			if (!args.quiet)
				cout<<"\tChromosome "<<chrom<<endl;
			unsigned long num_deletes = num_indels / 2;
				// (genome[chrom].size() * args.indel_rate) / 2;
			vector<unsigned long> del_indeces(num_deletes);
			for (unsigned long j = 0; j < del_indeces.size(); j++){
				unsigned long index = random() * (genome[chrom].size() - 1);
				del_indeces[j] = index;
			}
			sort(del_indeces.begin(), del_indeces.end());
			for (unsigned long j = 0; j < del_indeces.size(); j++){
				outfile<<chrom<<','<<genome[chrom][del_indeces[j]]<<',';
				if (del_indeces[j] >= genome[chrom].size())
					genome[chrom].erase(genome[chrom].end() - 1);
				else
					genome[chrom].erase(genome[chrom].begin() + del_indeces[j]);
				// genome[chrom][del_indeces[j]] = 'd';
				outfile<<genome[chrom][del_indeces[j]]<<','<<del_indeces[j]<<'\n';
			}

		}
	}
	outfile.close();
}*/

void generate_indels(vector<vector<char>>& genome){
	const unsigned short INDEL_MAX_LENGTH = 5;
	if (!args.quiet)
		cout<<"Generating Indels..."<<endl;
	
	// generate insert/delete indeces
	unsigned long num_indels = args.genome_size * args.indel_rate;
	// vector<pair<unsigned int, unsigned long>> inserts(num_indels / 2);
	// vector<pair<unsigned int, unsigned long>> deletes(num_indels / 2);
	// (chrom, index , length)
	typedef tuple<unsigned int, unsigned long, unsigned short> indel_tuple; 
	vector<indel_tuple> inserts(num_indels / 2);
	vector<indel_tuple> deletes(num_indels / 2);
	if (!args.quiet)
		cout<<"\tgenerating indel indeces..."<<endl;
	unsigned int chrom;
	unsigned long index;
	unsigned short indel_length;
	for (unsigned long i = 0; i < num_indels / 2; i++){
		chrom = rand() % genome.size();
		index = random() * (genome.size() - 1);
		indel_length = rand() % INDEL_MAX_LENGTH;
		// inserts[i] = pair<unsigned int, unsigned long> (chrom, index);
		inserts[i] = make_tuple(chrom, index, indel_length);
		index = random() * (genome.size() - 1);
		// deletes[i] = pair<unsigned int, unsigned long> (chrom, index);
		deletes[i] = make_tuple(chrom, index, indel_length);
	}

	// write random inserts/deletes to file
	ofstream outfile((char*)ans_file_name.c_str(), ios::app);
	if (outfile.is_open()){
		if (!args.quiet)
			cout<<"\twriting inserts..."<<endl;
		for (unsigned long i = 0; i < inserts.size(); i++){
			// unsigned int chrom = inserts[i].first;
			unsigned int chrom = get<0>(inserts[i]);
			// unsigned long index = inserts[i].second;
			unsigned long index = get<1>(inserts[i]);
			string alleles = random_alleles(get<2>(inserts[i]));
			// genome[chrom].insert(
			// 	genome[chrom].begin() + index, alleles.begin(), alleles.end());
			genome[chrom].insert(
				genome[chrom].begin() + index, alleles.begin(), alleles.end());
		}

		if (!args.quiet)
				cout<<"\twriting deletes..."<<endl;
		for (unsigned long i = 0; i < deletes.size(); i++){
			;
		}
	}
	outfile.close();
}

void generate_snps(vector<vector<char>>& genome){
	ofstream outfile((char*)ans_file_name.c_str(), ios::app);
	if (outfile.is_open()){
		if (!args.quiet)
			cout<<"Generating SNPs..."<<endl;
		outfile<<">SNP\n";
		for (unsigned int chrom = 0; chrom < genome.size(); chrom++){
			if (!args.quiet)
				cout<<"\twriting chrom "<<chrom<<endl;
			unsigned long num_snps = genome[chrom].size() * args.snp_rate;
			vector<unsigned long> rand_indeces(num_snps);
			for (unsigned long j = 0; j < rand_indeces.size(); j++){
				rand_indeces[j] = random() * (genome[chrom].size() - 1);
			}
			sort(rand_indeces.begin(), rand_indeces.end());
			for (unsigned long j = 0; j < rand_indeces.size(); j++){
				outfile<<chrom<<','<<genome[chrom][rand_indeces[j]]<<',';
				// genome[chrom][rand_indeces[j]] = 
				// 	random_snp(genome[chrom][rand_indeces[j]]);
				genome[chrom][rand_indeces[j]] = '#';
				outfile<<genome[chrom][rand_indeces[j]]<<','<<rand_indeces[j]<<'\n';
			}
		}
	}
	outfile.close();
}

string get_garbage_read(){
	string read;
	for (unsigned int i = 0; i < args.read_length; i++){
		read += random_allele();
	}
	return read;
}

string get_slice(int chromosome, 
	unsigned long start, unsigned long end, const vector<char>& genome)
{
	string slice = "";
	for (unsigned long i = start; i < end; i++){
		slice += genome[i];
	}
	return string(slice);
}

string random_alleles(unsigned int length){
	string alleles;
	for (unsigned int i = 0; i < length; i++){
		alleles += random_allele();
	}
	return alleles;
}

char random_allele(){
	switch(rand() % 4){
		case 0: return 'A'; break;
		case 1: return 'C'; break;
		case 2: return 'G'; break;
		case 3: return 'T'; break;
	}
}

char random_snp(char base){
	switch(base){
		case 'A':
			switch(rand() % 3){
				case 0: return 'C'; break;
				case 1: return 'G'; break;
				case 2: return 'T'; break;
			} 
			break;
		case 'C': 
			switch(rand() % 3){
				case 0: return 'A'; break;
				case 1: return 'G'; break;
				case 2: return 'T'; break;
			}
			break;
		case 'G':
			switch(rand() % 3){
				case 0: return 'A'; break;
				case 1: return 'C'; break;
				case 2: return 'T'; break;
			}
			break;
		case 'T': 
			switch(rand() % 3){
				case 0: return 'A'; break;
				case 1: return 'C'; break;
				case 2: return 'G'; break;
			}
			break;
		default: 
			throw invalid_argument(
				"Valid values: {'A', 'C', 'G', 'T'}. Given: " + base);
	}
}

double random(){
	return rand() / double(RAND_MAX);
}

}; //end namespace


int main(int argc, char** argv){	
  	parse_args(argc, argv);
  	srand(time(NULL));
	
  	vector<vector<char>> genome(args.num_chroms);
	for (unsigned int i = 0; i < args.num_chroms; i++){
		genome[i] = vector<char>(args.chrom_size);
	}

  	write_ref_genome(genome);
  	generate_mutations(genome);
  	write_private_genome(genome);
  	write_reads(genome);
  	cout<<"DONE"<<endl;

	return 0;
}