#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <random>

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
		cmd.parse(argc, argv);
		if (coverage.getValue() < 0)
			throw TCLAP::ArgParseException("value cannot be < 0", "coverage");
		if (read_error_rate.getValue() < 0 || read_error_rate.getValue() > 1)
			throw TCLAP::ArgParseException("out of valid range: [0, 1]", "read-error");
		if (read_garbage_rate.getValue() < 0 || read_garbage_rate.getValue() > 1)
			throw TCLAP::ArgParseException("out of valid range: [0, 1]", "read-garbage");

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
		// generate random errors for reads
		unsigned long num_errors = args.genome_size * args.read_error_rate;
		cout<<"#errors: "<<num_errors<<endl;
		for (unsigned long i = 0; i < num_errors; i++){
			unsigned int chrom_num = rand() % genome.size();
			unsigned long index = 
				(unsigned long)(random()*genome[chrom_num].size());
			genome[chrom_num][index] = random_snp(genome[chrom_num][index]);
			// genome[chrom_num][index] = '~';
		}
		// generate reads from genome
		unsigned long num_reads = 
			((args.coverage * args.genome_size) / args.read_length);
		for (unsigned long i = 0; i < num_reads / 2; i++){
			string read1, read2;
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
					chrom_num, index1, index1 + args.read_length, genome[chrom_num]);
			}
			if (random() <= args.read_garbage_rate){
				read2 = get_garbage_read();	
				// read2 = "garbage";
			}
			else{
				read2 = get_slice(
					chrom_num, index2, index2 + args.read_length, genome[chrom_num]);
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
	for (int i = 0; i < genome.size(); i++){
		if (!args.quiet)
			cout<<"Generating mutations in chromosome: "<<i<<endl;
		generate_insertions(genome[i]);
		generate_deletions(genome[i]);
		generate_snps(genome[i]);
	}
}

void generate_insertions(vector<char>& chromosome){
	// TODO
}

void generate_deletions(vector<char>& chromosome){
	// TODO
}

void generate_snps(vector<char>& chromosome){
	const float SNP_RATE = 0.003; // 0.3%
	const unsigned long NUM_SNPS = 
		(unsigned long)(chromosome.size() * SNP_RATE);
	if (!args.quiet)
		cout<<"\tGenerating SNPs..."<<endl;
	for (unsigned long i = 0; i < NUM_SNPS; i++){
		unsigned long index = random() * args.chrom_size;
		chromosome[index] = random_snp(chromosome[index]);
	}
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
  	// generate_mutations(genome);
  	write_private_genome(genome);
  	write_reads(genome);

  	// const unsigned long ITERATIONS = 1000000000;
  	// double rand_num;
  	// ranlux48_base engine;
  	// uniform_int_distribution<unsigned long> generator(0ul, 3000000000ul);
  	// // normal_distribution<double> generator(5, 3);
  	// // cout<<"mean: "<<generator.mean()<<endl;
  	// // cout<<"sdev: "<<generator.stddev()<<endl;
  	// for (unsigned long i = 0; i < ITERATIONS; i++){
  	// 	// rand_num = genome_gen::random() * 3000000000;
  	// 	rand_num = generator(engine);
  	// 	// cout<<rand_num<<endl;
  	// }

	return 0;
}