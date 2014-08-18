#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <stdexcept>
#include <vector>

#include "genomegen.hpp"
#include "tclap/CmdLine.h"

using namespace std;
using namespace genome_gen;

namespace genome_gen{

static const string REF_PRE = "ref_";
static const string PRIV_PRE = "private_";
static const string READS_PRE = "reads_";
static const string ANS_PRE = "ans_";

struct cmd_args{
	string genome_id;
	int num_chroms;
	unsigned long chrom_size;
	char scale;
	float error_rate;
	float garbage_rate;
	float coverage;
} args;

void print_args(){
	cout<<"\ngenome-ID:\t"<<args.genome_id<<endl;
	cout<<"num-chroms:\t"<<args.num_chroms<<endl;
	cout<<"chrom-size:\t"<<args.chrom_size<<endl;
	cout<<"scale:\t"<<args.scale<<endl;
}

void parse_args(int argc, char** argv){
	try{
		// Define Command line parser and description
		TCLAP::CmdLine cmd(
			"DESCRIPTION\nThis program generates a reference and mutated donor "
            "genome as a set of files. The following files are "
            "created:\n1) reference genome \'ref_<genome_id>.txt\'\n2) mutated "
            "donor genome \'private_<genome_id>.txt\'\n3) paired-end reads "
            "\'reads_<genome_id>.txt\' from donor genome\n4) mutation answer "
			"key \'ans_<genome_id>.txt\'", 
			' '
		);
		// define args
		TCLAP::UnlabeledValueArg<string> genome_id(
			"genome_id", 
			"The name or ID of this genome", 
			true, 
			"genome1",
			"genome_id(string)"
		);
		TCLAP::UnlabeledValueArg<int> num_chroms(
			"num_chromosomes",
			"The number of chromosomes to generate for the genome",
			true,
			1,
			"num_chroms(int)"
		);
		TCLAP::UnlabeledValueArg<unsigned long> chrom_size(
			"chromosome_size",
			"The size of each chromosome, by default in thousands of bp. "
			"Change scale with -s option",
			true,
			1,
			"chrom_size(int)"
		);
		TCLAP::ValueArg<char> scale(
			"s", "scale",
			"The amount to scale chromosome-size by. k: 1-thousand, m: "
			"1-million, b: 1-billion. By default, scale is 1-thousand (k).",
			false,
			'k',
			"k, m, b"
		);
		// bind args to parser and parse
		cmd.add(genome_id);
		cmd.add(num_chroms);
		cmd.add(chrom_size);
		cmd.add(scale);		
		cmd.parse(argc, argv);
		// set args struct and return
		args.genome_id = genome_id.getValue();
		args.num_chroms = num_chroms.getValue();
		args.chrom_size = chrom_size.getValue();
		args.scale = scale.getValue();
		args.error_rate = 0.0;
		args.garbage_rate = 0.0;
		// print arg values to user, if valid
		cout<<"\ngenome-ID:\t"<<args.genome_id<<endl;
    	cout<<"num-chroms:\t"<<args.num_chroms<<endl;
    	switch(args.scale){
    		case 'k': 
    			cout<<"chrom-size:\t"<<args.chrom_size<<" thousand"<<endl;
        		args.chrom_size *= 1000;
    			break;
    		case 'm': 
    			cout<<"chrom-size:\t"<<args.chrom_size<<" million"<<endl;
        		args.chrom_size *= 1000000;
    			break;
    		case 'b': 
    			cout<<"chrom-size:\t"<<args.chrom_size<<" billion"<<endl;
        		args.chrom_size *= 1000000000;
    			break;
    		default:
    			throw TCLAP::ArgException(
					"Invalid argument value. Choices <k, m, b>",
					"--scale"
				);
    	}
	}
	catch(TCLAP::ArgException &e){
		cout<<"ARG PARSE FAILURE:\n\t"<<e.what()<<endl;
		exit(1);
	}
}

void write_ref_genome(vector<vector<char>>& genome){
	string file_name = REF_PRE + args.genome_id + ".txt";
	ofstream outfile((char*)file_name.c_str());
	char base;
	if (outfile.is_open()){
		cout<<"Generating reference genome..."<<endl;
		outfile<<">" + args.genome_id;
    	for (int i = 0; i < args.num_chroms; i++){
    		string chrom_num = static_cast<ostringstream*>( 
    			&(ostringstream() << i) )->str();
    		outfile<<"\n>chromosome_" + chrom_num;
    		for (unsigned long j = 0; j < args.chrom_size; j++){
    			if (j % 80 == 0)
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
	string file_name = PRIV_PRE + args.genome_id + ".txt";
	ofstream outfile((char*)file_name.c_str());
	if (outfile.is_open()){
		cout<<"Writing private genome..."<<endl;
		outfile<<">" + args.genome_id;
		for (int i = 0; i < args.num_chroms; i++){
			string chrom_num = static_cast<ostringstream*>( 
    			&(ostringstream() << i) )->str();
			outfile<<"\n>chromosome_" + chrom_num;
			for (unsigned long j = 0; j < args.chrom_size; j++){
				if (j % 80 == 0)
					outfile<<'\n';
				outfile<<genome[i][j];
			}
		}
		outfile<<"\n";
		outfile.close();
	}
}

void write_reads(vector<vector<char>>& genome){
	string file_name = READS_PRE + args.genome_id + ".txt";
	ofstream outfile((char*)file_name.c_str());
	cout<<"Generating reads..."<<endl;
	if (outfile.is_open()){
		string read1 = "";
		string read2 = "";
		// TODO fix everything
		for (unsigned long i = 0; i < 50; i++){
			unsigned long index1 = 0;
			unsigned long index2 = 0;
			int gap = 0;
		}
	}
	outfile.close();
}


void generate_insertions(vector<vector<char>>& genome){
	// TODO
}

void generate_deletions(vector<vector<char>>& genome){
	// TODO
}

void generate_snps(vector<vector<char>>& genome){
	const float SNP_RATE = 0.003; // 0.3%
	const unsigned long NUM_SNPS = 
		(unsigned long)(args.num_chroms * args.chrom_size * SNP_RATE);
	cout<<"Generating SNPs..."<<endl;
	for (unsigned long i = 0; i < NUM_SNPS; i++){
		int chrom_num = rand_num() * args.num_chroms;
		unsigned long index = rand_num() * args.chrom_size;
		genome[chrom_num][index] = random_snp(genome[chrom_num][index]);
	}
}

// void generate_copies(char** genome){
// 	// TODO
// }

// void generate_inversions(char** genome){
// 	// TODO
// }


// void generate_alus(char** genome){
// 	// TODO
// }

// void generate_strs(char** genome){
// 	// TODO
// }

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
			string msg = "Valid values: {'A', 'C', 'G', 'T'}. Given: " + base;
			throw invalid_argument(msg);
			// throw invalid_argument("Valid values: {'A', 'C', 'G', 'T'}");
	}
}

double rand_num(){
	return rand() / double(RAND_MAX);
}

}; //end namespace

int main(int argc, char** argv){
	// initialize rand() and get args
  	parse_args(argc, argv);
  	srand(time(NULL));
	// generate reference genome and mutations
  	vector<vector<char>> genome(args.num_chroms);
	for (int i = 0; i < args.num_chroms; i++){
		genome[i] = vector<char>(args.chrom_size);
	}

  	write_ref_genome(genome);
  	generate_snps(genome);
  	write_private_genome(genome);
  	// write_reads(genome);
  	
	return 0;
}