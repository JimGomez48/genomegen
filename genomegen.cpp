#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <stdexcept>

#include "genomegen.hpp"
#include "tclap/CmdLine.h"

using namespace std;
using namespace genome_gen;

namespace genome_gen{

static const string REF_PRE = "ref_";
static const string PRIV_PRE = "private_";
static const string READS_PRE = "reads_";
static const string ANS_PRE = "ans_";

cmd_args parse_args(int argc, char** argv){
	struct cmd_args args;
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
			"The number of chromosomes to generate for the genome",
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
		// print arg values to user, if valid
		cout<<"\ngenome-ID:\t"<<args.genome_id<<endl;
    	cout<<"num-chroms:\t"<<args.num_chroms<<endl;
    	if (args.scale == 'k'){
        	cout<<"chrom-size:\t"<<args.chrom_size<<" thousand"<<endl;
        	args.chrom_size *= 1000;
        }
    	else if (args.scale == 'm'){
        	cout<<"chrom-size:\t"<<args.chrom_size<<" million"<<endl;
        	args.chrom_size *= 1000000;
        }
    	else if (args.scale == 'b'){
        	cout<<"chrom-size:\t"<<args.chrom_size<<" billion"<<endl;
        	args.chrom_size *= 1000000000;
        }
        else{
        	throw TCLAP::ArgException(
				"Invalid argument value. Choices <k, m, b>",
				"--scale"
			);
        }
		
		return args;
	}
	catch(TCLAP::ArgException &e){
		cout<<"ARG PARSE FAILURE:\n\t"<<e.what()<<endl;
		exit(1);
	}
}

void generate_ref_genome(string id, int num_chroms, unsigned long chrom_size, vector<char>& genome){
	string file_name = REF_PRE + id + ".txt";
	char *c_file_name = (char*)file_name.c_str();
	ofstream outfile(c_file_name);
	if (outfile.is_open()){
		cout<<"Generating reference genome..."<<endl;
		outfile<<">" + id;
    	for (int i = 0; i < num_chroms; i++){
    		string num = static_cast<ostringstream*>( 
    			&(ostringstream() << i) )->str();
    		outfile<<"\n>chromosome_" + num;
    		for (unsigned long j = 0; j < chrom_size; j++){
    			if (j % 80 == 0)
    				outfile<<"\n";
    			switch(rand() % 4){
    				case 0: outfile<<'A'; genome[i * chrom_size + j] = 'A'; break;
    				case 1: outfile<<'C'; genome[i * chrom_size + j] = 'A'; break;
    				case 2: outfile<<'G'; genome[i * chrom_size + j] = 'A'; break;
    				case 3: outfile<<'T'; genome[i * chrom_size + j] = 'A'; break;
    			}
    		}
    	}
    	outfile<<"\n";
    	outfile.close();
  	}
}

void generate_copies(vector<char>& genome){
	// TODO
}

void generate_inversions(vector<char>& genome){
	// TODO
}

void generate_insertions(vector<char>& genome){
	// TODO
}

void generate_deletions(vector<char>& genome){
	// TODO
}

void generate_snps(vector<char>& genome){
	const float SNP_RATE = 0.003; // 0.3%
	cout<<"Generating SNPs..."<<endl;
	for (unsigned long i = 0; i < genome.size() * SNP_RATE; i++){
		unsigned long index = rand_num() * 3000000000ul;
		genome[i] = random_snp(genome[i]);
	}
}

void generate_alus(vector<char>& genome){
	// TODO
}

void generate_strs(vector<char>& genome){
	// TODO
}

void generate_reads(string id){
	// TODO
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
			throw invalid_argument("Valid values: {'A', 'C', 'G', 'T'}");
	}
}

double rand_num(){
	return rand() / double(RAND_MAX);
}

}; //end namespace

int main(int argc, char** argv){
	// initialize rand() and get args
  	struct cmd_args args = parse_args(argc, argv);
  	srand(time(NULL));
  	// initialize container for in-memory genome
  	vector<char> genome;
	genome.reserve(args.num_chroms * args.chrom_size);
	// generate reference genome and mutations
  	generate_ref_genome(
  		args.genome_id, args.num_chroms, args.chrom_size, genome);
  	generate_snps(genome);

	return 0;
}