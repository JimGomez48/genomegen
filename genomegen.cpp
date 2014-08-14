#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include "genomegen.hpp"
#include "tclap/CmdLine.h"

using namespace std;


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
			"The number of chromosomes, by default in thousands of bp. "
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

string generate_ref_genome(string id, int num_chroms, unsigned long chrom_size){
	string file_name = "ref_" + id + ".txt";
	char *c_file_name = (char*)file_name.c_str();
	ofstream outfile;
	outfile.open(c_file_name);
	string line;
	if (outfile.is_open()){
		outfile<<">" + id;
    	for (int i = 0; i < num_chroms; i++){
    		string num = static_cast<ostringstream*>( 
    			&(ostringstream() << i) )->str();
    		outfile<<"\n>chromosome_" + num;
    		for (unsigned long j = 0; j < chrom_size; j++){
    			if (j % 80 == 0)
    				outfile<<"\n";
    			switch(rand() % 4){
    				case 0: outfile<<"A"; break;
    				case 1: outfile<<"C"; break;
    				case 2: outfile<<"G"; break;
    				case 3: outfile<<"T"; break;
    			}
    		}
    	}
    	outfile<<"\n";
    	outfile.close();
  	}

	return file_name;
}

int main(int argc, char** argv){
  	struct cmd_args args = parse_args(argc, argv);
  	string filename = generate_ref_genome(
  		args.genome_id, args.num_chroms, args.chrom_size);
  	cout<<filename<<endl;

	return 0;
}