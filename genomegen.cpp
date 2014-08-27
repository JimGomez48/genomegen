#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <tuple>

#include "genomegen.hpp"
#include "tclap/CmdLine.h"

using namespace std;
using namespace genome_gen;


namespace genome_gen {

static const string REF_PRE = "ref_";
static const string PRIV_PRE = "priv_";
static const string READS_PRE = "reads_";
static const string ANS_PRE = "ans_";
static const string SUFFIX = ".txt";
static const unsigned int LINE_WIDTH = 80;
static const char DEL_MARKER = '-';

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
	double insert_rate;
	double delete_rate;
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
			"vary by up to +/- 10 bp. Default: 100",
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
		TCLAP::ValueArg<double> insert_rate(
			"", "insert",
			"The rate at which inserts will be introduced into the mutated "
			"genome. Default: 0.001",
			false,
			0.001,
			"float"
		);
		TCLAP::ValueArg<double> delete_rate(
			"", "delete",
			"The rate at which deletes will be introduced into the mutated "
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
		// cmd.add(insert_rate);
		// cmd.add(delete_rate);
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
		// if (insert_rate.getValue() < 0 || insert_rate.getValue() > 1)
		// 	throw TCLAP::ArgParseException("out of valid range: [0, 1]", "insert");
		// if (delete_rate.getValue() < 0 || delete_rate.getValue() > 1)
		// 	throw TCLAP::ArgParseException("out of valid range: [0, 1]", "delete");

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
    	args.insert_rate 		= args.indel_rate / 2.0;
    	args.delete_rate 		= args.indel_rate / 2.0;

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
	// cout<<"insert-rate:\t"<<args.insert_rate<<endl;
	// cout<<"delete-rate:\t"<<args.delete_rate<<endl;
	cout<<"================================="<<endl;
}

void write_ref_genome(vector<vector<char>>& genome){
	if (args.genome_size){
		ofstream outfile((char*)ref_file_name.c_str());
		if (outfile.is_open()){
			if (!args.quiet)
				cout<<"Generating reference genome..."<<endl;
			outfile<<">" + args.genome_id;
	    	for (unsigned int i = 0; i < args.num_chroms; i++){
	    		string chrom_num = static_cast<ostringstream*>( 
	    			&(ostringstream() << i) )->str();
	    		outfile<<"\n>chromosome_" + chrom_num;
				char base;
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
}

void write_private_genome(vector<vector<char>>& genome){
	if (args.genome_size){
		ofstream outfile((char*)priv_file_name.c_str());
		if (outfile.is_open()){
			if (!args.quiet)
				cout<<"Writing private genome..."<<endl;
			outfile<<">" + args.genome_id;
			for (unsigned int chrom = 0; chrom < genome.size(); chrom++){
				string chrom_num = static_cast<ostringstream*>( 
	    			&(ostringstream() << chrom) )->str();
				outfile<<"\n>chromosome_" + chrom_num;
				unsigned long line_pos = 0;
				for (unsigned long j = 0; j < genome[chrom].size(); j++){
					if (line_pos % LINE_WIDTH == 0)
						outfile<<'\n';
					if (genome[chrom][j] != DEL_MARKER){
						outfile<<genome[chrom][j];
						line_pos++;
					}
					// outfile<<genome[chrom][j];
					// line_pos++;
				}
			}
			outfile<<"\n";
			outfile.close();
		}
	}
}

void write_reads(vector<vector<char>>& genome){
	if (args.coverage && args.genome_size){
		if (!args.quiet)
			cout<<"Generating reads..."<<endl;
		
		if (!args.quiet)
			cout<<"\tremoving place holders from mutated genome..."<<endl;
		for (unsigned int chrom = 0; chrom < genome.size(); chrom++){
			vector<char>::iterator i;
			for (i = genome[chrom].begin(); i < genome[chrom].end(); i++){
				if ((*i) == DEL_MARKER){
					i = genome[chrom].erase(i) - 1;
				}
			}
		}

		ofstream outfile((char*)reads_file_name.c_str());
		if (outfile.is_open()){
			if (!args.quiet)
				cout<<"\twriting reads from mutated genome..."<<endl;
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
				
				if (random() <= args.read_garbage_rate)
					read1 = get_garbage_read();	
				else{
					read1 = get_slice(
						chrom_num, index1, 
						index1 + args.read_length, 
						genome
					);
					// introduce errors
					for (unsigned int i = 0; i < read1.size(); i++){
						if (random() < args.read_error_rate)
							read1[i] = random_snp(read1[i]);
					}
				}
				
				if (random() <= args.read_garbage_rate)
					read2 = get_garbage_read();	
				else{
					read2 = get_slice(
						chrom_num, index2,
						index2 + args.read_length,
						genome
					);
					// introduce errors
					for (unsigned int i = 0; i < read2.size(); i++){
						if (random() < args.read_error_rate)
							read2[i] = random_snp(read2[i]);
					}
				}
				outfile<<read1<<','<<read2<<'\n';
			}
		}
		outfile.close();
	}
}

void generate_mutations(vector<vector<char>>& genome){
	if (args.genome_size){
		ofstream outfile((char*)ans_file_name.c_str());
		if (outfile.is_open()){
			outfile<<'>'<<args.genome_id<<'\n';
		}
		outfile.close();
	}
	generate_indels(genome);
	// generate_deletes(genome);
	// generate_inserts(genome);
	generate_snps(genome);
}

// 			  chrom_num,   	index, 		   alleles
typedef tuple<unsigned int, unsigned long, string> indel_tuple; 

bool indel_tuple_comparator(indel_tuple a, indel_tuple b){
	if (get<0>(a) < get<0>(b))
		return true;
	if (get<0>(a) == get<0>(b))
		return get<1>(a) < get<1>(b);
	return false;
}

void generate_indels(vector<vector<char>>& genome){
	if (args.indel_rate){
		const unsigned int INDEL_MAX_LENGTH = 5;
		if (!args.quiet)
			cout<<"Generating Indels"<<endl;

		const unsigned long num_indels = args.genome_size * args.indel_rate;
		const unsigned long num_inserts = num_indels / 2 + 1;
		const unsigned long num_deletes = num_indels / 2 + 1;
		vector<indel_tuple> inserts(num_inserts);
		vector<indel_tuple> deletes(num_deletes);
		if (!args.quiet)
			cout<<"\tmutating genome..."<<endl;
		for (unsigned long i = 0; i < num_indels; i++){
			const unsigned int length = (rand() % INDEL_MAX_LENGTH) + 1;
			unsigned int chrom_num = rand() % genome.size();
			vector<char> chrom = genome[chrom_num];
			// delete temp (replace with dummies)
			const unsigned long del_index = random() * (chrom.size() - length);
			const string del_alleles = get_slice(
				chrom_num, 
				del_index, 
				del_index + length,
				genome
			);
			vector<char>::iterator del_start = chrom.begin() + del_index;
			vector<char>::iterator del_end = chrom.begin() + del_index + length;
			chrom.erase(del_start, del_end);
			const string dummy = string(length, DEL_MARKER);
			chrom.insert(del_start, dummy.begin(), dummy.end());
			// insert
			unsigned long ins_index;
			cout<<"\tgetting insert index"<<endl;
			do{
				ins_index = random() * (chrom.size() - 1);
			} while (ins_index >= del_index && del_index + length >= ins_index);
			cout<<"\tdone"<<endl;
			const string ins_alleles = random_alleles(length);
			chrom.insert(
				chrom.begin() + ins_index, 
				ins_alleles.begin(), 
				ins_alleles.end()
			);
			// remove dummy placeholders
			vector<char>::iterator iter = chrom.begin() + del_index;
			if (ins_index < del_index){
				iter += length;
			}
			else if (ins_index > del_index + length){
				;
			}
			else{
				cout<<"Iter: "<<i<<endl;
				cout<<"Insert Index: "<<ins_index<<endl;
				cout<<"Delete range: ["<<del_index<<":"<<del_index + length<<"]"<<endl;
				throw runtime_error("insert index overlaps delete range");
			}
			vector<char>::iterator stop = iter + length;
			for (; iter < stop; iter++){
				if ((*iter) == DEL_MARKER){
					iter = chrom.erase(iter) - 1;
				}
			}
			inserts[i] = make_tuple(chrom_num, ins_index, ins_alleles);
			deletes[i] = make_tuple(chrom_num, del_index, del_alleles);
		}
		sort(inserts.begin(), inserts.end(), indel_tuple_comparator);
		sort(deletes.begin(), deletes.end(), indel_tuple_comparator);

		ofstream outfile((char*)ans_file_name.c_str(), ios::app);
		// write inserts to file
		if (!args.quiet)
			cout<<"\twriting insertions to file..."<<endl;
		outfile<<">INSERT\n";
		enum {CHROM, INDEX, ALLELES};
		vector<indel_tuple>::iterator i;
		for (i = inserts.begin(); i < inserts.end(); i++){
			const unsigned int chrom = get<CHROM>((*i));
			const unsigned long index = get<INDEX>((*i));
			const string alleles = get<ALLELES>((*i));
			outfile<<chrom<<','<<alleles<<','<<index<<'\n';
		}

		// write deletes to file
		if (!args.quiet)
			cout<<"\twriting deletions to file..."<<endl;
		outfile<<">DELETE\n";
		vector<indel_tuple>::iterator d;
		for (d = deletes.begin(); d < deletes.end(); d++){
			const unsigned int chrom = get<CHROM>((*d));
			const unsigned long index = get<INDEX>((*d));
			const string alleles = get<ALLELES>((*d));
			outfile<<chrom<<','<<alleles<<','<<index<<'\n';
		}
		outfile.close();
	}
}

// void generate_inserts(vector<vector<char>>& genome){
// 	if (args.insert_rate && args.genome_size){
// 		const unsigned short INS_MAX_LENGTH = 5;

// 		if (!args.quiet)
// 			cout<<"Generating Inserts..."<<endl;
		
// 		// generate inserts indeces and values
// 		unsigned long num_inserts = args.genome_size * args.insert_rate;
// 		vector<indel_tuple> inserts(num_inserts);
// 		if (!args.quiet)
// 			cout<<"\tgenerating random indeces..."<<endl;
// 		for (unsigned long i = 0; i < num_inserts; i++){
// 			unsigned int chrom = rand() % genome.size();
// 			unsigned long index = random() * (genome[chrom].size() - 1);
// 			unsigned short ins_length = (rand() % INS_MAX_LENGTH) + 1;
// 			const string alleles = random_alleles(ins_length);
// 			inserts[i] = make_tuple(chrom, index, alleles);
// 		}
// 		sort(inserts.begin(), inserts.end(), indel_tuple_comparator);

// 		// perform insertions
// 		if (!args.quiet)
// 			cout<<"\tinserting at generated indeces..."<<endl;
// 		enum {CHROM, INDEX, ALLELES};
// 		// end to begin to ensure correct indexing relative to reference genome
// 		vector<indel_tuple>::reverse_iterator i;
// 		for (i = inserts.rbegin(); i < inserts.rend(); i++){
// 			const unsigned int chrom = get<CHROM>((*i));
// 			const unsigned long index = get<INDEX>((*i));
// 			const string alleles = get<ALLELES>((*i));
// 			genome[chrom].insert(
// 				genome[chrom].begin() + index, 
// 				alleles.begin(), 
// 				alleles.end()
// 			);
// 		}

// 		// write inserts to file
// 		ofstream outfile((char*)ans_file_name.c_str(), ios::app);
// 		if (outfile.is_open()){
// 			if (!args.quiet)
// 				cout<<"\twriting inserts..."<<endl;
// 			outfile<<">INSERT\n";
// 			vector<indel_tuple>::iterator i;
// 			for (i = inserts.begin(); i < inserts.end(); i++){
// 				const unsigned int chrom = get<CHROM>((*i));
// 				const unsigned long index = get<INDEX>((*i));
// 				const string alleles = get<ALLELES>((*i));
// 				outfile<<chrom<<','<<alleles<<','<<index<<'\n';
// 			}
// 		}
// 		outfile.close();
// 	}
// }

// void generate_deletes(vector<vector<char>>& genome){
// 	if (args.delete_rate && args.genome_size){
// 		const unsigned short DEL_MAX_LENGTH = 5;

// 		if (!args.quiet)
// 			cout<<"Generating Deletes..."<<endl;
		
// 		// generate random deletion indeces 
// 		unsigned long num_dels = args.genome_size * args.delete_rate;
// 		vector<indel_tuple> deletes(num_dels);
// 		if (!args.quiet)
// 			cout<<"\tdeleting random slices..."<<endl;
// 		for (unsigned long i = 0; i < deletes.size(); i++){
// 			unsigned int chrom = rand() % genome.size();
// 			unsigned long index = random() * (genome[chrom].size() - 1);
// 			unsigned short length = (rand() % DEL_MAX_LENGTH) + 1;
// 			const string alleles = get_slice(
// 				chrom, 
// 				index, 
// 				index + length,
// 				genome
// 			);
// 			deletes[i] = make_tuple(chrom, index, alleles);
// 		}
// 		sort(deletes.begin(), deletes.end(), indel_tuple_comparator);

// 		// replace deletion slices with dummy sequences
// 		enum {CHROM, INDEX, ALLELES};
// 		vector<indel_tuple>::reverse_iterator i;
// 		for (i = deletes.rbegin(); i < deletes.rend(); i++){
// 			const unsigned int chrom = get<CHROM>((*i));
// 			const unsigned long index = get<INDEX>((*i));
// 			const unsigned short length = get<ALLELES>((*i)).size();
// 			vector<char>::iterator start = genome[chrom].begin() + index;
// 			vector<char>::iterator end = genome[chrom].begin() + index + length;
// 			genome[chrom].erase(start, end);
// 			const string dummy = string(length, DEL_MARKER);
// 			genome[chrom].insert(start, dummy.begin(), dummy.end());
// 		}

// 		// write deletes to file
// 		ofstream outfile((char*)ans_file_name.c_str(), ios::app);
// 		if (outfile.is_open()){
// 			if (!args.quiet)
// 					cout<<"\twriting deletes..."<<endl;
// 			outfile<<">DELETE\n";
// 			vector<indel_tuple>::iterator i;
// 			for (i = deletes.begin(); i < deletes.end(); i++){
// 				const unsigned int chrom = get<CHROM>((*i));
// 				const unsigned long index = get<INDEX>((*i));
// 				const string alleles = get<ALLELES>((*i));
// 				outfile<<chrom<<','<<alleles<<','<<index<<'\n';
// 			}
// 		}
// 		outfile.close();
// 	}
// }

// 			  index, 		 old   new
typedef tuple<unsigned long, char, char> snp_tuple; 

bool snp_tuple_comparator(snp_tuple a, snp_tuple b){
	if (get<0>(a) < get<0>(b))
		return true;
	return false;
}

void generate_snps(vector<vector<char>>& genome){
	if (args.snp_rate && args.genome_size){
		ofstream outfile((char*)ans_file_name.c_str(), ios::app);
		if (outfile.is_open()){
			if (!args.quiet)
				cout<<"Generating SNPs..."<<endl;
			outfile<<">SNP\n";
			for (unsigned int chrom = 0; chrom < genome.size(); chrom++){
				unsigned long num_snps = genome[chrom].size() * args.snp_rate;
				if (!args.quiet)
					cout<<"\tchrom "<<chrom<<endl;
				vector<snp_tuple> random_snps(num_snps);
				for (unsigned long j = 0; j < random_snps.size(); j++){
					while (1){
						unsigned long index = random() * (genome[chrom].size() - 1);
						// try until a valid genome position is found
						try{
							char snp = random_snp(genome[chrom][index]);
							random_snps[j] = make_tuple(
								index,
								genome[chrom][index],
								snp
							);
							genome[chrom][index] = snp;
							break;
						}
						catch (invalid_argument &e){ //thrown by random_snp()
							continue; 
						}
					}
				}
				sort(random_snps.begin(), random_snps.end(), snp_tuple_comparator);
				
				enum {INDEX, OLD, NEW};
				vector<snp_tuple>::iterator i;
				for (i = random_snps.begin(); i < random_snps.end(); i++){
					const char old = get<OLD>((*i));
					const char snp = get<NEW>((*i));
					const unsigned long index = get<INDEX>((*i));
					outfile<<chrom<<','<<old<<','<<snp<<','<<index<<'\n';
				}
			}
		}
		outfile.close();
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
	unsigned long start, unsigned long end, const vector<vector<char>>& genome)
{
	string slice = "";
	for (unsigned long i = start; i < end; i++){
		slice += genome[chromosome][i];
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

bool overlaps(const unsigned long start1, const unsigned int length1,
	const unsigned long start2, const unsigned int length2)
{
	unsigned long end1 = start1 + length1;
	unsigned long end2 = start2 + length2;
	return end1 > start2 && end2 > start1; 
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
  	if (!args.quiet)
  		cout<<"DONE"<<endl;

	return 0;
}