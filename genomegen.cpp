#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <tuple>
#include <iomanip>

#include "genomegen.hpp"
#include "tclap/CmdLine.h"

using namespace std;
using namespace genome_gen;


namespace genome_gen {

static const string REF_PRE = "ref_";
static const string PRIV_PRE = "donor_";
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
		// TCLAP::UnlabeledValueArg<unsigned int> num_chroms(
		// 	"num-chromosomes",
		// 	"Num-Chromosomes: The number of chromosomes to generate for the "
		// 	"genome",
		// 	true,
		// 	1,
		// 	"int"
		// );
		TCLAP::UnlabeledValueArg<double> chrom_size(
			"genome-size",
			"Genome-Size: The size of the genome, by default in "
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

		// bind args to parser and parse
		cmd.add(quiet);
		cmd.add(genome_id);
		// cmd.add(num_chroms);
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
		args.quiet 							= quiet.getValue();
		args.genome_id 					= genome_id.getValue();
		// args.num_chroms 				= num_chroms.getValue();
		args.num_chroms 				= 1;
		args.chrom_size 				= factor * chrom_size.getValue();
		args.scale 							= scale.getValue();
		args.coverage 					= coverage.getValue();
		args.read_length 				= read_length.getValue();
		args.read_gap 					= read_gap.getValue();
		args.read_error_rate 		= read_error_rate.getValue();
		args.read_garbage_rate 	= read_garbage_rate.getValue();
		args.genome_size 				= args.chrom_size * args.num_chroms;
		args.snp_rate 					= snp_rate.getValue();
		args.indel_rate 				= indel_rate.getValue();

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
	// cout<<"num-chroms:\t"<<args.num_chroms<<endl;
	// switch(args.scale){
	// 	case 'k': 
	// 	cout<<"chrom-size:\t"<<args.chrom_size / 1000.0<<"-thousand bp"
	// 	<<endl; break;
	// 	case 'm': 
	// 	cout<<"chrom-size:\t"<<args.chrom_size / 1000000.0<<"-million bp"
	// 	<<endl; break;
	// 	case 'b': 
	// 	cout<<"chrom-size:\t"<<args.chrom_size / 1000000000.0<<"-billion bp"
	// 	<<endl; break;
	// }
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
	if (args.genome_size){
		ofstream outfile((char*)ref_file_name.c_str());
		if (outfile.is_open()){
			if (!args.quiet)
				cout<<"Generating reference genome..."<<endl;
			outfile<<">" + args.genome_id;
			for (unsigned int i = 0; i < args.num_chroms; i++){
				// string chrom_num = static_cast<ostringstream*>( 
				// 	&(ostringstream() << i) )->str();
				// outfile<<"\n>chromosome_" + chrom_num;
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

void write_donor_genome(vector<vector<char>>& genome){
	if (args.genome_size){
		ofstream outfile((char*)priv_file_name.c_str());
		if (outfile.is_open()){
			if (!args.quiet)
				cout<<"Writing donor genome..."<<endl;
			outfile<<">" + args.genome_id;
			for (unsigned int chrom = 0; chrom < genome.size(); chrom++){
				// string chrom_num = static_cast<ostringstream*>( 
				// 	&(ostringstream() << chrom) )->str();
				// outfile<<"\n>chromosome_" + chrom_num;
				unsigned long line_pos = 0;
				for (unsigned long j = 0; j < genome[chrom].size(); j++){
					if (line_pos % LINE_WIDTH == 0)
						outfile<<'\n';
					outfile<<genome[chrom][j];
					line_pos++;
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
				unsigned long index1 = (unsigned long)(random() * 
					(genome[chrom_num].size() - gap - 2*args.read_length));
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
	generate_snps(genome);
}

// 			  chrom_num,   	index, 		   alleles
typedef tuple<unsigned int, unsigned long, string> indel_tuple; 

bool indel_tuple_compare(indel_tuple a, indel_tuple b){
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
		const unsigned long num_inserts = num_indels / 2;
		const unsigned long num_deletes = num_indels / 2;
		vector<indel_tuple> inserts(num_inserts);
		vector<indel_tuple> deletes(num_deletes);
		if (!args.quiet)
			cout<<"\tmutating genome...\n";
		for (unsigned long i = 0; i < num_indels / 2; i++){
			// cout<<'\r';
			const unsigned int length = (rand() % INDEL_MAX_LENGTH) + 1;
			unsigned int chrom = rand() % genome.size();
			// delete slice and replace with dummies
			const unsigned long del_index = 
				random() * (genome[chrom].size() - length);
			const string del_alleles = get_slice(
				chrom, 
				del_index, 
				del_index + length,
				genome
			);
			vector<char>::iterator del_start = 
				genome[chrom].begin() + del_index;
			vector<char>::iterator del_end = 
				genome[chrom].begin() + del_index + length;
			genome[chrom].erase(del_start, del_end);
			const string dummy = string(length, DEL_MARKER);
			genome[chrom].insert(del_start, dummy.begin(), dummy.end());
			// insert
			unsigned long ins_index;
			do{
				ins_index = random() * (genome[chrom].size() - 1);
			} while (ins_index >= del_index && del_index + length >= ins_index);
			const string ins_alleles = random_alleles(length);
			genome[chrom].insert(
				genome[chrom].begin() + ins_index, 
				ins_alleles.begin(), 
				ins_alleles.end()
			);
			// remove dummy placeholders
			vector<char>::iterator iter = genome[chrom].begin() + del_index;
			if (ins_index < del_index) iter += length;
			else if (ins_index > del_index + length);
			else{
				cout<<"Iter: "<<i<<endl;
				cout<<"Insert Index: "<<ins_index<<endl;
				cout<<"Delete range: ["
					<<del_index<<":"<<del_index + length<<"]"<<endl;
				throw runtime_error("insert index overlaps delete range");
			}
			vector<char>::iterator stop = iter + length;
			for (; iter < stop; iter++){
				if ((*iter) == DEL_MARKER)
					iter = genome[chrom].erase(iter) - 1;
			}
			inserts[i] = make_tuple(chrom, ins_index, ins_alleles);
			deletes[i] = make_tuple(chrom, del_index, del_alleles);
			cout<<setw(4)<<setprecision(3)<<
				"\r\t"<<(200 * (float)(i / (float)num_indels))<<"%"; 
		}
		cout<<endl;
		sort(inserts.begin(), inserts.end(), indel_tuple_compare);
		sort(deletes.begin(), deletes.end(), indel_tuple_compare);

		ofstream outfile((char*)ans_file_name.c_str(), ios::app);
		// write inserts to file
		if (!args.quiet)
			cout<<"\twriting insertions to file..."<<endl;
		outfile<<">INS\n";
		enum {CHROM, INDEX, ALLELES};
		vector<indel_tuple>::iterator i;
		for (i = inserts.begin(); i < inserts.end(); i++){
			const unsigned int chrom = get<CHROM>((*i));
			const unsigned long index = get<INDEX>((*i));
			const string alleles = get<ALLELES>((*i));
			outfile/*<<chrom<<','*/<<alleles<<','<<index<<'\n';
		}

		// write deletes to file
		if (!args.quiet)
			cout<<"\twriting deletions to file..."<<endl;
		outfile<<">DEL\n";
		vector<indel_tuple>::iterator d;
		for (d = deletes.begin(); d < deletes.end(); d++){
			const unsigned int chrom = get<CHROM>((*d));
			const unsigned long index = get<INDEX>((*d));
			const string alleles = get<ALLELES>((*d));
			outfile/*<<chrom<<','*/<<alleles<<','<<index<<'\n';
		}
		outfile.close();
	}
}

// 			  index, 		 old   new
typedef tuple<unsigned long, char, char> snp_tuple; 

bool snp_tuple_compare(snp_tuple a, snp_tuple b){
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
				// if (!args.quiet)
				// 	cout<<"\tchrom "<<chrom<<endl;
				vector<snp_tuple> snps(num_snps);
				for (unsigned long j = 0; j < snps.size(); j++){
					unsigned long index = random() * (genome[chrom].size() - 1);
					char snp = random_snp(genome[chrom][index]);
					snps[j] = make_tuple(index, genome[chrom][index], snp);
					genome[chrom][index] = snp;
				}
				sort(snps.begin(), snps.end(), snp_tuple_compare);
				
				enum {INDEX, OLD, NEW};
				vector<snp_tuple>::iterator i;
				for (i = snps.begin(); i < snps.end(); i++){
					const char old = get<OLD>((*i));
					const char snp = get<NEW>((*i));
					const unsigned long index = get<INDEX>((*i));
					outfile/*<<chrom<<','*/<<old<<','<<snp<<','<<index<<'\n';
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
	write_donor_genome(genome);
	write_reads(genome);
	if (!args.quiet)
		cout<<"DONE"<<endl;

	return 0;
}