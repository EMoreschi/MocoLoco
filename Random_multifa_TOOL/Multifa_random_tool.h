#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <iterator>
#include <list>
#include <sstream>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream>
#include <getopt.h>
#include <random>
#include <map>

using namespace std;


class multifasta_class{

	private:
		
		vector<string> headers;
		vector<string> sequences;

		void multifasta_map_creation();
		void multifasta_outfile(map<unsigned int,string>, string);

	public:

		multifasta_class(unsigned int length, unsigned int n_seq, unsigned int i){
			
			multifasta_map_creation();			
			multifasta_outfile(multifasta_map, "random_multifa_"+to_string(i+1)+".fasta");
		}
		map<unsigned int,string> multifasta_map;
		
};

class implanting_class{

	private:
		map<unsigned int,string> multifasta_map_implanted;
		vector<unsigned int> unique_rnd;
		vector<string> oligo_vector;

		void implanting_oligo(map<vector<unsigned int>, vector<vector<unsigned int>>>);
		void multifasta_outfile_2(map<unsigned int,string>, string);
		void unique_random_generator();
		void oligo_creation(map<vector<unsigned int>, vector<vector<unsigned int>>>::iterator, int);
		void print_debug_matrix(vector<vector<unsigned int>>);

	public:

		implanting_class(map<vector<unsigned int>, vector<vector<unsigned int>>> jaspar_map, map<unsigned int,string> multifasta_map, unsigned int i){
		
			multifasta_map_implanted = multifasta_map;	
			implanting_oligo(jaspar_map);		
			multifasta_outfile_2(multifasta_map_implanted, "random_multifa_implanted"+to_string(i+1)+".fasta");
		}
};

void generic_vector_creation(string, vector<unsigned int>&);

void command_line_parser(int, char **);
void read_input();
void display_help(); 
unsigned int random_number(unsigned int, unsigned int);
string reverse_complement(string);
char from_n_to_base(unsigned int);
void wobble_vector_creation(string);
void check_oligo_number();
void check_position_vector();
void check_wobble();
void check_overlapping(map<vector<unsigned int>, vector<vector<unsigned int>>>);
void check_input();
void freq_strand_plus_matrix(map<vector<unsigned int>, vector<vector<unsigned int>>>&);
void filling_jaspar_map(map<vector<unsigned int>, vector<vector<unsigned int>>>&, string, string);
void ordering_matrix_names();
vector<vector<unsigned int>> read_JASPAR(string,string&,string&);
void check_exceeding(map<vector<unsigned int>,vector<vector<unsigned int>>>);
void implanting_cycle(unsigned int);
void find_oligo_number();
bool is_file_exist(string, string);
