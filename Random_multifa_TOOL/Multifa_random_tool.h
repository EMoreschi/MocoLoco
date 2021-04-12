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
unsigned int length = 500;				//length of sequences
unsigned int n_seq = 200;				//number of sequences
vector<string> JASPAR_FILE_vector;			//Vector of Jaspar_file strings
string JASPAR_F;					
string oligo_perc = "80";				
vector<unsigned int> n_oligo_vector;			//vector of oligo_number we are implanting (calculated from percentages in input)
string position = "10";
string wobble = "0";
vector<unsigned int> position_vector;			//vector of implanting position (taken as input)
vector<unsigned int> wobble_vector;			//vector of wobble interval (taken as input)
map<unsigned int,string> position_jaspar_map;		//map of (implant_position | Jaspar_file_string)
unsigned int cycles = 1;

class matrix_class {

	private: 

		string matrix_name;
		string tf_name;
		vector<vector<unsigned int>> matrix;
		vector<vector<unsigned int>> matrix_sum;		

		void read_JASPAR(string);
		void oligo_creation(unsigned int);
		void check_oligo_number();

	public:
		matrix_class(string Jaspar_string, unsigned int p){
			
			check_oligo_number();				//checking that oligo number is < n_seq
			read_JASPAR(Jaspar_string);			//reading jaspar file and saving in matrix variable
			oligo_creation(p);				//creating a vector of random oligos coming from matrix frequences -> the size of oligos vector follows the input n_oligo_vector values for each matrix
			matrix_size = matrix[0].size();			//calculating the matrix size -> we need for implanting
		}

		unsigned int n_oligo;
		vector<string> oligo_vector;
		unsigned int matrix_size;
		void print_debug_matrix();
};

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

		void implanting_oligo(vector<matrix_class>);
		void multifasta_outfile_2(map<unsigned int,string>, string);
		void unique_random_generator();

	public:

		implanting_class(vector<matrix_class> MATRIX_VECTOR, map<unsigned int,string> multifasta_map, unsigned int i){
		
			multifasta_map_implanted = multifasta_map;	
			implanting_oligo(MATRIX_VECTOR);		
			multifasta_outfile_2(multifasta_map_implanted, "random_multifa_implanted"+to_string(i+1)+".fasta");
		}
};

void command_line_parser(int, char **);
void display_help(); 
unsigned int random_number(unsigned int, unsigned int);
string reverse_complement(string);
char from_n_to_base(unsigned int);
void position_vector_creation(string);
void wobble_vector_creation(string);
void n_oligo_vector_creation(string);
void check_overlapping(vector<matrix_class>);
void check_input();
void check_positions(vector<matrix_class>);
void check_jaspar_exist(unsigned int);
void find_oligo_number();
bool is_file_exist(string, string);
