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
unsigned int length = 500;
unsigned int n_seq = 200;
vector<string> JASPAR_FILE_vector;
string JASPAR_F;
unsigned int n_oligo = 80;
string position = "250";
vector<unsigned int> position_vector;
map<unsigned int,string> position_jaspar_map;
bool flag_JASPAR = 0;

class matrix_class {

	private: 

		string matrix_name;
		string tf_name;
		vector<vector<unsigned int>> matrix;
		vector<vector<unsigned int>> matrix_sum;		

		void read_JASPAR(string);
		void oligo_creation();
		void check_oligo_number();

	public:
		matrix_class(string Jaspar_string){
			
			check_oligo_number();
			read_JASPAR(Jaspar_string);
			oligo_creation();
			matrix_size = matrix[0].size();
		}

		vector<string> oligo_vector;
		unsigned int matrix_size;
		void print_debug_matrix();
		void print_oligo_vector();
};

class multifasta_class{

	private:
		
		vector<string> headers;
		vector<string> sequences;

		void multifasta_map_creation();
		void multifasta_outfile(map<unsigned int,string>, string);

	public:

		multifasta_class(unsigned int length, unsigned int n_seq){
			
			multifasta_map_creation();
			multifasta_outfile(multifasta_map, "random_multifa.fasta");
		}
		map<unsigned int,string> multifasta_map;
		
};

class implanting_class{

	private:
		map<unsigned int,string> multifasta_map_implanted;

		void implanting_oligo(vector<matrix_class>);
		void multifasta_outfile_2(map<unsigned int,string>, string);

	public:

		implanting_class(vector<matrix_class> MATRIX_VECTOR, map<unsigned int,string> multifasta_map){
		
			multifasta_map_implanted = multifasta_map;
			implanting_oligo(MATRIX_VECTOR);
			multifasta_outfile_2(multifasta_map_implanted, "random_multifa_implanted.fasta");
		}
};

void command_line_parser(int, char **);
void display_help(); 
unsigned int random_number(unsigned int, unsigned int);
char from_n_to_base(unsigned int);
void position_vector_creation(string);
void check_overlapping(vector<matrix_class>);
void check_input();
void check_positions(vector<matrix_class>);
void print_debug_matrixclass(vector<matrix_class>);
bool is_file_exist(string, string);
