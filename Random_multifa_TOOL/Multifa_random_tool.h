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
string JASPAR_FILE;
unsigned int n_oligo = 80;
string position = "250";
vector<int> position_vector;
bool flag_JASPAR = 0;
class matrix_class {

	private: //field definition

		string matrix_name;
		string tf_name;
		vector<vector<int>> matrix;
		vector<vector<int>> matrix_sum;		

		void read_JASPAR(string);
		void oligo_creation();
		void print_debug_matrix();
		void check_oligo_number();
		void print_oligo_vector();

	public:
		matrix_class(string JASPAR_FILE){
			
			check_oligo_number();
			read_JASPAR(JASPAR_FILE);
			print_debug_matrix();
			oligo_creation();
			print_oligo_vector();
			matrix_size = matrix[0].size();
		}

		vector<string> oligo_vector;
		int matrix_size;
};

class multifasta_class{

	private:
		
		vector<string> headers;
		vector<string> sequences;

		void multifasta_map_creation();
		void multifasta_outfile(map<int,string>, string);

	public:

		multifasta_class(int length, int n_seq){
			
			multifasta_map_creation();
			multifasta_outfile(multifasta_map, "random_multifa.fasta");
		}
		map<int,string> multifasta_map;
		
};

class implanting_class{

	private:
		map<int,string> multifasta_map_implanted;

		void implanting_oligo(vector<matrix_class>);
		void multifasta_outfile_2(map<int,string>, string);

	public:

		implanting_class(vector<matrix_class> MATRIX_VECTOR, map<int,string> multifasta_map){
		
			multifasta_map_implanted = multifasta_map;
			implanting_oligo(MATRIX_VECTOR);
			multifasta_outfile_2(multifasta_map_implanted, "random_multifa_implanted.fasta");
		}
};

void command_line_parser(int, char **);
void display_help(); 						//Display help function
int random_number(int, int);
char from_n_to_base(int);
void position_vector_creation(string);
bool is_file_exist(string, string);
