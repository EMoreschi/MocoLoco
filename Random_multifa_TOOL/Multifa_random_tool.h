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
string JASPAR_FILE;
unsigned int n_oligo = 80;
unsigned int position = 250;
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
//			print_debug_matrix();
			oligo_creation();
//			print_oligo_vector();
			matrix_size = matrix[0].size();
		}

		vector<string> oligo_vector;
		int matrix_size;
};

class multifasta_class{

	private:
		
		map<int,string> multifasta_map;
		map<int,string> multifasta_map_implanted;
		vector<string> oligo_vector;
		vector<string> headers;
		vector<string> sequences;
		int matrix_size;

		void multifasta_map_creation();
		void multifasta_outfile(map<int,string>, string);
		void implanting_oligo(int);

	public:

		multifasta_class(vector<string> oligo_v, int matrix_size){
			
			oligo_vector = oligo_v;
			multifasta_map_creation();
			multifasta_map_implanted = multifasta_map;
			multifasta_outfile(multifasta_map, "random_multifa.fasta");
			implanting_oligo(matrix_size);
			multifasta_outfile(multifasta_map_implanted, "random_multifa_implanted.fasta");
		}
		
};


void command_line_parser(int, char **);
void display_help(); 						//Display help function
int random_number(int, int);
char from_n_to_base(int);
bool is_file_exist(string, string);
