#include <cstdio>
#include <sys/stat.h>
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
int length = 500;
int n_seq = 200;
string JASPAR_FILE;


class matrix_class {

	private: //field definition

		string matrix_name;
		string tf_name;
		vector<vector<double>> matrix;
		vector<vector<int>> matrix_sum;		


		void read_JASPAR(string);
		void find_col_sum();
		void print_debug_matrix();

	public:
		matrix_class(string JASPAR_FILE){

			read_JASPAR(JASPAR_FILE);
			print_debug_matrix();
			find_col_sum();

		}
};


void command_line_parser(int, char **);
void display_help(); 						//Display help function
int random_number(int, int);
char from_n_to_base(int);
bool is_file_exist(string, string);
