#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <dirent.h>
#include <cmath>
#include <valarray>
#include <vector>
#include <algorithm>
#include <typeinfo> 
#include <sstream>
#include <numeric>
#include "./TwoBit/twobit.h"
#include "./TwoBit/twobit.c"

using namespace std;

string BED_FILE; 		//initializing const char variarible for Bed_file input reading
int parameter = 150; 		//default parameter 150
string TWOBIT_FILE;	//initializing const char variable for Twobit_file input reading
string JASPAR_FILE;
const int overhead = 25;
const double pseudoc = 0.01;

class genomic_position { //creation public class of genomic_position type        

	private:	//field definition

		string chr_coord;
		int start_coord;
		int end_coord;
		bool flag;
		string sequence;

		void centering_function(int, int, int, const int);
		void read_line(string);
		void flag_control(int, int);
		void extract_seq(TwoBit*, int);

	public:				//field definition
		genomic_position(){	//default constructor

			chr_coord = "";
			start_coord = 0;
			end_coord = 0;
			flag = 0;	

		}

		genomic_position(int p, string line, TwoBit* tb,int n_line){

			read_line(line);					//reading bed line
			flag_control(start_coord,end_coord);			//controlling coordinates
			centering_function(start_coord, end_coord, p, overhead);		//centering the coordinates
			extract_seq(tb, n_line);				//extracting the sequence

		}
		void print_debug_GEP(genomic_position);
		string return_sequence(genomic_position);

};

class matrix_class {

	private: //field definition

		string matrix_name;
		string tf;
		vector<vector<double>> matrix;
		vector<vector<double>> norm_matrix;
		vector<vector<double>> inverse_complement_matrix;
		vector<vector<double>> matrix_log;
		vector<double> local_mins;	
		vector<double> local_maxes;
		double global_min;
		double global_max;
		vector<double> col_sum;		


		void matrix_normalization_pseudoc(vector<vector<double>>, double);
		void matrix_normalization(vector<vector<double>>);
		void matrix_logarithmic(vector<vector<double>>);
		void read_JASPAR(string);
		void inverse_matrix(vector<vector<double>>);
		void find_minmax(vector<vector<double>>);
		vector<double> find_col_sum(vector<vector<double>>);

		

	public:
		matrix_class(string JASPAR_FILE){

			read_JASPAR(JASPAR_FILE);
			matrix_normalization_pseudoc(matrix, pseudoc);			//Calling matrix normalization function
			matrix_normalization(norm_matrix);
			matrix_logarithmic(norm_matrix);
			inverse_matrix(norm_matrix);
			find_minmax(matrix_log);		
		}
		void print_debug_matrix(vector<vector<double>>, string);
                void shifting(string seq, int p, int length, vector<double>&);
		vector<vector<double>> return_matrix(int);
		vector<vector<double>> return_norm_matrix(int);
		vector<vector<double>> return_inverse_matrix(int);
		vector<vector<double>> return_log_matrix(int);

};

void GEP_creation(string, string, vector<genomic_position>&);
void command_line_parser(int, char **);
void display_help();
bool exist_test0(string);
void is_file_exist(string fileName, string);
bool isDir(string);
