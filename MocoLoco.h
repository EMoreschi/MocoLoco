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
#include <unistd.h>
#include "./TwoBit/twobit.h"
#include "./TwoBit/twobit.c"

using namespace std;

string BED_FILE; 		//initializing const char variarible for Bed_file input reading
int parameter = 150; 		//default parameter 150
string TWOBIT_FILE;	//initializing const char variable for Twobit_file input reading
string JASPAR_FILE;
const int overhead = 25;
const double pseudoc = 0.01;

class bed_class { //creation public class of bed_class type        

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
		bed_class(){	//default constructor

			chr_coord = "";
			start_coord = 0;
			end_coord = 0;
			flag = 0;	

		}

		bed_class(int p, string line, TwoBit* tb,int n_line){

			read_line(line);					//reading bed line
			flag_control(start_coord,end_coord);			//controlling coordinates
			centering_function(start_coord, end_coord, p, overhead);		//centering the coordinates
			extract_seq(tb, n_line);				//extracting the sequence

		}
		void print_debug_GEP(bed_class);
		string return_sequence(bed_class);

};

class matrix_class {

	private: //field definition

		string matrix_name;
		string tf_name;
		vector<vector<double>> matrix;
		vector<vector<double>> norm_matrix;
		vector<vector<double>> inverse_complement_matrix;
		vector<vector<double>> matrix_log;
		vector<double> col_sum;		


		void matrix_normalization_pseudoc(vector<vector<double>>, double);
		void matrix_normalization(vector<vector<double>>);
		void matrix_logarithmic(vector<vector<double>>);
		void read_JASPAR(string);
		void inverse_matrix(vector<vector<double>>);
		vector<double> find_col_sum(vector<vector<double>>);

		

	public:
		matrix_class(string JASPAR_FILE){

			read_JASPAR(JASPAR_FILE);
			matrix_normalization_pseudoc(matrix, pseudoc);			//Calling matrix normalization function
			matrix_normalization(norm_matrix);
			matrix_logarithmic(norm_matrix);
			inverse_matrix(norm_matrix);
		}
		void print_debug_matrix(vector<vector<double>>, string);
                void shifting(string seq, int p, int length, vector<double>&);
		vector<vector<double>> return_matrix(int);
		vector<vector<double>> return_norm_matrix(int);
		vector<vector<double>> return_inverse_matrix(int);
		vector<vector<double>> return_log_matrix(int);

};


class oligo_class{

	private:
		
			
		vector<double> oligo_scores;	
		vector<double> o_matrix_mins;	
		vector<double> o_matrix_maxes;
		double min_possible_score;
		double max_possible_score;
		double best_score;
		string best_oligo;

		void find_minmax(vector<vector<double>>);
		void find_best_score(vector<double>);

	public:

		oligo_class(vector<vector<double>> matrix, string sequence){
		
			find_minmax(matrix);		
			shifting(matrix, sequence, 0);
			find_best_score(oligo_scores);
		}

		void shifting(vector<vector<double>>, string, int);
		vector<double> return_oligo_scores(int);
		double return_best_score(int);


};


void GEP_creation(string, string, vector<bed_class>&);
void command_line_parser(int, char **);
void display_help();
bool exist_test0(string);
bool is_file_exist(string fileName, string buf);
