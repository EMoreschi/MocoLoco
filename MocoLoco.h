#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <dirent.h>
#include <vector>
#include <typeinfo> 
#include <sstream>
#include "./TwoBit/twobit.h"
#include "./TwoBit/twobit.c"

using namespace std;

string BED_FILE; 		//initializing const char variarible for Bed_file input reading
int parameter = 150; 		//default parameter 150
string TWOBIT_FILE;	//initializing const char variable for Twobit_file input reading
string JASPAR_FILE;
const int overhead = 25;

class genomic_position { //creation public class of genomic_position type        

	private:	//field definition

		string chr_coord;
		int start_coord;
		int end_coord;
		bool flag;
		string sequence;

		void centering_function(int start, int end, int p, const int overhead);
		void read_line(string line);
		void flag_control(int start, int end);
		void extract_seq(TwoBit* tb, int n_line);

	public:
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

};

class matrix_class {

	private: //field definition

		string matrix_name;
		string tf;
		vector<vector<double>> matrix;
		void read_JASPAR(string JASPAR_FILE);
	public:
		matrix_class(string JASPAR_FILE){

			read_JASPAR(JASPAR_FILE);

		}
		void print_debug_matrix(matrix_class);
};

void GEP_creation(string, string, vector<genomic_position>&);
void command_line_parser(int, char **);
void display_help();
bool exist_test0(string);
bool is_file_exist(string fileName);
bool isDir(string);
