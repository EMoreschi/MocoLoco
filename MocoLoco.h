#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <typeinfo> 
#include <sstream>
#include "./TwoBit/twobit.h"
#include "./TwoBit/twobit.c"

using namespace std;

class genomic_position { //creation public class of genomic_position type        

	public:	//definition private field

		string chr_coord;
		int start_coord;
		int end_coord;
		bool flag;
		string sequence;
		void centering_function(int start, int end, int p);
		void read_line(string line);
		void flag_control(int start, int end);
                void extract_seq(TwoBit* tb, int n_line);


		//public:			//definition public field


		//CONSTRUCTORSstart_coord :
		genomic_position(){	//default constructor

			chr_coord = "";
			start_coord = 0;
			end_coord = 0;
			flag = 0;	

		}

		genomic_position(int p, string line, TwoBit* tb,int n_line){


			read_line(line);
			flag_control(start_coord,end_coord);
			centering_function(start_coord, end_coord, p); //function to center the coordinates
			extract_seq(tb, n_line);
				


		}

};

vector<genomic_position> GEP_creation(const char*, const char*);
void stamp_debug(vector<genomic_position> gep);
void command_line_parser(int, char **);
void display_help();
bool exist_test0(const char*);
bool is_file_exist(const char *fileName);
