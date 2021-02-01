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
const char * BED_FILE;
int parameter = 150; //default parameter 150
const char * TWOBIT_FILE;

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


		//public:			//definition public field


		//CONSTRUCTORSstart_coord :
		genomic_position(){	//default constructor

			chr_coord = "";
			start_coord = 0;
			end_coord = 0;
			flag = 0;	
		
		}

		genomic_position(int p, string line){
		
			read_line(line);
			flag_control(start_coord,end_coord);
			centering_function(start_coord, end_coord, p); //function to center the coordinates


		}
		
};

		void genomic_position::read_line(string line){
			
			istringstream mystream(line);
			mystream >> chr_coord >> start_coord >> end_coord;
	
		}

		void genomic_position::centering_function ( int start,  int end, int p){
	       
			int overhead = 25;
		       	int centro = (start + end)/2;
			start_coord = centro - p -overhead;
			end_coord = centro + p +overhead;
		}


		void genomic_position::flag_control( int start,  int end){ //function which controls that start coordinates are < then end coordinates

			if(start > end || start == end){		//if start coordinates are > or == then end coordinates, flag is setted to 0
				flag = 0;
			}
			else{ flag = 1;}
		}

void GEP_objects_creation(const char*, const char*);
void command_line_parser(int, char **);
void display_help();
bool exist_test0(const char*);
bool is_file_exist(const char *fileName);
