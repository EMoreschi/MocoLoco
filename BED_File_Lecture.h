#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <typeinfo> 
#include <sstream>

using namespace std;

string BED_FILE;
int parameter = 150; //default parameter 150

class genomic_position { //creation public class of genomic_position type        

	public:	//definition private field

		string chr_coord;
		unsigned int start_coord;
		unsigned int end_coord;
		bool flag;



		//public:			//definition public field


		//CONSTRUCTORS:
		genomic_position(){	//default constructor

			chr_coord = "";
			start_coord = 0;
			end_coord = 0;
			flag = 0;	
		} 


		bool flag_control(unsigned int st, unsigned int en){ //function which controls that start coordinates are < then end coordinates

			if(st > en){		//if start coordinates are > then end coordinates, flag is setted to 0
				return 0;
			}
			else{ return 1;}
		}	

		bool get_flag(){		//print end coordinate funcion
			return flag;
		}
};

void centering_function (unsigned int*, unsigned int*, int);
void command_line_parser(int, char **);
void display_help();
