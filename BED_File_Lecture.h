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
				
		string chr;
		unsigned int start;
		unsigned int end;
		bool flag;
		
		

		//public:			//definition public field
		

		//CONSTRUCTORS:
		genomic_position(){	//default constructor
		
		chr = "";
		start = 0;
		end = 0;
		flag = 0;	
		} 

		genomic_position(string chr_coord, unsigned int start_coord, unsigned int end_coord) {  //constructor given member values
			
			
			chr = chr_coord;
			start = start_coord;
			end = end_coord;
			flag = flag_control(start_coord, end_coord);
		}

		//FUNCTIONS:
		string get_char(){	//print char function
		return chr;
		}

		int get_start(){	//print start coordinate function
		return start;
		}

		int get_end(){		//print end coordinate funcion
		return end;
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
