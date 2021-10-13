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
#include "./../TwoBit/twobit.c"
#include "./../TwoBit/twobit.h"
using namespace std;

int half_length = 150;
unsigned int overhead = 25;

class bed_class {
	private:
		string chr_coord;
		unsigned int start_coord;
		unsigned int end_coord;
		bool flag;
		string sequence;

		void read_line(string);
		void flag_control();
	public:
		bed_class (string line, TwoBit* tb, unsigned int n_line){
			read_line(line);
			flag_control();
            		centering_function();
			extract_seq(tb, n_line);
		}
		string return_sequence();
		void centering_function();
		void extract_seq(TwoBit*, unsigned int);
		unsigned int return_start_coord();
		unsigned int return_end_coord();
		string return_chr_coord();
};
/*
class coordinator_class{                     //Coordinator class to connect Matrix to Bed and Oligos_vector

    private:

        void GEP_creation(vector<bed_class>&);

    public:

        coordinator_class(){

            //GEP (vector of bed class) creation. An empty GEP vector is passed by reference to be filled and saved in this class
            GEP_creation(GEP);

        }
        
        vector<bed_class> GEP;
        
};
*/
//Multifasta class: Creating a set of random seuences following input parameters
class multifasta_class{

	private:
		
		vector<string> headers;
		vector<string> sequences;

		void multifasta_map_creation();
		void multifasta_outfile(map<unsigned int,string>, string);

	public:

		//Constructor of the class --> passed length, number of seq and current cycle number (for the file name)
		multifasta_class(unsigned int length, unsigned int n_seq, unsigned int i){
			
			//Function to create the set of random sequences
			multifasta_map_creation();			

			//Function to handle the output file --> Printing the sequences on it (or them for more cycles)
			multifasta_outfile(multifasta_map, "random_multifa_"+to_string(i+1)+".fasta");
		}

		map<unsigned int,string> multifasta_map;
		
};

class implanting_class{

	private:
		map<unsigned int,string> multifasta_map_implanted;
		vector<unsigned int> unique_rnd;
		vector<string> oligo_vector;

		void implanting_oligo(map<vector<unsigned int>, vector<vector<unsigned int>>>);
		void multifasta_outfile_2(map<unsigned int,string>, string);
		void unique_random_generator();
		void oligo_creation(map<vector<unsigned int>, vector<vector<unsigned int>>>::iterator, int);
		void print_debug_matrix(vector<vector<unsigned int>>);

	public:

		implanting_class(map<vector<unsigned int>, vector<vector<unsigned int>>> jaspar_map, map<unsigned int,string> multifasta_map, unsigned int i){
		
			//Copying the random multifasta sequences into another map (where the iplanting are going to be performed)
			multifasta_map_implanted = multifasta_map;	

			//Function to implant oligos (generated following input parameters) into the sequences
			implanting_oligo(jaspar_map);		

			//Printing sequences + implants on an Output file (where i+1 is the number of current cycle)
			multifasta_outfile_2(multifasta_map_implanted, "random_multifa_implanted"+to_string(i+1)+".fasta");
		}
};

void generic_vector_creation(string, vector<unsigned int>&);
void command_line_parser(int, char **);
void pathway_multifasta();
void pathway_bed();
void read_input();
void display_help(); 
unsigned int random_number(unsigned int, unsigned int);
string reverse_complement(string);
char from_n_to_base(unsigned int);
void wobble_vector_creation(string);
void check_oligo_number();
void check_position_vector();
void check_frequence();
void check_wobble();
void check_overlapping(map<vector<unsigned int>, vector<vector<unsigned int>>>);
void check_input();
void freq_strand_plus_matrix(map<vector<unsigned int>, vector<vector<unsigned int>>>&);
void filling_jaspar_map(map<vector<unsigned int>, vector<vector<unsigned int>>>&, string, string);
void ordering_matrix_names();
vector<vector<unsigned int>> read_JASPAR(string,string&,string&);
void check_exceeding(map<vector<unsigned int>,vector<vector<unsigned int>>>);
void implanting_cycle(unsigned int);
void find_oligo_number();
bool is_file_exist(string, string);
void bed_class_creation(vector<bed_class>&);
void print_debug_bed(vector<bed_class>);
