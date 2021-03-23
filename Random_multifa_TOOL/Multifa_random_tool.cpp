#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iterator>
#include <list>
#include <vector>
#include <algorithm>
#include <map>
#include <typeinfo>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <sys/time.h>
#include <sstream>
#include <numeric>
#include <getopt.h>
#include <unordered_map>

using namespace std;
int length = 500;
int n_seq = 200;

void command_line_parser(int, char **);
void display_help(); 						//Display help function
char random_number();
char from_n_to_base(int);

int main(int argc, char *argv[]){

	command_line_parser(argc, argv);					//Parser function called to handle aguments
	string sequence;
	ofstream outfile;
	outfile.open("random_multifasta.fasta");

	for(int j=0; j<n_seq; j++){
		
		sequence.clear();
		
		for(int i=0; i<length; i++){

			char base = random_number();
			sequence = sequence + base;
		}

		outfile << ">random sequence " << j << " containing " << length << " bases." << endl;
		outfile << sequence << endl;
		outfile << endl;
	}
}

char random_number(){

	const gsl_rng_type * T;
	gsl_rng *r;

	gsl_rng_env_setup();
	struct timeval tv;		//Generate a seed depending on time to get every time different random number
	gettimeofday(&tv, 0);
	unsigned long mySeed = tv.tv_sec + tv.tv_usec;

	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, mySeed);
		
	int r_number = gsl_rng_uniform_int(r,4);
	
	char base = from_n_to_base(r_number);

	gsl_rng_free(r);

	return 	base;
}

char from_n_to_base (int n){

	char base;
			switch(n){

				case 0:

					base = 'A';
					break;

				case 1:

					base = 'T';
					break;

				case 2:

					base = 'C';
					break;

				case 3:

					base = 'G';
					break;

				default:
					
					cout << "Generation of random sequence failed. Please check the random generation code.!" << endl;
					exit(1);
					break;
			}
			return base;
}


void command_line_parser(int argc, char** argv){
	
	const char* const short_opts = "h:l:n";

	//Specifying the expected options
	const option long_opts[] ={
		{"help",      no_argument, nullptr,  'h' },
		{"length",      required_argument, nullptr,  'l' },
		{"nseq",   required_argument, nullptr,  'n' },
	};

	while (true)
	{
		const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

		if (-1 == opt)
			break;


		switch (opt) {
			case 'h' : display_help();
				   break;
			case 'l' : length = stoi(optarg); 
				   break;
			case 'n' : n_seq = stoi(optarg); 
				   break;
			case '?': // Unrecognized option
			default:
				   display_help();
				   break;
		}
	}
}

void display_help() 						//Display help function
{
	cerr << "\n --help || -h show this message" << endl;
	cerr << "\n --length || -l <number> to insert the length of Multifasta sequences (DEFAULT: 500)" << endl;
	cerr << "\n --nseq || -n <number> to insert the number of Multifasta sequences (DAFAULT: 200" << endl;
	cerr << endl;

	exit(EXIT_SUCCESS);
}
