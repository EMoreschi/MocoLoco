#include "MocoLoco.h"

const char * BED_FILE; 		//initializing const char variarible for Bed_file input reading
int parameter = 150; 		//default parameter 150
const char * TWOBIT_FILE;	//initializing const char variable for Twobit_file input reading
const char * JASPAR_FILE;

int main(int argc, char *argv[]){


	if(argc == 1){             //If arguments number is 1 means that no input file has been inserted - display help
		display_help();
	}

	command_line_parser(argc, argv);					//parser function called to handle aguments
	
	vector<genomic_position> GEP;
	GEP_creation(BED_FILE, TWOBIT_FILE, GEP); 			//function to read BED and 2Bit files and create GEP objects vector
	read_JASPAR(JASPAR_FILE);
	stamp_debug(GEP); 	 						//print vector function (debug only)


}

void genomic_position::read_line(string line){				//Read line function: it takes in input each line from BED file 

	istringstream mystream(line);					//Split the line word by word and extract chromosome coordinates (chr, start, end)
	mystream >> chr_coord >> start_coord >> end_coord;		

}

void genomic_position::centering_function ( int start,  int end, int p){	//Centering function: in takes start and end coordinate and re-sets them -
										//following an input parameter value (overhead added)
	int overhead = 25;
	int centro = (start + end)/2;						
	start_coord = centro - p;			//no overhead for start
	end_coord = centro + p +overhead;		//overhead for end
}


void genomic_position::flag_control( int start,  int end){ 	//Flag control function: start coordinates must be < then end coordinates

	if(start > end || start == end){		//if start coordinates are > or == then end coordinates, flag is setted to 0
		flag = 0;
	}
	else{ flag = 1;}
}

void GEP_creation(const char* Bed_file, const char* Twobit_file, vector<genomic_position> &GEP){		//Function to read BED and 2Bit files and create GEP object vector

	ifstream in(Bed_file); 						//Opening file in lecture mode
	TwoBit * tb;							//Creating a TwoBit* variable called tb
	tb = twobit_open(Twobit_file);					//Opening 2Bit file with twobit_open function and saved in tb 
	string line; 							//defining line string
	int n_line = 0;							//line counter initialization

	while(getline(in,line)){  					//reading input file line by line with getline function

		if(line.empty())   					//if line is empty or commented --> continue
			continue;
		if(line[0]=='#')
			continue;
		
		genomic_position new_class(parameter,line,tb, n_line);  //Called the object constructor passing the Bed line, parameter P, twobit file tb, and the line counter n_line
		GEP.emplace_back(new_class);				//Put the new object in the GEP vector with emplace function

		n_line = n_line + 1;					//pass to next line 

	}
}

void read_JASPAR(const char* file_jaspar){

	ifstream file(file_jaspar);
	string line;
	string testo;
	int riga = 0;
	string word;
	int col = col_number(file_jaspar);
	int number_words = 0;
	double matrix[4][col-3];

	while(getline(file,line)){

		if(line[0]=='>'){

			testo = line;

		}

	else{
		istringstream mystream(line);
	
		for(int i = 0; i <(col-1); i++){

			double a;	
			char inutili;
			
			if(i == 0){
				
				mystream >> inutili;
			}
			else if(i == 1){
				
				mystream >> inutili;
			}
	
			else {
				
				mystream >> a;
				matrix[riga-1][i-2]=a;

			}
		}
	}
	riga = riga + 1;

	
	}

	cout << testo << "\n";

	for(int i = 0; i<4; i++){
		for(int j=0; j<(col-3); j++){

			cout << matrix[i][j]<< " ";
		}
		cout << "\n";
	}
	file.close();
}

int col_number(string file_s){

	ifstream file(file_s);
	string line;
	string word;
	int number_words = 0;
	int colonne = 0;

	while(getline(file,line)){

		if(line[0]=='>')
			continue;
	
		istringstream mystream(line);
		
		while(mystream >> word){
		
			number_words++;
		}
		colonne = number_words;
		number_words = 0;
	}
	file.close();
	return colonne;

}

void stamp_debug( vector<genomic_position> GEP_print){			//Debug function: Print the GEP vector to control the working flow

	for (int i=0; i<GEP_print.size(); ++i){    			// from 0 to GEP vector length

		cout << ">" << GEP_print[i].chr_coord <<":"<< GEP_print[i].start_coord << "-" << GEP_print[i].end_coord << "\n";	//print chr,start,end
		cout << GEP_print[i].sequence<<"\n";											//print DNA sequence

	}

}

void genomic_position::extract_seq(TwoBit* tb, int n_line){			//Extract sequence function: Extract, from Twobit hg38 genome, the DNA sequence with (chr, start, end) coordinates -
										//extracted from Bed line
	if(flag == 1){								//CONTROL: if flag is 1 means that the current line has starting coordinate > end coordinate, so it is correct
		const char* chrom = chr_coord.c_str(); 				//Put in chrom the string of chr_coord
		sequence = twobit_sequence(tb,chrom,start_coord,end_coord-1); 	//Extract the sequence from the object with the twobit_sequence function
	}
	else {		
		cerr << "WARNING: the line " << n_line <<" is omitted because starting coordinates > end coordinates, please check your BED file!" << "\n";
		//if flag is not 1 means that the current line has starting coordinate < end coordinate: PRINT WARNING!		
	}
}

void command_line_parser(int argc, char **argv){

	int control_bed = 0;		
	int control_twobit = 0;
	int control_p = 0;

	for(int i = 1; i < argc; i++){

		string buf = argv[i];

		if(buf == "--help" || buf == "-h"){

			display_help();

		}

		if(buf == "--BED" || buf == "-B"){

			if(i < argc - 1){

				BED_FILE = argv[++i];
				control_bed = 1;

				bool bed_check = is_file_exist(BED_FILE);
				if(bed_check == 0){
					cout << "File BED does not exist, please insert a BED file as input. \n";
					cout << "FATAL ERROR \n";
					exit(EXIT_SUCCESS);
				}
				continue;

			}
		}

		if(buf == "--param" || buf == "-p"){

			if(i < argc - 1){

				parameter = stoi(argv[++i]);
				control_p = 1;
				continue;
			}
		}
		
		if(buf == "--j" || buf == "-J"){

			if(i < argc - 1){

				JASPAR_FILE = argv[++i];

				bool jaspar_check = is_file_exist(JASPAR_FILE);
				if(jaspar_check == 0){
					cout << "JASPAR matrix does not exist, please insert a JASPAR matrix as input. \n";
					cout << "FATAL ERROR \n";
					exit(EXIT_SUCCESS);
				}
				continue;
			}
		}
		
		if(buf == "--twobit" || buf == "-tb"){

			if(i < argc - 1){

				TWOBIT_FILE = argv[++i];
				control_twobit = 1;

				bool two_bit_check = is_file_exist(TWOBIT_FILE);
				if(two_bit_check == 0){
					cout << "File 2bit does not exist, please insert a 2bit file as input. \n";
					cout << "FATAL ERROR \n";
					exit(EXIT_SUCCESS);
				}
				continue;
			}
		}


	}

	if(control_bed == 0 && control_twobit == 0){

		cout << "BED file and TwoBit file missed or wrong parameters calling!\n";
		cout << "Please insert as input a BED file using -B or --BED annotation before it.\n";
		cout << "Please insert as input a Twobit file using -tb or --twobit annotation before it.\n";

		exit(EXIT_SUCCESS);
	}

	if(control_bed == 0){

		cout << "Wrong BED file immission or BED file missed!\n ";
		cout << "Please insert as input a BED file using -B or --BED annotation before it.\n";

		exit(EXIT_SUCCESS);
	}

	if(control_twobit == 0){

		cout << "Wrong Twobit file immission or Twobit file missed!\n";
		cout << "Please insert as input a Twobit file using -tb or --twobit annotation before it.\n";

		exit(EXIT_SUCCESS);

	}

	if(control_p == 0){

		cout << "WARNING: Sequence length parameter not inserted --> Default parameter is 150.\n";
		cout << "To put a parameter as input write -p or --param before parameter value.";

	}
}

bool is_file_exist(const char *fileName)		//Input files existence control
{
	std::ifstream infile(fileName);
	return infile.good();
}

void display_help() 						//Display help function
{
	cerr << "\n --help: show this message" << endl;
	cerr << "\n --BED -B <file_bed>: input bed file" << endl;
	cerr << endl;

	exit(EXIT_SUCCESS);
}
