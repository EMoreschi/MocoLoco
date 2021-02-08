#include "MocoLoco.h"

const char * BED_FILE; 		//initializing const char variarible for Bed_file input reading //FZ BASTA STRINGHE C, VI BUCO IL PALLONE SE NE VEDO ANCOORA //
int parameter = 150; 		//default parameter 150 //FZ MA NESSUN PARAMETRO VA CHIAMATO PARAMETRO :D PARAMETRO E' UN NOME GENERIICO. DARE UN NOME SENSATO. DOVREBBE ESSERE UNA VARIABILE DI CLASSE DI TIPO STATIC E NON UNA VARIABILE GLOBALE //
const char * TWOBIT_FILE;	//initializing const char variable for Twobit_file input reading
const char * JASPAR_FILE;

int main(int argc, char *argv[]){


	if(argc == 1){             //If arguments number is 1 means that no input file has been inserted - display help
		display_help();
	}

	command_line_parser(argc, argv);					//parser function called to handle aguments

	vector<genomic_position> GEP;
	GEP_creation(BED_FILE, TWOBIT_FILE, GEP); 			//function to read BED and 2Bit files and create GEP objects vector
	jaspar_PWM JASPAR_MTX(JASPAR_FILE);
	JASPAR_MTX.stamp_debug_matrix(JASPAR_MTX);
	stamp_debug(GEP);


}

void genomic_position::read_line(string line){				//Read line function: it takes in input each line from BED file 

	istringstream mystream(line);					//Split the line word by word and extract chromosome coordinates (chr, start, end)
	mystream >> chr_coord >> start_coord >> end_coord;		

}

void genomic_position::centering_function ( int start,  int end, int p){	//Centering function: in takes start and end coordinate and re-set/s them - 
	//following an input parameter value (overhead added)
	int overhead = 25;    //FZ L'OVERHEAD DEVE ESSERE UNA COSTANTE STATIC  DI CLASSE 
	int centro = (start + end)/2; //FZ EVITARE NOMI VARIABILI IN ITALIANO						
	start_coord = centro - p;			//no overhead for start
	end_coord = centro + p +overhead;		//overhead for end
}


void genomic_position::flag_control( int start,  int end){ 	//Flag control function: start coordinates must be < then end coordinates

	if(start > end || start == end){		//if start coordinates are > or == then end coordinates, flag is setted to 0 //FZ START == END E' LEGALE
		flag = 0;
	}
	else{ flag = 1;}
}

void GEP_creation(const char* Bed_file, const char* Twobit_file, vector<genomic_position> &GEP){		//Function to read BED and 2Bit files and create GEP object vector

	ifstream in(Bed_file); 						//Opening file in lecture mode
	TwoBit * tb;							//Creating a TwoBit* variable called tb
	tb = twobit_open(Twobit_file);					//Opening 2Bit file with twobit_open function and saved in tb 
	string line; 							//defining line string
	int n_line = 0;							//line counter initialization //FZ LA CONTA  DELLE LINEE DOVREBBE COMINCIARE DA 1 //

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

void jaspar_PWM::read_JASPAR(const char* file_jaspar){

	ifstream file(file_jaspar);
	string line;
	while(getline(file,line)){;

		if(line[0]=='>'){
			istringstream mystream(line);
			mystream >> matrix_name >> tf;
		}

		else{
			line.erase(0,line.find('[') +1);
			line.erase(line.find(']'));
			vector<double> baseQ;
			istringstream mystream(line);
			for (double num; mystream >> num;){
				baseQ.emplace_back(num);	
			}
			matrix.emplace_back(baseQ);

		}

	}
	file.close();
}

void jaspar_PWM::stamp_debug_matrix(jaspar_PWM){      

	cout << "\n" << matrix_name << "\n" << tf <<  "\n";
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix[i].size(); j++)
			cout << matrix[i][j] << " ";
		cout << endl;
	}
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
	ifstream infile(fileName);
	return infile.good();
}

void display_help() 						//Display help function //FZ AGGIORNARE E TENERE AGGIORNATO
{
	cerr << "\n --help: show this message" << endl;
	cerr << "\n --BED -B <file_bed>: input bed file" << endl;
	cerr << endl;

	exit(EXIT_SUCCESS);
}
