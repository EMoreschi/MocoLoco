#include "MocoLoco.h"

int main(int argc, char *argv[]){


	if(argc == 1){             //If arguments number is 1 means that no input file has been inserted - display help
		display_help();
	}

	command_line_parser(argc, argv);					//Parser function called to handle aguments
	string test;
	int p = 0;
	int l = 16;
	int end = 4;
	vector<double> oligo;
        genomic_position obj;
	vector<genomic_position> GEP;					//Initializing GEP --> vector of genomic_position classes
	GEP_creation(BED_FILE, TWOBIT_FILE, GEP); 			//Function to read BED and 2Bit files and create GEP objects vector
	matrix_class JASPAR_MATRIX(JASPAR_FILE);				//Function to read JASPAR PWM file, extract value from it and create a matrix class called JASPAR_MTX
	
	vector<vector<double>> matrix;
	matrix = JASPAR_MATRIX.return_matrix(1);
	JASPAR_MATRIX.print_debug_matrix(matrix, " ");

	matrix = JASPAR_MATRIX.return_norm_matrix(1);
	JASPAR_MATRIX.print_debug_matrix(matrix, " NORMALIZED");

	matrix = JASPAR_MATRIX.return_inverse_matrix(1);
	JASPAR_MATRIX.print_debug_matrix(matrix, " INVERSE COMPLEMENT");
	//Print the matrix for debugging
	matrix = JASPAR_MATRIX.return_log_matrix(1);
	JASPAR_MATRIX.print_debug_matrix(matrix, " LOGARITHMIC");

//	for(int i=0; i<GEP.size();i++){
//
//         test = GEP[i].return_sequence(GEP[i]);					//Print GEP vector for debugging
//	 JASPAR_MATRIX.shifting(test, p, l, oligo);
//    
//	}	
//	JASPAR_MATRIX.shifting(test, p, l, oligo); 
//        test = GEP[1].return_sequence(GEP[1]);					//Print GEP vector for debugging
//	 cout << test << "\n";
//	 JASPAR_MATRIX.shifting(test, p, l, oligo);


//        for(int i=0; i<oligo.size(); i++){
//	cout << oligo[i] << "\n";
// 	}
//	for(int i=0; i<GEP.size();i++){

//	GEP[i].print_debug_GEP(GEP[i]);					//Print GEP vector for debugging
//	}
}

void genomic_position::read_line(string line){				//Read line function: it takes in input each line from BED file 

	istringstream mystream(line);					//Split the line word by word and extract chromosome coordinates (chr, start, end)
	mystream >> chr_coord >> start_coord >> end_coord;		

}

void genomic_position::centering_function ( int start,  int end, int p, const int overhead){	//Centering function: in takes start and end coordinate and re-sets them -
	//following an input parameter value (overhead added)
	int center = (start + end)/2;						
	start_coord = center - p;			//No overhead for start
	end_coord = center + p +overhead;		//Overhead for end
}


void genomic_position::flag_control( int start,  int end){ 	//Flag control function: start coordinates must be < then end coordinates

	if(start > end){		//if start coordinates are > or == then end coordinates, flag is setted to 0
		flag = 0;
	}
	else{ flag = 1;}
}

void GEP_creation(string Bed_file, string Twobit_file, vector<genomic_position> &GEP){		//Function to read BED and 2Bit files and create GEP object vector

	ifstream in(Bed_file); 						//Opening file in lecture mode
	TwoBit * tb;				//Creating a TwoBit* variable called tb
	tb = twobit_open(Twobit_file.c_str());					//Opening 2Bit file with twobit_open function and saved in tb 
	string line; 							//defining line string

	int n_line = 1;							//line counter initialization

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

void matrix_class::shifting(string seq, int p, int l, vector<double> &oligo){
		
	double sum_oligo = 0;
	
	if(p <= seq.size() - matrix_log[0].size() ) {

	for(int i=0; i< matrix_log[0].size(); i++){

			switch(seq[i+p]){

				case 'A':
				       
					sum_oligo = sum_oligo + matrix_log[0][i];
					break;

				case 'C':
				       
					sum_oligo = sum_oligo + matrix_log[1][i];
					break;

				case 'G':
				       
					sum_oligo = sum_oligo + matrix_log[2][i];
					break;

				case 'T':
				       
					sum_oligo = sum_oligo + matrix_log[3][i];
					break;

		}
	}
	
	oligo.emplace_back(sum_oligo);
	shifting(seq, p+1, l, oligo);
	}

}

string genomic_position::return_sequence(genomic_position){ 
       return sequence;
}

void matrix_class::read_JASPAR(string JASPAR_FILE){			//Function to read JASPAR PWM file, extract values and create a matrix class

	ifstream file(JASPAR_FILE);					//opening JASPAR PWM file
	string line;							
	while(getline(file,line)){					//For each line of the file do:

		if(line[0]=='>'){					//If line start with ">"
			istringstream mystream(line);			
			mystream >> matrix_name >> tf;			//Extract the first two words and put into matrix_name string variable and tf string variable
		}

		else{							//Else, if line does not start with ">"
			line.erase(0,line.find('[') +1);		//Take line charachters after "["...
			line.erase(line.find(']'));			//...and line charachters before "]"
			vector<double> baseQ;				//Initializing baseQ vector of double
			istringstream mystream(line);			//Splitting the line in words
			for (double num; mystream >> num;){		//Put every word(number of matrix), ricorsively, in double variable num
				baseQ.emplace_back(num);		//Put every number(num) in baseQ vector
			}
			matrix.emplace_back(baseQ);			//Put baseQ vector (corrisponding to matrix line values) in our matrix

		}

	}
	file.close();						//Closing file
}

vector<double> matrix_class::find_col_sum(vector<vector<double>> matrix){

	vector<double> col_sum;						//Vector of columns sum
	double sum = 0;							//Sum initialized as 0

	for (int i = 0; i < matrix[0].size(); i++) {			//From 0 to number of columns of line 0
		for (int j = 0; j < 4; j++){				//From 0 to 4 (line number)

			sum = sum + matrix[j][i];			//Calculate the sum of columns
		}

		col_sum.emplace_back(sum);				//Put the column sum in vector col_sum
		sum = 0;						//Restore the sum to 0 for the next column
	}
	return col_sum;
}

void matrix_class::matrix_normalization_pseudoc(vector<vector<double>> matrix, double p){  
	
	double norm;							//Norm variable initialized
	vector<double> col_sum = find_col_sum(matrix);

	for (int i = 0; i < matrix.size(); i++) {		//From 0 to number of matrix lines

		vector<double> baseQ;				//baseQ vector to store the lines initialized
		for (int j = 0; j < matrix[i].size(); j++){	//From 0 to number of matrix columns

				norm = matrix[i][j]/col_sum[j];		//Put matrix value (divided for the corresponding column sum) into double variable norm
				baseQ.emplace_back(norm + p);		//Put norm value (with p added) in baseQ vector
		}

			norm_matrix.emplace_back(baseQ);	//Put baseQ vector (which carries line values) in norm_matrix
	}
}

void matrix_class::matrix_normalization(vector<vector<double>> matrix){

	vector<double> col_sum = find_col_sum(matrix);

	for (int i = 0; i < matrix.size(); i++) {		//From 0 to number of matrix lines

		vector<double> baseQ;				//baseQ vector to store the lines initialized
		for (int j = 0; j < matrix[i].size(); j++){	//From 0 to number of matrix columns

				norm_matrix[i][j] = matrix[i][j]/col_sum[j];	//Substitution of first normalized values with new normalized ones
		}
	}
}

void matrix_class::matrix_logarithmic(vector<vector<double>> matrix){
	
	for(int i=0; i < matrix.size(); i++){
		vector<double> baseQ;
		double value_log;

		for(int j=0; j < norm_matrix[i].size(); j++){
			
			value_log = log(norm_matrix[i][j]);
			baseQ.emplace_back(value_log);
		}
		matrix_log.emplace_back(baseQ);
	}
}
	

void matrix_class::inverse_matrix(vector<vector<double>> matrix){

	inverse_complement_matrix = norm_matrix;
	reverse(inverse_complement_matrix.begin(), inverse_complement_matrix.end());
	for (int i = 0; i < 4; i++) {		//From 0 to number of matrix lines
		vector<double> baseQ;
		reverse(inverse_complement_matrix[i].begin(), inverse_complement_matrix[i].end());
	}

}

void matrix_class::find_minmax(vector<vector<double>> matrix){

	for(int i=0; i < matrix[0].size(); i++){
		vector<double> colum;		   	
		for(int j=0; j < matrix.size(); j++){
			colum.emplace_back(matrix[j][i]);
		}
		local_mins.emplace_back(*min_element(colum.begin(),colum.end()));
		local_maxes.emplace_back(*max_element(colum.begin(),colum.end()));
	}
	global_min = accumulate(local_mins.begin(), local_mins.end(), 0.0);
	global_max = accumulate(local_maxes.begin(), local_maxes.end(), 0.0);

}	

vector<vector<double>> matrix_class::return_matrix(int i){

	return matrix;
}
vector<vector<double>> matrix_class::return_norm_matrix(int i){

	return norm_matrix;
}
vector<vector<double>> matrix_class::return_inverse_matrix(int i){

	return inverse_complement_matrix;
}
vector<vector<double>> matrix_class::return_log_matrix(int i){

	return matrix_log;
}

void matrix_class::print_debug_matrix(vector<vector<double>> matrix, string type){			//Debugging of matrix
	
	cout << "\n" << matrix_name << " " << tf << type << ":" << endl;

	for(int i=0; i < matrix.size(); i++){
		for(int j=0; j<matrix[i].size(); j++){

			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
	
	if(type == " LOGARITHMIC"){
	
		for(int i=0; i < local_maxes.size(); i++){
			
			cout  << local_maxes[i] << " ";
			cout  << local_mins[i] << " ";
		}
		
		cout << endl;
		cout << "\nThe global min is: " << global_min << endl;
		cout << "The global max is: " << global_max << endl;
	}
}


void genomic_position::print_debug_GEP(genomic_position){			//Debug function: Print the GEP vector to control the working flow

	cout << ">" << chr_coord << ":" << start_coord << " - " << end_coord << endl;	//Printing chr, start and end coordinates
	cout << sequence << endl;					//Printing sequence

}

void genomic_position::extract_seq(TwoBit* tb, int n_line){			//Extract sequence function: Extract, from Twobit hg38 genome, the DNA sequence with (chr, start, end) coordinates -
	//extracted from Bed line
	if(flag == 1){								//CONTROL: if flag is 1 means that the current line has starting coordinate > end coordinate, so it is correct
		string chrom = chr_coord; 				//Put in chrom the string of chr_coord
		sequence = twobit_sequence(tb,chrom.c_str(),start_coord,end_coord-1); 	//Extract the sequence from the object with the twobit_sequence function
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
				bool dir = isDir(BED_FILE);
				if(dir == 1){
					cout << "ERROR: BED file inserted is a directory!\nPlease insert a BED file.\n!";
					exit(EXIT_SUCCESS);
				}
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
				bool dir = isDir(JASPAR_FILE);
				if(dir == 1){
					cout << "ERROR: JASPAR file inserted is a directory!\nPlease insert a JASPAR file.\n!";
					exit(EXIT_SUCCESS);
				}
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
				bool dir = isDir(TWOBIT_FILE);
				if(dir == 1){
					cout << "ERROR: TWOBIT file inserted is a directory!\nPlease insert a TWOBIT file.\n!";
					exit(EXIT_SUCCESS);
				}
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

bool is_file_exist(string fileName)		//Input files existence control
{
	ifstream infile(fileName);
	if(!infile)
		return 0;
	else{
		return 1;
	}
}

bool isDir(string filename){

	DIR *pDir;
	bool exists = false;
	pDir = opendir(filename.c_str());
	if(pDir != 0){
		exists = true;
		(void)closedir(pDir);
	}
	return exists;
}


void display_help() 						//Display help function
{
	cerr << "\n --help: show this message" << endl;
	cerr << "\n --BED -B <file_bed>: input bed file" << endl;
	cerr << "\n --twobit -tb <file_twobit>: input twobit file" << endl;
	cerr << "\n --j -J <JASPAR_file>: input JASPAR file" << endl;
	cerr << "\n --param -p <parameter>: input parameter to select bases number to keep around the chip seq signal" << endl;
	cerr << endl;

	exit(EXIT_SUCCESS);
}
