#include "MocoLoco.h"

int main(int argc, char *argv[]){


	if(argc == 1){             //If arguments number is 1 means that no input file has been inserted - display help
		display_help();
	}

	command_line_parser(argc, argv);					//Parser function called to handle aguments
	vector<bed_class> GEP;					//Initializing GEP --> vector of bed_class classes
	vector<oligo_class> oligos_vector;
	GEP_creation(BED_FILE, TWOBIT_FILE, GEP); 			//Function to read BED and 2Bit files and create GEP objects vector
	matrix_class JASPAR_MATRIX(JASPAR_FILE);				//Function to read JASPAR PWM file, extract value from it and create a matrix class called JASPAR_MTX
	
	vector<vector<double>> matrix;
	vector<vector<double>> matrix_log;
	matrix = JASPAR_MATRIX.return_matrix(1);
	JASPAR_MATRIX.print_debug_matrix(matrix, " ");

	matrix = JASPAR_MATRIX.return_norm_matrix(1);
	JASPAR_MATRIX.print_debug_matrix(matrix, " NORMALIZED");

	matrix = JASPAR_MATRIX.return_inverse_matrix(1);
	JASPAR_MATRIX.print_debug_matrix(matrix, " INVERSE COMPLEMENT");
	//Print the matrix for debugging
	matrix_log = JASPAR_MATRIX.return_log_matrix(1);
	JASPAR_MATRIX.print_debug_matrix(matrix_log, " LOGARITHMIC");

	for(int i=0; i<5; i++){

	string sequence = GEP[i].return_sequence(GEP[i]);
	oligo_class SHIFTING(matrix_log, sequence);
	oligos_vector.emplace_back(SHIFTING);
	vector<double> oligo_scores = SHIFTING.return_oligo_scores(1);
	double best_score = SHIFTING.return_best_score(1);
	cout << endl;
	cout << sequence << endl;

	for(int j=0; j<oligo_scores.size(); j++){
		
		cout << oligo_scores[j] << " ";
	}
	cout << endl;
	cout << "The best oligo in sequence has a score of " << best_score << endl;
	}

//	for(int i=0; i<5; i++){
//
//		cout << SHIFTING.oligo_scores[i] << endl;
//	}
//	}
//	for(int i=0; i<GEP.size();i++){

//	GEP[i].print_debug_GEP(GEP[i]);					//Print GEP vector for debugging
//	}
 return 0;
}

void bed_class::read_line(string line){				//Read line function: it takes in input each line from BED file 

	istringstream mystream(line);					//Split the line word by word and extract chromosome coordinates (chr, start, end)
	mystream >> chr_coord >> start_coord >> end_coord;		

}

void bed_class::centering_function ( int start,  int end, int parameter, const int overhead){	//Centering function: in takes start and end coordinate and re-sets them -
	//following an input parameter value (overhead added)
	int center = (start + end)/2;						
	start_coord = center - parameter;			//No overhead for start
	end_coord = center + parameter +overhead;		//Overhead for end
}


void bed_class::flag_control( int start,  int end){ 	//Flag control function: start coordinates must be < then end coordinates

	if(start > end){		//if start coordinates are > or == then end coordinates, flag is setted to 0
		flag = 0;
	}
	else{ flag = 1;}
}

void GEP_creation(string Bed_file, string Twobit_file, vector<bed_class> &GEP){		//Function to read BED and 2Bit files and create GEP object vector

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

		bed_class new_class(parameter,line,tb, n_line);  //Called the object constructor passing the Bed line, parameter P, twobit file tb, and the line counter n_line
		GEP.emplace_back(new_class);				//Put the new object in the GEP vector with emplace function

		n_line = n_line + 1;					//pass to next line 

	}
}

void oligo_class::shifting(vector<vector<double>> matrix, string sequence, int s_iterator){
		
	double sum_oligo = 0;
	
	if(s_iterator <= sequence.size() - matrix[0].size() ) {

	for(int i=0; i< matrix[0].size(); i++){

			switch(sequence[i+s_iterator]){

				case 'A':
				       
					sum_oligo = sum_oligo + matrix[0][i];
					break;

				case 'C':
				       
					sum_oligo = sum_oligo + matrix[1][i];
					break;

				case 'G':
				       
					sum_oligo = sum_oligo + matrix[2][i];
					break;

				case 'T':
				       
					sum_oligo = sum_oligo + matrix[3][i];
					break;
				
				default:
				       
					sum_oligo = sum_oligo + o_matrix_mins[i];
					break;

		}
	}
	
	oligo_scores.emplace_back(sum_oligo);
	shifting(matrix, sequence, s_iterator+1);

	}

}

string bed_class::return_sequence(bed_class){ 
       return sequence;
}

void matrix_class::read_JASPAR(string JASPAR_FILE){			//Function to read JASPAR PWM file, extract values and create a matrix class

	ifstream file(JASPAR_FILE);					//opening JASPAR PWM file
	string line;							
	while(getline(file,line)){					//For each line of the file do:

		if(line[0]=='>'){					//If line start with ">"
			istringstream mystream(line);			
			mystream >> matrix_name >> tf_name;			//Extract the first two words and put into matrix_name string variable and tf_name string variable
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

void oligo_class::find_minmax(vector<vector<double>> matrix){

	for(int i=0; i < matrix[0].size(); i++){
		vector<double> colum;		   	
		for(int j=0; j < matrix.size(); j++){
			colum.emplace_back(matrix[j][i]);
		}
		o_matrix_mins.emplace_back(*min_element(colum.begin(),colum.end()));
		o_matrix_maxes.emplace_back(*max_element(colum.begin(),colum.end()));
	}
	min_possible_score = accumulate(o_matrix_mins.begin(), o_matrix_mins.end(), 0.0);
	max_possible_score = accumulate(o_matrix_maxes.begin(), o_matrix_maxes.end(), 0.0);

}	

void oligo_class::find_best_score(vector<double> oligo_scores){

	best_score = *max_element(oligo_scores.begin(), oligo_scores.end());

}

void bed_class::extract_seq(TwoBit* tb, int n_line){			//Extract sequence function: Extract, from Twobit hg38 genome, the DNA sequence with (chr, start, end) coordinates -
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

/////DEBUG/////////////////////////////////////////////////////////

vector<double> oligo_class::return_oligo_scores(int i){

	return oligo_scores;
}

double oligo_class::return_best_score(int i){

	return best_score;
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
	
	cout << "\n" << matrix_name << " " << tf_name << type << ":" << endl;

	for(int i=0; i < matrix.size(); i++){
		for(int j=0; j<matrix[i].size(); j++){

			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
	
}


void bed_class::print_debug_GEP(bed_class){			//Debug function: Print the GEP vector to control the working flow

	cout << ">" << chr_coord << ":" << start_coord << " - " << end_coord << endl;	//Printing chr, start and end coordinates
	cout << sequence << endl;					//Printing sequence

}
///////////////////////////////////////////////////////////////////////////////////////////////


////////////////////PARSER////////////////////////////////////////////////////////////////////
void command_line_parser(int argc, char **argv){

	for(int i = 1; i < argc; i++){

		string buf = argv[i];

		if(buf == "--help" || buf == "-h"){

			display_help();

		}

		else if(buf == "--BED" || buf == "-B"){

			if(i < argc - 1){

				BED_FILE = argv[++i];
				is_file_exist(BED_FILE, buf);
				continue;

			}
		}

		else if(buf == "--param" || buf == "-p"){

			if(i < argc - 1){
				try{

				parameter = stoi(argv[++i]);}
				catch(exception &err)
				{
					cerr<< buf <<" parameter is not a number"<<endl;
				}
				continue;
			}
		}

		else if(buf == "--jaspar" || buf == "-J"){

			if(i < argc - 1){

				JASPAR_FILE = argv[++i];
				is_file_exist(JASPAR_FILE, buf);
				continue;

			}
		}

		else if(buf == "--twobit" || buf == "-tb"){

			if(i < argc - 1){

				TWOBIT_FILE = argv[++i];
				is_file_exist(TWOBIT_FILE, buf);
				continue;
			}
		}
		else
		{
			cerr << "Unknown option: " << buf <<  endl;
			exit(1);
		}

	}
}

bool is_file_exist(string fileName, string buf)		//Input files existence control
{
	DIR *Dir;
	int result;
	const char * C_fileName = fileName.c_str();
	result = access (C_fileName, W_OK);
	Dir = opendir(C_fileName);
	if (result == 0 && Dir == 0 ) {
		return result;
	}
	else
	{ 
		cerr <<"ERROR: "<< buf << " parameter has wrong argument"<< endl;
		exit(1);	
	}
}


void display_help() 						//Display help function
{
	cerr << "\n --help show this message" << endl;
	cerr << "\n --BED -B <file_bed>: input bed file" << endl;
	cerr << "\n --twobit -tb <file_twobit>: input twobit file" << endl;
	cerr << "\n --jaspar -J <JASPAR_file>: input JASPAR file" << endl;
	cerr << "\n --param -p <parameter>: input parameter to select bases number to keep around the chip seq signal" << endl;
	cerr << endl;

	exit(EXIT_SUCCESS);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
