#include "Multifa_random_tool.h"

int main(int argc, char *argv[]){

	command_line_parser(argc, argv);					//Parser function called to handle aguments
	if(argc == 1){             //If arguments number is 1 means that no input file has been inserted - display help
		display_help();
	}
	
	vector<matrix_class> MATRIX_VECTOR;
	position_vector_creation(position);

	cout << "Il numero di JASPAR inserite è: " << JASPAR_FILE_vector.size() << endl;
	for(unsigned int i=0; i<JASPAR_FILE_vector.size(); i++){
	
		cout << JASPAR_FILE_vector[i];
		cout << endl;
		matrix_class NEW_MATRIX(JASPAR_FILE_vector[i]);
		MATRIX_VECTOR.emplace_back(NEW_MATRIX);
	}

	multifasta_class MULTIFA(length,n_seq);
	implanting_class IMPLANTED(MATRIX_VECTOR, MULTIFA.multifasta_map);

}

void multifasta_class::multifasta_map_creation(){

	string sequence;

	for(unsigned int j=0; j<n_seq; j++){
		
		sequence.clear();
		
		for(unsigned int i=0; i<length; i++){

			int random_int = random_number(0,3);
			char base = from_n_to_base(random_int);
			sequence = sequence + base;
		}
	
		multifasta_map.insert(pair<int,string>(j+1, sequence));
	}
	
}

int random_number(int range_begin, int range_end){

        mt19937_64 generator (clock());	
        uniform_int_distribution<int> dis(range_begin, range_end);	
        int r_number = dis(generator);;
	return 	r_number;
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

void matrix_class::read_JASPAR(string JASPAR_FILE){			//Function to read JASPAR PWM file, extract values and create a matrix class

	cout << "Reading JASPAR MATRIX file and extracting values...\n";

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
			vector<int> baseQ;				//Initializing baseQ vector of double
			istringstream mystream(line);			//Splitting the line in words
			for (int num; mystream >> num;){		//Put every word(number of matrix), ricorsively, in double variable num
				baseQ.emplace_back(num);		//Put every number(num) in baseQ vector
			}
			matrix.emplace_back(baseQ);			//Put baseQ vector (corrisponding to matrix line values) in our matrix

		}

	}
	file.close();						//Closing file
}

void matrix_class::oligo_creation(){
	
	string oligo;
	for(unsigned int j=0; j<n_oligo; j++){
	
		oligo.clear();

	for (unsigned int i = 0; i < matrix[0].size(); i++) {			//From 0 to number of columns of line 0

		int somma = matrix[0][i] + matrix[1][i] + matrix[2][i]+ matrix[3][i];
		int random_score = random_number(1,somma);
		
		if(random_score <= matrix[0][i]){
			oligo = oligo + 'A';
		}
		else if(random_score > matrix[0][i] && random_score <= (matrix[0][i] + matrix[1][i])){
			oligo = oligo + 'C';
		}
		else if(random_score > (matrix[0][i]+matrix[1][i]) && random_score <= (matrix[0][i] + matrix[1][i] + matrix[2][i])){
			oligo = oligo + 'G';
		}
		else{
			oligo = oligo + 'T';
		}
	}
	oligo_vector.emplace_back(oligo);
	}
}

void matrix_class::check_oligo_number(){

	if(n_oligo > n_seq){
		
		cerr << "The number of oligo can't be > then n_seq.";
		exit(1);
	}
}

void implanting_class::implanting_oligo(vector<matrix_class> MATRIX_VECTOR){
	
	for(int j=0; j<MATRIX_VECTOR.size(); j++){
		int i=0;
		for(map<int,string>::iterator it = multifasta_map_implanted.begin(); it->first <= MATRIX_VECTOR[j].oligo_vector.size() ; it++, i++){
			it->second.replace(position_vector[j], MATRIX_VECTOR[j].matrix_size, MATRIX_VECTOR[j].oligo_vector[i]);
		}
	}
}

void position_vector_creation(string position){

	int index;
	
	while(index != -1){
		index = position.find(",");
		position_vector.emplace_back(stoi(position.substr(0,index)));
		position.erase(0,index+1);
	}
}
/////////////////////////////////////// DEBUG ////////////////////////////////////


void matrix_class::print_debug_matrix(){		//Print matrix function

	cout << "\n" << matrix_name << " " << tf_name << ":" << endl;

	for(unsigned int i=0; i < matrix.size(); i++){
		for(unsigned int j=0; j<matrix[i].size(); j++){

			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
}

void matrix_class::print_oligo_vector(){
	
	cout << endl;
	for(unsigned int i=0; i<oligo_vector.size();i++){

		cout << oligo_vector[i] << endl;
	}
	cout << "\n\n----------------------------------------------------------------------------"<< "\n\n";
}

void multifasta_class::multifasta_outfile(map<int,string> multifasta_map, string filename){

	ofstream outfile;
	outfile.open(filename);
	
	for(map<int,string>::iterator it = multifasta_map.begin(); it != multifasta_map.end(); it++){

		
		outfile << ">random multifasta sequence number " + to_string(it->first) + " containing "+ to_string(length) +" bases:"<<endl;
		outfile << it->second << endl;
		outfile << endl;
	}
	outfile.close();
}

void implanting_class::multifasta_outfile_2(map<int,string> multifasta_map, string filename){

	ofstream outfile;
	outfile.open(filename);
	
	for(map<int,string>::iterator it = multifasta_map.begin(); it != multifasta_map.end(); it++){

		
		outfile << ">random multifasta sequence number " + to_string(it->first) + " containing "+ to_string(length) +" bases:"<<endl;
		outfile << it->second << endl;
		outfile << endl;
	}
	outfile.close();
}
/////////////////////////////////////// PARSER ///////////////////////////////////

void command_line_parser(int argc, char** argv){
	
	const char* const short_opts = "hl:n:j:o:p:";

	//Specifying the expected options
	const option long_opts[] ={
		{"help",      no_argument, nullptr,  'h' },
		{"length",      required_argument, nullptr,  'l' },
		{"jaspar",   required_argument, nullptr,  'j' },
		{"nseq",   required_argument, nullptr,  'n' },
		{"noligo",   required_argument, nullptr,  'o' },
		{"position",   required_argument, nullptr,  'p' },
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
			case 'o' : n_oligo = stoi(optarg); 
				   break;
			case 'p' : position = string(optarg); 
				   break;
			case 'j' : JASPAR_FILE = (string(optarg));
				   is_file_exist(JASPAR_FILE, ("--jaspar || -j number 1"));
				   JASPAR_FILE_vector.emplace_back(JASPAR_FILE);

				   if(argv[optind-1] !=  argv[argc-1]){
					   
					   if(string(argv[optind])[0] != '-'){
						   
						   JASPAR_FILE = (string(argv[optind]));
						   is_file_exist(JASPAR_FILE, ("--jaspar || -j number 2"));
						   JASPAR_FILE_vector.emplace_back(JASPAR_FILE);
						   
						   if(argv[optind]!=  argv[argc-1]){
							   
							   if(string(argv[optind+1])[0] != '-'){
								   
								   JASPAR_FILE = (string(argv[optind+1]));
								   is_file_exist(JASPAR_FILE, ("--jaspar || -j number 3"));
								   JASPAR_FILE_vector.emplace_back(JASPAR_FILE);
							   }
						   }
					   }
				   }
				   break;
			case '?': // Unrecognized option
			default:
				   display_help();
				   break;
		}
	}
}

bool is_file_exist(string fileName, string buf){		//Input files existence control

	struct stat check;
	int regular_check, existing_check;
	const char * C_fileName = fileName.c_str();
	existing_check = stat(C_fileName, &check );

	regular_check = S_ISREG( check.st_mode );

	if ( regular_check == 0 || existing_check != 0) {
		cerr <<"ERROR: "<< buf << " file does not exist!\n"<< endl;
		display_help();
		exit(1);	
	}
	return 0;
}

void display_help() 						//Display help function
{
	cerr << "\n --help || -h show this message" << endl;
	cerr << "\n --length || -l <number> to insert the length of Multifasta sequences (DEFAULT: 500)" << endl;
	cerr << "\n --nseq || -n <number> to insert the number of Multifasta sequences (DAFAULT: 200)" << endl;
	cerr << "\n --noligo || -o <number> to insert the number of oligo to put in Multifasta sequences (DAFAULT: 80)" << endl;
	cerr << "\n --position || -p <number> to insert the position in Multifasta sequences where you want to implant the oligos (DAFAULT: 80)" << endl;
	cerr << endl;

	exit(EXIT_SUCCESS);
}
