#include "Multifa_random_tool.h"

int main(int argc, char *argv[]){

	command_line_parser(argc, argv);					//Parser function called to handle aguments
	if(argc == 1){             //If arguments number is 1 means that no input file has been inserted - display help
		display_help();
	}

	n_oligo_vector_creation(oligo_perc);	
	find_oligo_number();
	for(unsigned int i=0; i<cycles; i++){
	check_jaspar_exist(i);
	}
}	

void find_oligo_number(){

	for(unsigned int i=0; i<n_oligo_vector.size(); i++){
		
		if(n_oligo_vector[i] > 100){

			cout << "ERROR: The percentage need to be from 0 to 100!" << endl;
			exit(1);
		}
		n_oligo_vector[i] = (n_oligo_vector[i] * n_seq)/100;
	}

}

void check_jaspar_exist(unsigned int i){

	if(position.size() != 0){
		
		position_vector.clear();
		position_vector_creation(position);	//position vector creation from input positions, passed as a string
		wobble_vector_creation(wobble);
		vector<matrix_class> MATRIX_VECTOR;	//MATRIX_VECTOR initialization

		cout << "\nThe number of JASPAR matrices is: " << JASPAR_FILE_vector.size() << endl; //Debug for a beautiful and meaningful output -> here you can check if you have put the input that you wanted
		cout << "The number of implanting position is: " << position_vector.size() << endl;
		cout << "The length of multifasta random sequences is: " << length << endl;
		cout << "The number of multifasta random sequences generated is: " << n_seq << endl;
		cout << "The number of oligo randomly generated for each Jaspar matrix is: ";
		
		for(unsigned int i=0; i<n_oligo_vector.size(); i++){
		cout <<n_oligo_vector[i] << " ";
		}
		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------" << endl;

		check_input();		//controlling that the number of Jaspar matrices in input are = to number of -p (implanting position) in input. If the control is positive, the map<position,jaspar_name> start to be filled	
		unsigned int p=0;
		for(map<unsigned int,string>::iterator it = position_jaspar_map.begin(); it != position_jaspar_map.end(); it++, p++){	//for each element in the map

			matrix_class NEW_MATRIX(it->second, p);		//a new matrix_class is created starting from the jaspar_name string
			MATRIX_VECTOR.emplace_back(NEW_MATRIX);		//and a vector of matrix_class is filled
		}
		
		if(position_vector.size() > 1){		//if the jaspar input are more then 1
			check_overlapping(MATRIX_VECTOR);	//checking if -p implanting position in input are different and if the implanting does not overlap
		}

		check_positions(MATRIX_VECTOR);		//checking if -p implanting position don't bring the oligos to exceed from the sequences length 
		multifasta_class MULTIFA(length,n_seq,i); 	//generating a random multifasta_class
		implanting_class IMPLANTED(MATRIX_VECTOR, MULTIFA.multifasta_map,i);	//implanting the oligos in the position -p gave as input on the multifasta sequences from multifasta_class previouly generated
	
		MATRIX_VECTOR.clear();
	}

	else{
	
		cout << "\nWARNING: No Jaspar matrices and implanting position given as input." << endl;	
		cout << "Generating a Random Multifasta file of " << n_seq << " sequences of " << length << " bases length...\n" << endl;
		cout << "\nThe number of JASPAR matrices is: " << JASPAR_FILE_vector.size() << endl; //Debug for a beautiful and meaningful output -> here you can check if you have put the input that you wanted
		cout << "The number of implanting position is: " << position_vector.size() << endl;
		cout << "The length of multifasta random sequences is: " << length << endl;
		cout << "The number of multifasta random sequences generated is: " << n_seq << endl;
		cout << "The number of oligo randomly generated for each Jaspar matrix is: ";
		
		for(unsigned int i=0; i<n_oligo_vector.size(); i++){
		cout <<n_oligo_vector[i] << " ";
		}
		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------" << endl;
		multifasta_class MULTIFA(length,n_seq,i); 	//generating a random multifasta_class
	}

}

void multifasta_class::multifasta_map_creation(){

	string sequence;

	for(unsigned int j=0; j<n_seq; j++){
		
		sequence.clear();
		
		for(unsigned int i=0; i<length; i++){

			unsigned int random_int = random_number(0,3);
			char base = from_n_to_base(random_int);
			sequence = sequence + base;
		}
	
		multifasta_map.insert(pair<unsigned int,string>(j+1, sequence));
	}
	
}

unsigned int random_number(unsigned int range_begin, unsigned int range_end){

        mt19937_64 generator (clock());	
        uniform_int_distribution<unsigned int> dis(range_begin, range_end);	
        unsigned int r_number = dis(generator);;
	return 	r_number;
}

string reverse_complement(string oligo){

	string reverse_oligo;
	for(unsigned int i=0; i<oligo.size(); i++){

		char base;
		base = oligo[i];
		switch (base) {

			case 'A' : reverse_oligo.append("T"); 
				   break;
			case 'T' : reverse_oligo.append("A"); 
				   break;
			case 'G' : reverse_oligo.append("C"); 
				   break;
			case 'C' : reverse_oligo.append("G"); 
				   break;
			case 'N' : reverse_oligo.append("N"); 
				   break;
		}
	}

	reverse(reverse_oligo.begin(), reverse_oligo.end());
	return reverse_oligo;
}


char from_n_to_base (unsigned int n){

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
					
					cerr << "ERROR: Generation of random sequence failed. Please check the random generation code.!" << endl;
					exit(1);
					break;
			}
			return base;
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
			vector<unsigned int> baseQ;				//Initializing baseQ vector of double
			istringstream mystream(line);			//Splitting the line in words
			for (unsigned int num; mystream >> num;){		//Put every word(number of matrix), ricorsively, in double variable num
				baseQ.emplace_back(num);		//Put every number(num) in baseQ vector
			}
			matrix.emplace_back(baseQ);			//Put baseQ vector (corrisponding to matrix line values) in our matrix

		}

	}
	file.close();						//Closing file
}

void matrix_class::oligo_creation(unsigned int p){
	
	unsigned int strand;
	string oligo;
	
	for(unsigned int j=0; j<n_oligo_vector[p]; j++){
	
		strand = random_number(0,1);
		oligo.clear();

	for (unsigned int i = 0; i < matrix[0].size(); i++) {			//From 0 to number of columns of line 0

		unsigned int somma = matrix[0][i] + matrix[1][i] + matrix[2][i]+ matrix[3][i];
		unsigned int random_score = random_number(1,somma);
		
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
	
	if(strand == 0){
	oligo_vector.emplace_back(oligo);
	}
	else{
		oligo = reverse_complement(oligo); 
		oligo_vector.emplace_back(oligo);
	}

	}
}

void matrix_class::check_oligo_number(){

	for(unsigned int i=0; i<n_oligo_vector.size(); i++){
		if(n_oligo_vector[i] > n_seq){

			cerr << "\nERROR: The number of oligo can't be > then n_seq.";
			exit(1);
		}
	}
}

void check_input(){

	if(JASPAR_FILE_vector.size() != position_vector.size()){

		cerr << "\nERROR: Please insert the rigth number of Jaspar files and implanting position.\nTheir number need to be equal to make a correct implanting!" << endl;
		exit(1);
	}

		for(unsigned int i=0; i<position_vector.size(); i++){	
			position_jaspar_map.insert(pair<unsigned int,string>(position_vector[i], JASPAR_FILE_vector[i]));
		}
}

void check_positions(vector<matrix_class> MATRIX_VECTOR){

	for(unsigned int i=0; i<position_vector.size(); i++){

		if(position_vector[i] + MATRIX_VECTOR[i].matrix_size > length){

			cerr << "\nERROR: Position in input lead the oligos to exceed the length of the sequences.\nPlease check your implanting position!" << endl;
			exit(1);
		}
	}
}

void check_overlapping(vector<matrix_class> MATRIX_VECTOR){

	map<unsigned int, string>::iterator it_before = position_jaspar_map.begin();
	map<unsigned int, string>::iterator it_after = ++position_jaspar_map.begin();
	
	if(position_vector.size() != position_jaspar_map.size()){

		cerr << "\nERROR: You have inserted 2 or more equal implanting positions!\nPlease check your -p input." << endl;
		exit(1);
	}

	for(unsigned int i=0; i<MATRIX_VECTOR.size(); i++, it_after++, it_before++){
		
		if(it_after->first >= it_before->first && it_after->first < (it_before->first + MATRIX_VECTOR[i].matrix_size)){

				cerr << "\nERROR: The oligos coming from matrix " << i+1 << " and " << i+2 << " that you are trying to implant overlap!"<< endl;
				exit(1);
		}
	}
}

void implanting_class::implanting_oligo(vector<matrix_class> MATRIX_VECTOR){

		map<unsigned int,string>::iterator pos_it = position_jaspar_map.begin(); 

		for(unsigned int j=0; j<MATRIX_VECTOR.size(); j++){		//for every matrix_class created
		
			unique_random_generator();				//generating a vector of unique random numbers from 1 to n_seq (length n_seq)

			vector<unsigned int> index_vec;
			
			for(unsigned int w=0; w<MATRIX_VECTOR[j].oligo_vector.size(); w++){	//Generating a vector of new indexes starting from input position and adding wobble

				unsigned int index = random_number((pos_it->first - wobble_vector[j]), (pos_it->first + wobble_vector[j]));
				index_vec.emplace_back(index);		//index_vec is vector of new idexes (wobble variation added)	
			}

			map<unsigned int, string>::iterator it;
		
			cout << "\nStarting implant position: "<< pos_it->first << endl;
			cout << "Wobble inserted: " << wobble_vector[j] << endl;
			cout << "Matrix: " << pos_it->second << endl;
			MATRIX_VECTOR[j].print_debug_matrix();
			cout << endl;

			for(unsigned int i=0; i<MATRIX_VECTOR[j].oligo_vector.size(); i++){ 	//from 0 to oligo_vector_size() (for every matrix class)
				
				it = multifasta_map_implanted.find(unique_rnd[i]);	//find in multifasta map the unique_rnd[i](random number) sequence
				it->second.replace(index_vec[i], MATRIX_VECTOR[j].matrix_size, MATRIX_VECTOR[j].oligo_vector[i]);		//Implant the oligo generated from jaspar in the right position in the right sequence

				cout << "Implanted string " << MATRIX_VECTOR[j].oligo_vector[i] << " in sequence number " << unique_rnd[i] << " in position number " << index_vec[i] << "." << endl;
			
			}
				
			cout << endl << "---------------------------------------------------------------------------------------" << endl;	
			++pos_it;
			
			index_vec.clear();
			unique_rnd.clear();
		}
}

void implanting_class::unique_random_generator(){

	mt19937 eng{random_device{}()};
	
	for(unsigned int i=1; i<=n_seq; i++){

		unique_rnd.emplace_back(i);
	}
	shuffle(unique_rnd.begin(), unique_rnd.end(), eng);
}

void n_oligo_vector_creation(string position){

	int index;
	
	while(index != -1){
		index = oligo_perc.find(",");
		n_oligo_vector.emplace_back(stoi(oligo_perc.substr(0,index)));
		oligo_perc.erase(0,index+1);
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

void wobble_vector_creation(string wobble){

	int index;
	
	while(index != -1){
		index = wobble.find(",");
		wobble_vector.emplace_back(stoi(wobble.substr(0,index)));
		wobble.erase(0,index+1);
	}

	if(position_vector.size() < wobble_vector.size()){

		cout << "WARNING: wobbles parameters inserted as input are more then position parameters!" << endl;
		cout << "Please check your input data. " << endl;
	}	

	if(position_vector.size() > wobble_vector.size()){

		int difference = position_vector.size() - wobble_vector.size();

		for(int i=0; i<difference; i++){

			wobble_vector.emplace_back(0);
		}
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

void multifasta_class::multifasta_outfile(map<unsigned int,string> multifasta_map, string filename){

	ofstream outfile;
	outfile.open(filename);
	
	for(map<unsigned int,string>::iterator it = multifasta_map.begin(); it != multifasta_map.end(); it++){

		
		outfile << ">random multifasta sequence number " + to_string(it->first) + " containing "+ to_string(length) +" bases:"<<endl;
		outfile << it->second << endl;
		outfile << endl;
	}
	outfile.close();
}

void implanting_class::multifasta_outfile_2(map<unsigned int,string> multifasta_map, string filename){

	ofstream outfile;
	outfile.open(filename);
	
	for(map<unsigned int,string>::iterator it = multifasta_map.begin(); it != multifasta_map.end(); it++){

		
		outfile << ">random multifasta sequence number " + to_string(it->first) + " containing "+ to_string(length) +" bases:"<<endl;
		outfile << it->second << endl;
		outfile << endl;
	}
	outfile.close();
}

/////////////////////////////////////// PARSER ///////////////////////////////////

void command_line_parser(int argc, char** argv){
	
	const char* const short_opts = "hl:n:c:j:o:p:w:";

	//Specifying the expected options
	const option long_opts[] ={
		{"help",      no_argument, nullptr,  'h' },
		{"length",      required_argument, nullptr,  'l' },
		{"jaspar",   required_argument, nullptr,  'j' },
		{"nseq",   required_argument, nullptr,  'n' },
		{"cycles",   required_argument, nullptr,  'c' },
		{"oligop",   required_argument, nullptr,  'o' },
		{"position",   required_argument, nullptr,  'p' },
		{"wobble",   required_argument, nullptr,  'w' },
		{nullptr, no_argument, nullptr,  0   }
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
			case 'c' : cycles = stoi(optarg); 
				   break;
			case 'o' : oligo_perc = string(optarg); 
				   break;
			case 'p' : position = string(optarg); 
				   break;
			case 'w' : wobble = string(optarg); 
				   break;
			case 'j' : JASPAR_F = (string(optarg));
				   is_file_exist(JASPAR_F, ("--jaspar || -j number 1"));
				   JASPAR_FILE_vector.emplace_back(JASPAR_F);
				   for (;optind < argc && *argv[optind] != '-';optind++){
					   JASPAR_F = (string(argv[optind]));
					   is_file_exist(JASPAR_F, ("--jaspar || -j one of files do not exist"));
					   JASPAR_FILE_vector.emplace_back(JASPAR_F);
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
	cerr << "\n --jaspar || -j <JASPAR_FILE_1, JASPAR_FILE_2, ..., JASPAR_FILE_N> to import the jaspar matrices from which the oligos will be generated." << endl;
	cerr << "\n --length || -l <number> to insert the length of Multifasta sequences (DEFAULT: 500)" << endl;
	cerr << "\n --nseq || -n <number> to insert the number of Multifasta sequences (DEFAULT: 200)" << endl;
	cerr << "\n --oligop || -o <number> to insert the percentage of sequences in which the oligos will be implanted (DEFAULT: 80%) ---- NB: The number of sequences extracted from the percentage will be rounded down" << endl;
	cerr << "\n --position || -p <n1,n2..,nN> to insert the position in Multifasta sequences where you want to implant the oligos (DEFAULT: 10)" << endl;
	cerr << "\n --cycles || -c <number> to choose how many random multifasta files (and implanted also) this tool will produce" << endl;
	cerr << endl;

	exit(EXIT_SUCCESS);
}
