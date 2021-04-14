#include "Multifa_random_tool.h"

unsigned int length = 500;
unsigned int n_seq = 200;
string JASPAR_F;
string oligo_perc = "0";
string position;
string wobble = "0";
vector<string> JASPAR_FILE_vector;
vector<string> matrix_n;
vector<string> matrix_tf;
vector<unsigned int>n_oligo_vector;
vector<unsigned int> position_vector;
vector<unsigned int> wobble_vector;
unsigned int cycles = 1;

int main(int argc, char *argv[]){

	command_line_parser(argc, argv);					//Parser function called to handle aguments

	if(argc != 1){  

		read_input();	
	}

	for(unsigned int i=0; i<cycles; i++){

		implanting_cycle(i);
	}
}

void read_input(){

	generic_vector_creation(oligo_perc, n_oligo_vector);	
	find_oligo_number();
	generic_vector_creation(position, position_vector);	//position vector creation from input positions, passed as a string
	check_oligo_number();
	check_position_vector();
	wobble_vector_creation(wobble);
	check_wobble();
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

void implanting_cycle(unsigned int i){

	if(position.size() != 0){


		map<vector<unsigned int>, vector<vector<unsigned int>>> jaspar_map;
		string matrix_name;
		string tf_name;

		cout << "\nThe number of JASPAR matrices is: " << JASPAR_FILE_vector.size() << endl; //Debug for a meaningful output -> here you can check if you have put the input that you wanted
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
		filling_jaspar_map(jaspar_map,matrix_name,tf_name);
		ordering_matrix_names();
		multifasta_class MULTIFA(length,n_seq,i); 	//generating a random multifasta_class
		
		if(position_vector.size() > 1){		//if the jaspar input are more then 1
			
			check_overlapping(jaspar_map);	//checking if -p implanting position in input are different and if the implanting does not overlap
		}

		check_exceeding(jaspar_map);		//checking if -p implanting position don't bring the oligos to exceed from the sequences length 
		implanting_class IMPLANTED(jaspar_map, MULTIFA.multifasta_map,i);	//implanting the oligos in the position -p gave as input on the multifasta sequences from multifasta_class previouly generated
		matrix_n.clear();
		matrix_tf.clear();	
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

void filling_jaspar_map(map<vector<unsigned int>, vector<vector<unsigned int>>>& jaspar_map, string matrix_name, string tf_name){
	

	for(unsigned int i=0; i<JASPAR_FILE_vector.size(); i++){

		vector<vector<unsigned int>> matrix = read_JASPAR(JASPAR_FILE_vector[i],matrix_name,tf_name);
		matrix_n.emplace_back(matrix_name);		//saving matrix name and tf name in a vector	
		matrix_tf.emplace_back(tf_name);		//for future printing

		vector<unsigned int> parameters;
		parameters.emplace_back(position_vector[i]);
		parameters.emplace_back(wobble_vector[i]);
		parameters.emplace_back(n_oligo_vector[i]);
		
		jaspar_map.insert({parameters,matrix});
		parameters.clear();

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

vector<vector<unsigned int>> read_JASPAR(string JASPAR_FILE, string& matrix_name, string& tf_name){			//Function to read JASPAR PWM file, extract values and create a matrix class

	ifstream file(JASPAR_FILE);					//opening JASPAR PWM file
	string line;
	vector<vector<unsigned int>> matrix;
	
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
	file.close();
	return matrix;
}

void ordering_matrix_names(){

		map<unsigned int,pair<string,string>> names;

		for(unsigned int i=0; i<position_vector.size(); i++){

			pair<string,string> name_tf_pair;
			name_tf_pair.first = matrix_n[i];
			name_tf_pair.second = matrix_tf[i];
			names.insert({position_vector[i],name_tf_pair});
		}
		matrix_n.clear();
		matrix_tf.clear();

		for(map<unsigned int,pair<string,string>>::iterator it = names.begin(); it!=names.end(); it++){

			matrix_n.emplace_back(it->second.first);
			matrix_tf.emplace_back(it->second.second);
		}
}

void implanting_class::oligo_creation(map<vector<unsigned int>,vector<vector<unsigned int>>>::iterator it){

		unsigned int strand;
		string oligo;

		for(unsigned int j=0; j<it->first[2]; j++){

			strand = random_number(0,1);
			oligo.clear();

			for (unsigned int i = 0; i < it->second[0].size(); i++) {			//From 0 to number of columns of line 0

				unsigned int somma = it->second[0][i] + it->second[1][i] + it->second[2][i]+ it->second[3][i];
			unsigned int random_score = random_number(1,somma);

			if(random_score <= it->second[0][i]){
				oligo = oligo + 'A';
			}
			else if(random_score > it->second[0][i] && random_score <= (it->second[0][i] + it->second[1][i])){
				oligo = oligo + 'C';
			}
			else if(random_score > (it->second[0][i]+it->second[1][i]) && random_score <= (it->second[0][i] + it->second[1][i] + it->second[2][i])){
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

void check_position_vector(){

	unsigned int counter = 0;

	for(unsigned int i=0; i<position_vector.size(); i++){
		for(unsigned int j=0; j<position_vector.size(); j++){

			if(position_vector[i] == position_vector[j]){
				
				counter++;
			}
		}
		if(counter != 1){
			cerr << "\nERROR: You have inserted 2 or more equal implanting positions!\nPlease check your -p input." << endl;
			exit(1);
		}
		counter = 0;
	}
}

void check_oligo_number(){

	for(unsigned int i=0; i<n_oligo_vector.size(); i++){
		if(n_oligo_vector[i] > n_seq){

			cerr << "\nERROR: The number of oligo can't be > then n_seq.";
			exit(1);
		}
		if(n_oligo_vector[i] == 0 || n_oligo_vector.size() != position_vector.size()){

			cerr << "WARNING: There is one or more 0% oligo generation frequence" << endl;
		}
	}
}

void check_wobble(){

	for(unsigned int i=0; i<wobble_vector.size(); i++){

		if(wobble_vector[i] > 10){

			cerr << "ERROR: The maximum value allowed for wobble parameters is 10.\nPlease check the -w parameters inserted!" << endl;
			exit(1);
		}
	}
}

void check_input(){

	if(JASPAR_FILE_vector.size() != position_vector.size()){

		cerr << "\nERROR: Please insert the rigth number of Jaspar files and implanting position.\nTheir number need to be equal to make a correct implanting!" << endl;
		exit(1);
	}
}

void check_exceeding(map<vector<unsigned int>,vector<vector<unsigned int>>> jaspar_map){

	for(map<vector<unsigned int>,vector<vector<unsigned int>>>::iterator it = jaspar_map.begin(); it != jaspar_map.end(); it++){
		int is_less_zero = it->first[0] - it->first[1];
			
		if((it->first[0] + it->first[1] + it->second.size()) > length){

			cerr << "\nERROR: Position in input lead the oligos to exceed the length of the sequences.\nImplanting position > " << length << " try to be generated.\nPlease check your implanting position!" << endl;
			exit(1);
		}
		if(is_less_zero < 0){

			cerr << "\nERROR: Position in input lead the oligos to exceed the length of the sequences.\nImplanting position < 0 try to be generated.\nPlease check your implanting position!" << endl;
			exit(1);
		}
	}
}

void check_overlapping(map<vector<unsigned int>,vector<vector<unsigned int>>> jaspar_map){

	map<vector<unsigned int>, vector<vector<unsigned int>>>::iterator it_before = jaspar_map.begin();
	map<vector<unsigned int>, vector<vector<unsigned int>>>::iterator it_after = ++jaspar_map.begin();

	for(unsigned int i=0; i<(position_vector.size()-1); i++, it_after++, it_before++){
		
		if(it_after->first[0] <= (it_before->first[0] + it_before->second[0].size())){
			cerr << "\nERROR: The oligos coming from matrix " << matrix_n[i] << " and " << matrix_n[i+1] << " that you are trying to implant overlap ---> wrong positions inserted!"<< endl;
			exit(1);

			}

		else if((it_after->first[0] - it_after->first[1]) <= (it_before->first[0] + it_before->second[0].size() + it_before->first[1])){

			cerr << "\nERROR: The oligos coming from matrix " << matrix_n[i] << " and " << matrix_n[i+1] << " that you are trying to implant overlap ---> wrong wobbles inserted!"<< endl;
			exit(1);
		}
	}
}

void implanting_class::implanting_oligo(map<vector<unsigned int>, vector<vector<unsigned int>>> jaspar_map){

		int j=0;
		for(map<vector<unsigned int>,vector<vector<unsigned int>>>::iterator it = jaspar_map.begin(); it!=jaspar_map.end(); it++,j++){

			unique_random_generator();				//generating a vector of unique random numbers from 1 to n_seq (length n_seq)
			oligo_creation(it);				//creating a vector of random oligos coming from matrix frequences -> the size of oligos vector follows the input n_oligo_vector values for each matrix
			vector<unsigned int> index_vec;
			
			for(unsigned int w=0; w < it->first[2]; w++){	//Generating a vector of new indexes starting from input position and adding wobble

				unsigned int index = random_number((it->first[0] - it->first[1]), (it->first[0] + it->first[1]));
				index_vec.emplace_back(index);		//index_vec is vector of new idexes (wobble variation added)	
			}

			cout << "\nStarting implant position: "<< it->first[0] << endl;
			cout << "Wobble inserted: " << it->first[1] << endl<<endl;
			cout << matrix_n[j] << ": " << matrix_tf[j] << endl;
			print_debug_matrix(it->second);
			cout << endl;

			map<unsigned int, string>::iterator implant_it;

			for(unsigned int i=0; i<it->first[2]; i++){ 	//from 0 to oligo_vector_size() (for every matrix class)
				
				implant_it = multifasta_map_implanted.find(unique_rnd[i]);	//find in multifasta map the unique_rnd[i](random number) sequence
				implant_it->second.replace(index_vec[i], it->second[0].size(), oligo_vector[i]);		//Implant the oligo generated from jaspar in the right position in the right sequence

				cout << "Implanted string " << oligo_vector[i] << " in sequence number " << unique_rnd[i] << " in position number " << index_vec[i] << "." << endl;
			
			}
				
			cout << endl << "---------------------------------------------------------------------------------------" << endl;	
			
			index_vec.clear();
			unique_rnd.clear();
			oligo_vector.clear();
		}
}

void implanting_class::unique_random_generator(){

	mt19937 eng{random_device{}()};
	
	for(unsigned int i=1; i<=n_seq; i++){

		unique_rnd.emplace_back(i);
	}
	shuffle(unique_rnd.begin(), unique_rnd.end(), eng);
}


void generic_vector_creation(string oligo_perc, vector<unsigned int> &n_oligo_vector){

	int index;
	
	while(index != -1){
		index = oligo_perc.find(",");
		n_oligo_vector.emplace_back(stoi(oligo_perc.substr(0,index)));
		oligo_perc.erase(0,index+1);
	}
}


void wobble_vector_creation(string wobble){

	generic_vector_creation(wobble,wobble_vector);

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


void implanting_class::print_debug_matrix(vector<vector<unsigned int>> matrix){		//Print matrix function

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

void display_help(){

	cerr << "\n --help || -h show this message" << endl;
	cerr << "\n --jaspar || -j <JASPAR_FILE_1, JASPAR_FILE_2, ..., JASPAR_FILE_N> to import the jaspar matrices from which the oligos will be generated." << endl;
	cerr << "\n --length || -l <number> to insert the length of Multifasta sequences (DEFAULT: 500)." << endl;
	cerr << "\n --nseq || -n <number> to insert the number of Multifasta sequences (DEFAULT: 200)." << endl;
	cerr << "\n --oligop || -o <n1,n2,...,nN> to insert the percentage of sequences in which the oligos will be implanted (DEFAULT: 0%) ---- NB: The number of sequences extracted from the percentage will be rounded down." << endl;
	cerr << "\n --position || -p <n1,n2,...,nN> to insert the position in Multifasta sequences where you want to implant the oligos." << endl;
	cerr << "\n --wobble || -w <n1,n2,...,nN> to set the wobble parameter for every implanting position. The implanting position, for every oligo, will be randomly choosen between p-w and p+w interval. (DEFAULT: 0)" << endl;
	cerr << "\n --cycles || -c <number> to choose how many random multifasta files (and implanted also) this tool will produce. (DEFAULT: 1)" << endl;
	cerr << endl;
	exit(EXIT_SUCCESS);
}

bool is_file_exist(string fileName, string buf){		//Input files existence control

	struct stat check;
	int regular_check, existing_check;
	const char * C_fileName = fileName.c_str();
	existing_check = stat(C_fileName, &check );

	regular_check = S_ISREG( check.st_mode );

	if(regular_check == 0 || existing_check != 0){
		
		cerr<<"ERROR: "<<buf<<" file does not exist!\n"<<endl;
		display_help();
		exit(1);
	}
	return 0;
}

