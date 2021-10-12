#include "Prova_BED_tool.h"

/////////////////// GLOBAL VARIABLES ////////////////////////////////////////////////////////////////

unsigned int length = 500;
unsigned int n_seq = 200;
string BED_FILE;   
string TWOBIT_FILE;
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
string freq_strand_plus;
vector<unsigned int> freq_strand_plus_vector;
vector<bool> plus_minus;
vector<vector<bool>> plus_minus_matrix;

/////////////////////// MAIN FUNCTION //////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]){
	
	//Parser function to handle input arguments
	command_line_parser(argc, argv);
	
	if(BED_FILE.size() == 0){
		pathway_multifasta();
	} else{
		pathway_bed();
	}
}

void pathway_multifasta(){

	if(position.size() != 0){  

		read_input();	
	
	}
	//For evry cycle -c inserted as input do inplanting_cycle function --> creating a random multifa + implanting (if any implaning position is inserted as input)
	for(unsigned int i=0; i<cycles; i++){

		implanting_cycle(i);
	}
}	

void pathway_bed(){
	cout << BED_FILE << endl; 
}
//////////////////////// INPUT READING AND CHECKING ////////////////////////////////////////////////////////////

//Function to read the input passed on command line and fill the vector to store the input information
void read_input(){

	//Calling of generic vector creation function to transform a string into an unsigned int vector (For -o parameters)
	generic_vector_creation(oligo_perc, n_oligo_vector);	

	//Check if oligo implanting percentage is between 0 and 100
	find_oligo_number();
	
	//Calling of generic vector creation function to transform a string into an unsigned int vector (For -p parameters)
	generic_vector_creation(position, position_vector);

	//Cheking if to any positions (-p) passed as input correspond an oligo implanting percentage (-o) --> Print a Warning
	check_oligo_number();

	//Checking if there are two or more equal implanting position --> If any print error and exit
	check_position_vector();
	
	//Calling of wobble-specific vector creation function to transform a string into an unsigned int vector (For -w parameters)
	wobble_vector_creation(wobble);

	//Checking if wobble parameters passed as input are not greater than 10 --> if any, print an error
	check_wobble();

	//Check if any FWD/REV frequency is inserted --> if any, call the generic_vector_creation function		
	if(freq_strand_plus.size() != 0){
			
		//Calling of generic vector creation function to transform a string into an unsigned int vector (For -f parameters)
		generic_vector_creation(freq_strand_plus, freq_strand_plus_vector);
	}
	
	check_frequence();
}

//Function to convert a string in a vector of unsigned int
void generic_vector_creation(string oligo_perc, vector<unsigned int> &n_oligo_vector){

	int index;
	
	//While index points to something belonging to the string (if not his value is -1)
	while(index != -1){
		
		//Find the index of the first "," in the string
		index = oligo_perc.find(",");

		//Save all the string charachter from 0 to index(,)
		n_oligo_vector.emplace_back(stoi(oligo_perc.substr(0,index)));

		//Delete from the string the value just saved and restart the process with a new cycle
		oligo_perc.erase(0,index+1);
	}
}

//Checking if the implanting percentage overcomes 100% or if it is lower than 0% --> Exit and print error if it happens
void find_oligo_number(){

	for(unsigned int i=0; i<n_oligo_vector.size(); i++){
		
		if(n_oligo_vector[i] > 100 || n_oligo_vector[i] < 0){

			cout << "ERROR: The percentage need to be from 0 to 100!" << endl;
			exit(1);
		}
		
		//Transforming the implanting percentage into the real number of implants that need to be done
		n_oligo_vector[i] = (n_oligo_vector[i] * n_seq)/100;
	}

}

//Cheking if to any positions (-p) passed as input correspond an oligo implanting percentage (-o) --> Print a Warning
void check_oligo_number(){

	for(unsigned int i=0; i<n_oligo_vector.size(); i++){
		
		//If there is a oligo percentage equal to 0 --> print a warning
		if(n_oligo_vector[i] == 0){

			cerr << "WARNING: There is one or more 0% oligo generation frequence" << endl;
		}

		//If number of oligo implanting percentages and positions in input are different --> print a warning
		if(n_oligo_vector.size() != position_vector.size()){

			cerr << "WARNING: The number of -o (oligo implanting percentage) and -p (position) parameters are different. Some parameter can be setted as default. Please check your input command." << endl;
		}
	}
}

//Checking if there are two or more equal implanting position --> If any print error and exit
void check_position_vector(){

	unsigned int counter = 0;

	//Scrolling the position vector with i variable
	for(unsigned int i=0; i<position_vector.size(); i++){
		
		//Scrolling again the same vector with j variable
		for(unsigned int j=0; j<position_vector.size(); j++){
			
			//Increment the counter if there are some equal values
			if(position_vector[i] == position_vector[j]){
				
				counter++;
			}
		}
		
		//For a correct position vector input the counter must be 1 every cycle --> if not, print an error
		if(counter != 1){
			cerr << "\nERROR: You have inserted 2 or more equal implanting positions!\nPlease check your -p input." << endl;
			exit(1);
		}
		
		//Re-setting the counter for the next cycle
		counter = 0;
	}
}

//Function (specific for wobble string) to create a vector of wobble and then re-set the positions
void wobble_vector_creation(string wobble){

	//Calling the generic vector creation function to transform the wobble input string into an unsigned int vector
	generic_vector_creation(wobble,wobble_vector);

	//If the positions inserted are less than wobbles --> print a warning
	if(position_vector.size() < wobble_vector.size()){

		cout << "WARNING: wobbles parameters inserted as input are more then position parameters!" << endl;
		cout << "Please check your input data. " << endl;
	}	
	
	//If positions inserted are more than wobbles set the missing wobble values to 0
	if(position_vector.size() > wobble_vector.size()){

		int difference = position_vector.size() - wobble_vector.size();

		for(int i=0; i<difference; i++){

			wobble_vector.emplace_back(0);
		}
	}
}

//Checking if wobble parameters passed as input are not greater than 10 --> if any, print an error
void check_wobble(){

	for(unsigned int i=0; i<wobble_vector.size(); i++){

		if(wobble_vector[i] > 10){

			cerr << "ERROR: The maximum value allowed for wobble parameters is 10.\nPlease check the -w parameters inserted!" << endl;
			exit(1);
		}
	}
}

//Check and compare the frequencies number to the jaspar matrices number
void check_frequence(){

	unsigned int freq_number = freq_strand_plus_vector.size();

	//If freq number is 0 and there is a Jaspar matrix that is going to be implanted
	if(freq_number == 0 && JASPAR_FILE_vector.size() > 0){

		cerr << "WARNING: No frequence inserted, the tool will set all the frequences on 50%" << endl;
		
		//Set the frequency at 50% in frequency (strand plus) vector
		for(unsigned int i=0; i<JASPAR_FILE_vector.size(); i++){

			freq_strand_plus_vector.emplace_back(50);
		}
	}
	
	//Else if number of frequences insterted are more than 0 but less than the matrices that are going to be implanted
	else if(freq_number > 0 && freq_number < JASPAR_FILE_vector.size()){
		
		cerr << "WARNING: Frequence number inserted is less than matrix inserted, the tool will set on 50% the frequences not specified" << endl;
		
		//Set on 50% the frequences not specified until the size of freq are equal to jaspar matrices number
		while(freq_strand_plus_vector.size() != JASPAR_FILE_vector.size()){

			freq_strand_plus_vector.emplace_back(50);
			}
	}
	
	//Else if there are more frequencies than the matrices in input --> print a warning and delete the excess frequencies
	else if(freq_number > JASPAR_FILE_vector.size()){

		cerr << "WARNING: Frequence number inserted is more than matrix inserted, Please check your input parameters " << endl;
	
		//Erase function to delete the excess frequencies
		freq_strand_plus_vector.erase(freq_strand_plus_vector.begin() + JASPAR_FILE_vector.size(), freq_strand_plus_vector.end());

	}

}

//////////////////////////////// IMPLANTING CYCLES ////////////////////////////////////////////////////////////////

//Function to handle random multifasta creation and any implants
void implanting_cycle(unsigned int i){

	//If at least one implanting position is inserted 
	if(position.size() != 0){


		map<vector<unsigned int>, vector<vector<unsigned int>>> jaspar_map;
		string matrix_name;
		string tf_name;

		//Debug for a meaningful output -> here you can check if the inputs inserted are correct
		cout << "\nThe number of JASPAR matrices is: " << JASPAR_FILE_vector.size() << endl;
		cout << "The number of implanting position is: " << position_vector.size() << endl;
		cout << "The length of multifasta random sequences is: " << length << endl;
		cout << "The number of multifasta random sequences generated is: " << n_seq << endl;
		cout << "The Jaspar matrices inserted as input are: ";
		
		//Printing of jaspar matrices inserted as input
		for(unsigned int i=0; i<JASPAR_FILE_vector.size(); i++){
			
			cout <<JASPAR_FILE_vector[i] << " ";
		}

		cout << endl;

		cout << "The implanting positions for each Jaspar matrix are: ";
		
		//Printing of the implanting positions inserted as input
		for(unsigned int i=0; i<position_vector.size(); i++){
			
			cout <<position_vector[i] << " ";
		}

		cout << endl;

		cout << "The number of oligos randomly generated for each Jaspar matrix are: ";

		//Printing of oligo implanting percentages inserted as input
		for(unsigned int i=0; i<n_oligo_vector.size(); i++){

			cout <<n_oligo_vector[i] << " ";
		}

		cout << endl;

		cout << "The frequences of fwd strand oligo generation for each Jaspar matrix are: ";

		//Printing of FWD strance frequencies inserted as input
		for(unsigned int i=0; i<freq_strand_plus_vector.size(); i++){
			
			cout <<freq_strand_plus_vector[i] << "% ";
		}

		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------" << endl;
		
		
		//Controlling that the number of Jaspar matrices in input are equal to number of -p (implanting position) in input. If the control is positive, the tool associate each matirx to an implanting position --> the map <position,jaspar_name> starts to be filled	
		check_input();	

		//Function to associate each Jaspar matrix to his parameters correctly
		filling_jaspar_map(jaspar_map,matrix_name,tf_name);

		//Function to order matrix name and tf name to mantain their association with the matrix
		ordering_matrix_names();

		//Function to generate a random distribution of FDW and REV matrix implanting
		freq_strand_plus_matrix(jaspar_map);

		//Generation of a random multifasta class
		multifasta_class MULTIFA(length,n_seq,i);

		//if the jaspar input are more then 1
		if(position_vector.size() > 1){
			
			//checking if the implanting positions don't bring to any overlap
			check_overlapping(jaspar_map);
		}

		//checking if -p implanting positions don't bring the oligos to exceed from the sequences length 
		check_exceeding(jaspar_map);
 
		//implanting class constructor calling --> The map matrix + parameters, the random multifa, and the current cycle number are passed to constructor
		implanting_class IMPLANTED(jaspar_map, MULTIFA.multifasta_map,i);
		
		matrix_n.clear();
		matrix_tf.clear();	
	}
	
	//Else if there are no implanting positions passed as input
	else{
		
		//Printing of warning messages
		cout << "\nWARNING: No Jaspar matrices and implanting position given as input." << endl;	
		cout << "Generating a Random Multifasta file of " << n_seq << " sequences of " << length << " bases length...\n" << endl;
		//Debug for a meaningful output -> here you can check if the inputs inserted are correct
		cout << "\nThe number of JASPAR matrices is: " << JASPAR_FILE_vector.size() << endl;
		cout << "The number of implanting position is: " << position_vector.size() << endl;
		cout << "The length of multifasta random sequences is: " << length << endl;
		cout << "The number of multifasta random sequences generated is: " << n_seq << endl;
		cout << "The number of oligo randomly generated for each Jaspar matrix is: ";
		cout << endl << endl;
		cout << "-----------------------------------------------------------------------------------" << endl;
		
		//Generating a simple random multifasta file
		multifasta_class MULTIFA(length,n_seq,i);
	}

}

//Checking if Jaspar matreces number is equal to positions number inserted as input --> if not print an error and exit
void check_input(){

	if(JASPAR_FILE_vector.size() != position_vector.size()){

		cerr << "\nERROR: Please insert the rigth number of Jaspar files and implanting position.\nTheir number need to be equal to make a correct implanting!" << endl;
		exit(1);
	}
}

//**************************************************************************************************************************
//Function to create a map which associate, to each jaspar, a vector containing parameters (position, wobble, implants number, FWD strand frequency)
void filling_jaspar_map(map<vector<unsigned int>, vector<vector<unsigned int>>>& jaspar_map, string matrix_name, string tf_name){
	
	//For every Jaspar matrix passed
	for(unsigned int i=0; i<JASPAR_FILE_vector.size(); i++){
		
		//Extract the matrix values with read_JASPAR function and save it in a vector of vector variable
		vector<vector<unsigned int>> matrix = read_JASPAR(JASPAR_FILE_vector[i],matrix_name,tf_name);

		//saving matrix name and tf name in a vector for future printing
		matrix_n.emplace_back(matrix_name);	
		matrix_tf.emplace_back(tf_name);

		//Filling the parameters vector with (1) position, (2) wobble, (3) implants number, (4) FWD freq
		vector<unsigned int> parameters;
		parameters.emplace_back(position_vector[i]);
		parameters.emplace_back(wobble_vector[i]);
		parameters.emplace_back(n_oligo_vector[i]);
		parameters.emplace_back(freq_strand_plus_vector[i]);
		
		//Insert the matrix and the parameters into a map (Ordered by the first parameters value --> by positions)
		jaspar_map.insert({parameters,matrix});

		//Clear the vector to avoid interferences in the next cycle
		parameters.clear();

	}	
}	

//******************************************************************************************************************
//Funtion to order the matrix name and tf name (strings) following the jaspar matrix values positioning into the map
void ordering_matrix_names(){
		
		//Making a map of position and name + tf strings --> to order them by positions
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

//***************************************************************************
//Function to generate a random distribution of FDW and REV matrix implanting
void freq_strand_plus_matrix(map<vector<unsigned int>, vector<vector<unsigned int>>>& jaspar_map){
	
	//For each matrix a vector of shuffled boolean (of implants number size) is created
	for(map<vector<unsigned int>, vector<vector<unsigned int>>>::iterator it = jaspar_map.begin(); it != jaspar_map.end(); it++){
		
		//Initializing a random device of mt19937 type
		mt19937 eng{random_device{}()};

		//Exact number of FWD implanting
		unsigned int number_plus = (it->first[3]*it->first[2]/100);
		
		//For each implants to be done
		for(unsigned int j=1; j<=it->first[2]; j++){
			
			if(j<=number_plus){

				plus_minus.emplace_back(0);
			}
			else{
				plus_minus.emplace_back(1);
			}
		}

		//Random shuffling of the 0/1 plus and minus vector and save it into a matrix of boolean (plus minus matrix)
		shuffle(plus_minus.begin(), plus_minus.end(),eng);	
		plus_minus_matrix.emplace_back(plus_minus);
		plus_minus.clear();	
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

void implanting_class::oligo_creation(map<vector<unsigned int>,vector<vector<unsigned int>>>::iterator it, int n){

		unsigned int strand;
		string oligo;

		for(unsigned int j=0; j<it->first[2]; j++){

			strand = plus_minus_matrix[n][j];
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

//****************************************************************************************??
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
			oligo_creation(it,j);				//creating a vector of random oligos coming from matrix frequences -> the size of oligos vector follows the input n_oligo_vector values for each matrix
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
				
				unsigned int half;

				if(it->second[0].size()%2 == 0){
					half = it->second[0].size()/2;
				}
				else{
					half = it->second[0].size()/2 +1;
				}

				implant_it = multifasta_map_implanted.find(unique_rnd[i]);	//find in multifasta map the unique_rnd[i](random number) sequence
				implant_it->second.replace((index_vec[i]-half), it->second[0].size(), oligo_vector[i]);		//Implant the oligo generated from jaspar in the right position in the right sequence

				cout << "Implanted string " << oligo_vector[i] << " in sequence number " << unique_rnd[i] << " centered in position " << index_vec[i] << "." << endl;
			
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
/*
//FUNCTION FOR BED FILE CLASS

void bed_class::read_line(string line){

    //Split the line word by word and extract chromosome coordinates (chr, start, end)
    istringstream mystream(line);
    mystream >> chr_coord >> start_coord >> end_coord;

}

//Flag control function: start coordinates must be < then end coordinates
void bed_class::flag_control( unsigned int start,  unsigned int end){

    //if start coordinates are >  end coordinates flag is setted to 0 --> WARNING printed to warn users
    if(start > end){

        flag = 0;
    }

    else{
        flag = 1;
    }
}

void bed_class::centering_function ( unsigned int start,  unsigned int end, unsigned int half_length, const unsigned int overhead){

    unsigned int center = (start + end)/2;

    //No overhead for start coordinates but overhead added to end coordinates
    start_coord = center - half_length;
    end_coord = center + half_length +overhead;

}

//Extract sequence function: Extract, from Twobit hg38 genome, the DNA sequence with (chr, start, end) coordinates extracted from Bed line
void bed_class::extract_seq(TwoBit* tb, unsigned int n_line){

    //CONTROL: if flag is 1 means that the current line has starting coordinate > end coordinate, so it is correct
    if(flag == 1){
        
        string chrom = chr_coord;

        //Extract the sequence from the object with the twobit_sequence function
        sequence = twobit_sequence(tb,chrom.c_str(),start_coord,end_coord-1);
    }

    //if flag is not 1 means that the current line has starting coordinate < end coordinate: PRINT WARNING!
    else {
        cerr << "WARNING: the line " << n_line <<" is omitted because starting coordinates > end coordinates, please check your BED file!" << "\n";
    }
}

void coordinator_class::GEP_creation(vector<bed_class> &GEP){

    cout << "\n- [1] Extract bed coordinate sequences from reference genome  \n";

    ifstream in(BED_FILE);
    TwoBit * tb;

    //Opening 2Bit file with twobit_open function from andrelmartens code and saved in tb variable
    tb = twobit_open(TWOBIT_FILE.c_str());

    string line;

    //Line counter initialization
    unsigned int n_line = 1;

    //For each line in BED file create a bed class called bed_line in which chr, start and end coordinate are ridden, centered and saved
    //Then the Fasta sequence is extracted from Twobit genome following the coordinated and saved into a string variable
    while(getline(in,line)){


        //if line is empty or commented --> continue
        if(line.empty() || line[0] == '#'){

            continue;
        }

        bed_class bed_line(line,tb, n_line);

        //For each line a bed class is created --> All the bed classes are saved in GEP vector (vector og bed class)
        GEP.emplace_back(bed_line);

        n_line = n_line + 1;

    }
}

*/
/////////////////////////////////////// PARSER ///////////////////////////////////

void command_line_parser(int argc, char** argv){
	
	const char* const short_opts = "hl:n:b:t:c:j:o:p:w:f:";

	//Specifying the expected options
	const option long_opts[] ={
		{"help",      no_argument, nullptr,  'h' },
		{"bed",	required_argument, nullptr, 'b'},         //PROVA
		{"twobit",	required_argument, nullptr, 't'},         //PROVA
		{"length",      required_argument, nullptr,  'l' },
		{"jaspar",   required_argument, nullptr,  'j' },
		{"nseq",   required_argument, nullptr,  'n' },
		{"cycles",   required_argument, nullptr,  'c' },
		{"oligop",   required_argument, nullptr,  'o' },
		{"position",   required_argument, nullptr,  'p' },
		{"wobble",   required_argument, nullptr,  'w' },
		{"freq",   required_argument, nullptr,  'f' },
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
			case 'b' : BED_FILE = string(optarg);
                                   is_file_exist(BED_FILE, "--bed || -b");	
				   break;
			case 't' : TWOBIT_FILE = string(optarg);		
                                   is_file_exist(TWOBIT_FILE, "--twobit || -t");
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
			case 'f' : freq_strand_plus = string(optarg); 
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
	cerr << "\n --jaspar || -j <JASPAR_FILE_1 JASPAR_FILE_2 ... JASPAR_FILE_N> to import the jaspar matrices from which the oligos will be generated." << endl;
	cerr << "\n --bed || -b <file_bed>: input bed file" << endl;	
	cerr << "\n --twobit || -t <file_twobit>: input TwoBit file" << endl;
	cerr << "\n --length || -l <number> to insert the length of Multifasta sequences (DEFAULT: 500)." << endl;
	cerr << "\n --nseq || -n <number> to insert the number of Multifasta sequences (DEFAULT: 200)." << endl;
	cerr << "\n --oligop || -o <n1,n2,...,nN> to insert the percentage of sequences in which the oligos will be implanted (DEFAULT: 0%) ---- NB: The number of sequences extracted from the percentage will be rounded down." << endl;
	cerr << "\n --position || -p <n1,n2,...,nN> to insert the position in Multifasta sequences where you want to implant the oligos." << endl;
	cerr << "\n --wobble || -w <n1,n2,...,nN> to set the wobble parameter for every implanting position. The implanting position, for every oligo, will be randomly choosen between p-w and p+w interval. (DEFAULT: 0)" << endl;
	cerr << "\n --cycles || -c <number> to choose how many random multifasta files (and implanted also) this tool will produce. (DEFAULT: 1)" << endl;
	cerr << "\n --freq || -f <number> to choose the percentage of fwd/rev strand to select for the implants. (DEFAULT: 50)" << endl;
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

