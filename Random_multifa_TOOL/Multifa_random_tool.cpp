#include "Multifa_random_tool.h"

/////////////////// GLOBAL VARIABLES ////////////////////////////////////////////////////////////////

vector<bed_class> GEP;
string var;
map<unsigned int, string> multiBED_map;
map<unsigned int, string> BED_map;
vector<string> header;
vector<unsigned int> unique_rnd;
unsigned int half_length = 250;
//unsigned int overhead = 25;
unsigned int length = 500;
unsigned int n_seq = 0;
string JASPAR_F;
string TWOBIT_FILE;
string oligo_perc = "0";
string position;
string wobble = "0";
vector<string> BED_FILE_vector;
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
map<unsigned int, string> pre_multiBED_map;

/////////////////////// MAIN FUNCTION //////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]){

	//Parser function to handle input arguments
	command_line_parser(argc, argv);

	//If any implanting position is inserted read all the other inputs
	if(position.size() != 0){  

		read_input();	
	}
	else{
		if (BED_FILE.size() != 0){
		
			bed_class_creation(GEP);

			//Storage of variables useful for fasta file (header and sequences)
			GEP_parameters(GEP);

			//Now we create a bed map like the multifasta map
			multiBED_map_creation(GEP);
		}
	}
	
	//For evry cycle -c inserted as input do inplanting_cycle function --> creating a random multifa + implanting (if any implaning position is inserted as input)
	for(unsigned int i=0; i<cycles; i++){
		implanting_cycle(i);
	}
	
}
//////////////////////// INPUT READING AND CHECKING ////////////////////////////////////////////////////////////

//Function to read the input passed on command line and fill the vector to store the input information
void read_input(){

	if (BED_FILE.size() != 0){
		
		bed_class_creation(GEP);

		//Storage of variables useful for fasta file (header and sequences)
		GEP_parameters(GEP);

		//Now we create a bed map like the multifasta map
		multiBED_map_creation(GEP);
	}
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

	if(BED_FILE.size() == 0){
		for(unsigned int i=0; i<n_oligo_vector.size(); i++){
		
			if(n_oligo_vector[i] > 100 || n_oligo_vector[i] < 0){

				cout << "ERROR: -o parameter need to be from 0 to 100!" << endl;
				exit(1);
			}
			//Transforming the implanting percentage into the real number of implants that need to be done

			n_oligo_vector[i] = (n_oligo_vector[i] * n_seq)/100;
		}
	}
	
	else{
		if (n_seq == 0){
			for(unsigned int i=0; i<n_oligo_vector.size(); i++){
		
				if(n_oligo_vector[i] > 100 || n_oligo_vector[i] < 0){

					cout << "ERROR: -o parameter need to be from 0 to 100!" << endl;
					exit(1);
				}
				//Transforming the implanting percentage into the real number of implants that need to be done

				n_oligo_vector[i] = (n_oligo_vector[i] * header.size())/100;
			}
		}
		else{
			for(unsigned int i=0; i<n_oligo_vector.size(); i++){
		
				if(n_oligo_vector[i] > 100 || n_oligo_vector[i] < 0){

					cout << "ERROR: -o parameter need to be from 0 to 100!" << endl;
					exit(1);
				}
				//Transforming the implanting percentage into the real number of implants that need to be done

				n_oligo_vector[i] = (n_oligo_vector[i] * n_seq)/100;
			}
		}
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

		cerr << "\nWARNING: No frequence inserted, the tool will set all the frequences on 50%" << endl;
		
		//Set the frequency at 50% in frequency (strand plus) vector
		for(unsigned int i=0; i<JASPAR_FILE_vector.size(); i++){

			freq_strand_plus_vector.emplace_back(50);
		}
	}
	
	//Else if number of frequences insterted are more than 0 but less than the matrices that are going to be implanted
	else if(freq_number > 0 && freq_number < JASPAR_FILE_vector.size()){
		
		cerr << "\nWARNING: Frequence number inserted is less than matrix inserted, the tool will set on 50% the frequences not specified" << endl;
		
		//Set on 50% the frequences not specified until the size of freq are equal to jaspar matrices number
		while(freq_strand_plus_vector.size() != JASPAR_FILE_vector.size()){

			freq_strand_plus_vector.emplace_back(50);
			}
	}
	
	//Else if there are more frequencies than the matrices in input --> print a warning and delete the excess frequencies
	else if(freq_number > JASPAR_FILE_vector.size()){

		cerr << "\nWARNING: Frequence number inserted is more than matrix inserted, Please check your input parameters " << endl;
	
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
        
		debug_implanting();
		
		//Controlling that the number of Jaspar matrices in input are equal to number of -p (implanting position) in input. If the control is positive, the tool associate each matirx to an implanting position --> the map <position,jaspar_name> starts to be filled	
		check_input();	

		//Function to associate each Jaspar matrix to his parameters correctly
		filling_jaspar_map(jaspar_map,matrix_name,tf_name);

		//Function to order matrix name and tf name to mantain their association with the matrix
		ordering_matrix_names();

		//Function to generate a random distribution of FDW and REV matrix implanting
		freq_strand_plus_matrix(jaspar_map);
        if (BED_FILE.size() == 0){
			if (n_seq == 0){
				cerr << "WARNING! The -n (--nseq) parameter is not present, the program will use a default number of sequences of 200!" << endl;
				n_seq = 200;
			}
        	multifasta_class MULTIFA(length,n_seq,i);
			//if the jaspar input are more then 1
			if(position_vector.size() > 1){
			
				//checking if the implanting positions don't bring to any overlap
				check_overlapping(jaspar_map);
			}

			//checking if -p implanting positions don't bring the oligos to exceed from the sequences length 
			//check_exceeding(jaspar_map);
			
			
			//implanting class constructor calling --> The map composed by matrix + parameters, the random multifa map, and the current cycle number are passed to constructor
			implanting_class (jaspar_map, MULTIFA.multifasta_map,i);
		
			matrix_n.clear();
			matrix_tf.clear();
 
        }
        else{

			if (n_seq == 0){
				multiBED_map.insert(pre_multiBED_map.begin(), pre_multiBED_map.end());	
			}
			else{
				casual_map_filtering(pre_multiBED_map);
			}
			implanting_class (jaspar_map, multiBED_map,i);

			
			
			if(position_vector.size() > 1){	
				//checking if the implanting positions don't bring to any overlap
				check_overlapping(jaspar_map);
			}
					
			//implanting class constructor calling --> The map composed by matrix + parameters, the random multifa map, and the current cycle number are passed to constructor
			matrix_n.clear();
			matrix_tf.clear();
			multiBED_map.clear();
        }
	
	}
	
	//Else if there are no implanting positions passed as input
	else{
		if (BED_FILE.size() == 0){
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

			multifasta_outfile(MULTIFA.multifasta_map, "random_multifa_"+to_string(i+1)+".fasta");
		}
		else{
			//Printing of warning messages
			
			cout << "\nWARNING: No Jaspar matrices and implanting position given as input." << endl;	
			if (n_seq == 0){
				cout <<"\n Generating a random BED file of " << header.size() << " sequences of " << length << " bases length\n";
			}
			else{
				cout <<"\n Generating a random BED file of " << n_seq << " sequences of " << length << " bases length\n";
			}
			if (n_seq == 0){
				multiBED_map.insert(pre_multiBED_map.begin(), pre_multiBED_map.end());	
			}
			else{
				casual_map_filtering(pre_multiBED_map);
			}
			
			multifasta_outfile(multiBED_map, "BED_"+to_string(i+1)+".fasta");
		
		}
	}

}


//Checking if Jaspar matreces number is equal to positions number inserted as input --> if not print an error and exit
void check_input(){

	if(JASPAR_FILE_vector.size() != position_vector.size()){

		cerr << "\nERROR: Please insert the rigth number of Jaspar files and implanting position.\nTheir number need to be equal to make a correct implanting!" << endl;
		exit(1);
	}
}

void bed_class_creation(vector<bed_class> &GEP){
	
	ifstream in(BED_FILE);
		
	TwoBit * tb;

	tb = twobit_open(TWOBIT_FILE.c_str());

	string line;

	unsigned int n_line = 1;
	
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

void bed_class::read_line(string line){ 

	//Split the line word by word and extract chromosome coordinates (chr, start, end)
	istringstream mystream(line);
	mystream >> chr_coord >> start_coord >> end_coord;		

}

//Flag control function: start coordinates must be < then end coordinates
void bed_class::flag_control(){

	//if start coordinates are >  end coordinates flag is setted to 0 --> WARNING printed to warn users
	switch (start_coord > end_coord){
		case 1:
			flag = 0;
			break;
		case 0:
			flag = 1;
			break;
	}
}

void bed_class::centering_function(){

	unsigned int center = (start_coord + end_coord)/2;	
	half_length = length/2;					
	//No overhead for start coordinates but overhead added to end coordinates
	start_coord = center - half_length;
	end_coord = center + half_length;

}

void bed_class::extract_seq(TwoBit* tb, unsigned int n_line){
	switch (flag){
		case 1:
			//Extract the sequence from the object with the twobit_sequence function
			sequence = twobit_sequence(tb,chr_coord.c_str(),start_coord,end_coord-1);
			break;
		case 0: 
			//if flag is not 1 means that the current line has starting coordinate > end coordinate: PRINT WARNING!		
			cerr << "ERROR: the line " << n_line <<" is omitted because starting coordinates > end coordinates, please check your BED file!" << "\n";
			exit(1);
	}
}

void GEP_parameters(vector<bed_class> GEP){

	for(unsigned int element = 0; element < GEP.size(); element++){
		string chrom = GEP[element].return_chr_coord();
		unsigned int start = GEP[element].return_start_coord();
		unsigned int end = GEP[element].return_end_coord();
		string seq = GEP[element].return_sequence();
		stringstream ss;
    	ss << ">" << chrom << "-" << start << ":" << end;
		string s = ss.str();
		header.emplace_back(s);

	}
	if (n_seq > header.size()){
		cerr << "\nWARNING! The -n parameter is higher than the sequences in BED file, the programm will consider all the sequences\n" << endl;
	}
}

string bed_class::return_chr_coord(){

	return chr_coord;
}

unsigned int bed_class::return_start_coord(){

	return start_coord;
}

unsigned int bed_class::return_end_coord(){

	return end_coord;
}

string bed_class::return_sequence(){

	return sequence;
}

void multiBED_map_creation(vector<bed_class> GEP){
	//For each sequence
    for(unsigned int j=0; j<(GEP.size()); j++){
  
        pre_multiBED_map.insert({j, GEP[j].return_sequence()});
    }
	
}
  
void casual_map_filtering(map<unsigned int, string> pre_multiBED_map){
  
	vector<unsigned int> shuffled_vector;
	mt19937 eng{random_device{}()};
  
    for(unsigned int i=0; i<header.size(); i++){
  
        shuffled_vector.emplace_back(i);
    }
  
    shuffle(shuffled_vector.begin(), shuffled_vector.end(), eng);
    map<unsigned int, string>::iterator itr = pre_multiBED_map.begin();
  
    for(unsigned int i=0; i<n_seq; i++){
  
        itr = pre_multiBED_map.find(shuffled_vector[i]);
        multiBED_map.insert({itr->first, itr->second});
    }
	shuffled_vector.clear();
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

		//Filling the parameters vector with (0) position, (1) wobble, (2) implants number, (3) FWD freq
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

//////////////////////////////////////// MULTIFASTA RANDOM GENERATION ///////////////////////////////////////

//Function to create the random multifasta sequences
void multifasta_class::multifasta_map_creation(){

	string sequence;
	
	//For each sequence
	for(unsigned int j=0; j<n_seq; j++){

		//Clearing the string sequence to avoid interference with the one of the previous cycle
		sequence.clear();
		
		//For each base from 0 to seq_length
		for(unsigned int i=0; i<length; i++){

			//A random number between 0 and 3 is generated (A=0, C=1, G=2, T=3)
			unsigned int random_int = random_number(0,3);

			//Conversion from the number randomly generated to an effective base (0123 -> ACGT)
			char base = from_n_to_base(random_int);

			//Concatenate the base just created to the under construction sequence
			sequence = sequence + base;
		}
		
		//Insert the the pair number of current seq + sequence (for ex 1-ATCGGGA, 2-CCTCCA...) in a map
		multifasta_map.insert(pair<unsigned int,string>(j+1, sequence));
	}
	
}

//Function to generate a random number between a start number/end number range (in this case between 0 and 3)
unsigned int random_number(unsigned int range_begin, unsigned int range_end){
	
	//Unsing the mt19937_64 random generator (for more search it on internet)
        mt19937_64 generator (clock());	

	//A distribution of values (between 0 and 3) is generated.
        uniform_int_distribution<unsigned int> dis(range_begin, range_end);	

	//From this number distribution a random number is generated and then returned
        unsigned int r_number = dis(generator);;
	return 	r_number;
}

//Function to convert a random number from 0 to 3 into a base (If a number < 0 or > 3 is passed an error is generated)
char from_n_to_base (unsigned int n){

	char base;
			switch(n){

				case 0:

					base = 'A';
					break;

				case 1:

					base = 'C';
					break;

				case 2:

					base = 'G';
					break;

				case 3:

					base = 'T';
					break;

				default:
					
					cerr << "ERROR: Generation of random sequence failed. Please check the random generation code.!" << endl;
					exit(1);
					break;
			}
			return base;
}

//Function to read JASPAR PWM file, extract values and create a matrix class
vector<vector<unsigned int>> read_JASPAR(string JASPAR_FILE, string& matrix_name, string& tf_name){		

	ifstream file(JASPAR_FILE);
	string line;
	vector<vector<unsigned int>> matrix;
	
	//For each line of the JASPAR file	
	while(getline(file,line)){

		//If the line start with '>' character save the words into matrix_name string and into tf_name string
		if(line[0]=='>'){

			istringstream mystream(line);			
			mystream >> matrix_name >> tf_name;
		}

		//If the line does not start with '>' delete the '[',']' and 'A,T,C,G' characters, extract the scores and save them into a scores matrix
		else{	

			//Deleting from line the first character (A,T,C,G), the '[' and the ']' characters
			line.erase(0,line.find('[') +1);
			line.erase(line.find(']'));

			vector<unsigned int> scores_line;
			istringstream mystream(line);

			//For each words (number) in line put the current number in num variables
			for (unsigned int num; mystream >> num;){

				scores_line.emplace_back(num);
			}

			matrix.emplace_back(scores_line);			
		}

	}
	file.close();
	return matrix;
}

/////////////////////////////////// OVERLAP & EXCEED CHECKING ////////////////////////////////////////////////////

//Check, if more than 1 implanting pos -p is inserted, if the implanting positions, wobbles and jaspar lengths, can bring to an implant overlapping
void check_overlapping(map<vector<unsigned int>,vector<vector<unsigned int>>> jaspar_map){
	
	//Since the jaspars are ordered in the map by their implanting positions, it is possible to compare a jaspar only to the next jaspar to understand if an overlap can happens
	//Iterator before for the current jaspar in the map | Iterator after for the next jaspar in the map
	map<vector<unsigned int>, vector<vector<unsigned int>>>::iterator it_before = jaspar_map.begin();
	map<vector<unsigned int>, vector<vector<unsigned int>>>::iterator it_after = ++jaspar_map.begin();

	//From all the position -p in the position vector (-1 because the last can't be compared to any position) 
	for(unsigned int i=0; i<(position_vector.size()-1); i++, it_after++, it_before++){
		
		//If "after" position is < or = to "before" position + "before" jaspar length --> wrong positions inserted 
		if(it_after->first[0] <= (it_before->first[0] + it_before->second[0].size())){

			cerr << "\nERROR: The oligos coming from matrix " << matrix_n[i] << " and " << matrix_n[i+1] << " that you are trying to implant overlap ---> wrong positions inserted!"<< endl;
			exit(1);
			}

		//If "after" position + "after" wobble is < or = to "before" position + "before" wobble + "before" jaspar length --> wrong wobble inserted 
		else if((it_after->first[0] - it_after->first[1]) <= (it_before->first[0] + it_before->second[0].size() + it_before->first[1])){

			cerr << "\nERROR: The oligos coming from matrix " << matrix_n[i] << " and " << matrix_n[i+1] << " that you are trying to implant overlap ---> wrong wobbles inserted!"<< endl;
			exit(1);
		}
	}
}

/////////////////////////////// IMPLANTING /////////////////////////////////////////////////////////////////////////////

//******************************************************************************************************
//Function to implant the oligo generated from jaspars into the random sequences (in correct position)
void implanting_class::implanting_oligo(map<vector<unsigned int>, vector<vector<unsigned int>>> jaspar_map){

		int j=0;
		
		//For each jaspar inserted as input
		for(map<vector<unsigned int>,vector<vector<unsigned int>>>::iterator jaspar_it = jaspar_map.begin(); jaspar_it!=jaspar_map.end(); jaspar_it++,j++){
			
			//generating a vector of unique random numbers from 1 to n_seq (length n_seq)
			unique_random_generator();

			//creating a vector of random oligos coming from jaspar matrix scores -> the number of oligo generated follows the input n_oligo_vector values for each jaspar matrix
			oligo_creation(jaspar_it,j);
			
			//An empty vector of index (new implanting position -p generated by wobble addition) is initialized
			vector<unsigned int> index_vec;
			
			//From 0 to -w wobble value inserted as input
			for(unsigned int w=0; w < jaspar_it->first[2]; w++){
				
				//An index (new implanting position) is generated randomly following -p and -w inputs
				//For example if p=150 and w=3 the new implating position is a number randomly generated between 147 (150-3) and 153 (150+3)
				int index = random_number((jaspar_it->first[0] - jaspar_it->first[1]), (jaspar_it->first[0] + jaspar_it->first[1]));				
				//index_vec is vector of new idexes (implanting positions) with wobble variation added
				index_vec.emplace_back(index);
			}

			//Some control printing in cout to check if the implanting is performing correctly
			cout << "\nStarting implant position: "<< jaspar_it->first[0] << endl;
			cout << "Wobble inserted: " << jaspar_it->first[1] << endl<<endl;
			cout << matrix_n[j] << ": " << matrix_tf[j] << endl;
			print_debug_matrix(jaspar_it->second);
			cout << endl;
			
			check_exceeding(jaspar_map);
			
			//An iterator for the multifasta sequences map is initialized
			map<unsigned int, string>::iterator multifa_it;
			//From 0 to -o (number of oligos that must be implanted)
			for(unsigned int i=0; i<jaspar_it->first[2]; i++){ 	
				
				//Here the half length of jaspar matrix is calculated --> to center the oligos in implanting position -p (and not make the implants starting in -p)
				unsigned int half;

				//If the jaspar length is even
				if(jaspar_it->second[0].size()%2 == 0){
					half = jaspar_it->second[0].size()/2;
				}
				//If the jaspar length is odd
				else{
					half = jaspar_it->second[0].size()/2 +1;
				}

				//Select the correct sequence (in multifa sequences map) for the implant, following the shuffling unique_rnd vector previously generated)
				multifa_it = multifasta_map_implanted.find(unique_rnd[i]);

				//Replace in the correct sequence the bases from starting implant pos (implanting index - half_matrix length) to end implant pos (start pos + jaspar length) with the current oligo in oligo vector
				multifa_it->second.replace((index_vec[i]-half), jaspar_it->second[0].size(), oligo_vector[i]);

				//Debug printing to control if the implants are performed correctly
				if(BED_FILE.size() == 0){
					cout << "Implanted string " << oligo_vector[i] << " in sequence number " << unique_rnd[i] << " centered in position " << index_vec[i] << "." << endl;
				}
				else{
					cout << "Implanted string " << oligo_vector[i] << " in sequence " << header[unique_rnd[i]] << " centered in position " << index_vec[i] << "." << endl;
				}
			}

			cout << endl << "---------------------------------------------------------------------------------------" << endl;	
			
			//Clearing of index vector, unique random vector and oligo vector for next jaspar implanting
			index_vec.clear();
			unique_rnd.clear();
			oligo_vector.clear();
		}
}

//Check if there are some possible out-of-sequence implant --> analyzing the wobble -w, implanting pos -p and jaspar length
void check_exceeding(map<vector<unsigned int>,vector<vector<unsigned int>>> jaspar_map){
	
	//For each jaspar inserted as input
	for(map<vector<unsigned int>,vector<vector<unsigned int>>>::iterator jaspar_it = jaspar_map.begin(); jaspar_it != jaspar_map.end(); jaspar_it++){

		unsigned int matrix_half_length = jaspar_it->second[0].size()/2;

		for(unsigned int i = 0; i <= position_vector.size(); i++){
			
			if ((matrix_half_length + jaspar_it -> first[1]) >= jaspar_it->first[0] || (jaspar_it->first[0] + jaspar_it->first[1] + matrix_half_length) > length){
				cerr << "\nERROR: You've inserted position " << jaspar_it->first[0] << ", this position implies that the implant exceed from the sequence of length " << length << "!\nPlease check your implanting position or length of the sequence!" << endl;
				exit(1);
			}

		}

	}
}
//Function to generate a vector of number from 1 to number of seq --> then the function shuffles the vector to generate a 0-n_seq number distribution (for ex: from 1-2-3-4-5-6 --> 3-5-1-4-2-6)
void unique_random_generator(){

	//using the mt19937 random eng generator fot shuffling
	mt19937 eng{random_device{}()};
	
	//Put the number from 1 to n_seq in an ordered vector
	if(BED_FILE.size() == 0){
		for(unsigned int i=1; i<=n_seq; i++){

			unique_rnd.emplace_back(i);
		}

		//shuffling the vector
		shuffle(unique_rnd.begin(), unique_rnd.end(), eng);
	}
	else{
		for(map<unsigned int, string>::iterator itr = multiBED_map.begin(); itr != multiBED_map.end(); itr++){
			
			unique_rnd.emplace_back(itr -> first);
		}

		//shuffling the vector
		shuffle(unique_rnd.begin(), unique_rnd.end(), eng);
	}
}

//*************************************************************************************************
//Generating (for each matrix) a vector containing -o oligos coming from jaspar scores frequences
void implanting_class::oligo_creation(map<vector<unsigned int>,vector<vector<unsigned int>>>::iterator it, int n){
		
		unsigned int strand;
		string oligo;
		
		//from 0 to -o parameter (for each oligo that must be generated)
		for(unsigned int j=0; j<it->first[2]; j++){

			//Deciding the strand where the oligo is going to be implant scrolling the plus_minus_matrix previously generated (followig the -f frequences parameter) --> (for ex a line of the matrix can be: 0010001100 if -f is 70%)
			//(for ex: if "+" is selected --> implanting of ATTCAA | if "-" selected --> implanting of TTGAAT)
			strand = plus_minus_matrix[n][j];
			
			//Clearing oligo vector (filled in previous cycle) to avoid interferences
			oligo.clear();

			//For each jaspar matrix column
			for (unsigned int i = 0; i < it->second[0].size(); i++) {

				//Sum the column scores
				unsigned int col_sum = it->second[0][i] + it->second[1][i] + it->second[2][i]+ it->second[3][i];

			//A random number from 1 to col sum is generated (using random_number function)
			unsigned int random_score = random_number(1,col_sum);
			
			//If the number is from 1 to A-score --> A is assigned to oligo
			if(random_score <= it->second[0][i]){
				oligo = oligo + 'A';
			}
			
			//Else if the number is from A-score to A-score + C-score --> C is assigned to oligo
			else if(random_score > it->second[0][i] && random_score <= (it->second[0][i] + it->second[1][i])){
				oligo = oligo + 'C';
			}

			//Else if the number is from A-score + C-score to A-score + C-score + G-score --> G is assigned to oligo
			else if(random_score > (it->second[0][i]+it->second[1][i]) && random_score <= (it->second[0][i] + it->second[1][i] + it->second[2][i])){
				oligo = oligo + 'G';
			}
			
			//Else --> T is assigned to oligo
			else{
				oligo = oligo + 'T';
			}
		}
		
		//If strand is equal to 0 (means that the strand selected is "+") --> put the oligo in oligo vector
		if(strand == 0){
			oligo_vector.emplace_back(oligo);
		}

		//Else means that the strand is 1 (strand selected is "-") --> Generate the reverse complement and put the reverse complement into the oligo vector
		else{
			oligo = reverse_complement(oligo); 
			oligo_vector.emplace_back(oligo);
		}

	}
}

//Function to generate, from a string in input, its reverse complement (ATAAC --> GTTAT)
string reverse_complement(string oligo){

	//Empty string "reverse" is initialized
	string reverse_oligo;

	//For each base in input string
	for(unsigned int i=0; i<oligo.size(); i++){
		
		char base;
		base = oligo[i];
		switch (base) {
			
			//Using the append function to insert a char into a string
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
	
	//Function reverse to reverse the string just created (form TATTG --> GTTAT)
	reverse(reverse_oligo.begin(), reverse_oligo.end());

	//Return the reverse complement
	return reverse_oligo;
}

/////////////////////////////////////// DEBUGGING ////////////////////////////////////////////////////////////////////////
void debug_implanting(){
	if(BED_FILE.size() == 0){
    	//Debug for a meaningful output -> here you can check if the inputs inserted are correct
		cout << "\nThe number of JASPAR matrices is: " << JASPAR_FILE_vector.size() << endl;
		cout << "The number of implanting position is: " << position_vector.size() << endl;
		cout << "The length of multifasta random sequences is: " << length << endl;
		cout << "The number of multifasta random sequences generated is: " << n_seq << endl;
		cout << "The Jaspar matrices inserted as input are: ";
	}
	else{
		cout << "\nThe number of JASPAR matrices is: " << JASPAR_FILE_vector.size() << endl;
		cout << "The number of implanting position is: " << position_vector.size() << endl;
		cout << "The length of BED sequences is: " << length << endl;
		if (n_seq != 0){
			if (n_seq >= header.size()){
				cout << "You choose to pick all the sequences for your BED file." << endl;	
			}
			else{
				cout << "You choose to pick " << n_seq << " sequences for your BED file." << endl;
			}
		}
		cout << "The number of BED sequences is: " << header.size() << endl;
		cout << "The BED files inserted as input are: ";
		for(unsigned int i=0; i<BED_FILE_vector.size(); i++){
			cout <<BED_FILE_vector[i] << " ";
		}
		cout << endl;
		cout << "The number of BED file is: " << BED_FILE_vector.size() << endl; 
		cout << "The Jaspar matrices inserted as input are: ";		
	}
	
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
}

void implanting_class::print_debug_matrix(vector<vector<unsigned int>> matrix){		//Print matrix function

	for(unsigned int i=0; i < matrix.size(); i++){
		for(unsigned int j=0; j<matrix[i].size(); j++){

			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
}

//Function to print the random multifasta sequences generated into a file called random_multifasta_(i).fasta
void multifasta_outfile(map<unsigned int,string> multifasta_map, string filename){

	ofstream outfile;
	outfile.open(filename);
	unsigned int i = 0;
	//For each sequence in multifasta map
	for(map<unsigned int,string>::iterator it = multifasta_map.begin(); it != multifasta_map.end(); it++){
		if(BED_FILE.size() == 0){
			outfile << ">random multifasta sequence number " + to_string(it->first) + " containing "+ to_string(length) +" bases:"<<endl;

			//Printing the sequence
			outfile << it->second << endl;
			outfile << endl;
		}
		else{
			outfile << header[it->first] << endl;
			//Printing the sequence
			outfile << it->second << endl;
			outfile << endl;	
			i++;
		}

	}
	outfile.close();
}

/////////////////////////////////////// PARSER //////////////////////////////////////////////////////////////////

void command_line_parser(int argc, char** argv){
	
	const char* const short_opts = "hl:n:b:t:c:j:o:p:w:f:";

	//Specifying the expected options
	const option long_opts[] ={
		{"help",      no_argument, nullptr,  'h' },
		{"length",      required_argument, nullptr,  'l' },
        {"bed",	required_argument, nullptr, 'b'}, 
		{"twobit", required_argument, nullptr, 't'},
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
			case 'l' : length = stoi(optarg); 
				   break;
			case 'b' : BED_FILE = string(optarg);
                    is_file_exist(BED_FILE, "--bed || -b");
					BED_FILE_vector.emplace_back(BED_FILE);
					for (;optind < argc && *argv[optind] != '-';optind++){
						BED_FILE = (string(argv[optind]));
						is_file_exist(BED_FILE, ("--bed || -b one of files do not exist"));
						BED_FILE_vector.emplace_back(BED_FILE);
				   	}	
			        break;
            case 't' : TWOBIT_FILE = string(optarg);		
                    is_file_exist(TWOBIT_FILE, "--twobit || -t");
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

////////////////////////////////// DISPLAY HELP /////////////////////////////////////////////////////////////

void display_help(){

	cerr << "\n --help || -h show this message" << endl;
	cerr << "\n --jaspar || -j <JASPAR_FILE_1 JASPAR_FILE_2 ... JASPAR_FILE_N> to import the jaspar matrices from which the oligos will be generated." << endl;
	cerr << "\n --length || -l <number> to insert the length of Multifasta sequences (DEFAULT: 500)." << endl;
    cerr << "\n --bed || -b <file_bed>: input bed file" << endl;
	cerr << "\n --twobit || -t <twobit_file> to import the genome. (This parameter is needed only if is present BED file)" << endl;	
	cerr << "\n --nseq || -n <number> to insert the number of Multifasta sequences (DEFAULT: 200) or in the case of BED input the number of sequences you want to take into account (DEFAULT: all the seuqneces in BED file)." << endl;
	cerr << "\n --oligop || -o <n1,n2,...,nN> to insert the percentage of sequences in which the oligos will be implanted (DEFAULT: 0%) ---- NB: The number of sequences extracted from the percentage will be rounded down." << endl;
	cerr << "\n --position || -p <n1,n2,...,nN> to insert the position in Multifasta sequences where you want to implant the oligos." << endl;
	cerr << "\n --wobble || -w <n1,n2,...,nN> to set the wobble parameter for every implanting position. The implanting position, for every oligo, will be randomly choosen between p-w and p+w interval. (DEFAULT: 0)" << endl;
	cerr << "\n --cycles || -c <number> to choose how many random multifasta files (and implanted also) this tool will produce. (DEFAULT: 1)" << endl;
	cerr << "\n --freq || -f <number> to choose the percentage of fwd/rev strand to select for the implants. (DEFAULT: 50)" << endl;
	cerr << endl;
	exit(EXIT_SUCCESS);
}
