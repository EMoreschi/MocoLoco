#include "MocoLoco.h" 
//#include "Profiling.h"
#include <sys/resource.h>


int main(int argc, char *argv[]){
	Timer timer;
	//Instrumentor::Get().BeginSession("MocoLoco");
	{

		//If arguments number is 1 means that no input file has been inserted - display help
		if(argc == 1){

			display_help();
		}
		command_line_parser(argc, argv);
	

		GEP_path();
	
		return 0;
	}
	//Instrumentor::Get().EndSession();
}

//Function to analyse the RAM usage of the tool, it returns the maximum amount of memory allocated by the program
void RAM_usage(){
	int who = RUSAGE_SELF;
    struct rusage usage;
    int ret;
    ret = getrusage(who, &usage);
	cout << ret << endl;
    cout << endl << "Maximum resident set size: " << usage.ru_maxrss/1000 << " Mb" << endl << "User CPU time used: " << usage.ru_utime.tv_sec << " s" << endl << "System CPU time used: " << usage.ru_stime.tv_usec << " micros" << endl;

}

//Function to choose the pathway to follow. 2 input options:
//1) Bed-Twobit-Jaspar input
//2) Multifasta input
void  GEP_path(){
	//PROFILE_FUNCTION();

	//if the input is Bed-Twobit-Jaspar	
	if(MFASTA_FILE.size() == 0){	
		tb = twobit_open(TWOBIT_FILE.c_str()); 
		coordinator_class C; 
		twobit_close(tb);
		//Create a .fasta file to check if the coordinates and the sequences extracted are correct
		C.print_GEP(C.GEP);
		
		//Creating map class: input are GEP vector created from bed-twobit analysis, kmers, hamming distance
		map_class MAP(C.GEP,kmers,dist);
		
		RAM_usage();
	}

	//else if the input is a Multifasta file
	else{

		multifasta_class MULTIFA(MFASTA_FILE);

		//Creating a Map class: input are GEP vector created from multifasta file analysis, kmers, hamming distance
		map_class MAP(MULTIFA.GEP,kmers,dist);
		RAM_usage();
	}

}

void bed_class::read_line(string line){ 
	//PROFILE_FUNCTION();
	//Split the line word by word and extract chromosome coordinates (chr, start, end)
	istringstream mystream(line);
	mystream >> chr_coord >> start_coord >> end_coord;		
}
/*
//Flag control function: start coordinates must be < then end coordinates
void bed_class::flag_control( unsigned int start,  unsigned int end){
	//PROFILE_FUNCTION();
	//if start coordinates are >  end coordinates flag is setted to 0 --> WARNING printed to warn users
	if(start > end){

		flag = 0;
	}

	else{ 
		flag = 1;
	}
}
*/
void bed_class::centering_function (unsigned int start,  unsigned int end, int half_length, const unsigned int overhead){
	//PROFILE_FUNCTION();
	unsigned int center = (start + end)/2;						
	if(start > end){

		flag = 0;
	}

	else{ 
		flag = 1;
	}
	//No overhead for start coordinates but overhead added to end coordinates
	start_coord = center - half_length;
	end_coord = center + half_length +overhead;
}

//Extract sequence function: Extract, from Twobit hg38 genome, the DNA sequence with (chr, start, end) coordinates extracted from Bed line
void bed_class::extract_seq(TwoBit* tb, unsigned int n_line){
	//PROFILE_FUNCTION();
	//CONTROL: if flag is 1 means that the current line has starting coordinate > end coordinate, so it is correct
	if(flag == 1){	
		
		string chrom = chr_coord;

		//Extract the sequence from the object with the twobit_sequence function
		sequence = twobit_sequence(tb,chrom.c_str(),start_coord,end_coord-1);
	}

	//if flag is not 1 means that the current line has starting coordinate < end coordinate: PRINT WARNING!		
	else {		
		err = true;
		cerr << "WARNING: the line " << n_line << " (" << chr_coord << ":" << start_coord << "-" << end_coord << ")" << " is omitted because starting coordinates > end coordinates, please check your BED file!" << "\n";
	}
}

//Function useful to normalize matrix scores and adding a pseudocount to them
void matrix_class::matrix_normalization_pseudoc(vector<vector<double>> &matrix){  						//CHANGE: Possible use of reference
	//PROFILE_FUNCTION();
	//double normalized_score;

	//Calculate and save the scores column sum into a vector to perform a faster normalization step	
	vector<double> col_sum = find_col_sum(matrix);

	for (unsigned int i = 0; i < matrix.size(); i++) {

		vector<double> normalized_matrix_line;
		for (unsigned int j = 0; j < matrix[i].size(); j++){

			//normalized_score = matrix[i][j]/col_sum[j];
			normalized_matrix_line.emplace_back(matrix[i][j]/col_sum[j] + pseudoc);
		}

		norm_matrix.emplace_back(normalized_matrix_line);
	}
}

//Function which saves into a vector called col_sum all the score column sums --> This is made to perform the next Normalization step faster
vector<double> matrix_class::find_col_sum(vector<vector<double>> &matrix){								//CHANGE: Possible use of reference
	//PROFILE_FUNCTION();
	vector<double> col_sum;						
	double sum = 0;								
	for (unsigned int i = 0; i < matrix[0].size(); i++){	
		for (unsigned int j = 0; j < 4; j++){			
			
			sum += matrix[j][i];			
		}

		col_sum.emplace_back(sum);				
		sum = 0;						
	}

	return col_sum;
}

//Function to perform a second normalization on matrix scores (without a pseudocount addition)
void matrix_class::matrix_normalization(vector<vector<double>> &matrix){									//CHANGE: Possible use of reference
	//PROFILE_FUNCTION();
	//Calculate and save again the scores column sum into a vector to perform a faster normalization step
	vector<double> col_sum = find_col_sum(matrix);

	for (unsigned int i = 0; i < matrix.size(); i++) {

		for (unsigned int j = 0; j < matrix[i].size(); j++){

			//Substitution of first normalized values with new normalized ones
			norm_matrix[i][j] = matrix[i][j]/col_sum[j];
		}
	}
}

//Function to calculate, from the normalized matrix, the logarithmic values of the scores and creates a new matrix called matrix_log
void matrix_class::matrix_logarithmic(vector<vector<double>> &matrix){									//CHANGE: Possible use of reference
	//PROFILE_FUNCTION();
	for(unsigned int i=0; i < matrix.size(); i++){
		
		vector<double> log_matrix_line;
		//double log_scores;

		for(unsigned int j=0; j < norm_matrix[i].size(); j++){

			//log_scores = log(norm_matrix[i][j]);
			log_matrix_line.emplace_back(log(norm_matrix[i][j]));
		}

		matrix_log.emplace_back(log_matrix_line);
	}
}


//Function which return the Transposed matrix from a matrix in input
vector<vector<double>> matrix_class::reverse_matrix(vector<vector<double>> &matrix){		
	//PROFILE_FUNCTION();
	vector<vector<double>> rev_matrix = matrix;
	reverse(rev_matrix.begin(), rev_matrix.end());

	for (int i = 0; i < 4; i++) {

		reverse(rev_matrix[i].begin(), rev_matrix[i].end());
	}

	return rev_matrix;

}

//Finding the best and the worst score that an oligo can reach based on current JASPAR matrix
void oligo_class::find_minmax(vector<vector<double>> &matrix){
	//PROFILE_FUNCTION();
	//Extract the mins and the maxes from each columns, saved into vectors and their total sum will be the best and the worst score that an oligo can reach
	for(unsigned int i=0; i < matrix[0].size(); i++){

		vector<double> colum;		   	
		for(unsigned int j=0; j < matrix.size(); j++){

			colum.emplace_back(matrix[j][i]);
		}

		o_matrix_mins.emplace_back(*min_element(colum.begin(),colum.end()));
		o_matrix_maxes.emplace_back(*max_element(colum.begin(),colum.end()));
	}

	min_possible_score = accumulate(o_matrix_mins.begin(), o_matrix_mins.end(), 0.0);
	max_possible_score = accumulate(o_matrix_maxes.begin(), o_matrix_maxes.end(), 0.0);

}	

//Function to calculate the score of a general oligo against a JASPAR matrix
void oligo_class::shifting(vector<vector<double>> &matrix, string &sequence/*, unsigned int s_iterator*/){
	//PROFILE_FUNCTION();
	unsigned int max;
	max = sequence.size() - matrix[0].size();
	for (unsigned int s_iterator = 0; s_iterator <= max; s_iterator++){
		double sum_scores = 0;
		//For each oligo in the current sequence a score is calculated
	//	if(s_iterator <= sequence.size() - matrix[0].size()) {
		for(unsigned int i=0; i< matrix[0].size(); i++){

			switch(sequence[i+s_iterator]){

				case 'A':

					sum_scores += matrix[0][i];
					break;

				case 'C':

					sum_scores += matrix[1][i];
					break;

				case 'G':

					sum_scores += matrix[2][i];
					break;

				case 'T':

					sum_scores += matrix[3][i];
					break;

				default:				//Case if there is N

					sum_scores += o_matrix_mins[i];
					break;
			}
		}
		
		//The total score of an oligo is saved into an oligo_scores vector
		oligo_scores.emplace_back(sum_scores);

		//Then the function is recalled recursively shifting on the next oligo thanks to the iterator progression
		//shifting(matrix, sequence, s_iterator+1);

	}

}

//Best score normalization with normalization formula (The parameter to normalize have already been calculated and saved into the class)
void oligo_class::scores_normalization(){
	//PROFILE_FUNCTION();
	double score_normalized;
	vector<double> oligo_scores_normalized;

	for(unsigned int score=0; score < oligo_scores.size(); score++){

	score_normalized = 1 + ((oligo_scores[score] - max_possible_score)/(max_possible_score - min_possible_score));
	oligo_scores_normalized.emplace_back(score_normalized);
	}
	
	oligo_scores = oligo_scores_normalized;
}

//Function to find the best oligo score. From every sequence from GEP the best oligo is calculated and both oligo and position in the window are saved
unsigned int oligo_class::find_best_score(){
	//PROFILE_FUNCTION();
	//Extracting the best score from oligo_scores with function max_element
	best_score = *max_element(oligo_scores.begin(), oligo_scores.end());

	vector<int> positions;
	vector<unsigned int> dist_center;
	unsigned int matches = 0;
	unsigned int min_distance;
	vector<unsigned int>::iterator itr;
	
	//Check if there are more than one oligo with the best score --> if any count their numbers and save their position in sequence into a position vector
	for(unsigned int i=0; i < oligo_scores.size(); i++){

		if(oligo_scores[i] == best_score){

			matches = matches + 1;
			positions.emplace_back(i);
		}
	}

	//If more than one oligo carries the best score calculate distances and select the nearest to the center
	if(matches > 1){ 

		for (int& p: positions){

			unsigned int distance;

			//The distance from the center need to be calculated as an absolute value (no negative values)
			distance = abs( p - half_length); 
			dist_center.emplace_back(distance);
		}

		//Once distances have been calculated select as best the one who is the nearest to the sequence center --> select the min element from distance vector
		min_distance = *min_element(dist_center.begin(), dist_center.end());

		//Find min_distance on dist_centre vector and save its index
		itr = find(dist_center.begin(),dist_center.end(), min_distance);
		unsigned int index = distance(dist_center.begin(), itr);

		//Index of min_distance on distance vector == the corrisponding position in the positions vector
		return positions[index];
	

	}

	//If just one best oligo score has been found return its position (it is the first and the only element in positions vector)
	return positions[0];

}

//The oligo which has generated the best score is extracted from fasta sequence and saved into best_oligo_seq variable
void oligo_class::find_best_sequence(string sequence, unsigned int length){
	//PROFILE_FUNCTION();
	best_oligo_seq = sequence.substr(local_position,length);

}

//Best oligo coordinates are saved
void oligo_class::find_coordinate( unsigned int length, string chr_coord_GEP, unsigned int start_coord_GEP){
	//PROFILE_FUNCTION();
	chr_coord_oligo = chr_coord_GEP;
	start_coord_oligo = start_coord_GEP + local_position;
	end_coord_oligo = start_coord_oligo + length;
	
}
	
//Function to read BED and 2Bit files and create GEP (vector of bed class)
void coordinator_class::GEP_creation(vector<bed_class> &GEP){
	//PROFILE_FUNCTION();
	//RAM_usage();
	cout << "\n- [1] Extract bed coordinate sequences from reference genome  \n";

	ifstream in(BED_FILE); 					
	//TwoBit * tb;

	//Opening 2Bit file with twobit_open function from andrelmartens code and saved in tb variable 
	//tb = twobit_open(TWOBIT_FILE.c_str()); 

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
	if (err == true){
		exit(1);
	}
	
}

//Function to read JASPAR PWM file, extract values and create a matrix class
vector<vector<double>> coordinator_class::read_JASPAR(){
	//PROFILE_FUNCTION();
	//RAM_usage();
	cout << "- [2] Reading JASPAR MATRIX file and extracting values\n";

	ifstream file(JASPAR_FILE);
	string line;		

	//For each line of the JASPAR file	
	while(getline(file,line)){
		
		//If the line start with '>' character save the words into matrix_name string and into tf_name string
		if(line[0]=='>'){

			istringstream mystream(line);			
			//mystream >> matrix_name >> tf_name;
		}

		//If the line does not start with '>' delete the '[',']' and 'A,T,C,G' characters, extract the scores and save them into a scores matrix
		else{
			
			//Deleting from line the first character (A,T,C,G), the '[' and the ']' characters
			line.erase(0,line.find('[') +1);
			line.erase(line.find(']'));

			vector<double> scores_line;	
			istringstream mystream(line);
			
			//For each words (number) in line put the current number in num variables
			for (double num; mystream >> num;){

				scores_line.emplace_back(num);
			}

			matrix.emplace_back(scores_line);
		}
	}
    
	file.close();
	//RAM_usage();
	//Cout of step 3/4 here because the normalization and reverse function will be re-utilized during the workflow
	cout << "- [3] Jaspar Matrix normalization\n";
	cout << "- [4] Jaspar Matrix reverse complement determination to analize the reverse strand\n";

	return matrix;
}

//Function to create oligos_vector (a oligo class vector)
void coordinator_class::oligos_vector_creation(vector<oligo_class> &oligos_vector, vector<vector<double>> &matrix_log, vector<vector<double>> &matrix_log_inverse, vector<bed_class> &GEP){
	//PROFILE_FUNCTION();
	//RAM_usage();
	cout << "- [5] Analyzing sequences using Jaspar matrix\n";

	//For every sequences into GEP vector
	for(unsigned int i=0; i<GEP.size(); i++){

		//Calling the oligo_class constructor to analyze the shifting of the sequence on log_matrix (FWD strand analysis)
		oligo_class SHIFTING(matrix_log, GEP[i].sequence, GEP[i].chr_coord, GEP[i].start_coord, '+');

		//The oligo class just created is saved into oligos_vector (oligo_class vector)
		oligos_vector.emplace_back(SHIFTING);

		//If the analysis is on DS calling the oligo_class constructor to analyze the shifting of sequence on inverse_log_matrix (REVERSE strande analysis)
		if(DS == 1){

			oligo_class SHIFTING(matrix_log_inverse, GEP[i].sequence, GEP[i].chr_coord, GEP[i].start_coord, '-');
			oligos_vector.emplace_back(SHIFTING);
		}
	}	
	//RAM_usage();
	cout << "- [6] Selecting the best Jaspar's oligo for each sequence \n";
}

//Function useful, if the analysis is performed on DS, to choose from the best FWD strand oligo and the best REV strand oligo the best one to keep as "Best oligo" --> The oligos vector is divided in half and only the best strand for each sequence is kept
void coordinator_class::best_strand(){
	//PROFILE_FUNCTION();
	vector<oligo_class> comparison;
	
	for(unsigned int i=0; i<oligos_vector.size(); i+=2){
	
		//The comparison is made by oligo_class in i position against the oligo class in i+1 position (The fwd and rev strand of the same sequence, which are consecutive into the oligos_vector)
		double best_score_norm_positive = oligos_vector[i].best_score;
		double best_score_norm_negative = oligos_vector[i+1].best_score;

		if(best_score_norm_positive >= best_score_norm_negative){
			
			comparison.emplace_back(oligos_vector[i]);

		}

		else{
			comparison.emplace_back(oligos_vector[i+1]);
			
		}
	}

		//The new oligos_vector is replaced by comparison vector, which contains only the best strand
	oligos_vector.clear();
	
	//oligos_vector.shrink_to_fit();
	oligos_vector = comparison;
}

//Function to re-set the genomic coordinates and the sequences window --> centered on the best oligo found for each sequence
void coordinator_class::centering_oligo(){
	//PROFILE_FUNCTION();

	int center_oligo ;

	//To center on the best oligo, centering_function and extract_seq functions from bed_class need to be recalled with updated input parameters
	for(unsigned int i=0; i<oligos_vector.size(); i++){

		//The center of the window is exactly on the center of the best oligo (which length depends to the JASPAR matrix size)
		center_oligo = oligos_vector[i].start_coord_oligo + matrix_log[0].size()/2;
		GEP[i].centering_function(center_oligo,center_oligo,half_length,0);
		GEP[i].extract_seq(tb,0);
	}
}


//Function able to convert a string (containing numbers separated by ",") into a vector of unsigned int
vector<unsigned int> map_class::generic_vector_creation(string numbers){
	//PROFILE_FUNCTION();
	int index;
	vector<unsigned int> vec;

	//When index is == -1 means that it is pointin to the end character of the string
	while(index != -1){

		//Find the first "," character into the string
		index = numbers.find(",");

		//Put everything before "," into a vector
		vec.emplace_back(stoi(numbers.substr(0,index)));

		//Erase from the string the number already inserted and the ","
		numbers.erase(0,index+1);
	}

	return vec;
}

//Function to extract fasta sequences from a Multifasta file and save them into a vector of string
void multifasta_class::extract_sequences(){
	//PROFILE_FUNCTION();
	cout << "\n- [1] Extracting sequences from MultiFasta file \n";

	ifstream file(MFASTA_FILE);
	string line;
	string current_sequence;

	//First line is flagged to 1 to ignore it (because in MF format sometimes there is an header first line)
	bool first_line = 1;
	
	//The while cycle follows a complex reasoning: 
	//if in line there is the first line --> all ignored --> then first line flag set to 0
	//if in line there is the FASTA sequence --> line is saved in current_sequence variable --> then pass to next line
	//if in line there is the header (but not the first line) --> current_sequence previously filled by FASTA is saved in sequences vector (a vector of strings) --> then current_sequence is clean
	while(getline(file,line)){

		if(line[0] == '>' && !first_line){
			
			
			sequences.emplace_back(current_sequence);
			current_sequence.clear();
			
		}

		else if (!first_line){

			if(line[0] != ' ' && line.size() != 0){	

				//Before to save the sequence in current_sequence variable the FASTA characters are capitalized
				transform(line.begin(), line.end(), line.begin(), ::toupper);	
				current_sequence = current_sequence + line; 
			}
		}
		
		//After the first cicle the first line flag is set to 0
		first_line = 0;	
	}
	sequences.emplace_back(current_sequence);
}

//Controlling that all the MF sequences have the same langth
void multifasta_class::length_control(vector<string> sequences){
	//PROFILE_FUNCTION();
	cout << "- [2] Multifasta Sequences length check\n";

	//unsigned int size = sequences[0].size();

	for(unsigned int i=0; i<sequences.size(); i++){
		
		//If only one sequence in the vector are longer an error is generated
		if(sequences[i].size() != sequences[0].size()){

			cerr << "Sequences are not of the same length!\n";
			exit(1);
		}
	}
}

//Creating alias file following the multifasta input filename
void multifasta_class::alias_output_filename(){
	//PROFILE_FUNCTION();
	alias_file = alias_file + MFASTA_FILE.erase(0,MFASTA_FILE.find_last_of("_")+1); 
	alias_file = MFASTA_FILE.erase(MFASTA_FILE.find_last_of("."), MFASTA_FILE.size()) + "_";
}

void multifasta_class::GEP_creation_MF(vector<string> sequences){
	//PROFILE_FUNCTION();
	cout << "- [3] Sorting Multifasta sequences\n";

	for(unsigned int i=0; i<sequences.size(); i++){
		
		//A bed class is created for each sequence
		bed_class BED_MULTIFASTA(sequences[i]);

		//The classes are stored into a GEP vector
		GEP.emplace_back(BED_MULTIFASTA);
	}

}

//Checking function to control if, for any k-mers inserted as input, there is a distance parameter
void map_class::check_kmer_dist(){
	//PROFILE_FUNCTION();
	
	if(kmers_vector.size() != distance_vector.size()){
		
		//If the number is not equal --> ERROR printed and help visualized
		cerr << "\nERROR: Please insert an equal number of k-mers and distance parameters!\n";
		display_help();
		exit(1);
	}
	for(unsigned int i = 0; i < kmers_vector.size(); i++){
		if(distance_vector[i] > kmers_vector[i]){
			cerr << "\nERROR: Please distance parmater must be lower or equal to k-mers length!\n";
		display_help();
		exit(1);
		}
	}
}

//Function to create a map of oligos + their occurrences along all the sequences
void map_class::table_creation_orizzontal(vector<bed_class> &GEP){ 
	//PROFILE_FUNCTION();

	if (MFASTA_FILE.size() == 0 ){
		//RAM_usage();
		cout << "- [7] Counting all k-mers occurrences for sequences and making orizzontal maps  \n";
	}
	else{
		//RAM_usage();
		cout << "- [4] Counting all k-mers occurrences for sequences and making orizzontal maps  \n";
	}
	//A map is created for each k-mer inserted as input
	for(unsigned int k=0; k<kmers_vector.size(); k++){
		
		//For every sequence contained into GEP vector
		for(unsigned int j=0; j<GEP.size(); j++){
			
			//Extracted and analyzed all words of length k that are found by scrolling through the sequence
			for(unsigned int i=0; i < (GEP[j].sequence.size() - kmers_vector[k] + 1); i++){
				
				//The current k-length oligo is saved into bases string
				string bases = GEP[j].sequence.substr(i,kmers_vector[k]);
				//kmer_oligo.emplace_back(bases);
				
				//Function to fill orizzontal plus/minus maps --> current oligo "bases" is passed as parameter to be inserted into the maps
				or_ver_kmer_count(bases,orizzontal_plus,orizzontal_minus);
			}
		}

		//Once map is created it is saved into a vector of map --> then an analysis with new k value (if any) is performed
		orizzontal_plus_debug.emplace_back(orizzontal_plus);
		orizzontal_minus_debug.emplace_back(orizzontal_minus);

		//Once maps are saved they need to be cleaned to avoid any interference between two different k-mers analysis
		orizzontal_plus.clear();
		orizzontal_minus.clear();
	}
	//RAM_usage();
}

//Function to fill orizzontal plus/minus maps --> current oligo "bases" is passed as parameter to be inserted into the maps
void map_class::or_ver_kmer_count(string bases,unordered_map<string,unsigned int> &plus, unordered_map<string,unsigned int> &minus){
	//PROFILE_FUNCTION();
	unordered_map<string,unsigned int>::iterator it_plus;
	unordered_map<string,unsigned int>::iterator it_minus;

	//Finding the oligo "bases" into the orizzontal plus map (FWD strand)
	it_plus = plus.find(bases);

	//Check if the current oligo "bases" is palindrome --> this function has also the aim to generate "reverse bases", the reverse complement of the current oligo
	check_palindrome(bases, reverse_bases);

	//Finding the oligo "bases" into the orizzontal minus map (REV strand)
	it_minus = minus.find(reverse_bases);
	
	//If the current oligo "bases" is already present into the map --> increase the occurrences by 1
	if(it_plus!=plus.end()){

		it_plus->second++;
		it_minus->second++;
	}
	
	//If the oligo is not present yet --> insert it into the map with an occurrences value of 1
	else{

		plus.insert({bases,1});
		minus.insert({reverse_bases,1});
	}

	//Clear the string bases and reverse bases to avoid interference from different oligos
	bases.clear();
	//bases.shrink_to_fit();
	reverse_bases.clear();
	//reverse_bases.shrink_to_fit();
}

//Function to create maps to count oligos occurrences for each sequence position (in this function also frequences are calculated)
void map_class::table_creation_vertical(vector<bed_class> &GEP){
	//PROFILE_FUNCTION();
	if (MFASTA_FILE.size() ==0){
		//RAM_usage();
		cout << "- [8] Counting all k-mers positional occurrences and making vertical maps  \n";
	}
	else{
		//RAM_usage();
		cout << "- [5] Counting all k-mers positional occurrences and making vertical maps  \n";
	}

	//A vector of map is created for each k-mer inserted as input
	for(unsigned int k=0; k<kmers_vector.size(); k++){

		//vector to store the total number of possible oligo per position, useful to calculate frequences
		vector<unsigned int> tot_freq_vec;


		//Extracted and analyzed the oligo in position "i"
		for(unsigned int i=0; i < (GEP[0].sequence.size() - kmers_vector[k] + 1); i++){

			unsigned int tot_freq = 0;
			
			//Make the analysis of all the sequences' oligo in position "i" (vertical analisys)
			for(unsigned int j=0; j<GEP.size(); j++){
				
				//string sequence = GEP[j].sequence;
				string bases = GEP[j].sequence.substr(i,kmers_vector[k]);
				
				//Calling of function to count oligo occurrences and to create and fill the maps
				vertical_kmer_count(bases, vertical_plus, tot_freq);
			}
			
			if(DS==1){

				//Into vertical plus map we have (oligo + RC) pair but also (RC + oligo) pair --> here we select the best one to keep (it will be the one who has the most occurrences in FWD strand)
				select_best(vertical_plus);
			}

			//The maps just created is saved into a vector of maps (one for each position)
			maps_vector_positions_plus.emplace_back(vertical_plus);

			//Then the maps are cleaved to make a new analysis on a new position
			vertical_plus.clear();
			
			//the counting of all the possible oligos is saved into a vector too
			tot_freq_vec.emplace_back(tot_freq);
		}
		
		//The vector of maps created is stored into a vector to create a matrix in which on the dimension 1 we have the k-mers used and on dimension 2 the positions on sequences
		vector_kmers_maps_plus.emplace_back(maps_vector_positions_plus);

		//Also the all possible oligos to calculate frequences are stored in a matrix (k-mers x position)
		tot_freq_matrix.emplace_back(tot_freq_vec);
		
		//The maps vector are cleaved to make a new analysis with new k
		maps_vector_positions_plus.clear();

		//RAM_usage();
	}
}

//Function to count oligo occurrences and to create and fill the maps. It also count all the total oligo present in position (taking into account to the palindrome oligos) --> this count will be useful to calculate the frequency
void map_class::vertical_kmer_count(string bases,map<pair<string,string>,pair<unsigned int, unsigned int>>&plus, unsigned int& tot_freq){
	//PROFILE_FUNCTION();
	map<pair<string,string>,pair<unsigned int, unsigned int>>::iterator it_plus;
	map<pair<string,string>,pair<unsigned int, unsigned int>>::iterator it_plus_rev;

	//Check if the current oligo is palindrome
	bool pal = check_palindrome(bases, reverse_bases);

	//Check if the analisys is in double strand
	if(DS == 1){	

		//If the current oligo is not palindrome and the analysis is on DS --> total number of possible oligo per position (to calculate frequences) ++2, else if it is palindrome ++1
		if(!pal){
			tot_freq += 2;
		}
		else{
			tot_freq++;
		}
	}

	//If analysis is on SS --> possible oligos per position ++1 (it will be equal to the number of sequences)
	else{
		tot_freq++;
	}
	/*
	//Creation of a pair of strings containing current oligo (first) and his reverse complement (second)
	pair<string,string> pair_bases;
	pair_bases.first= bases;
	pair_bases.second= reverse_bases;
	
	//Creation of a pair of strings containing reverse complement (first) and current oligo (second)
	pair<string,string> pair_bases_reverse;
	pair_bases_reverse.first= reverse_bases;
	pair_bases_reverse.second= bases;
	*/
	//Try to find the two pairs just created into plus map
	it_plus = plus.find(make_pair(bases, reverse_bases));
	it_plus_rev = plus.find(make_pair(reverse_bases, bases));

	//If the pair oligo + RC is already present in the plus map --> increase occurrences by 1
	if(it_plus!=plus.end()){

		it_plus->second.first++;
	}	
	
	//Else If the pair RC + oligo is already present in the plus map --> increase occurrences by 1
	else if (it_plus==plus.end() && it_plus_rev != plus.end()) {
	
		it_plus_rev->second.second++;

	}
	
	//Else insert the pair as new into plus and minus map
	else{
		
		plus.insert({{bases,reverse_bases},{1,0}});
	}

	//Clear bases and reverse bases to avoid interference in the next oligo analysis 
	bases.clear();
	reverse_bases.clear();
	
}

//Function to select the best pair (Oligo + RC or RC + Oligo) into vertical plus map
void map_class::select_best(map<pair<string,string>,pair<unsigned int,unsigned int>>& vertical_plus){
	//PROFILE_FUNCTION();
	//Initialization of an empty map called copy
	map<pair<string,string>,pair<unsigned int,unsigned int>> copy;
	
	//For every element into the map a comparison is performed
	for(map<pair<string,string>,pair<unsigned int,unsigned int>>::iterator it = vertical_plus.begin(); it!=vertical_plus.end(); it++){

		//If Oligo occurrences are less then RC occurrences in FWD strand
		if(it->second.first < it->second.second){
			
			string oligo1 = it->first.second;	
			string oligo2 = it->first.first;
			unsigned int occ1 = it->second.second;
			unsigned int occ2 = it->second.first;
			
			//Insert the pair RC + Oligo in the map "copy"
			//swap(it->second.first, it->second.second);
			//swap(oligo1, oligo2);
			copy.insert({{oligo1,oligo2},{occ1,occ2}});		
		}

		//Else if Oligo occurrences are higher then RC occurrences is FWD strand
		else{

			//Insert the pair Oligo + RC in the map "copy"
			copy.insert({{it->first.first, it->first.second},{it->second.first,it->second.second}});		
		}
	}

	//Clearing the map vertical plus to be substituted by "copy", which contains only the best pairs
	vertical_plus.clear();
	vertical_plus = copy;
}

//Function to create a matrix of p_value class (1 dimension = different k-mers, 2 dimension = positions)
void map_class::P_VALUE_MATRIX_creation(){
	//PROFILE_FUNCTION();
	if (MFASTA_FILE.size() ==0){
		//RAM_usage();
		cout << "- [9] Calculating oligos p_value and flling P_Value class matrix  \n";
	}
	else{
		//RAM_usage();
		cout << "- [6] Calculating oligos p_value and flling P_Value class matrix  \n";
	}
	//For each k-mers inserted as input
	for(unsigned int j=0; j<vector_kmers_maps_plus.size(); j++){
		
		//Outfile_header function to handle Output file
		ofstream outfile = outfile_header(j);
		
		//For each position in sequence
		for(unsigned int i=0; i<vector_kmers_maps_plus[j].size(); i++){
			
			//Create a p_value_class object P passing the vertical and orizzontal maps, the number of sequences T, i to fix the current position and the outfile to print into the Output file
			p_value_class P(vector_kmers_maps_plus[j][i], orizzontal_plus_debug[j], sequences_number_T, i, outfile, tot_freq_matrix[j][i]);

			//A vector of p_value classes is created 
			P_VALUE_VECTOR.emplace_back(P); 
			
			//The sum of the first N oligos occurrences are saved into tot_sum_vector
			tot_sum_vector.emplace_back(P.sum_top_N);
			
		}
		
		//The vector of p_value_classes is stored into a bigger P_VALUE_MATRIX
		P_VALUE_MATRIX.emplace_back(P_VALUE_VECTOR);

		//Vector is cleaned up to avoid interference on the next cycle
		P_VALUE_VECTOR.clear();

		//Sums vector of the first N oligo occurrences saved into a matrix and then is cleaned up
		tot_sum_matrix.emplace_back(tot_sum_vector);
		tot_sum_vector.clear();
		
		outfile.close();

	}
}

//Function to create a matrix of hamming class (1 dimension = different k-mers, 2 dimension = positions)
void map_class::HAMMING_MATRIX_creation(vector<bed_class> &GEP){
	//PROFILE_FUNCTION();
	if (MFASTA_FILE.size() ==0){
		//RAM_usage();
		cout << "- [10] Calculating best oligos hamming neighbours and filling Hamming matrix  \n";
	}
	else{
		//RAM_usage();
		cout << "- [7] Calculating best oligos hamming neighbours and filling Hamming matrix  \n";
	}
	//For every kmer in input	
	for(unsigned int j=0; j<P_VALUE_MATRIX.size(); j++){
		
		//Outfile_header function to handle Output file
		ofstream outfile = outfile_header_hamming(j);
		
		//For each position in sequence
		for(unsigned int i=0; i<P_VALUE_MATRIX[j].size(); i++){
			
			//Create a hamming_class H passing the vertical multimap, the distance inserted as input, the current position i, the number of different oligos (contained in tot_freq_matrix), the orizzontal matrix, the outfile to print the Output file and finally the GEP (which contains the sequences) --> the ordering p or not distinguishes the multimap passed as input
			
			hamming_class H(P_VALUE_MATRIX[j][i].vertical_multimap, vector_kmers_maps_plus[j][i], P_VALUE_MATRIX[j][i].p_value_sort,distance_vector[j],i,tot_freq_matrix[j][i],orizzontal_plus_debug[j], orizzontal_minus_debug[j], outfile, GEP, kmers_vector);
			
			//A vector of hamming classes is created 
			HAMMING_VECTOR.emplace_back(H);


		}

		//The vector of hamming_classes is stored into a bigger P_VALUE_MATRIX
		HAMMING_MATRIX.emplace_back(HAMMING_VECTOR);
		
		//Vector is cleaned up to avoid interference on the next cycle
		HAMMING_VECTOR.clear();
		//HAMMING_VECTOR.shrink_to_fit();

		outfile.close();
	}

	if (MFASTA_FILE.size() ==0){
		//RAM_usage();
		cout << "- [11] PWM matrices from hit positions calculated \n";
	}
	else{
		//RAM_usage();
		cout << "- [8] PWM matrices from hit positions calculated \n";
	}
}

//Function to create a matrix of z_test class
void map_class::Z_TEST_MATRIX_creation(vector<bed_class> &GEP){
	//PROFILE_FUNCTION();
	bool local_max = 1;
	
	if (MFASTA_FILE.size() ==0){
		//RAM_usage();
		cout << "- [12] Calculating Z-test from PWM matrices shifing and filling Z-test matrix  \n";
	}
	else{
		//RAM_usage();
		cout << "- [9] Calculating Z-test from PWM matrices shifing and filling Z-test matrix  \n";
	}
	//For every kmer in input	
	for(unsigned int i=0; i<HAMMING_MATRIX.size(); i++){
		
		//For each position in sequence
		for (unsigned int j=0; j<HAMMING_MATRIX[i].size(); j++){
			
			//Return the PWM matrix calculated from the corresponding hamming class
			//vector<vector<double>> PWM_matrix = HAMMING_MATRIX[i][j].PWM_hamming;

			//Extract the frequence_1 from the corresponding hamming class
			//double Frequence_1 = HAMMING_MATRIX[i][j].FREQUENCE_1;
			
			//If the local_maxima filtering is enabled
			if(local_maxima_grouping == true){
				
				
				double Frequence_1_prev = 0;
				double Frequence_1_post = 0;
				
				//if the position is not the first -> return the pos-1 freq_1 -> else the pos-1 freq_1 value remain 0
				if(j != 0){
					Frequence_1_prev = HAMMING_MATRIX[i][j-1].FREQUENCE_1;
				}
				
				//if the position is not the last -> return the pos+1 freq_1 -> else the pos+1 freq_1 value remain 0
				if(j != HAMMING_MATRIX[i].size()){
					Frequence_1_post = HAMMING_MATRIX[i][j+1].FREQUENCE_1;
				}
				
				//Analyze if the current pos is a local max
				local_max = find_local_max(HAMMING_MATRIX[i][j].FREQUENCE_1,Frequence_1_prev,Frequence_1_post);
			}
			
			//If it is a local max and its freq_1 value overcome the threshold build a z_test_class Z and then save it into a vector and finally into a matrix (as done with hamming and p_value classes)
			if(HAMMING_MATRIX[i][j].FREQUENCE_1 >= freq_treshold && local_max == 1){

				z_test_class Z(HAMMING_MATRIX[i][j].PWM_hamming, GEP,j+1,kmers_vector,HAMMING_MATRIX);
				Z_TEST_VECTOR.emplace_back(Z);
			}
		}
		Z_TEST_MATRIX.emplace_back(Z_TEST_VECTOR);
		Z_TEST_VECTOR.clear();
		
	}
}

//Return 1 only if the frequence in position is higher than frequences in pos-1 and in pos+1
bool map_class::find_local_max(double center, double prev, double post){
	//PROFILE_FUNCTION();
	
	if(center > prev && center >= post){

		return 1;
	}
	
	else{
		return 0;
	}

}

//Transforming a map of pair into a multimap, ordered by oligo occurrences
void p_value_class::multimap_creation(map<pair<string,string>, pair<unsigned int,unsigned int>> pair_map){
	//PROFILE_FUNCTION();
	for(map<pair<string,string>,pair<unsigned int, unsigned int>>::iterator it = pair_map.begin(); it != pair_map.end(); it++){ 	

		vertical_multimap.insert({it->second,it->first});
	}
}

//Finding N2 parameters thanks to the accumulate function, which allows to sum all the values present in a map -> accumulate function uses a lambda function (https://en.cppreference.com/w/cpp/algorithm/accumulate) 
void p_value_class::N2_calculation(unordered_map<string,unsigned int> &orizzontal_map){
	//PROFILE_FUNCTION();
	total_oligo_N2 = 0;
	total_oligo_N2 = accumulate(begin(orizzontal_map), end(orizzontal_map), 0, [] (unsigned int val, const unordered_map<string,int>::value_type& p) {return val + p.second;});
}

//Function to calculate and store into vectors K,N1,N2 parameters used to find oligos' p-values
void p_value_class::filling_KNT_vectors(unordered_map<string,unsigned int> &orizzontal_map){
	//PROFILE_FUNCTION();
	unsigned int K;
	unsigned int N1;
	unsigned int N2;

	//For each oligo present in positional vertical multimap (in the current position)
	for(multimap<pair<unsigned int, unsigned int>,pair<string,string>>::reverse_iterator it_rev = vertical_multimap.rbegin(); it_rev != vertical_multimap.rend(); it_rev++){

		//Check if the oligo is palindrome -> K value is differentially calculated
		bool pal = check_palindrome(it_rev->second.first, reverse_bases);

		if(pal == 0 && DS == 1){

			K = it_rev->first.first + it_rev->first.second;
		}
		else{

			K = it_rev->first.first;
		}
		
		//searching the current oligo in orizzontal map to find N1 value (oligo occurrences along all the sequences positions)
		it_N1_plus = orizzontal_map.find(it_rev->second.first);
		it_N1_minus = orizzontal_map.find(it_rev->second.second);
		
		//If the oligo RC is present in the orizzontal map (and analysis is in DS) the N1 parameter is calculated as the sum of the oligo + oligo RC occurrences in orizzontal map. Else, only the oligo occurrences can be used for the N1 definition.
		if(it_N1_minus != orizzontal_map.end() && DS == 1){

			N1 = it_N1_plus->second + it_N1_minus->second;
		}
		else{
			N1 = it_N1_plus->second;
		}
		
		//Using the total_oligo_N2 (calculated before in N2_calculation function) defining the N2 value as definition (Total number of k-mers - N1)
		N2 = total_oligo_N2 - N1;
		
		//Storing values into vectors
		K_vec.emplace_back(K);
		N1_vec.emplace_back(N1);
		N2_vec.emplace_back(N2);

	}
}

//Using gsl library hypergeometric_Q function and the parameters found before (K,N1,N2,T) calculate each oligo p_value and store it into a vector
void p_value_class::calculating_p_value(){
	//PROFILE_FUNCTION();
	for(unsigned int i=0; i<K_vec.size(); i++){

		double p_value =  gsl_cdf_hypergeometric_Q(K_vec[i],N1_vec[i],N2_vec[i],T);

		//Checking if p_value is too low and therefore is automatically rounded to 0
		p_value = check_p_value(p_value);
		p_value_vec.emplace_back(p_value);	
	}
}

//Function to sort all the oligos in each position following the lowest p_value scores
void p_value_class::sorting_p_value(){
	//PROFILE_FUNCTION();
	vector<unsigned int> KNT;
	unsigned int i=0;

	//Creating 2 multimaps:
	//1. p_value_sort = a multimap containing p_value and the pair oligo + oligo's RC to order oligos by lowest p_value 
	//2. p_value_KNT = a multimap composed by p_value at first position to guarantee the order and oligo + K,N1,N2,T vector at second position
	for(multimap<pair<unsigned int, unsigned int>, pair<string, string>>::reverse_iterator it_rev = vertical_multimap.rbegin(); it_rev!=vertical_multimap.rend(); it_rev++, i++){
		
		//Fill the KNT vector with the right K,N1,N2,T values
		KNT.emplace_back(K_vec[i]);	
		KNT.emplace_back(N1_vec[i]);	
		KNT.emplace_back(N2_vec[i]);	
		KNT.emplace_back(T);	

		p_value_sort.insert({p_value_vec[i],{it_rev->second.first, it_rev->second.second}});
		p_value_KNT.insert({p_value_vec[i],{it_rev->second.first, KNT}}); 
		
		//Clear the KNT vector to mantain its size (need to be formed by 4 elements(K,N1,N2,T))
		KNT.clear();
	}
}

//Check and select if oligos are going to be printed in output file ordered by occurrences or p_values
void p_value_class::checking_ordering(map<pair<string,string>,pair<unsigned int,unsigned int>> &pair_map, unsigned int position, ofstream &outfile, unsigned int freq){
	//PROFILE_FUNCTION();
	//If p_value ordering has been selected
	if(ordering == "p"){
		
		//Select for Double strand analysis or Single strand analysis
		if(DS==1){
			
			print_debug_p_value_DS(pair_map, position, outfile, freq);	
		}	

		else{

			print_debug_p_value_SS(pair_map, position, outfile, freq);	
		}
	}

	else{
		if(DS == 1){

			print_debug_occurrences_DS(pair_map, position, outfile, freq, p_value_vec);	
		}

		else{

			print_debug_occurrences_SS(pair_map, position, outfile, freq, p_value_vec);	
		}
	}
}


multimap<unsigned int,pair<string,string>> hamming_class::creating_sum_occurrences_multimap(multimap<pair<unsigned int,unsigned int>,pair<string,string>>& pair_map){

	multimap<unsigned int,pair<string,string>> sum_occurrences_multimap;

	for(multimap<pair<unsigned int,unsigned int>,pair<string,string>>::reverse_iterator rev_it = pair_map.rbegin(); rev_it != pair_map.rend(); rev_it ++){

		sum_occurrences_multimap.insert({(rev_it->first.first + rev_it->first.second),make_pair(rev_it->second.first,rev_it->second.second)});
	}

	return sum_occurrences_multimap;
}

//Scrolling the vertical positional multimap find the best oligo for occurrences. If more than one is present selecting the oligo which has more hamming neighbours.
void hamming_class::find_best_oligos(multimap<pair<unsigned int,unsigned int>, pair<string,string>>& vertical_multimap, map<pair<string,string>,pair<unsigned int,unsigned int>>& vertical_map, multimap<unsigned int,pair<string,string>>& sum_occurrences_multimap){


	//PROFILE_FUNCTION();
	multimap<unsigned int, pair<string,string>>::reverse_iterator it_rev_DS = sum_occurrences_multimap.rbegin();
	multimap<pair<unsigned int,unsigned int>, pair<string,string>>::reverse_iterator it_rev_SS = vertical_multimap.rbegin();

	//Saving the best oligo occurrences extracting from the last multimap element (higher first value) its occurrences -> working with multimap allows to have always the max occurrences in the last position of the map
	if(DS == 1){
		real_best_oligo_occurrences = it_rev_DS->first;
	
	}

	else{
		
		real_best_oligo_occurrences = it_rev_SS->first.first;
	}

	//Flag and counter to control if the function does more cycle than multimap size (to avoid infinite loop as happened)
	bool flag = 1;
	unsigned int counter = 1;

	//If all the sequences have the same oligo (100%) in a specific position --> vertical_multimap.size() == 1 --> This control is made to avoid an infinite while cycle
	//This control is specific for DS analysis
	if(DS==1){
		while(it_rev_DS->first == real_best_oligo_occurrences && flag == 1){

			if(counter == vertical_multimap.size()){
				flag = 0;
			}
			
			//If another oligo has the same occurrences as the best -> save it into a vector of best_oligos
			best_oligos.emplace_back(it_rev_DS->second.first);
			it_rev_DS++;
			counter++;
		}
	}

	//If all the sequences have the same oligo (100%) in a specific position --> vertical_multimap.size() == 1 --> This control is made to avoid an infinite while cycle
	//This control is specific for SS analysis
	else{
		while(it_rev_SS->first.first == real_best_oligo_occurrences && flag == 1){

			if(counter == vertical_multimap.size()){
				flag = 0;
			}

			//If another oligo has the same occurrences as the best -> save it into a vector of best_oligos
			best_oligos.emplace_back(it_rev_SS->second.first);
			it_rev_SS++;
			counter++;
		}
	}
}

//Checking if there is one or more best oligos
void hamming_class::checking_best_oligo(unsigned int distance, multimap<pair<unsigned int,unsigned int>, pair<string,string>> &vertical_multimap){
	//PROFILE_FUNCTION();
	//If there is only one best oligo for occurrences
	if(best_oligos.size() == 1){
		
		//Set that the real best oligo is it
		real_best_oligo = best_oligos[0];

		//Proceed to find his hamming distance neighbours
		find_distanced_oligos(real_best_oligo,distance, vertical_multimap);	
	}
	
	//else means that there are more than one best oligo for occurrences
	else{
		
		//Set the real best oligo after a selection function
		real_best_oligo = select_real_best_oligo(distance, vertical_multimap);

		//Proceed to find his hamming distance neighbours
		find_distanced_oligos(real_best_oligo,distance, vertical_multimap);	
	}

}

//Function for find the real_best_oligo (selecting the one who has more neighbours than the others) if more than one best oligo has been found
string hamming_class::select_real_best_oligo(unsigned int distance, multimap<pair<unsigned int,unsigned int>, pair<string,string>> &vertical_multimap){
	//PROFILE_FUNCTION();
	unsigned int max_similarity;
	unsigned int index;
	
	//For each best oligo found
	for(unsigned int i=0; i<best_oligos.size(); i++){
		
		//Find all of its hamming neighbours
		find_distanced_oligos(best_oligos[i], distance, vertical_multimap);
		
		//If the analysis is on the first best oligo -> the max similarity (that is the max number of neighbour) corresponds to its similar_oligos number
		if(i==0){

			max_similarity = similar_oligos.size();
			index = 0;
		}	
		
		//If an oligo has a number of neighbour (similar_oligo.size()) which exceeds the best partial so far (max_similarity) -> its neighbour number becames the best and the index of oligo (into the best_oligos vector) is saved.
		if(similar_oligos.size() > max_similarity){

			max_similarity = similar_oligos.size();
			index = i;
		}
		
		//Similar oligos vector need to be cleaned up for the next cycle
		similar_oligos.clear(); 
		//similar_oligos.shrink_to_fit();
		//And also similar_oligos_occurrences vector needs to be cleaned up to avoid subsequently interferences (because the function find_distanced_oligos calculates also the neighbours occurrences and saves them into similar_oligo_occurrences vector)
		similar_oligos_occurrences.clear();  
		//similar_oligos_occurrences.shrink_to_fit();
	}
	
	//Return the real best oligo
	return best_oligos[index];
}

//Function to perform a secondary hamming --> find of distanced d hamming from the similar oligos found
void hamming_class::find_secondary_hamming(unsigned int distance, unsigned int number_first_hamming,multimap<pair<unsigned int,unsigned int>, pair<string,string>> &vertical_multimap){
	//PROFILE_FUNCTION();
	//For each similar oligo
	for(unsigned int neighbour=0; neighbour < number_first_hamming; neighbour++){

		find_distanced_oligos(similar_oligos[neighbour],distance, vertical_multimap);
	}
}

//Function to find oligo's hamming neighbours
void hamming_class::find_distanced_oligos(string best, unsigned int distance, multimap<pair<unsigned int,unsigned int>, pair<string,string>> &vertical_multimap){
	//PROFILE_FUNCTION();
	bool is_similar;

	//Scrolling all the positional vertical multimap
	for(multimap<pair<unsigned int,unsigned int>, pair<string,string>>::reverse_iterator it_rev = vertical_multimap.rbegin(); it_rev != vertical_multimap.rend(); it_rev++){

		//If the oligo in analysis is not the real_best_oligo (to avoid the comparison to himself and to avoid to be added to similar oligos vector in case of secondary hamming)
		if(it_rev->second.first != best && it_rev->second.first != real_best_oligo){
			
			//Call the function is_similar_oligo to compare the real_best_oligo to the current oligo. The function returns 1 if it is, otherwise 0
			is_similar = is_similar_oligo(best, it_rev->second.first, distance);

			//If real_best_oligo and the current oligo are similar they can be considered neighbours and the current oligo (from FWD strand) and its occurrences can be saved into vectors
			if(is_similar == 1){

				//If -r option is active check if the neighbour found is already present in similar oligos vector
				if(refining_matrix == 1){

					bool is_present = checking_neighbour_presence(it_rev->second.first);

					//If the oligo is not present --> add it
					if(is_present == 0){

						similar_oligos.emplace_back(it_rev->second.first);
						similar_oligos_occurrences.emplace_back(it_rev->first.first);
					}
				}

				else{

					similar_oligos.emplace_back(it_rev->second.first);
					similar_oligos_occurrences.emplace_back(it_rev->first.first);

				}
			}

			//If analysis is in DS -> call the function is_similar_oligo to compare the real_best_oligo to the current oligo's Reverse Complement (RC). The function returns 1 if it is, otherwise 0
			if(DS==1){

				is_similar = is_similar_oligo(best, it_rev->second.second, distance);
				//If they are similar add the oligo's RC and its occurrences to neighbours vectors 
				if(is_similar == 1){

					//If -r option is active check if the neighbour found is already present in similar oligos vector
					if(refining_matrix == 1){

						bool is_present = checking_neighbour_presence(it_rev->second.second);

						//If the oligo is not present --> add it
						if(is_present == 0){

							similar_oligos.emplace_back(it_rev->second.second);
							similar_oligos_occurrences.emplace_back(it_rev->first.second);
						}
					}
					else{

						similar_oligos.emplace_back(it_rev->second.second);
						similar_oligos_occurrences.emplace_back(it_rev->first.second);
					}
				}
			}
		}	
	}
}

//Function to check if an oligo is already present in similar_oligos vector
bool hamming_class::checking_neighbour_presence(string oligo){
	//PROFILE_FUNCTION();
	for(unsigned int neigh = 0; neigh < similar_oligos.size(); neigh++){
		
		if(oligo == similar_oligos[neigh]){

			return 1;
		}
	}
	
	return 0;
}

//Function to compare two oligos and find out how many charachters are apart (distanced)
bool hamming_class::is_similar_oligo(string oligo_1, string oligo_2, unsigned int distance){
	//PROFILE_FUNCTION();
	//Counter to count the character differences 
	unsigned int counter = 0;
	
	//For each character of the oligos (which have the same dimension)
	for(unsigned int i = 0; i<oligo_1.size() && counter <= distance; i++){
		
		//If the letter is different -> increase the counter
		if(oligo_1[i] != oligo_2[i]){
			
			counter++;	
		}
	}
	
	//Return 1 if counter is lower or equal to distance, otherwise 0
	return(counter<=distance);
}

//Function to calculate the frequence_1 (total of similar occurrences / total of possible oligos in the position)
double hamming_class::frequence_1_calculation(unsigned int freq){
	//PROFILE_FUNCTION();
	tot_similar_occurrences = 0;
	//Sum of the total occurrences (real best oligo + his neighbours)	
	for(unsigned int i=0; i<similar_oligos_occurrences.size(); i++){

		tot_similar_occurrences += similar_oligos_occurrences[i];

	}

	FREQUENCE_1 = tot_similar_occurrences/freq;
	return FREQUENCE_1;

}

//Function to calculate the frequence_2 (total of similar occurrences / total number of best+hamming occurrences in sequences)
double hamming_class::frequence_2_calculation(unordered_map<string,unsigned int> &orizzontal_map_plus, unordered_map<string,unsigned int> &orizzontal_map_minus, unsigned int position){
	//PROFILE_FUNCTION();

	unsigned int total_orizzontal_occurrences = finding_orizzontal_occurrences(orizzontal_map_plus, orizzontal_map_minus);
	FREQUENCE_2 = tot_similar_occurrences/total_orizzontal_occurrences;
	return FREQUENCE_2;
}

//Function to find total number of best+hamming occurrences in sequences --> useful to calculate Freq_2
unsigned int hamming_class::finding_orizzontal_occurrences(unordered_map<string,unsigned int> &orizzontal_map_plus, unordered_map<string,unsigned int> &orizzontal_map_minus){
	//PROFILE_FUNCTION();
	//Search the best oligo in orizzontal map plus 
	unordered_map<string,unsigned int>::iterator it;
	unsigned int total_orizz_occ = 0;
	
	//search orizzontal occurrences for all the hamming neighbours	
	for(unsigned int i=0; i<similar_oligos.size(); i++){
		
		it = orizzontal_map_plus.find(similar_oligos[i]);
		if(it == orizzontal_map_plus.end()){
			
			it = orizzontal_map_minus.find(similar_oligos[i]);
		}
		total_orizz_occ = (total_orizz_occ + it->second);
	}

	return total_orizz_occ;
}

//Starting from the best oligo and his hamming neighbours the function builds a PWM_matrices following their sequences and occurrences
void hamming_class::PWM_hamming_creation(){
	//PROFILE_FUNCTION();
	//Vector for each bases initialized to count, position by position, the occurrences of that base
	vector<double> vec_A, vec_C, vec_G, vec_T;
	
	//For each bases (position) of oligo
	for(unsigned int character = 0; character < similar_oligos[0].size(); character++){

		double counter_A = 0;
		double counter_C = 0;
		double counter_G = 0;
		double counter_T = 0;
		
		//For each oligo in similar oligos vector
		for(unsigned int oligo = 0; oligo < similar_oligos.size(); oligo++){
			string oligo_map = similar_oligos[oligo];
			similar_oligos_map.insert(pair<string, double>(oligo_map, similar_oligos_occurrences[oligo]));
			switch(similar_oligos[oligo][character]){

				//Increment base counters of the oligo occurrences
				case 'A' : counter_A += similar_oligos_occurrences[oligo]; 
					   break;
				case 'C' : counter_C += similar_oligos_occurrences[oligo]; 
					   break;
				case 'G' : counter_G += similar_oligos_occurrences[oligo];
					   break;
				case 'T' :  counter_T += similar_oligos_occurrences[oligo]; 
					   break;
			}
		}
		
		//Fill base vectors
		vec_A.emplace_back(counter_A);
		vec_C.emplace_back(counter_C);
		vec_G.emplace_back(counter_G);
		vec_T.emplace_back(counter_T);
	}
	
	//Build the PWM matrix
	PWM_hamming.emplace_back(vec_A);
	PWM_hamming.emplace_back(vec_C);
	PWM_hamming.emplace_back(vec_G);
	PWM_hamming.emplace_back(vec_T);
}
void hamming_class::print_PWM(string name, vector<vector<double>> PWM_hamming){
	cout << name << endl;
	for (unsigned short int i = 0; i<PWM_hamming.size(); i++){
		for (unsigned short int j = 0; j<PWM_hamming[i].size(); j++){
			cout << PWM_hamming[i][j] << "\t";
		}
		cout << endl;
	}
	cout << endl;
}
//In this function we transfrom the position frequency matrix into position probability matrix
void hamming_class::EM_Ipwm(vector<vector<double>> &PWM_hamming,vector<bed_class> &GEP) 
{
	//PROFILE_FUNCTION();

	//matrix::normalize()

	double sum = 0;
	double corr = 0;

	for(unsigned int j = 0; j < PWM_hamming.size(); j++){
		sum = sum + PWM_hamming[j][0];
	}

	corr = sqrt(sum);

	for(unsigned int x = 0; x < PWM_hamming.size(); x++){
		for(unsigned int y = 0; y < PWM_hamming[0].size(); y++){
			PWM_hamming[x][y] = PWM_hamming[x][y] + corr;
		}
	}

	sum = 0;
	for(unsigned int j = 0; j < PWM_hamming.size(); j++){
		sum += PWM_hamming[j][0];
	}
	for (unsigned short int i = 0; i<PWM_hamming.size(); i++){
		for (unsigned short int j = 0; j<PWM_hamming[i].size(); j++){
			PWM_hamming[i][j] = PWM_hamming[i][j]/sum;
		}
	}	
}

/*
This function is about the expectation step of the expectation-maximization algorithm
where we obtain the likelihood ratio for each oligo in the vertical map
*/

void hamming_class::EM_Epart(vector<bed_class> &GEP, unordered_map<string,unsigned int>& orizzontal_map_minus, unsigned int position,unordered_map<string,unsigned int>& orizzontal_map_plus) 
{
	//PROFILE_FUNCTION();
	/*
	cout << "---------------------------------" << endl;
	cout << "#POSITION " << position + 1 << " before Epart"<< endl;
	cout << "---------------------------------" << endl;

	print_PWM("", PWM_hamming);
	*/
	
	unsigned int sum_hor = 0;

	for(map<string, double>::iterator it = similar_oligos_map.begin(); it != similar_oligos_map.end(); it++){
		unordered_map<string,unsigned int>::iterator occ_oligo_it = orizzontal_map_plus.begin();
		unordered_map<string,unsigned int>::iterator occ_oligo_it_rev = orizzontal_map_minus.begin();
	
		occ_oligo_it = orizzontal_map_plus.find(it->first);
		occ_oligo_it_rev = orizzontal_map_minus.find(it->first);
		
		if(occ_oligo_it != orizzontal_map_plus.end()){
			sum_hor = sum_hor + occ_oligo_it->second;
		}
		if(occ_oligo_it_rev != orizzontal_map_minus.end()){
			sum_hor = sum_hor + occ_oligo_it_rev->second;
		}
	}

	vector<double> LR;
	vector<string> oligo;


	//In this cycle for each element in vertical map we calculate the probability that this oligo is present in the hamming matrix 
	for(map<string,double>::iterator it = similar_oligos_map.begin(); it != similar_oligos_map.end(); it++){
		unordered_map<string,unsigned int>::iterator occ_oligo_it = orizzontal_map_plus.begin();
		unordered_map<string,unsigned int>::iterator occ_oligo_it_rev = orizzontal_map_minus.begin();
		string similar_oligo = it->first;
		double P_oligo = 1;
		double P_bg = 0;
		double horizontal_occurences = 0;
		double likelihood_ratio = 0;

		occ_oligo_it = orizzontal_map_plus.find(similar_oligo);
		occ_oligo_it_rev = orizzontal_map_minus.find(similar_oligo);

		if(occ_oligo_it != orizzontal_map_plus.end()){
			horizontal_occurences = horizontal_occurences + occ_oligo_it->second;
		}
		if(occ_oligo_it_rev != orizzontal_map_minus.end()){
			horizontal_occurences = horizontal_occurences + occ_oligo_it_rev->second;
		}
		
		P_bg = horizontal_occurences/sum_hor;

		for (unsigned int k = 0; k < PWM_hamming[0].size(); k++){
			switch(it->first[k]){

				case 'A':			
					P_oligo *= PWM_hamming[0][k];
					break;

				case 'C':
					P_oligo *= PWM_hamming[1][k];
					break;

				case 'G':
					P_oligo *= PWM_hamming[2][k];
					break;

				case 'T':
					P_oligo *= PWM_hamming[3][k];
					break;

				default:				//Case if there is N
					P_oligo *=1;
					break;
			}
		}
		
		likelihood_ratio = P_oligo/P_bg;
		LR.emplace_back(likelihood_ratio);
		oligo.emplace_back(similar_oligo);
		//cout << "Oligo: "<< similar_oligo << "\tProbabiliy oligo: " << P_oligo << "\tProbability background: " << P_bg << "\tLR: " << likelihood_ratio << endl;
		
		/*
		 *The like_ratio_map is a map where for each oligo present in the vertical map
		 *we couple the likelihood ratio previously calculated with the ratio between the probability
		 *to have the oligo in the PWM_hamming and the background probability
		 */
	}

	double sum = 0;
	double Nsites = ceil((similar_oligos_map.size()/2) + 1);

	for (unsigned int i = 0; i < similar_oligos_map.size(); i++){
		sum += LR[i];
	}

	for (unsigned int i = 0; i < similar_oligos_map.size(); i++){
		LR[i] = LR[i]/sum;
		LR[i] = LR[i] * Nsites;
		like_ratio_map.insert(pair<string, double>(oligo[i], LR[i]));
	}


	//matrix::squash()

	double total = 0;
	
	bool renorm = true;

	for(map<string, double >::iterator it = like_ratio_map.begin();it != like_ratio_map.end(); ++it){
		total += it->second;
	}
	while(renorm)
    {
    	renorm = false;

        double norm = total/Nsites;
        total = 0;

        for(map<string, double >::iterator it = like_ratio_map.begin();it != like_ratio_map.end(); ++it){
            
			double p = it->second;

            if(p < 1){
				p /= norm; 
			}                            

            if(p > 1)
            {
                p = 1;
                Nsites--;
                renorm = true;
            }

            it->second = p;

            if(p <= 1){
				total += p;
			}              
        }
    }
	
	LR.clear();
	oligo.clear();
}

//In this function there is the second part of the EM algorithm and this is the maximization part
void hamming_class::EM_Mpart(unsigned int position,unordered_map<string,unsigned int> &orizzontal_map_plus){
	//PROFILE_FUNCTION();
	/*
	cout << "---------------------------------" << endl;
	cout << "#POSITION " << position + 1 << " before Mpart"<< endl;
	cout << "---------------------------------" << endl;
	*/

	//matrix::get_p()
	
	for(unsigned int x = 0; x < PWM_hamming.size(); x++)
	{
		for(unsigned int y = 0; y < PWM_hamming[0].size(); y++)
		{
			PWM_hamming[x][y] = 0;
		}
	}
	for (unsigned int b = 0; b < PWM_hamming[0].size(); b++){
		for(map<string, double >::iterator it = like_ratio_map.begin();it != like_ratio_map.end(); ++it){
			switch(it->first[b]){
				case 'A':			
					PWM_hamming[0][b] = PWM_hamming[0][b] + it->second;
					break;

				case 'C':
					PWM_hamming[1][b] = PWM_hamming[1][b] + it->second;	
					break;

				case 'G':
					PWM_hamming[2][b] = PWM_hamming[2][b] + it->second;
					break;

				case 'T':
					PWM_hamming[3][b] = PWM_hamming[3][b] + it->second;
					break;

				default:				//Case if there is N
					break;
			}
		}
	}

	double sum = 0;

	for(unsigned int j = 0; j < PWM_hamming.size(); j++){
		sum += PWM_hamming[j][0];
	}

	for (unsigned short int i = 0; i<PWM_hamming.size(); i++){
		for (unsigned short int j = 0; j<PWM_hamming[i].size(); j++){
			PWM_hamming[i][j] = PWM_hamming[i][j]/sum;
		}
	}

/*
	for (unsigned int i = 0; i < PWM_hamming.size(); i++){
		for (unsigned int j = 0; j < PWM_hamming[i].size(); j++){
			PWM_hamming[i][j] = 0;
		}
	}
	//For each position of the PWM_hamming we add the likelihood ratios to modify the previous PWM	
	for (unsigned int b = 0; b < PWM_hamming[0].size(); b++){
		double sum = 0;
		for(map<string, double >::iterator it = like_ratio_map.begin();it != like_ratio_map.end(); ++it){
			switch(it->first[b]){
				case 'A':			
					PWM_hamming[0][b] = PWM_hamming[0][b] + it->second;
					break;

				case 'C':
					PWM_hamming[1][b] = PWM_hamming[1][b] + it->second;	
					break;

				case 'G':
					PWM_hamming[2][b] = PWM_hamming[2][b] + it->second;
					break;

				case 'T':
					PWM_hamming[3][b] = PWM_hamming[3][b] + it->second;
					break;

				default:				//Case if there is N
					break;
			}
		}
	}
	*/

	/*
	//And here we divide each position for the sum of all the likelihood ratio
	for (unsigned int i = 0; i<PWM_hamming.size(); i++){
		for (unsigned int j = 0; j < PWM_hamming[i].size(); j++){
	         PWM_hamming[i][j] = PWM_hamming[i][j]/sum_vect[j];
		}
	}
	
	cout << "---------------------------------" << endl;
	cout << "#POSITION " << position + 1 << " after Mpart"<< endl;
	cout << "---------------------------------" << endl;
	
	print_PWM("", PWM_hamming);
	*/
}
// This function is made to check if the EM_cycle reaches convergence
bool hamming_class::EM_convergence(vector<vector<double>>& PWM_old, vector<vector<double>> &PWM_hamming, bool conv){

	vector<vector<double>> PWM_old_conv;
	vector<vector<double>> PWM_hamming_conv;

	PWM_old_conv = PWM_old;
	PWM_hamming_conv = PWM_hamming;

	conv = false;
	for (unsigned int i = 0; i < PWM_hamming.size(); i++)
	{
        for (unsigned int j = 0; j < PWM_hamming[0].size(); j++)
		{
			// if (PWM_old_conv[i][j] < 0.00001){
			// 	PWM_old_conv[i][j] = 0;
			// }
			// if (PWM_hamming_conv[i][j] < 0.00001){
			// 	PWM_old_conv[i][j] = 0;
			// }
            if (abs(PWM_old_conv[i][j] - PWM_hamming_conv[i][j])>0.001){
                conv = true;
				break;
			}
		}
	}
	
	return conv;
}

void hamming_class::EM_cycle(vector<bed_class> &GEP, unordered_map<string,unsigned int>& orizzontal_map_minus, unsigned int position, unordered_map<string,unsigned int> &orizzontal_map_plus){
	//PROFILE_FUNCTION();
	bool conv = true;
	int i = 0;

	vector<vector<double>> PWM_old;
	
	//In this cycle we repeat the EM until the convergence is reached
	//for (unsigned int i = 0; i < exp_max; i++){ 
	if (exp_max == "c"){

		while(conv && i < 200){
			
			PWM_old = PWM_hamming;

			EM_Epart(GEP, orizzontal_map_minus, position, orizzontal_map_plus);
			EM_Mpart(position, orizzontal_map_plus);
		
			like_ratio_map.clear();

			conv = EM_convergence(PWM_old, PWM_hamming, conv);
			PWM_old = PWM_hamming;
			i++;
		}

	}
	else{

		double em = stod(exp_max);
		for(unsigned int i = 0; i < em; i++){
			EM_Epart(GEP, orizzontal_map_minus, position, orizzontal_map_plus);
			EM_Mpart(position, orizzontal_map_plus);
		
			like_ratio_map.clear();
		}

	}
}

//Shifting the PWM_matrix on the sequences and calculate local scores (from positon where the matrix has been generated), and the global scores from each sequences positions
void z_test_class::oligos_vector_creation_PWM(vector<bed_class> &GEP){
	//PROFILE_FUNCTION();
	//For every sequence
	for(unsigned int i=0; i<GEP.size(); i++){
		
		//string sequence = GEP[i].sequence;
		
		//Calling oligo class to accede to all functions useful to shift a matrix on sequences --> Shifting on FWD strand
		oligo_class SHIFTING_PWM(matrix_log, GEP[i].sequence);

		//Return oligo scores calculated from previous shifting
		oligo_scores_orizzontal_FWD = SHIFTING_PWM.oligo_scores;
	
		//If analysis is in Double strand
		if(DS == 1){
			
			//Make the shifting also on reverse strand, putting as input the inverse_log_matrix --> Shifting on REV strand
			oligo_class SHIFTING_PWM_2(inverse_matrix_log, GEP[i].sequence);

			//Retrun oligo scores from previous shifting
			oligo_scores_orizzontal_REV = SHIFTING_PWM_2.oligo_scores;

			//Select the best scores between FWD and REV strand (for each position)
			check_best_strand_oligo();	

			//Fill the local scores vector with scores found in position where the matrix has been generated
			all_local_scores.emplace_back(oligo_scores_orizzontal_BEST[local_pos-1]);			
			//Fill the global scores vector with all scores generated from shifting and selected from check_best_strand_oligo function
			all_global_scores.insert(all_global_scores.end(), oligo_scores_orizzontal_BEST.begin(), oligo_scores_orizzontal_BEST.end());
		}

		//If analysis is in Single strand
		else{
			
			//Local best scores are all from FWD strand
			all_local_scores.emplace_back(oligo_scores_orizzontal_FWD[local_pos-1]);	
			all_global_scores.insert(all_global_scores.end(), oligo_scores_orizzontal_FWD.begin(), oligo_scores_orizzontal_FWD.end());
		}
		
		//Clearing of orizzontal best score for the next sequence cycle
		oligo_scores_orizzontal_BEST.clear();	
		//oligo_scores_orizzontal_BEST.shrink_to_fit();
	}	
	
}

//Function to select the best scores between FWD and REV strand (for each sequence position)
void z_test_class::check_best_strand_oligo(){
	//PROFILE_FUNCTION();
	//For all oligo scores 
	for(unsigned int oligo = 0; oligo < oligo_scores_orizzontal_FWD.size(); oligo ++){
		
		//If FWD is better than REV
		if(oligo_scores_orizzontal_FWD[oligo] >= oligo_scores_orizzontal_REV[oligo]){
			
			oligo_scores_orizzontal_BEST.emplace_back(oligo_scores_orizzontal_FWD[oligo]);
		}
		
		//If REV is better than FWD
		else{

			oligo_scores_orizzontal_BEST.emplace_back(oligo_scores_orizzontal_REV[oligo]);
		}
	}
}

//Function to calculate all the parameters useul to z-score calculation
void z_test_class::z_score_parameters_calculation(){
	//PROFILE_FUNCTION();
	double local_sum = accumulate(all_local_scores.begin(), all_local_scores.end(),0.0);	
	double global_sum = accumulate(all_global_scores.begin(), all_global_scores.end(), 0.0);
	double tot_sq_sum_global = inner_product(all_global_scores.begin(), all_global_scores.end(), all_global_scores.begin(), 0.0);

	double tot_sq_sum_local = inner_product(all_local_scores.begin(), all_local_scores.end(), all_local_scores.begin(), 0.0);
	global_mean = global_sum/all_global_scores.size();
	local_mean = local_sum/all_local_scores.size();
	global_dev_std = sqrt(tot_sq_sum_global/all_global_scores.size() - global_mean * global_mean);
	local_dev_std = sqrt(tot_sq_sum_local/all_local_scores.size() - local_mean * local_mean);
	
}

//Z-score calculation function
void z_test_class::z_score_calculation(){
	//PROFILE_FUNCTION();

	z_score = ((global_mean - local_mean)/ (global_dev_std / sqrt(all_local_scores.size()))); 

	const double Z  = z_score;
	Zpvalue = gsl_cdf_ugaussian_P(Z);

	Zpvalue = check_p_value(Zpvalue);


}

//Function to check, given an oligo as input, if this oligo is palindrome or not
bool check_palindrome(string bases,string& reverse_bases){
	reverse_bases.clear();
	//For any character of the string insert into another string (called reverse bases) the complementary character
	for(int i=bases.length()-1; i>=0; i--){

		char base;
		base = bases[i];
		switch (base) {

			case 'A' : reverse_bases.append("T"); 
				   break;
			case 'T' : reverse_bases.append("A"); 
				   break;
			case 'G' : reverse_bases.append("C"); 
				   break;
			case 'C' : reverse_bases.append("G"); 
				   break;
			case 'N' : reverse_bases.append("N"); 
				   break;
		}
	}
	
	//If they are equal --> it means that the oligo "bases" is palindrome
	if (reverse_bases == bases){
		return true;
	}
	else {
		return false;
	}
	
}

//If the p value is rounded to 0 assigne it a standar low value of 1.000001e-300 to avoid possible future errors
double check_p_value(double p){
	if(p == 0){

		p = 1.000001e-300;
	}

	return p;
}

/////DEBUG/////////////////////////////////////////////////////////

vector<vector<double>> matrix_class::return_inverse_log_matrix(){
	//PROFILE_FUNCTION();
	return inverse_matrix_log;
}

vector<vector<double>> matrix_class::return_log_matrix(){
	//PROFILE_FUNCTION();
	return matrix_log;
}

vector<vector<double>> matrix_class::return_norm_matrix(){
	//PROFILE_FUNCTION();
	return norm_matrix;
}

//Debug function: Print sequences and coordinates from GEP vector into a .fasta file to check if the sequences extraction is correct
void coordinator_class::print_GEP(vector<bed_class> &GEP){
	//PROFILE_FUNCTION();
	//Twobit_JASPAR_Bed used to create GEP vector saved into alias file to name the outputs	
	alias_file = (TWOBIT_FILE.erase(0,TWOBIT_FILE.find_last_of("/")+1)+"_"+ JASPAR_FILE.erase(0,JASPAR_FILE.find_last_of("/")+1)+"_"+ BED_FILE.erase(0,BED_FILE.find_last_of("/")+1));

	//Output file .bed carrying the centered coordinates
	ofstream outfile;	
	JASPAR_FILE = JASPAR_FILE.erase(JASPAR_FILE.find_last_of("."), JASPAR_FILE.size());
	outfile.open(alias_file);

	for(unsigned int i=0; i<GEP.size(); i++){

		outfile << GEP[i].chr_coord << ":" << GEP[i].start_coord << "-" << GEP[i].end_coord << endl;
	}

	outfile.close();

	//Output file .fasta carrying the centered coordinates and the sequences extracted
	BED_FILE = BED_FILE.erase(BED_FILE.find_last_of("."), BED_FILE.size());
	outfile.open(alias_file+".fasta");

	for(unsigned int i=0; i<GEP.size(); i++){

		outfile << ">" << GEP[i].chr_coord << ":" << GEP[i].start_coord << "-" << GEP[i].end_coord << endl;	
		outfile << GEP[i].sequence << endl;

	}

	outfile.close();
}

//If the analysis is on Double Strand and oligos need to be ordered by p_value this is the function which writes on the output file
void p_value_class::print_debug_p_value_DS(map<pair<string,string>,pair<unsigned int,unsigned int>> &pair_map, unsigned int position, ofstream& outfile, unsigned int freq){
	//PROFILE_FUNCTION();
	unsigned int c=0;
	sum_top_N = 0;

	multimap<pair<string,string>,pair<unsigned int, unsigned int>>::iterator it_multi;
	
	//Scrollig the p_value_sort multimap from first element to -nth element assign the rigth value to each feature and print them
	for(multimap<double,pair<string,string>>::iterator it_pair = p_value_sort.begin(); it_pair!=p_value_sort.end() && c<top_N; it_pair++, c++){

		double FREQ, Sum_Occ_Oligo;
		it_multi = pair_map.find(it_pair->second);
		string Oligo = it_multi->first.first;
		string Oligo_RC = it_multi->first.second;
		unsigned int Num_Occ_FWD = it_multi->second.first;
		unsigned int Num_Occ_REV = it_multi->second.second;
		bool pal = check_palindrome(Oligo, reverse_bases);
		unsigned int Rank = c;
		string PAL;
		unsigned int Num_Occ_RC_FWD, Num_Occ_RC_REV, Sum_Occ_RC;
		double P_VAL = it_pair->first;

		if(pal == 1){

			PAL = "TRUE";
			Num_Occ_REV = Sum_Occ_Oligo = Num_Occ_RC_FWD = Num_Occ_RC_REV = Sum_Occ_RC = Num_Occ_FWD;
			FREQ = Sum_Occ_Oligo/freq;
		}

		else{
			PAL = "FALSE";
			Sum_Occ_Oligo = Num_Occ_FWD + Num_Occ_REV;
			Num_Occ_RC_FWD = Num_Occ_REV;
			Num_Occ_RC_REV = Num_Occ_FWD;
			Sum_Occ_RC = Num_Occ_RC_FWD + Num_Occ_RC_REV;
			FREQ = Sum_Occ_Oligo/freq;
		}

		outfile << position+1 << "\t" << Rank+1 << "\t";
		outfile << Oligo << "\t" << Num_Occ_FWD << "\t" << Num_Occ_REV << "\t" << Sum_Occ_Oligo << "\t";
		outfile << Oligo_RC << "\t" << Num_Occ_RC_FWD << "\t" << Num_Occ_RC_REV << "\t" << Sum_Occ_RC << "\t";
		outfile << PAL << "\t" << Sum_Occ_Oligo << "\t" << FREQ << "\t" << P_VAL << endl;	
		sum_top_N = sum_top_N + Sum_Occ_Oligo;

	}
}

//If the analysis is on Single Strand and oligos need to be ordered by p_value this is the function which writes on the output file
void p_value_class::print_debug_p_value_SS(map<pair<string,string>,pair<unsigned int,unsigned int>> &pair_map, unsigned int position, ofstream& outfile, unsigned int freq){
	//PROFILE_FUNCTION();
	unsigned int c = 0;
	sum_top_N = 0;

	multimap<pair<string,string>,pair<unsigned int, unsigned int>>::iterator it_multi;

	//Scrollig the p_value_sort multimap from first element to -nth element assign the rigth value to each feature and print them
	for(multimap<double,pair<string,string>>::iterator it_pair = p_value_sort.begin(); it_pair!=p_value_sort.end() && c<top_N; it_pair++, c++){

		double FREQ, Num_Occ_FWD;
		it_multi = pair_map.find(it_pair->second);
		string Oligo = it_multi->first.first;
		unsigned int Num_Occ_Oligo = it_multi->second.first;
		bool pal = check_palindrome(Oligo, reverse_bases);
		unsigned int Rank = c;
		string PAL;
		double P_VAL = it_pair->first;
		Num_Occ_FWD = Num_Occ_Oligo;
		FREQ = Num_Occ_FWD/freq;

		if(pal == 0){

			PAL = "FALSE";
		}
		else{	PAL = "TRUE";}

		outfile << position+1 << "\t" << Rank+1 << "\t";
		outfile << Oligo << "\t" << Num_Occ_FWD  << "\t";
		outfile << PAL << "\t" << FREQ << "\t" << P_VAL << endl;	
		sum_top_N = sum_top_N + Num_Occ_FWD;
		
	}
}

//If the analysis is on Double Strand and oligos need to be ordered by occurrences this is the function which writes on the output file
void p_value_class::print_debug_occurrences_DS(map<pair<string,string>,pair<unsigned int,unsigned int>> &pair_map, unsigned int position, ofstream& outfile, unsigned int freq, vector<double> p_value_vec){
	//PROFILE_FUNCTION();
	unsigned int c=0;
	sum_top_N = 0;
	multimap<unsigned int,pair<string,string>> sum_occurrences_multimap;
	multimap<pair<string,string>,pair<unsigned int,unsigned int>>::iterator pair_map_it = pair_map.begin();
	
	sum_occurrences_multimap = creating_sum_occurrences_multimap(vertical_multimap);

	//Scrollig the p_value_sort multimap from first element to -nth element assign the rigth value to each feature and print them
	for(multimap<unsigned int,pair<string,string>>::reverse_iterator it_rev = sum_occurrences_multimap.rbegin(); it_rev!=sum_occurrences_multimap.rend() && c<top_N; it_rev++, c++){
		
		pair_map_it = pair_map.find(make_pair(it_rev->second.first, it_rev->second.second));

		double FREQ, Sum_Occ_Oligo;
		string Oligo = pair_map_it->first.first;
		string Oligo_RC = pair_map_it->first.second;
		unsigned int Num_Occ_FWD = pair_map_it->second.first;
		unsigned int Num_Occ_REV = pair_map_it->second.second;
		bool pal = check_palindrome(Oligo, reverse_bases);
		unsigned int Rank = c;
		string PAL;
		unsigned int Num_Occ_RC_FWD, Num_Occ_RC_REV, Sum_Occ_RC;
		double P_VAL = p_value_vec[c];

		if(pal == 1){

			PAL = "TRUE";
			Num_Occ_REV = Sum_Occ_Oligo = Num_Occ_RC_FWD = Num_Occ_RC_REV = Sum_Occ_RC = Num_Occ_FWD;
			FREQ = Sum_Occ_Oligo/freq;
		}

		else{
			PAL = "FALSE";
			Sum_Occ_Oligo = Num_Occ_FWD + Num_Occ_REV;
			Num_Occ_RC_FWD = Num_Occ_REV;
			Num_Occ_RC_REV = Num_Occ_FWD;
			Sum_Occ_RC = Num_Occ_RC_FWD + Num_Occ_RC_REV;
			FREQ = Sum_Occ_Oligo/freq;
		}

		outfile << position+1 << "\t" << Rank+1 << "\t";
		outfile << Oligo << "\t" << Num_Occ_FWD << "\t" << Num_Occ_REV << "\t" << Sum_Occ_Oligo << "\t";
		outfile << Oligo_RC << "\t" << Num_Occ_RC_FWD << "\t" << Num_Occ_RC_REV << "\t" << Sum_Occ_RC << "\t";
		outfile << PAL << "\t" << Sum_Occ_Oligo << "\t" << FREQ << "\t" << P_VAL << endl;	
		sum_top_N = sum_top_N + Sum_Occ_Oligo;

	}
}


multimap<unsigned int,pair<string,string>> p_value_class::creating_sum_occurrences_multimap(multimap<pair<unsigned int,unsigned int>,pair<string,string>>& pair_map){

	multimap<unsigned int,pair<string,string>> sum_occurrences_multimap;

	for(multimap<pair<unsigned int,unsigned int>,pair<string,string>>::reverse_iterator rev_it = pair_map.rbegin(); rev_it != pair_map.rend(); rev_it ++){

		sum_occurrences_multimap.insert({(rev_it->first.first + rev_it->first.second),make_pair(rev_it->second.first,rev_it->second.second)});
	}

	return sum_occurrences_multimap;
}

//If the analysis is on Single Strand and oligos need to be ordered by occurrences this is the function which writes on the output file
void p_value_class::print_debug_occurrences_SS(map<pair<string,string>,pair<unsigned int,unsigned int>> &pair_map, unsigned int position, ofstream& outfile, unsigned int freq, vector<double> p_value_vec){
	//PROFILE_FUNCTION();
	unsigned int c = 0;
	sum_top_N = 0;

	//Scrollig the p_value_sort multimap from first element to -nth element assign the rigth value to each feature and print them
	for(multimap<pair<unsigned int,unsigned int>,pair<string,string>>::reverse_iterator it_rev = vertical_multimap.rbegin(); it_rev!=vertical_multimap.rend() && c<top_N; it_rev++, c++){

		double FREQ, Num_Occ_FWD;
		string Oligo = it_rev->second.first;
		unsigned int Num_Occ_Oligo = it_rev->first.first;
		bool pal = check_palindrome(Oligo,reverse_bases);
		unsigned int Rank = c;
		string PAL;
		double P_VAL = p_value_vec[c];
		Num_Occ_FWD = Num_Occ_Oligo;
		FREQ = Num_Occ_FWD/freq;

		if(pal == 0){

			PAL = "FALSE";
		}
		else{	PAL = "TRUE";}

		outfile << position+1 << "\t" << Rank+1 << "\t";
		outfile << Oligo << "\t" << Num_Occ_FWD  << "\t";
		outfile << PAL << "\t" << FREQ << "\t" << P_VAL << endl;	
		sum_top_N = sum_top_N + Num_Occ_FWD;
		
	}
}

//Function to create an output file of orizzontal map. Oligos are ranked by their occurrences
void map_class::print_debug_orizzontal(){
	//PROFILE_FUNCTION();
	for(unsigned int i=0; i<orizzontal_plus_debug.size(); i++){

		ofstream outfile;

		if(DS==1){
			outfile.open(to_string(kmers_vector[i])+"-mers_occurrences_"+alias_file+"_DS.txt");	
		}

		else{
			outfile.open(to_string(kmers_vector[i])+"-mers_occurrences_"+alias_file+"_SS.txt");	
		}

		multimap<unsigned int,string> orizzontal_output;

		for (unordered_map<string,unsigned int>::iterator it = orizzontal_plus_debug[i].begin() ; it != orizzontal_plus_debug[i].end(); it++ ){
			orizzontal_output.insert({it->second, it->first});	
		}	      

		for (multimap<unsigned int,string>::reverse_iterator it_rev = orizzontal_output.rbegin(); it_rev!=orizzontal_output.rend(); it_rev++){

			if(DS==1){
				reverse_bases.clear();	
				//reverse_bases.shrink_to_fit();

				bool palindrome = check_palindrome(it_rev->second, reverse_bases);

				if(!palindrome){

					unordered_map<string,unsigned int>::iterator find_RC = orizzontal_minus_debug[i].find(reverse_bases);
					outfile << it_rev->second << "\t" << it_rev->first << "\t" << find_RC->first << "\t" << find_RC->second << "\t\n";

				}

				else{
					outfile << it_rev->second << "\t" << it_rev->first <<  endl;
				}
			}

			else{
				outfile << it_rev->second << "\t" << it_rev->first <<  endl;
			}
		}
		outfile.close();
	}
}

//Function to write the correct header in output positional file, following analysis parameters
ofstream map_class::outfile_header(unsigned int j){
	//PROFILE_FUNCTION();

	ofstream outfile;

	if(DS==1){
		
		if(ordering == "p"){
		
			outfile.open(to_string(kmers_vector[j])+"-mers_positional_occurrences_"+alias_file+"DS_p_val.txt");
		}
		else{
		
			outfile.open(to_string(kmers_vector[j])+"-mers_positional_occurrences_"+alias_file+"DS_occ.txt");
		}
		
		outfile << "#Maps vector with kmers occurences (Double Strand) counted for positions in sequence (for k = " << kmers_vector[j] << "):\n";
		outfile << "#Position" << "\t" << "Rank" << "\t" << "Oligo" << "\t" << "Num_Occ_FWD" << "\t" << "Num_Occ_REV" << "\t" << "Sum_Occ_Oligo" << "\t" << "Oligo_RC" << "\t" << "Num_Occ_RC_FWD" << "\t" << "Num_Occ_RC_REV" << "\t" << "Sum_Occ_RC" << "\t" << "PAL" << "\t" << "Tot_Occ" << "\t" << "FREQ" << "\t" << "P_VALUE\n";

	}

	else{
		
		if(ordering == "p"){
		
			outfile.open(to_string(kmers_vector[j])+"-mers_positional_occurrences_"+alias_file+"SS_p_val.txt");
		}
		else{
		
			outfile.open(to_string(kmers_vector[j])+"-mers_positional_occurrences_"+alias_file+"SS_occ.txt");
		}

		outfile << "#Maps vector with kmers occurences (Single Strand) counted for positions in sequence (for k = " << kmers_vector[j] << "):\n";
		outfile << "#Position" << "\t" << "Rank" << "\t" << "Oligo" << "\t" << "Num_Occ_Oligo" << "\t" << "PAL" << "\t" << "FREQ" << "\t" << "P_VALUE\n";

	}
	return outfile;	
}

//Function to print Top N sum and their frequences agains total oligo number in position in Output file
void map_class::TopN_sum_and_freq(){
	//PROFILE_FUNCTION();
	ofstream outfile;

	for(unsigned int i=0; i<tot_sum_matrix.size(); i++){

		if(DS==1 && ordering == "p"){	
			outfile.open(to_string(kmers_vector[i])+"-mers_Top"+to_string(top_N)+"_sum_and_frequence_DS_p_val.txt");
		}
		
		else if(DS==1 && ordering != "p"){	
			outfile.open(to_string(kmers_vector[i])+"-mers_Top"+to_string(top_N)+"_sum_and_frequence_DS_occ.txt");
		}

		else if(DS==0 && ordering == "p"){	
			outfile.open(to_string(kmers_vector[i])+"-mers_Top"+to_string(top_N)+"_sum_and_frequence_SS_p_val.txt");
		}

		else{
			outfile.open(to_string(kmers_vector[i])+"-mers_Top"+to_string(top_N)+"_sum_and_frequence_SS_occ.txt");
		}

		outfile << "###Top " << top_N << " occurrences sum with k = " << kmers_vector[i] << ":\n"; 
		outfile << "Position" << "\t" << "Sum" << "\t" << "Frequences\n"; 

		for(unsigned int j=0; j<tot_sum_matrix[i].size(); j++){

			double sum = tot_sum_matrix[i][j];		//Put as double to dont loose the precision	
			double frequence = sum/tot_freq_matrix[i][j];	
			outfile << j+1 << "\t" << tot_sum_matrix[i][j] << "\t" << frequence << endl; 

		}
		outfile.close();
	}
}

//Function to print the parameters used to calculate p-values. For each oligo ranked all its parameters are visualized
void map_class::p_value_parameters_debug_p_val(){
	//PROFILE_FUNCTION();
	ofstream outfile;
	
	//For each k-mer in input create an output file
	for(unsigned int j = 0; j<P_VALUE_MATRIX.size(); j++){
		
		//Different names to differentiate if the analysis has been made on Double strand or Single strand
		if(DS == 1){
		outfile.open(to_string(kmers_vector[j])+"-mers_p_value_parameters_control_p_val_"+BED_FILE+"DS.txt");
		}
		else{
		outfile.open(to_string(kmers_vector[j])+"-mers_p_value_parameters_control_p_val_"+BED_FILE+"SS.txt");
		}
		
		//Defining the file header 
		outfile << "#Parameters used to calculate p_value for each oligo positionally ranked\n";
		outfile << "#Position\tRank\tOligo\tK\tN1\tN2\tT\tP_VALUE\tP_VALUE_LOG10\n";

		//For each position
		for(unsigned int i = 0; i<P_VALUE_MATRIX[j].size(); i++){
			
			//Setting the rank to 0
			unsigned int rank = 0;

			//Printing into the file the parameters used to calculate each oligo p-value in the top N ranking
			for(multimap<double,pair<string,vector<unsigned int>>>::iterator it = P_VALUE_MATRIX[j][i].p_value_KNT.begin(); it != P_VALUE_MATRIX[j][i].p_value_KNT.end() && rank<top_N; it++, rank++){

				outfile << i+1 << "\t";
				outfile << rank+1 << "\t";
				outfile << it->second.first << "\t";
				outfile << it->second.second[0] << "\t";
				outfile << it->second.second[1] << "\t";
				outfile << it->second.second[2] << "\t";
				outfile << it->second.second[3] << "\t";
				outfile << it->first << "\t";
				outfile << log10(it->first)*-1 << endl;
			}

		}
		outfile.close();
	}
}

void map_class::p_value_parameters_debug_occ(){
	//PROFILE_FUNCTION();
	ofstream outfile;
	
	for(unsigned int j = 0; j<P_VALUE_MATRIX.size(); j++){
		
		if(DS == 1){
		outfile.open(to_string(kmers_vector[j])+"-mers_p_value_parameters_control_occ_"+BED_FILE+"DS.txt");
		}
		else{
		outfile.open(to_string(kmers_vector[j])+"-mers_p_value_parameters_control_occ_"+BED_FILE+"SS.txt");
		}

		outfile << "#Parameters used to calculate p_value for each oligo positionally ranked\n";
		outfile << "#Position" << "\t" << "Rank" << "\t" << "Oligo" << "\t" << "K" << "\t" << "N1" << "\t" << "N2" << "\t" << "T" << "\t" << "P_VALUE" << "\t" <<"P_VALUE_LOG10\n";

		for(unsigned int i = 0; i<P_VALUE_MATRIX[j].size(); i++){

			unsigned int rank = 0;

			for(multimap<pair<unsigned int,unsigned int>,pair<string,string>>::reverse_iterator it_rev = P_VALUE_MATRIX[j][i].vertical_multimap.rbegin(); it_rev != P_VALUE_MATRIX[j][i].vertical_multimap.rend() && rank<top_N; it_rev++, rank++){

				outfile << i+1 << "\t";
				outfile << rank+1 << "\t";
				outfile << it_rev->second.first << "\t";
				outfile << P_VALUE_MATRIX[j][i].K_vec[rank] << "\t";
				outfile << P_VALUE_MATRIX[j][i].N1_vec[rank] << "\t";
				outfile << P_VALUE_MATRIX[j][i].N2_vec[rank] << "\t";
				outfile << P_VALUE_MATRIX[j][i].T << "\t";
				outfile << P_VALUE_MATRIX[j][i].p_value_vec[rank] << "\t";
				outfile << log10(P_VALUE_MATRIX[j][i].p_value_vec[rank])*-1 << endl;
			}

		}
		outfile.close();
	}
}

//Function to open and print the header of hamming output file, following the analysis parameters
ofstream map_class::outfile_header_hamming(unsigned int j){
	//PROFILE_FUNCTION();

	ofstream outfile;

	if(DS==1){
		outfile.open(to_string(kmers_vector[j])+"-mers_hamming_"+alias_file+"DS.txt");
	}
	else{
		outfile.open(to_string(kmers_vector[j])+"-mers_hamming_"+alias_file+"SS.txt");
	}
			
	outfile << "#For each best oligo in positions, a table representing hamming distance oligos and 2 different frequences(FREQUENCE_1 = N.occurrences/N.sequences | FREQUENCE_2 = N.occurrences/Orizzontal_occurrences) (for k = " << kmers_vector[j] << "):\n";
	outfile << "#Position" << "\t" << "Best_oligo" << "\t" << "Num_Occ_Best_oligo" << "\t" << "Humming_found_oligo_number" << "\t" << "Total_occurrences(best+hamming)" << "\t" << "FREQUENCE_1" << "\t" << "FREQUENCE_2\n";

	return outfile;	
}


//Print debug hamming features for each position in hamming output file passed as reference
void hamming_class::print_debug_hamming(unsigned int position, ofstream& outfile){
	//PROFILE_FUNCTION();
	outfile << position+1 << "\t" << real_best_oligo << "\t" << real_best_oligo_occurrences << "\t" << similar_oligos.size()-1 << "\t" << tot_similar_occurrences << "\t" << FREQUENCE_1 << "\t" << FREQUENCE_2 << endl;

}

//Print debug for PWM_hamming outfile -> Selection of output filename
void map_class::Outfile_PWM_matrices(){
	//PROFILE_FUNCTION();
	ofstream outfile;
	
	for(unsigned int j=0; j<kmers_vector.size(); j++){

		if(DS == 1){
			
			outfile.open(to_string(kmers_vector[j])+"-mers_PWM_hamming_matrices_"+alias_file+"DS.txt");
			if (tomtom){
				print_debug_PWM_hamming_tomtom(outfile, j, kmers_vector[j]);
			}
			else
			{
				print_debug_PWM_hamming(outfile, j, kmers_vector[j]);
			}			
			outfile.close();
		}

		else{
			outfile.open(to_string(kmers_vector[j])+"-mers_PWM_hamming_matrices_"+alias_file+"SS.txt");

			if (tomtom){
				print_debug_PWM_hamming_tomtom(outfile, j, kmers_vector[j]);
			}
			else
			{
				print_debug_PWM_hamming(outfile, j, kmers_vector[j]);
			}			
			outfile.close();
		}

	}
}

void map_class::print_debug_PWM_hamming_tomtom(ofstream& outfile, unsigned int j, unsigned int k){
	//PROFILE_FUNCTION();

	vector<vector<double>> PWM_hamming;
	string ACGT = "ACGT";
	for(unsigned int position = 0; position < Z_TEST_MATRIX[j].size(); position++){
	
		PWM_hamming = HAMMING_MATRIX[j][Z_TEST_MATRIX[j][position].local_pos-1].PWM_hamming;
		
		outfile << ">Position" << Z_TEST_MATRIX[j][position].local_pos << " " << Z_TEST_MATRIX[j][position].Zpvalue <<endl; 


		for(unsigned int i = 0; i< PWM_hamming.size(); i++){
			outfile << ACGT[i] << "\t" << "[" << "\t";
			for(unsigned int j = 0; j<PWM_hamming[i].size(); j++){
				outfile << PWM_hamming[i][j] << "\t";
			}
			outfile << "]\n";
		}
		
		//outfile << endl << endl;
	}
	
}
//PWM_matrices, parameters to calculate z-score, z-score and p-value printing
void map_class::print_debug_PWM_hamming(ofstream& outfile, unsigned int j, unsigned int k){
	//PROFILE_FUNCTION();
	outfile << "#PWM Matrices calculated from the best oligo for each position and his hamming distanced oligos - k = " << k << endl << endl;

	vector<vector<double>> PWM_hamming;
	string ACGT = "ACGT";
	
	for(unsigned int position = 0; position < Z_TEST_MATRIX[j].size(); position++){
	
		PWM_hamming = HAMMING_MATRIX[j][Z_TEST_MATRIX[j][position].local_pos-1].PWM_hamming;
		
		outfile << "#Position " << Z_TEST_MATRIX[j][position].local_pos << ": \n#PWM calculated from oligo " << 
		HAMMING_MATRIX[j][Z_TEST_MATRIX[j][position].local_pos-1].real_best_oligo << " and his " <<
		HAMMING_MATRIX[j][Z_TEST_MATRIX[j][position].local_pos-1].similar_oligos.size()-1 << " hamming distanced neighbours.\n\n";

		for(unsigned int i = 0; i< PWM_hamming.size(); i++){
			
			outfile << ACGT[i] << "\t" << "[" << "\t";

			for(unsigned int j = 0; j<PWM_hamming[i].size(); j++){

				outfile << PWM_hamming[i][j] << "\t";
			}
			outfile << "]\n";
		}
		
		outfile << endl;
		outfile << "The global mean is: " << Z_TEST_MATRIX[j][position].global_mean << endl;
		outfile << "The global standard deviation is: " << Z_TEST_MATRIX[j][position].global_dev_std << endl;
		outfile << "The local mean is: " << Z_TEST_MATRIX[j][position].local_mean << endl;
		outfile << "The local standard deviation is: " << Z_TEST_MATRIX[j][position].local_dev_std << endl << endl;
		outfile << "The zscore calculated is: " << Z_TEST_MATRIX[j][position].z_score<< endl << endl;
		outfile << "The pvalue calculated from the Z score is: " << Z_TEST_MATRIX[j][position].Zpvalue<< endl << endl;
		outfile << "-------------------------------------------------------------------" << endl;	
	}
	
}

//Print debug for PWM_hamming outfile -> Selection of output filename
void map_class::Outfile_Z_score_values(){
	//PROFILE_FUNCTION();
	ofstream outfile;
	
	for(unsigned int j=0; j<kmers_vector.size(); j++){

		if(DS == 1){

			outfile.open(to_string(kmers_vector[j])+"-mers_Z_scores_"+alias_file+"DS.txt");
			print_debug_Z_scores(outfile, j, kmers_vector[j]);
			outfile.close();
		}

		else{

			outfile.open(to_string(kmers_vector[j])+"-mers_Z_scores_"+alias_file+"SS.txt");

			print_debug_Z_scores(outfile, j, kmers_vector[j]);
			outfile.close();
		}

	}
}

//PWM_matrices, parameters to calculate z-score, z-score and p-value printing
void map_class::print_debug_Z_scores(ofstream& outfile, unsigned int j, unsigned int k){
	//PROFILE_FUNCTION();
	outfile << "#Z_score parameters and p-value for hit positions - k = " << k << endl << endl;
	string best_oligo;
	
	outfile << "#Position" << "\t" << "best_oligo" << "\t" << "Local_mean" << "\t" << "Global_mean" << "\t" << "Local_std_dev" << "\t" << "Global_std_dev" << "\t" << "Z_score" << "\t" << "P-value" << "\t" << "P-value_Log10" << "\t" << "Bonferroni P-value" << "\t" << "Bonferroni_Log10\n";
	
	for(unsigned int position = 0; position < Z_TEST_MATRIX[j].size(); position++){
	
		double Zpvalue_Log10  = abs(log10(Z_TEST_MATRIX[j][position].Zpvalue));
		double bonferroni = Z_TEST_MATRIX[j][position].Zpvalue * ((half_length*2) - k);
		double bonferroni_Log10 = abs(log10(bonferroni));

		//best_oligo = HAMMING_MATRIX[j][local_pos-1].real_best_oligo;	

		outfile << Z_TEST_MATRIX[j][position].local_pos << "\t" << 
		HAMMING_MATRIX[j][Z_TEST_MATRIX[j][position].local_pos-1].real_best_oligo << 
		"\t" << Z_TEST_MATRIX[j][position].local_mean << "\t" << 
		Z_TEST_MATRIX[j][position].global_mean << "\t" << 
		Z_TEST_MATRIX[j][position].local_dev_std << 
		"\t" << Z_TEST_MATRIX[j][position].global_dev_std << "\t" << 
		Z_TEST_MATRIX[j][position].z_score << "\t" << Z_TEST_MATRIX[j][position].Zpvalue <<
		 "\t" << Zpvalue_Log10 << "\t" << bonferroni << "\t" << bonferroni_Log10 << endl;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////

////////////////////PARSER////////////////////////////////////////////////////////////////////

void command_line_parser(int argc, char** argv){
	
	const char* const short_opts = "hp:k:b:j:m:d:o:f:alrt:n:e:s";

	//Specifying the expected options
	const option long_opts[] ={
		{"help",      no_argument, nullptr,  'h' },
		{"param",      required_argument, nullptr,  'p' },
		{"ntop",      required_argument, nullptr,  'n' },
		{"kmer",   required_argument, nullptr,  'k' },
		{"all",   no_argument, nullptr,  'a' },
		{"tomtom", no_argument, nullptr, 'l'},
		{"refine",   no_argument, nullptr,  'r' },
		{"freq",   required_argument, nullptr,  'f' },
		{"distance",   required_argument, nullptr,  'd' },
		{"bed",    required_argument, nullptr,  'b' },
		{"ordering",    required_argument, nullptr,  'o' },
		{"jaspar",   required_argument, nullptr,  'j' },
		{"mf",   required_argument, nullptr,  'm' },
		{"twobit",   required_argument, nullptr,  't' },
		{"ss",   no_argument, nullptr,  's' },
		{"exp_maximization", required_argument, nullptr, 'e'},
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
			case 'p' : half_length = stoi(optarg); 
				   break;
			case 'n' : top_N = stoi(optarg); 
				   break;
			case 'b' : BED_FILE = string(optarg);
				   is_file_exist(BED_FILE, "--bed || -b");
				   break;
			case 'j' : JASPAR_FILE = string(optarg);
				   is_file_exist(JASPAR_FILE, "--jaspar || -j");
				   break;
			case 't' : TWOBIT_FILE = string(optarg);
				   is_file_exist(TWOBIT_FILE, "--twobit || -t ");
				   break;
			case 'o' : ordering = string(optarg);
				   if(ordering != "p"){
					   cerr << "ERROR: Wrong -o parameter inserted.\n\n\n";
					   display_help();
					   exit(1);}
				   break;
			case 'k' : kmers.clear();
				   kmers = string(optarg);
				   break;
			case 'a' : local_maxima_grouping = false;
				   break;
			case 'l' : 
				tomtom = true;
				break;
			case 'r' : refining_matrix = 1;
				   break;
			case 'e' : exp_max = string(optarg);
				   break;
			case 'f' : freq_treshold = stod(optarg);
				   if(freq_treshold == 0){

					   cout << "WARNING: frequency threshold 0 inserted\n";
				   }

				   if(freq_treshold < 0 || freq_treshold >= 1){
					   cerr << "ERROR: please insert a frequency treshold between 0 and 1.\n\n\n";
					   display_help();
					   exit(1);
				   }
				   break; 
			case 'd' : dist.clear();
				   dist = string(optarg);
				   break;
			case 's' : DS = 0;
				   break;
			case 'm' : MFASTA_FILE = string(optarg);
				   is_file_exist(MFASTA_FILE, "--mf || -m ");
				   break;
			case '?': // Unrecognized option
			default:
				   display_help();
				   break;
		}
	}
	check_input_file();
}

bool is_file_exist(string fileName, string buf){		//Input files existence control

	struct stat check;
	int regular_check, existing_check;
	const char * C_fileName = fileName.c_str();
	existing_check = stat(C_fileName, &check );

	regular_check = S_ISREG( check.st_mode );

	if ( regular_check == 0 || existing_check != 0) {
		cerr <<"ERROR: "<< buf << " file does not exist!\n\n";
		display_help();
		exit(1);	
	}
	return 0;
}

void check_input_file(){

	if(MFASTA_FILE.size() != 0 && (BED_FILE.size() != 0 || JASPAR_FILE.size() != 0 || TWOBIT_FILE.size() != 0)){

		cerr << "FATAL ERROR: Too many input arguments!\nPlease insert Multifasta file or Bed, Twobit and Jaspar file.\n\n";
	        display_help();
		exit(1);
	}
	if ((TWOBIT_FILE.size() == 0 ||  JASPAR_FILE.size() == 0 || BED_FILE.size() == 0) && MFASTA_FILE.size() == 0){
		cerr << "FATAL ERROR: some arguments needed \n"<<endl;	
	        display_help();
		exit(1);
	}
	
}

void display_help(){
	cerr << "\n --help || -h show this message\n";
	cerr << "\n --bed || -b <file_bed>: input bed file\n";
	cerr << "\n --kmer || -k <n1,n2,..,nN>: to select k-mers length for the analysis(DEFAULT: 6,8,10)\n" ;
	cerr << "\n --twobit || -t <file_twobit>: input twobit file\n";
	cerr << "\n --jaspar || -j <JASPAR_file>: input JASPAR file\n";
	cerr << "\n --param || -p <half_length>: half_length to select bases number to keep around the chip seq signal (DEFAULT: 150) \n";
	cerr << "\n --ntop || -n <number>: to decide the top n oligos to classify in positional sequence occurrences (DEFAULT: 10) \n";
	cerr << "\n --mf || -m <multifasta-file>: use multifasta instead of bed file [ -j,-b,-t,-p options not needed ]\n";
	cerr << "\n -s || --ss as input to make the analysis along the single strand. (DEFAULT: double strand)\n";
	cerr << "\n -o || --ordering 'p' to order the top N oligos by p-value and not by occurrences. (DEFAULT: ordering by occurrences)\n";
	cerr << "\n --distance || -d <n1,n2,...,nN> to select the hamming distances. (DEFAULT: 1,2,3)\n";
	cerr << "\n --freq || -f <n1> to set the frequence treshold to calculate the z_score. (DEFAULT: 0.02)\n";
	cerr << "\n --all || -a to disable the local maxima filtering\n"; 
	cerr << "\n --refine || -r to refine PWM matrices with secondary hamming\n"; 
	cerr << "\n --exp_maximization || -e to refine PWM matrices with the expectation maximization method, if you type a number this will be the number of cycles but if you want to reach convergence you can type just 'c'\n\n";
	cerr << "\n --tomtom || -t will give as output a format of matrices adapted for tomtom analysis\n\n";
	exit(EXIT_SUCCESS);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
