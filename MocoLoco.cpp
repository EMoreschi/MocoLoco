#include "MocoLoco.h" 

int main(int argc, char *argv[]){


	//If arguments number is 1 means that no input file has been inserted - display help
	if(argc == 1){

		display_help();
	}
	command_line_parser(argc, argv);
	GEP_path();

	return 0;
}

//Function to choose the pathway to follow. 2 input options:
//1) Bed-Twobit-Jaspar input
//2) Multifasta input
void  GEP_path(){

	//if the input is Bed-Twobit-Jaspar	
	if(MFASTA_FILE.size() == 0){	

		coordinator_class C; 
		
		//Create a .fasta file to check if the coordinates and the sequences extracted are correct
		C.print_debug_GEP(C.GEP);
		
		//Creating map class: input are GEP vector created from bed-twobit analysis, kmers, hamming distance
		map_class MAP(C.GEP,kmers,dist);
	}

	//else if the input is a Multifasta file
	else{

		multifasta_class MULTIFA(MFASTA_FILE);
		
		//Creating a Map class: input are GEP vector created from multifasta file analysis, kmers, hamming distance
		map_class MAP(MULTIFA.GEP,kmers,dist);
	}		
}

//Checking function to control if, for any k-mers inserted as input, there is a distance parameter
void map_class::check_kmer_dist(){

	if(kmers_vector.size() != distance_vector.size()){
	
		//If the number is not equal --> ERROR printed and help visualized
		cerr << "\nERROR: Please insert an equal number of k-mers and distance parameters!" << endl;
		display_help();
		exit(1);
	}
}

void bed_class::read_line(string line){ 

	//Split the line word by word and extract chromosome coordinates (chr, start, end)
	istringstream mystream(line);
	mystream >> chr_coord >> start_coord >> end_coord;		

}

void bed_class::centering_function ( unsigned int start,  unsigned int end, unsigned int half_length, const unsigned int overhead){

	unsigned int center = (start + end)/2;						

	//No overhead for start coordinates but overhead added to end coordinates
	start_coord = center - half_length;
	end_coord = center + half_length +overhead;

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

//Function to read BED and 2Bit files and create GEP (vector of bed class)
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

//Function to create oligos_vector (a oligo class vector)
void coordinator_class::oligos_vector_creation(vector<oligo_class> &oligos_vector, vector<vector<double>> matrix_log, vector<vector<double>> matrix_log_inverse, vector<bed_class> GEP){

	cout << "- [5] Analyzing sequences using Jaspar matrix\n";

	//For every sequences into GEP vector
	for(unsigned int i=0; i<GEP.size(); i++){

		string sequence = GEP[i].return_sequence(GEP[i]);
		string chr_coord = GEP[i].return_chr_coord_GEP();
		unsigned int start_coord = GEP[i].return_start_coord_GEP();

		//Calling the oligo_class constructor to analyze the shifting of the sequence on log_matrix (FWD strand analysis)
		oligo_class SHIFTING(matrix_log, sequence, chr_coord, start_coord, '+');

		//The oligo class just created is saved into oligos_vector (oligo_class vector)
		oligos_vector.emplace_back(SHIFTING);

		//If the analysis is on DS calling the oligo_class constructor to analyze the shifting of sequence on inverse_log_matrix (REVERSE strande analysis)
		if(DS == 1){

			oligo_class SHIFTING(matrix_log_inverse, sequence, chr_coord, start_coord, '-');
			oligos_vector.emplace_back(SHIFTING);
		}
	}	

	cout << "- [6] Selecting the best Jaspar's oligo for each sequence \n";
}

//Function useful, if the analysis is performed on DS, to choose from the best FWD strand oligo and the best REV strand oligo the best one to keep as "Best oligo" --> The oligos vector is divided in half and only the best strand for each sequence is kept
void coordinator_class::best_strand(){
	
	//Only if the analysis is on Double Strand
	if(DS == 1){

		vector<oligo_class> comparison;
		
		for(unsigned int i=0; i<oligos_vector.size(); i+=2){
			
			//The comparison is made by oligo_class in i position against the oligo class in i+1 position (The fwd and rev strand of the same sequence, which are consecutive into the oligos_vector)
			double best_score_norm_positive = oligos_vector[i].return_best_score_normalized();
			double best_score_norm_negative = oligos_vector[i+1].return_best_score_normalized();

			if(best_score_norm_positive >= best_score_norm_negative){

				comparison.emplace_back(oligos_vector[i]);
			}

			else{
				comparison.emplace_back(oligos_vector[i+1]);
			}
		}

		//The new oligos_vector is replaced by comparison vector, which contains only the best strand
		oligos_vector = comparison;
	}
}

//Function to calculate the score of a general oligo against a JASPAR matrix
void oligo_class::shifting(vector<vector<double>> matrix, string sequence, unsigned int s_iterator){

	double sum_scores = 0;
	
	//For each oligo in the current sequence a score is calculated
	if(s_iterator < sequence.size() - matrix[0].size() ) {

		for(unsigned int i=0; i< matrix[0].size(); i++){

			switch(sequence[i+s_iterator]){

				case 'A':

					sum_scores = sum_scores + matrix[0][i];
					break;

				case 'C':

					sum_scores = sum_scores + matrix[1][i];
					break;

				case 'G':

					sum_scores = sum_scores + matrix[2][i];
					break;

				case 'T':

					sum_scores = sum_scores + matrix[3][i];
					break;

				default:

					sum_scores = sum_scores + o_matrix_mins[i];
					break;
			}
		}
		
		//The total score of an oligo is saved into an oligo_scores vector
		oligo_scores.emplace_back(sum_scores);

		//Then the function is recalled recursively shifting on the next oligo thanks to the iterator progression
		shifting(matrix, sequence, s_iterator+1);
	}
}

//Function to read JASPAR PWM file, extract values and create a matrix class
vector<vector<double>> coordinator_class::read_JASPAR(){

	cout << "- [2] Reading JASPAR MATRIX file and extracting values\n";

	ifstream file(JASPAR_FILE);
	string line;		

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

	//Cout of step 3/4 here because the normalization and reverse function will be re-utilized during the workflow
	cout << "- [3] Jaspar Matrix normalization\n";
	cout << "- [4] Jaspar Matrix reverse complement determination to analize the reverse strand\n";

	return matrix;
}

//Function which saves into a vector called col_sum all the score column sums --> This is made to perform the next Normalization step faster
vector<double> matrix_class::find_col_sum(vector<vector<double>> matrix){

	vector<double> col_sum;						
	double sum = 0;							

	for (unsigned int i = 0; i < matrix[0].size(); i++) {		
		for (unsigned int j = 0; j < 4; j++){			

			sum = sum + matrix[j][i];			
		}

		col_sum.emplace_back(sum);				
		sum = 0;						
	}

	return col_sum;
}

//Function useful to normalize matrix scores and adding a pseudocount to them
void matrix_class::matrix_normalization_pseudoc(vector<vector<double>> matrix){  

	double normalized_score;

	//Calculate and save the scores column sum into a vector to perform a faster normalization step	
	vector<double> col_sum = find_col_sum(matrix);

	for (unsigned int i = 0; i < matrix.size(); i++) {		

		vector<double> normalized_matrix_line;
		for (unsigned int j = 0; j < matrix[i].size(); j++){

			normalized_score = matrix[i][j]/col_sum[j];
			normalized_matrix_line.emplace_back(normalized_score + pseudoc);
		}

		norm_matrix.emplace_back(normalized_matrix_line);
	}
}

//Function to perform a second normalization on matrix scores (without a pseudocount addition)
void matrix_class::matrix_normalization(vector<vector<double>> matrix){

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
void matrix_class::matrix_logarithmic(vector<vector<double>> matrix){

	for(unsigned int i=0; i < matrix.size(); i++){
		
		vector<double> log_matrix_line;
		double log_scores;

		for(unsigned int j=0; j < norm_matrix[i].size(); j++){

			log_scores = log(norm_matrix[i][j]);
			log_matrix_line.emplace_back(log_scores);
		}

		matrix_log.emplace_back(log_matrix_line);
	}
}

//Function which return the Transposed matrix from a matrix in input
vector<vector<double>> matrix_class::reverse_matrix(vector<vector<double>> matrix){

	vector<vector<double>> rev_matrix = matrix;
	reverse(rev_matrix.begin(), rev_matrix.end());

	for (int i = 0; i < 4; i++) {

		reverse(rev_matrix[i].begin(), rev_matrix[i].end());
	}

	return rev_matrix;

}

//Finding the best and the worst score that an oligo can reach based on current JASPAR matrix
void oligo_class::find_minmax(vector<vector<double>> matrix){

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

//Function to find the best oligo score. From every sequence from GEP the best oligo is calculated and both oligo and position in the window are saved
unsigned int oligo_class::find_best_score(){

	//Extracting the best score from oligo_scores with function max_element
	best_score = *max_element(oligo_scores.begin(), oligo_scores.end());

	vector<int> positions;
	vector<int> dist_center;
	unsigned int matches = 0;
	int min_distance;
	vector<int>::iterator itr;

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

			int distance;

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

//Best score normalization with normalization formula (The parameter to normalize have already been calculated and saved into the class)
void oligo_class::best_score_normalization(){
	
	double score_normalized;
	vector<double> oligo_scores_normalized;

	for(unsigned int score=0; score < oligo_scores.size(); score++){

	score_normalized = 1 + ((oligo_scores[score] - max_possible_score)/(max_possible_score - min_possible_score));
	oligo_scores_normalized.emplace_back(score_normalized);
	}
	
	oligo_scores = oligo_scores_normalized;
}

//The oligo which has generated the best score is extracted from fasta sequence and saved into best_oligo_seq variable
void oligo_class::find_best_sequence(string sequence, unsigned int length){

	best_oligo_seq = sequence.substr(local_position,length);
}

//Best oligo coordinates are saved
void oligo_class::find_coordinate( unsigned int length, string chr_coord_GEP, unsigned int start_coord_GEP){

	chr_coord_oligo = chr_coord_GEP;
	start_coord_oligo = start_coord_GEP + local_position;
	end_coord_oligo = start_coord_oligo + length;

}
	
//Function to re-set the genomic coordinates and the sequences window --> centered on the best oligo found for each sequence
void coordinator_class::centering_oligo(){

	TwoBit * tb;
	tb = twobit_open(TWOBIT_FILE.c_str());
	int center_oligo ;
	
	//To center on the best oligo, centering_function and extract_seq functions from bed_class need to be recalled with updated input parameters
	for(unsigned int i=0; i<oligos_vector.size(); i++){

		//The center of the window is exactly on the center of the best oligo (which length depends to the JASPAR matrix size)
		center_oligo = oligos_vector[i].return_start_coord_oligo() + matrix_log[0].size()/2;
		GEP[i].centering_function(center_oligo,center_oligo,half_length,0);
		GEP[i].extract_seq(tb,0);
	}
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

//Function able to convert a string (containing numbers separated by ",") into a vector of unsigned int
vector<unsigned int> map_class::generic_vector_creation(string numbers){

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

//Function to create a map of oligos + their occurrences along all the sequences
void map_class::table_creation_orizzontal(vector<bed_class> GEP){ 


	if (MFASTA_FILE.size() ==0) 
		cout << "- [7] Counting all k-mers occurrences for sequence and positions  \n";
	else
		cout << "- [4] Counting all k-mers occurrences for sequence and positions  \n";
	
	//A map is created for each k-mer inserted as input
	for(unsigned int k=0; k<kmers_vector.size(); k++){
		
		//For every sequence contained into GEP vector
		for(unsigned int j=0; j<GEP.size(); j++){
			
			//Extract the FASTA sequence from each bed class in GEP
			string sequence = GEP[j].return_sequence(GEP[j]);
			
			//Extracted and analyzed all words of length k that are found by scrolling through the sequence
			for(unsigned int i=0; i < (sequence.size() - kmers_vector[k] + 1); i++){
				
				//The current k-length oligo is saved into bases string
				string bases = sequence.substr(i,kmers_vector[k]);

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
}

//Function to create maps to count oligos occurrences for each sequence position (in this function also frequences are calculated)
void map_class::table_creation_vertical(vector<bed_class> GEP){
	
	//Return the fisrt sequence to know the sequences length
	string seq_length = GEP[0].return_sequence(GEP[0]);

	//A vector of map is created for each k-mer inserted as input
	for(unsigned int k=0; k<kmers_vector.size(); k++){

		//vector to store the total number of possible oligo per position, useful to calculate frequences
		vector<unsigned int> tot_freq_vec;


		//Extracted and analyzed the oligo in position "i"
		for(unsigned int i=0; i < (seq_length.size() - kmers_vector[k] + 1); i++){

			unsigned int tot_freq = 0;
			
			//Make the analysis of all the sequences' oligo in position "i" (vertical analisys)
			for(unsigned int j=0; j<GEP.size(); j++){
				
				string sequence = GEP[j].return_sequence(GEP[j]);
				string bases = sequence.substr(i,kmers_vector[k]);
				
				//Calling of function to count oligo occurrnces and to create and fill the maps
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
	}
}

//Function to fill orizzontal plus/minus maps --> current oligo "bases" is passed as parameter to be inserted into the maps
void map_class::or_ver_kmer_count(string bases,unordered_map<string,unsigned int> &plus, unordered_map<string,unsigned int> &minus){

	unordered_map<string,unsigned int>::iterator it_plus;
	unordered_map<string,unsigned int>::iterator it_minus;
	
	//Finding the oligo "bases" into the orizzontal plus map (FWD strand)
	it_plus = plus.find(bases);

	//Check if the current oligo "bases" is palindrome --> this function has also the aim to generate "reverse bases", the reverse complement of the current oligo
	check_palindrome(bases);

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
	reverse_bases.clear();
}

void map_class::vertical_kmer_count(string bases,map<pair<string,string>,pair<unsigned int, unsigned int>>&plus, unsigned int& tot_freq){

	map<pair<string,string>,pair<unsigned int, unsigned int>>::iterator it_plus;
	map<pair<string,string>,pair<unsigned int, unsigned int>>::iterator it_plus_rev;
	map<pair<string,string>,pair<unsigned int, unsigned int>>::iterator it_minus;
	map<pair<string,string>,pair<unsigned int, unsigned int>>::iterator it_minus_rev;

	//Check if the current oligo is palindrome
	bool pal = check_palindrome(bases);

	//Check if the analisys is in double strand
	if(DS == 1){	

		//If the current oligo is not palindrome and the analysis is on DS --> total number of possible oligo per position (to calculate frequences) ++2, else if it is palindrome ++1
		if(!pal){
			tot_freq = tot_freq+2;
		}
		else{
			tot_freq++;
		}
	}

	//If analysis is on SS --> possible oligos per position ++1 (it will be equal to the number of sequences)
	else{
		tot_freq++;
	}
	
	//Creation of a pair of strings containing current oligo (first) and his reverse complement (second)
	pair<string,string> pair_bases;
	pair_bases.first= bases;
	pair_bases.second= reverse_bases;
	
	//Creation of a pair of strings containing reverse complement (first) and current oligo (second)
	pair<string,string> pair_bases_reverse;
	pair_bases_reverse.first= reverse_bases;
	pair_bases_reverse.second= bases;
	
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

void p_value_class::N2_calculation(unordered_map<string,unsigned int> orizzontal_map){

	total_oligo_N2 = 0;
	total_oligo_N2 = accumulate(begin(orizzontal_map), end(orizzontal_map), 0, [] (unsigned int val, const unordered_map<string,int>::value_type& p) {return val + p.second;});
}

//Function to check, given an oligo as input, if this oligo is palindrome or not
bool map_class::check_palindrome(string bases){

	//For any character of the string insert into another string (called reverse bases) the complementary character
	for(unsigned int i=0; i<bases.size(); i++){

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
	
	//Then reverse the string 
	reverse(reverse_bases.begin(), reverse_bases.end());
	
	//If they are equal --> it means that the oligo "bases" is palindrome
	if (reverse_bases == bases){
		return true;
	}
	else {return false;}

}

bool p_value_class::check_palindrome2(string bases){

	reverse_bases.clear();
	for(unsigned int i=0; i<bases.size(); i++){

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

	reverse(reverse_bases.begin(), reverse_bases.end());
	if (reverse_bases == bases){
		return true;
	}
	else {return false;}

}


void multifasta_class::length_control(vector<string> sequences){

	cout << "- [2] Multifasta Sequences length check\n";

	unsigned int size = sequences[0].size();

	for(unsigned int i=0; i<sequences.size(); i++){
		
		//If only one sequence in the vector are longer an error is generated
		if(sequences[i].size() != size){

			cerr << "Sequences are not of the same length!" << endl;
			exit(1);
		}
	}
}

//Function to extract fasta sequences from a Multifasta file and save them into a vector of string
void multifasta_class::extract_sequences(){

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

void multifasta_class::GEP_creation_MF(vector<string> sequences){

	cout << "- [3] Sorting Multifasta sequences\n";

	for(unsigned int i=0; i<sequences.size(); i++){
		
		//A bed class is created for each sequence
		bed_class BED_MULTIFASTA(sequences[i]);

		//The classes are stored into a GEP vector
		GEP.emplace_back(BED_MULTIFASTA);
	}

}

void map_class::P_VALUE_MATRIX_creation(){

	for(unsigned int j=0; j<vector_kmers_maps_plus.size(); j++){

		ofstream outfile = outfile_header(j);

		for(unsigned int i=0; i<vector_kmers_maps_plus[j].size(); i++){

			p_value_class P(vector_kmers_maps_plus[j][i], orizzontal_plus_debug[j], sequences_number_T, i, outfile, tot_freq_matrix[j][i]);
			P_VALUE_VECTOR.emplace_back(P); //creating a vector for every position
			tot_sum_vector.emplace_back(P.return_sum_top_N());
		}
		P_VALUE_MATRIX.emplace_back(P_VALUE_VECTOR); //creating a vector for every k
		P_VALUE_VECTOR.clear();
		tot_sum_matrix.emplace_back(tot_sum_vector);
		tot_sum_vector.clear();
		outfile.close();

	}
}


multimap<pair<unsigned int, unsigned int>,pair<string,string>> p_value_class::multimap_creation(map<pair<string,string>, pair<unsigned int,unsigned int>> pair_map){

	for(map<pair<string,string>,pair<unsigned int, unsigned int>>::iterator it = pair_map.begin(); it != pair_map.end(); it++){ 	

		vertical_multimap.insert({it->second,it->first});
	}

	return vertical_multimap;
}

void p_value_class::filling_KNT_vectors(unordered_map<string,unsigned int> orizzontal_map){

	unsigned int K;
	unsigned int N1;
	unsigned int N2;

	for(multimap<pair<unsigned int, unsigned int>,pair<string,string>>::reverse_iterator it_rev = vertical_multimap.rbegin(); it_rev != vertical_multimap.rend(); it_rev++){

		bool pal = check_palindrome2(it_rev->second.first);

		if(pal == 0 && DS == 1){

			K = it_rev->first.first + it_rev->first.second;
		}
		else{

			K = it_rev->first.first;
		}
		it_N1_plus = orizzontal_map.find(it_rev->second.first);
		it_N1_minus = orizzontal_map.find(it_rev->second.second);

		if(it_N1_minus != orizzontal_map.end() && DS == 1){

			N1 = it_N1_plus->second + it_N1_minus->second;
		}
		else{
			N1 = it_N1_plus->second;
		}

		N2 = total_oligo_N2 - N1;

		K_vec.emplace_back(K);
		N1_vec.emplace_back(N1);
		N2_vec.emplace_back(N2);

	}
}

void p_value_class::calculating_p_value(){

	for(unsigned int i=0; i<K_vec.size(); i++){

		double p_value =  gsl_cdf_hypergeometric_Q(K_vec[i],N1_vec[i],N2_vec[i],T);
		p_value = check_p_value(p_value);
		p_value_vec.emplace_back(p_value);	
	}
}

double p_value_class::check_p_value(double p){

	if(p == 0){

		p = 1.000001e-300;
	}

	return p;
}

void p_value_class::sorting_p_value(){

	vector<unsigned int> KNT;
	unsigned int i=0;
	for(multimap<pair<unsigned int, unsigned int>, pair<string, string>>::reverse_iterator it_rev = vertical_multimap.rbegin(); it_rev!=vertical_multimap.rend(); it_rev++){

		KNT.emplace_back(K_vec[i]);	
		KNT.emplace_back(N1_vec[i]);	
		KNT.emplace_back(N2_vec[i]);	
		KNT.emplace_back(T);	

		p_value_sort.insert({p_value_vec[i],{it_rev->second.first, it_rev->second.second}});
		p_value_KNT.insert({p_value_vec[i],{it_rev->second.first, KNT}}); 
		i = i+1;
		KNT.clear();
	}
}

void p_value_class::checking_ordering(map<pair<string,string>,pair<unsigned int,unsigned int>> pair_map, unsigned int position, ofstream &outfile, unsigned int freq){

	if(ordering == "p"){

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

void p_value_class::print_debug_p_value_DS(map<pair<string,string>,pair<unsigned int,unsigned int>> pair_map, unsigned int position, ofstream& outfile, unsigned int freq){

	unsigned int c=0;
	sum_top_N = 0;

	multimap<pair<string,string>,pair<unsigned int, unsigned int>>::iterator it_multi;

	for(multimap<double,pair<string,string>>::iterator it_pair = p_value_sort.begin(); it_pair!=p_value_sort.end() && c<top_N; it_pair++, c++){

		double FREQ, Sum_Occ_Oligo;
		it_multi = pair_map.find(it_pair->second);
		string Oligo = it_multi->first.first;
		string Oligo_RC = it_multi->first.second;
		unsigned int Num_Occ_FWD = it_multi->second.first;
		unsigned int Num_Occ_REV = it_multi->second.second;
		bool pal = check_palindrome2(Oligo);
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

void p_value_class::print_debug_p_value_SS(map<pair<string,string>,pair<unsigned int,unsigned int>> pair_map, unsigned int position, ofstream& outfile, unsigned int freq){

	unsigned int c = 0;
	sum_top_N = 0;

	multimap<pair<string,string>,pair<unsigned int, unsigned int>>::iterator it_multi;

	for(multimap<double,pair<string,string>>::iterator it_pair = p_value_sort.begin(); it_pair!=p_value_sort.end() && c<top_N; it_pair++, c++){

		double FREQ, Num_Occ_FWD;
		it_multi = pair_map.find(it_pair->second);
		string Oligo = it_multi->first.first;
		unsigned int Num_Occ_Oligo = it_multi->second.first;
		bool pal = check_palindrome2(Oligo);
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

void p_value_class::print_debug_occurrences_DS(map<pair<string,string>,pair<unsigned int,unsigned int>> pair_map, unsigned int position, ofstream& outfile, unsigned int freq, vector<double> p_value_vec){

	unsigned int c=0;
	sum_top_N = 0;

	for(multimap<pair<unsigned int,unsigned int>,pair<string,string>>::reverse_iterator it_rev = vertical_multimap.rbegin(); it_rev!=vertical_multimap.rend() && c<top_N; it_rev++, c++){

		double FREQ, Sum_Occ_Oligo;
		string Oligo = it_rev->second.first;
		string Oligo_RC = it_rev->second.second;
		unsigned int Num_Occ_FWD = it_rev->first.first;
		unsigned int Num_Occ_REV = it_rev->first.second;
		bool pal = check_palindrome2(Oligo);
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

void p_value_class::print_debug_occurrences_SS(map<pair<string,string>,pair<unsigned int,unsigned int>> pair_map, unsigned int position, ofstream& outfile, unsigned int freq, vector<double> p_value_vec){

	unsigned int c = 0;
	sum_top_N = 0;


	for(multimap<pair<unsigned int,unsigned int>,pair<string,string>>::reverse_iterator it_rev = vertical_multimap.rbegin(); it_rev!=vertical_multimap.rend() && c<top_N; it_rev++, c++){

		double FREQ, Num_Occ_FWD;
		string Oligo = it_rev->second.first;
		unsigned int Num_Occ_Oligo = it_rev->first.first;
		bool pal = check_palindrome2(Oligo);
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

void map_class::HAMMING_MATRIX_creation(vector<bed_class> GEP){
	
	for(unsigned int j=0; j<P_VALUE_MATRIX.size(); j++){

		ofstream outfile = outfile_header_hamming(j);
		
		for(unsigned int i=0; i<P_VALUE_MATRIX[j].size(); i++){

			multimap<pair<unsigned int, unsigned int>, pair<string,string>> vertical_multimap = P_VALUE_MATRIX[j][i].return_vertical_multimap();

			hamming_class H(vertical_multimap,distance_vector[j],i,tot_freq_matrix[j][i],orizzontal_plus_debug[j], orizzontal_minus_debug[j], outfile, GEP);
			HAMMING_VECTOR.emplace_back(H);
		}

		HAMMING_MATRIX.emplace_back(HAMMING_VECTOR);
	//	Z_TEST_MATRIX.emplace_back(Z_TEST_VECTOR);
		HAMMING_VECTOR.clear();
	//	Z_TEST_VECTOR.clear();
		outfile.close();
	}
}

void hamming_class::find_best_oligos(){

	multimap<pair<unsigned int,unsigned int>, pair<string,string>>::reverse_iterator it_rev = vertical_multimap.rbegin();
	real_best_oligo_occurrences = (it_rev->first.first + it_rev->first.second);
	unsigned int verical_size = vertical_multimap.size();
	bool flag = 1;
	unsigned int counter = 1;
	
	//If all the sequences have the same oligo (100%) in a specific position --> vertical_multimap.size() == 1 --> This control is made to avoid an infinite while cycling
	while(it_rev->first.first + it_rev->first.second == real_best_oligo_occurrences && flag == 1){
		
		if(counter == verical_size){
			flag = 0;
		}
		best_oligos.emplace_back(it_rev->second.first);
		it_rev++;
		counter++;
	}
}

void hamming_class::checking_best_oligo(unsigned int distance){
	
	if(best_oligos.size() == 1){
		
		real_best_oligo = best_oligos[0];
		find_distanced_oligos(real_best_oligo,distance);	
	}

	else{

		real_best_oligo = select_real_best_oligo(distance);
		find_distanced_oligos(real_best_oligo,distance);	
	}

}

string hamming_class::select_real_best_oligo(unsigned int distance){
	
	unsigned int max_similarity;
	unsigned int index;

	for(unsigned int i=0; i<best_oligos.size(); i++){

		find_distanced_oligos(best_oligos[i], distance);
		
		if(i==0){
			max_similarity = similar_oligos.size();
			index = 0;
		}	

		if(similar_oligos.size() > max_similarity){

			max_similarity = similar_oligos.size();
			index = i;
		}

		similar_oligos.clear(); 
		similar_oligos_occurrences.clear();  

	}
	
	return best_oligos[index];
}

void hamming_class::find_distanced_oligos(string best, unsigned int distance){

	bool is_similar;

	for(multimap<pair<unsigned int,unsigned int>, pair<string,string>>::reverse_iterator it_rev = vertical_multimap.rbegin(); it_rev != vertical_multimap.rend(); it_rev++){

		if(it_rev->second.first != best){

			is_similar = is_similar_oligo(best, it_rev->second.first, distance);

			if(is_similar == 1){

				similar_oligos.emplace_back(it_rev->second.first);
				similar_oligos_occurrences.emplace_back(it_rev->first.first);
			}

			if(DS==1){

				is_similar = is_similar_oligo(best, it_rev->second.second, distance);

				if(is_similar == 1){

					similar_oligos.emplace_back(it_rev->second.second);
					similar_oligos_occurrences.emplace_back(it_rev->first.second);
				}
			}
		}	

	}
}

bool hamming_class::is_similar_oligo(string oligo_1, string oligo_2, unsigned int distance){
	
	unsigned int counter = 0;

	for(unsigned int i = 0; i<oligo_1.size() && counter <= distance; i++){

		if(oligo_1[i] != oligo_2[i]){
			
			counter++;	
		}
	}
	
	return(counter<=distance);
}

double hamming_class::frquence_1_calculation(unsigned int freq){

	tot_similar_occurrences = real_best_oligo_occurrences;
	
	for(unsigned int i=0; i<similar_oligos_occurrences.size(); i++){

		tot_similar_occurrences = tot_similar_occurrences + similar_oligos_occurrences[i];

	}

	double FREQ_1 = tot_similar_occurrences/freq;

	return FREQ_1;
}

double hamming_class::frquence_2_calculation(unordered_map<string,unsigned int> orizzontal_map_plus, unordered_map<string,unsigned int> orizzontal_map_minus, unsigned int position){
	

	unsigned int total_orizzontal_occurrences = finding_orizzontal_occurrences(orizzontal_map_plus, orizzontal_map_minus);
	double FREQ_2 = tot_similar_occurrences/total_orizzontal_occurrences;

	return FREQ_2;
}

unsigned int hamming_class::finding_orizzontal_occurrences(unordered_map<string,unsigned int> orizzontal_map_plus, unordered_map<string,unsigned int> orizzontal_map_minus){

	unordered_map<string,unsigned int>::iterator it = orizzontal_map_plus.find(real_best_oligo);
	unsigned int total_orizz_occ = it->second;
	
	for(unsigned int i=0; i<similar_oligos.size(); i++){
		
		it = orizzontal_map_plus.find(similar_oligos[i]);
		if(it == orizzontal_map_plus.end()){
			
			it = orizzontal_map_minus.find(similar_oligos[i]);
		}
		total_orizz_occ = (total_orizz_occ + it->second);
	}

	return total_orizz_occ;
}

void hamming_class::PWM_hamming_creation(){

	double counter_A = 0;
	double counter_C = 0;
	double counter_G = 0;
	double counter_T = 0;
	vector<double> vec_A, vec_C, vec_G, vec_T;
	
	for(unsigned int character = 0; character < similar_oligos[0].size(); character++){

		counter_A = 0;
		counter_C = 0;
		counter_G = 0;
		counter_T = 0;

		for(unsigned int oligo = 0; oligo < similar_oligos.size(); oligo++){

			switch(similar_oligos[oligo][character]){


				case 'A' : counter_A = counter_A + similar_oligos_occurrences[oligo]; 
					   break;
				case 'C' : counter_C = counter_C + similar_oligos_occurrences[oligo]; 
					   break;
				case 'G' : counter_G = counter_G + similar_oligos_occurrences[oligo];
					   break;
				case 'T' :  counter_T = counter_T + similar_oligos_occurrences[oligo]; 
					   break;
			}
		}

		vec_A.emplace_back(counter_A);
		vec_C.emplace_back(counter_C);
		vec_G.emplace_back(counter_G);
		vec_T.emplace_back(counter_T);
	}
	
	PWM_hamming.emplace_back(vec_A);
	PWM_hamming.emplace_back(vec_C);
	PWM_hamming.emplace_back(vec_G);
	PWM_hamming.emplace_back(vec_T);
}

void map_class::Z_TEST_MATRIX_creation(vector<bed_class> GEP){

	for(unsigned int i=0; i<HAMMING_MATRIX.size(); i++){
		for (unsigned int j=0; j<HAMMING_MATRIX[i].size(); j++){

			vector<vector<double>> PWM_matrix = HAMMING_MATRIX[i][j].return_PWM_hamming();
			double Frequence_1 = HAMMING_MATRIX[i][j].return_FREQUENCE_1();

			if(Frequence_1 >= freq_treshold){

				z_test_class Z(PWM_matrix, GEP,j+1,kmers_vector,HAMMING_MATRIX);
				Z_TEST_VECTOR.emplace_back(Z);
			}
		}

		Z_TEST_MATRIX.emplace_back(Z_TEST_VECTOR);
		Z_TEST_VECTOR.clear();
	}
}

void z_test_class::oligos_vector_creation_PWM(vector<bed_class> GEP){

	for(unsigned int i=0; i<GEP.size(); i++){
	
		string sequence = GEP[i].return_sequence(GEP[i]);
		oligo_class SHIFTING_PWM(matrix_log, sequence);
		oligo_scores_orizzontal_FWD = SHIFTING_PWM.return_oligo_scores();
	
		if(DS == 1){

			oligo_class SHIFTING_PWM_2(inverse_matrix_log, sequence);
			oligo_scores_orizzontal_REV = SHIFTING_PWM_2.return_oligo_scores();
			check_best_strand_oligo();			
			all_local_scores.emplace_back(oligo_scores_orizzontal_BEST[local_pos-1]);	
			all_global_scores.insert(all_global_scores.end(), oligo_scores_orizzontal_BEST.begin(), oligo_scores_orizzontal_BEST.end());
		}

		else{
			
			all_local_scores.emplace_back(oligo_scores_orizzontal_FWD[local_pos-1]);	
			all_global_scores.insert(all_global_scores.end(), oligo_scores_orizzontal_FWD.begin(), oligo_scores_orizzontal_FWD.end());
		}
	
		oligo_scores_orizzontal_BEST.clear();	
	}	
	
}

void z_test_class::check_best_strand_oligo(){

	for(unsigned int oligo = 0; oligo < oligo_scores_orizzontal_FWD.size(); oligo ++){

		if(oligo_scores_orizzontal_FWD[oligo] >= oligo_scores_orizzontal_REV[oligo]){

			oligo_scores_orizzontal_BEST.emplace_back(oligo_scores_orizzontal_FWD[oligo]);
		}

		else{

			oligo_scores_orizzontal_BEST.emplace_back(oligo_scores_orizzontal_REV[oligo]);
		}
	}
}

void z_test_class::global_mean_calculation(){

	double local_sum = accumulate(all_local_scores.begin(), all_local_scores.end(),0.0);	
	double global_sum = accumulate(all_global_scores.begin(), all_global_scores.end(), 0.0);
	double tot_sq_sum_global = inner_product(all_global_scores.begin(), all_global_scores.end(), all_global_scores.begin(), 0.0);
	double tot_sq_sum_local = inner_product(all_local_scores.begin(), all_local_scores.end(), all_local_scores.begin(), 0.0);
	global_mean = global_sum/all_global_scores.size();
	local_mean = local_sum/all_local_scores.size();
	global_dev_std = sqrt(tot_sq_sum_global/all_global_scores.size() - global_mean * global_mean);
	local_dev_std = sqrt(tot_sq_sum_local/all_local_scores.size() - local_mean * local_mean);
	
}
void z_test_class::z_score_calculation(){


	z_score = ((global_mean - local_mean)/ (local_dev_std / sqrt(all_local_scores.size()))); 

	const double Z  = z_score;
	Zpvalue = gsl_cdf_ugaussian_P(Z);


}

/////DEBUG/////////////////////////////////////////////////////////

unsigned int oligo_class::return_start_coord_oligo(){

	return start_coord_oligo;
}

double oligo_class::return_best_score_normalized(){

	return best_score_normalized;
}

string bed_class::return_chr_coord_GEP(){

	return chr_coord;
}

unsigned int bed_class::return_start_coord_GEP(){

	return start_coord;
}

vector<vector<double>> matrix_class::return_inverse_log_matrix(){

	return inverse_matrix_log;
}

vector<vector<double>> matrix_class::return_log_matrix(){

	return matrix_log;
}

string bed_class::return_sequence(bed_class){

	return sequence;
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

multimap<double,pair<string,vector<unsigned int>>> p_value_class::return_p_value_KNT(){

	return p_value_KNT;
}

multimap<pair<unsigned int,unsigned int>, pair<string,string>> p_value_class::return_vertical_multimap(){

	return vertical_multimap;
}

unsigned int p_value_class::return_sum_top_N(){

	return sum_top_N;
}

vector<unsigned int> p_value_class::return_K_vec(){

	return K_vec;
}

vector<unsigned int> p_value_class::return_N1_vec(){

	return N1_vec;
}

vector<unsigned int> p_value_class::return_N2_vec(){

	return N2_vec;
}

unsigned int p_value_class::return_T(){

	return T;
}

vector<double> p_value_class::return_p_value_vec(){

	return p_value_vec;
}

string hamming_class::return_real_best_oligo(){

	return real_best_oligo;
}

unsigned int hamming_class::return_similar_oligo_size(){

	return similar_oligos.size();
}

vector<vector<double>> hamming_class::return_PWM_hamming(){

	return PWM_hamming;
}

vector<double> oligo_class::return_oligo_scores(){

	return oligo_scores;
}

double hamming_class::return_FREQUENCE_1(){

	return FREQUENCE_1;
}

unsigned int z_test_class::return_local_pos(){

	return local_pos;
}

double z_test_class::return_local_mean(){

	return local_mean;
}

double z_test_class::return_global_mean(){

	return global_mean;
}

double z_test_class::return_local_std_dev(){

	return local_dev_std;
}

double z_test_class::return_global_std_dev(){

	return global_dev_std;
}
double z_test_class::return_Zpvalue(){

	return Zpvalue;
}
double z_test_class::return_z_score(){

	return z_score;
}

//Function to matrix debugging --> it prints the scores extracted from JASPAR file, the normalized scores, the logarithmic scores and the logarithmic score of transposed matrix
void matrix_class::debug_matrix(matrix_class M){

	M.print_debug_matrix(matrix, " ");
	M.print_debug_matrix(norm_matrix, " NORMALIZED");
	M.print_debug_matrix(matrix_log, " LOGARITHMIC MATRIX");
	M.print_debug_matrix(inverse_matrix_log, " INVERSE LOGARITHMIC MATRIX");
}

//Function to matrix debugging --> it prints the scores extracted from JASPAR file, the matrix name and the tf name
void matrix_class::print_debug_matrix(vector<vector<double>> matrix, string type){

	cout << "\n" << matrix_name << " " << tf_name << type << ":" << endl;

	for(unsigned int i=0; i < matrix.size(); i++){
		for(unsigned int j=0; j<matrix[i].size(); j++){

			cout << matrix[i][j] << " ";
		}

		cout << endl;
	}
}

//Debug function: Print sequences and coordinates from GEP vector into a .fasta file to check if the sequences extraction is correct
void coordinator_class::print_debug_GEP(vector<bed_class> GEP){

	//Twobit_JASPAR_Bed used to create GEP vector saved into alias file to name the outputs	
	alias_file = (TWOBIT_FILE.erase(0,TWOBIT_FILE.find_last_of("/")+1)+"_"+ JASPAR_FILE.erase(0,JASPAR_FILE.find_last_of("/")+1)+"_"+ BED_FILE.erase(0,BED_FILE.find_last_of("/")+1));

	//Output file .bed carrying the centered coordinates
	ofstream outfile;	
	JASPAR_FILE = JASPAR_FILE.erase(JASPAR_FILE.find_last_of("."), JASPAR_FILE.size());
	outfile.open(alias_file);

	for(unsigned int i=0; i<GEP.size(); i++){

		string chr_coord = GEP[i].return_chr_coord();
		unsigned int start_coord = GEP[i].return_start_coord();
		unsigned int end_coord = GEP[i].return_end_coord();
		outfile << chr_coord << "\t" << start_coord << "\t" << end_coord << endl;
	}

	outfile.close();

	//Output file .fasta carrying the centered coordinates and the sequences extracted
	BED_FILE = BED_FILE.erase(BED_FILE.find_last_of("."), BED_FILE.size());
	outfile.open(alias_file+".fasta");

	for(unsigned int i=0; i<GEP.size(); i++){

		string chr_coord = GEP[i].return_chr_coord();
		unsigned int start_coord = GEP[i].return_start_coord();
		unsigned int end_coord = GEP[i].return_end_coord();
		string sequence = GEP[i].return_sequence(GEP[i]);				

		outfile << ">" << chr_coord << ":" << start_coord << "-" << end_coord << endl;	
		outfile << sequence << endl;
	}

	outfile.close();
}

void oligo_class::oligos_vector_debug(){	//Debug function to print the best oligo features

	cout << endl;
	cout << "Sequence: " << global_sequence << endl;
	cout << "The hit position is " << local_position << endl;
	cout << "The genomic coordinates are:\n> " << chr_coord_oligo << ": " << start_coord_oligo << " - " << end_coord_oligo << endl;
	cout << "The best score is " << best_score << endl;
	cout << "The best score normalized is " << best_score_normalized << endl;
	cout << "The best oligo sequence is " << best_oligo_seq << endl;
	cout << "Strand  " << strand << endl;
	cout << endl;
}

void map_class::print_debug_orizzontal(){

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
				bool palindrome = check_palindrome(it_rev->second);

				if(!palindrome){

					unordered_map<string,unsigned int>::iterator find_RC = orizzontal_minus_debug[i].find(reverse_bases);
					outfile << it_rev->second << "\t" << it_rev->first << "\t" << find_RC->first << "\t" << find_RC->second << "\t" << endl;

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

ofstream map_class::outfile_header(unsigned int j){


	ofstream outfile;

	if(DS==1){
		
		if(ordering == "p"){
		
			outfile.open(to_string(kmers_vector[j])+"-mers_positional_occurrences_"+alias_file+"DS_p_val.txt");
		}
		else{
		
			outfile.open(to_string(kmers_vector[j])+"-mers_positional_occurrences_"+alias_file+"DS_occ.txt");
		}
		
		outfile << "#Maps vector with kmers occurences (Double Strand) counted for positions in sequence (for k = " << kmers_vector[j] << "):" << endl;
		outfile << "#Position" << "\t" << "Rank" << "\t" << "Oligo" << "\t" << "Num_Occ_FWD" << "\t" << "Num_Occ_REV" << "\t" << "Sum_Occ_Oligo" << "\t" << "Oligo_RC" << "\t" << "Num_Occ_RC_FWD" << "\t" << "Num_Occ_RC_REV" << "\t" << "Sum_Occ_RC" << "\t" << "PAL" << "\t" << "Tot_Occ" << "\t" << "FREQ" << "\t" << "P_VALUE" << endl;

	}

	else{
		
		if(ordering == "p"){
		
			outfile.open(to_string(kmers_vector[j])+"-mers_positional_occurrences_"+alias_file+"SS_p_val.txt");
		}
		else{
		
			outfile.open(to_string(kmers_vector[j])+"-mers_positional_occurrences_"+alias_file+"SS_occ.txt");
		}

		outfile << "#Maps vector with kmers occurences (Single Strand) counted for positions in sequence (for k = " << kmers_vector[j] << "):" << endl;
		outfile << "#Position" << "\t" << "Rank" << "\t" << "Oligo" << "\t" << "Num_Occ_Oligo" << "\t" << "PAL" << "\t" << "FREQ" << "\t" << "P_VALUE" << endl;

	}
	return outfile;	
}


void map_class::TopN_sum_and_freq(){

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

		outfile << "###Top " << top_N << " occurrences sum with k = " << kmers_vector[i] << ":" << endl; 
		outfile << "Position" << "\t" << "Sum" << "\t" << "Frequences" << endl; 

		for(unsigned int j=0; j<tot_sum_matrix[i].size(); j++){

			double sum = tot_sum_matrix[i][j];		//Put as double to dont loose the precision	
			double frequence = sum/tot_freq_matrix[i][j];	
			outfile << j+1 << "\t" << tot_sum_matrix[i][j] << "\t" << frequence << endl; 

		}
		outfile.close();
	}
}

void map_class::p_value_parameters_debug_p_val(){

	ofstream outfile;

	for(unsigned int j = 0; j<P_VALUE_MATRIX.size(); j++){
		
		if(DS == 1){
		outfile.open(to_string(kmers_vector[j])+"-mers_p_value_parameters_control_p_val_"+BED_FILE+"DS.txt");
		}
		else{
		outfile.open(to_string(kmers_vector[j])+"-mers_p_value_parameters_control_p_val_"+BED_FILE+"SS.txt");
		}

		outfile << "#Parameters used to calculate p_value for each oligo positionally ranked" << endl;
		outfile << "#Position" << "\t" << "Rank" << "\t" << "Oligo" << "\t" << "K" << "\t" << "N1" << "\t" << "N2" << "\t" << "T" << "\t" << "P_VALUE" << "\t" <<"P_VALUE_LOG10" << endl;

		for(unsigned int i = 0; i<P_VALUE_MATRIX[j].size(); i++){

			multimap<double,pair<string,vector<unsigned int>>> KNT_multimap = P_VALUE_MATRIX[j][i].return_p_value_KNT();
			unsigned int rank = 0;

			for(multimap<double,pair<string,vector<unsigned int>>>::iterator it = KNT_multimap.begin(); it != KNT_multimap.end() && rank<top_N; it++, rank++){

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

	ofstream outfile;

	for(unsigned int j = 0; j<P_VALUE_MATRIX.size(); j++){
		
		if(DS == 1){
		outfile.open(to_string(kmers_vector[j])+"-mers_p_value_parameters_control_occ_"+BED_FILE+"DS.txt");
		}
		else{
		outfile.open(to_string(kmers_vector[j])+"-mers_p_value_parameters_control_occ_"+BED_FILE+"SS.txt");
		}

		outfile << "#Parameters used to calculate p_value for each oligo positionally ranked" << endl;
		outfile << "#Position" << "\t" << "Rank" << "\t" << "Oligo" << "\t" << "K" << "\t" << "N1" << "\t" << "N2" << "\t" << "T" << "\t" << "P_VALUE" << "\t" <<"P_VALUE_LOG10" << endl;

		for(unsigned int i = 0; i<P_VALUE_MATRIX[j].size(); i++){

			multimap<pair<unsigned int,unsigned int>, pair<string,string>> vertical_multimap = P_VALUE_MATRIX[j][i].return_vertical_multimap();
			unsigned int rank = 0;
			vector<unsigned int> K_vec = P_VALUE_MATRIX[j][i].return_K_vec();
			vector<unsigned int> N1_vec = P_VALUE_MATRIX[j][i].return_N1_vec();
			vector<unsigned int> N2_vec = P_VALUE_MATRIX[j][i].return_N2_vec();
			unsigned int T = P_VALUE_MATRIX[j][i].return_T();
			vector<double> p_value_vec = P_VALUE_MATRIX[j][i].return_p_value_vec();

			for(multimap<pair<unsigned int,unsigned int>,pair<string,string>>::reverse_iterator it_rev = vertical_multimap.rbegin(); it_rev != vertical_multimap.rend() && rank<top_N; it_rev++, rank++){

				outfile << i+1 << "\t";
				outfile << rank+1 << "\t";
				outfile << it_rev->second.first << "\t";
				outfile << K_vec[rank] << "\t";
				outfile << N1_vec[rank] << "\t";
				outfile << N2_vec[rank] << "\t";
				outfile << T << "\t";
				outfile << p_value_vec[rank] << "\t";
				outfile << log10(p_value_vec[rank])*-1 << endl;
			}

		}
		outfile.close();
	}
}

ofstream map_class::outfile_header_hamming(unsigned int j){


	ofstream outfile;

	if(DS==1){
		
		outfile.open(to_string(kmers_vector[j])+"-mers_hamming_"+alias_file+"DS.txt");
		
		outfile << "#For each best oligo in positions, a table representing hamming distance oligos and 2 different frequences(FREQUENCE_1 = N.occurrences/N.sequences | FREQUENCE_2 = N.occurrences/Orizzontal_occurrences) (for k = " << kmers_vector[j] << "):" << endl;
		outfile << "#Position" << "\t" << "Best_oligo" << "\t" << "Num_Occ_Best_oligo" << "\t" << "Humming_found_oligo_number" << "\t" << "Total_occurrences(best+hamming)" << "\t" << "FREQUENCE_1" << "\t" << "FREQUENCE_2" << endl;

	}
	else{
		outfile.open(to_string(kmers_vector[j])+"-mers_hamming_"+alias_file+"SS.txt");
		
		outfile << "#For each best oligo in positions, a table representing hamming distance oligos and 2 different frequences(FREQUENCE_1 = N.occurrences/N.sequences | FREQUENCE_2 = N.occurrences/Orizzontal_occurrences) (for k = " << kmers_vector[j] << "):" << endl;
		outfile << "#Position" << "\t" << "Best_oligo" << "\t" << "Num_Occ_Best_oligo" << "\t" << "Humming_found_oligo_number" << "\t" << "Total_occurrences(best+hamming)" << "\t" << "FREQUENCE_1" << "\t" << "FREQUENCE_2" << endl;

	}
		
	return outfile;	
}



void hamming_class::print_debug_hamming(unsigned int position, ofstream& outfile){

	outfile << position+1 << "\t" << real_best_oligo << "\t" << real_best_oligo_occurrences << "\t" << similar_oligos.size() << "\t" << tot_similar_occurrences << "\t" << FREQUENCE_1 << "\t" << FREQUENCE_2 << endl;

}

void map_class::Outfile_PWM_hamming(){

	ofstream outfile;

	for(unsigned int j=0; j<kmers_vector.size(); j++){

		if(DS == 1){

			outfile.open(to_string(kmers_vector[j])+"-mers_PWM_hamming_matrices_"+alias_file+"DS.txt");
			print_debug_PWM_hamming(outfile, j, kmers_vector[j]);
			outfile.close();
		}

		else{

			outfile.open(to_string(kmers_vector[j])+"-mers_PWM_hamming_matrices_"+alias_file+"SS.txt");

			print_debug_PWM_hamming(outfile, j, kmers_vector[j]);
			outfile.close();
		}
	}
}

void map_class::print_debug_PWM_hamming(ofstream& outfile, unsigned int j, unsigned int k){
	
	outfile << "#PWM Matrices calculated from the best oligo for each position and his hamming distanced oligos - k = " << k << endl << endl;
	string best_oligo;
	unsigned int neighbour_numb;
	vector<vector<double>> PWM_hamming;
	string ACGT = "ACGT";

	for(unsigned int position = 0; position < Z_TEST_MATRIX[j].size(); position++){

		unsigned int local_pos = Z_TEST_MATRIX[j][position].return_local_pos();
		double global_mean = Z_TEST_MATRIX[j][position].return_global_mean();
		double local_mean = Z_TEST_MATRIX[j][position].return_local_mean();
		double local_dev_std = Z_TEST_MATRIX[j][position].return_local_std_dev();
		double global_dev_std = Z_TEST_MATRIX[j][position].return_global_std_dev();
		double Zpvalue  = Z_TEST_MATRIX[j][position].return_Zpvalue();
		double z_score  = Z_TEST_MATRIX[j][position].return_z_score();

		best_oligo = HAMMING_MATRIX[j][local_pos].return_real_best_oligo();	
		neighbour_numb = HAMMING_MATRIX[j][local_pos].return_similar_oligo_size();
		PWM_hamming = HAMMING_MATRIX[j][local_pos].return_PWM_hamming();
		
		outfile << "#Position " << local_pos << ": \n#PWM calculated from oligo " << best_oligo << " and his " << neighbour_numb << " hamming distanced neighbours. " << endl << endl;

		for(unsigned int i = 0; i< PWM_hamming.size(); i++){
			
			outfile << ACGT[i] << "\t" << "[" << "\t";

			for(unsigned int j = 0; j<PWM_hamming[i].size(); j++){

				outfile << PWM_hamming[i][j] << "\t";
			}
			outfile << "]" << endl;
		}
		
		outfile << endl;
		outfile << "The global mean is: " << global_mean << endl;
		outfile << "The global standard deviation is: " << global_dev_std << endl;
		outfile << "The local mean is: " << local_mean << endl;
		outfile << "The local standard deviation is: " << local_dev_std << endl << endl;
		outfile << "The zscore calculated is: " << z_score<< endl << endl;
		outfile << "The pvalue calculated from the Z score is: " << Zpvalue<< endl << endl;
		outfile << "-------------------------------------------------------------------" << endl;	
	}
}

//void z_test_class::print_debug_oligo_vec(vector<vector<double>> PWM_hamming){
//
//	cout << endl << endl;
//	
//	for(int i=0; i<PWM_hamming.size(); i++){
//		for(int j=0; j<PWM_hamming[0].size(); j++){
//
//			cout << PWM_hamming[i][j] << " ";
//		}
//		cout << endl;
//	}
//	cout << endl;
//
//	for(int i=0; i<matrix_log.size(); i++){
//		for(int j=0; j<matrix_log[0].size(); j++){
//
//			cout << matrix_log[i][j] << " ";
//		}
//		cout << endl;
//	}
//	cout << endl;
//	
//	for(int i=0; i<inverse_matrix_log.size(); i++){
//		for(int j=0; j<inverse_matrix_log[0].size(); j++){
//
//			cout << inverse_matrix_log[i][j] << " ";
//		}
//		cout << endl;
//	}
//	cout << "-------------------------------------------------------" << endl;
//}

///////////////////////////////////////////////////////////////////////////////////////////////

////////////////////PARSER////////////////////////////////////////////////////////////////////

void command_line_parser(int argc, char** argv){
	
	const char* const short_opts = "hp:k:b:j:m:d:o:f:t:n:s";

	//Specifying the expected options
	const option long_opts[] ={
		{"help",      no_argument, nullptr,  'h' },
		{"param",      required_argument, nullptr,  'p' },
		{"ntop",      required_argument, nullptr,  'n' },
		{"kmer",   required_argument, nullptr,  'k' },
		{"freq",   required_argument, nullptr,  'f' },
		{"distance",   required_argument, nullptr,  'd' },
		{"bed",    required_argument, nullptr,  'b' },
		{"ordering",    required_argument, nullptr,  'o' },
		{"jaspar",   required_argument, nullptr,  'j' },
		{"mf",   required_argument, nullptr,  'm' },
		{"twobit",   required_argument, nullptr,  't' },
		{"ss",   no_argument, nullptr,  's' },
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
					   cerr << "ERROR: Wrong -o parameter inserted.\n\n" << endl;
					   display_help();
					   exit(1);}
				   break;
			case 'k' : kmers.clear();
				   kmers = string(optarg);
				   break;
			case 'f' : freq_treshold = stod(optarg);
				   if(freq_treshold <= 0 || freq_treshold >= 1){
					   cerr << "ERROR: please insert a frequency treshold between 0 and 1.\n\n" << endl;
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
		cerr <<"ERROR: "<< buf << " file does not exist!\n"<< endl;
		display_help();
		exit(1);	
	}
	return 0;
}

void check_input_file(){

	if(MFASTA_FILE.size() != 0 && (BED_FILE.size() != 0 || JASPAR_FILE.size() != 0 || TWOBIT_FILE.size() != 0)){

		cerr << "FATAL ERROR: Too many input arguments!\nPlease insert Multifasta file or Bed, Twobit and Jaspar file.\n" << endl;
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
	cerr << "\n --help || -h show this message" << endl;
	cerr << "\n --bed || -b <file_bed>: input bed file" << endl;
	cerr << "\n --kmer || -k <n1,n2,..,nN>: to select k-mers length for the analysis(DEFAULT: 6,8,10) " << endl;
	cerr << "\n --twobit || -t <file_twobit>: input twobit file" << endl;
	cerr << "\n --jaspar || -j <JASPAR_file>: input JASPAR file" << endl;
	cerr << "\n --param || -p <half_length>: half_length to select bases number to keep around the chip seq signal (DEFAULT: 150) " << endl;
	cerr << "\n --ntop || -n <number>: to decide the top n oligos to classify in positional sequence occurrences (DEFAULT: 10) " << endl;
	cerr << "\n --mf || -m <multifasta-file>: use multifasta instead of bed file [ -j,-b,-t,-p options not needed ]" << endl;
	cerr << "\n -s || --ss as input to make the analysis along the single strand. (DEFAULT: double strand)" << endl;
	cerr << "\n -o || --ordering 'p' to order the top N oligos by p-value and not by occurrences. (DEFAULT: ordering by occurrences)" << endl;
	cerr << "\n --distance || -d <n1,n2,...,nN> to select the hamming distances. (DEFAULT: 1,2,3)" << endl;
	cerr << "\n --freq || -f <n1> to set the frequence treshold to calculate the z_score. (DEFAULT: 0.02)" << endl;
	cerr << endl;

	exit(EXIT_SUCCESS);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
