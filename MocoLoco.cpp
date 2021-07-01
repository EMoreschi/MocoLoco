#include "MocoLoco.h"

int main(int argc, char *argv[]){


	if(argc == 1){             //If arguments number is 1 means that no input file has been inserted - display help
		display_help();
	}

	command_line_parser(argc, argv);					//Parser function called to handle aguments
	GEP_path();							//Calling to GEP pathway

	return 0;
}

void  GEP_path(){

	if(MFASTA_FILE.size() == 0){	

		coordinator_class C;
		C.print_debug_GEP(C.GEP);
		map_class MAP(C.GEP,kmers,dist);
	}

	else{

		multifasta_class MULTI(MFASTA_FILE);
		map_class MAP(MULTI.GEP,kmers,dist);
	}		
}

void map_class::check_kmer_dist(){
	
	if(kmers_vector.size() != distance_vector.size()){
		
		cerr << "\nERROR: Please insert an equal number of k-mers and distance parameters!" << endl;
		display_help();
		exit(1);
	}

}

void bed_class::read_line(string line){				//Read line function: it takes in input each line from BED file 

istringstream mystream(line);					//Split the line word by word and extract chromosome coordinates (chr, start, end)
mystream >> chr_coord >> start_coord >> end_coord;		

}

void bed_class::centering_function ( unsigned int start,  unsigned int end, unsigned int half_length, const unsigned int overhead){	//Centering function: in takes start and end coordinate and re-sets them -
//following an input half_length value (overhead added to the end)
unsigned int center = (start + end)/2;						
start_coord = center - half_length;			//No overhead for start
end_coord = center + half_length +overhead;		//Overhead for end
}

void bed_class::flag_control( unsigned int start,  unsigned int end){ 	//Flag control function: start coordinates must be < then end coordinates

if(start > end){		//if start coordinates are > or == then end coordinates, flag is setted to 0
	flag = 0;
}
else{ flag = 1;}
}

void coordinator_class::GEP_creation(string Bed_file, string Twobit_file, vector<bed_class> &GEP){		//Function to read BED and 2Bit files and create GEP object vector

cout << "\n- [1] Extract bed coordinate sequences from reference genome  \n";

ifstream in(Bed_file); 						//Opening file in lecture mode
TwoBit * tb;				//Creating a TwoBit* variable called tb
tb = twobit_open(Twobit_file.c_str());					//Opening 2Bit file with twobit_open function and saved in tb 
string line; 							//defining line string

unsigned int n_line = 1;							//line counter initialization

while(getline(in,line)){  					//reading input file line by line with getline function

	if(line.empty())   					//if line is empty or commented --> continue
		continue;
	if(line[0]=='#')
		continue;

	bed_class new_class(half_length,line,tb, n_line);  //Called the object constructor passing the Bed line, half_length P, twobit file tb, and the line counter n_line
	GEP.emplace_back(new_class);				//Put the new object in the GEP vector with emplace function

	n_line = n_line + 1;					//pass to next line 

}
}

void coordinator_class::oligos_vector_creation(vector<oligo_class> &oligos_vector, vector<vector<double>> matrix_log, vector<vector<double>> matrix_log_inverse, vector<bed_class> GEP){

cout << "- [5] Analyzing sequences using Jaspar matrix\n";

for(unsigned int i=0; i<GEP.size(); i++){
	string sequence = GEP[i].return_sequence(GEP[i]);
	string chr_coord = GEP[i].return_chr_coord_GEP();
	unsigned int start_coord = GEP[i].return_start_coord_GEP();


	oligo_class SHIFTING(matrix_log, sequence, chr_coord, start_coord, '+');
	oligos_vector.emplace_back(SHIFTING);

	if(DS == 1){

		oligo_class SHIFTING(matrix_log_inverse, sequence, chr_coord, start_coord, '-');
		oligos_vector.emplace_back(SHIFTING);
	}
}	

cout << "- [6] Selecting the best Jaspar's oligo for each sequence \n";
}

void coordinator_class::best_strand(vector<oligo_class> oligos_vec){

if(DS == 1){
	vector<oligo_class> comparison;
	for(unsigned int i=0; i<oligos_vec.size(); i+=2){

		double best_score_norm_positive = oligos_vec[i].return_best_score_normalized();
		double best_score_norm_negative = oligos_vec[i+1].return_best_score_normalized(); 
		if(best_score_norm_positive >= best_score_norm_negative){

			comparison.emplace_back(oligos_vec[i]);
		}
		else{
			comparison.emplace_back(oligos_vec[i+1]);
		}
	}
	oligos_vector = comparison;
}
}

void oligo_class::shifting(vector<vector<double>> matrix, string sequence, unsigned int s_iterator){

double sum_scores = 0;

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

	oligo_scores.emplace_back(sum_scores);
	shifting(matrix, sequence, s_iterator+1);

}

}

void matrix_class::read_JASPAR(string JASPAR_FILE){			//Function to read JASPAR PWM file, extract values and create a matrix class

cout << "- [2] Reading JASPAR MATRIX file and extracting values\n";

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

for (unsigned int i = 0; i < matrix[0].size(); i++) {			//From 0 to number of columns of line 0
	for (unsigned int j = 0; j < 4; j++){				//From 0 to 4 (line number)

		sum = sum + matrix[j][i];			//Calculate the sum of columns
	}

	col_sum.emplace_back(sum);				//Put the column sum in vector col_sum
	sum = 0;						//Restore the sum to 0 for the next column
}
return col_sum;
}

void matrix_class::matrix_normalization_pseudoc(vector<vector<double>> matrix, double p){  

cout << "- [3] Jaspar Matrix normalization\n";

double norm;							//Norm variable initialized
vector<double> col_sum = find_col_sum(matrix);

for (unsigned int i = 0; i < matrix.size(); i++) {		//From 0 to number of matrix lines

	vector<double> baseQ;				//baseQ vector to store the lines initialized
	for (unsigned int j = 0; j < matrix[i].size(); j++){	//From 0 to number of matrix columns

		norm = matrix[i][j]/col_sum[j];		//Put matrix value (divided for the corresponding column sum) into double variable norm
		baseQ.emplace_back(norm + p);		//Put norm value (with p added) in baseQ vector
	}

	norm_matrix.emplace_back(baseQ);	//Put baseQ vector (which carries line values) in norm_matrix
}
}

void matrix_class::matrix_normalization(vector<vector<double>> matrix){

//	cout << "Second Matrix normalization...\n";

vector<double> col_sum = find_col_sum(matrix);

for (unsigned int i = 0; i < matrix.size(); i++) {		//From 0 to number of matrix lines

	vector<double> baseQ;				//baseQ vector to store the lines initialized
	for (unsigned int j = 0; j < matrix[i].size(); j++){	//From 0 to number of matrix columns

		norm_matrix[i][j] = matrix[i][j]/col_sum[j];	//Substitution of first normalized values with new normalized ones
	}
}
}

void matrix_class::matrix_logarithmic(vector<vector<double>> matrix){

for(unsigned int i=0; i < matrix.size(); i++){
	vector<double> baseQ;
	double value_log;

	for(unsigned int j=0; j < norm_matrix[i].size(); j++){

		value_log = log(norm_matrix[i][j]);
		baseQ.emplace_back(value_log);
	}
	matrix_log.emplace_back(baseQ);
}

cout << "- [4] Jaspar Matrix reverse complement determination to analize the reverse strand\n";
}


vector<vector<double>> matrix_class::reverse_matrix(vector<vector<double>> matrix){

vector<vector<double>> inv_matrix = matrix;
reverse(inv_matrix.begin(), inv_matrix.end());
for (int i = 0; i < 4; i++) {		//From 0 to number of matrix lines
	vector<double> baseQ;
	reverse(inv_matrix[i].begin(), inv_matrix[i].end());
}
return inv_matrix;

}

void oligo_class::find_minmax(vector<vector<double>> matrix){

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

unsigned int oligo_class::find_best_score(vector<double> oligo_scores){

best_score = *max_element(oligo_scores.begin(), oligo_scores.end());

vector<int> positions;
vector<int> dist_center;
unsigned int matches = 0;
int min_distance;
vector<int>::iterator itr;

for(unsigned int i=0; i < oligo_scores.size(); i++){

	if(oligo_scores[i] == best_score){

		matches = matches + 1;
		positions.emplace_back(i);
	}
}
if(matches > 1){ 

	for (int& p: positions){
		int distance;
		distance = abs( p - half_length); 
		dist_center.emplace_back(distance);
	}

	min_distance = *min_element(dist_center.begin(), dist_center.end());
	itr = find(dist_center.begin(),dist_center.end(), min_distance);
	unsigned int index = distance(dist_center.begin(), itr);
	return positions[index];

}
return positions[0];
}

void oligo_class::best_score_normalization(){

	best_score_normalized = 1 + ((best_score - max_possible_score)/(max_possible_score - min_possible_score));

}

void oligo_class::find_best_sequence(string sequence, unsigned int local_position, unsigned int length){

	best_oligo_seq = sequence.substr(local_position,length);
}

void oligo_class::find_coordinate(unsigned int local_position, unsigned int length, string chr_coord_GEP, unsigned int start_coord_GEP){

	chr_coord_oligo = chr_coord_GEP;
	start_coord_oligo = start_coord_GEP + local_position;
	end_coord_oligo = start_coord_oligo + length;

}

void coordinator_class::centering_oligo(){

	TwoBit * tb;
	tb = twobit_open(TWOBIT_FILE.c_str());
	int center_oligo ;

	for(unsigned int i=0; i<oligos_vector.size(); i++){
		center_oligo = oligos_vector[i].return_start_coord_oligo() + matrix_log[0].size()/2;
		GEP[i].centering_function(center_oligo,center_oligo,half_length,0);
		GEP[i].extract_seq(tb,0);
	}
}

void bed_class::extract_seq(TwoBit* tb, unsigned int n_line){			//Extract sequence function: Extract, from Twobit hg38 genome, the DNA sequence with (chr, start, end) coordinates -
	//extracted from Bed line
	if(flag == 1){								//CONTROL: if flag is 1 means that the current line has starting coordinate > end coordinate, so it is correct
		string chrom = chr_coord;		//Put in chrom the string of chr_coord
		sequence = twobit_sequence(tb,chrom.c_str(),start_coord,end_coord-1); 	//Extract the sequence from the object with the twobit_sequence function
	}
	else {		
		cerr << "WARNING: the line " << n_line <<" is omitted because starting coordinates > end coordinates, please check your BED file!" << "\n";
		//if flag is not 1 means that the current line has starting coordinate < end coordinate: PRINT WARNING!		
	}
}

vector<unsigned int> map_class::generic_vector_creation(string numbers){

	int index;
	vector<unsigned int> vec;

	while(index != -1){
		index = numbers.find(",");
		vec.emplace_back(stoi(numbers.substr(0,index)));
		numbers.erase(0,index+1);
	}

	return vec;
}

void map_class::table_creation_orizzontal(vector<bed_class> GEP){ 

	if (MFASTA_FILE.size() ==0) 
		cout << "- [7] Counting all k-mers occurrences for sequence and positions  \n";
	else
		cout << "- [4] Counting all k-mers occurrences for sequence and positions  \n";

	for(unsigned int k=0; k<kmers_vector.size(); k++){

		for(unsigned int j=0; j<GEP.size(); j++){

			string sequence = GEP[j].return_sequence(GEP[j]);

			for(unsigned int i=0; i < (sequence.size() - kmers_vector[k] + 1); i++){

				string bases = sequence.substr(i,kmers_vector[k]);
				or_ver_kmer_count(bases,orizzontal_plus,orizzontal_minus);
			}
		}
		orizzontal_plus_debug.emplace_back(orizzontal_plus);
		orizzontal_minus_debug.emplace_back(orizzontal_minus);
		orizzontal_plus.clear();
		orizzontal_minus.clear();
	}
}

void map_class::table_creation_vertical(vector<bed_class> GEP){

	string seq_length = GEP[0].return_sequence(GEP[0]);

	for(unsigned int k=0; k<kmers_vector.size(); k++){

		vector<unsigned int> tot_freq_vec;

		for(unsigned int i=0; i < (seq_length.size() - kmers_vector[k] + 1); i++){

			unsigned int tot_freq = 0;

			for(unsigned int j=0; j<GEP.size(); j++){

				string sequence = GEP[j].return_sequence(GEP[j]);
				string bases = sequence.substr(i,kmers_vector[k]);
				vertical_kmer_count(bases, vertical_plus,vertical_minus, tot_freq);
			}

			if(DS==1){
				select_best(vertical_plus);
			}
			maps_vector_positions_plus.emplace_back(vertical_plus);
			maps_vector_positions_minus.emplace_back(vertical_minus);
			vertical_plus.clear();
			vertical_minus.clear();
			tot_freq_vec.emplace_back(tot_freq);
		}

		vector_kmers_maps_plus.emplace_back(maps_vector_positions_plus);
		vector_kmers_maps_minus.emplace_back(maps_vector_positions_minus);
		tot_freq_matrix.emplace_back(tot_freq_vec);
		maps_vector_positions_plus.clear();
		maps_vector_positions_minus.clear();
	}
}

void map_class::or_ver_kmer_count(string bases,unordered_map<string,unsigned int> &plus, unordered_map<string,unsigned int> &minus){

	unordered_map<string,unsigned int>::iterator it_plus;
	unordered_map<string,unsigned int>::iterator it_minus;
	it_plus = plus.find(bases);
	check_palindrome(bases);
	it_minus = minus.find(reverse_bases);

	if(it_plus!=plus.end()){

		it_plus->second++;
		it_minus->second++;
	}

	else{

		plus.insert({bases,1});
		minus.insert({reverse_bases,1});
	}


	bases.clear();
	reverse_bases.clear();
}

void map_class::vertical_kmer_count(string bases,map<pair<string,string>,pair<unsigned int, unsigned int>>&plus, map<pair<string,string>,pair<unsigned int, unsigned int>> &minus, unsigned int& tot_freq){



	map<pair<string,string>,pair<unsigned int, unsigned int>>::iterator it_plus;
	map<pair<string,string>,pair<unsigned int, unsigned int>>::iterator it_plus_rev;
	map<pair<string,string>,pair<unsigned int, unsigned int>>::iterator it_minus;
	map<pair<string,string>,pair<unsigned int, unsigned int>>::iterator it_minus_rev;

	bool pal = check_palindrome(bases);

	if(DS == 1){	
		if(!pal){
			tot_freq = tot_freq+2;
		}
		else{
			tot_freq++;
		}
	}
	else{
		tot_freq++;
	}

	pair<string,string> pair_bases;
	pair_bases.first= bases;
	pair_bases.second= reverse_bases;

	pair<string,string> pair_bases_reverse;
	pair_bases_reverse.first= reverse_bases;
	pair_bases_reverse.second= bases;

	it_plus = plus.find(make_pair(bases, reverse_bases));
	it_plus_rev = plus.find(make_pair(reverse_bases, bases));

	it_minus = minus.find(pair_bases);
	it_minus_rev = minus.find(pair_bases_reverse);

	if(it_plus!=plus.end()){

		it_plus->second.first++;
		it_minus->second.first++;
	}	
	if (it_plus==plus.end() && it_plus_rev != plus.end()) {
		it_plus_rev->second.second++;
		it_minus_rev->second.second++;

	}

	else{

		plus.insert({{bases,reverse_bases},{1,0}});
		minus.insert({{bases,reverse_bases},{1,0}});
	}


	bases.clear();
	reverse_bases.clear();
}

void map_class::select_best(map<pair<string,string>,pair<unsigned int,unsigned int>>& vertical_plus){

	map<pair<string,string>,pair<unsigned int,unsigned int>> copy;

	for(map<pair<string,string>,pair<unsigned int,unsigned int>>::iterator it = vertical_plus.begin(); it!=vertical_plus.end(); it++){

		if(it->second.first < it->second.second){

			string oligo1 = it->first.second;	
			string oligo2 = it->first.first;
			unsigned int occ1 = it->second.second;
			unsigned int occ2 = it->second.first;

			copy.insert({{oligo1,oligo2},{occ1,occ2}});		
		}
		else{

			copy.insert({{it->first.first, it->first.second},{it->second.first,it->second.second}});		
		}
	}
	vertical_plus.clear();
	vertical_plus = copy;
}

void p_value_class::N2_calculation(unordered_map<string,unsigned int> orizzontal_map){

	total_oligo_N2 = 0;
	total_oligo_N2 = accumulate(begin(orizzontal_map), end(orizzontal_map), 0, [] (unsigned int val, const unordered_map<string,int>::value_type& p) {return val + p.second;});
}

bool map_class::check_palindrome(string bases){

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

		if(sequences[i].size() != size){

			cerr << "Sequences are not of the same length!" << endl;
			exit(1);
		}
	}
}

void multifasta_class::extract_sequences(string MFasta_file){

	cout << "\n- [1] Extracting sequences from MultiFasta file \n";

	ifstream file(MFasta_file);
	string line;
	string current_sequence;
	bool first_line = 1;

	while(getline(file,line)){

		if(line[0] == '>' && !first_line){

			sequences.emplace_back(current_sequence);
			current_sequence.clear();

		}

		else if (!first_line){

			if(line[0] != ' ' && line.size() != 0){	
				transform(line.begin(), line.end(), line.begin(), ::toupper);	
				current_sequence = current_sequence + line; 
			}
		}

		first_line = 0;	
	}
	sequences.emplace_back(current_sequence);
}

void multifasta_class::GEP_creation_MF(vector<string> sequences){

	cout << "- [3] Sorting Multifasta sequences\n";

	for(unsigned int i=0; i<sequences.size(); i++){

		bed_class new_class(sequences[i]);
		GEP.emplace_back(new_class);
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

void map_class::HUMMING_MATRIX_creation(){
	
	for(unsigned int j=0; j<P_VALUE_MATRIX.size(); j++){

		ofstream outfile = outfile_header_hamming(j);
		
		for(unsigned int i=0; i<P_VALUE_MATRIX[j].size(); i++){

			multimap<pair<unsigned int, unsigned int>, pair<string,string>> vertical_multimap = P_VALUE_MATRIX[j][i].return_vertical_multimap();

			hamming_class H(vertical_multimap,distance_vector[j],i,tot_freq_matrix[j][i],orizzontal_plus_debug[j], orizzontal_minus_debug[j], outfile);
			HUMMING_VECTOR.emplace_back(H);
		}

		HUMMING_MATRIX.emplace_back(HUMMING_VECTOR);
		HUMMING_VECTOR.clear();
		outfile.close();
	}
}

void hamming_class::find_best_oligos(){
	
	multimap<pair<unsigned int,unsigned int>, pair<string,string>>::reverse_iterator it_rev = vertical_multimap.rbegin();
	real_best_oligo_occurrences = (it_rev->first.first + it_rev->first.second);

	while(it_rev->first.first + it_rev->first.second == real_best_oligo_occurrences){

		best_oligos.emplace_back(it_rev->second.first);
		it_rev++;
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

double hamming_class::frquence_2_calculation(unordered_map<string,unsigned int> orizzontal_map_plus, unordered_map<string,unsigned int> orizzontal_map_minus){
	

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

	unsigned int counter_A = 0;
	unsigned int counter_C = 0;
	unsigned int counter_G = 0;
	unsigned int counter_T = 0;
	vector<unsigned int> vec_A, vec_C, vec_G, vec_T;
	
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

vector<vector<unsigned int>> hamming_class::return_PWM_hamming(){

	return PWM_hamming;
}

void matrix_class::debug_matrix(matrix_class M){		//Debugging of matrices: calling print matrix function

	M.print_debug_matrix(matrix, " ");
	M.print_debug_matrix(norm_matrix, " NORMALIZED");
	M.print_debug_matrix(inverse_norm_matrix, " INVERSE NORMALIZED MATRIX");
	M.print_debug_matrix(matrix_log, " LOGARITHMIC MATRIX");
	M.print_debug_matrix(inverse_matrix_log, " INVERSE LOGARITHMIC MATRIX");
}

void matrix_class::print_debug_matrix(vector<vector<double>> matrix, string type){		//Print matrix function

	cout << "\n" << matrix_name << " " << tf_name << type << ":" << endl;

	for(unsigned int i=0; i < matrix.size(); i++){
		for(unsigned int j=0; j<matrix[i].size(); j++){

			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}

}


void coordinator_class::print_debug_GEP(vector<bed_class> GEP){			//Debug function: Print the GEP vector to control the working flow
	
	alias_file = (TWOBIT_FILE.erase(0,TWOBIT_FILE.find_last_of("/")+1)+"_"+ JASPAR_FILE.erase(0,JASPAR_FILE.find_last_of("/")+1)+"_"+ BED_FILE.erase(0,BED_FILE.find_last_of("/")+1));
	ofstream outfile;	
	JASPAR_FILE = JASPAR_FILE.erase(JASPAR_FILE.find_last_of("."), JASPAR_FILE.size());
	outfile.open(alias_file);
	
	for(unsigned int i=0; i<GEP.size(); i++){
		string chr_coord = GEP[i].return_chr_coord();
		unsigned int start_coord = GEP[i].return_start_coord();
		unsigned int end_coord = GEP[i].return_end_coord();
		outfile << chr_coord << "\t" << start_coord << "\t" << end_coord << endl;	//Printing chr, start and end coordinates
	}
	outfile.close();
	
	BED_FILE = BED_FILE.erase(BED_FILE.find_last_of("."), BED_FILE.size());
	outfile.open(alias_file+".fasta");
	for(unsigned int i=0; i<GEP.size(); i++){
		string chr_coord = GEP[i].return_chr_coord();
		unsigned int start_coord = GEP[i].return_start_coord();
		unsigned int end_coord = GEP[i].return_end_coord();
		string sequence = GEP[i].return_sequence(GEP[i]);					//Printing sequence
		outfile << ">" << chr_coord << ":" << start_coord << "-" << end_coord << endl;	//Printing chr, start and end coordinates
		outfile << sequence << endl;
	}
	outfile.close();

}

void oligo_class::oligos_vector_debug(oligo_class oligos_vector){	//Debug function to print the best oligo features

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

			print_debug_PWM_hamming(outfile, j);
		}

		else{

			outfile.open(to_string(kmers_vector[j])+"-mers_PWM_hamming_matrices_"+alias_file+"DS.txt");

			print_debug_PWM_hamming(outfile, j);
		}
	}
}

void map_class::print_debug_PWM_hamming(ofstream& outfile, unsigned int j){

	outfile << "#PWM Matrices calculated from the best oligo for each position and his hamming distanced oligos" << endl << endl;
	string best_oligo;
	unsigned int neighbour_numb;
	vector<vector<unsigned int>> PWM_hamming;
	string ATCG = "ATCG";

	for(unsigned int position = 1; position <= HUMMING_MATRIX[0].size(); position++){

		best_oligo = HUMMING_MATRIX[j][position-1].return_real_best_oligo();	
		neighbour_numb = HUMMING_MATRIX[j][position-1].return_similar_oligo_size();
		PWM_hamming = HUMMING_MATRIX[j][position-1].return_PWM_hamming();

		outfile << "#Position " << position << ": \n#PWM calculated from oligo " << best_oligo << " and his " << neighbour_numb << " hamming distanced neighbours. " << endl << endl;

		for(unsigned int i = 0; i< PWM_hamming.size(); i++){
			
			outfile << ATCG[i] << "\t" << "[" << "\t";

			for(unsigned int j = 0; j<PWM_hamming[i].size(); j++){

				outfile << PWM_hamming[i][j] << "\t";
			}
			outfile << "]" << endl;
		}
		outfile << endl;
		outfile << "-------------------------------------------------------------------" << endl;	
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////

////////////////////PARSER////////////////////////////////////////////////////////////////////

void command_line_parser(int argc, char** argv){
	
	const char* const short_opts = "hp:k:b:j:m:d:o:t:n:s";

	//Specifying the expected options
	const option long_opts[] ={
		{"help",      no_argument, nullptr,  'h' },
		{"param",      required_argument, nullptr,  'p' },
		{"ntop",      required_argument, nullptr,  'n' },
		{"kmer",   required_argument, nullptr,  'k' },
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
	cerr << endl;

	exit(EXIT_SUCCESS);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
