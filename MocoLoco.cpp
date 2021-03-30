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
		map_class MAP(C.GEP,kmers);
	}

	else{

		multifasta_class MULTI(MFASTA_FILE);
		map_class MAP(MULTI.GEP,kmers);
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

//	cout << "Calculating the log matrix...\n";

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

void map_class::kmers_vector_creation(string kmers){

	int index;
	
	while(index != -1){
		index = kmers.find(",");
		kmers_vector.emplace_back(stoi(kmers.substr(0,index)));
		kmers.erase(0,index+1);
	}
}

void map_class::table_creation(unordered_map<string,int> moco_table, vector<int> kmers_vector, vector<bed_class> GEP){ 
	if (MFASTA_FILE.size() ==0) 
		cout << "- [7] Counting all k-mers occurrences for sequence and positions  \n";
	else
		cout << "- [4] Counting all k-mers occurrences for sequence and positions  \n";

	for(unsigned int k=0; k<kmers_vector.size(); k++){

		maps_vector_positions.clear();

		for(unsigned int j=0; j<GEP.size(); j++){

			string sequence = GEP[j].return_sequence(GEP[j]);
			unordered_map<string,int>::iterator it;

			for(unsigned int i=0; i < (sequence.size() - kmers_vector[k] + 1); i++){

				string bases = sequence.substr(i,kmers_vector[k]);
				it = moco_table.find(bases);
				bool palindrome = check_palindrome(bases);

				if(j==0){

					maps_vector_positions.emplace_back(moco_pos);
					moco_pos.clear();

				}

				unordered_map<string,int>::iterator it_pos;
				it_pos = maps_vector_positions[i].find(bases);

				if (!palindrome && DS){
					unordered_map<string, int>::iterator it_rev = moco_table.find(reverse_bases);
					unordered_map<string, int>::iterator it_rev_pos = maps_vector_positions[i].find(reverse_bases);

					if (it != moco_table.end()){
						it->second++;
						it_rev->second++;
					}
					else{

						moco_table.insert( pair<string,int>(bases,1) );
						moco_table.insert( pair<string,int>(reverse_bases,1) );

					}


					if(it_pos != maps_vector_positions[i].end()){

						it_pos->second++;
						it_rev_pos->second++;
					}

					else {

						maps_vector_positions[i].insert(pair<string,int>(bases,1));
						maps_vector_positions[i].insert(pair<string,int>(reverse_bases,1));
					}
				}

				else{
					if (it != moco_table.end()){
						it->second++;
					}
					else{
						moco_table.insert( pair<string,int>(bases,1) );

					}


					if(it_pos != maps_vector_positions[i].end()){

						it_pos->second++;
					}
					else{

						maps_vector_positions[i].insert(pair<string,int>(bases,1));						
					}
				}

				bases.clear();
				reverse_bases.clear();
			}
		}

		maps_vector_debug.emplace_back(moco_table);
		moco_table.clear();

		v_v_maps.emplace_back(maps_vector_positions);
	}
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

void map_class::print_debug_maps(vector<unordered_map<string,int>> maps_vector_debug, vector<int> kmers_vector, string direction){
	for(unsigned int i=0; i<maps_vector_debug.size(); i++){
		
		ofstream outfile;
		outfile.open(to_string(kmers_vector[i])+"-mer_ordered_map_"+alias_file+direction+".txt");	
			multimap<int,string> moco_multimap;
		for (unordered_map<string,int>::iterator it = maps_vector_debug[i].begin(); it !=maps_vector_debug[i].end(); it++) {
			moco_multimap.insert( { it->second, it->first });
		}
		for (multimap<int,string>::reverse_iterator rev_it = moco_multimap.rbegin(); rev_it != moco_multimap.rend(); rev_it++) {
			outfile << rev_it->second <<"\t"<< rev_it->first <<"\n";
		}
		outfile.close();
	}
}

void map_class::print_debug_maps_positions(){

	ofstream outfile;
	for(unsigned int j=0; j<v_v_maps.size(); j++){

		outfile.open("vertical_ordered_map_"+to_string(kmers_vector[j])+"kmers"+alias_file+".txt");

		outfile << "# Maps vector with kmers occurences counted for positions in sequence (for k = " << kmers_vector[j] << "):" << endl;

		for(unsigned int i=0; i<v_v_maps[j].size(); i++){

			outfile << "### kmers occurred in position " << i << ":" << endl;
			multimap<int,string> moco_multimap;

			for (unordered_map<string,int>::iterator it = v_v_maps[j][i].begin(); it !=v_v_maps[j][i].end(); it++) {
				moco_multimap.insert({it->second, it->first});
			}
			multimap<int,string>::reverse_iterator rev_it = moco_multimap.rbegin();

			for_each(rev_it, next(rev_it,10),[&outfile](pair<int,string> element){

					outfile << element.second << "\t" << element.first << "\n";
					});

		}

		outfile.close();
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////

////////////////////PARSER////////////////////////////////////////////////////////////////////

void command_line_parser(int argc, char** argv){
	
	const char* const short_opts = "hp:k:b:j:m:t:s";

	//Specifying the expected options
	const option long_opts[] ={
		{"help",      no_argument, nullptr,  'h' },
		{"param",      required_argument, nullptr,  'p' },
		{"kmer",   required_argument, nullptr,  'k' },
		{"bed",    required_argument, nullptr,  'b' },
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
			case 'b' : BED_FILE = string(optarg);
				   is_file_exist(BED_FILE, "--bed || -b");
				   break;
			case 'j' : JASPAR_FILE = string(optarg);
				   is_file_exist(JASPAR_FILE, "--jaspar || -j");
				   break;
			case 't' : TWOBIT_FILE = string(optarg);
				   is_file_exist(TWOBIT_FILE, "--twobit || -t ");
				   break;
			case 'k' : kmers.clear();
				   kmers = string(optarg);
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

void display_help() 						//Display help function
{
	cerr << "\n --help || -h show this message" << endl;
	cerr << "\n --bed || -b <file_bed>: input bed file" << endl;
	cerr << "\n --kmer || -k <n1,n2,..,nX>:input at least one k-mer length (DEFAULT: 6,8,10) " << endl;
	cerr << "\n --twobit || -t <file_twobit>: input twobit file" << endl;
	cerr << "\n --jaspar || -j <JASPAR_file>: input JASPAR file" << endl;
	cerr << "\n --param || -p <half_length>: input half_length to select bases number to keep around the chip seq signal (DEFAULT: 150) " << endl;
	cerr << "\n --mf || -m <multifasta-file>: use multifasta instead of bed file [ -j,-b,-t,-p options not needed ]" << endl;
	cerr << "\n -s || --ss as input to make the analysis along the single strand. (DEFAULT: double strand)" << endl;
	cerr << endl;

	exit(EXIT_SUCCESS);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
