#include "MocoLoco.h"
//#include "Profiling.h"
#include <sys/resource.h>

int main(int argc, char *argv[]) {
  Timer timer;
  //Instrumentor::Get().BeginSession("MocoLoco");
  {

    // If arguments number is 1 means that no input file has been inserted -
    // display help
    if (argc == 1) {

      display_help();
    }
    // Collect all the parameters given as input
    command_line_parser(argc, argv);

    // This line can be read as "If MFASTA_FILE is empty do BED_pat() else do MULTIFA_path()"
    (MFASTA_FILE.empty()) ? (BED_path()) : (MULTIFA_path());

    return 0;
  }
  //Instrumentor::Get().EndSession();
}

// Function to choose the pathway to follow. 2 input options:
// 1) Bed-Twobit-Jaspar input
// 2) Multifasta input
void BED_path() {
  // PROFILE_FUNCTION();

  tb = twobit_open(TWOBIT_FILE.c_str());
  coordinator_class C;
  twobit_close(tb);
  // Create a .fasta file to check if the coordinates and the sequences
  // extracted are correct
  C.print_GEP(C.GEP);
  // Reading k-mers, distance and freq in input and saving them into vectors
  kmers_vector = generic_vector_creation(kmers);
  distance_vector = generic_vector_creation(dist);
  freq_vector = freq_vector_creation(freq_threshold);
  
  if(distance_vector.size() != kmers_vector.size() ||
      distance_vector.size() != freq_vector.size() ||
      kmers_vector.size() != freq_vector.size()){
    cerr << endl 
    << "ERROR! The number of kmers must be equal to the number of hamming distances and frequencies" << endl << endl;
    exit(1);
  }
  // Vector len contains all the lengths of sequences for each kmer 
  for(unsigned int i = 0; i < kmers_vector.size(); i++){
    len.emplace_back(C.GEP[0].sequence.size() - kmers_vector[i] + 1);
  }
  //In this class horizontal and vertical maps are created, these maps will  
  //be useful later on 
  MapClass M(C.GEP);
  //For each kmer
  for (unsigned int i = 0, d = 0; i < kmers_vector.size(); i++, d++) {
    vector<PvalueClass> P_vector;
    vector<string> seed_oligo;
    //For each position in the sequence
    for (unsigned int j = 0; j < len[i]; j++) {
      double Pval = 0;
      unsigned int counter = 0;
      vector<string> pos_oligo_vec;
      // Loop for eventually other matrices that are hidden by the best motif
      while(counter < max_matrix) {

        cout << "Position: " << j << endl;
          
        //For each oligo present in the vertical map
        for (multimap<int, string>::iterator it =
        M.vector_positions_occurrences[i][j].begin();
        it != M.vector_positions_occurrences[i][j].end(); it++) {
          //In the PvalueClass to each oligo is associated its pvalue
          PvalueClass P(C.GEP, it, M.vector_map_hor[i], i);
          P_vector.push_back(P);
        }
        //The element in P_vector are ordered on the basis of pvalues
        //otherwise the normal ordering (by occurrences in vertical map)
        //is maintained
        sort(begin(P_vector), end(P_vector), comp_occ);

        if (ordering == "p") {
          sort(begin(P_vector), end(P_vector), comp);
        }
        // Debug for PValueClass
        // DVector(P_vector, j);
          
        //Creation of clusters of oligos at hamming distance
        //and creation of PWM for each position
        HammingClass H(P_vector[0].oligo,
                        M.vector_map_ver[i],
                        M.vector_positions_occurrences[i][j],
                        M.vector_map_hor[i], j, d);
        // Perform the expectation maximization algorithm 
        // (https://en.wikipedia.org/wiki/Expectation–maximization_algorithm)
        if (!exp_max.empty()){
          EMClass E(H.cluster_map, H.PWM_hamming, 
                    M.vector_map_hor[i]);
        }
        //If the frequence of seed oligo is higher than a threshold
        //the z_score is calculated

        if(H.freq1 >= freq_vector[i]){
          z_test_class Z(H.PWM_hamming, C.GEP, 
                          j + 1, len[i]);
          Pval = Z.Zpvalue_bonf;
            
          //If it is the first cycle of while loop or if the pval is lower 
          //than a certain threshold the z_score and PWM are calculated
          if(Pval <= (z_pval_threshold * len[i])){
              seed_oligo.emplace_back(P_vector[0].oligo);
              Z_TEST_VECTOR.emplace_back(Z);
              H_HAMMING_VECTOR.emplace_back(H);
          }
        }  
        P_vector.clear(); 
        counter++;
      }
      pos_oligo_vec.clear();
    }
    Z_TEST_MATRIX.emplace_back(Z_TEST_VECTOR);
    H_HAMMING_MATRIX.emplace_back(H_HAMMING_VECTOR);
    Z_TEST_VECTOR.clear();
    H_HAMMING_VECTOR.clear();
    
    //Outfile functions 
    Outfile_PWM_matrices(i, seed_oligo);
    Outfile_Z_score_values(i, seed_oligo);
  }
  RAM_usage();
}
  
// else if the input is a Multifasta file
void MULTIFA_path(){

  multifasta_class MULTIFA(MFASTA_FILE);
  // Reading k-mers, dist and freq in input and saving them into vectors
  kmers_vector = generic_vector_creation(kmers);
  distance_vector = generic_vector_creation(dist);
  freq_vector = freq_vector_creation(freq_threshold);
  
  if(distance_vector.size() != kmers_vector.size() ||
      distance_vector.size() != freq_vector.size() ||
      kmers_vector.size() != freq_vector.size()){
    cerr << endl 
    << "ERROR! The number of kmers must be equal to the number of hamming distances and frequencies" << endl << endl;
    exit(1);
  }
  // Vector len contains all the lengths of sequences for each kmer 
  for(unsigned int i = 0; i < kmers_vector.size(); i++){
    len.emplace_back(MULTIFA.GEP[0].sequence.size() - kmers_vector[i] + 1);
  }
  //In this class horizontal and vertical maps are created, these maps will  
  //be useful later on 
  MapClass M(MULTIFA.GEP);
  //For each kmer
  for (unsigned int i = 0, d = 0; i < kmers_vector.size(); i++, d++) {
    vector<PvalueClass> P_vector;
    vector<string> seed_oligo;
    //For each position in the sequence
    for (unsigned int j = 0; j < len[i]; j++) {
      double Pval = 0;
      unsigned int counter = 0;
      vector<string> pos_oligo_vec;
      // Loop for eventually other matrices that are hidden by the best motif
      while(counter < max_matrix) {

        cout << "Position: " << j << endl;
          
        //For each oligo present in the vertical map
        for (multimap<int, string>::iterator it =
        M.vector_positions_occurrences[i][j].begin();
        it != M.vector_positions_occurrences[i][j].end(); it++) {
          //In the PvalueClass to each oligo is associated its pvalue
          PvalueClass P(MULTIFA.GEP, it, M.vector_map_hor[i], i);
          P_vector.push_back(P);
        }
        //The element in P_vector are ordered on the basis of pvalues
        //otherwise the normal ordering (by occurrences in vertical map)
        //is maintained

        sort(begin(P_vector), end(P_vector), comp_occ);

        if (ordering == "p") {
          sort(begin(P_vector), end(P_vector), comp);
        }
        // Debug for PValueClass
        // DVector(P_vector, j);
          
        //Creation of clusters of oligos at hamming distance
        //and creation of PWM for each position
        HammingClass H(P_vector[0].oligo,
                        M.vector_map_ver[i],
                        M.vector_positions_occurrences[i][j],
                        M.vector_map_hor[i], j, d);
        // Perform the expectation maximization algorithm 
        // (https://en.wikipedia.org/wiki/Expectation–maximization_algorithm)
        if (!exp_max.empty()){
          EMClass E(H.cluster_map, H.PWM_hamming, 
                    M.vector_map_hor[i]);
        }
        
        //If the frequence of seed oligo is higher than a threshold
        //the z_score is calculated
        if(H.freq1 >= freq_vector[i]){
          z_test_class Z(H.PWM_hamming, MULTIFA.GEP, 
                          j + 1, len[i]);
          Pval = Z.Zpvalue_bonf;
            
          //If it is the first cycle of while loop or if the pval is lower 
          //than a certain threshold the z_score and PWM are calculated
          if(Pval <= (z_pval_threshold*len[i])){

              seed_oligo.emplace_back(P_vector[0].oligo);
              Z_TEST_VECTOR.emplace_back(Z);
              H_HAMMING_VECTOR.emplace_back(H);
          }
        }  
        P_vector.clear(); 
        counter++;
      }
      pos_oligo_vec.clear();
    }
    Z_TEST_MATRIX.emplace_back(Z_TEST_VECTOR);
    H_HAMMING_MATRIX.emplace_back(H_HAMMING_VECTOR);
    Z_TEST_VECTOR.clear();
    H_HAMMING_VECTOR.clear();
    
    //Outfile functions 
    Outfile_PWM_matrices(i, seed_oligo);
    Outfile_Z_score_values(i, seed_oligo);
  }

  RAM_usage();
}


void bed_class::read_line(string line, char ** result) {
  // PROFILE_FUNCTION();
  // Split the line word by word and extract chromosome coordinates (chr, start,
  // end)
  istringstream mystream(line);
  mystream >> chr_coord >> start_coord >> end_coord;
  
  bool found = false;
  // Check out if the format of the BED file chromosomes is the same of the 2bit genome
	for (int i = 0; result[i] != NULL; i++) {
		if (result[i] == chr_coord) {
			found = true;
      continue;
		}
	}

	if (!found) {
		cerr << endl << "Check out your BED file!" <<endl;
    cerr << endl << "In the BED file there is chromosome: " << chr_coord << endl;
    cerr << "This chromosome is not present in 2bit genome." << endl << endl;
    
    cerr << "Here a list of all the chromosomes present in .2bit file: " << endl << endl;
    for (int i = 0; result[i] != NULL; i++) {
      cerr << result[i] << endl;
    }
    cerr << endl;
    cerr << endl << "Be sure that the chromosomes in the bed file are of the same format" << endl;
    exit(1);
	}
}

void bed_class::centering_function(unsigned int start, unsigned int end,
                                   int half_length,
                                   const unsigned int overhead) {
  // PROFILE_FUNCTION();
  unsigned int center = (start + end) / 2;
  if (start > end) {

    flag = false;
  }

  else {
    flag = true;
  }
  // No overhead for start coordinates but overhead added to end coordinates
  start_coord = center - half_length;
  end_coord = center + half_length + overhead;
}

// Extract sequence function: Extract, from Twobit hg38 genome, the DNA sequence
// with (chr, start, end) coordinates extracted from Bed line
void bed_class::extract_seq(TwoBit *tb, unsigned int n_line) {
  // PROFILE_FUNCTION();
  // CONTROL: if flag is 1 means that the current line has starting coordinate >
  // end coordinate, so it is correct
  if (flag) {
    
    // Extract the sequence from the object with the twobit_sequence function
    sequence =
        twobit_sequence(tb, chr_coord.c_str(), start_coord, end_coord - 1);
  }

  // if flag is not 1 means that the current line has starting coordinate < end
  // coordinate: PRINT WARNING!
  else {
    err = true;
    cerr << "WARNING: the line " << n_line << " (" << chr_coord << ":"
         << start_coord << "-" << end_coord << ")"
         << " is omitted because starting coordinates > end coordinates, "
            "please check your BED file!"
         << "\n";
  }
}

// Function useful to normalize matrix scores and adding a pseudocount to them
void matrix_class::matrix_normalization_pseudoc(
    vector<vector<double>> &matrix) { // CHANGE: Possible use of reference
  // PROFILE_FUNCTION();

  // Calculate and save the scores column sum into a vector to perform a faster
  // normalization step
  vector<double> col_sum = find_col_sum(matrix);
  vector<double> normalized_matrix_line;
  for (unsigned int i = 0; i < matrix.size(); i++) {

    for (unsigned int j = 0; j < matrix[i].size(); j++) {

      // normalized_score = matrix[i][j]/col_sum[j];
      normalized_matrix_line.emplace_back(matrix[i][j] / col_sum[j] + pseudoc);
    }

    norm_matrix.emplace_back(normalized_matrix_line);
    normalized_matrix_line.clear();
  }
}

// Function which saves into a vector called col_sum all the score column sums
// --> This is made to perform the next Normalization step faster
vector<double> matrix_class::find_col_sum(
    vector<vector<double>> &matrix) { // CHANGE: Possible use of reference
  // PROFILE_FUNCTION();
  vector<double> col_sum;
  double sum = 0;
  for (unsigned int i = 0; i < matrix[0].size(); i++) {
    for (unsigned int j = 0; j < 4; j++) {

      sum += matrix[j][i];
    }

    col_sum.emplace_back(sum);
    sum = 0;
  }

  return col_sum;
}

// Function to perform a second normalization on matrix scores (without a
// pseudocount addition)
void matrix_class::matrix_normalization(
    vector<vector<double>> &matrix) { // CHANGE: Possible use of reference
  // PROFILE_FUNCTION();
  // Calculate and save again the scores column sum into a vector to perform a
  // faster normalization step
  vector<double> col_sum = find_col_sum(matrix);

  for (unsigned int i = 0; i < matrix.size(); i++) {

    for (unsigned int j = 0; j < matrix[i].size(); j++) {

      // Substitution of first normalized values with new normalized ones
      norm_matrix[i][j] = matrix[i][j] / col_sum[j];
    }
  }
}

// Function to calculate, from the normalized matrix, the logarithmic values of
// the scores and creates a new matrix called matrix_log
void matrix_class::matrix_logarithmic(
    vector<vector<double>> &matrix) { // CHANGE: Possible use of reference
  // PROFILE_FUNCTION();
  vector<double> log_matrix_line;
  for (unsigned int i = 0; i < matrix.size(); i++) {

    for (unsigned int j = 0; j < norm_matrix[i].size(); j++) {

      // log_scores = log(norm_matrix[i][j]);
      log_matrix_line.emplace_back(log(norm_matrix[i][j]));
    }

    matrix_log.emplace_back(log_matrix_line);
    log_matrix_line.clear();
  }
}

// Function which return the Transposed matrix from a matrix in input
vector<vector<double>>
matrix_class::reverse_matrix(vector<vector<double>> &matrix) {
  // PROFILE_FUNCTION();
  vector<vector<double>> rev_matrix = matrix;
  reverse(rev_matrix.begin(), rev_matrix.end());

  for (int i = 0; i < 4; i++) {

    reverse(rev_matrix[i].begin(), rev_matrix[i].end());
  }

  return rev_matrix;
}

// Finding the best and the worst score that an oligo can reach based on current
// JASPAR matrix
void oligo_class::find_minmax(vector<vector<double>> &matrix) {
  // PROFILE_FUNCTION();
  vector<double> column;
  vector<double> o_matrix_maxes;


  // Extract the mins and the maxes from each columns, saved into vectors and
  // their total sum will be the best and the worst score that an oligo can
  // reach
  for (unsigned int i = 0; i < matrix[0].size(); i++) {

    for (unsigned int j = 0; j < matrix.size(); j++) {

      column.emplace_back(matrix[j][i]);
    }

    o_matrix_mins.emplace_back(*min_element(column.begin(), column.end()));
    o_matrix_maxes.emplace_back(*max_element(column.begin(), column.end()));
    column.clear();
  }

  min_possible_score =
      accumulate(o_matrix_mins.begin(), o_matrix_mins.end(), 0.0);
  max_possible_score =
      accumulate(o_matrix_maxes.begin(), o_matrix_maxes.end(), 0.0);
}

// Function to calculate the score of a general oligo against a JASPAR matrix
void oligo_class::shifting(vector<vector<double>> &matrix, string &sequence) {
  // PROFILE_FUNCTION();
  unsigned int max = 0;
  max = sequence.size() - matrix[0].size();
  for (unsigned int s_iterator = 0; s_iterator <= max; s_iterator++) {
    double sum_scores = 0;
    // For each oligo in the current sequence a score is calculated
    for (unsigned int i = 0; i < matrix[0].size(); i++) {

      switch (sequence[i + s_iterator]) {

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

      default: // Case if there is N

        sum_scores += o_matrix_mins[i];
        break;
      }
    }

    // The total score of an oligo is saved into an oligo_scores vector
    oligo_scores.emplace_back(sum_scores);
  }
}

// Best score normalization with normalization formula (The parameter to
// normalize have already been calculated and saved into the class)
void oligo_class::scores_normalization() {
  // PROFILE_FUNCTION();

  for (unsigned int i = 0; i < oligo_scores.size(); i++) {

    oligo_scores[i] = 1 + ((oligo_scores[i] - max_possible_score) /
                           (max_possible_score - min_possible_score));
  }
}

// Function to find the best oligo score. From every sequence from GEP the best
// oligo is calculated and both oligo and position in the window are saved
unsigned int oligo_class::find_best_score() {
  // PROFILE_FUNCTION();
  // Extracting the best score from oligo_scores with function max_element
  best_score = *max_element(oligo_scores.begin(), oligo_scores.end());

  vector<int> positions;
  unsigned int matches = 0;

  // Check if there are more than one oligo with the best score --> if any count
  // their numbers and save their position in sequence into a position vector
  for (unsigned int i = 0; i < oligo_scores.size(); i++) {

    if (oligo_scores[i] == best_score) {

      matches = matches + 1;
      positions.emplace_back(i);
    }
  }

  // If more than one oligo carries the best score calculate distances and
  // select the nearest to the center
  if (matches > 1) {
    vector<unsigned int> dist_center;
    unsigned int min_distance = 0;
    vector<unsigned int>::iterator itr;
    for (int &p : positions) {

      unsigned int distance;

      // The distance from the center need to be calculated as an absolute value
      // (no negative values)
      distance = abs(p - half_length);
      dist_center.emplace_back(distance);
    }

    // Once distances have been calculated select as best the one who is the
    // nearest to the sequence center --> select the min element from distance
    // vector
    min_distance = *min_element(dist_center.begin(), dist_center.end());

    // Find min_distance on dist_centre vector and save its index
    itr = find(dist_center.begin(), dist_center.end(), min_distance);
    unsigned int index = distance(dist_center.begin(), itr);

    // Index of min_distance on distance vector == the corresponding position in
    // the positions vector
    return positions[index];
  }

  // If just one best oligo score has been found return its position (it is the
  // first and the only element in positions vector)
  return positions[0];
}


// Function to read BED and 2Bit files and create GEP (vector of bed class)
void coordinator_class::GEP_creation(vector<bed_class> &GEP, char ** result) {
  // PROFILE_FUNCTION();
  // RAM_usage();
  cout << "\n- [1] Extract bed coordinate sequences from reference genome  \n";

  ifstream in(BED_FILE);

  string line;

  // Line counter initialization
  unsigned int n_line = 1;

  // For each line in BED file create a bed class called bed_line in which chr,
  // start and end coordinate are ridden, centered and saved Then the Fasta
  // sequence is extracted from Twobit genome following the coordinated and
  // saved into a string variable
  while (getline(in, line)) {

    // if line is empty or commented --> continue
    if (line.empty() || line[0] == '#') {

      continue;
    }

    bed_class bed_line(line, tb, n_line, result);

    // For each line a bed class is created --> All the bed classes are saved in
    // GEP vector (vector og bed class)
    GEP.emplace_back(bed_line);

    n_line = n_line + 1;
  }
  if (err == true) {
    exit(1);
  }
}

// Function to read JASPAR PWM file, extract values and create a matrix class
vector<vector<double>> coordinator_class::read_JASPAR() {
  // PROFILE_FUNCTION();
  // RAM_usage();
  cout << "- [2] Reading JASPAR MATRIX file and extracting values\n";

  ifstream file(JASPAR_FILE);
  string line;

  // For each line of the JASPAR file
  while (getline(file, line)) {

    // If the line start with '>' character save the words into matrix_name
    // string and into tf_name string
    if (line[0] == '>') {

      istringstream mystream(line);
      // mystream >> matrix_name >> tf_name;
    }

    // If the line does not start with '>' delete the '[',']' and 'A,T,C,G'
    // characters, extract the scores and save them into a scores matrix
    else {

      // Deleting from line the first character (A,T,C,G), the '[' and the ']'
      // characters
      line.erase(0, line.find('[') + 1);
      line.erase(line.find(']'));

      vector<double> scores_line;
      istringstream mystream(line);

      // For each words (number) in line put the current number in num variables
      for (double num; mystream >> num;) {

        scores_line.emplace_back(num);
      }

      matrix.emplace_back(scores_line);
    }
  }

  file.close();
  // RAM_usage();
  // Cout of step 3/4 here because the normalization and reverse function will
  // be re-utilized during the workflow
  cout << "- [3] Jaspar Matrix normalization\n";
  cout << "- [4] Jaspar Matrix reverse complement determination to analize the "
          "reverse strand\n";

  return matrix;
}

// Function to create oligos_vector (a oligo class vector)
void coordinator_class::oligos_vector_creation(
    vector<oligo_class> &oligos_vector, vector<vector<double>> &matrix_log,
    vector<vector<double>> &matrix_log_inverse, vector<bed_class> &GEP) {
  // PROFILE_FUNCTION();

  cout << "- [5] Analyzing sequences using Jaspar matrix\n";

  // For every sequences into GEP vector
  for (unsigned int i = 0; i < GEP.size(); i++) {

    // Calling the oligo_class constructor to analyze the shifting of the
    // sequence on log_matrix (FWD strand analysis)
    oligo_class SHIFTING(matrix_log, GEP[i].sequence, GEP[i].start_coord);

    // The oligo class just created is saved into oligos_vector (oligo_class
    // vector)
    oligos_vector.emplace_back(SHIFTING);

    // If the analysis is on DS calling the oligo_class constructor to analyze
    // the shifting of sequence on inverse_log_matrix (REVERSE strande analysis)
    if (DS) {

      oligo_class SHIFTING(matrix_log_inverse, GEP[i].sequence, GEP[i].start_coord);
      oligos_vector.emplace_back(SHIFTING);
    }
  }
  // RAM_usage();
  cout << "- [6] Selecting the best Jaspar's oligo for each sequence \n";
}

// Function useful, if the analysis is performed on DS, to choose from the best
// FWD strand oligo and the best REV strand oligo the best one to keep as "Best
// oligo" --> The oligos vector is divided in half and only the best strand for
// each sequence is kept
void coordinator_class::best_strand() {
  // PROFILE_FUNCTION();
  vector<oligo_class> comparison;

  for (unsigned int i = 0; i < oligos_vector.size(); i += 2) {

    // The comparison is made by oligo_class in i position against the oligo
    // class in i+1 position (The fwd and rev strand of the same sequence, which
    // are consecutive into the oligos_vector)
    double best_score_norm_positive = oligos_vector[i].best_score;
    double best_score_norm_negative = oligos_vector[i + 1].best_score;

    if (best_score_norm_positive >= best_score_norm_negative) {

      comparison.emplace_back(oligos_vector[i]);
      rev.emplace_back(false);
    }

    else {
      rev.emplace_back(true);
      comparison.emplace_back(oligos_vector[i + 1]);
    }
  }

  // The new oligos_vector is replaced by comparison vector, which contains only
  // the best strand
  oligos_vector.clear();

  oligos_vector = comparison;
}

// Function to re-set the genomic coordinates and the sequences window -->
// centered on the best oligo found for each sequence
void coordinator_class::centering_oligo() {
  // PROFILE_FUNCTION();

  int center_oligo;

  // To center on the best oligo, centering_function and extract_seq functions
  // from bed_class need to be recalled with updated input parameters
  for (unsigned int i = 0; i < oligos_vector.size(); i++) {

    // The center of the window is exactly on the center of the best oligo
    // (which length depends to the JASPAR matrix size)
    center_oligo =
        oligos_vector[i].start_coord_oligo + matrix_log[0].size() / 2;
    GEP[i].centering_function(center_oligo, center_oligo, half_length, 0);
    GEP[i].extract_seq(tb, 0);
    if(DS){
      if (rev[i] && direction) {
        check_palindrome(GEP[i].sequence, reverse_bases);
        GEP[i].sequence = reverse_bases;
        if (matrix_log[0].size() % 2 == 0) {
          center_oligo =
              (oligos_vector[i].start_coord_oligo + matrix_log[0].size() / 2);
          GEP[i].centering_function(center_oligo, center_oligo, half_length, 0);
        } else {
          center_oligo =
              (oligos_vector[i].start_coord_oligo + matrix_log[0].size() / 2) + 1;
          GEP[i].centering_function(center_oligo, center_oligo, half_length, 0);
        }
        GEP[i].extract_seq(tb, 0);
        check_palindrome(GEP[i].sequence, reverse_bases);
        GEP[i].sequence = reverse_bases;
      }
    }
  }
}

// Function able to convert a string (containing numbers separated by ",") into
// a vector of unsigned int
vector<unsigned int> generic_vector_creation(string numbers) {
  // PROFILE_FUNCTION();
  int index = 0;
  vector<unsigned int> vec;

  // When index is == -1 means that it is pointing to the end character of the
  // string
  while (index != -1) {

    // Find the first "," character into the string
    index = numbers.find(",");

    // Put everything before "," into a vector
    vec.emplace_back(stoi(numbers.substr(0, index)));

    // Erase from the string the number already inserted and the ","
    numbers.erase(0, index + 1);
  }

  return vec;
}

// Function to extract fasta sequences from a Multifasta file and save them into
// a vector of string
void multifasta_class::extract_sequences() {
  // PROFILE_FUNCTION();
  cout << "\n- [1] Extracting sequences from MultiFasta file \n";

  ifstream file(MFASTA_FILE);
  string line;
  string current_sequence;

  // First line is flagged to 1 to ignore it (because in MF format sometimes
  // there is an header first line)
  bool first_line = true;

  // The while cycle follows a complex reasoning:
  // if in line there is the first line --> all ignored --> then first line flag
  // set to 0 if in line there is the FASTA sequence --> line is saved in
  // current_sequence variable --> then pass to next line if in line there is
  // the header (but not the first line) --> current_sequence previously filled
  // by FASTA is saved in sequences vector (a vector of strings) --> then
  // current_sequence is clean
  while (getline(file, line)) {

    if (line[0] == '>' && !first_line) {

      sequences.emplace_back(current_sequence);
      current_sequence.clear();

    }

    else if (!first_line) {

      if (line[0] != ' ' && line.size() != 0) {

        // Before to save the sequence in current_sequence variable the FASTA
        // characters are capitalized
        transform(line.begin(), line.end(), line.begin(), ::toupper);
        current_sequence = current_sequence + line;
      }
    }

    // After the first cicle the first line flag is set to 0
    first_line = false;
  }
  sequences.emplace_back(current_sequence);
}

// Controlling that all the MF sequences have the same langth
void multifasta_class::length_control(vector<string> sequences) {
  // PROFILE_FUNCTION();
  cout << "- [2] Multifasta Sequences length check\n";

  // unsigned int size = sequences[0].size();

  for (unsigned int i = 0; i < sequences.size(); i++) {

    // If only one sequence in the vector are longer an error is generated
    if (sequences[i].size() != sequences[0].size()) {

      cerr << "Sequences are not of the same length!\n";
      exit(1);
    }
  }
}

// Creating alias file following the multifasta input filename
void multifasta_class::alias_output_filename() {
  // PROFILE_FUNCTION();
  alias_file =
      alias_file + MFASTA_FILE.erase(0, MFASTA_FILE.find_last_of("_") + 1);
  alias_file =
      MFASTA_FILE.erase(MFASTA_FILE.find_last_of("."), MFASTA_FILE.size()) +
      "_";
}

void multifasta_class::GEP_creation_MF(vector<string> sequences) {
  // PROFILE_FUNCTION();
  cout << "- [3] Sorting Multifasta sequences\n";

  for (unsigned int i = 0; i < sequences.size(); i++) {

    // A bed class is created for each sequence
    bed_class BED_MULTIFASTA(sequences[i]);

    // The classes are stored into a GEP vector
    GEP.emplace_back(BED_MULTIFASTA);
  }
}

// Function used to order the elements of a class (in this case it orders PValueClass elements)
bool comp(const PvalueClass &P1, const PvalueClass &P2) {
  // PROFILE_FUNCTION();
  return (P1.pvalue == P2.pvalue) ? (P1.K > P2.K) : (P1.pvalue < P2.pvalue);
}

bool comp_occ(const PvalueClass &P1, const PvalueClass &P2) {
  // PROFILE_FUNCTION();
  return (P1.K == P2.K) ? ((P1.N1 == P2.N1)? (P1.oligo < P2.oligo) : (P1.N1 < P2.N1)) : (P1.K > P2.K);
}

void HammingClass::PWMHammingCalc() {
  // PROFILE_FUNCTION();
  // Vector for each bases initialized to count, position by position, the
  // occurrences of that base
  vector<double> vec_A, vec_C, vec_G, vec_T;


  // For each bases (position) of oligo
  for (unsigned int i = 0; i < hamming_seed[0].size();i++){

    double counter_A = 0;
    double counter_C = 0;
    double counter_G = 0;
    double counter_T = 0;

    // For each oligo in similar oligos vector
    for (unsigned int j = 0; j < hamming_seed.size(); j++) {

      cluster_map.insert(pair<string, double>(
        hamming_seed[j], vert_vector[j]));

      switch (hamming_seed[j][i]) {

      // Increment base counters of the oligo occurrences
      case 'A':
        counter_A += vert_vector[j];
        break;
      case 'C':
        counter_C += vert_vector[j];
        break;
      case 'G':
        counter_G += vert_vector[j];
        break;
      case 'T':
        counter_T += vert_vector[j];
        break;
      }
    }

    // Fill base vectors
    vec_A.emplace_back(counter_A);
    vec_C.emplace_back(counter_C);
    vec_G.emplace_back(counter_G);
    vec_T.emplace_back(counter_T);
  }

  // Build the PWM matrix
  PWM_hamming.emplace_back(vec_A);
  PWM_hamming.emplace_back(vec_C);
  PWM_hamming.emplace_back(vec_G);
  PWM_hamming.emplace_back(vec_T);
}

void HammingClass::DPWMHamming(vector<vector<double>> &PWM_hamming){
  // PROFILE_FUNCTION();
  for (unsigned short int i = 0; i < PWM_hamming.size(); i++) {
    for (unsigned short int j = 0; j < PWM_hamming[i].size(); j++) {
      cout << PWM_hamming[i][j] << "\t";
    }
    cout << endl;
  }
  cout << endl;
}

void HammingClass::CheckSeed(string seed,
                             unordered_map<string, VerticalClass> &map_vertical,
                             multimap<int, string, greater<int>> &pos,
                             unsigned int position, unsigned int d) {
  // PROFILE_FUNCTION();
  for (multimap<int, string>::iterator it = pos.begin(); it != pos.end();
       it++) {
    hamming_v_occ = 0;
    tot_freq += it->first;
    string oligo = it->second;
    if (oligo == seed){
      seed_vertical = it->first;
    }
    unsigned int i = 0, count = 0;
    // Here the hamming distance is calculated
    while (seed[i] != '\0') {
      if (seed[i] != oligo[i]){
        count++;
      }
      i++;
    }
    // If the hamming distance is less or equal the distance by the user 
    // the oligo is added to the cluster and all the occurrences of the oligos present in the 
    // cluster are set to zero (for secondary matrix analysis)
    if (count <= distance_vector[d]) {

      string reverse_o = reverse_oligo(oligo);
      unordered_map<string, VerticalClass>::iterator it_ver = map_vertical.find(oligo);
      unordered_map<string, VerticalClass>::iterator it_ver_rc = map_vertical.find(reverse_o);
      if (it_ver != map_vertical.end()) {
        hamming_v_occ += it_ver->second.vertical_count_FWD[position];
      } else {
        hamming_v_occ += it_ver_rc->second.vertical_count_rc_FWD[position];
      }        
      if(DS){
        if (it_ver != map_vertical.end()) {
          hamming_v_occ += it_ver->second.vertical_count_REV[position];
        } else {
          hamming_v_occ += it_ver_rc->second.vertical_count_FWD[position];
        } 
      }
      vert_vector.push_back(hamming_v_occ);
      hamming_seed.push_back(oligo);
    }
  }
}

void HammingClass::Freq1Calc() {
  // PROFILE_FUNCTION();
  freq1 = static_cast<double>(seed_vertical) / static_cast<double>(tot_freq);
  cout << "FREQ: " << freq1 << endl;
  // double freq_right = static_cast<double>(tot_freq)/static_cast<double>(frequence);
  // cout << frequence << endl;
  // cout << "OLD FREQ: " << freq_right << endl;
}

void HammingClass::HoccCalc(unordered_map<string, HorizontalClass> &map_horizontal) {
 // PROFILE_FUNCTION();
  hamming_H_occ = 0;
  for (unsigned int i = 0; i < hamming_seed.size(); i++) {
    string oligo = hamming_seed[i];
    string reverse_o = reverse_oligo(oligo);
    unordered_map<string, HorizontalClass>::iterator it = map_horizontal.find(oligo);
    unordered_map<string, HorizontalClass>::iterator itrc = map_horizontal.find(reverse_o);
    if (it != map_horizontal.end()) {
      hamming_H_occ += it->second.horizontal_count;
    } else {
      hamming_H_occ += itrc->second.horizontal_count_rc;
    }
  }
}

void HammingClass::Freq2Calc() {
  // PROFILE_FUNCTION();
  freq2 = static_cast<double>(seed_vertical) / static_cast<double>(hamming_H_occ);
}

void HammingClass::ClearVertical(multimap<int, string, greater<int>> &pos,
                                  unordered_map<string, VerticalClass> &map_vertical, unsigned int j){

  for(unsigned int i = 0; i < hamming_seed.size(); i++){
    
    for(multimap<int, string, greater<int>>::iterator it = pos.begin();
      it != pos.end(); ++it){
      
      string oligo = it->second;
    
      if(oligo == hamming_seed[i]){
        pos.erase(it);
        break;
      }
    }
  }
}

void DVector(vector<PvalueClass> &P_vector, unsigned int j) {
  // PROFILE_FUNCTION();
  for (unsigned int c = 0; c < P_vector.size() && c < 10; c++) {
    cout << j + 1 << "\t";
    cout << c + 1 << "\t";
    cout << P_vector[c].oligo << "\t";
    cout << P_vector[c].K << "\t";
    cout << P_vector[c].N1 << "\t";
    cout << P_vector[c].N2 << "\t";
    cout << P_vector[c].pvalue << "\t";
    cout << log10(P_vector[c].pvalue) * -1 << endl;
  }
}

// Function where K, N1, N2 and T are calculated in order to obtain the p value
void PvalueClass::TKN1Calc(vector<bed_class> &GEP,
                           multimap<int, string>::iterator &it,
                           unordered_map<string, HorizontalClass> &vector_map_hor,
                           unsigned int i) {
  // PROFILE_FUNCTION();
  // T is the number of sequences
  unsigned int T = GEP.size();

  // For each oligo in the multimap of vertical occurrences T, N1, N2 and K are
  // calculated.

  // Remember N1 is the number of horizontal occurrences of oligo, T is the
  // total number of sequences, N2 is the total number of oligos in all
  // sequences minus N1 and K is the vertical occurrences of the oligo.

  N1 = 0;
  K = it->first;
  oligo = it->second;

  unordered_map<string, HorizontalClass>::iterator itBigMap =
      vector_map_hor.find(oligo);
  unordered_map<string, HorizontalClass>::iterator itBigMap_rc =
      vector_map_hor.find(reverse_oligo(oligo));
  if (itBigMap == vector_map_hor.end()) {

    if (!itBigMap_rc->second.palindrome) {
      N1 = itBigMap_rc->second.horizontal_count_rc;
    }
  } else {

    N1 = itBigMap->second.horizontal_count;
  }
  // Calculation of total number of oligos in the multifasta
  tot_oligos = T * len[i];

  // Calculation of N2
  N2 = tot_oligos - N1;

  // Using the gsl library for the hypergeometric p_value
  pvalue = gsl_cdf_hypergeometric_Q(K, N1, N2, T);
  if (pvalue == 0) {
    pvalue = 1e-300;
  }
}

// Function for counting the horizontal and vertical occurrences for each oligo
// and putting them inside a map.
void MapClass::CountOccurrencesHor(string sequence, int k) {
  // PROFILE_FUNCTION();
  for (unsigned int i = 0; i < (sequence.size() - k + 1); i++) {
    string oligo = sequence.substr(i, k);
    string oligo_rc = reverse_oligo(oligo);
    unordered_map<string, HorizontalClass>::iterator it = horizontal_map.find(oligo);
    unordered_map<string, HorizontalClass>::iterator it_rc = horizontal_map.find(oligo_rc);
    if (it != horizontal_map.end()) {
      it->second.horizontal_count++;
      it->second.horizontal_count_FWD++;
      if (DS) {
        if(it->second.palindrome){
          it->second.horizontal_count_REV++;
          it->second.horizontal_count_rc_FWD++;
          it->second.horizontal_count++;
          it->second.horizontal_count_rc++;
        }
        it->second.horizontal_count_rc++;
        it->second.horizontal_count_rc_REV++;
      }
    } else if (it_rc != horizontal_map.end()) {
      if (!it_rc->second.palindrome) {
        it_rc->second.horizontal_count_rc++;
        it_rc->second.horizontal_count_rc_FWD++;
        if (DS) {
          it_rc->second.horizontal_count_REV++;
          it_rc->second.horizontal_count++;
        }
      }

    } else {
      HorizontalClass Hor;
      Hor.oligo = oligo, Hor.oligo_rc = oligo_rc;
      Hor.palindrome = bool(oligo == oligo_rc);
      Hor.horizontal_count = 1, Hor.horizontal_count_rc = 0;
      Hor.horizontal_count_FWD = 1, Hor.horizontal_count_rc_FWD = 0;
      Hor.horizontal_count_REV = 0, Hor.horizontal_count_rc_REV = 0;
      if (DS) {
        Hor.horizontal_count_FWD = 1, Hor.horizontal_count_rc_FWD = 0;
        Hor.horizontal_count_REV = 0, Hor.horizontal_count_rc_REV = 1;
        Hor.horizontal_count_rc = 1, Hor.horizontal_count = 1;
        if(oligo == oligo_rc){
          Hor.horizontal_count_rc_FWD = 1;
          Hor.horizontal_count_REV = 1;
          Hor.horizontal_count = 2, Hor.horizontal_count_rc = 2;
        }
      }
    
      horizontal_map.emplace(oligo, Hor);
    }
  }
}

void MapClass::CountOccurrencesVer(string sequence, int k) {
  // PROFILE_FUNCTION();
  for (unsigned int i = 0; i < (sequence.size() - k + 1); i++) {
    // unsigned int tot_freq = 0;
    string oligo = sequence.substr(i, k);
    string oligo_rc = reverse_oligo(oligo);
    unordered_map<string, VerticalClass>::iterator it = vertical_map.find(oligo);
    unordered_map<string, VerticalClass>::iterator it_rc = vertical_map.find(oligo_rc);
    if (it != vertical_map.end()) {
      it->second.vertical_count[i]++;
      it->second.vertical_count_FWD[i]++;
      it->second.vertical_count_rc_REV[i]++;
      it->second.vertical_count_rc[i]++;
      // tot_freq++;
      if (DS) {
        if(it->second.palindrome){         
          it->second.vertical_count[i]++;
          it->second.vertical_count_rc[i]++;
          it->second.vertical_count_REV[i]++;
          it->second.vertical_count_rc_FWD[i]++;
        }
      //   tot_freq++;
        
      }
    } else if (it_rc != vertical_map.end()) {
      if (!it_rc->second.palindrome) {
        it_rc->second.vertical_count_rc[i]++;
        it_rc->second.vertical_count_rc_FWD[i]++;
        if (DS) {
          it_rc->second.vertical_count[i]++;
          it_rc->second.vertical_count_REV[i]++;
        }
      }

    } else {
      VerticalClass Ver;
      Ver.oligo = oligo, Ver.oligo_rc = oligo_rc;
      Ver.palindrome = bool(oligo == oligo_rc);
      Ver.vertical_count.resize(((sequence.size() - k) + 1), 0);
      Ver.vertical_count_rc.resize(((sequence.size() - k) + 1), 0);
      Ver.vertical_count_FWD.resize(((sequence.size() - k) + 1), 0);
      Ver.vertical_count_rc_FWD.resize(((sequence.size() - k) + 1), 0);
      Ver.vertical_count_REV.resize(((sequence.size() - k) + 1), 0);
      Ver.vertical_count_rc_REV.resize(((sequence.size() - k) + 1), 0);
      Ver.vertical_count_FWD[i] = 1;
      Ver.vertical_count[i] = 1;
      if (DS) {
        Ver.vertical_count_rc[i] = 1;
        Ver.vertical_count_FWD[i] = 1, Ver.vertical_count_rc_REV[i] = 1;
        Ver.vertical_count_REV[i] = 0, Ver.vertical_count_rc_FWD[i] = 0;
        if(oligo == oligo_rc){
          Ver.vertical_count_rc_FWD[i] = 1;
          Ver.vertical_count_REV[i] = 1;
          Ver.vertical_count[i] = 2;
          Ver.vertical_count_rc[i] = 2;
        }
      }
    
      vertical_map.emplace(oligo, Ver);
    }
  }
}

string reverse_oligo(string bases) {
  // PROFILE_FUNCTION();
  string reverse_string;
  for (int i = bases.length() - 1; i >= 0; i--) {

    switch (bases[i]) {

    case 'A':
      reverse_string.append("T");
      break;
    case 'T':
      reverse_string.append("A");
      break;
    case 'G':
      reverse_string.append("C");
      break;
    case 'C':
      reverse_string.append("G");
      break;
    case 'N':
      reverse_string.append("N");
      break;
    }
  }
  return reverse_string;
}

void MapClass::MainMapVector(vector<bed_class> &GEP) {
  // PROFILE_FUNCTION();
  for (unsigned int i = 0; i < kmers_vector.size(); i++) {
    for (unsigned int j = 0; j < GEP.size(); j++) {
      CountOccurrencesHor(GEP[j].sequence, kmers_vector[i]);
      CountOccurrencesVer(GEP[j].sequence, kmers_vector[i]);
    }
    vector_map_hor.push_back(horizontal_map);
    horizontal_map.clear();
    vector_map_ver.push_back(vertical_map);
    vertical_map.clear();
  }
}

void MapClass::VerticalMapVector() {
  // PROFILE_FUNCTION();
  for (unsigned int i = 0; i < kmers_vector.size(); i++) {
    for (unsigned int j = 0; j < len[i]; j++) {
      multimap<int, string, greater<int>> pos;
      for (unordered_map<string, VerticalClass>::iterator it =
               vector_map_ver[i].begin();
           it != vector_map_ver[i].end(); it++) {
        if (it->second.vertical_count[j] > 0) {
          pos.emplace(it->second.vertical_count[j], it->first);
        }
        if (!it->second.palindrome && it->second.vertical_count_rc[j] > 0) {
          pos.emplace(it->second.vertical_count_rc[j], it->second.oligo_rc);
        }
      }
      positions_occurrences.push_back(pos);
      pos.clear();
    }
    vector_positions_occurrences.push_back(positions_occurrences);
    positions_occurrences.clear();
  }
}

// Debug function
void MapClass::DVerticalMapVector() {
  // PROFILE_FUNCTION();
  
  for (unsigned int i = 0; i < kmers_vector.size(); i++) {
    cout << "I: " << i << endl;
    for (unsigned int j = 0; j < vector_positions_occurrences[i].size(); j++) {
      cout << "J: " << j << endl;
      int sum = 0;
      for (multimap<int, string>::iterator it =
               vector_positions_occurrences[i][j].begin();
           it != vector_positions_occurrences[i][j].end(); it++) {
        cout << it->first << "   " << it->second << endl;
        sum = sum + it->first;
      }
      cout << "SOMMA TOTALE VERTICAL: " << sum << endl << endl;
    }
  }
}

// Debug function
void MapClass::DMainMapVectorDS() {
  // PROFILE_FUNCTION();
  for (unsigned int i = 0; i < kmers_vector.size(); i++) {
    for (unordered_map<string, HorizontalClass>::iterator it =
             vector_map_hor[i].begin();
         it != vector_map_hor[i].end(); it++) {
      cout << it->first << " " << it->second.horizontal_count_FWD << " " 
           << it->second.horizontal_count_REV << " " << it->second.horizontal_count
           << "\t" << it->second.oligo_rc << " " << it->second.horizontal_count_rc_FWD
           << " " << it->second.horizontal_count_rc_REV << " "
           << it->second.horizontal_count_rc << " "
           << " " << boolalpha << it->second.palindrome << noboolalpha << endl;
    }
    cout << "\n\n VERTICAL MAP: \n\n" << endl;
    for (unordered_map<string, VerticalClass>::iterator it =
             vector_map_ver[i].begin();
         it != vector_map_ver[i].end(); it++) {
      cout << it->first << " " << it->second.vertical_count_FWD[21] << " " 
           << it->second.vertical_count_REV[21] << " " << it->second.vertical_count[21]
           << "\t" << it->second.oligo_rc << " " << it->second.vertical_count_rc_FWD[21]
           << " " << it->second.vertical_count_rc_REV[21] << " "
           << it->second.vertical_count_rc[21] << " "
           << " " << boolalpha << it->second.palindrome << noboolalpha << endl;
    }
  }
}

void MapClass::DMainMapVectorSS() {
  // PROFILE_FUNCTION();
  for (unsigned int i = 0; i < kmers_vector.size(); i++) {
    for (unordered_map<string, HorizontalClass>::iterator it =
             vector_map_hor[i].begin();
         it != vector_map_hor[i].end(); it++) {
      cout << it->first << " " << it->second.horizontal_count << "\t" 
           << it->second.oligo_rc << " " << it->second.horizontal_count_rc << " "
           << " " << boolalpha << it->second.palindrome << noboolalpha << endl;
    }
        for (unordered_map<string, VerticalClass>::iterator it =
             vector_map_ver[i].begin();
         it != vector_map_ver[i].end(); it++) {
      cout << it->first << " " << it->second.vertical_count_FWD[21] << " " 
           << it->second.vertical_count_REV[21] << " " << it->second.vertical_count[21]
           << "\t" << it->second.oligo_rc << " " << it->second.vertical_count_rc_FWD[21]
           << " " << it->second.vertical_count_rc_REV[21] << " "
           << it->second.vertical_count_rc[21] << " "
           << " " << boolalpha << it->second.palindrome << noboolalpha << endl;
    }
  }
}

void EMClass::print_PWM(string name, vector<vector<double>> &PWM_hamming) {
  // PROFILE_FUNCTION();
  cout << name << endl;
  for (unsigned short int i = 0; i < PWM_hamming.size(); i++) {
    for (unsigned short int j = 0; j < PWM_hamming[i].size(); j++) {
      cout << PWM_hamming[i][j] << "\t";
    }
    cout << endl;
  }
  cout << endl;
}

void EMClass::EM_Ipwm(vector<vector<double>> &PWM_hamming) {
  // PROFILE_FUNCTION();

  // matrix::normalize()

  double sum = 0;
  double corr = 0;

  for (unsigned int j = 0; j < PWM_hamming.size(); j++) {
    sum = sum + PWM_hamming[j][0];
  }

  corr = sqrt(sum);

  for (unsigned int x = 0; x < PWM_hamming.size(); x++) {
    for (unsigned int y = 0; y < PWM_hamming[0].size(); y++) {
      PWM_hamming[x][y] = PWM_hamming[x][y] + corr;
    }
  }

  sum = 0;
  for (unsigned int j = 0; j < PWM_hamming.size(); j++) {
    sum += PWM_hamming[j][0];
  }
  for (unsigned short int i = 0; i < PWM_hamming.size(); i++) {
    for (unsigned short int j = 0; j < PWM_hamming[i].size(); j++) {
      PWM_hamming[i][j] = PWM_hamming[i][j] / sum;
    }
  }
}

/*
This function is about the expectation step of the expectation-maximization
algorithm where we obtain the likelihood ratio for each oligo in the vertical
map
*/
void EMClass::EM_Epart(map<string,double> &cluster_map, vector<vector<double>> &PWM_hamming,
    unordered_map<string, HorizontalClass> &map_horizontal) {
  // PROFILE_FUNCTION();
  
  unsigned int sum_hor = 0;

  for (map<string, double>::iterator it = cluster_map.begin();
       it != cluster_map.end(); it++) {
    unordered_map<string, HorizontalClass>::iterator occ_oligo_it =
        map_horizontal.begin();
    unordered_map<string, HorizontalClass>::iterator occ_oligo_it_rev =
        map_horizontal.begin();
    occ_oligo_it = map_horizontal.find(it->first);
    occ_oligo_it_rev = map_horizontal.find(reverse_oligo(it->first));

    if (occ_oligo_it != map_horizontal.end()) {
      sum_hor = sum_hor + occ_oligo_it->second.horizontal_count;
    }
    if (occ_oligo_it_rev != map_horizontal.end()) {
      sum_hor = sum_hor + occ_oligo_it_rev->second.horizontal_count_rc;
    }
  }

  vector<double> LR;
  vector<string> oligo;

  // In this cycle for each element in vertical map we calculate the probability
  // that this oligo is present in the hamming matrix
  for (map<string, double>::iterator it = cluster_map.begin();
       it != cluster_map.end(); it++) {
    unordered_map<string, HorizontalClass>::iterator occ_oligo_it =
        map_horizontal.begin();
    unordered_map<string, HorizontalClass>::iterator occ_oligo_it_rev =
        map_horizontal.begin();
    string similar_oligo = it->first;
    double P_oligo = 1;
    double P_bg = 0;
    double horizontal_occurences = 0;
    double likelihood_ratio = 0;

    occ_oligo_it = map_horizontal.find(it->first);
    occ_oligo_it_rev = map_horizontal.find(reverse_oligo(it->first));

    if (occ_oligo_it != map_horizontal.end()) {
      horizontal_occurences = horizontal_occurences + occ_oligo_it->second.horizontal_count;
    }
    if (occ_oligo_it_rev != map_horizontal.end()) {
      horizontal_occurences = horizontal_occurences + occ_oligo_it_rev->second.horizontal_count_rc;
    }

    P_bg = horizontal_occurences / sum_hor;

    for (unsigned int k = 0; k < PWM_hamming[0].size(); k++) {
      switch (it->first[k]) {

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

      default: // Case if there is N
        P_oligo *= 1;
        break;
      }
    }

    likelihood_ratio = P_oligo / P_bg;
    LR.emplace_back(likelihood_ratio);
    oligo.emplace_back(similar_oligo);

    /*
     *The like_ratio_map is a map where for each oligo present in the vertical
     *map we couple the likelihood ratio previously calculated with the ratio
     *between the probability to have the oligo in the PWM_hamming and the
     *background probability
     */
  }

  double sum = 0;
  double Nsites = ceil((cluster_map.size() / 2) + 1);

  for (unsigned int i = 0; i < cluster_map.size(); i++) {
    sum += LR[i];
  }

  for (unsigned int i = 0; i < cluster_map.size(); i++) {
    LR[i] = LR[i] / sum;
    LR[i] = LR[i] * Nsites;
    // cout << "Oligo: "<< oligo[i] << "\tLR: " << LR[i] << endl;
    like_ratio_map.insert(pair<string, double>(oligo[i], LR[i]));
  }

  // matrix::squash()

  double total = 0;

  bool renorm = true;

  for (map<string, double>::iterator it = like_ratio_map.begin();
       it != like_ratio_map.end(); ++it) {
    total += it->second;
  }

  while (renorm) {
    renorm = false;

    double norm = total / Nsites;
    total = 0;

    for (map<string, double>::iterator it = like_ratio_map.begin();
         it != like_ratio_map.end(); ++it) {
      double p = it->second;
      if (p < 1) {
        p /= norm;
      }

      if (p > 1) {
        p = 1;
        Nsites--;
        renorm = true;
      }

      it->second = p;

      if (p <= 1) {
        total += p;
      }
    }
  }

  LR.clear();
  oligo.clear();
}

// In this function there is the second part of the EM algorithm and this is the
// maximization part
void EMClass::EM_Mpart(vector<vector<double>> &PWM_hamming) {
  // PROFILE_FUNCTION();

  // matrix::get_p()

  for (unsigned int x = 0; x < PWM_hamming.size(); x++) {
    for (unsigned int y = 0; y < PWM_hamming[0].size(); y++) {
      PWM_hamming[x][y] = 0;
    }
  }
  for (unsigned int b = 0; b < PWM_hamming[0].size(); b++) {
    for (map<string, double>::iterator it = like_ratio_map.begin();
         it != like_ratio_map.end(); ++it) {
      switch (it->first[b]) {
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

      default: // Case if there is N
        break;
      }
    }
  }

  double sum = 0;

  for (unsigned int j = 0; j < PWM_hamming.size(); j++) {
    sum += PWM_hamming[j][0];
  }

  for (unsigned short int i = 0; i < PWM_hamming.size(); i++) {
    for (unsigned short int j = 0; j < PWM_hamming[i].size(); j++) {
      PWM_hamming[i][j] = PWM_hamming[i][j] / sum;
    }
  }
}
// This function is made to check if the EM_cycle reaches convergence
bool EMClass::EM_convergence(vector<vector<double>> &PWM_old,
                                   vector<vector<double>> &PWM_hamming,
                                   bool conv) {
  // PROFILE_FUNCTION();
  vector<vector<double>> PWM_old_conv;
  vector<vector<double>> PWM_hamming_conv;

  PWM_old_conv = PWM_old;
  PWM_hamming_conv = PWM_hamming;

  conv = false;
  for (unsigned int i = 0; i < PWM_hamming.size(); i++) {
    for (unsigned int j = 0; j < PWM_hamming[0].size(); j++) {
      if (abs(PWM_old_conv[i][j] - PWM_hamming_conv[i][j]) < sim_tresh) {
        conv = true;
        break;
      }
    }
  }

  return conv;
}

void EMClass::EM_cycle(map<string,double> &cluster_map,
    vector<vector<double>> &PWM_hamming,
    unordered_map<string, HorizontalClass> &map_horizontal) {
  // PROFILE_FUNCTION();
  bool conv = true;
  int i = 0;

  vector<vector<double>> PWM_old;

  // In this cycle we repeat the EM until the convergence is reached
  if (exp_max == "c") {

    while (conv && i < 200) {

      PWM_old = PWM_hamming;

      EM_Epart(cluster_map, PWM_hamming, map_horizontal);
      EM_Mpart(PWM_hamming);

      like_ratio_map.clear();

      conv = EM_convergence(PWM_old, PWM_hamming, conv);
      PWM_old = PWM_hamming;
      i++;
    }

  } else {

    double em = stod(exp_max);
    for (unsigned int i = 0; i < em; i++) {
      EM_Epart(cluster_map, PWM_hamming, map_horizontal);
      EM_Mpart(PWM_hamming);

      like_ratio_map.clear();
      
    }
  }
}

// Shifting the PWM_matrix on the sequences and calculate local scores (from
// positon where the matrix has been generated), and the global scores from each
// sequences positions
void z_test_class::oligos_vector_creation_PWM(vector<bed_class> &GEP) {
  // PROFILE_FUNCTION();
  // For every sequence
  for (unsigned int i = 0; i < GEP.size(); i++) {

    // Calling oligo class to accede to all functions useful to shift a matrix
    // on sequences --> Shifting on FWD strand
    oligo_class SHIFTING_PWM(matrix_log, GEP[i].sequence);

    // Return oligo scores calculated from previous shifting
    oligo_scores_horizontal_FWD = SHIFTING_PWM.oligo_scores;

    // If analysis is in Double strand
    if (DS) {

      // Make the shifting also on reverse strand, putting as input the
      // inverse_log_matrix --> Shifting on REV strand
      oligo_class SHIFTING_PWM_2(inverse_matrix_log, GEP[i].sequence);

      // Retrun oligo scores from previous shifting
      oligo_scores_horizontal_REV = SHIFTING_PWM_2.oligo_scores;

      // Select the best scores between FWD and REV strand (for each position)
      check_best_strand_oligo();

      // Fill the local scores vector with scores found in position where the
      // matrix has been generated
      all_local_scores.emplace_back(
          oligo_scores_horizontal_BEST[local_pos - 1]);
      // Fill the global scores vector with all scores generated from shifting
      // and selected from check_best_strand_oligo function
      all_global_scores.insert(all_global_scores.end(),
                               oligo_scores_horizontal_BEST.begin(),
                               oligo_scores_horizontal_BEST.end());
    }

    // If analysis is in Single strand
    else {

      // Local best scores are all from FWD strand
      all_local_scores.emplace_back(oligo_scores_horizontal_FWD[local_pos - 1]);
      all_global_scores.insert(all_global_scores.end(),
                               oligo_scores_horizontal_FWD.begin(),
                               oligo_scores_horizontal_FWD.end());
    }

    // Clearing of horizontal best score for the next sequence cycle
    oligo_scores_horizontal_BEST.clear();
  }
}

// Function to select the best scores between FWD and REV strand (for each
// sequence position)
void z_test_class::check_best_strand_oligo() {
  // PROFILE_FUNCTION();
  // For all oligo scores
  for (unsigned int oligo = 0; oligo < oligo_scores_horizontal_FWD.size();
       oligo++) {

    // If FWD is better than REV
    if (oligo_scores_horizontal_FWD[oligo] >=
        oligo_scores_horizontal_REV[oligo]) {

      oligo_scores_horizontal_BEST.emplace_back(
          oligo_scores_horizontal_FWD[oligo]);
    }

    // If REV is better than FWD
    else {

      oligo_scores_horizontal_BEST.emplace_back(
          oligo_scores_horizontal_REV[oligo]);
    }
  }
}

// Function to calculate all the parameters useul to z-score calculation
void z_test_class::z_score_parameters_calculation() {
  // PROFILE_FUNCTION();
  double local_sum =
      accumulate(all_local_scores.begin(), all_local_scores.end(), 0.0);
  double global_sum =
      accumulate(all_global_scores.begin(), all_global_scores.end(), 0.0);
  double tot_sq_sum_global =
      inner_product(all_global_scores.begin(), all_global_scores.end(),
                    all_global_scores.begin(), 0.0);

  double tot_sq_sum_local =
      inner_product(all_local_scores.begin(), all_local_scores.end(),
                    all_local_scores.begin(), 0.0);
  global_mean = global_sum / all_global_scores.size();
  local_mean = local_sum / all_local_scores.size();
  global_dev_std = sqrt(tot_sq_sum_global / all_global_scores.size() -
                        global_mean * global_mean);
  local_dev_std = sqrt(tot_sq_sum_local / all_local_scores.size() -
                       local_mean * local_mean);
}

// Z-score calculation function
void z_test_class::z_score_calculation(unsigned int len) {
  // PROFILE_FUNCTION();

  z_score = ((global_mean - local_mean) /
             (global_dev_std / sqrt(all_local_scores.size())));

  const double Z = z_score;
  Zpvalue = gsl_cdf_ugaussian_P(Z);

  Zpvalue = check_p_value(Zpvalue, "");

  Zpvalue_bonf = Zpvalue * len;
}


// Function to check, given an oligo as input, if this oligo is palindrome or
// not
bool check_palindrome(string bases, string &reverse_bases) {
  reverse_bases.clear();
  // For any character of the string insert into another string (called reverse
  // bases) the complementary character
  for (int i = bases.length() - 1; i >= 0; i--) {

    char base;
    base = bases[i];
    switch (base) {

    case 'A':
      reverse_bases.append("T");
      break;
    case 'T':
      reverse_bases.append("A");
      break;
    case 'G':
      reverse_bases.append("C");
      break;
    case 'C':
      reverse_bases.append("G");
      break;
    case 'N':
      reverse_bases.append("N");
      break;
    }
  }

  // If they are equal --> it means that the oligo "bases" is palindrome
  if (reverse_bases == bases) {
    return true;
  } else {
    return false;
  }
}

// If the p value is rounded to 0 assigne it a standar low value
// of 1.000001e-300 to avoid possible future errors
double check_p_value(double p, string oligo) {
  if (p == 0 || p <= 1.000001e-300) {

    p = 1.000001e-300;
  }
  if (oligo.find("N") != string::npos) {
    p = 1;
  }
  return p;
}

/////DEBUG/////////////////////////////////////////////////////////

vector<vector<double>> matrix_class::return_inverse_log_matrix() {
  // PROFILE_FUNCTION();
  return inverse_matrix_log;
}

vector<vector<double>> matrix_class::return_log_matrix() {
  // PROFILE_FUNCTION();
  return matrix_log;
}

// Function to analyse the RAM usage of the tool, it returns the maximum amount
// of memory allocated by the program
void RAM_usage() {
  int who = RUSAGE_SELF;
  struct rusage usage;
  int ret;
  ret = getrusage(who, &usage);
  cout << ret << endl;
  cout << endl
       << "Maximum resident set size: " << usage.ru_maxrss / 1000 << " Mb"
       << endl
       << "User CPU time used: " << usage.ru_utime.tv_sec << " s" << endl
       << "System CPU time used: " << usage.ru_stime.tv_usec << " micros"
       << endl;
}

// Debug function: Print sequences and coordinates from GEP vector into a .fasta
// file to check if the sequences extraction is correct
void coordinator_class::print_GEP(vector<bed_class> &GEP) {
  // PROFILE_FUNCTION();
  // Twobit_JASPAR_Bed used to create GEP vector saved into alias file to name
  // the outputs
  alias_file = (TWOBIT_FILE.erase(0, TWOBIT_FILE.find_last_of("/") + 1) + "_" +
                JASPAR_FILE.erase(0, JASPAR_FILE.find_last_of("/") + 1) + "_" +
                BED_FILE.erase(0, BED_FILE.find_last_of("/") + 1));

  // Output file .bed carrying the centered coordinates
  ofstream outfile;
  JASPAR_FILE =
      JASPAR_FILE.erase(JASPAR_FILE.find_last_of("."), JASPAR_FILE.size());
  outfile.open(alias_file);

  for (unsigned int i = 0; i < GEP.size(); i++) {

    outfile << GEP[i].chr_coord << ":" << GEP[i].start_coord << "-"
            << GEP[i].end_coord << endl;
  }

  outfile.close();

  // Output file .fasta carrying the centered coordinates and the sequences
  // extracted
  BED_FILE = BED_FILE.erase(BED_FILE.find_last_of("."), BED_FILE.size());
  outfile.open(alias_file + ".fasta");

  for (unsigned int i = 0; i < GEP.size(); i++) {

    outfile << ">" << GEP[i].chr_coord << ":" << GEP[i].start_coord << "-"
            << GEP[i].end_coord << endl;
    outfile << GEP[i].sequence << endl;
  }

  outfile.close();
}

void Outfile_PWM_matrices(unsigned int j, vector<string> &seed_oligo) {
  // // PROFILE_FUNCTION();
  ofstream outfile;
   if (DS) {
     outfile.open(to_string(kmers_vector[j]) + "-mers_PWM_hamming_matrices_" +
                   alias_file + "DS.txt");
    if (tomtom) {
      print_debug_PWM_hamming_tomtom(outfile, j, kmers_vector[j]);
    } else {
      print_debug_PWM_hamming(outfile, j, kmers_vector[j], seed_oligo);
    }
    outfile.close();
  }

  else {
    outfile.open(to_string(kmers_vector[j]) + "-mers_PWM_hamming_matrices_" +
                 alias_file + "SS.txt");
    if (tomtom) {
      print_debug_PWM_hamming_tomtom(outfile, j, kmers_vector[j]);
    } else {
      print_debug_PWM_hamming(outfile, j, kmers_vector[j], seed_oligo);
    }
    outfile.close();
  }
}

void print_debug_PWM_hamming_tomtom(ofstream &outfile, unsigned int j, 
                                  unsigned int k) {
  // PROFILE_FUNCTION();
  vector<vector<double>> PWM_hamming;
  string ACGT = "ACGT";
  for (unsigned int position = 0; position < Z_TEST_MATRIX[j].size();
       position++) {
    PWM_hamming =
        H_HAMMING_MATRIX[j][position].PWM_hamming;
    outfile << ">Position" << Z_TEST_MATRIX[j][position].local_pos << "("
            << position << ")" << 
            " " << Z_TEST_MATRIX[j][position].Zpvalue << endl;
    for (unsigned int i = 0; i < PWM_hamming.size(); i++) {
      outfile << ACGT[i] << "\t"
              << "["
              << "\t";
      for (unsigned int j = 0; j < PWM_hamming[i].size(); j++) {
        if (!exp_max.empty()){
        outfile << round(PWM_hamming[i][j] * 100) << "\t";
        }
        else{
          outfile << PWM_hamming[i][j] << "\t";
        }
      }
      outfile << "]\n";
  
    }
  }
}
// PWM_matrices, parameters to calculate z-score, z-score and p-value printing
void print_debug_PWM_hamming(ofstream &outfile, unsigned int j, 
                              unsigned int k, vector<string> &seed_oligo) {
  // PROFILE_FUNCTION();
  outfile << "#PWM Matrices calculated from the best oligo for each position "
             "and his hamming distanced oligos - k = "
          << k << endl
          << endl;

  vector<vector<double>> PWM_hamming;
  string ACGT = "ACGT";

  for (unsigned int position = 0; position < Z_TEST_MATRIX[j].size();
       position++) {

    PWM_hamming =
        H_HAMMING_MATRIX[j][position].PWM_hamming;

    outfile << "#Position " << Z_TEST_MATRIX[j][position].local_pos
            << ": \n#PWM calculated from oligo "
            << seed_oligo[position]
            << "\n\n";

    for (unsigned int i = 0; i < PWM_hamming.size(); i++) {

      outfile << ACGT[i] << "\t"
              << "["
              << "\t";

      for (unsigned int j = 0; j < PWM_hamming[i].size(); j++) {

        outfile << PWM_hamming[i][j] << "\t";
      }
      outfile << "]\n";
    }

    outfile << endl;
    outfile << "The global mean is: " << Z_TEST_MATRIX[j][position].global_mean
            << endl;
    outfile << "The global standard deviation is: "
            << Z_TEST_MATRIX[j][position].global_dev_std << endl;
    outfile << "The local mean is: " << Z_TEST_MATRIX[j][position].local_mean
            << endl;
    outfile << "The local standard deviation is: "
            << Z_TEST_MATRIX[j][position].local_dev_std << endl
            << endl;
    outfile << "The zscore calculated is: "
            << Z_TEST_MATRIX[j][position].z_score << endl
            << endl;
    outfile << "The pvalue calculated from the Z score is: "
            << Z_TEST_MATRIX[j][position].Zpvalue << endl
            << endl;
    outfile
        << "-------------------------------------------------------------------"
        << endl;
  }
}

// // Print debug for PWM_hamming outfile -> Selection of output filename
void Outfile_Z_score_values(unsigned int j, vector<string> &seed_oligo) {
  // PROFILE_FUNCTION();
  ofstream outfile;
    if (DS) {

      outfile.open(to_string(kmers_vector[j]) + "-mers_Z_scores_" + alias_file +
                   "DS.txt");
      print_debug_Z_scores(outfile, j, kmers_vector[j], seed_oligo);
      outfile.close();
    }

    else {

      outfile.open(to_string(kmers_vector[j]) + "-mers_Z_scores_" + alias_file +
                   "SS.txt");

      print_debug_Z_scores(outfile, j, kmers_vector[j], seed_oligo);
      outfile.close();
    }
}

// // PWM_matrices, parameters to calculate z-score, z-score and p-value printing
void print_debug_Z_scores(ofstream &outfile, unsigned int j,
                                     unsigned int k, vector<string> &seed_oligo) {
  // PROFILE_FUNCTION();
  outfile << "#Z_score parameters and p-value for hit positions - k = " << k
          << " and -f = " << freq_vector[j] << endl
          << endl;
  string best_oligo;

  outfile << "#Position"
          << "\t"
          << "best_oligo"
          << "\t"
          << "Local_mean"
          << "\t"
          << "Global_mean"
          << "\t"
          << "Local_std_dev"
          << "\t"
          << "Global_std_dev"
          << "\t"
          << "Z_score"
          << "\t"
          << "P-value"
          << "\t"
          << "P-value_Log10"
          << "\t"
          << "Bonferroni P-value"
          << "\t"
          << "Bonferroni_Log10\n";

  for (unsigned int position = 0; position < Z_TEST_MATRIX[j].size();
       position++) {

    double Zpvalue_Log10 = abs(log10(Z_TEST_MATRIX[j][position].Zpvalue));
    double bonferroni_Log10 = abs(log10(Z_TEST_MATRIX[j][position].Zpvalue_bonf));

    // best_oligo = HAMMING_MATRIX[j][local_pos-1].real_best_oligo;

    outfile << Z_TEST_MATRIX[j][position].local_pos << "\t"
            << seed_oligo[position] << "\t" 
            << Z_TEST_MATRIX[j][position].local_mean << "\t"
            << Z_TEST_MATRIX[j][position].global_mean << "\t"
            << Z_TEST_MATRIX[j][position].local_dev_std << "\t"
            << Z_TEST_MATRIX[j][position].global_dev_std << "\t"
            << Z_TEST_MATRIX[j][position].z_score << "\t"
            << Z_TEST_MATRIX[j][position].Zpvalue << "\t" << Zpvalue_Log10
            << "\t" << Z_TEST_MATRIX[j][position].Zpvalue_bonf << "\t" 
            << bonferroni_Log10 << endl;
  }
}

vector<double> freq_vector_creation(string numbers){
   // PROFILE_FUNCTION();
  int index = 0;
  vector<double> vec;

  // When index is == -1 means that it is pointing to the end character of the
  // string
  while (index != -1) {

    // Find the first "," character into the string
    index = numbers.find(",");
    double freq = stod(numbers.substr(0, index));
     if (freq == 0) {
      cout << "WARNING: frequency threshold 0 inserted\n";
    }

    if (freq < 0 || freq >= 1) {
      cerr << "ERROR: please insert a frequency treshold between 0 and "
              "1.\n\n\n";
      display_help();
      exit(1);
    }
    // Put everything before "," into a vector
    vec.emplace_back(freq);

    // Erase from the string the number already inserted and the ","
    numbers.erase(0, index + 1);
  }

  return vec;
  
  for(unsigned int i = 0; i < freq_vector.size(); i++){
    if (freq_vector[i] == 0) {
      cout << "WARNING: frequency threshold 0 inserted\n";
    }

    if (freq_vector[i] < 0 || freq_vector[i] >= 1) {
      cerr << "ERROR: please insert a frequency treshold between 0 and "
              "1.\n\n\n";
      display_help();
      exit(1);
    }
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////

////////////////////PARSER////////////////////////////////////////////////////////////////////
void command_line_parser(int argc, char **argv) {

  const char *const short_opts = "hp:k:b:j:m:ud:o:f:lr:t:e:z:s";

  // Specifying the expected options
  const option long_opts[] = {
      {"help", no_argument, nullptr, 'h'},
      {"param", required_argument, nullptr, 'p'},
      // {"ntop", required_argument, nullptr, 'n'},
      {"kmer", required_argument, nullptr, 'k'},
      // {"all", no_argument, nullptr, 'a'},
      {"tomtom", no_argument, nullptr, 'l'},
      {"unidirection", no_argument, nullptr, 'u'},
      // {"refine", no_argument, nullptr, 'r'},
      {"freq", required_argument, nullptr, 'f'},
      {"distance", required_argument, nullptr, 'd'},
      {"bed", required_argument, nullptr, 'b'},
      {"ordering", required_argument, nullptr, 'o'},
      {"jaspar", required_argument, nullptr, 'j'},
      {"mf", required_argument, nullptr, 'm'},
      {"twobit", required_argument, nullptr, 't'},
      {"ss", no_argument, nullptr, 's'},
      {"exp_maximization", required_argument, nullptr, 'e'},
      {"seconday_matrices", required_argument, nullptr, 'r'},
      {"z_pval_threshold", required_argument, nullptr, 'z'},
      {nullptr, no_argument, nullptr, 0}};

  while (true) {
    const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

    if (-1 == opt)
      break;

    switch (opt) {
    case 'h':
      display_help();
      break;
    case 'p':
      half_length = stoi(optarg);
      break;
    // case 'n':
    //   top_N = stoi(optarg);
    //   break;
    case 'b':
      BED_FILE = string(optarg);
      is_file_exist(BED_FILE, "--bed || -b");
      break;
    case 'j':
      JASPAR_FILE = string(optarg);
      is_file_exist(JASPAR_FILE, "--jaspar || -j");
      break;
    case 't':
      TWOBIT_FILE = string(optarg);
      is_file_exist(TWOBIT_FILE, "--twobit || -t ");
      break;
    case 'o':
      ordering = string(optarg);
      if (ordering != "p") {
        cerr << "ERROR: Wrong -o parameter inserted.\n\n\n";
        display_help();
        exit(1);
      }
      break;
    case 'k':
      kmers.clear();
      kmers = string(optarg);
      break;
    case 'a':
      local_maxima_grouping = true;
      break;
    case 'l':
      tomtom = true;
      break;
    case 'u':
      direction = true;
      break;
    // case 'r':
    //   refining_matrix = true;
    //   break;
    case 'e':
      exp_max = string(optarg);
      break;
    case 'f':
      freq_threshold.clear();
      freq_threshold = string(optarg);
      break;
    case 'd':
      dist.clear();
      dist = string(optarg);
      break;
    case 's':
      DS = false;
      break;
    case 'm':
      MFASTA_FILE = string(optarg);
      is_file_exist(MFASTA_FILE, "--mf || -m ");
      break;
    case 'r':
      max_matrix = stoi(optarg);
      break;
    case 'z':
      z_pval_threshold = stod(optarg); 
      break;
    case '?': // Unrecognized option
    default:
      display_help();
      break;
    }
  }
  check_input_file();
}

bool is_file_exist(string fileName,
                   string buf) { // Input files existence control

  struct stat check;
  int regular_check, existing_check;
  const char *C_fileName = fileName.c_str();
  existing_check = stat(C_fileName, &check);

  regular_check = S_ISREG(check.st_mode);

  if (regular_check == 0 || existing_check != 0) {
    cerr << "ERROR: " << buf << " file does not exist!\n\n";
    display_help();
    exit(1);
  }
  return 0;
}

void check_input_file() {

  if (!MFASTA_FILE.empty() &&
      (!BED_FILE.empty() || !JASPAR_FILE.empty() ||
       !TWOBIT_FILE.empty())) {

    cerr << "FATAL ERROR: Too many input arguments!\nPlease insert Multifasta "
            "file or Bed, Twobit and Jaspar file.\n\n";
    display_help();
    exit(1);
  }
  if ((TWOBIT_FILE.empty() || JASPAR_FILE.empty() ||
       BED_FILE.empty()) &&
      MFASTA_FILE.empty()) {
    cerr << "FATAL ERROR: some arguments needed \n" << endl;
    display_help();
    exit(1);
  }
}

void display_help() {
  cerr << "\n --help || -h show this message\n";
  cerr << "\n --bed || -b <file_bed>: input bed file\n";
  cerr << "\n --kmer || -k <n1,n2,..,nN>: to select k-mers length for the "
          "analysis(DEFAULT: 6,8,10)\n";
  cerr << "\n --twobit || -t <file_twobit>: input twobit file\n";
  cerr << "\n --jaspar || -j <JASPAR_file>: input JASPAR file\n";
  cerr << "\n --param || -p <half_length>: half_length to select bases number "
          "to keep around the chip seq signal (DEFAULT: 150) \n";
  cerr << "\n --mf || -m <multifasta-file>: use multifasta instead of bed file "
          "[ -j,-b,-t,-p options not needed ]\n";
  cerr << "\n -s || --ss as input to make the analysis along the single "
          "strand. (DEFAULT: double strand)\n";
  cerr << "\n -o || --ordering 'p' to order the top N oligos by p-value and "
          "not by occurrences. (DEFAULT: ordering by occurrences)\n";
  cerr << "\n --distance || -d <n1,n2,...,nN> to select the hamming distances. "
          "(DEFAULT: 1,2,3)\n";
  cerr << "\n --freq || -f <n1,n2,...,nN> to set the frequence treshold to "
          "calculate the z_score. (DEFAULT: 0.006, 0.004, 0.003)\n";
  cerr << "\n --exp_maximization || -e to refine PWM matrices with the "
          "expectation maximization method, if you type a number this will be "
          "the number of cycles but if you want to reach convergence you can "
          "type just 'c'\n\n";
  cerr << "\n --tomtom || -l will give as output a format of matrices adapted "
          "for tomtom analysis\n\n";
  cerr << "\n --unidirection || -u parameter orders the sequences based on the "
          "matrix direction \n\n";
  cerr << "\n --secondary_matrices || -r parameter for secondary matrices \n\n";
  cerr << "\n --z_pval_threshold || -z parameter to set a threshold for the PWM's "
          "Z_pvalue (DEFAULT: 1)\n\n";
  exit(EXIT_SUCCESS);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
