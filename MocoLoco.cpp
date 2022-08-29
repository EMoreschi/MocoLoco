#include "MocoLoco.h"
#include "Profiling.h"
#include <sys/resource.h>

int main(int argc, char *argv[]) {
  Timer timer;
  // Instrumentor::Get().BeginSession("MocoLoco");
  // {
  vector<double> Best_vector;
  // If arguments number is 1 means that no input file has been inserted -
  // display help
  if (argc == 1) {
    display_help();
  }
  // Collect all the parameters given as input
  command_line_parser(argc, argv);

  //If multifasta is not provided the input file is BED file and twobit file 
  //is opened, otherwise if the input is multifasta the tb variable has value 0
  (MFASTA_FILE.empty()) ? tb = twobit_open(TWOBIT_FILE.c_str()) : tb = 0;
  //In BedClass is created a vector of structures with sequence and coordinates
  BedClass B(BED_FILE, MFASTA_FILE, tb);

  if(MFASTA_FILE.empty()){
    JasparClass J(JASPAR_FILE);
    MatrixClass mat(J.ReturnJaspar());
    for (unsigned int i = 0; i < B.bed_v.size(); i++) {
    
      ScoreClass score(mat.ReturnMatrix(), B.bed_v[i].Sequence, 
                      B.bed_v[i].Start);
      if(DS == false){ 
        Best_vector.emplace_back(score.MaxScore);
    	  Centering(B.bed_v[i], score.CenteredStart, 
                    mat.ReturnMatrix()[0].size());
      }
      else{
        ScoreClass score_rev(mat.ReturnInverseMatrix(), B.bed_v[i].Sequence, 
                        B.bed_v[i].Start);
        if(score_rev.MaxScore > score.MaxScore){
          Best_vector.emplace_back(score_rev.MaxScore);
          Centering(B.bed_v[i], score_rev.CenteredStart, 
                     mat.ReturnInverseMatrix()[0].size());
          if(direction){
            ReCentering(score_rev.CenteredStart, B.bed_v[i],
                   mat.ReturnInverseMatrix()[0].size());
          }
        }
        else{
          Best_vector.emplace_back(score.MaxScore);
          Centering(B.bed_v[i], score.CenteredStart, 
                    mat.ReturnMatrix()[0].size());
        }
      }
    }
    //Function to delete sequences with low primary motif scores 
    ClearingGEP(B.bed_v, Best_vector);
    twobit_close(tb);
    print_GEP(B.bed_v);
  }
  // Vector len contains all the lengths of sequences for each kmer 
  for(unsigned int i = 0; i < kmers_vector.size(); i++){
    len.emplace_back(B.bed_v[0].Sequence.size() - kmers_vector[i] + 1);
  }
  
  MapClass M(B.bed_v);
  //For each kmer
  for (unsigned int i = 0; i < kmers_vector.size(); i++) {
    vector<PvalueClass> P_vector;
    vector<string> seed_oligo;
    //For each position in the sequence
    for (unsigned int j = 0; j < len[i]; j++) {
      double Pval = 0;
      unsigned int counter = 0;
      // Loop for eventually other matrices that are hidden by the best motif
      while(counter < max_matrix) {

        cout << "Position: " << j << endl;
          
        //For each oligo present in the vertical map
        for (multimap<int, string>::iterator it =
        M.vector_positions_occurrences[i][j].begin();
        it != M.vector_positions_occurrences[i][j].end(); it++) {
          //In the PvalueClass to each oligo is associated its pvalue
          PvalueClass P(B.bed_v, it, M.vector_map_hor[i], i);
          P_vector.push_back(P);
        }
        //The element in P_vector are ordered on the basis of pvalues
        //otherwise the normal ordering (by occurrences in vertical map)
        //is maintained
        (ordering == "p") ? sort(begin(P_vector), end(P_vector), comp) :
          sort(begin(P_vector), end(P_vector), comp_occ);

        // Debug for PValueClass
        // DVector(P_vector, j);
          
        //Creation of clusters of oligos at hamming distance
        //and creation of PWM for each position
        HammingClass H(P_vector[0].oligo,
                        M.vector_map_ver[i],
                        M.vector_positions_occurrences[i][j],
                        M.vector_map_hor[i], j, i);
        // Perform the expectation maximization algorithm 
        // (https://en.wikipedia.org/wiki/Expectation–maximization_algorithm)
        if (!exp_max.empty()){
          EMClass E(H.cluster_map, H.PWM_hamming, 
                    M.vector_map_hor[i]);
        }
        //If the frequence of seed oligo is higher than a threshold
        //the z_score is calculated

        if(H.freq1 >= freq_vector[i]){
          z_test_class Z(H.PWM_hamming, B.bed_v, 
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
    }
    Z_TEST_MATRIX.emplace_back(Z_TEST_VECTOR);
    H_HAMMING_MATRIX.emplace_back(H_HAMMING_VECTOR);
    Z_TEST_VECTOR.clear();
    H_HAMMING_VECTOR.clear();

    //Outfile functions 
    Outfile_PWM_matrices(i, seed_oligo);
    Outfile_Z_score_values(i, seed_oligo);
    seed_oligo.clear();
  }
  RAM_usage();
  return 0;
// }
//   Instrumentor::Get().EndSession();
}

void BedClass::ReadBed(string BED_FILE, string MFASTA_FILE, TwoBit *tb) {
  // PROFILE_FUNCTION();
  string line;
  vector<string> sequences;
  if(MFASTA_FILE.empty()){
    ifstream in(BED_FILE);
    while (getline(in, line)) {
      if (line.empty() || line[0] == '#') {
        continue;
      }
      istringstream mystream(line);
      bed_s bed_in;
        mystream >> bed_in.Chromosome >> bed_in.Start >> bed_in.End;
        if (bed_in.Start > bed_in.End) {
          cout << "error " << endl;
        } 
        else {
          int center = (bed_in.Start + bed_in.End) / 2;
          bed_in.Start = center - half_length;
          bed_in.End = center + half_length + overhead;
          bed_in.Sequence =
            twobit_sequence(tb, bed_in.Chromosome.c_str(), bed_in.Start, 
                            bed_in.End - 1);

          bed_v.push_back(bed_in);
        }
    }
  }
  else{
    string file_name = MFASTA_FILE;
    alias_file =
      alias_file + file_name.erase(0, file_name.find_last_of("_") + 1);
    alias_file =
      file_name.erase(file_name.find_last_of("."), file_name.size()) +
      "_";
    ifstream file(MFASTA_FILE);
    string current_sequence;
    
    bool first_line = true;

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
  for (unsigned int i = 0; i < sequences.size(); i ++){
    if (sequences[0].size() != sequences[i].size()){
      cerr << "ERROR! The length of fasta sequences isn't equal in all the file!\n\nCheck it out!\n\n";
    }
    bed_s bed_in;
    bed_in.Chromosome = "MULTIFASTA";
    bed_in.Start = 0;
    bed_in.End = 0;
    bed_in.Sequence = sequences[i];
    bed_v.push_back(bed_in);
  }
}

void JasparClass::ReadJaspar(string JASPAR_FILE) {
  // PROFILE_FUNCTION();
  ifstream file(JASPAR_FILE);
  string line;

  // For each line of the JASPAR file
  while (getline(file, line)) {

    // If the line Start with '>' character save the words into matrix_name
    // string and into tf_name string
    if (line[0] == '>') {

      istringstream mystream(line);
      // mystream >> matrix_name >> tf_name;
    }

    // If the line does not Start with '>' delete the '[',']' and 'A,T,C,G'
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
      mJasparMatrix.emplace_back(scores_line);
    }
  }
  file.close();
}

vector<vector<double>> JasparClass::ReturnJaspar() { return mJasparMatrix; }
vector<vector<double>>& MatrixClass::ReturnMatrix() { return mLogMatrix; }
vector<vector<double>>& MatrixClass::ReturnInverseMatrix() { return mInverseLogMatrix; }

// void matrix_c::Norm(std::vector<double> col_sum,
// std::vector<std::vector<double>>& ma){
void MatrixClass::MatrixNormalization(double psdcount, vector<double> col_sum,
                    vector<vector<double>> &mJasparMatrix) {
  // PROFILE_FUNCTION();
  for (unsigned int x = 0; x < mJasparMatrix.size(); x++) {
    for (unsigned int i = 0; i < mJasparMatrix[x].size(); i++) {
      double *ptr = &mJasparMatrix[x][i];
      *ptr = (*ptr / col_sum[i] + psdcount);
    }
  }
}

vector<double> MatrixClass::ColSum(vector<vector<double>> &mJasparMatrix) {
  // PROFILE_FUNCTION();
  vector<double> col_sum;
  for (unsigned int i = 0; i < mJasparMatrix[0].size(); i++) {
    double sum = 0;
    for (unsigned int j = 0; j < 4; j++) {
      sum += mJasparMatrix[j][i];
    }
    col_sum.emplace_back(sum);
  }
  return col_sum;
}

void MatrixClass::MatrixLog(vector<vector<double>> &mJasparMatrix) {
  // PROFILE_FUNCTION();
  for (unsigned int x = 0; x < mJasparMatrix.size(); x++) {
    for (unsigned int i = 0; i < mJasparMatrix[x].size(); i++) {
      double *ptr = &mJasparMatrix[x][i];
      *ptr = log(mJasparMatrix[x][i]);
    }
  }
  mLogMatrix = mJasparMatrix;
}

void ScoreClass::FindMinMax(vector<vector<double>> &matrix) {
  // PROFILE_FUNCTION();
  vector<double> max_sum;
  for (unsigned int i = 0; i < matrix[0].size(); i++) {
    vector<double> column;
    for (unsigned int j = 0; j < matrix.size(); j++) {
      column.emplace_back(matrix[j][i]);
      // cout << matrix[j][i] << endl;
    }
    mMinColumnSum.emplace_back(*min_element(column.begin(), column.end()));
    max_sum.emplace_back(*max_element(column.begin(), column.end()));
  }

  double min = accumulate(mMinColumnSum.begin(), mMinColumnSum.end(), 0.0);
  double max = accumulate(max_sum.begin(), max_sum.end(), 0.0);
  mVectorMinMax.emplace_back(min);
  mVectorMinMax.emplace_back(max);
}

vector<double> ScoreClass::Shifting(vector<vector<double>> &matrix,
                                 string &sequence) {
  // PROFILE_FUNCTION();
  int max = 0;
  max = sequence.size() - matrix[0].size();
  vector<double> seq_scores;

  for (int s_iterator = 0; s_iterator <= max; s_iterator++) {
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

        sum_scores += mMinColumnSum[i];
        break;
      }
    }

    // Best score normalization with normalization formula
    sum_scores = 1 + (sum_scores - mVectorMinMax[1]) / (mVectorMinMax[1] - mVectorMinMax[0]);
    seq_scores.emplace_back(sum_scores);
    // cout << seq_scores[s_iterator];
    //  cout << *seq_scores[s_iterator] << endl;
    //  The total score of an oligo is saved into an oligo_scores vector
  }
  return seq_scores;
}

void MatrixClass::InversemLogMatrix() {
  // PROFILE_FUNCTION();
  mInverseLogMatrix = mLogMatrix;
  reverse(mInverseLogMatrix.begin(), mInverseLogMatrix.end());

  for (int i = 0; i < 4; i++) {

    reverse(mInverseLogMatrix[i].begin(), mInverseLogMatrix[i].end());
  }
}

unsigned int ScoreClass::BestScore(vector<double> &ScoreVector){
  MaxScore = *max_element(ScoreVector.begin(), ScoreVector.end());
  int match = -25;
  for (int j = 0; j < ScoreVector.size(); j++) {
    if (ScoreVector[j] == MaxScore) {
      if (abs(j-half_length) < abs(match-half_length)){
        match = j; 
      }
    }
  }
  return match;
}

void ReCentering(unsigned int center, BedClass::bed_s &GEP, 
                    unsigned int matrix_size){
  // PROFILE_FUNCTION();
  ReverseString(GEP.Sequence, reverse_bases);
  GEP.Sequence = reverse_bases;
  if(matrix_size % 2 == 0){
    unsigned int center_oligo =
            (center + matrix_size / 2);
    GEP.Start = center_oligo - half_length;
    GEP.End = center_oligo + half_length;
  }
  else{
    unsigned int center_oligo =
            (center + matrix_size / 2) + 1;
    GEP.Start = center_oligo - half_length;
    GEP.End = center_oligo + half_length; 
  }
  GEP.Sequence =
          twobit_sequence(tb, GEP.Chromosome.c_str(), GEP.Start, 
                            GEP.End - 1);
  ReverseString(GEP.Sequence, reverse_bases);
  GEP.Sequence = reverse_bases;
}

void Centering(BedClass::bed_s &GEP, unsigned int start_coord, 
                unsigned int matrix_size){
  // PROFILE_FUNCTION();
  unsigned int center_oligo =
        start_coord + (matrix_size / 2);
  GEP.Start = center_oligo - half_length;
  GEP.End = center_oligo + half_length;
  GEP.Sequence =
        twobit_sequence(tb, GEP.Chromosome.c_str(), GEP.Start, 
                          GEP.End - 1);
}

void ReverseString(string bases, string &reverse_bases) {
  // PROFILE_FUNCTION();
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
}

void MapClass::MainMapVector(vector<BedClass::bed_s> &GEP) {
  // PROFILE_FUNCTION();
  for (unsigned int i = 0; i < kmers_vector.size(); i++) {
    for (unsigned int j = 0; j < GEP.size(); j++) {
      CountOccurrencesHor(GEP[j].Sequence, kmers_vector[i]);
      CountOccurrencesVer(GEP[j].Sequence, kmers_vector[i]);
    }
    vector_map_hor.push_back(horizontal_map);
    horizontal_map.clear();
    vector_map_ver.push_back(vertical_map);
    vertical_map.clear();
  }
}

void MapClass::CountOccurrencesHor(string &sequence, unsigned int k) {
  // PROFILE_FUNCTION();
  for (unsigned int i = 0; i < (sequence.size() - k + 1); i++) {
    string oligo = sequence.substr(i, k);
    string oligo_rc;
    ReverseString(oligo, oligo_rc);
    unordered_map<string, HorizontalClass>::iterator it = horizontal_map.find(oligo);
    unordered_map<string, HorizontalClass>::iterator it_rc = horizontal_map.find(oligo_rc);
    if (it != horizontal_map.end()) {
      it->second.horizontal_count++;
      if (DS) {
        it->second.horizontal_count_rc++;
      }
    } else if (it_rc != horizontal_map.end()) {
      if (!it_rc->second.palindrome) {
        it_rc->second.horizontal_count_rc++;
        if (DS) {
          it_rc->second.horizontal_count++;
        }
      }

    } else {
      HorizontalClass Hor;
      Hor.oligo = oligo, Hor.oligo_rc = oligo_rc;
      Hor.palindrome = bool(oligo == oligo_rc);
      Hor.horizontal_count = 1, Hor.horizontal_count_rc = 0;
      if (DS) {
        Hor.horizontal_count_rc = 1, Hor.horizontal_count = 1;
      }
    
      horizontal_map.emplace(oligo, Hor);
    }
  }
}

void MapClass::CountOccurrencesVer(string &sequence, unsigned int k) {
  // PROFILE_FUNCTION();
  for (unsigned int i = 0; i < (sequence.size() - k + 1); i++) {
    // unsigned int tot_freq = 0;
    string oligo = sequence.substr(i, k);
    string oligo_rc;
    ReverseString(oligo, oligo_rc);
    unordered_map<string, VerticalClass>::iterator it = vertical_map.find(oligo);
    unordered_map<string, VerticalClass>::iterator it_rc = vertical_map.find(oligo_rc);
    if (it != vertical_map.end()) {
      it->second.vertical_count[i]++;
      if(DS){
        it->second.vertical_count_rc[i]++;
      }
    } else if (it_rc != vertical_map.end()) {
      if (!it_rc->second.palindrome) {
        it_rc->second.vertical_count_rc[i]++;
        if (DS) {
          it_rc->second.vertical_count[i]++;
        }
      }
    } else {
      VerticalClass Ver;
      Ver.oligo = oligo, Ver.oligo_rc = oligo_rc;
      Ver.palindrome = bool(oligo == oligo_rc);
      Ver.vertical_count.resize(((sequence.size() - k) + 1), 0);
      Ver.vertical_count_rc.resize(((sequence.size() - k) + 1), 0);
      Ver.vertical_count[i] = 1;
      if (DS) {
        Ver.vertical_count_rc[i] = 1;        
      }
    
      vertical_map.emplace(oligo, Ver);
    }
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
  // // PROFILE_FUNCTION();
  
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
  // // PROFILE_FUNCTION();
  for (unsigned int i = 0; i < kmers_vector.size(); i++) {

    for (unordered_map<string, HorizontalClass>::iterator it =
             vector_map_hor[i].begin();
         it != vector_map_hor[i].end(); it++) {
      cout << it->first << " " << it->second.horizontal_count
           << "\t" << it->second.oligo_rc
           << " " << it->second.horizontal_count_rc << " "
           << " " << boolalpha << it->second.palindrome << noboolalpha << endl;
    }
    cout << "\n\n VERTICAL MAP: \n\n" << endl;
    for (unordered_map<string, VerticalClass>::iterator it =
             vector_map_ver[i].begin();
         it != vector_map_ver[i].end(); it++) {
      cout << it->first << " " << it->second.vertical_count[0]
           << "\t" << it->second.oligo_rc << " "
           << it->second.vertical_count_rc[0] << " "
           << " " << boolalpha << it->second.palindrome << noboolalpha << endl;
    }
  }
}

void MapClass::DMainMapVectorSS() {
  // // PROFILE_FUNCTION();
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
      cout << it->first << " " << it->second.vertical_count[21]
           << "\t" << it->second.oligo_rc << " "
           << it->second.vertical_count_rc[21] << " "
           << " " << boolalpha << it->second.palindrome << noboolalpha << endl;
    }
  }
}

// Function where K, N1, N2 and T are calculated in order to obtain the p value
void PvalueClass::TKN1Calc(vector<BedClass::bed_s> &GEP,
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
  string oligo_rc;
  ReverseString(oligo, oligo_rc);
  unordered_map<string, HorizontalClass>::iterator itBigMap =
      vector_map_hor.find(oligo);
  unordered_map<string, HorizontalClass>::iterator itBigMap_rc =
      vector_map_hor.find(oligo_rc);
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
      // cout << "Hamming seed: " << hamming_seed[j] << "\tOcc: " << vert_vector[j] << endl;
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
  // // PROFILE_FUNCTION();
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
      
      string reverse_o;
      ReverseString(oligo, reverse_o);
      unordered_map<string, VerticalClass>::iterator it_ver = map_vertical.find(oligo);
      unordered_map<string, VerticalClass>::iterator it_ver_rc = map_vertical.find(reverse_o);
      if (it_ver != map_vertical.end()) {
        hamming_v_occ += it_ver->second.vertical_count[position];
      } else {
        hamming_v_occ += it_ver_rc->second.vertical_count_rc[position];
      }        
      vert_vector.push_back(hamming_v_occ);
      hamming_seed.push_back(oligo);
      hamming_seed_rev.push_back(reverse_o);
    }
  }
}

void HammingClass::Freq1Calc() {
  // PROFILE_FUNCTION();
  double cluster_occ = accumulate(vert_vector.begin(), vert_vector.end(),
                                decltype(vert_vector)::value_type(0));
  
  freq1 = static_cast<double>(cluster_occ) / static_cast<double>(tot_freq);
}

void HammingClass::HoccCalc(unordered_map<string, HorizontalClass> &map_horizontal) {
 // PROFILE_FUNCTION();
  hamming_H_occ = 0;
  for (unsigned int i = 0; i < hamming_seed.size(); i++) {
    string oligo = hamming_seed[i];
    string reverse_o;
    ReverseString(oligo, reverse_o);
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
  // PROFILE_FUNCTION();
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
  for(unsigned int i = 0; i < hamming_seed_rev.size(); i++){
    
    for(multimap<int, string, greater<int>>::iterator it = pos.begin();
      it != pos.end(); ++it){
      
      string oligo = it->second;

      if(oligo == hamming_seed_rev[i]){
        pos.erase(it);
        break;
      }
    }
  }
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
    string oligo_rc;
    ReverseString(it->first, oligo_rc);
    occ_oligo_it_rev = map_horizontal.find(oligo_rc);

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
    string oligo_rc;
    ReverseString(it->first, oligo_rc);
    occ_oligo_it_rev = map_horizontal.find(oligo_rc);

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

void ClearingGEP(vector<BedClass::bed_s> &GEP, vector<double> &ScoreVector){
  int count = 0;
  //Copy of GEP to keep just the sequences with good scores
  vector<BedClass::bed_s> ClearedGEP;
  //Calculation of standard deviation
  double SumScores =
      accumulate(ScoreVector.begin(), ScoreVector.end(), 0.0);
  double Mean = SumScores / ScoreVector.size();
  double SquareSumScores =
      inner_product(ScoreVector.begin(), ScoreVector.end(),
                    ScoreVector.begin(), 0.0);
  double StdDev = sqrt(SquareSumScores / ScoreVector.size() -
                        Mean * Mean);
  //This is the value used as threshold equal to the mean of the scores 
  //minus two times the standard deviation
  double Comparison = Mean - (2 * StdDev);
  //For each sequence
  for (unsigned int i = 0; i < GEP.size(); i++){
    //If the score of the sequence is higher than the threshold the BedClass structure is loaded in the new vector
    if(ScoreVector[i] > Comparison){
      ClearedGEP.emplace_back(GEP[i]);
    }
    else{
      cerr << "The sequence " << i + 1 << " is not taken into account because low score \n\n"; 
      count += 1;
    }
    cout << Mean << "\t" << StdDev << endl << endl;
    
  }
  cerr << count << " sequences are eliminated \n\n";
  //The old GEP is cleared and replaced with the new GEP
  GEP.clear();
  GEP = ClearedGEP;
}

void DVector(vector<PvalueClass> &P_vector, unsigned int j) {
  // // PROFILE_FUNCTION();
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

// Function used to order the elements of a class (in this case it orders PValueClass elements)
bool comp(const PvalueClass &P1, const PvalueClass &P2) {
  // PROFILE_FUNCTION();
  return (P1.pvalue == P2.pvalue) ? (P1.K > P2.K) : (P1.pvalue < P2.pvalue);
}

bool comp_occ(const PvalueClass &P1, const PvalueClass &P2) {
  // PROFILE_FUNCTION();
  return (P1.K == P2.K) ? ((P1.N1 == P2.N1)? (P1.oligo < P2.oligo) : (P1.N1 < P2.N1)) : (P1.K > P2.K);
}

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

void z_test_class::oligos_vector_creation_PWM(vector<BedClass::bed_s> &GEP) {
  // PROFILE_FUNCTION();
  // For every sequence
  for (unsigned int i = 0; i < GEP.size(); i++) {

    // Calling oligo class to accede to all functions useful to shift a matrix
    // on sequences --> Shifting on FWD strand
    ScoreClass SHIFTING_PWM(mMatrixLog, GEP[i].Sequence);

    // Return oligo scores calculated from previous shifting
    oligo_scores_horizontal_FWD = SHIFTING_PWM.SeqScores;

    // If analysis is in Double strand
    if (DS) {

      // Make the shifting also on reverse strand, putting as input the
      // inverse_log_matrix --> Shifting on REV strand
      ScoreClass SHIFTING_PWM_2(mInverseMatrixLog, GEP[i].Sequence);

      // Retrun oligo scores from previous shifting
      oligo_scores_horizontal_REV = SHIFTING_PWM_2.SeqScores;

      // Select the best scores between FWD and REV strand (for each position)
      check_best_strand_oligo();

      // Fill the local scores vector with scores found in position where the
      // matrix has been generated
      all_local_scores.emplace_back(oligo_scores_horizontal_BEST[LocalPosition - 1]);
      // Fill the global scores vector with all scores generated from shifting
      all_global_scores.insert(all_global_scores.end(),
                               oligo_scores_horizontal_BEST.begin(),
                               oligo_scores_horizontal_BEST.end());
    }

    // If analysis is in Single strand
    else {

      // Local best scores are all from FWD strand
      all_local_scores.emplace_back(oligo_scores_horizontal_FWD[i]);
      all_global_scores.insert(all_global_scores.end(),
                               oligo_scores_horizontal_FWD.begin(),
                               oligo_scores_horizontal_FWD.end());
    }

    // Clearing of horizontal best score for the next sequence cycle
    oligo_scores_horizontal_BEST.clear();
  }
}

void z_test_class::check_best_strand_oligo() {
  // PROFILE_FUNCTION();
  // For all oligo scores
  for (unsigned int oligo = 0; oligo < oligo_scores_horizontal_FWD.size();
       oligo++) {

    (oligo_scores_horizontal_FWD[oligo] >= oligo_scores_horizontal_REV[oligo]) ? 
        oligo_scores_horizontal_BEST.emplace_back(oligo_scores_horizontal_FWD[oligo]) : 
        oligo_scores_horizontal_BEST.emplace_back(oligo_scores_horizontal_REV[oligo]);
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

void Outfile_PWM_matrices(unsigned int j, vector<string> &seed_oligo) {
  // PROFILE_FUNCTION();
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
    outfile << ">Position" << Z_TEST_MATRIX[j][position].LocalPosition << "("
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

    outfile << "#Position " << Z_TEST_MATRIX[j][position].LocalPosition
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

  outfile << "#Position\tbest_oligo\tLocal_mean\tGlobal_mean\tLocal_std_dev"
          << "\tGlobal_std_dev\tZ_score\tP-value\tP-value_Log10"
          << "\tBonferroni P-value\tBonferroni_Log10\n";

  for (unsigned int position = 0; position < Z_TEST_MATRIX[j].size();
       position++) {

    double Zpvalue_Log10 = abs(log10(Z_TEST_MATRIX[j][position].Zpvalue));
    double bonferroni_Log10 = abs(log10(Z_TEST_MATRIX[j][position].Zpvalue_bonf));

    // best_oligo = HAMMING_MATRIX[j][LocalPosition-1].real_best_oligo;
    outfile << Z_TEST_MATRIX[j][position].LocalPosition << "\t"
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

void print_GEP(vector<BedClass::bed_s> &GEP) {
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

    outfile << GEP[i].Chromosome << ":" << GEP[i].Start << "-"
            << GEP[i].End << endl;
  }

  outfile.close();

  // Output file .fasta carrying the centered coordinates and the sequences
  // extracted
  BED_FILE = BED_FILE.erase(BED_FILE.find_last_of("."), BED_FILE.size());
  outfile.open(alias_file + ".fasta");

  for (unsigned int i = 0; i < GEP.size(); i++) {

    outfile << ">" << GEP[i].Chromosome << ":" << GEP[i].Start << "-"
            << GEP[i].End << endl;
    outfile << GEP[i].Sequence << endl;
  }

  outfile.close();
}

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
      kmers_vector = generic_vector_creation(kmers);
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
    case 'e':
      exp_max = string(optarg);
      break;
    case 'f':
      freq_threshold.clear();
      freq_threshold = string(optarg);
      freq_vector = freq_vector_creation(freq_threshold);
      if(kmers_vector.size() != freq_vector.size()){
        cerr << endl 
            << "ERROR! The number of kmers must be equal to the number of hamming distances and frequencies" << endl << endl;
        exit(1);
      }
      break;
    case 'd':
      dist.clear();
      dist = string(optarg);
      distance_vector = generic_vector_creation(dist);
      if(kmers_vector.size() != distance_vector.size()){
    cerr << endl 
      << "ERROR! The number of kmers must be equal to the number of hamming distances and frequencies" << endl << endl;
    exit(1);
  }
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
