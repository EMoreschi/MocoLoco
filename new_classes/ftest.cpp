#include "ftest.h"
#include <sys/resource.h>

int main(int argc, char *argv[]) {
  Timer timer;

  vector<vector<double>> Gscore;
  vector<double> Best_vector;

  // If arguments number is 1 means that no input file has been inserted -
  // display help
  if (argc == 1) {
    display_help();
  }
  // Collect all the parameters given as input
  command_line_parser(argc, argv);
  
  tb = twobit_open(TWOBIT_FILE.c_str());
  BedClass B(BED_FILE, tb);
  JasparClass J(JASPAR_FILE);
  MatrixClass mat(J.ReturnJaspar());

  for (unsigned int i = 0; i < B.bed_v.size(); i++) {
    
    ScoreClass score(mat.ReturnMatrix(), B.bed_v[i].Sequence, 
                      B.bed_v[i].Start);
    
    Centering(B.bed_v[i], score.start_coord_oligo, 
                    mat.ReturnMatrix()[0].size());
    // if(DS){
    //   ScoreClass score_rev(mat.ReturnInverseMatrix(), B.bed_v[i].Sequence, 
    //                     B.bed_v[i].Start);
          
    // if(score_rev.MaxScore > score.MaxScore){
    //   cout << score_rev.MaxScore << endl;
    //   cout << score_rev.start_coord_oligo << endl;
    //   Centering(B.bed_v[i], score_rev.start_coord_oligo, 
    //                 mat.ReturnInverseMatrix()[0].size());
    //     ReverseCentering(score_rev.start_coord_oligo, B.bed_v[i],
    //             mat.ReturnInverseMatrix()[0].size());
    // }
    // else{
    //   cout << score.MaxScore << endl;
    //   cout << score.start_coord_oligo << endl;
    //     Centering(B.bed_v[i], score.start_coord_oligo, 
    //                 mat.ReturnMatrix()[0].size());
    // }
    // }
    // cout << B.bed_v[i].Sequence << endl;
  }
  twobit_close(tb);

  // Vector len contains all the lengths of sequences for each kmer 
  kmers_vector = generic_vector_creation(kmers);
  distance_vector = generic_vector_creation(dist);
  freq_vector = freq_vector_creation(freq_threshold);
  if(kmers_vector.size() != freq_vector.size()){
    cerr << endl 
      << "ERROR! The number of kmers must be equal to the number of hamming distances and frequencies" << endl << endl;
    exit(1);
  }
  for(unsigned int i = 0; i < kmers_vector.size(); i++){
    len.emplace_back(B.bed_v[0].Sequence.size() - kmers_vector[i] + 1);
  }
  
  MapClass M(B.bed_v);
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
          PvalueClass P(B.bed_v, it, M.vector_map_hor[i], i);
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

        // if(H.freq1 >= freq_vector[i]){
        //   z_test_class Z(H.PWM_hamming, C.GEP, 
        //                   j + 1, len[i]);
        //   Pval = Z.Zpvalue_bonf;
            
        //   //If it is the first cycle of while loop or if the pval is lower 
        //   //than a certain threshold the z_score and PWM are calculated
        //   if(Pval <= (z_pval_threshold * len[i])){
        //       seed_oligo.emplace_back(P_vector[0].oligo);
        //       Z_TEST_VECTOR.emplace_back(Z);
        //       H_HAMMING_VECTOR.emplace_back(H);
        //   }
        // }  
        P_vector.clear(); 
        counter++;
      }
      pos_oligo_vec.clear();
    }
    // Z_TEST_MATRIX.emplace_back(Z_TEST_VECTOR);
    H_HAMMING_MATRIX.emplace_back(H_HAMMING_VECTOR);
    // Z_TEST_VECTOR.clear();
    H_HAMMING_VECTOR.clear();
    
    //Outfile functions 
    // Outfile_PWM_matrices(i, seed_oligo);
    // Outfile_Z_score_values(i, seed_oligo);
  }
  RAM_usage();
  return 0;
}

void BedClass::ReadBed(string BED_FILE, TwoBit *tb) {
  ifstream in(BED_FILE);
  string line;
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

void JasparClass::ReadJaspar(string JASPAR_FILE) {

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
vector<vector<double>> MatrixClass::ReturnMatrix() { return mLogMatrix; }
vector<vector<double>> MatrixClass::ReturnInverseMatrix() { return mInverseLogMatrix; }

// void matrix_c::Norm(std::vector<double> col_sum,
// std::vector<std::vector<double>>& ma){
void MatrixClass::MatrixNormalization(double psdcount, vector<double> col_sum,
                    vector<vector<double>> &mJasparMatrix) {
  for (unsigned int x = 0; x < mJasparMatrix.size(); x++) {
    for (unsigned int i = 0; i < mJasparMatrix[x].size(); i++) {
      double *ptr = &mJasparMatrix[x][i];
      *ptr = (*ptr / col_sum[i] + psdcount);
    }
  }
}

vector<double> MatrixClass::ColSum(vector<vector<double>> &mJasparMatrix) {
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
  for (unsigned int x = 0; x < mJasparMatrix.size(); x++) {
    for (unsigned int i = 0; i < mJasparMatrix[x].size(); i++) {
      double *ptr = &mJasparMatrix[x][i];
      *ptr = log(mJasparMatrix[x][i]);
    }
  }
  mLogMatrix = mJasparMatrix;
}

void ScoreClass::FindMinMax(vector<vector<double>> &matrix) {
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

unsigned int ScoreClass::BestScore(vector<double> &ScoreMatrix){
  MaxScore = *max_element(ScoreMatrix.begin(), ScoreMatrix.end());
  int match = -25;
  for (int j = 0; j < ScoreMatrix.size(); j++) {
    if (ScoreMatrix[j] == MaxScore) {
      if (abs(j-half_length) < abs(match-half_length)){
        match = j; 
      }
    }
  }
  return match;
}

void ReverseCentering(unsigned int center, BedClass::bed_s &GEP, 
                    unsigned int matrix_size){
  
  check_palindrome(GEP.Sequence, reverse_bases);
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
  check_palindrome(GEP.Sequence, reverse_bases);
  GEP.Sequence = reverse_bases;
}

void Centering(BedClass::bed_s &GEP, unsigned int start_coord, 
                unsigned int matrix_size){
  unsigned int center_oligo =
        start_coord + (matrix_size / 2);
  GEP.Start = center_oligo - half_length;
  GEP.End = center_oligo + half_length;
  GEP.Sequence =
        twobit_sequence(tb, GEP.Chromosome.c_str(), GEP.Start, 
                          GEP.End - 1);
}

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
  double cluster_occ = accumulate(vert_vector.begin(), vert_vector.end(),
                                decltype(vert_vector)::value_type(0));
  
  freq1 = static_cast<double>(cluster_occ) / static_cast<double>(tot_freq);
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

// Function used to order the elements of a class (in this case it orders PValueClass elements)
bool comp(const PvalueClass &P1, const PvalueClass &P2) {
  // PROFILE_FUNCTION();
  return (P1.pvalue == P2.pvalue) ? (P1.K > P2.K) : (P1.pvalue < P2.pvalue);
}

bool comp_occ(const PvalueClass &P1, const PvalueClass &P2) {
  // PROFILE_FUNCTION();
  return (P1.K == P2.K) ? ((P1.N1 == P2.N1)? (P1.oligo < P2.oligo) : (P1.N1 < P2.N1)) : (P1.K > P2.K);
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
