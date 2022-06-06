#include "./TwoBit/twobit.c"
#include "./TwoBit/twobit.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <vector>

using namespace std;

#pragma once
#include <algorithm>
#include <thread>

#define PROFILING 1
#if PROFILING
#define PROFILE_SCOPE(name) InstrumentationTimer timer##__LINE__(name)
#define PROFILE_FUNCTION() PROFILE_SCOPE(__PRETTY_FUNCTION__)
#else
#define PROFILE_SCOPE(name)
#endif

////////////////////////GLOBAL VARIABLES////////////////////////
const unsigned int overhead = 25;
const double pseudoc = 0.01;
vector<bool> rev;
TwoBit *tb;
double sim_tresh = 0.001;
string reverse_bases;
vector<unsigned int> kmers_vector;
vector<unsigned int> distance_vector;
vector<unsigned int> len;

////////////////////////PARSER VARIABLES////////////////////////
string BED_FILE;
int half_length = 150;
string TWOBIT_FILE;
string JASPAR_FILE;
string alias_file = "multifasta_";
string MFASTA_FILE;
string ordering;
bool DS = true;
string kmers = "6,8,10";
string dist = "1,2,3";
unsigned int top_N = 10;
double freq_treshold = 0.02;
bool local_maxima_grouping = false;
bool refining_matrix = false;
string exp_max;
bool tomtom = false;
bool err = false;
bool direction = false;
unsigned int seed_vertical = 0;

double pval_threshold = 10e-30;
unsigned int max_matrix = 2;

class Timer {
  public:
    Timer() { m_StartTimepoint = chrono::high_resolution_clock::now(); }
    ~Timer() { Stop(); }

    void Stop() {
      auto endTimepoint = chrono::high_resolution_clock::now();

      auto start = chrono::time_point_cast<chrono::microseconds>(m_StartTimepoint)
                     .time_since_epoch()
                     .count();
      auto end = chrono::time_point_cast<chrono::microseconds>(endTimepoint)
                   .time_since_epoch()
                   .count();

      auto duration = end - start;
      double ms = duration * 0.000001;

      cout << duration << "us (" << ms << "s)\n";
    }

  private:
    chrono::time_point<chrono::high_resolution_clock> m_StartTimepoint;
};

class bed_class {

private:
  string chr_coord;
  unsigned int start_coord;
  unsigned int end_coord;

  bool flag;

  friend class coordinator_class;
  friend class z_test_class;


  friend class MapClass;
  void read_line(string, char **);

public:
  string sequence;
  // Bed class constructor if input is a Multifasta file
  bed_class(string seq) {

    chr_coord = "MULTIFASTA";
    start_coord = 0;
    end_coord = 0;
    sequence = seq;
  }

  // Bed class constructor if input are Bed-Twobit-Jaspar files
  bed_class(string line, TwoBit *tb, unsigned int n_line, char ** result) {

    // Take line from bed file and extract chr_coord, start_coord and end_coord
    read_line(line, result);

    // Set the new start and end coordinates following p (half_length) input and
    // add overhead to end coordinates
    centering_function(start_coord, end_coord, half_length, overhead);

    // Extract from twobit genome the sequence following chr, start and end
    // coordinates
    extract_seq(tb, n_line);
  }

  void centering_function(unsigned int, unsigned int, int, const unsigned int);
  void extract_seq(TwoBit *, unsigned int);
};

class matrix_class {

private:
  vector<vector<double>> matrix;
  vector<vector<double>> norm_matrix;
  vector<vector<double>> inverse_norm_matrix;
  vector<vector<double>> matrix_log;
  vector<vector<double>> inverse_matrix_log;

  void matrix_normalization_pseudoc(vector<vector<double>> &);
  void matrix_normalization(vector<vector<double>> &);
  void matrix_logarithmic(vector<vector<double>> &);
  vector<vector<double>> reverse_matrix(vector<vector<double>> &);
  vector<double> find_col_sum(vector<vector<double>> &);
  void print_debug_matrix(vector<vector<double>>, string);

public:
  matrix_class(vector<vector<double>> &mat) {

    // Function to normalize the matrix scores and add a pseudocount
    matrix_normalization_pseudoc(mat);

    // Function to normalize again the matrix after pseudocount addition
    matrix_normalization(norm_matrix);

    // Function to calculate logarithmic matrix from normalized matrix
    matrix_logarithmic(norm_matrix);

    // Function to reverse the logarithmic normalized matrix to read the oligo
    // in reverse strand
    inverse_matrix_log = reverse_matrix(matrix_log);
  }


  vector<vector<double>> return_log_matrix();
  vector<vector<double>> return_inverse_log_matrix();

};

class oligo_class {

private:
  vector<double> o_matrix_mins;
  double min_possible_score;
  double max_possible_score;
  double best_score;
  unsigned int local_position;
  unsigned int start_coord_oligo;

  void find_minmax(vector<vector<double>> &);
  unsigned int find_best_score();
  void scores_normalization();
  friend class coordinator_class;
  friend class z_test_class;

public:
  vector<double> oligo_scores;

  oligo_class(vector<vector<double>> &matrix, string &sequence,
              unsigned int start_coord_GEP) {

    // Function to annotate in min_possible_score and max_possible_score the
    // worst and the best score that an oligo can reach based on the current
    // jaspar matrix
    find_minmax(matrix);

    // for each oligo in the current sequence a total score of similarity is
    // calculated against the JASPAR matrix
    shifting(matrix, sequence);

    // Function to normalize the scores with the normalization formula
    scores_normalization();

    // Find the best score position and save it into local_position variable (If
    // find more than one select as best the nearest to the center
    local_position = find_best_score();

    start_coord_oligo = start_coord_GEP + local_position;
  }

  oligo_class(vector<vector<double>> &matrix, string &sequence) {

    find_minmax(matrix);

    shifting(matrix, sequence);

    scores_normalization();
  }

  void shifting(vector<vector<double>> &, string &);
};

class coordinator_class { // Coordinator class to connect Matrix to Bed and
                          // Oligos_vector

private:
  vector<vector<double>> matrix;
  vector<vector<double>> matrix_log;
  vector<vector<double>> inverse_matrix_log;
  vector<oligo_class> oligos_vector;
  char ** result;

  void centering_oligo();
  void best_strand();
  void GEP_creation(vector<bed_class> &, char **);
  void oligos_vector_creation(vector<oligo_class> &, vector<vector<double>> &,
                              vector<vector<double>> &, vector<bed_class> &);
  vector<vector<double>> read_JASPAR();

public:
  coordinator_class() {
    result = twobit_sequence_names(tb);
    // GEP (vector of bed class) creation. An empty GEP vector is passed by
    // reference to be filled and saved in this class
    GEP_creation(GEP, result);

    // reading Jaspar file and returning scores as a matrix of double, saved in
    // a variable called matrix
    read_JASPAR();

    // Creating matrix class: input matrix scores, name and tf matrix name
    matrix_class M(matrix);

    // matrix_log and inverse_matrix_log calculated are returned from Matrix
    // class to be saved here --> These are the two matrices on which the
    // analysis will be performed
    matrix_log = M.return_log_matrix();
    inverse_matrix_log = M.return_inverse_log_matrix();

    // The sequences are shifting on the matrix and, for each position, an oligo
    // score is calculated, based on log and inverse_log matrices An oligo_class
    // is created for each sequence shifting to analyze all oligo scores. All
    // the information are saved on oligos_vector, which is an oligo_class
    // vector passed by reference to this function and filled sequence by
    // sequence.
    oligos_vector_creation(oligos_vector, matrix_log, inverse_matrix_log, GEP);

    // If the analysis is performed on Double Strand the function best_strand is
    // useful to select if an oligo has the best score on FWD or REV strand,
    // discarding the other
    if (DS) {
      best_strand();
    }

    // The best oligo selected for each sequence becames the new center of the
    // window, re-setting the GEP coordinates
    centering_oligo();
  }

  vector<bed_class> GEP;

  void print_GEP(vector<bed_class> &);
};

class multifasta_class {

private:
  vector<string> sequences;

  void length_control(vector<string>);
  void extract_sequences();
  void GEP_creation_MF(vector<string>);
  void alias_output_filename();

public:
  vector<bed_class> GEP;

  multifasta_class(string MFASTA_FILE) {

    // Firstly the fasta sequences from multifasta file are extracted and saved
    // into a vector of strings
    extract_sequences();

    // Then the length control is performed --> All the MF sequences must be of
    // the same langth
    length_control(sequences);

    // Function to handle output names
    alias_output_filename();

    // For every sequence in vector "sequences" a bed class is created to store
    // the FASTA seqinto a GEP vector
    GEP_creation_MF(sequences);
  }
};

class z_test_class {

private:
  // zscores and PWM pvalues
  friend class map_class;

  vector<vector<double>> matrix_log;
  vector<vector<double>> inverse_matrix_log;
  vector<double> oligo_scores_horizontal_FWD;
  vector<double> oligo_scores_horizontal_REV;
  vector<double> oligo_scores_horizontal_BEST;
  vector<double> all_global_scores;
  vector<double> all_local_scores;

  void print_debug_oligo_vec(vector<vector<double>> &);
  void oligos_vector_creation_PWM(vector<bed_class> &);
  void z_score_parameters_calculation();
  void z_score_calculation();
  void check_best_strand_oligo();

public:

  double z_score;
  double Zpvalue;


  double global_mean;
  double global_dev_std;
  double local_mean;
  double local_dev_std;
  unsigned int local_pos;
  z_test_class(vector<vector<double>> &PWM_hamming, vector<bed_class> &GEP,
               unsigned int local_p, vector<unsigned int> kmers_vector) {

    local_pos = local_p;

    // Calling matrix class constructor passing PWM_hamming matrix
    matrix_class PWM_hamming_mat(PWM_hamming);

    // Return from matrix class the log_PWM_hamming matrix
    matrix_log = PWM_hamming_mat.return_log_matrix();

    if (DS) {

      // if analysis is in DS return from matrix class the
      // inverse_log_PWM_hamming matrix to shift the reverse strand
      inverse_matrix_log = PWM_hamming_mat.return_inverse_log_matrix();
    }

    // Function to shift PWM_matrix on sequences using oligo class functions
    oligos_vector_creation_PWM(GEP);

    // From local and global scores calculated finding all the parameters useful
    // to z-score calculation
    z_score_parameters_calculation();

    // Calculating z-score and p-value from it
    z_score_calculation();

    all_local_scores.clear();
    all_global_scores.clear();
    oligo_scores_horizontal_FWD.clear();
    oligo_scores_horizontal_REV.clear();
  }
};

vector<z_test_class> Z_TEST_VECTOR;
vector<vector<z_test_class>> Z_TEST_MATRIX;

// Class taking into account all occurrences of the oligos in the input file
class HorizontalClass {
  friend class MapClass;
  friend class PvalueClass;
  friend class HammingClass;
  friend class EMClass;
  private:
    string oligo, oligo_rc;                    // TTGCAT - ATGCAA
    int horizontal_count, horizontal_count_rc; // in occorrenze reverse complement
    int horizontal_count_REV, horizontal_count_rc_REV;
    int horizontal_count_FWD, horizontal_count_rc_FWD;
    bool palindrome;
};

// Class taking into account the occurrences of oligos per columns (pos 1, pos 2, etc.)
class VerticalClass {
  friend class MapClass;
  friend class PvalueClass;
  friend class HammingClass;
  private:
    string oligo, oligo_rc;  
    vector<int> vertical_count,
      vertical_count_rc;
    vector<int> vertical_count_FWD,
      vertical_count_rc_FWD;
    vector<int> vertical_count_REV,
      vertical_count_rc_REV;
    bool palindrome;
};

class MapClass {

  friend class PvalueClass;
  friend class HammingClass;

public:
  unordered_map<string, HorizontalClass> horizontal_map;
  vector<unordered_map<string, HorizontalClass>> vector_map_hor;
  unordered_map<string, VerticalClass> vertical_map;
  vector<unordered_map<string, VerticalClass>> vector_map_ver;
  vector<multimap<int, string, greater<int>>> positions_occurrences;
  vector<vector<multimap<int, string, greater<int>>>>
      vector_positions_occurrences;
  void CountOccurrencesHor(string, int);
  void CountOccurrencesVer(string, int);
  void MainMapVector(vector<bed_class> &);
  void VerticalMapVector();
  void DMainMapVectorDS();
  void DMainMapVectorSS();
  void DVerticalMapVector();

  MapClass(vector<bed_class> &GEP) {
    // Main function where all the maps are constructed 
    MainMapVector(GEP);
    // FIltered vertical map with only the oligo that are effectively present at
    // certain position 
    VerticalMapVector();
    // Functions that starts with D are debug functions
    if(DS){
      // DMainMapVectorDS();
    }
    else{
      // DMainMapVectorSS();
    }
    // DVerticalMapVector();
  }
};

class PvalueClass {

  // K vertical occurences
  // N2 general number of occurrences
  // K = numero occorrenze oligo nella posizione, per ogni elemento del vettore
  // position occurrence: oligo  .
  // N1 = numero occorrenze oligo in tutte le sequenze --> oligo
  // N2 = numero totale di tutti gli oligo di tutte le sequenze - N1 labda su
  // elemento vettore vector_map_hor T = Tot oligo della colonna (numero
  // sequenze)position_occurrences[i].size N1
  unsigned int tot_oligos;
  unsigned int vertical_occurrences;
  void TKN1Calc(vector<bed_class> &, multimap<int, string>::iterator &,
                unordered_map<string, HorizontalClass> &, unsigned int);

public:
  unsigned int K, N1, N2;
  double pvalue;
  string oligo;


  PvalueClass(vector<bed_class> &GEP, multimap<int, string>::iterator &it,
              unordered_map<string, HorizontalClass> &vector_map_hor, unsigned int i) {
    // Here all the paramters for pvalue are calculated
    TKN1Calc(GEP, it, vector_map_hor, i);

  }
};

class HammingClass {

private:
  string seed;
  vector<string> hamming_seed;
  vector<unsigned int> vert_vector;
  unsigned int hamming_v_occ;
  unsigned int hamming_H_occ;
  double freq2;
  unsigned int tot_freq = 0;
  unsigned int pos;

  void HoccCalc(unordered_map<string, HorizontalClass> &);
  void CheckSeed(string, unordered_map<string, VerticalClass> &, multimap<int, string, greater<int>> &,
                  unsigned int, unsigned int);
  void Freq1Calc();
  void Freq2Calc();
  void PWMHammingCalc();
  void DPWMHamming(vector<vector<double>>&);
  void ClearVertical(multimap<int, string, greater<int>> &,
                      unordered_map<string, VerticalClass> &, unsigned int);

public:
  double freq1;
  vector<vector<double>> PWM_hamming;
  map<string, double> cluster_map;

  HammingClass(string seed,
               unordered_map<string, VerticalClass> &map_vertical,
               multimap<int, string, greater<int>> &position_occurrences,
               unordered_map<string, HorizontalClass> &map_horizontal, 
               unsigned int position, unsigned int d) {
    pos = position;
    // Creation of oligo cluster 
    CheckSeed(seed, map_vertical, position_occurrences, position, d);
    // Calculation of freq 1
    Freq1Calc();
    HoccCalc(map_horizontal);
    // Calculation of freq 2
    Freq2Calc();
    PWMHammingCalc();
    ClearVertical(position_occurrences, map_vertical, pos);
    // DPWMHamming(PWM_hamming);
  }
};

vector<HammingClass> H_HAMMING_VECTOR;
vector<vector<HammingClass>> H_HAMMING_MATRIX;

class EMClass{
  private:
    
    map<string, double> like_ratio_map;
    void EM_Ipwm(vector<vector<double>> &);
    void EM_Epart(map<string,double> &, vector<vector<double>> &,unordered_map<string, HorizontalClass> &);
    void EM_Mpart(vector<vector<double>> &);
    bool EM_convergence(vector<vector<double>> &, vector<vector<double>> &, bool);
    void EM_cycle(map<string,double> &, vector<vector<double>> &, unordered_map<string, HorizontalClass> &);
    void print_PWM(string, vector<vector<double>> &);
  
  public:
    
    EMClass(map<string,double> &cluster_map, vector<vector<double>> &PWM_hamming, 
            unordered_map<string, HorizontalClass> &map_horizontal){
      EM_Ipwm(PWM_hamming);
      EM_cycle(cluster_map, PWM_hamming, map_horizontal);
    }
};

string reverse_oligo(string);

// Debug pvalue vector
void DVector(vector<PvalueClass> &, unsigned int);

// Comparison function
bool comp(const PvalueClass &, const PvalueClass &);

void BED_path();
void MULTIFA_path();
void command_line_parser(int, char **);
void display_help();
bool is_file_exist(string fileName, string buf);
void check_input_file();
bool check_palindrome(string, string &);
double check_p_value(double, string);
void RAM_usage();
vector<unsigned int> generic_vector_creation(string);

void Outfile_PWM_matrices(unsigned int, vector<string> &);

void print_debug_PWM_hamming(ofstream &, unsigned int, 
                              unsigned int, vector<string> &);
void print_debug_PWM_hamming_tomtom(ofstream &, unsigned int, 
                                      unsigned int);

void Outfile_Z_score_values(unsigned int, vector<string> &);
void print_debug_Z_scores(ofstream &, unsigned int, unsigned int, 
                            vector<string> &);

