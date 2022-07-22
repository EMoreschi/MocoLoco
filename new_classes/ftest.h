#include "../TwoBit/twobit.c"
#include "../TwoBit/twobit.h"
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

vector<bool> Reverse;
string reverse_bases;
const unsigned int overhead = 25;
TwoBit *tb;
vector<unsigned int> kmers_vector;
vector<unsigned int> distance_vector;
vector<double> freq_vector;
vector<unsigned int> len;
unsigned int seed_vertical = 0;
double sim_tresh = 0.001;
/////////////////////////////PARSER VARIABLES/////////////////////////////
string BED_FILE;
int half_length = 150;
string TWOBIT_FILE;
string JASPAR_FILE;
string alias_file = "multifasta_";
string MFASTA_FILE;
string ordering;
bool DS = false;
string kmers = "6,8,10";
string dist = "1,2,3";
unsigned int top_N = 10;
string freq_threshold = "0.006, 0.004, 0.003";
bool local_maxima_grouping = false;
bool refining_matrix = false;
string exp_max;
bool tomtom = false;
bool err = false;
bool direction = false;
double z_pval_threshold = 1;
unsigned int max_matrix = 1;

template <typename T> void PrintVector(const T &t) {
  copy(t.cbegin(), t.cend(),
       ostream_iterator<typename T::value_type>(cout, ", "));
  cout << '\n';
}

template <typename T> void PrintMatrix(const T &t) {
  for_each(t.cbegin(), t.cend(), PrintVector<typename T::value_type>);
}

class BedClass {

public:
  struct bed_s {
    string Chromosome, Sequence;
    unsigned int Start, End;
  };
  vector<bed_s> bed_v;
  void ReadBed(string, TwoBit *);
  BedClass(string BED_FILE, TwoBit *tb) { ReadBed(BED_FILE, tb); };
};

class JasparClass {
private:
  vector<vector<double>> mJasparMatrix;
  void ReadJaspar(string);

public:
  vector<vector<double>> ReturnJaspar();

  JasparClass(string line) { ReadJaspar(line); };
};

class MatrixClass {
  vector<vector<double>> mLogMatrix;
  vector<vector<double>> mInverseLogMatrix;
  void MatrixNormalization(double, vector<double>, 
                            vector<vector<double>> &);
  vector<double> ColSum(vector<vector<double>> &);
  void MatrixLog(vector<vector<double>> &);
  void InversemLogMatrix();

public:
  vector<vector<double>> ReturnMatrix();
  vector<vector<double>> ReturnInverseMatrix();
  MatrixClass(vector<vector<double>> mJasparMatrix) {
    const double pseudoc = 0.01;
    MatrixNormalization(pseudoc, ColSum(mJasparMatrix), mJasparMatrix);
    // PrintMatrix(mJasparMatrix);
    // cout << '\n';
    MatrixNormalization(0, ColSum(mJasparMatrix), mJasparMatrix);
    // PrintMatrix(mJasparMatrix);
    // cout << '\n';
    MatrixLog(mJasparMatrix);
    // PrintMatrix(mLogMatrix);
    // cout << '\n';
    InversemLogMatrix();
    // PrintMatrix(mInverseLogMatrix);
  }
};

class ScoreClass {
  vector<double> mMinColumnSum;
  vector<double> mVectorMinMax;
  unsigned int mCentralPosition;
  
public:
  double MaxScore;
  unsigned int start_coord_oligo;
  vector<double> seq_scores;
  void FindMinMax(vector<vector<double>> &);
  vector<double> Shifting(vector<vector<double>> &, string &);
  unsigned int BestScore(vector<double> &);

  ScoreClass(vector<vector<double>> log_matrix, string seq, unsigned int start_coord_GEP) {
    FindMinMax(log_matrix);
    seq_scores = Shifting(log_matrix, seq);
    mCentralPosition = BestScore(seq_scores);
    start_coord_oligo = start_coord_GEP + mCentralPosition;
  };
};

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
  void MainMapVector(vector<BedClass::bed_s> &);
  void VerticalMapVector();
  void DMainMapVectorDS();
  void DMainMapVectorSS();
  void DVerticalMapVector();

  MapClass(vector<BedClass::bed_s> &GEP) {
        

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
  void TKN1Calc(vector<BedClass::bed_s> &, multimap<int, string>::iterator &,
                unordered_map<string, HorizontalClass> &, unsigned int);

public:
  unsigned int K, N1, N2;
  double pvalue;
  string oligo;


  PvalueClass(vector<BedClass::bed_s> &GEP, multimap<int, string>::iterator &it,
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
    // cout << "Size cluster: " << hamming_seed.size() << endl;
    // cout << "Occurrences: " << accumulate(vert_vector.begin(),vert_vector.end(),0) << endl;
    // cout << "Seed oligo: " << seed << endl;
    // DPWMHamming(PWM_hamming);
    ClearVertical(position_occurrences, map_vertical, pos);
    
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

void RAM_usage();
void ReverseCentering(unsigned int, BedClass::bed_s &,
                  unsigned int);
void Centering(BedClass::bed_s &, unsigned int, unsigned int);
bool check_palindrome(string, string &);
string reverse_oligo(string);
// Debug pvalue vector
void DVector(vector<PvalueClass> &, unsigned int);
// Comparison function
bool comp(const PvalueClass &, const PvalueClass &);
bool comp_occ(const PvalueClass &, const PvalueClass &);
//Fare classe debug/error

//Parser functions
void command_line_parser(int, char **);
void display_help();
bool is_file_exist(string fileName, string buf);
void check_input_file();
vector<unsigned int> generic_vector_creation(string);
vector<double> freq_vector_creation(string);