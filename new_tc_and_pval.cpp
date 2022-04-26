// This C++ script takes as input a multifasta file and it does the equivalent
// of MocoLoco table_creation_horizontal and vertical functions and the
// p_value_class

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
#include <thread>
#include <typeinfo>
#include <unordered_map>
#include <vector>

using namespace std;

unsigned int secondary = 1;
unsigned int dist = 1;
// Variable for the ordination of pvalues
bool ordering = true;
// Double strand condition
bool DS = false;
// Input file
string MFASTA_FILE;
// k-mer vector
vector<unsigned int> kv = {6};
// sequence length
vector<unsigned int> len;
// sequences vector
vector<string> sequences;

string reverse_oligo(string);

void extract_sequences();

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

class HorizontalClass {
  friend class MapClass;
  friend class PvalueClass;
  friend class HammingClass;
  private:
    string oligo, oligo_rc;                    // TTGCAT - ATGCAA
    int horizontal_count, horizontal_count_rc; // in occorrenze reverse complement
    int horizontal_count_REV, horizontal_count_rc_REV;
    int horizontal_count_FWD, horizontal_count_rc_FWD;
    bool palindrome;
};

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

// class KmerClass {
//   friend class MapClass;
//   friend class PvalueClass;
//   friend class HammingClass;

//   private:
//     string oligo, oligo_rc;                    // TTGCAT - ATGCAA
//     int horizontal_count, horizontal_count_rc; // in occorrenze reverse complement
//     int horizontal_count_REV, horizontal_count_rc_REV;
//     int horizontal_count_FWD, horizontal_count_rc_FWD;
//     vector<int> vertical_count,
//       vertical_count_rc; // vettore di 0 grande come seq -k +1
//     vector<int> vertical_count_FWD,
//       vertical_count_rc_FWD;
//     vector<int> vertical_count_REV,
//       vertical_count_rc_REV;
//     bool palindrome;       // T or F

// };

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
  void MainMapVector();
  void VerticalMapVector();
  void DMainMapVectorDS();
  void DMainMapVectorSS();
  void DVerticalMapVector();

  MapClass() {
    MainMapVector();
    VerticalMapVector();
    // Functions that starts with D are debug functions
    
    if(DS){
      DMainMapVectorDS();
    }
    else{
      DMainMapVectorSS();
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

public:
  unsigned int K, N1, N2;
  double pvalue;
  string oligo;
  unsigned int vertical_occurrences;

  void TKN1Calc(multimap<int, string>::iterator &,
                unordered_map<string, HorizontalClass> &, unsigned int);
  // int N2Calc(unordered_map<string, KmerClass> &);
  void DTKN1Calc();
  void DN2Calc();
  void Dpvalues();
  // void N2vFill(unordered_map<string, KmerClass>);

  PvalueClass(multimap<int, string>::iterator &it,
              unordered_map<string, HorizontalClass> &vector_map_hor, unsigned int i) {
    // N2Calc(vector_map_hor);
    TKN1Calc(it, vector_map_hor, i);
    // Dpvalues();
    // DN2Calc();
    // DTKN1Calc();
  }
};

class SeedClass {
private:
public:
  vector<string> seed_vector; // questo deve diventare mappa posizione seed
  void FillSeed(vector<PvalueClass> &);

  SeedClass(vector<PvalueClass> &P_vector) { FillSeed(P_vector); }
};

class HammingClass {

private:
public:
  string seed;
  vector<string> hamming_seed;
  // vector<vector<string>> hamming_oligos;
  vector<unsigned int> vert_vector;
  unsigned int hamming_v_occ;
  unsigned int hamming_H_occ;
  double freq1, freq2;
  unsigned int tot_freq = 0;
  vector<vector<double>> PWM_hamming;

  void HoccCalc(unordered_map<string, HorizontalClass>);
  void CheckSeed(string, unordered_map<string, VerticalClass>, multimap<int, string, greater<int>>, unsigned int);
  void Freq1Calc();
  void Freq2Calc();
  void PWMHammingCalc();
  void DPWMHamming(vector<vector<double>>&);

  HammingClass(string seed,
               unordered_map<string, VerticalClass> map_vertical,
               multimap<int, string, greater<int>> position_occurrences,
               unordered_map<string, HorizontalClass> map_horizontal, 
               unsigned int position) {
    
    CheckSeed(seed, map_vertical, position_occurrences, position);
    Freq1Calc();
    HoccCalc(map_horizontal);
    Freq2Calc();
    PWMHammingCalc();
    DPWMHamming(PWM_hamming);
  }
};

class EMClass{
  private:

  public:
    
};

// Debug pvalue vector
void DVector(vector<PvalueClass> &);

// Comparison function
bool comp(const PvalueClass &, const PvalueClass &);

// Function from MocoLoco to read a multifasta file and store the sequences in a
// vector
void extract_sequences() {

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

bool comp(const PvalueClass &P1, const PvalueClass &P2) {
  return P1.pvalue < P2.pvalue;
}

void HammingClass::PWMHammingCalc() {
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

      // similar_oligos_map.insert(pair<string, double>(
      //     hamming_seed[j], vert_vector[j]));

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
  for (unsigned short int i = 0; i < PWM_hamming.size(); i++) {
    for (unsigned short int j = 0; j < PWM_hamming[i].size(); j++) {
      cout << PWM_hamming[i][j] << "\t";
    }
    cout << endl;
  }
  cout << endl;
}

void HammingClass::CheckSeed(string seed,
                             unordered_map<string, VerticalClass> map_vertical,
                             multimap<int, string, greater<int>> pos,
                             unsigned int position) {
  
  for (multimap<int, string>::iterator it = pos.begin(); it != pos.end();
       it++) {
    hamming_v_occ = 0;
    tot_freq += it->first;
    string oligo = it->second;
    unsigned int i = 0, count = 0;
    while (seed[i] != '\0') {
      if (seed[i] != oligo[i])
        count++;
      i++;
    }

    if (count <= dist) {
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
      cout << oligo << " " << hamming_v_occ << endl;
    }
  }
  // hamming_oligos.push_back(hamming_seed);
  cout << "Best oligo: " << seed << endl;
  for (unsigned int i = 0; i < hamming_seed.size(); i++) {
    cout << hamming_seed[i] << endl;
  }
  // hamming_seed.clear();
}

void HammingClass::Freq1Calc() {

  freq1 = static_cast<double>(hamming_v_occ) / static_cast<double>(tot_freq);
  cout << "Freq1: " << freq1 << endl;
}

void HammingClass::HoccCalc(unordered_map<string, HorizontalClass> map_horizontal) {
 
  hamming_H_occ = 0;
  for (unsigned int i = 0; i < hamming_seed.size(); i++) {
    string oligo = hamming_seed[i];
    string reverse_o = reverse_oligo(oligo);
    unordered_map<string, HorizontalClass>::iterator it = map_horizontal.find(oligo);
    unordered_map<string, HorizontalClass>::iterator itrc = map_horizontal.find(reverse_o);
    if (it != map_horizontal.end()) {
      hamming_H_occ += it->second.horizontal_count;
      cout << it->first << " " << it->second.horizontal_count << endl;
    } else {
      hamming_H_occ += itrc->second.horizontal_count_rc;
      cout << itrc->second.oligo_rc << " " << itrc->second.horizontal_count_rc << endl;
    }
  }
}

void HammingClass::Freq2Calc() {

  freq2 = static_cast<double>(hamming_v_occ) / static_cast<double>(hamming_H_occ);
  cout << "Freq2: " << freq2 << endl;
}

void SeedClass::FillSeed(vector<PvalueClass> &P_vector) {

  double pval = 0;
  unsigned int vertical_occ = 0;
  for (unsigned int i = 0; i < P_vector.size() && i < secondary; i++) {
    if (ordering) {
      if (P_vector[i].pvalue == pval || pval == 0) {
        pval = P_vector[i].pvalue;
        seed_vector.emplace_back(P_vector[i].oligo);
      } else {
        break;
      }
      cout << "Seed vector size: " << seed_vector.size() << endl;
    } else {
      if (P_vector[i].vertical_occurrences == vertical_occ ||
          vertical_occ == 0) {
        vertical_occ = P_vector[i].vertical_occurrences;
        seed_vector.emplace_back(P_vector[i].oligo);
      } else {
        break;
      }
      cout << "Seed vector size: " << seed_vector.size() << endl;
    }
  }
}
void DVector(vector<PvalueClass> &P_vector) {
  for (unsigned int i = 0; i < P_vector.size(); i++) {
    cout << "Oligo: " << P_vector[i].oligo << endl;
    cout << "K: " << P_vector[i].K << " N1: " << P_vector[i].N1
         << " N2: " << P_vector[i].N2 << endl;
    cout << "Ver_occ: " << P_vector[i].vertical_occurrences
         << " Pval: " << P_vector[i].pvalue << endl;
  }
}

// Function where K, N1, N2 and T are calculated in order to obtain the p value
void PvalueClass::TKN1Calc(multimap<int, string>::iterator &it,
                           unordered_map<string, HorizontalClass> &vector_map_hor,
                           unsigned int i) {

  // T is the number of sequences
  int T = sequences.size();

  // For each oligo in the multimap of vertical occurrences T, N1, N2 and K are
  // calculated.

  // Remember N1 is the number of horizontal occurrences of oligo, T is the
  // total number of sequences, N2 is the total number of oligos in all
  // sequences minus N1 and K is the vertical occurrences of the oligo.
  //  for (multimap<int, string>::iterator it = positions_occurrences.begin();
  //         it != positions_occurrences.end(); it++) {
  N1 = 0;
  vertical_occurrences = it->first;
  K = it->first;
  // Kv.emplace_back(K);
  oligo = it->second;
  // cout << "Oligo: " << oligos << endl;
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
  tot_oligos = sequences.size() * len[i];

  // Calculation of N2
  N2 = tot_oligos - N1;

  // Using the gsl library for the hypergeometric p_value
  pvalue = gsl_cdf_hypergeometric_Q(K, N1, N2, T);
  if (pvalue == 0) {
    pvalue = 1e-300;
  }
  // cout << pval << endl;
  // All the p_values are inserted in a vector

  // N1v.emplace_back(N1);
  // }
}
// N2 calculation with lambda (good for ss but not for ds)
//  int PvalueClass::N2Calc(unordered_map<string, KmerClass> &mappa) {
//    int N2 = accumulate(begin(mappa), end(mappa), 0,
//                    [](unsigned int val,
//                       const unordered_map<string, KmerClass>::value_type &p)
//                       {
//                      return val + (p.second.horizontal_count +
//                                    p.second.horizontal_count_rc);
//                    });
//    return N2;
//  }

// Debug output of pvalues
void PvalueClass::Dpvalues() {

  // for (unsigned int j = 0; j < pvalues.size(); j++){
  //   // cout << "P_value: " << pvalues[j] << endl;
  // }
}

// Fil the vector of N2s
//  void PvalueClass::N2vFill(unordered_map<string, KmerClass> map_vector) {
//    N2v.push_back(N2Calc(map_vector));
//  }

// Other debug functions

// void PvalueClass::DN2Calc() {
//   for (unsigned int i = 0; i < N2v.size(); i++) {
//     cout << "N2 " << i << "  " << N2v[i] << endl;
//   }
// }
//
// void PvalueClass::DTKN1Calc() {
//   for (unsigned int i = 0; i < Kv.size(); i++) {
//     // cout << "N1 "
//     //     << "posizione vettore " << i << "  " << N1v[i] << endl;
//     cout << "K "
//          << "posizione vettore " << i << "  " << Kv[i] << endl;
//   }
// }

// Function for counting the horizontal and vertical occurrences for each oligo
// and putting them inside a map.
void MapClass::CountOccurrencesHor(string sequence, int k) {

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
        it->second.horizontal_count_FWD++;
      }
    } else if (it_rc != horizontal_map.end()) {
      if (!it_rc->second.palindrome) {
        it_rc->second.horizontal_count_rc++;
        it_rc->second.horizontal_count_rc_FWD++;
        if (DS) {
          it_rc->second.horizontal_count_REV++;
          it_rc->second.horizontal_count_rc_FWD++;
          it_rc->second.horizontal_count++;
        }
      }

    } else {
      HorizontalClass Hor;
      Hor.oligo = oligo, Hor.oligo_rc = oligo_rc;
      Hor.palindrome = bool(oligo == oligo_rc);
      Hor.horizontal_count = 1, Hor.horizontal_count_rc = 0;
      if (DS) {
        Hor.horizontal_count_FWD = 1, Hor.horizontal_count_rc_FWD = 0;
        Hor.horizontal_count_REV = 0, Hor.horizontal_count_rc_REV = 1;
        Hor.horizontal_count_rc = 1;
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

  for (unsigned int i = 0; i < (sequence.size() - k + 1); i++) {
    unsigned int tot_freq = 0;
    string oligo = sequence.substr(i, k);
    string oligo_rc = reverse_oligo(oligo);
    unordered_map<string, VerticalClass>::iterator it = vertical_map.find(oligo);
    unordered_map<string, VerticalClass>::iterator it_rc = vertical_map.find(oligo_rc);
    if (it != vertical_map.end()) {
      it->second.vertical_count[i]++;
      it->second.vertical_count_FWD[i]++;
      tot_freq++;
      if (DS) {
        if(it->second.palindrome){         
          it->second.vertical_count[i]++;
          it->second.vertical_count_rc[i]++;
          it->second.vertical_count_REV[i]++;
          it->second.vertical_count_rc_FWD[i]++;
        }
        tot_freq++;
        it->second.vertical_count_FWD[i]++;
        it->second.vertical_count_rc_REV[i]++;
        it->second.vertical_count_rc[i]++;
      }
    } else if (it_rc != vertical_map.end()) {
      if (!it_rc->second.palindrome) {
        it_rc->second.vertical_count_rc[i]++;
        it_rc->second.vertical_count_rc_FWD[i]++;
        if (DS) {
          it_rc->second.vertical_count[i]++;
          it_rc->second.vertical_count_REV[i]++;
          it_rc->second.vertical_count_rc_FWD[i]++;
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

void MapClass::MainMapVector() {
  for (unsigned int i = 0; i < kv.size(); i++) {
    for (unsigned int j = 0; j < sequences.size(); j++) {
      CountOccurrencesHor(sequences[j], kv[i]);
      CountOccurrencesVer(sequences[j], kv[i]);
    }
    vector_map_hor.push_back(horizontal_map);
    horizontal_map.clear();
    vector_map_ver.push_back(vertical_map);
    vertical_map.clear();
  }
}

void MapClass::VerticalMapVector() {

  for (unsigned int i = 0; i < kv.size(); i++) {

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
  for (unsigned int i = 0; i < kv.size(); i++) {
    cout << "I: " << i << endl;
    for (unsigned int j = 0; j < vector_positions_occurrences[i].size(); j++) {
      cout << "J: " << j << endl;
      for (multimap<int, string>::iterator it =
               vector_positions_occurrences[i][j].begin();
           it != vector_positions_occurrences[i][j].end(); it++) {
        cout << it->first << "   " << it->second << endl;
      }
    }
  }
}

// Debug function
void MapClass::DMainMapVectorDS() {
  for (unsigned int i = 0; i < kv.size(); i++) {
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
    for (unordered_map<string, VerticalClass>::iterator it =
             vector_map_ver[i].begin();
         it != vector_map_ver[i].end(); it++) {
      cout << it->first << " " << it->second.vertical_count_FWD[0] << " " 
           << it->second.vertical_count_REV[0] << " " << it->second.vertical_count[0]
           << "\t" << it->second.oligo_rc << " " << it->second.vertical_count_rc_FWD[0]
           << " " << it->second.vertical_count_rc_REV[0] << " "
           << it->second.vertical_count_rc[0] << " "
           << " " << boolalpha << it->second.palindrome << noboolalpha << endl;
    }
  }
}

void MapClass::DMainMapVectorSS() {
  for (unsigned int i = 0; i < kv.size(); i++) {
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
      cout << it->first << " " << it->second.vertical_count_FWD[0] << " " 
           << it->second.vertical_count_REV[0] << " " << it->second.vertical_count[0]
           << "\t" << it->second.oligo_rc << " " << it->second.vertical_count_rc_FWD[0]
           << " " << it->second.vertical_count_rc_REV[0] << " "
           << it->second.vertical_count_rc[0] << " "
           << " " << boolalpha << it->second.palindrome << noboolalpha << endl;
    }
  }
}

int main(int argc, char **argv) {

  MFASTA_FILE = argv[1];
  extract_sequences();
  Timer timer;
  for (unsigned int i = 0; i < kv.size(); i++) {
    len.emplace_back(sequences[0].size() - kv[i] + 1);
  }
  MapClass M;
  vector<PvalueClass> P_vector;
  for (unsigned int i = 0; i < kv.size(); i++) {
    for (unsigned int j = 0; j < len[i]; j++) {
      cout << "Position: " << j << endl;
      for (multimap<int, string>::iterator it =
               M.vector_positions_occurrences[i][j].begin();
           it != M.vector_positions_occurrences[i][j].end(); it++) {
        PvalueClass P(it, M.vector_map_hor[i], i);
        P_vector.push_back(P);
      }
      if (ordering) {
        sort(begin(P_vector), end(P_vector), comp);
      }
      //DVector(P_vector);
      SeedClass S(P_vector);
      for (unsigned int x = 0; x < S.seed_vector.size(); x++) {
        HammingClass H(S.seed_vector[x],
                        M.vector_map_ver[i],
                        M.vector_positions_occurrences[i][j],
                        M.vector_map_hor[i], j);

             cout << H.seed << endl
                   << "vertical occurrences :" << H.hamming_v_occ << endl;
        cout << H.seed << endl
             << "Horizontal occurrences :" << H.hamming_H_occ << endl;
        //      cout << H.freq1 << endl << "freq1" << endl;
        //}
      }
      P_vector.clear();
    }
  }
}
// unordered map <string, KmerClass>
