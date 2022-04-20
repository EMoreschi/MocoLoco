//This C++ script takes as input a multifasta file and it does the equivalent of 
//MocoLoco table_creation_horizontal and vertical functions and the p_value_class

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
#include <algorithm>
#include <thread>

using namespace std;

//Variable for the ordination of pvalues
bool ordering = false;
//Double strand condition
bool DS = true;
//Input file
string MFASTA_FILE;
//k-mer vector
vector<unsigned int> kv = {6};
//sequence length
vector<unsigned int> len;
//sequences vector
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


class KmerClass {
  friend class MapClass;
  friend class PvalueClass;

private:
  string oligo, oligo_rc;                    // TTGCAT - ATGCAA
  int orizzontal_count, orizzontal_count_rc; // in occorrenze reverse complement
  vector<int> vertical_count, vertical_count_rc; // vettore di 0 grande come seq -k +1
  bool palindrome;  //T or F
};

class MapClass {

  friend class PvalueClass;

public:
  unordered_map<string, KmerClass> mappa;
  vector<unordered_map<string, KmerClass>> mappa_vector;
  vector<multimap<int, string, greater<int> >> positions_occurrences;
  vector<vector<multimap<int, string, greater<int>>>> vector_positions_occurrences;
  void CountOccurrences(string, int);
  void MainMapVector();
  void VerticalMapVector();
  void DMainMapVector();
  void DVerticalMapVector();

  MapClass() {
    MainMapVector();
    VerticalMapVector();
    // Functions that starts with D are debug functions
    // DMainMapVector();
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
  // elemento vettore mappa_vector T = Tot oligo della colonna (numero
  // sequenze)position_occurrences[i].size N1
  unsigned int tot_oligos;

  
public:
  unsigned int K, N1, N2;
  double pvalue;
  string oligo;
  unsigned int vertical_occurrences;

  void TKN1Calc(multimap<int, string>::iterator &,
                unordered_map<string, KmerClass> &, unsigned int);
  // int N2Calc(unordered_map<string, KmerClass> &);
  void DTKN1Calc();
  void DN2Calc();
  void Dpvalues();
  // void N2vFill(unordered_map<string, KmerClass>);



  PvalueClass(multimap<int, string>::iterator &it,
              unordered_map<string, KmerClass> &mappa_vector, unsigned int i) {
    //N2Calc(mappa_vector);
    TKN1Calc(it, mappa_vector, i);
    // Dpvalues();
    // DN2Calc();
    // DTKN1Calc();

  }
};

//Debug pvalue vector
void DVector(vector<PvalueClass>&);

//Comparison function
bool comp(const PvalueClass& , const PvalueClass& );

//Function from MocoLoco to read a multifasta file and store the sequences in a vector
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


bool comp(const PvalueClass& P1, const PvalueClass& P2)
{
  return P1.pvalue < P2.pvalue;
}

void DVector(vector<PvalueClass> &P_vector){
  for (unsigned int i = 0; i < P_vector.size(); i++){
    cout << "Oligo: " << P_vector[i].oligo << endl;
    cout << "K: " << P_vector[i].K << " N1: " << P_vector[i].N1 << " N2: " << P_vector[i].N2 << endl;
    cout << "Ver_occ: " << P_vector[i].vertical_occurrences <<
    " Pval: " << P_vector[i].pvalue << endl;
  }
}

//Function where K, N1, N2 and T are calculated in order to obtain the p value
void PvalueClass::TKN1Calc(
    multimap<int, string>::iterator &it,
    unordered_map<string, KmerClass> &mappa_vector, unsigned int i) {
  
  //T is the number of sequences
  int T = sequences.size();

  //For each oligo in the multimap of vertical occurrences T, N1, N2 and K are calculated.

  //Remember N1 is the number of horizontal occurrences of oligo, T is the total number of sequences, N2 is the 
  //total number of oligos in all sequences minus N1 and K is the vertical occurrences of the oligo.
  // for (multimap<int, string>::iterator it = positions_occurrences.begin();
  //        it != positions_occurrences.end(); it++) {
      N1 = 0;
      vertical_occurrences = it->first;
      K = it->first;
      // Kv.emplace_back(K);
      oligo = it->second;
      // cout << "Oligo: " << oligos << endl;
      unordered_map<string, KmerClass>::iterator itBigMap =
          mappa_vector.find(oligo);
      unordered_map<string, KmerClass>::iterator itBigMap_rc =
          mappa_vector.find(reverse_oligo(oligo));
      if (itBigMap == mappa_vector.end()) {

        if (!itBigMap_rc->second.palindrome) {
          N1 = itBigMap_rc->second.orizzontal_count_rc;
        }
      } else {

        N1 = itBigMap->second.orizzontal_count;

      }
      //Calculation of total number of oligos in the multifasta
      tot_oligos = sequences.size() * len[i];

      //Calculation of N2
      N2 = tot_oligos - N1;

      //Using the gsl library for the hypergeometric p_value
      pvalue = gsl_cdf_hypergeometric_Q(K, N1, N2, T);
      if (pvalue == 0){
        pvalue = 1e-300;
      }
      // cout << pval << endl;
      //All the p_values are inserted in a vector

      // N1v.emplace_back(N1);
    // }
}
//N2 calculation with lambda (good for ss but not for ds)
// int PvalueClass::N2Calc(unordered_map<string, KmerClass> &mappa) {
//   int N2 = accumulate(begin(mappa), end(mappa), 0,
//                   [](unsigned int val,
//                      const unordered_map<string, KmerClass>::value_type &p) {
//                     return val + (p.second.orizzontal_count +
//                                   p.second.orizzontal_count_rc);
//                   });
//   return N2;
// }

//Debug output of pvalues
void PvalueClass::Dpvalues() {

    // for (unsigned int j = 0; j < pvalues.size(); j++){
    //   // cout << "P_value: " << pvalues[j] << endl;
    // }
}

//Fil the vector of N2s
// void PvalueClass::N2vFill(unordered_map<string, KmerClass> map_vector) {
//   N2v.push_back(N2Calc(map_vector));
// }

//Other debug functions

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

//Function for counting the horizontal and vertical occurrences for each oligo and 
//putting them inside a map.
void MapClass::CountOccurrences(string sequence, int k) {

  for (unsigned int i = 0; i < (sequence.size() - k + 1); i++) {
    string oligo = sequence.substr(i, k);
    string oligo_rc = reverse_oligo(oligo);
    unordered_map<string, KmerClass>::iterator it = mappa.find(oligo);
    unordered_map<string, KmerClass>::iterator it_rc = mappa.find(oligo_rc);
    if (it != mappa.end()) {

      it->second.orizzontal_count++;
      it->second.vertical_count[i]++;
      if(DS){
          it->second.orizzontal_count_rc++;
          it->second.vertical_count_rc[i]++;
      }
    } else if (it_rc != mappa.end()) {
      if(!it_rc->second.palindrome){
        it_rc->second.orizzontal_count_rc++;
        it_rc->second.vertical_count_rc[i]++;
        if(DS){
          it_rc->second.orizzontal_count++;
          it_rc->second.vertical_count[i]++;
        }
      }

    } else {
      KmerClass M;
      M.oligo = oligo, M.oligo_rc = oligo_rc;
      M.palindrome = bool(oligo == oligo_rc);
      M.orizzontal_count = 1, M.orizzontal_count_rc = 0;
      M.vertical_count.resize(((sequence.size() - k) + 1), 0);
      M.vertical_count_rc.resize(((sequence.size() - k) + 1), 0);
      M.vertical_count[i] = 1;
      if(DS){
          M.vertical_count_rc[i] = 1;
          M.orizzontal_count_rc = 1;
      }
      mappa.emplace(oligo, M);
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
      CountOccurrences(sequences[j], kv[i]);
    }
    mappa_vector.push_back(mappa);
    mappa.clear();
  }
}
void MapClass::VerticalMapVector() {
  
  for (unsigned int i = 0; i < kv.size(); i++) {

    for (unsigned int j = 0; j < len[i] ; j++) {
    multimap<int, string, greater<int> > pos;
    for (unordered_map<string, KmerClass>::iterator it = mappa_vector[i].begin();
         it != mappa_vector[i].end(); it++) {
        if(it->second.vertical_count[j] > 0){
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

//Debug function
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

//Debug function
void MapClass::DMainMapVector() {
  for (unsigned int i = 0; i < kv.size(); i++) {
    int sum = 0;
    for (unordered_map<string, KmerClass>::iterator it =
             mappa_vector[i].begin();
         it != mappa_vector[i].end(); it++) {
      cout << it->first << " " << it->second.vertical_count[0] << " "
           << it->second.orizzontal_count << "\t"
           << it->second.oligo_rc << " " << it->second.vertical_count_rc[0] << " " 
           << it->second.orizzontal_count_rc << " " 
           << " " << boolalpha << it->second.palindrome << noboolalpha << endl;
      sum = sum + it->second.vertical_count_rc[0] + it->second.vertical_count[0];
    }
  }
}

int main(int argc,  char **argv){
  
  MFASTA_FILE = argv[1];
  extract_sequences();
  Timer timer;
  for(unsigned int i = 0; i < kv.size(); i++){
    len.emplace_back(sequences[0].size() - kv[i] + 1);
  }
  MapClass M;
  vector<PvalueClass> P_vector;
  for (unsigned int i = 0; i < kv.size(); i++){
    for (unsigned int j = 0; j < len[i]; j++){
      cout << "Position: " << j << endl;
      for (multimap<int, string>::iterator it = M.vector_positions_occurrences[i][j].begin();
         it != M.vector_positions_occurrences[i][j].end(); it++) {
            PvalueClass P(it, M.mappa_vector[i], i);
            P_vector.push_back(P);
          }
          if (ordering){
            sort(begin(P_vector), end(P_vector), comp);
          }
          DVector(P_vector);
          P_vector.clear();
    }
  }
}
// unordered map <string, KmerClass>
