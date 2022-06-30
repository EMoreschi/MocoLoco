#include "../TwoBit/twobit.c"
#include "../TwoBit/twobit.h"
#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <math.h>
#include <numeric>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

using namespace std;

template <typename T> void PrintVector(const T &t) {
  copy(t.cbegin(), t.cend(),
       ostream_iterator<typename T::value_type>(cout, ", "));
  cout << '\n';
}

template <typename T> void PrintMatrix(const T &t) {
  for_each(t.cbegin(), t.cend(), PrintVector<typename T::value_type>);
}

class bed_c {

  struct bed_s {
    string chr, seq;
    unsigned int start, end;
  };

public:
  vector<bed_s> bed_v;
  vector<bed_s> Rbed_v();
  void ReadB(string, TwoBit *);
  bed_c(string BED_FILE, TwoBit *tb) { ReadB(BED_FILE, tb); };
};

class JR_c {
private:
  vector<vector<double>> JM;
  void ReadJ(string);

public:
  vector<vector<double>> RJM();

  JR_c(string line) { ReadJ(line); };
};

class matrix_c {
  vector<vector<double>> LogMatrix;
  vector<vector<double>> ICLogMatrix;
  void Norm(double, vector<double>, vector<vector<double>> &);
  vector<double> ColSum(vector<vector<double>> &);
  void MakeLog(vector<vector<double>> &);
  void InverseMa();

public:
  vector<vector<double>> RMatrix();
  vector<vector<double>> RICMatrix();
  matrix_c(vector<vector<double>> JM) {
    const double pseudoc = 0.01;
    Norm(pseudoc, ColSum(JM), JM);
    PrintMatrix(JM);
    cout << '\n';
    Norm(0, ColSum(JM), JM);
    PrintMatrix(JM);
    cout << '\n';
    MakeLog(JM);
    PrintMatrix(LogMatrix);
    cout << '\n';
    InverseMa();
    PrintMatrix(ICLogMatrix);
  }
};

class score_c {
  vector<unsigned int> score_centered_pos;
  double global_scores;
  vector<double> min_sum;
  vector<double> minmax_v;

public:
  vector<double> seq_scores;
  void FindMinMax(vector<vector<double>> &);
  vector<double> Shifting(vector<vector<double>> &, string &);
  score_c(vector<vector<double>> log_matrix, string seq) {
    FindMinMax(log_matrix);
    seq_scores = Shifting(log_matrix, seq);
    // cout << *seq_scores;
  }
};
