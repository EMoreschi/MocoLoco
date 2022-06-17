#include "../TwoBit/twobit.c"
#include "../TwoBit/twobit.h"
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
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
    string chr;
    unsigned int start, end;
  };
  vector<bed_s> bed_v;

public:
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
  vector<vector<int>> NLogMatrix;
  vector<vector<int>> INLogMatrix;
  void Norm(int, vector<double>, vector<vector<double>>);
  vector<double> ColSum(vector<vector<double>>);
  void MakeLog();

public:
  matrix_c(vector<vector<double>> JM, int psdcount) {
    vector<double> col_sum;
    col_sum = ColSum(JM);
    Norm(0, col_sum, JM);
  }
};
