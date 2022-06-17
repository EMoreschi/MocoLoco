#include "ftest.h"
#include <fstream>
#include <unistd.h>

unsigned int half_length = 150, overhead = 25;

int main() {
  /// READING BED
  string BED_FILE = "../Test_Bed/nfy_k562_hg19.bed";
  string TWOBIT_FILE = "../Genomes/hg19/hg19.2bit";
  string JASPAR_FILE = "../Jaspar_2020/MA0060.1.jaspar";
  TwoBit *tb;
  char **result;
  tb = twobit_open(TWOBIT_FILE.c_str());
  bed_c B(BED_FILE, tb);
  JR_c J(JASPAR_FILE);
  matrix_c mat(J.RJM(), 0);
  return 0;
}

void bed_c::ReadB(string BED_FILE, TwoBit *tb) {
  ifstream in(BED_FILE);
  string line;
  while (getline(in, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    istringstream mystream(line);
    bed_s bed_in;
    mystream >> bed_in.chr >> bed_in.start >> bed_in.end;
    if (bed_in.start > bed_in.end) {
      cout << "error " << endl;
    } else {
      unsigned int center = (bed_in.start + bed_in.end) / 2;
      bed_in.start = center - half_length;
      bed_in.end = center + half_length + overhead;
      bed_v.push_back(bed_in);
    }
  }
  //  for (unsigned int i = 0; i < bed_v.size(); i++) {
  //    cout << bed_v[i].chr << '\t' << bed_v[i].start << '\t' << bed_v[i].end
  //         << endl;
  //  }
}
void JR_c::ReadJ(string JASPAR_FILE) {

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
      JM.emplace_back(scores_line);
    }
  }
  file.close();
}

vector<vector<double>> JR_c::RJM() { return JM; }

// void matrix_c::Norm(std::vector<double> col_sum,
// std::vector<std::vector<double>>& ma){
void matrix_c::Norm(double psdcount, vector<double> col_sum,
                    vector<vector<double>> JM) {
  for (unsigned int x = 0; x < JM.size(); x++) {
    for (unsigned int i = 0; i < JM[x].size(); i++) {
      JM[x][i] = JM[x][i] / col_sum[i] + psdcount;
    }
  }
}

vector<double> matrix_c::ColSum(vector<vector<double>> JM) {
  vector<double> col_sum;
  double sum = 0;
  for (unsigned int i = 0; i < JM[0].size(); i++) {
    for (unsigned int j = 0; j < 2; j++) {
      sum += JM[j][i];
    }
    col_sum.emplace_back(sum);
    sum = 0;
  }
  return col_sum;
}

void matrix_c::MakeLog() {}
