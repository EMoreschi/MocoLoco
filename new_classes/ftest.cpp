#include "ftest.h"

unsigned int half_length = 150, overhead = 25;

int main() {
  /// READING BED
  string BED_FILE = "../Test_Bed/nfy_k562_hg19.bed";
  string TWOBIT_FILE = "../Genomes/hg19/hg19.2bit";
  string JASPAR_FILE = "../Jaspar_2020/MA0006.1.jaspar";
  TwoBit *tb;
  char **result;
  tb = twobit_open(TWOBIT_FILE.c_str());
  bed_c B(BED_FILE, tb);
  JR_c J(JASPAR_FILE);
  matrix_c mat(J.RJM());
  vector<vector<double>> Gscore;

  for (unsigned int i = 0; i < B.bed_v.size(); i++) {
    score_c score(mat.RMatrix(), B.bed_v[i].seq);
    Gscore.emplace_back(score.seq_scores);
  }
  for (unsigned int i = 0; i < B.bed_v.size(); i++) {
    PrintVector(Gscore[i]);
  }
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
      bed_in.seq =
          twobit_sequence(tb, bed_in.chr.c_str(), bed_in.start, bed_in.end - 1);

      bed_v.push_back(bed_in);
    }
  }
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
vector<vector<double>> matrix_c::RMatrix() { return LogMatrix; }
vector<vector<double>> matrix_c::RICMatrix() { return ICLogMatrix; }

// void matrix_c::Norm(std::vector<double> col_sum,
// std::vector<std::vector<double>>& ma){
void matrix_c::Norm(double psdcount, vector<double> col_sum,
                    vector<vector<double>> &JM) {
  for (unsigned int x = 0; x < JM.size(); x++) {
    for (unsigned int i = 0; i < JM[x].size(); i++) {
      double *ptr = &JM[x][i];
      *ptr = (*ptr / col_sum[i] + psdcount);
    }
  }
}

vector<double> matrix_c::ColSum(vector<vector<double>> &JM) {
  vector<double> col_sum;
  for (unsigned int i = 0; i < JM[0].size(); i++) {
    double sum = 0;
    for (unsigned int j = 0; j < 4; j++) {
      sum += JM[j][i];
    }
    col_sum.emplace_back(sum);
  }
  return col_sum;
}

void matrix_c::MakeLog(vector<vector<double>> &JM) {
  for (unsigned int x = 0; x < JM.size(); x++) {
    for (unsigned int i = 0; i < JM[x].size(); i++) {
      double *ptr = &JM[x][i];
      *ptr = log(JM[x][i]);
    }
  }
  LogMatrix = JM;
}

void score_c::FindMinMax(vector<vector<double>> &matrix) {
  vector<double> max_sum;
  for (unsigned int i = 0; i < matrix[0].size(); i++) {
    vector<double> column;
    for (unsigned int j = 0; j < matrix.size(); j++) {
      column.emplace_back(matrix[j][i]);
      // cout << matrix[j][i] << endl;
    }
    min_sum.emplace_back(*min_element(column.begin(), column.end()));
    max_sum.emplace_back(*max_element(column.begin(), column.end()));
  }

  double min = accumulate(min_sum.begin(), min_sum.end(), 0.0);
  double max = accumulate(max_sum.begin(), max_sum.end(), 0.0);
  minmax_v.emplace_back(min);
  minmax_v.emplace_back(max);
}

vector<double> score_c::Shifting(vector<vector<double>> &matrix,
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

        sum_scores += min_sum[i];
        break;
      }
    }

    // Best score normalization with normalization formula
    sum_scores = 1 + (sum_scores - minmax_v[1]) / (minmax_v[1] - minmax_v[0]);
    seq_scores.emplace_back(sum_scores);
    // cout << seq_scores[s_iterator];
    //  cout << *seq_scores[s_iterator] << endl;
    //  The total score of an oligo is saved into an oligo_scores vector
  }
  return seq_scores;
}

void matrix_c::InverseMa() {
  // PROFILE_FUNCTION();
  ICLogMatrix = LogMatrix;
  reverse(ICLogMatrix.begin(), ICLogMatrix.end());

  for (int i = 0; i < 4; i++) {

    reverse(ICLogMatrix[i].begin(), ICLogMatrix[i].end());
  }
}
