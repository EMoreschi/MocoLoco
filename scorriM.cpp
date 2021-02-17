#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <dirent.h>
#include <cmath>
#include <valarray>
#include <vector>
#include <algorithm>
#include <typeinfo> 
#include <sstream>

using namespace std;

void scorrimento(vector<vector<int>> mat, string seq, int p, int length, vector<int>&);

int main(){
	
	vector<int> oligo;
//	int matrix [4][4] = {3,1,4,0,1,3,8,2,0,3,6,1,2,2,7,1};
	vector<vector<int>> matrix = {{3,1,4,0},{1,3,8,2},{0,3,6,1},{2,2,7,1}};
	string seq = "ATAGTCCACGATCACTGGGA";
	int p = 0;
	int l = 16;
	int end = 4;

	

	scorrimento(matrix, seq, p, l, oligo);

	for(int i=0; i<oligo.size(); i++){
		cout << oligo[i] << " ";
	}
}

void scorrimento(vector<vector<int>> mat, string seq, int p, int l, vector<int> &oligo){
		
	int sum_oligo = 0;
	
	if(p <= seq.size() -4 ) {

	for(int i=0; i<4; i++){

			switch(seq[i+p]){

				case 'A':
				       
					sum_oligo = sum_oligo + mat[0][i];
					break;

				case 'C':
				       
					sum_oligo = sum_oligo + mat[1][i];
					break;

				case 'G':
				       
					sum_oligo = sum_oligo + mat[2][i];
					break;

				case 'T':
				       
					sum_oligo = sum_oligo + mat[3][i];
					break;

		}
	}
	
	oligo.emplace_back(sum_oligo);
	scorrimento(mat, seq, p+1, l, oligo);
	}

}



