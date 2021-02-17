#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <typeinfo> 
#include <sstream>

using namespace std;

void scorrimento(int mat[4][4], char seq[20], int p, int length, vector<int>&);

int main(){
	
	vector<int> oligo;
	int matrix [4][4] = {3,1,4,0,1,3,8,2,0,3,6,1,2,2,7,1};
	char seq [20] = {'A','T','A','G','T','C','C','A','C','G','A','T','C','A','C','T','G','G','G','A'};
	int p = 0;
	int l = 16;
	int end = 4;

	scorrimento(matrix, seq, p, l, oligo);

	for(int i=0; i<oligo.size(); i++){
		cout << oligo[i] << " ";
	}
}

void scorrimento(int mat[4][4], char seq[20], int p, int l, vector<int> &oligo){
		
	int sum_oligo = 0;
	
	if(p <= l){

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



