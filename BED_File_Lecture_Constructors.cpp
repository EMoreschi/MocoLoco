#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <typeinfo> 
#include <sstream>
using namespace std;

class genomic_position { //creation public class of genomic_position type        
	        
		private:	//definition private field
		string chr;
		unsigned int start;
		unsigned int end;
		

		public:			//definition public field
		genomic_position(){	//default constructor
		
		chr = "";
		start = 0;
		end = 0;
		
			
		} 

		genomic_position(string chr_coord, unsigned int start_coord, unsigned int end_coord) {  //constructor given member values

			chr = chr_coord;
			start = start_coord;
			end = end_coord;
		
		}

		string get_char(){	//print char function
		return chr;
		}

		int get_start(){	//print start coordinate function
		return start;
		}

		int get_end(){		//print end coordinate funcion
		return end;
		}
};



int main(){
	ifstream myfile ("BED_di_prova.BED"); //Opening file in lecture mode// it must not be hard-coded!!!!
	vector<genomic_position> GEP;	 //defining vector of genomic_position datas
	vector<string> x;		//defining string vector x	
	string line; 			//defining line string
	string token;			//defining token string
	genomic_position prova();	//initialization of class prova of genomic_position type using the default constructor
 	int i = 0;
	int n_arg = 3;			//number of columns in file (feature, args)

        while(getline(myfile,line)){  //reading input file line by line with getline function
  	     if (line.empty())		   //CONTROL: if line is empty pass to next line
		     continue;
	      if (line[0] == '#')	//CONTROL: if line starts with # (possible headers or comments) pass to next line
		     continue;

		istringstream my_stream(line); //istringstream function to split each line word by word	
	      
		while(my_stream >> token){	//put every word in token string while words in the line are not finished
		
		x.push_back(string{token});	//put every word in string vector called x until the words in the line are finished	
		}
		unsigned int s = stoul(x[i+1]);  //The word corrisponding to start coordinate converted from string to unsigned int
		unsigned int e = stoul(x[i+2]);	//The word corrisponding to end coordinate converted from string to unsigned int
		genomic_position prova(x[i],s,e);  //modifiyng genomic_position class prova (inizialized as default before)
						   //with the corrisponding data from the current file line
		
		GEP.push_back(genomic_position{prova});	//put the class prova in GAP (vector of classes of type genomic_position)
		
		i = i + n_arg;	//increase i of narg value to slide in the correct way the string vector x 
	}

	for (int i=0; i<GEP.size(); ++i){    // from 0 to GEP vector length
		cout<< GEP[i].get_char() << "\n" << GEP[i].get_start() << "\n" << GEP[i].get_end() << "\n\n"; //control print
		}

}



