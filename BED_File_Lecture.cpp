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
		bool flag;
		
		

		public:			//definition public field
		

		//CONSTRUCTORS:
		genomic_position(){	//default constructor
		
		chr = "";
		start = 0;
		end = 0;
		flag = 0;	
		} 

		genomic_position(string chr_coord, unsigned int start_coord, unsigned int end_coord) {  //constructor given member values
			
			
			chr = chr_coord;
			start = start_coord;
			end = end_coord;
			flag = flag_control(start_coord, end_coord);
		}

		//FUNCTIONS:
		string get_char(){	//print char function
		return chr;
		}

		int get_start(){	//print start coordinate function
		return start;
		}

		int get_end(){		//print end coordinate funcion
		return end;
		}
		
		bool flag_control(unsigned int st, unsigned int en){ //function which controls that start coordinates are < then end coordinates

		if(st > en){		//if start coordinates are > then end coordinates, flag is setted to 0
			return 0;
		}
		else{ return 1;}
		}	

		bool get_flag(){		//print end coordinate funcion
		return flag;
		}
};



int main(){
	ifstream myfile ("BED_di_prova.BED"); //Opening file in lecture mode// it must not be hard-coded!!!!
	vector<genomic_position> GEP;	 //defining vector of genomic_position datas
	vector<string> x;		//defining string vector x	
	string line; 			//defining line string
	string token;			//defining token string
	genomic_position prova();	//initialization of class prova of genomic_position type using the default constructor 	

	int n_line = 0;			//line counter initialization
        while(getline(myfile,line)){  //reading input file line by line with getline function
  	     if (line.empty())		   //CONTROL: if line is empty pass to next line
		     continue;
	      if (line[0] == '#')	//CONTROL: if line starts with # (possible headers or comments) pass to next line
		     continue;

		istringstream my_stream(line); //istringstream function to split each line word by word	
	      
		while(my_stream >> token){	//put every word in token string while words in the line are not finished
		
		x.push_back(string{token});	//put every word in string vector called x until the words in the line are finished	
		}
		unsigned int s = stoul(x[1]);  //The word corrisponding to start coordinate converted from string to unsigned int
		unsigned int e = stoul(x[2]);	//The word corrisponding to end coordinate converted from string to unsigned int
		genomic_position prova(x[0],s,e);  //modifiyng genomic_position class prova (inizialized as default before)
						   //with the corrisponding data from the current file line
		
		if(prova.get_flag() == 1){	//CONTROL: if flag is 1 means that the current line has starting coordinate > end coordinate, so it is correct

			GEP.push_back(genomic_position{prova});	//put the class prova in GAP (vector of classes of type genomic_position)
		}
		else {		
			cout << "WARNING: the line " << n_line << " is omitted because starting coordinates > end coordinates, please check your BED file!" << "\n";
			//if flag is not 1 means that the current line has starting coordinate < end coordinate: PRINT WARNING!		
		}

		x.clear();				//Restore empty x vector
		n_line = n_line + 1;			//pass to next line 
	}

	for (int i=0; i<GEP.size(); ++i){    // from 0 to GEP vector length
		cout<< GEP[i].get_char() << "\n" << GEP[i].get_start() << "\n" << GEP[i].get_end() << "\n" << GEP[i].get_flag() << "\n\n"; //control print
		}

}


