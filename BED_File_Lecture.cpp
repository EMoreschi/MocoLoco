#include "BED_File_Lecture.h"


int main(int argc, char *argv[]){
	
	
	ifstream myfile (argv[1]); //Opening file in lecture mode// it must not be hard-coded!!!!
	vector<genomic_position> GEP;	 //defining vector of genomic_position datas
	vector<string> bed_words;		//defining string vector x	
	string line; 			//defining line string
	string token;			//defining token string
	genomic_position bed_line();	//initialization of class bed_line of genomic_position type using the default constructor 	

	int n_line = 0;			//line counter initialization
        while(getline(myfile,line)){  //reading input file line by line with getline function
  	     if (line.empty())		   //CONTROL: if line is empty pass to next line
		     continue;
	      if (line[0] == '#')	//CONTROL: if line starts with # (possible headers or comments) pass to next line
		     continue;

		istringstream my_stream(line); //istringstream function to split each line word by word	
	      
		while(my_stream >> token){	//put every word in token string while words in the line are not finished
		
		bed_words.push_back(string{token});	//put every word in string vector called x until the words in the line are finished	
		}
		unsigned int s = stoul(bed_words[1]);  //The word corrisponding to start coordinate converted from string to unsigned int
		unsigned int e = stoul(bed_words[2]);	//The word corrisponding to end coordinate converted from string to unsigned int
		genomic_position bed_line(bed_words[0],s,e);  //modifiyng genomic_position class bed_line (inizialized as default before)
						   //with the corrisponding data from the current file line
		
		if(bed_line.get_flag() == 1){	//CONTROL: if flag is 1 means that the current line has starting coordinate > end coordinate, so it is correct

			GEP.push_back(genomic_position{bed_line});	//put the class bed_line in GAP (vector of classes of type genomic_position)
		}
		else {		
			cerr << "WARNING: the line " << n_line << " is omitted because starting coordinates > end coordinates, please check your BED file!" << "\n";
			//if flag is not 1 means that the current line has starting coordinate < end coordinate: PRINT WARNING!		
		}

		bed_words.clear();				//Restore empty x vector
		n_line = n_line + 1;			//pass to next line 
	}

	for (int i=0; i<GEP.size(); ++i){    // from 0 to GEP vector length
		cout<< GEP[i].get_char() << "\n" << GEP[i].get_start() << "\n" << GEP[i].get_end() << "\n" << GEP[i].get_flag() << "\n\n"; //control print
		}

}


