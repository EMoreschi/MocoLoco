#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>


class genomic_position { //creation public class of genomic_position type        
	        public:
		std::string chr;
		int start;
		int end;
};

int main(){
	std::ifstream myfile ("BED_di_prova.BED"); //open file in lecture mode// it must not to be hard-coded!!!!
	std::vector<genomic_position> GEP; // defining vector
	std::string chr; //defining string
        unsigned int start;  // defining start coordinate variable
        unsigned int end;   // defining end coordinate variable
	std::string line; //defining line
 
        while(std::getline(myfile,line)){  //reading input file line by line with getline function
  	     if (line.empty())		   //CONTROL: if line is empty pass to next line
		     continue;
	      if (line[0] == '#')	//CONTROL: if line start with # (possible headers or comments) pass to next line
		     continue;

              std::istringstream my_stream(line); //istringstream function to separe each lines word by word	
	      my_stream	>> chr >> start >> end;  
	      GEP.push_back(genomic_position{chr, start, end}); //allocating in GEP vector a space to put in object of genomic_position class

			}

//	while(myfile >> chr >> start >> end) {   // until it reads object
//	       GEP.push_back(genomic_position{chr,start,end});	
//	}
	for (int i=0; i<GEP.size(); ++i)    // from 0 to GEP vector length
	       std::cout<< GEP[i].end <<" "<< GEP[i].start << " " << GEP[i].chr << "\n"; //control print 
	return 0;
}




