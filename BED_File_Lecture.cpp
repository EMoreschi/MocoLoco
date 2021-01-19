#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>


class genomic_position { //creazione classe pubblica di tipo genomic position          
	        public:
		std::string chr;
		int start;
		int end;
};

int main(){
	std::ifstream myfile ("imp_file.txt"); //apro il file in lettura// non deve essere hard-coded!!!!
	std::vector<genomic_position> GEP; //definisco vettore
	std::string chr; //definisco stringa
        unsigned int start;  // definisco in start
        unsigned int end;   // definisco end	
	std::string line; //definisco line
 
        while(std::getline(myfile,line)){  //leggo il file in input riga per riga con funzione getline
  	     if (line.empty())		   //se linea vuota passo alla successiva
		     continue;
	      if (line[0] == '#')	//se linea inizia con # (possibili header o commenti) passo alla linea successiva
		     continue;

              std::istringstream my_stream(line); //funzione istringstream per serparare, per ogni linea, parola per parola 	
	      my_stream	>> chr >> start >> end;  
	      GEP.push_back(genomic_position{chr, start, end}); //Alloco nel vettore GEP uno spazio per inserire oggetti di classe genomic_position

			}

//	while(myfile >> chr >> start >> end) {   // fino a quando leggo gli oggetti
//	       GEP.push_back(genomic_position{chr,start,end});	
//	}
	for (int i=0; i<GEP.size(); ++i)    // da 0 fino a tutta la lunghezza del vettore GEP
	       std::cout<< GEP[i].end <<" "<< GEP[i].start << " " << GEP[i].chr << "\n"; //stampa per controllo 
	return 0;
}




