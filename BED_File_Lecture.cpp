#include "BED_File_Lecture.h"


int main(int argc, char *argv[]){

	if(argc == 1){             //If arguments number is 1 means that no input file has been inserted - display help
		display_help();
	}

	command_line_parser(argc, argv); //parser function called to handle aguments
	
	
	ifstream myfile (BED_FILE); //Opening file in lecture mode// it must not be hard-coded!!!!
	vector<genomic_position> GEP;	 //defining vector of genomic_position datas
	string line; 			//defining line string
	string token;			//defining token string
	genomic_position prova();	//initialization of class prova of genomic_position type using the default constructor 	

	int n_line = 0;			//line counter initialization
        while(getline(myfile,line)){  //reading input file line by line with getline function

	 vector<string> x;		//defining string vector x

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
			cerr << "WARNING: the line " << n_line << " is omitted because starting coordinates > end coordinates, please check your BED file!" << "\n";
			//if flag is not 1 means that the current line has starting coordinate < end coordinate: PRINT WARNING!		
		}

		n_line = n_line + 1;			//pass to next line 
	}

	for (int i=0; i<GEP.size(); ++i){    // from 0 to GEP vector length
		cout<< GEP[i].get_char() << "\n" << GEP[i].get_start() << "\n" << GEP[i].get_end() << "\n" << GEP[i].get_flag() << "\n\n"; //control print
		}

}


void command_line_parser(int argc, char **argv){


for(int i = 1; i < argc; i++){
	
	string buf = argv[i];

	if(buf == "--help" || buf == "-h"){

	display_help();

	}

	if(buf == "--BED" || buf == "-B"){

		if(i < argc - 1){
			
			BED_FILE = argv[++i];

		continue;
        	}
	}
	
	

}
	

}

void display_help()
{
	cerr << "\n --help: show this message" << endl;
	cerr << "\n --BED -B <file_bed>: input bed file" << endl;
	cerr << endl;
	
	exit(EXIT_SUCCESS);
}


