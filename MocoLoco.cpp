#include "MocoLoco.h"

int main(int argc, char *argv[]){


if(argc == 1){             //If arguments number is 1 means that no input file has been inserted - display help
	display_help();
}

command_line_parser(argc, argv); //parser function called to handle aguments

GEP_objects_creation(BED_FILE, TWOBIT_FILE);


}

void GEP_objects_creation(const char* Bed_file, const char* Twobit_file){
	
	ifstream myfile (Bed_file); //Opening file in lecture mode// it must not be hard-coded!!!!
	const char * chrom;
	TwoBit * tb;
	tb = twobit_open(Twobit_file);
	
	//if (tb == NULL) {
	//	
	//	fprintf(stderr, "Failed to open: %s\n", TWOBIT_FILE);
	//	return EXIT_FAILURE;
	//}

	vector<genomic_position> GEP;	 //defining vector of genomic_position datas
	string line; 			//defining line string
	string token;			//defining token string
	genomic_position new_class();	//initialization of class prova of genomic_position type using the default constructor 	
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

		int s = stoul(x[1])-1;  //The word corrisponding to start coordinate converted from string to  int
		int e = stoul(x[2])-1;	//The word corrisponding to end coordinate converted from string to  int -1 because ucsc count from 

		genomic_position new_class(x[0],s,e,parameter);
		
		if(new_class.flag == 1){	//CONTROL: if flag is 1 means that the current line has starting coordinate > end coordinate, so it is correct
		        chrom = &*new_class.chr_coord.begin();
			new_class.sequence = twobit_sequence(tb,chrom,new_class.start_coord,new_class.end_coord-1);
			GEP.push_back(genomic_position{new_class});	//put the class prova in GAP (vector of classes of type genomic_position)
		}

		else {		
			cerr << "WARNING: the line " << n_line << " is omitted because starting coordinates > end coordinates, please check your BED file!" << "\n";
			//if flag is not 1 means that the current line has starting coordinate < end coordinate: PRINT WARNING!		
		}

		n_line = n_line + 1;			//pass to next line 
	}



	for (int i=0; i<GEP.size(); ++i){    // from 0 to GEP vector length
                 cout<< ">" << GEP[i].chr_coord <<":"<< GEP[i].start_coord << "-" << GEP[i].end_coord << "\n";
                 //cout<< GEP[i].chr_coord <<"\t"<< GEP[i].start_coord << "\t" << GEP[i].end_coord << "\n";
		 cout << GEP[i].sequence<<"\n";
		//cout<< GEP[i].sequence; //control print
	}

}

void command_line_parser(int argc, char **argv){

	int control_bed = 0;
	int control_twobit = 0;
	int control_p = 0;

	for(int i = 1; i < argc; i++){

		string buf = argv[i];

		if(buf == "--help" || buf == "-h"){

			display_help();

		}

		if(buf == "--BED" || buf == "-B"){

			if(i < argc - 1){

				BED_FILE = argv[++i];
				control_bed = 1;
				continue;
			
			}
		}

		if(buf == "--param" || buf == "-p"){

			if(i < argc - 1){

				parameter = stoi(argv[++i]);
				control_p = 1;
				continue;
			}
		}
		if(buf == "--twobit" || buf == "-tb"){

			if(i < argc - 1){

				TWOBIT_FILE = argv[++i];
				control_twobit = 1;
				continue;
			}
		}


	}

	if(control_bed == 0 && control_twobit == 0){
		
		cout << "BED file and TwoBit file missed or wrong parameters calling!\n";
		cout << "Please insert as input a BED file using -B or --BED annotation before it.\n";
		cout << "Please insert as input a Twobit file using -tb or --twobit annotation before it.\n";
		
		exit(EXIT_SUCCESS);
	}
		
	if(control_bed == 0){
		
		cout << "Wrong BED file immission or BED file missed!\n ";
		cout << "Please insert as input a BED file using -B or --BED annotation before it.\n";
		
		exit(EXIT_SUCCESS);
	}

	if(control_twobit == 0){

		cout << "Wrong Twobit file immission or Twobit file missed!\n";
		cout << "Please insert as input a Twobit file using -tb or --twobit annotation before it.\n";

		exit(EXIT_SUCCESS);

	}

	if(control_p == 0){

		cout << "WARNING: Sequence length parameter not inserted --> Default parameter is 150.\n";
		cout << "To put a parameter as input write -p or --param before parameter value.";
	
	}
}



void display_help()
{
	cerr << "\n --help: show this message" << endl;
	cerr << "\n --BED -B <file_bed>: input bed file" << endl;
	cerr << endl;

	exit(EXIT_SUCCESS);
}

