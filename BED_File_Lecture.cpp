#include "BED_File_Lecture.h.cpp"


int main(int argc, char *argv[]){

	if(argc == 1){             //If arguments number is 1 means that no input file has been inserted - display help
		display_help();
	}

	command_line_parser(argc, argv); //parser function called to handle aguments
	
	
	ifstream myfile (argv[1]); //Opening file in lecture mode// it must not be hard-coded!!!!
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
	cerr << "\nSynopsis: aScan --rna rnaseq_bamfile --vcf vcf_file --gtf gtf_file [-nt transcript/gene_correspondence_file] [-p thread_num] [--filter filter_key_word]" << endl;
	cerr << "\n--rna rnaseq_bamfile\nBAM file containing RNA-Seq reads mapped to a reference genome. The BAM file does not need to be sorted. Please notice that aScan does not currently support reads mapped to the reference transcriptome." << endl;
	cerr << "\n--vcf vcf_file\nVCF file with genomic variants from the same individual as the RNA-Seq data. Please notice that aScan currently does not support multi-VCF files. aScan considers only single nucleotide\nsubstitutions. aScan can work both with phased or unphased VCF files." << endl;
	cerr << "\n--gtf gtf_file GTF file with reference transcripts annotation.\nPlease notice that aScan will consider only \"exon\" named features. The attribute column of the GTF file must include a transcript_id value for each exon." << endl;
	cerr << "\n-nt, --nametable transcript/gene_correspondence_file.\nA tab-separated tabular file associating transcript IDs (first column) to gene names or gene IDs (second column)." << endl;
	cerr << "\n-p threadnum\nThe number of threads to use (default 1). While aScan supports multithreading, the execution speed bottleneck is usually represented by the drive read bandwidth,\nmeaning that normally the benefits of using more than 2 or 3 threads are minimal." << endl;
	cerr << "\n--filter filter_key_word\nUsing this option, all the variants flagged with filter_key_word in the FILTER field of the VCF file will be ignored. Use this option to discard low-quality variant calls." << endl;
	
	cerr << endl;
	
	exit(EXIT_SUCCESS);
}


