## **INSTRUCTION** and **REPORT** MocoLoco.cpp and MocoLoco.h files project
<br>  
<br>  

### The program has **two different ways to be executed**.

   1. BED, TwoBit and Jaspar files in input.
   2. Multifasta file in input.

<br>  

### To **compile** the program the command needed in both cases is:

     	g++ -o test -Wall Mocoloco.cpp -lgsl -lgslclbas
<br>  

### To **execute** it 2 different ways have been implemented:

1. BED,TwoBit,Jaspar input:
       
       ./test -b <file.bed>  -p <number>(optional) -t <file.2bit> -j <file.JASPAR> -k <number, number, ...>(optional) -n <number> (optional) -ss(optional)`

2. Multifasta input:
       
       ./test -m <file_multifasta.fa>  -p <number>(optional) -k <number, number, ...>(optional) -ss(optional)
<br>  

The output is composed by 3 different files .txt for each k inserted as input. We can have different file compositions, following if the analysis has made on Double Strand or Single Strand.
		
    A. The first one is the k-mers count in all the sequences extractred from bed coordinates inserted, ordered from the most relevant to the less.
       		
       A1. If we are in DS analysis the output file contains oligo sequence + its occurrences in FWD strand and RC sequence + its occurrences in reverse strand (the two occurrences must be equal).
       The name of the file is k-mers_occurrences_twobit_jaspar_bed_DS.txt.

       A2. If we are in SS analysis the output file contains oligo sequence + its occurrences.
       The name of the file is k-mers_occurrences_twobit_jaspar_bed_SS.txt.
       
    B. The second one is a file with the k-mers count for each position of sequences.
       The file contains information, for each position, of the top -n kmers ranked by occurrences. 
       
       B1. For DS analysis these information are:
       
    • The position.
    • The rank.
    • The oligo occurrences in FWD strand.
    • The reverse complement occurrences in reverse strand.
    • The sum of their occurrences.
    • The occurrences of reverse complement in forward strand.
    • The occurrences of the oligo in reverse strand.
    • The sum of their occurrences.
    • TRUE if the oligo is palindrome, FALSE if not.
    • The frequency of occurrences at that position.
       The filename is k-mers_positional_occurrences_twobit_jaspar_bed_DS.txt
       
       B2. For SS analysis these information are:
       
    • The position.
    • The rank.
    • The oligo occurrences in FWD strand.
    • The occurrences of reverse complement on FWD strand.
    • TRUE if the oligo is palindrome, FALSE if not.
    • The frequency of oligo occurrences at that position.
       The filename is k-mers_positional_occurrences_twobit_jaspar_bed_SS.txt
       	
    C. The third is a file with the sum of the occurrences of the first top -n oligos for each position and the frequences of these sums for each position.
       The output structure for DS and SS analysis is the same but obviously the data has come from different processes.
       The name of the file is k-mers_Topn_sum_and_frequence_DS.txt (For DS)
       The name of the file is k-mers_Topn_sum_and_frequence_SS.txt (For SS)
	
There are also two others output file, one in .bed format and the other in .fasta format.
These output are created only if the program is executed with bed/twobit/jaspar input, and they are:

    D. A .bed file with the genomic coordinates of the sequences analized by MocoLoco.cpp tool. The file contains, for each sequence, information about:

    • The Chromosome
    • The Starting Coordinate
    • The Ending Coordinate
	The filename is twobit_jaspar_bed.bed.

    E. A .fasta file which contains the sequences analyzed by MocoLoco tool in fasta format.
       For every sequence it contains:

    • Header (fasta format) with genomic coordinates.
    • The ATCG fasta sequence extracted from TwoBit file, centered following the coordinates bed in input and as long as indicated by the parameter p (length = p x 2).
	The filename is twobit_jaspar_bed.fasta.

Finally, there is another supplementary output file called k-mers_p_value_control_parameters.txt (one for each k inserted as input), which contains all the parameter values used to calculate Oligos P-values.
This file contains information about:
	
    • Position in sequence.
    • Rank in the position.
    • Oligo string.
    • K parameter (Occurrences of the oligo in that position).
    • N1 parameter (Occurrences of the oligo in every position).
    • N2 parameter (Total possible oligos - N1).
    • T parameter (Number of sequences).
    • P_value.

MOCOLOCO WORKFLOW

Firstly the file MocoLoco.h is included with all the libraries, structures, global variables and functions in it.
The funcrion Main read the arguments in input and with the function command_line_parser assigns input parameters to the right global variables. This function is also able to control the validity and the existance of parameters and it is is trained to report errors in case of incorrect arguments, followed by help to assist users.
Once the global variables are assigned correctly, the tool is ready to process the informations.
The first function is GEP_path, which is able to distinguish if input are in bed/twobit/jaspar or multifasta format and, cosequently, it chooses the right path.


BED/TWOBIT/JASPAR Pathway:

If MocoLoco recognizes a bed/twobit/jaspar input it immediately creates a coordinator class C, which handles and connects the arguments.
this class is the connecting hub between the bed, twobit and jaspar files. It allows the creation of bed classes by extracting genomic sequences from the twobit file; it also allows to create and manage a jaspar matrix class and finally to connect them to each other thanks to creation of an additional oligo class.
The first coordinator class function is GEP_creation:

	This function has the task of opening the input bed file, reading its contents and finally, for 	each genomic coordinate readed in the bed file, creating a bed class.
	Each bed class created is then inserted into a vector of bed classes called GEP (GEnomic 	Position).  
	An object of bed class is created for each genomic coordinate in the bed file. The class 	constructor provides that the object can be constructed in two different ways: one for a 		bed/twobit/jaspar input and one for a multifasta input. Now the first case is analyzed, so the 	constructor used is the one with passage of 4 parameters (p = half window legth, line = 	current bed line tb = Twobit file, n_line = current bed line number).
	Once the constructor calling and the parameters passing, the bed class is able to handle the 	data and it does it calling 4 functions:

    1. read_line:
       This function has the task of taking the input line (from bed file) and splits it in words, extracting and saving information about chromosome, starting coordinate and end coordinate.
    2. flag_control:
       This is a controlling function. Its job is to check if end coordinate is greater than start coordinate as it should be. If this is confirmed the flag is set to 1, otherwise the flag is set to 0.
    3. centering_function:
       It has the task of centering the genomic coordinates in the middle of start/end window.
       Then the function reassigns the start coordinate calculating it making 
       
