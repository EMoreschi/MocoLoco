# **MOCOLOCO** tool instruction and report <br>
## A tool to analyze a TF chip-seq data and to find nearby TFBS to discover or validate co-activity between different TFs in transcription activation process
<br>  

## **MOCOLOCO INPUT/OUTPUT**
<br>

### The program has **supports multiple inputs**.
<ul><br>
1. BED, TwoBit and Jaspar files in input.<br>
<ul><li> A BED file format is a text file used to store genomic regions as coordinates (and eventually also associated annotations).
<li> A TwoBit is a binary file format used for storing big biological data (Genomes,chromosomes, etc..). It is useful because of his reading speed, which is an order of magnitude faster than a FASTA file format.
<li> A JASPAR matrix is a file format coming from JASPAR database, in which TF binding profiles are stored as position frequency matrices (PWMs).
<br><br></ul>
2. Multifasta file in input.
<ul><li> A Multifasta file format contains multiple FASTA sequences, which must be of the same length.
</ul><br>  

### To **compile** the program the command needed in both cases is:

      g++ -o test -Wall Mocoloco.cpp -lgsl (-lgslcblas command is need to be added in certain systems)
<br>  

### The tool can be **executed** in different ways:
<br><ul>
1.**BED,TwoBit,Jaspar input**:
       
       ./test -b <file.bed>  -p <number>(optional) -t <file.2bit> -j <file.JASPAR> -k <number, number, ...>(optional) -n <number> (optional) -ss(optional) -a(optional) -r(optional) -e(optional)`
       
2.**Multifasta input**:
       
       ./test -m <file_multifasta.fa>  -p <number>(optional) -k <number, number, ...>(optional) -ss(optional) -a(optional) -r(optional) -e(optional)
</ul><br>  

### For each k inserted as input the **output** is composed by different files in .txt format. The program can be run in Double Strand, to analyze both the DNA helices, or in Single Strand, to perform the analysis only on the forward strand.
	
   <br>
   
   ### A. Output file which contains the total k-mers occurrences in all the sequences, ordered from the most numerous to the less.<br> 
   <ul> 		
   <li> If the analysis is in DS, the output file contains oligo occurrences in forward strand and Reverse Complement occurrences in reverse strand (the two occurrences must be equal).

   <li> If the analysis is in SS the output file contains just oligo occurrences.
   </ul><br>

   ### <br> B. Output file which contains the rank of the best -n oligo for each position. 
   <br>

         #Position       Rank    Oligo   Num_Occ_FWD     Num_Occ_REV     Sum_Occ_Oligo   Oligo_RC        Num_Occ_RC_FWD  Num_Occ_RC_REV  Sum_Occ_RC      PAL     Tot_Occ FREQ    P_VALUE

   <ul>
   <li> For DS analysis the tool reports:
   <ul>    
   <li> The position.
   <li> The rank.
   <li> The oligo occurrences in FWD strand.
   <li> The reverse complement occurrences in reverse strand.
   <li> The sum of their occurrences.
   <li> The occurrences of reverse complement in forward strand.
   <li> The occurrences of the oligo in reverse strand.
   <li> The sum of their occurrences.
   <li> TRUE if the oligo is palindrome, FALSE if not.
   <li> The frequency of occurrences at that position.
   </ul>
       
   <br>
   <li> B2. For SS analysis:
   <ul>    
   <li> The position.
   <li> The rank.
   <li> The oligo occurrences in FWD strand.
   <li> The occurrences of reverse complement on FWD strand.
   <li> TRUE if the oligo is palindrome, FALSE if not.
   <li> The frequency of oligo occurrences at that position.
   </ul>
   </ul>

  ### <br> C. Output which contains the sum and frequences of the top -n oligos for each position.
<br>

### D. Output which contains, for each oligo, his p-value (and all the parameters used to calculate it):
<br>
<ul>	
<li> Position in sequence.
<li> Rank in the position.
<li> Oligo string.
<li> K parameter (Occurrences of the oligo in that position).
<li> N1 parameter (Occurrences of the oligo in every position).
<li> N2 parameter (Total possible oligos - N1).
<li> T parameter (Number of sequences).
<li> P_value.
</ul> <br>

### There are also **two supplementary output file**, one in .bed format and the other in .fasta format. <br>These output are created only if the program is executed with **bed/twobit/jaspar input**. They are:
<br>

### D. A .bed file with the genomic coordinates of the sequences analized by MocoLoco.cpp tool. The file contains, for each sequence, information about:
<br>
<ul>
   <li> The Chromosome
   <li> The Starting Coordinate
   <li> The Ending Coordinate
</ul>
<br>

### E. A .fasta file which contains the sequences analyzed by MocoLoco tool in fasta format.
   <br>
   <ul>
   For every sequence it contains:

   <li> Header (fasta format) with genomic coordinates.
   <li> The fasta sequence extracted from TwoBit file. The sequences are centered where the PWM JASPAR matrix fits better (probably on the TFBS) and a window (of -p x 2 length) is generated.
</ul>
<br>

### All the output file names contains a brief description of their content. Furthermore they are generated using the bed name, the twobit genome name, the JASPAR matrix name, (eventually the multifasta name), the k-mer value and finally DS/SS annotation. If multiple analysis are performed at the same time this method guarantees to always be able to trace the origin of the data.
<br>

## **MOCOLOCO WORKFLOW**
<br>
Firstly all libraries, global variables, classes and function declarations, which are written in MocoLoco.h file, are included.
Then the <i><u>command_line_parser</i></u> function handles the input parameter assigning them to the rignt variables. This function is also able to control the validity and the existance of parameters and it is trained to report errors in case of incorrect arguments.<br>
Once global variables are assigned correctly, the tool is ready to process the informations.<br>
The function <i><u>GEP_path</i></u> analyze the input and choose the right workflow.<br>
<br>

### **BED/TWOBIT/JASPAR Pathway**:

If MocoLoco recognizes a bed/twobit/jaspar input, it immediately creates a <i><u>coordinator class C</i></u>, which handles and connects the arguments.<br><br>
### **COORDINATOR CLASS** <br>
**This class is the connecting hub between the bed, twobit and jaspar files**. It allows the **creation of bed classes** by extracting genomic sequences from the twobit file. Then it guides the **creation and managing of a jaspar matrix class**. Finally the Coordinator class **connects them to each other thanks to the creation of an additional oligo class**.<br>
The first coordinator class function is <i><u>GEP creation</i></u>:
<br>This function has the task of **opening the input bed file**,** reading its contents** and finally, for each genomic coordinate readed in the bed file, **creating a bed class**. Each bed class created is then inserted into a **vector of bed classes called GEP** (GEnomic Position).<br><br>
<ul>

### **BED CLASS** <br>
For each genomic coordinate in the bed file an object of <i><u>bed class</i></u> is created.<br> 
Two different class constructors has been implemented to handle different inputs. <br>
In this case (Bed/twobit/jaspar input) **the constructor used is the one with the passage of 4 parameters** (**p** = half window legth, **line** = current bed line **tb** = Twobit file, **n_line** = current bed line number).<br>
Once the constructor is called and the parameters are passed, subsequently the data are managed through 4 bed class functions:


1. <i><u>read_line</u></i>:<br>
Takes the input line (from bed file) and **splits word by word**, **extracting and saving information about chromosome**, **starting coordinate and end coordinate**.<br>

2. <i><u>flag_control</i></u>: <br>
**Checks if end coordinate is greater than start coordinate** as it should be. If this is confirmed the flag is set to 1, otherwise the flag is set to 0.<br>

3. <i><u>centering_function</i></u>:<br>
**Centers the genomic coordinates** in the middle of start/end coordinated provided by Bed file. 
Subsequently, if for example -p is 150, a window of 300 bases is created and an overhead of 25(default value) is added to end coordinates.
<br>

4. <i><u>extract_seq</i></u>:<br>
Takes the Twobit file passed as input by user and, **analyzing class variables**, **extracts from the 2bit genome the sequence corresponding to the class coordinates**, saving it in a string.
If the flag has been set to 0 the function can not extract any sequence and it returns a warning message.
<br><br>
</ul>
</ul>

When <i><u>GEP_creation</i></u> function ends, a **GEP vector** has been created. This vector carries informations about the genomic coordinates but, more importantly, it **carries the sequences extracted from the genome**.
<br>
The next step of Coordinator class is to call <i><u>matrix_class</i></u> contstructor, which takes as input the Jaspar matrix provided by user.<br><br>
<ul>

### **MATRIX CLASS**<br>
The matrix class has the function of **reading the jaspar matrix** entered by the user, **storing it** and finally **modifying its values ​​with normalization and logarithmic operations**.
The class constructor is made by 6 functions able to do the operations described before.<br>
These functions are:
<br><br>
1. <i><u>read_JASPAR</i></u>:<br>
This function is the hub of the constructor. It **takes the Jaspar file inserted in input and extracts the numeric values**, saving them in a double type vector of vectors. Then the header informations about matrix_name and tf_name are stored in two string variables.<br>

2. <i><u>matrix_normalization_prseudoc</i></u>:<br>
Here the **matrix values are normalized** for the first time. They are normalized **by column**, thus a vector of column sums is firstly calculated thanks to the function <i><u>find_col_sum</i></u>. Each value of the matrix is normalized by dividing it by the sum of its column. Then **a pseudocount (0,01) is added** to avoid to have some 0 into the matrix.<br>

3. <i><u>matrix_normalization</i></u>:<br>
A **second normalization is done** by this function to trim the values after the psudocout adding.
Also this time the <i><u>find_col_sum</i></u> function is used to recalculate the column sums.

4. <i><u> matrix_logarithmic</i></u>:<br>
After the two normalization another filtering is required. With this function the tool takes matrix values and **substitutes them with their natural logarithms**.<br>

5. <i><u> reverse_matrix</i></u>:<br>
A simple function called reverse matrix is used to reverse the matrix values. The elements of the matrix are rotated along two axis. For example the first value (0,0) becomes the last (n,n) and viceversa. The aim of that is to create a jaspar matrix which is able to analyze the reverse strand in double strand workflow.<br>

All these functions are used to read and work on the jaspar matrix in input. All the matrices obtained during normalizing and logarithmic processes are stored in double vector of vector varibles.
<br><br>
</ul>
Returning into the <i><u>coordinator_class</i></u> the program has created a vector of <i><u>bed_class</i></u> containing the genome sequences extracted following the bed coordinates in input and also a <i><u>matrix_class</i></u> containing the jaspar matrix provided by user normalized and filtered as described before.<br>
The coordinator class puts in communication these two classed with the function <i><u>oligos_vector_creation</i></u>:<br>
This function has the aim to connect GEP vector to jaspar matrix and to do that it calls the costructor of the new <i><u>oligo_class</i></u>, passing the matrix, the sequence, the start/end coordinates and the strand. If the analysis is on SS just a call (for forward strand) is required, otherwise the oligo_class constructor is called two times for the same sequence:<br>
1. Passing matrix_log and strand + sign to perform the analysis on the forward strand.<br>
2. Passing the matrix_log_inverse and strand - sign to perform the analysis on the reverse strand.<br>
The oligo_class contructor is called for each sequence into the GEP vector.<br> 

**The goal of each class is to scroll the matrix along each sequence and define a score for each oligo**. The score will be saved in a **score vector** and each sequence will have its own score vector.<br> 
All oligo classes created, for which a score vector has been generated and stored, are in turn **saved in an oligo_class vector**.<br><br>
<ul>

### **OLIGO CLASS** 

As said before this class takes in input 5 parameters: Jaspar matrix, sequence, start/end coordinates, strand. This class is composed principally by 6 fuctions, which aim is to analyze sequence and jaspar matrix connections. <br>
These functions are:<br>

1. <i><u>find_minmax</i></u>:<br>
This function calculates the maximum and the minimum possible score that an oligo can have and it saves them into two double variables.

2. <i><u>shifting</i></u>:<br>
This is the hub function of the class. It allows to **shift the sequence on the matrix** and, for each oligo, to **calculate the score following matrix values**. Not only, this function **stores the oligo score in a oligo_vector** to keep the information saved and recalls itself recoursively until the end of the sequence.

3. <i><u>find_best_score</i></u>:<br>
**Here the best match, for each sequence, between oligo and jaspar matrix is found**.<br> If the best oligo is only one then, once identified, its **position** is extracted, returned by the function and finally **saved into local_position variable**.<br>
On the other hand, if the best score is found in more then one oligo, the function ensures that the position returned is the one closest to the centre.

4. <i><u>best_score_normalization</i></u>:<br>
This function simply **normalizes the best score extracted using the normalization formula**:

$$ BestscoreNorm =  1 + \frac{best score - maxpossiblescore}{maxpossiblescore - minpossiblescore} $$ 
<br>

5. <i><u>find_best_sequence</i></u>:<br>
This function **extracts from the sequence the substring** corresponding to the oligo which gave the best score.

6. <i><u>find_coordinate</i></u>:<br>
Finally with find_coordinate function the tool is able to **find the coordinates of the oligo which gave the best score**. The chromosome, the oligo's starting/ending coordinates are therefore saved in 3 variables.
<br>

The process described takes place for each sequence contained in the GEP vector (if the analysis is on DS twice per sequence) and, in the end, what will be obtained will be **an <i><u>oligo_class vector</i></u>, in which the informations about each oligo class are stored** (best score, best score position, coordinates, strand, etc..).<br><br></ul>
Then, returning into <i><u>coordinator_class</i></u> the tool still has to perform two tasks: **select the best strand** (if a DS analysis is being performed) and then **center the sequences** of the GEP vector on the oligo that gave the best score.<br>
The first task is performed by the <i><u>best_strand</i></u> function, **which can choose**, comparing the best score obtained by forward strand with that obtained by reverse strand, **which helix to keep and which one to discard**.<br>
Once the strand has been chosen the last step is made by <i><u>centering_oligo</i></u> function:<br> This function scrolls the oligo vector element by element and, for each one, **extracts the coordinate of the best oligo center**.<br>
This coordinate will be the starting point for **centering the whole sequence right on the oligo center**, thus **having, for each genomic sequence, the best oligo in a central position**. To do this, the tool reuses two functions already described: <i><u>the centering function</i></u> and the <i><u>extract sequence</i></u>.<br>
Now the GEP vector carries, for each bed class, the genomic sequence exactly centered in the middle of the best score oligo.<br><br>
### **MULTIFASTA PATHWAY**
If the input inserted is a multifasta file the pathway followed by MocoLoco is different. In fact, returning on <i><u>GEP_path</i></u> function, **if the tool recognize a multifasta input, it proceeds to create a <i><u>multifasta_class</i></u>**.<br>
The multifasta class is a short and simple class which aims to **extract from .fasta file the sequences and store them into a GEP vector**.<br>
This class is composed by 3 functions and its constructor takes in input just the multifasta file string. These functions are:

1. <i><u>extract_sequences</i></u>:<br>
It **allows to extract, from a multifasta file, all the ATCG sequences deleting headers and other useless stuff**. The sequences extracted are stored into a string vector.

2. <i><u>length_control</i></u>:<br>
This is a control function. Its goal is to **check if all the sequences extracted from multifasta file have the same length**. If this is not the case the tool prints an error and exits.

3. <i><u>GEP_creation_MF</i></u>:<br>
This function has the aim to **create a GEP vector from the sequences extracted**. For each sequence, a bed_class constructor is called.<br>
This constructor has some differences from the one called by the other pathway; for example it assign to chromosome coordinate "Multifasta" and start/end coord are set to 0.<br>The main important operation is **the storing of the sequence** and, consequentely, **the creation of a GEP vector made by all these bed classes created**.
<br><br>

Once the GEP vector has been created the workflow proceeds along only one branch. **The two different pathways are in fact used only for create a GEP vector from two different input files**, but now they are merged together into a single branch.<br><br>

### **MAP CLASS** <br>

To analyze, classify and print out the result the tool uses a new class called <i><u>map_class</i></u>. 
**This is the main class of the project and it allows**, taking the GEP vector and the number of k parameters inserted by user, **to make a complete and successfull sequence analysis**.


