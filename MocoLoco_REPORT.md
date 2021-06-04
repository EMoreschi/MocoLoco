# **INSTRUCTION** and **REPORT** MocoLoco.cpp and MocoLoco.h files project
<br>  

## **MOCOLOCO INPUT/OUTPUT**
<br>

### The program has **two different ways to be executed**.
<ul><br>
1. BED, TwoBit and Jaspar files in input.<br>
2. Multifasta file in input.
</ul>
<br>  

### To **compile** the program the command needed in both cases is:

      g++ -o test -Wall Mocoloco.cpp -lgsl -lgslclbas
<br>  

### To **execute** it 2 different ways have been implemented:
<br><ul>
1.**BED,TwoBit,Jaspar input**:
       
       ./test -b <file.bed>  -p <number>(optional) -t <file.2bit> -j <file.JASPAR> -k <number, number, ...>(optional) -n <number> (optional) -ss(optional)`
       
2.**Multifasta input**:
       
       ./test -m <file_multifasta.fa>  -p <number>(optional) -k <number, number, ...>(optional) -ss(optional)
</ul><br>  

### The **output** is composed by **3 different files .txt** for each k inserted as input. The analysis can be done in Double Strand or in Single Strand.
	
   <br>
   
   ### A. The first one is the k-mers count in all the sequences extractred from bed coordinates inserted, ordered from the most relevant to the less.

   <br> 
   <ul> 		
   <li> If the analysis is in DS the output file contains oligo sequence + its occurrences in forward (FWD) strand and Reverse Complement (RC) sequence + its occurrences in reverse strand (the two occurrences must be equal). <br> The name of the file is <i> k-mers_occurrences_twobit_jaspar_bed_DS.txt.</i>

   <li> If the analysis is in SS the output file contains oligo sequence + its occurrences.
      <br>The name of the file is <i>k-mers_occurrences_twobit_jaspar_bed_SS.txt.</i>
   </ul>

   ### <br> B. The second one is a file with the k-mers count for each position of sequences. The file contains information, for each position, of the top -n kmers ranked by occurrences. 
   <br>
   <ul>
   <li> For DS analysis these information are:
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
   
   <br> The filename is <i>k-mers_positional_occurrences_twobit_jaspar_bed_DS.txt</i>
   </ul>
       
   <li> B2. For SS analysis these information are:
   <ul>    
   <li> The position.
   <li> The rank.
   <li> The oligo occurrences in FWD strand.
   <li> The occurrences of reverse complement on FWD strand.
   <li> TRUE if the oligo is palindrome, FALSE if not.
   <li> The frequency of oligo occurrences at that position.
   <br> The filename is <i>k-mers_positional_occurrences_twobit_jaspar_bed_SS.txt</i>
   </ul>
   </ul>

  ### <br> C. The third one is a file with the sum of the occurrences of the first top -n oligos for each position and the frequences of these sums for each position.
  <ul>
  <br> The output structure for DS and SS analysis is the same but obviously the data has come from different processes.
  <br>The name of the file is <i>k-mers_Topn_sum_and_frequence_DS.txt (For DS)</i>.
  <br> The name of the file is <i>k-mers_Topn_sum_and_frequence_SS.txt (For SS)</i>.
  </ul>
<br>

### There are also **two supplementary output file**, one in .bed format and the other in .fasta format. <br>These output are created only if the program is executed with **bed/twobit/jaspar input**. They are:
<br>

### D. A .bed file with the genomic coordinates of the sequences analized by MocoLoco.cpp tool. The file contains, for each sequence, information about:
<br>
<ul>
   <li> The Chromosome
   <li> The Starting Coordinate
   <li> The Ending Coordinate
	The filename is <i>twobit_jaspar_bed.bed.</i>
</ul>
<br>

### E. A .fasta file which contains the sequences analyzed by MocoLoco tool in fasta format.
   <br>
   <ul>
   For every sequence it contains:

   <li> Header (fasta format) with genomic coordinates.
   <li> The ATCG fasta sequence extracted from TwoBit file, centered following the coordinates bed in input and as long as indicated by the parameter p (length = p x 2).
	The filename is <i>twobit_jaspar_bed.fasta.</i>
</ul>
<br>

### Finally, there is **another supplementary output** file called <i>k-mers_p_value_control_parameters.txt</i> (one for each k inserted as input), which contains **all the parameter values used to calculate Oligos P-values**. <br>This file contains information about:
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

## **MOCOLOCO WORKFLOW**
<br>
Firstly the file MocoLoco.h is included with all the libraries, structures, global variables and functions in it.
The function <i><u>Main</i></u> read the arguments in input and, thanks to the function <i><u>command_line_parser</i></u>, it is able to assign input parameters to the right global variables. <br>
The <i><u>command_line_parser</i></u> function is also able to control the validity and the existance of parameters and it is trained to report errors in case of incorrect arguments, followed by the printing of <i><u>display_help</i></u> to assist users.<br>
Once global variables are assigned correctly, the tool is ready to process the informations.
The first function is <i><u>GEP_path</i></u>, which is able to distinguish if input is in bed/twobit/jaspar or multifasta format and, cosequently, choose the right path.<br>
<br>

### **BED/TWOBIT/JASPAR Pathway**:

If MocoLoco recognizes a bed/twobit/jaspar input, it immediately creates a <i><u>coordinator class C</i></u>, which handles and connects the arguments.<br><br>
### **COORDINATOR CLASS** <br>
**This class is the connecting hub between the bed, twobit and jaspar files**. It allows the **creation of bed classes** by extracting genomic sequences from the twobit file, the **creation and managing of a jaspar matrix class** and finally to **connect them to each other thanks to the creation of an additional oligo class**.<br>
The first coordinator class function is <i><u>GEP creation</i></u>:
<br>This function has the task of **opening the input bed file**,** reading its contents** and finally, for each genomic coordinate readed in the bed file, **creating a bed class**. Each bed class created is then inserted into a **vector of bed classes called GEP** (GEnomic Position).<br><br>
<ul>

### **BED CLASS** <br>
in object of <i><u>bed class</i></u> is created for each genomic coordinate in the bed file.<br> 
The class constructor provides that the object can be built in two different ways: one for a bed/twobit/jaspar input and one for a multifasta input. Now the first case is analyzed, so **the constructor used is the one with the passage of 4 parameters** (**p** = half window legth, **line** = current bed line **tb** = Twobit file, **n_line** = current bed line number).<br>
Once the constructor is called and the parameters are passed, the bed class is able to handle the data calling 4 functions:


1. <i><u>read_line</u></i>:<br>
This function has the task of taking the input line (from bed file) and **splits word by word**, **extracting and saving information about chromosome**, **starting coordinate and end coordinate**.<br>

2. <i><u>flag_control</i></u>: <br>
**This is a controlling function**. Its job is to **check if end coordinate is greater than start coordinate** as it should be. If this is confirmed the flag is set to 1, otherwise the flag is set to 0.<br>

3. <i><u>centering_function</i></u>:<br>
It has the task of **centering the genomic coordinates** in the middle of start/end window. The function reassigns the start and end coordinates following the half length parameter inserted as input by user and adds an overhead (default 25) to end coordinate.
The genomic window is centered on the basis of bed file informations and its length corresponds to that desired by user.<br>

4. <i><u>extract_seq</i></u>:<br>
This function takes the Twobit file passed as input by user and, **analyzing class variables**, **extracts from the 2bit genome the sequence corresponding to the class coordinates**; then the sequence is saved in a string.
If the flag has been set to 0 the function can not extract any sequence and it returns a warning message.
<br><br>
</ul>
</ul>

When <i><u>GEP_creation</i></u> function ends a **GEP vector** has been created. This vector carries informations about the genomic coordinates but more importantly it **carries the sequences extracted from the genome**, following those coordinates.
<br>
The next Coordinator class contructor step is the calling of <i><u>matrix_class</i></u> contstructor, which takes as input the Jaspar matrix provided by user.<br><br>
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

