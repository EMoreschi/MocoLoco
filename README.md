# **MocoLoco**

## MocoLoco project

MocoLoco aims to identify positionally-constrained conserved motifs within the genomic regions surrounding TF binding sites identified through ChIP-Seq or similar assays.

If you want to know how the tool works: [MocoLoco_REPORT](https://github.com/EMoreschi/MocoLoco/blob/main/MocoLoco_REPORT.md).

# **How to run MocoLoco**

### **How to compile:**
`g++ -Wall -o MOCO  MocoLoco.cpp -lgsl -lgslcblas -O3` 
### **Gsl-2.6 library**

The current stable version of GSL is always available from ftp.gnu.org
in the directory /pub/gnu/gsl.

A list of mirror sites can be found at http://www.gnu.org/order/ftp.html

The project homepage is http://www.gnu.org/software/gsl/
### **How to run (2 possibilities):**

``` bash
./MOCO -b <file.bed>  -p <number>(optional) -t <file.2bit> -j <file.JASPAR> -k <number, number, ...>(optional) -d <number, number, ...>(optional) -o p(optional) -e <number>(optional) -f <number, number, ...>(optional) -s(optional) -u(optional) -l(optional) -z <number>(optional) -r <number>(optional)

./MOCO -m <file_multifasta.fa>  -p <number>(optional) -k <number, number, ...>(optional) -d <number, number, ...>(optional) -o p(optional) -e <number>(optional) -f <number, number, ...>(optional) -s(optional) -l(optional) -z <number>(optional) -r <number>(optional)
```

|Option |Parameter  |Description                                        |
|:-----:|:---------:|:-------------------------------------------------:|
| -b    | --bed     |Input bed file                                     |
| -m    | --mf      |Input multifasta file                                                      |
| -t    | --twobit  |It is used with bed file to know from which genome extract the sequences   |
| -j    | --jaspar  |Specify the jaspar file used as primary motif (used only if bed file as input file)        |
| -k    | --kmer    |Length of k-mers in the analysis (DEFAULT: 6,8,10)                         |
| -d    | --distance|It refers to the Hamming distance for the creation of oligos cluster (DEFAULT: 1,2,3)   |
| -p    | --param   |It is the half-length of the sequences (DEFAULT: 150)                      |
| -s    | --ss      |To run the anlysis as single strand (DEFAULT: double strand)               |
| -o    | --ordering|It defines if the best oligos are ordered by pvalue (```-o p```) or by occurrences (DEFAULT)|
| -f    | --freq    |Set the frequence treshold to calculate the Zscore                         |
| -e    | --exp_max |Refine the PWM matrices using the expectation maximization algorithm       |
| -l    | --tomtom  |Set the output format suitable for Tomtom from MEME                        |
| -u    | --unidirection|It orders the sequences in the same direction                          |
| -r    | --secondary|Set how secondary matrices I need in the analysis                         |
| -z    | --z_pval_threshold|Set the Zscore p-value threshold for eventually secondary matrices |
| -a    | --cleaning|Disable the cleaning of sequences with low primary motif scores (DEFAULT: ON)|
# **Dependencies**

### **TwoBit**
This script uses the twobit.c and twobit.h file from [this repository](https://github.com/andrelmartins/TwoBit), you can find these files in TwoBit folder. 

# **Repository folders:**

In this repository there are some useful folders for the usage of MocoLoco:
* **Genomes** where there are some TwoBit files of most important genomes
* **Jaspar_2020** contains most of the matrices present in [Jaspar database](https://jaspar.genereg.net)
* **Test_Bed** is a folder with some bed files from ENCODE that can be used as test files
* **Random_multifa_TOOL** contains an important tool used during MocoLoco testing, this tool is explained in the next chapter.
* **TESTING_MOCOLOCO** is a folder with some bash scripts used during MocoLoco testing.

# **Testing tools**
## **Multifa random tool**

It is a tool able to generate multifasta files with user-defined TFBSs motifs implanted at specific positions.

#### How to compile

`g++ -o RMC Multifa_random_tool.cpp`

#### How to run (2 possibilities)

##### To generate random multifasta files and make implants:

``` bash
./RMC -j <JASPAR_FILE_1> <JASPAR_FILE_2> ... <JASPAR_FILE_N> -p <number, number, ....> -l <number> -n <number> -o <number, number,...> -w <number, number,...> -c <number> -f <number, number,..> 
```
##### To generate random multifasta files:
``` bash
./RMC -l <number> -n <number> -c <number> 
```
|Option |Parameter  |Description                                        |
|:-----:|:---------:|:-------------------------------------------------:|
| -j    | --jaspar  |Specify the jaspar file used as primary motif for implants |
| -l    | --length  |Length of Multifasta sequences (DEFAULT: 500)              |
| -n    | --nseq    |Number of Multifasta sequences (DEFAULT: 200)              |
| -o    | --oligop  |Percentage of sequences in which the oligos will be implanted (DEFAULT: 0%)|
| -p    | --position|Position of the implants in the sequences                  |
| -w    | --wobble  |The implanting position, for every oligo, will be randomly choosen between *p-w* and *p+w* interval (DEFAULT: 0)|
| -c    | --cycles  |How many random multifasta files this tool will produce (DEFAULT: 1)|
| -f    | --freq    |Percentage of fwd/rev strand to select for the implants. (DEFAULT: 50)|
