# **MocoLoco**

## MocoLoco project

MocoLoco aims to identify positionally-constrained conserved motifs within the genomic regions surrounding TF binding sites identified through ChIP-Seq or similar assays.

If you want to know how the tool works: [MocoLoco_REPORT](https://github.com/EMoreschi/MocoLoco/blob/main/MocoLoco_REPORT.md).

# **How to run MocoLoco**

### **How to compile:**
`g++ -Wall -o MOCO  MocoLoco.cpp -lgsl -lgslcblas -O3` 

### **How to run (2 possibilities):**

``` bash
./MOCO -b <file.bed>  -p <number>(optional) -t <file.2bit> -j <file.JASPAR> -k <number, number, ...>(optional) -d <number, number, ...>(optional) -o p(optional) -e <number>(optional) -f <number, number, ...>(optional) -s(optional) -u(optional) -l(optional) -z <number>(optional) -r <number>(optional)

./MOCO -m <file_multifasta.fa>  -p <number>(optional) -k <number, number, ...>(optional) -d <number, number, ...>(optional) -o p(optional) -e <number>(optional) -f <number, number, ...>(optional) -s(optional) -l(optional) -z <number>(optional) -r <number>(optional)
```
# **Dependencies**

### **TwoBit**
This script uses the twobit.c and twobit.h file from [this repository](https://github.com/andrelmartins/TwoBit), you can find these files in TwoBit folder. 

### **Gsl-2.6 library**

The current stable version of GSL is always available from ftp.gnu.org
in the directory /pub/gnu/gsl.

A list of mirror sites can be found at http://www.gnu.org/order/ftp.html

The project homepage is http://www.gnu.org/software/gsl/

# **Repository folders:**

In this repository there are some useful folders for the usage of MocoLoco:
* **Genomes** where there are some TwoBit files of most important genomes
* **Jaspar_2020** contains most of the matrices present in [Jaspar database](https://jaspar.genereg.net)
* **Test_Bed** is a folder with some bed files from ENCODE that can be used as test files
* **Random_multifa_TOOL** contains an important tool used during MocoLoco testing, this tool is explained in the next chapter.
* **TESTING_MOCOLOCO** is a folder with some bash scripts used during MocoLoco testing, as for Multifa random tool also these scripts are explained in the next chapter.

# **Testing tools**
## **Multifa random tool**

The tool creates a set of random multifasta files.

#### How to compile

`g++ -Wall -o RMC Multifa_random_tool.cpp`

#### How to run (2 possibilities)

##### To generate random multifasta files and make implants:

``` bash
./RMC -j <JASPAR_FILE_1> <JASPAR_FILE_2> ... <JASPAR_FILE_N> -p <number, number, ....> -l <number> (default 500) -n <number> (default 200) -o <number, number,...> -w <number, number,...> -c <number> -f <number, number,..> 
```
##### To generate random multifasta files:
``` bash
./RMC -l <number> (default 500) -n <number> (default 200) -c <number> 
```

## **Test multifasta**

## **Test Bed**