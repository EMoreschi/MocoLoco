# **MocoLoco**

## MocoLoco project

#### How to compile
`g++ -Wall -o MOCO  MocoLoco.cpp -lgsl -lgslcblas` 
#### How to run (2 possibilities)

``` bash
./MOCO -b <file.bed>  -p <number>(optional) -t <file.2bit> -j <file.JASPAR> -k <number, number, ...>(optional) -n <number> (optional) -ss(optional)

./MOCO -m <file_multifasta.fa>  -p <number>(optional) -k <number, number, ...>(optional) -ss(optional)
```
## Multifa random tool

The tool is contained in the Random_multifa_TOOL directory

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

### **This script uses the twobit.c and twobit.h file from**
https://github.com/andrelmartins/TwoBit 

### **This script uses the gsl-2.6 library**

The current stable version of GSL is always available from ftp.gnu.org
in the directory /pub/gnu/gsl.

A list of mirror sites can be found at http://www.gnu.org/order/ftp.html

The project homepage is http://www.gnu.org/software/gsl/

