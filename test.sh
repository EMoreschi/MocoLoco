#!/bin/bash

g++ -o test MocoLoco.cpp;

touch $1
rm $1

Jaspar="Jaspar_2020/MA0502.2.jaspar Jaspar_2020/MA0002.2.jaspar Jaspar_2020/MA0042.2.jaspar Jaspar_2020/MA0071.1.jaspar Jaspar_2020/MA0093.3.jaspar"
Twobit="Genomes/hg38/hg38.2bit Genomes/hg19/hg19.2bit"
BED="Test_Bed/baitedregionsb38.bed Test_Bed/nfy_k562_hg19.bed"

for a in $BED 
do
	for b in $Twobit 
	do
		for c in $Jaspar 
		do
			g++ -o test MocoLoco.cpp;
			./test --BED ${a} -tb ${b} -J ${c} -DS >> $1;
			echo "----------------------------------------------------------------------------------------------------------------------------" >> $1;
			
		done
	done
done
