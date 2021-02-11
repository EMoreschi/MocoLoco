#!/bin/bash

g++ -o test MocoLoco.cpp;

Jaspar="JASPAR_MA/MA0002.2.jaspar JASPAR_MA/MA0042.2.jaspar JASPAR_MA/MA0071.1.jaspar JASPAR_MA/MA0093.3.jaspar"
Twobit="hg38.2bit hg19.2bit"
BED="baitedreagionsb38.bed nfy_hg19.bed"

for a in $BED 
do
	for b in $Twobit 
	do
		for c in $Jaspar 
		do
			g++ -o test MocoLoco.cpp;
			./test --BED ${a} -tb ${b} -J ${c} >> $1;
			echo "----------------------------------------------------------------------------------------------------------------------------" >> $1;
			
		done
	done
done
