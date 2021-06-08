#!/bin/bash

Matrix=$1;
Position=$2;
N_seq=$3;
Seq_length=300;
Wobble=0;
Frequence_FWD=50;
K=8;
Output=$4;

Line=$(cat $Matrix | sed -n '2p' | gawk '{ print NF}');
Half_line=$(($(($Line - 3)) / 2));
Initial_position=$(($Position - $Half_line));

echo $Line;
echo $Half_line;
echo $Initial_position;

i=1

while [ $i -ne 100 ]
do
	
	g++ -o test -Wall ./Random_multifa_TOOL/Multifa_random_tool.cpp;
	./test -n $N_seq -l $Seq_length -j $Matrix -p $Position -o $i -w $Wobble -f $Frequence_FWD;

	g++ -o test -Wall ./MocoLoco.cpp -lgsl -lgslcblas &&
	./test -m random_multifa_implanted1.fasta -k $K;
	
	touch $Output;
	
	awk -v pos=$Initial_position '$1 == pos && $2 == 1 {print $line}' 8-mers_positional_occurrences__DS.txt >> $Output;
	
	i=$(($i+1));
done






















