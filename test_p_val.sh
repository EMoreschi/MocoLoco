#!/bin/bash

Matrix=$1;
Position=$2;
N_seq=$3;
Seq_length=300;
Wobble=0;
Frequence_FWD=50;
K=$4;
Output=$5;
Filename="-mers_positional_occurrences__DS.txt";

Line=$(cat $Matrix | sed -n '2p' | gawk '{ print NF}');
Half_line=$(($(($Line - 3)) / 2));
Initial_position=$(($Position - $Half_line));

#if matrix length is even, then $Line%2 != 0 
if (( $Line%2 != 0 )) 
then
	Initial_position=$(($Initial_position+1));
fi

i=1;

while [ $i -ne 100 ]
do
	
	g++ -o test -Wall ./Random_multifa_TOOL/Multifa_random_tool.cpp;
	./test -n $N_seq -l $Seq_length -j $Matrix -p $Position -o $i -w $Wobble -f $Frequence_FWD;

	g++ -o test -Wall ./MocoLoco.cpp -lgsl -lgslcblas &&
	./test -m random_multifa_implanted1.fasta -k $K;
	
	touch $Output;
	
	awk -v pos=$Initial_position '$1 == pos && $2 == 1 {print $line}' $K$Filename >> $Output;
	
	i=$(($i+1));
done






















