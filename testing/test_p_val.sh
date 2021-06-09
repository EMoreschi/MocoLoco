#!/bin/bash

Matrix=$1;
Position=$2;
N_seq=$3;
Seq_length=$4;
Wobble=0;
Frequence_FWD=50;
K=$5;
Output=$6;
Filename="-mers_positional_occurrences__DS.txt";

Line_char_number=$(cat $Matrix | sed -n '2p' | gawk '{ print NF}');
Half_line=$(($(($Line_char_number - 3)) / 2));
Initial_position=$(($Position - $Half_line));

#if matrix length is even, then $Line_char_number%2 != 0 -> made to ensure it extracts from the starting implant position
if (( $Line_char_number%2 != 0 )) 
then
	Initial_position=$(($Initial_position+1));
fi

i=1;

echo "IMPLANTING_FREQ	POSITION	RANK	OLIGO	OLIGO_FWD_OCC	OLIGO_REV_OCC	SUM_OLIGO_OCC	OLIGO_RC	OLIGO_RC_FWD_OCC	OLIGO_RC_REV_OCC	SUM_OLIGO_RC_OCC	PALINDROME	TOT_SUM	FREQ	P_VALUE" >> $Output;

while [ $i -ne 101 ]
do
	
	g++ -o test -Wall ../Random_multifa_TOOL/Multifa_random_tool.cpp;
	./test -n $N_seq -l $Seq_length -j $Matrix -p $Position -o $i -w $Wobble -f $Frequence_FWD;

	g++ -o test -Wall ../MocoLoco.cpp -lgsl -lgslcblas &&
	./test -m random_multifa_implanted1.fasta -k $K;
	
	touch $Output;
	
	awk -v pos=$Initial_position -v perc=$i '$1 == pos && $2 == 1 {print perc "%	" $line}' $K$Filename >> $Output;
	
	i=$(($i+1));
done

rm *mer*;
rm random_multifa_1.fasta;
rm random_multifa_implanted1.fasta;




















