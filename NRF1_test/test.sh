#!/bin/bash

touch $2;
RMC=$(realpath RMC)
MOCO=$(realpath MOCO)
path_out=$(realpath $2)

mkdir $1_Out;
echo "TESTING MOCOLOCOC">$path_out;
frequenze=(100 80 60 40 30 20 10 9 8 7 6 5 4) 

for i in ${frequenze[@]}
do 
  cd $1_Out;
  mkdir $i;
  cd $i;
	$RMC -n $3 -l $4 -j ../../$1 -p $5 -o $i && $MOCO -m random_multifa_implanted1.fasta  -k $6 -d $7 &
  cd ../..;
done

wait

 for x in ${frequenze[@]}
do
 cd $1_Out/$x;
        echo "# Freq:" $x >> $path_out;        
	a=$(awk  '/#Position/ {sub(/:/ ,"" ); sub(/#Position/,""); {ORS=" "}; print }' *mers_PWM_hamming_matrices_multifasta_DS.txt );
        echo "Hit : "$a >> $path_out;
        array=();
       for j in $a 
        do
        k=$(awk -v j="$j" '{if ($1 == j) print $9;}' *mers_p_value_parameters_control_occ_DS.txt | head -n 1)
       array+=($k)
        done  
       echo "Pval-log10: "${array[@]} >> $path_out; 
cd ../..;
done
