#!/bin/bash
#	$RMC -n $N -l $L -j ../../${J} -p $P -o $i && $MOCO -m random_multifa_implanted1.fasta  -k $K -d $D &
usage() { echo "Usage: $0 -j <JASPAR_MATRIX> -f <Hit.txt> -n <n> -l <l> -p <p> -k <k> -d <d> " 1>&2; exit 1; }

while getopts ":j:f:n:l:p:k:t:d:ra" o; do
    case "${o}" in
        j)
            J=${OPTARG}
            ;;
        f)
            F=${OPTARG}
            ;;
        n)
            N=${OPTARG}
            ;;
        t)
            T=${OPTARG}
            ;;
        l)
            L=${OPTARG}
            ;;
        p)
            P=${OPTARG}
            ;;
        k)
            K=${OPTARG}
            ;;
        d)
            D=${OPTARG}
            ;;
        r) 
           Refine="-r"
           ;;
        a) 
           all="-a"
           ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${J}" ] || [ -z "${F}" ] || [ -z "${N}" ] || [ -z "${T}" ] || [ -z "${L}" ] || [ -z "${P}" ] || [ -z "${K}" ] || [ -z "${D}" ] || [ -z "${F}" ]; then
    usage
fi

touch $F;
RMC=$(realpath RMC)
MOCO=$(realpath MOCO)
path_out=$(realpath $F)

Out_dir=${J#../*/};
mkdir ${Out_dir}_test;
echo "#TESTING MOCOLOCO">$path_out;
echo -e "FREQ \t HIT \t PVAL \t PVAL-LOG10">>$path_out; 
frequenze=(100 50 40 30 20 15 10 9 8 7 6 5) 

for i in ${frequenze[@]}
do 
  cd ${Out_dir}_test;
  mkdir $i;
  cd $i;
	$RMC -n $N -l $L -j ../../${J} -p $P -o $i && $MOCO -m random_multifa_implanted1.fasta  $Refine -k $K -d $D -f $T $all &
  cd ../..;
done

wait

 for x in ${frequenze[@]}
do
 cd ${Out_dir}_test/$x;
	a=$(awk  '/#Position/ {sub(/:/ ,"" ); sub(/#Position/,""); {ORS=" "}; print }' *mers_PWM_hamming_matrices_multifasta_DS.txt );
        array=();
       for j in $a 
        do
        pval=$(awk -v j="$j" '{if ($1 == j) print $9;}' *mers_p_value_parameters_control_occ_DS.txt | head -n 1)
        array+=($pval)
        logpval=$(awk -v j="$j" '{if ($1 == j) print $9;}' *Z_scores* | head -n 1) 
	array+=($logpval)
        echo -e $x "\t" $j "\t" $pval "\t" $logpval >> $path_out; 
        done  
cd ../..;
done
