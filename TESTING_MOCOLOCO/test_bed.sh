#!/bin/bash
#	$RMC -n $N -l $L -j ../../${J} -p $P -o $i && $MOCO -m random_multifa_implanted1.fasta  -k $K -d $D &
usage() { echo "Usage: $0 -j <JASPAR_MATRIX> -f <Hit.txt> -n <n> -l <l> -p <p> -k <k> -d <d> " 1>&2; exit 1; }

while getopts ":j:f:b:t:p:k:s:d:ra" o; do
    case "${o}" in
        j)
            J=${OPTARG}
            ;;
        f)
            F=${OPTARG}
            ;;
        b)
            B=${OPTARG}
            ;;
        t)
            T=${OPTARG}
            ;;
        s)
            S=${OPTARG}
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

if [ -z "${J}" ] || [ -z "${F}" ] || [ -z "${B}" ] || [ -z "${S}" ] || [ -z "${T}" ] || [ -z "${K}" ] || [ -z "${D}" ] || [ -z "${F}" ]; then
    usage
fi

touch $F;
###RMC=$(realpath RMC)
MOCO=$(realpath MOCO)
path_out=$(realpath $F)

##mkdir ${J}_Out;
echo "#TESTING MOCOLOCO">$path_out;
echo -e "HIT \t OLIGO \t PVAL \t PVAL-LOG10">>$path_out; 
##frequenze=(100 50 40 30 20 15 10 9 8 7 6 5) 

##for i in ${frequenze[@]}
##do 
 ## cd ${J}_Out;
 Out_dir=${B#../*/};
  mkdir ${Out_dir}_test_k${K};
  cd ${Out_dir}_test_k${K};
	$MOCO -b ../${B} -t ../${T} -j ../${J} -k $K -d $D $Refine $all &
	###$RMC -n $N -l $L -j ../../${J} -p $P -o $i && $MOCO -m random_multifa_implanted1.fasta  $Refine -k $K -d $D -f $S $all &
  cd ..;
##done

wait

 ##for x in ${frequenze[@]}
 ##do
 cd ${Out_dir}_test_k${K};
	a=$(awk  '/#Position/ {sub(/:/ ,"" ); sub(/#Position/,""); {ORS=" "}; print }' *mers_PWM_hamming_matrices* );
        array=();
        for j in $a 
        do
        pval=$(awk -v j="$j" '{if ($1 == j) print $9;}' *mers_p_value_parameters_control* | head -n 1) 
	array+=($pval)
        logpval=$(awk -v j="$j" '{if ($1 == j) print $9;}' *Z_scores* | head -n 1) 
	array+=($logpval)
        bestoligo=$(awk -v j="$j" '{if ($1 == j) print $2;}' *Z_scores* | head -n 1) 
	array+=($bestoligo)
        
	echo -e $j "\t" $bestoligo "\t" $pval "\t" $logpval>> $path_out; 
        done  
cd ..;
##done
