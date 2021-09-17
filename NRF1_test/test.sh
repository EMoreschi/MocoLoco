#!/bin/bash
#	$RMC -n $N -l $L -j ../../${J} -p $P -o $i && $MOCO -m random_multifa_implanted1.fasta  -k $K -d $D &
usage() { echo "Usage: $0 -j <JASPAR_MATRIX> -f <Hit.txt> -n <n> -l <l> -p <p> -k <k> -d <d> " 1>&2; exit 1; }

while getopts ":j:f:n:l:p:k:d:" o; do
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
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${J}" ] || [ -z "${F}" ] || [ -z "${N}" ] || [ -z "${L}" ] || [ -z "${P}" ] || [ -z "${K}" ] || [ -z "${D}" ] || [ -z "${F}" ]; then
    usage
fi

touch $F;
RMC=$(realpath RMC)
MOCO=$(realpath MOCO)
path_out=$(realpath $F)

mkdir ${J}_Out;
echo "TESTING MOCOLOCOC">$path_out;
frequenze=(100 80 60 40 30 20 10 9 8 7 6 5 4) 

for i in ${frequenze[@]}
do 
  cd ${J}_Out;
  mkdir $i;
  cd $i;
	$RMC -n $N -l $L -j ../../${J} -p $P -o $i && $MOCO -m random_multifa_implanted1.fasta  -k $K -d $D &
  cd ../..;
done

wait

 for x in ${frequenze[@]}
do
 cd ${J}_Out/$x;
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
