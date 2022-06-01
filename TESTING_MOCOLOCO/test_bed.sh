#!/bin/bash
#	$RMC -n $N -l $L -j ../../${J} -p $P -o $i && $MOCO -m random_multifa_implanted1.fasta  -k $K -d $D &
usage() { echo "Usage: $0 -j <JASPAR_MATRIX> -f <Hit.txt> -n <n> -l <l> -p <p> -k <k> -d <d> -e <e> " 1>&2; exit 1; }

while getopts ":b:t:j:z:o:k:d:f:v:e:ra" o; do
    case "${o}" in
        
	b)
            B=${OPTARG}
            ;;
        t)
            T=${OPTARG}
            ;;
        j)
            J=${OPTARG}
            ;;
        z)
            Z=${OPTARG}
            ;;
        o)
            O=${OPTARG}
            ;;
        k)
            K=${OPTARG}
            ;;
        d)
            D=${OPTARG}
            ;;
        f)
            F=${OPTARG}
            ;;
        v)
            V=${OPTARG}
            ;;
	e)
	    E=${OPTARG}
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

if [ -z "${B}" ] || [ -z "${T}" ] || [ -z "${J}" ] || [ -z "${Z}" ] || [ -z "${K}" ] || [ -z "${D}" ];
then

    usage
fi

if [ -z "$E" ]
then
	E=0
fi

if [ -z "$V" ]
then
	V=1
fi

if [ -z "$F" ]
then
	F=0.02
fi

MOCO=$(realpath MOCO)
path_out=$(realpath $Z)
path_out_tot=${path_out::-4}_tot.txt
Out_dir=${B#../*/};

if [ -z "${Refine}" ]
then

	mkdir ${Out_dir}_test_k${K}_f${F}_em${E};

else	

	mkdir ${Out_dir}_test_k${K}_f${F}_r;

fi 

echo "#TESTING MOCOLOCO">$path_out_tot;
echo -e "HIT\tPVAL\tPVAL-LOG10\tBONF-PVAL\tLOG10BONF">>$path_out_tot; 

if [ -z "${Refine}" ]
then

	cd ${Out_dir}_test_k${K}_f${F}_em${E};

else	

	cd ${Out_dir}_test_k${K}_f${F}_r;

fi 

if [ -z "${O}" ]
then
	$MOCO -b ../${B} -t ../${T} -j ../${J} -k $K -d $D -f $F -e $E $Refine $all &
else
	$MOCO -b ../${B} -t ../${T} -j ../${J} -k $K -d $D -f $F -e $E $Refine $all -o $O &
fi

wait

	#Extraction of all the pvalue from the Z_scores_implanted files
	awk -v p_val=$V  '!/^#|^$/ { if($8<p_val) print $1"\t"$8"\t"$9"\t"$10"\t"$11 } '  *Z_scores_* >> $path_out_tot;

cd ..;	
#	a=$(awk  '/#Position/ {sub(/:/ ,"" ); sub(/#Position/,""); {ORS=" "}; print }' *mers_PWM_hamming_matrices* );
#        array=();
#        for j in $a 
#        do
#        pval=$(awk -v j="$j" '{if ($1 == j) print $9;}' *mers_p_value_parameters_control* | head -n 1) 
#	array+=($pval)
#        logpval=$(awk -v j="$j" '{if ($1 == j) print $9;}' *Z_scores* | head -n 1) 
#	array+=($logpval)
#        bestoligo=$(awk -v j="$j" '{if ($1 == j) print $2;}' *Z_scores* | head -n 1) 
#	array+=($bestoligo)
        
#	echo -e $j "\t" $bestoligo "\t" $pval "\t" $logpval>> $path_out; 
#        done  
cd ..;
