#!/bin/bash
#$RMC -n $N -l $L -j ../../${J} -p $P -o $i && $MOCO -m random_multifa_implanted1.fasta  -k $K -d $D &

#-------USAGE----------------------------------------------------------------------------------------------

usage() { echo "Usage: $0 -j <JASPAR_MATRIX> -f <Output_file> -n <Sequences_number> -l <Sequences_length> -p <Implanting_position> -k <kmers_analyzed> -d <Hamming_distance> -c <Cycle_number> " 1>&2; exit 1; }

#--------PARSER----------------------------------------------------------------------------------------------

while getopts ":j:f:n:l:c:p:k:t:d:ra" o; do
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
        c)
            C=${OPTARG}
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

if [ -z "${J}" ] || [ -z "${F}" ] || [ -z "${N}" ] || [ -z "${T}" ] || [ -z "${L}" ] || [ -z "${P}" ] || [ -z "${K}" ] || [ -z "${D}" ] || [ -z "${C}" ] || [ -z "${F}" ]; then
    usage
fi

#----------DIRECTORIES PREPARATION-------------------------------------------------------------------------
                                                                                                             
touch $F;
RMC=$(realpath RMC)
MOCO=$(realpath MOCO)
path_out=$(realpath $F)

#Defining output directory name (<Jaspar_name>_test_k<k>)
Out_dir=${J#../*/};
mkdir ${Out_dir}_test_k${K};

#Initializing headers in Output file (Outside directory just created before)
echo "#TESTING MOCOLOCO">$path_out;
echo -e "FREQ \t HIT \t OLIGO \t PVAL \t PVAL-LOG10">>$path_out; 

#Defining frequences for analysis
frequenze=(75 65 55 45 35 25 20 15 10 5);

cd ${Out_dir}_test_k${K};

for freq in ${frequenze[@]}
do 
	mkdir $freq;
done

#--------RANDOM & RANDOM IMPLANTED FILE CREATION--------------------------------------------------------------
	
for freq in ${frequenze[@]}
do 
	cd $freq;
		
		$RMC -n $N -l $L -j ../../${J} -p $P -o $freq -c $C &
  
	cd ..;
done
wait

#-------RUNNING MOCOLOCO ON IMPLANTED MULTIFASTA-------------------------------------------------------------

nc=10
(
for freq in ${frequenze[@]}
do 
	cd $freq;

	for j in $(seq 1 $C);
	do
        	((i=i%nc)); ((i++==0)) && wait
		echo "random_multifa_implanted${j}"
		$MOCO -m random_multifa_implanted${j}.fasta  $Refine -k $K -d $D -f $T $all &

		
	done
	cd ..;
done
)

