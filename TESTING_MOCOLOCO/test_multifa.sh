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

if [ -z "${J}" ] || [ -z "${F}" ] || [ -z "${N}" ] || [ -z "${T}" ] || [ -z "${L}" ] || [ -z "${P}" ] || [ -z "${K}" ] || [ -z "${D}" ] || [ -z "${C}" ] || [ -z "${F}" ] ; then
    usage
fi

#----------DIRECTORIES PREPARATION-------------------------------------------------------------------------
                                                                                                             
touch $F;
RMC=$(realpath RMC)
MOCO=$(realpath MOCO)
#Creation of the path for the output file with the best pvalue for each cycle in each frequence
path_out=$(realpath $F).txt
#Creation of the path for the output file with all the pvalue obtained by our script
path_out_tot=$(realpath $F)_tot.txt
#Defining output directory name (<Jaspar_name>_test_k<k>)
Out_dir=${J#../*/};

if [ -z "${Refine}" ]
then

	mkdir ${Out_dir}_test_k${K}_c${C};

else	

	mkdir ${Out_dir}_test_k${K}_c${C}_r;

fi 

#Initializing headers in Output file (Outside directory just created before) and in the Output file with best pvalue
echo "#TESTING MOCOLOCO">$path_out;
echo -e "#FREQ\tHIT\tPVAL-LOG10">>$path_out; 
echo "#TESTING MOCOLOCO">$path_out_tot;
echo -e "#FREQ\tHIT\tPVAL-LOG10">>$path_out_tot; 

#Defining frequences for analysis
frequenze=(75 65 55 45 35 25 20 15 10 5);

if [ -z "${Refine}" ]
then

	cd ${Out_dir}_test_k${K}_c${C};

else	

	cd ${Out_dir}_test_k${K}_c${C}_r;

fi 

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

multi_thread=10
(
for freq in ${frequenze[@]}
do 
	cd $freq;

	for j in $(seq 1 $C);
	do
        	((i=i%multi_thread)); ((i++==0)) && wait
		echo "random_multifa_implanted${j}"
		$MOCO -m random_multifa_implanted${j}.fasta  $Refine -k $K -d $D -f $T $all &
		
		
	done
	wait
	cd ..;
done
)
wait
#-------EXTRACTING P-VALUES---------------------------------------------------------------------------------

#With this construct I obtain just the higher pvalue for each implanted file
for freq in ${frequenze[@]}
do
	cd $freq;

	#Extraction of all the pvalue from the Z_scores_implanted files
	awk -v fre=$freq  '!/^#|^$/ { print fre "\t"$1"\t"$9}' *Z_scores_* >> $path_out_tot;
	#Extraction of the best pvalue for each cycle in each frequence
	for f in *Z_scores* ; do awk -v fr=$freq '!/^#|^$/ { print fr"\t"$1"\t"$9}' $f | sort -g -k3| tail -n1 >>$path_out ; done
	
	cd ..;
done

