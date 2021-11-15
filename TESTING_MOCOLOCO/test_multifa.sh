#!/bin/bash

#-------USAGE----------------------------------------------------------------------------------------------

# To run the test.sh file run one of those commands:

# 1) If you want to make implants at 10 different frequencies:

# > bash test_multifa.sh -j <JASPAR_MATRIX> -f <filename.txt> -n <n_seq> -l <length_seq> -p <implanting_pos> -k <k-mers_analysis> -d <hamming_distance> -t <freq1_threshold> -v <Log10p_value_threshold>

# Where:

# -j = Jaspar matrix implanted
# -f = Hitted position output filename
# -n = Number of random multifasta sequences created
# -l = Length of random sequences created
# -p = Position of implants
# -k = K-mers analyzed
# -h = Hamming distance
# -t= Frequency1 threshold 
# -v = Filtering hits on Log10 p-values

# 2) If you want to make a set of random multifasta files:

# > bash test_multifa.sh -f <filename.txt> -n <n_seq> -l <length_seq> -k <k-mers_analysis> -d <hamming_distance> -t <freq1_threshold> -v <Log10p_value_threshold>

# -f = Hitted position output filename
# -n = Number of random multifasta sequences created
# -l = Length of random sequences created
# -k = K-mers analyzed
# -h = Hamming distance
# -t = Frequency1 threshold 
# -v = Filtering hits on Log10 p-values

usage() { echo "Usage: $0 -j <JASPAR_MATRIX> -f <Output_file> -n <Sequences_number> -l <Sequences_length> -p <Implanting_position> -k <kmers_analyzed> -d <Hamming_distance> -c <Cycle_number> -t <frequency1 threshold> -v <Log10 p-value threshold> -b <Bed_file> -tb <Twobit_file>" 1>&2; exit 1; }

#--------PARSER----------------------------------------------------------------------------------------------

while getopts ":j:f:n:l:c:p:k:t:b:x:v:d:ra" o; do
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
        b)
            B=${OPTARG}
            ;;
        x)
            X=${OPTARG}
            ;;
        t)
            T=${OPTARG}
            ;;
        v)
            V=${OPTARG}
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

if [ -z "${F}" ] || [ -z "${N}" ] || [ -z "${T}" ] || [ -z "${L}" ] || [ -z "${K}" ] || [ -z "${D}" ] || [ -z "${C}" ]; then
	
	usage
fi

#----------DIRECTORIES PREPARATION-------------------------------------------------------------------------
                                                                                                             
RMC=$(realpath RMC)
MOCO=$(realpath MOCO)
#Creation of the path for the output file with the best pvalue for each cycle in each frequence
path_out=$(realpath $F)
#Creation of the path for the output file with all the pvalue obtained by our script
path_out_tot=${path_out::-4}_tot.txt
#Defining output directory name (<Jaspar_name>_test_k<k>)

if [ -z "$J" ]
then
	Out_dir="MIXED_BED_noImplant"

else

	Out_dir=${J#../*/};
fi

if [ -z "${Refine}" ]
then

	mkdir ${Out_dir}_test_k${K}_c${C}_f${T};

else	

	mkdir ${Out_dir}_test_k${K}_c${C}_f${T}_r;

fi 

#Initializing headers in Output file (Outside directory just created before) and in the Output file with best pvalue
echo "#TESTING MOCOLOCO">$path_out;
echo -e "FREQ\tHIT\tPVAL-LOG10\tPVAL">>$path_out; 
echo "#TESTING MOCOLOCO">$path_out_tot;
echo -e "FREQ\tHIT\tPVAL-LOG10\tPVAL">>$path_out_tot; 

#Defining frequences for analysis
frequenze=(75 65 55 45 35 25 20 15 10 5);

if [ -z "${Refine}" ]
then

	cd ${Out_dir}_test_k${K}_c${C}_f${T};

else	

	cd ${Out_dir}_test_k${K}_c${C}_f${T}_r;

fi 

for freq in ${frequenze[@]}
do 
	mkdir $freq;
done

#--------RANDOM & RANDOM IMPLANTED FILE CREATION--------------------------------------------------------------
	
for freq in ${frequenze[@]}
do 
	cd $freq;
		
		if [ -z "$J" ]
		then		
		
			$RMC -n $N -l $L -o $freq -c $C -b ../../${B} -t ../../${X} &
		
		else
			
			$RMC -b ../../${B} -t ../../${X} -n $N -l $L -j ../../${J} -p $P -o $freq -c $C &
		fi
  			
	cd ..;
done
wait

#-------RUNNING MOCOLOCO ON IMPLANTED MULTIFASTA-------------------------------------------------------------

multi_thread=20
(
for freq in ${frequenze[@]}
do 
	cd $freq;

	for j in $(seq 1 $C);
	do
        	((i=i%multi_thread)); ((i++==0)) && wait
		
		if [ -z "$J" ]
		then		
		
			$MOCO -m BED_${j}.fasta  $Refine -k $K -d $D -f $T $all &
		
		else

			$MOCO -m BED_${j}.fasta  $Refine -k $K -d $D -f $T $all &
                
		fi
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
	awk -v fre=$freq -v p_val=$V  '!/^#|^$/ { if($9>p_val) print fre "\t"$1"\t"$9"\t"$8 } '  *Z_scores_* >> $path_out_tot;
	#Extraction of the best pvalue for each cycle in each frequence
	for f in *Z_scores* ; do awk -v fr=$freq -v p_val=$V '!/^#|^$/ { if($9>p_val) print fr"\t"$1"\t"$9"\t"$8 }' $f | sort -g -k3| tail -n1 >>$path_out ; done
	
	cd ..;
done

