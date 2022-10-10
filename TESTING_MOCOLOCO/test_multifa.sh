#!/bin/bash

# To run the test_multifa.sh file run one of these 4 commands:

#########    1    #########################################################################################

# TO MAKE IMPLANTS AT DIFFERENT FREQUENCIES (75 | 65 | 55 | 45 | 35 | 25 | 20 | 15 | 10 | 5) INTO A SET OF RANDOM MULTIFASTA FILES:

# > ./test_multifa.sh -j <JASPAR_MATRIX> -o <filename.txt> -n <n_seq> -l <length_seq> -p <implanting_pos> -k <k-mers_analysis> -d <hamming_distance> -f <freq1_threshold> -s <strand_+_implanting_frequency> -c <cycles_number> -v <p_value_threshold> -r <optional> -a <optional>

# Where:

# -j (JASPAR) = Jaspar matrix implanted
# -o (OUTPUT FILE) = Hitted position output filename
# -n (N SEQ) = Number of random multifasta sequences created
# -l (L SEQ) = Length of random sequences created
# -p (POSITION) = Position of implants
# -k (K-MERS) = K-mers analyzed 
# -d (HAMMING DISTANCE) = Hamming distance
# -f (FREQ1 THRESHOLD) = Frequency1 threshold
# -s (STRAND +) = Strand + implanting frequency (DEFAULT 50%)
# -c (CYCLES) = Number of multifasta created for each implanting frequency
# -v (P-VALUE FILTER) = Filtering hits on p-values (it takes pval < of v | DEFAULT = 1)
# -r (REFINE) = To refine each hit's PWM MATRIX created
# -a (ALL) = No Local maxima selection

#########    2    #########################################################################################

# TO MAKE A SET OF RANDOM MULTIFASTA FILES (Without any implants)

# > ./test_multifa.sh -o <filename.txt> -n <n_seq> -l <length_seq> -k <k-mers_analysis> -d <hamming_distance> -f <freq1_threshold> -c <cycles_number> -v <p_value_threshold> -r <optional> -a <optional>

# -o (OUTPUT FILE) = Hitted position output filenamec
# -n (N SEQ) = Number of random multifasta sequences created
# -l (L SEQ) = Length of random sequences created
# -k (K-MERS) = K-mers analyzed
# -d (HAMMING DISTANCE) = Hamming distance
# -f (FREQ1 THRESHOLD) = Frequency1 threshold
# -c (CYCLES) = Number of multifasta created for each implanting frequency
# -v (P-VALUE FILTER) = Filtering hits on p-values (it takes pval < of v | DEFAULT = 1)
# -r (REFINE) = To refine each hit's PWM MATRIX created
# -a (ALL) = No Local maxima selection

#########    3    #########################################################################################

# TO MAKE IMPLANTS AT DIFFERENT FREQUENCIES (75 | 65 | 55 | 45 | 35 | 25 | 20 | 15 | 10 | 5) INTO A SET OF SEQUENCES COMING FROM BED COORDINATES:

# > ./test_multifa.sh -b <BED_FILE> -t <TWOBIT_FILE> -j <JASPAR_MATRIX> -o <filename.txt> -n <n_seq> -l <length_seq> -p <implanting_pos> -k <k-mers_analysis> -d <hamming_distance> -f <freq1_threshold> -s <strand_+_implanting_frequency> -c <cycles_number> -v <p_value_threshold> -r <optional> -a <optional>

# Where:

# -b (BED) = Bed file
# -t (TWOBIT) = Twobit file
# -j (JASPAR) = Jaspar matrix implanted
# -o (OUTPUT FILE) = Hitted position output filename
# -n (N SEQ) = Number of random multifasta sequences created
# -l (L SEQ) = Length of random sequences created
# -p (POSITION) = Position of implants
# -k (K-MERS) = K-mers analyzed
# -d (HAMMING DISTANCE) = Hamming distance
# -f (FREQ1 THRESHOLD) = Frequency1 threshold 
# -s (STRAND +) = Strand + implanting frequency (DEFAULT 50%)
# -c (CYCLES) = Number of multifasta created for each implanting frequency
# -v (P-VALUE FILTER) = Filtering hits on p-values (it takes pval < of v | DEFAULT 1)
# -r (REFINE) = To refine each hit's PWM MATRIX created
# -a (ALL) = No Local maxima selection

#########    4    #########################################################################################

# TO MAKE A SET OF MULTIFASTA FILES TAKEN FROM BED COORDINATES (Without any implants)

# > ./test_multifa.sh -b <BED_FILE> -t <TWOBIT_FILE> -o <filename.txt> -n <n_seq> -l <length_seq> -k <k-mers_analysis> -d <hamming_distance> -f <freq1_threshold> -c <cycles_number> -v <p_value_threshold> -r <optional> -a <optional>

# -b (BED) = Bed file
# -t (TWOBIT) = Twobit file
# -o (OUTPUT FILE) = Hitted position output filename
# -n (N SEQ) = Number of random multifasta sequences created
# -l (L SEQ) = Length of random sequences created
# -k (K-MERS) = K-mers analyzed
# -d (HAMMING DISTANCE) = Hamming distance
# -f (FREQ1 THRESHOLD) = Frequency1 threshold 
# -c (CYCLES) = Number of multifasta created for each implanting frequency
# -v (P-VALUE FILTER) = Filtering hits on p-values (it takes pval < of v | DEFAULT 1)
# -r (REFINE) = To refine each hit's PWM MATRIX created
# -a (ALL) = No Local maxima selection

############################################################################################################

#--------USAGE----------------------------------------------------------------------------------------------

usage() { 

echo -e "\nWRONG PARAMETERS INSERTION!"
echo -e "Run one of these 4 commands"

echo -e "\n1) $0 -j <JASPAR_MATRIX> -o <Output_file> -n <Sequences_number> -l <Sequences_length> -p <Implanting_position> -k <kmers_analyzed> -d <Hamming_distance> -f <frequency1 threshold> -s <strand_+_implant_freq> -c <cycles_number> -v <p-value threshold> -r <optional> -a <optional>" 1>&2;

echo -e "\n2) $0 -o <Output_file> -n <Sequences_number> -l <Sequences_length> -k <kmers_analyzed> -d <Hamming_distance> -f <frequency1 threshold> -c <cycles_number> -v <p-value threshold> -r <optional> -a <optional>" 1>&2; 

echo -e "\n3) $0 -b <BED_FILE> -t <TWOBIT_FILE> -j <JASPAR_MATRIX> -f <Output_file> -n <Sequences_number> -l <Sequences_length> -p <Implanting_position> -k <kmers_analyzed> -d <Hamming_distance> -f <frequency1 threshold> -s <strand_+_implant_freq> -c <cycles_number> -v <p-value threshold> -r <optional> -a <optional>" 1>&2; 

echo -e "\n4) $0 -b <BED_FILE> -t <TWOBIT_FILE> -f <Output_file> -n <Sequences_number> -l <Sequences_length> -k <kmers_analyzed> -d <Hamming_distance> -f <frequency1 threshold> -c <cycles_number> -v <p-value threshold> -r <optional> -a <optional> \n" 1>&2; exit 1; 

}

#--------PARSER----------------------------------------------------------------------------------------------

while getopts ":b:t:j:o:n:l:p:k:d:f:s:c:v:e:ra" o; do
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
        o)
            O=${OPTARG}
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
        f)
            F=${OPTARG}
            ;;
        s)
            S=${OPTARG}
            ;;
        c)
            C=${OPTARG}
            ;;
		e)
			E=${OPTARG}
			;;
        v)
            V=${OPTARG}
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

if [ -z "${O}" ] || [ -z "${N}" ] || [ -z "${L}" ] || [ -z "${K}" ] || [ -z "${D}" ] || [ -z "${F}" ] || [ -z "${C}" ]; then
	
	usage
fi

#----------DIRECTORIES PREPARATION & INPUT PARAMETERS CONTROL-------------------------------------------------------------------------   

RMC=$(realpath RMC)
MOCO=$(realpath NewMoco)
#Creation of the path for the output file with the best pvalue for each cycle in each frequence
path_out=$(realpath $O)
#Creation of the path for the output file with all the pvalue obtained by our script
path_out_tot=${path_out::-4}_tot.txt
#Defining output directory name (<Jaspar_name>_test_k<k>)

if [ -z "$V" ]
then
	V=1
fi

if [ -z "$S" ]
then
	S=50
fi

if [ -z "$B" ]
then
	if [ -z "$J" ]
	then
		Out_dir="RANDOM_noImplant"

	else

		Out_dir=${J#../*/};
	fi
else
	if [ -z "$J" ]
	then
		Out_dir="MIXED_BED_noImplant"
	
	else

		Out_dir="MIXED_BED"${J#../*/};
	fi
fi

if [ -z "$E" ]
then

	mkdir ${Out_dir}_test_k${K}_c${C}_f${F};

else	

	mkdir ${Out_dir}_test_k${K}_c${C}_f${F}_e${E};

fi 
#Initializing headers in Output file (Outside directory just created before) and in the Output file with best pvalue
echo "#TESTING MOCOLOCO">$path_out;
echo -e "FREQ\tHIT\tPVAL\tPVAL-LOG10\tBONF-PVAL\tLOG10BONF">>$path_out; 
echo "#TESTING MOCOLOCO">$path_out_tot;
echo -e "FREQ\tHIT\tPVAL\tPVAL-LOG10\tBONF-PVAL\tLOG10BONF">>$path_out_tot; 

#Defining frequences for analysis
frequenze=(75 65 55 45 35 25 20 15 10 5);

#frequenze=(9 8 7 6);

if [ -z "$E" ]
then

	cd ${Out_dir}_test_k${K}_c${C}_f${F};

else	

	cd ${Out_dir}_test_k${K}_c${C}_f${F}_e${E};

fi 

for freq in ${frequenze[@]}
do 
	mkdir $freq;
done

#--------RANDOM & RANDOM IMPLANTED FILE CREATION--------------------------------------------------------------
	
for freq in ${frequenze[@]}
do 
	cd $freq;
		
		if [ -z "$B" ]
		then		
			if [ -z "$J" ]
			then

				$RMC -n $N -l $L -o $freq -c $C &
		
			else
			
				$RMC -j ../../${J} -n ${N} -l ${L} -p $P -o $freq -c $C -f $S &
			fi
		else
			if [ -z "$J" ]
			then		
		
				$RMC -b ../../${B} -t ../../${T} -n $N -l $L -o $freq -c $C &
		
			else
			
				$RMC -b ../../${B} -t ../../${T} -j ../../${J} -n ${N} -l ${L} -p $P -o $freq -c $C -f $S &
			fi
  		fi
	cd ..;
done
wait

#-------RUNNING MOCOLOCO ON IMPLANTED MULTIFASTA-------------------------------------------------------------

multi_thread=40
(
for freq in ${frequenze[@]}
do 
	cd $freq;

	for j in $(seq 1 $C);
	do
        	((i=i%multi_thread)); ((i++==0)) && wait
		
		if [ -z "$B" ]
		then		
			if [ -z "$J" ]
			then
				$MOCO -m random_multifa_${j}.fasta -k $K -d $D -f $F $all &
		
			else

				$MOCO -m random_multifa_implanted${j}.fasta -k $K -d $D -f $F $all &
                
			fi
		else
		
			if [ -z "$E" ]
			then	
				if [ -z "$J" ]
				then	
		
					$MOCO -m BED_${j}.fasta -k $K -d $D -f $F $all &

				else

					$MOCO -m BED_implanted${j}.fasta -k $K -d $D -f $F $all &
				fi	
			else

				if [ -z "$J" ]
				then		
		
					$MOCO -m BED_${j}.fasta -k $K -d $D -f $F $all -e $E &
		
				else

					$MOCO -m BED_implanted${j}.fasta -k $K -d $D -f $F $all -e $E &
                		fi
			fi
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
	

	if [ -z "$V" ]
	then
		V=1
	fi

	#Extraction of all the pvalue from the Z_scores_implanted files
	awk -v fre=$freq -v p_val=$V  '!/^#|^$/ { if($8<p_val) print fre "\t"$1"\t"$8"\t"$9"\t"$10"\t"$11 } '  *Z_scores_* >> $path_out_tot;
	#Extraction of the best pvalue for each cycle in each frequence
	for f in *Z_scores* ; do awk -v fr=$freq -v p_val=$V '!/^#|^$/ { if($8<p_val) print fr"\t"$1"\t"$8"\t"$9"\t"$10"\t"$11 }' $f | sort -g -k4| tail -n1 >>$path_out ; done
	
	cd ..;
done

