# **P_VALUE TEST**

## Testing P_value on different implanting frequencies
<br>

### **How to run**

    bash test_p_val.sh <Jaspar Matrix> <Implanting position> <Sequences number> <Sequences length> <K-mer length> <Output_filename>

Example:

    bash test_p_val.sh MA0060.1.jaspar 50 500 300 8 Output_example.txt

### **Output**

The script generates an output files containing:<br> 
1. the implanting frequency percentage<br>
2. the first ranked line in the implant starting position. The line has been extracted from <i><u>k-mers_positional_occurrenced__DS.txt</i></u> file, produced by <i><u>MocoLoco.cpp</i></u>.<br><br>

Output lines example:<br>

##### <li> 14%     43      1       TGCATCGG        2       0       2               CCGATGCA        0       2       2       FALSE   2       0.002   7.84345e-07 <br>
##### <li> 15%     43      1       GTTCTGAT        2       0       2               ATCAGAAC        0       2       2       FALSE   2       0.00200401      3.93173e-07 <br>
</ul>
<br>

### **Notes**
The script <i><u>test_p_val.sh</i></u> implants sequences from Jaspar Matrix Oligo generation using <i><u>Multifa_random_TOOL.cpp</i></u> tool.<br>
The implanting position provided as input indicates where implanted oligos are centered. <br>
For example if a Jaspar M. generates 11 bases length oligos, the position provided is 50 and K is 7, Multifa_random_TOOL.cpp ensures analysis starts on 47th position.<br>
