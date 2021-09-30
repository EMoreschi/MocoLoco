# MOCOLOCO TEST REPORT
<br>
To validate and analyze the MocoLoco tool functioning we performed a test on 3 different Transcription Factor binding sites (coming from chip-seq analysis) and on a Multifasta random sequences with specific positional implants, generated from <i>MocoLoco_random_tool</i>.<br><br>
Transcription Factors bed:
<ul>
<li> NFY-A
<li> NRF1
<li> ADNP
<br>
</ul>
Multifasta:
<ul>
<li> 12 files composed by 3000 random sequences of 300 base length
</ul><br>

The bed files analized come from chip-seq peak experiments on k562 cell line. The DNA bed coordinates refer to hg19(for NFY) and hg38(for NRF1 and ADNP) genome assembly. The sequences (of 300 base length) have been extracted from hg19 and hg38 twobit format, following bed coordinates. <br>
The multifasta file has been generated from <i>MocoLoco_Random_tool</i>. The tool has been run 12 times with the following parameters:
<ul>
- Sequences number: 3000 <br>
- Sequences length: 300 <br>
- Implanting position: 150 <br>
- Implanted Jaspar: MA000.0_prova.jaspar (which produces GGGGGGGG string 100% times) <br>
- FWD/REV Jaspar frequences: 50% GGGGGGGG and 50% CCCCCCCC <br>
- Implanting frequences: Different for each files:<br>
<ul>100% | 50% | 40% | 30% | 20% | 15% | 10% | 9% | 8% | 7% | 6% | 5% <br>
</ul></ul><br>

The test has been performed running the bash script _test_bed.sh_ for bed files input and the script _test_multifa.sh_ for multifasta files input. <br><br>
For each file in input 6 different test have been performed:
<ul>
1. No -r and No -a options with k = 6 and d = 1 <br>
2. No -r and No -a options with k = 8 and d = 2 <br>
3. -r and No -a options with k = 6 and d = 1 <br>
4. -r and No -a options with k = 8 and d = 2 <br>
5. -r and -a options with k = 6 and d = 1 <br>
6. -r and -a options with k = 8 and d = 2 <br><br>
<b>A frequence1 threshold of 0.02 has been set for all the tests.<br>
Each position has a frequence1 associated with it and it is an index of how much a set of oligos is present with respect to the others within all sequences.</b>

$$ Frequence1 = \frac{best.oligo.occurrences + hamming.neighbours.occurrences}{Tot.sequences.number} $$
<br>

</ul>
MocoLoco tool produced 6 output files for each bed in analysis (18 bed output files) and 6 output files for each different implanting frequences of random multifasta (72 multifasta output files). <br><br>
Each output file contains information about:
<ul>
<li> Implanting frequence (for Multifasta)
<li> Hit position
<li> Best oligo of the position
<li> Log10 of p-value calculated from K,N1,N2,T parameters
<li> Log10 of p-value calculated from Z-test
</ul>
Then plots of the 3 bed files and 3 selected multifasta frequences (100% | 20% | 5%) has been generated using <i>R Studio</i>. <br>
To generate the graphs the Log10 p-values coming from Z-test have been plotted on the Hit positions.<br><br>

## NO -r NO -a analysis with k = 6 and d = 1
<br>
<img src='./noRnoA_k6.png'>
<i>Image 1</i>
<br><br>
The <i>Figure 1</i> represents the result of NFY hits. We can observe clearly that, without -r and -a options, the hits which overcome the threshold of 0.02 for frequence1 are very few and they correspond to central positions, where there are NFY chip-seq peaks. <br>
The same thing happens for NRF1, ADNP Transcription Factors and also for Multifasta. <br>
P-values associated to this hit are obviously high because of their location around TFs binding sites and around the implanting position for multifasta. Only the <i>Figure 6</i> shows lower p-value but this is correlated to the low (5%) implanting frequence.<br><br>

## NO -r NO -a analysis with k = 8 and d = 2
<br>
<img src='./noRnoA_k8.png'>
<i>Image 2</i>
<br><br>
Observing the analysis without -r and without -a on k = 8 we see how the behaviour is more or less the same.<br>
The TF graphs have hits only on peak regions (NFR1 does not even show any hit) and for multifasta the hits are in implanting regions and their p-value follows the implanting frequences decrease.<br><br>
We can conclude that, to see some possible TFBSs, we can not analyze only the chip-seq peak parts of the sequences but we need to select more hits. To do that a possible solution is to set the threshold on a lower value or, as we did, compute tests with the -r refine option addition. <br><br>

## YES -r NO -a analysis with k = 6 and d = 1
<br>
<img src='./RnoA_k6.png'>
<i>Image 3.1</i>
<br><br>
Using the option refine -r the hamming neighbours group, for each position, is composed by more oligos and then also the frequence1 value increases. This allows to select more hit than the previous tests. <br>
Since the option -a is still disabled, the hit represented in graph are the local maxima, calculated on frequence1 values.<br>
In <i>Figure 1</i> we can see the NFY peak surrounded by a lot of smaller peaks, meaning that this genomic region is pretty conserved in all the sequences. This high p-values distribution along all the genomic window means that this region is full of putative TF binding sites even if there are no peak which clearly dominates the others (the best ones are around pos.250 and pos.270).<br> Near the peak corresponding to the NFY binding site the p-values are as high as NFY. This happens because in these position there is the CCAAT box, very important for NFY binding.<br>
Observing the <i>Figure 2</i> it is appreciable the high peak corresponding to NRF1 binding site. Differently from <i>Figure 1</i> the surrounding region has very low p-value with minimal peaks around position 10,70,240,290.<br>
The reason of that can be the high ripetitiveness of the region, which is observable pretty good looking at the best oligo found in these positions. <br>
In <i>Figure 3</i>, where ADNP binding region is analyzed, another case is shown. This is a combination of a region full of repetition and two beautiful and appreciable peaks (around pos.100 and 200).<br>
To validate the strength of these peaks we extract the best oligo from the peaks (TGTGACCT) and we put as query motif on <i>Factorbook</i>. The results show a beautiful correlation with the TF CHD4 binding site. <br><br>
<img src='./CHD4_ADNP_correlation.png'>
<i>Image 3.2</i>
<br><br>
Going deeper we tryed to find a biological correlation between the two TFs and we found out that, as written in the article <i>Activity-dependent neuroprotective protein recruits HP1 and CHD4 to control ineage-specifying genes</i> (<A>https://pubmed.ncbi.nlm.nih.gov/29795351/</A>), ADNP, via the recruitment of HP1 and CHD4, regulates the expression of genes that are crucial for maintaining distinct cellular states and assures accurate cell fate decisions upon external signals.<br>
This correlation means that the tool works good, highlighting, near TF chip-seq peak binding site like ADNP, a conserved region which is probably associated to CHD4 binding.<br>
Another proof of the good behaviour of the program can be seen looking at <i>Figure 4</i>, <i>Figure 5</i> and <i>Figure 6</i> of <i>Image 3.1</i>. These analysis are performed on random multifasta and it is appreciable how the implantig site peaks emerge on the random regions (which has not any small peaks, no false positives), also with a low implanting frequency like 5%. <br>
<br>

## YES -r NO -a analysis with k = 8 and d = 2
<br><br>
<img src='./RnoA_k8.png'>
<i>Image 4</i>
<br><br>
Setting the k parameter to 8 the results are pretty similar to the previous ones. <br>
For NFY the behaviour is almost the same, but for NRF1 (<i>Figure 2</i>) we can observe an increase of peaks in position 10,50,240,290.
But looking at oligos which generate these peaks we saw that they are composed only by repetitions (AAAAAAAA or TTTTTTTT) and their high presence generates this strong false positives. <br>
Also looking at <i>Figure 3</i> we found some differences. The most important one is the lowering of the CHD4 peaks, even if they mantain a good distinguishable p-value.
<br><br>

## YES -r YES -a analysis with k = 6 and d = 1
<br><br>
<img src='./RA_k6.png'>
<i>Image 5</i>
<br><br>
This test has been made with both -r and -a parameters activated. The -a allows to consider also the hits which overcome the threshold and are not local maxima. In this case, setting a low threshold like 0.02, it is possible to observe the complete shapes of p-values along the positions.<br>
The results are similar to ones described in <i>Image 3</i> but the can be seen in a more precise way.
<br><br>

## YES -r YES -a analysis with k = 6 and d = 1
<br><br>
<img src='./RA_k8.png'>
<i>Image 6</i>
<br><br>
Finally this image represents the complete p-values shapes of the <i>Image 4</i>.<br>
Also in this case the results are similar to the ones described before. <br><br>

## CONCLUSIONS <br>
The tool generally works good, as proved by the random multifasta behaviour in complitely different states (100% implanting frequency, 20% and 5%). Analyzing the MF result it is possible to assume that the program can distinguish pretty good a random region from a conserved region (also with low conservation) without any false positive case.<br>

The main problem of the tool is the risk to be "trapped" in repetitions. This happens when a repetitive region (for example AAAAAAAA) is present more in a specific part of the genomic window than in the other positions. This situation brings the p-value to get up suddenly because of the high local conservation of the repeats in comparison to its distribution along the sequences.<br>
A way to overcome this problem can be an additional function implementation. My idea is a function that can distinguish if an oligo is a repeat (following criteria previously defined) and, if a certain best oligo for a certain position is classified as a repeat, this function can select as best the following (not repeat) oligo in the positional rank.<br>
Another easier (but i think non so much effective) solution can be to stretch the genomic window for the analysis. It is possible that the high frequency of the local repeats can be shielded increasing the probability to find them in the global region. By the way i think this method does not work, indeed the result can become worst for the amplification of the repetition presence in a specific position of a larger window.<br>
In any case i think that some tests with larger genomic window (600bp) has to be done, maybe also to explore a larger region and find out more interesting putative binding sites. However, before to do that, it is important to reduce the computatonal time of the tool, appling some threading funcions.



