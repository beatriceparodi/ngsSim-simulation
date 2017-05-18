## ngsSim simulation
### A tutorial for ngsSim data simulation 

ngsSim will allow you to simulate NEXT-GENERATION SEQUENCING (NGS) data for up to 3 populations, taking into account different variables (number of sites, sequencing error rates, sequencing or individual depths, FST value).
Also, using the command -expand you will be able to simulate population expansion.

ngsSim will outcomes different files:
* true genotype files for the whole population and n files for each population 
* reads files 
* genotypes likelihoods

#### Settings directories

1. First of all, **make sure that you have set your directories**. This will allow you to easily understand where your data and your results are. For instance, my directories are: 


`mkdir Data`

`mkdir Results`


2. Now we need to set directories for all programs that we are going to use in this simulation. 
Again, **make sure that you know where these programs are installed** in your computer. My paths are, for example:


`NGSTOOLS=~/Software/ngsTools`

`SAMTOOLS=~/Software/samtools-1.4.1/samtools`

`ANGSD=~/Software/angsd`


3. Now, what kind of data do we want to simulate?

Here my parameters:

* 3 populations
* 10 individuals per population
* 10000 indipendent sites
* for each population, 50% 5x coverage depth samples and 50% 2x coverage depth samples
* FST (fixation index) where population A and B are more closely related (0.1) and  population C is the more diverse (0.3)
* 0.01 error rate
* 0.005 minimun allele frequency
* seed number 12345
* 0.25 frequencies per each base (A,C,G,T)
* output files with prefix pop


#### Start the simulation:
  
1. First, set the sites number

`NSITES=10000`


2. Second, what about the depth?
If you want to generate data with all the same depth, for instance 4x, you`ll have -depth 4. Otherwise, if you want to generate data with individual depths per line, here what you need to do:

  * Set a path for your depths with 

`DEPTH=Data/depths.txt`

  * Create a file with individual depths per line. In this example, I want to generate a file with repetition of 2x and 5x per 30 individuals.

`declare -a seq1=(2 5 2 5 2 5)`

`declare -a repeat=(5 5 5 5 5 5)`

`tLen=${#seq1[@]}`



 `for l in `seq 0 $tLen`

`do`

	`for i in `seq 1 ${repeat[l]}`

	`do`

		echo ${seq1[l]} >> $DEPTH

	`done`

`done`

*Alternatively, for a small amount of data, you can simply create a .txt file and manually insert the depths per individual, one each line*


3. Run the simulation with ngsSim, and **make sure to remember your directories**

`$NGSTOOLS/ngsSim/ngsSim -outfiles Data/pop -npop 3 -nind 10 10 10 -nsites $NSITES -depth $DEPTH -errate 0.01 -pvar 0.10 -mfreq 0.005 -F 0.1 0.3 -model 1 -base_freq 0.25 0.25 0.25 0.25 -seed 12345 2> /dev/null`


*You can find all the parameters that you might want to use for your simulation in the following table:*

Parameter | Usage
----------| -----
`-outfiles POP` | Prefix to outfiles
`-npop INT` | Number of populations (**MUST** be before -nind)
`-nind INT` | Number of individuals per population
`-nsites INT` | Number of sites
`-errate INT` | Sequencing error rate
`-depth INT or FILE` | sequencing depth or file
`-pvar FLOAT` | Probability that a site is variable in a population
`-mfreq FLOAT` | Minimum population frequency
`-F FLOAT` | FST value
`-model INT` | 0 (fixed rate) or 1 (variable rate)
`-simpleRand INT` | Boolean, binary digits 0 (false) or 1 (true)
`-base_FREQ FLOAT FLOAT FLOAT` | Allele frequency for bases
`-seed INT` | Number (or vector) 
`-expansion INT` | Simulate population expansion, binarty digits 0 (false) or 1 (true)

 
4. Create a simulate sequence with an index

`perl -s -e 'print(">chrSIM\n".("A"x$n)."\n");' -- -n=$NSITES > Data/ref.fasta`

`$SAMTOOLS faidx Data/ref.fasta`

5. Now that we have a dataset, we want to investigate it using `ANGSD`

`$ANGSD/angsd -glf Data/pop.glf.gz -fai Data/ref.fasta.fai -nInd 30 -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 1 -doSaf 1 -anc data/ref.fasta -out testA`

*Check the table below for parameters*


Parameter | Usage
----------| -----
`-nInd`| Number of individuals
`-doMajorMinor`| Specification of how to assign the major and minor alleles
`-doMaf`| Estimate allele frequencies
`-doPost` | Calculate posterior probability of genotypes
`-doGeno` | Assign genotype probabilities at each site for each individual
`-doSaf` | Sample allele frequency based analysis


## Principal Component Analysis (PCA)
### How to perform a PCA with your data
PCA is a useful tool when analysing genetic low-depth data.
This statistical method will allow you to reduce your measurements into few principal components (PCs) that will explain the main pattern in your dataset (Reich, Prince & Patterson, 2008).

*You`ll find more detailed info about PCA in this link* [http://www.nature.com/ng/journal/v40/n5/full/ng0508-491.html]

1. To perform a PCA we first need a set of data.
You can use data that you have collected or you can simulate your own dataset with ngsSim.
Recalling what we have done so far, here some more in-depth detailes about ANGSD commands.


* `-doGeno` will assign genotype probabilities at each site for each individuals.
Different options are available for this:
          
 * `-doGeno1` print major and minor alleles
 * `-doGeno2` print called genotype encoded (-1,0,1,2)
 * `-doGeno4` print called genotype directly (AA,AC,AG)
 * `-doGeno8` posterior probability of all possible genotypes
 * `-doGeno16` posterior probability of called genotype
 * `-doGeno32` posterior probability of called genotype as binary

*You can combine the output by summing the numbers. For instance,* `-doGeno9`*(1+8) will gives you the major and minor alleles and the probability for the genotypes [(MM), (Mm), (mm)]*


* `-doPost` will calculate the posterior probability of genotypes defined by a choosen model. Again different options are available:

 * `-doPost1` will use frequency as a prior
 * `-doPost2` will use uniform prior


* `-doMaf` will estimate  minor allele frequencies. You can estimate allele frequencies by:

 * `-doMaf1` calculate the frequency of major and minor allele, based on EM     algorithm with genotypes likelihood
 * `-doMaf2` calculate frequency with only fixed major and unknown minor allele
    `-doMaf4` calculate the frequencies directly from genotype posterior probabilities
    `-doMaf8` calculate frequencies based in allele counts


* `-doMajorMinor` will nallow you to assign the major and minor alleles. Options are available:

 * `-doMajorMinor1` major and minor inferred from GL
 * `-doMajorMinor2` major and minor inferred from allele counts
 * `-doMajorMinor3` use major and minor specified in a file (sites- FILE)
 * `-doMajorMinor4` reference allele as major (ref- FILE)
 * `-doMajorMInor5` ancestral allele as major (anc- FILE)
       

* `-GL` specify the genotype likelihood model that you are interested in. In this case, we are going to use SAMtools, with `-GL1`
    
2. Before running the analyses remember to unzip the files that you want to use (in my case the file with the whole population) with `gunzip`: 

`gunzip Data/pop.glf.gz`

`$ANGSD/angsd -P 3 -b Data/pop.glf -ref Data/ref.fasta.fai -out results/all  -r chrSIM:1-100 -GL 1 -doMajorMinor 1 -doMaf 1 -doGeno 1 -doPost 1`

















