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

`echo $NSITES`

2. Second, what about the depth?
If you want to generate data with all the same depth, for instance 4x, you`ll have -depth 4. Otherwise, if you want to generate data with individual depths per line, here what you need to do:

  * Set a path for your depths with 

`DEPTH=Data/depths.txt`

`echo $DEPTH`

  * Create a file with individual depths per line. In this example, I want to generate a file with repetition of 2x and 5x per 30 individuals.


`declare -a seq1=(2 5 2 5 2 5)`

`declare -a repeat=(5 5 5 5 5 5)`

`tLen=${#seq1[@]}`



 `for l in `seq 0 $tLen`

`do`

	for i in `seq 1 ${repeat[l]}`

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

5. Now that we have the data, we want to investigate them using `ANGSD`

`$ANGSD/angsd -glf $SIM_DATA/testA.glf.gz -fai $SIM_DATA/testAF.ANC.fas.fai -nInd 30 -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 1 -doSaf 1 -anc $SIM_DATA/testAF.ANC.fas -out testA`

*Check the table below for the parameters*


Parameter | Usage
----------| -----
`-nInd`| Number of individuals
`-doMajorMinor`| Specification of how to assign the major and minor alleles
`-doMaf`| Estimate allele frequencies
`-doPost` | Calculate posterior probability of genotypes
`-doGeno` | Assign genotype probabilities at each site for each individual
`-doSaf` | Sample allele frequency based analysis






















