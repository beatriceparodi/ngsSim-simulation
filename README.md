## ngsSim simulation
### A tutorial for ngsSim data simulation 

ngsSim will allow you to simulate NEXT-GENERATION SEQUENCING (NGS) data for up to 3 populations, taking into account different variables (number of sites, sequencing error rates, sequencing or individual depths, FST value).
Also, using the command -expand you will be able to simulate population expansion.

ngsSim will outputs different files:
* true genotype files for the whole population and n files for each population 
* reads files 
* genotypes likelihoods

#### Settings directories

1. First of all, **make sure that you have set your directories**. This will allow you to easily understand where data and results are. For instance, my directories are: 


`mkdir Data`

`mkdir results`


2. Now we need to set directories for all programs that we are going to use in this simulation. 
Again, **make sure that you know where these programs are installed** in your computer. My paths are, for example:


`NGSTOOLS=~/Software/ngsTools`

`SAMTOOLS=~/Software/samtools-1.4.1/samtools`

`ANGSD=~/Software/angsd`

`NGSLD=~/Software/ngsLD`

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
If you want to generate data with all the same depth, for instance 4x, you`ll have `-depth 4` . Otherwise, if you want to generate data with individual depths per line, here what you need to do:

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
This statistical method will allow you to reduce your measurements into few principal components (PCs) that will explain the main pattern in your dataset. It is a linear transformation that chooses a new system of coordinates from your original dataset, and creates new axis called principal components.

You`ll find more detailed info about PCA [here](http://www.nature.com/ng/journal/v40/n5/full/ng0508-491.html).

1. To perform a PCA we first need a set of data.
You can use data that you have collected or you can simulate your own dataset with ngsSim.
Recalling what we have done so far, here some more in-depth detailes about ANGSD commands.

`-doGeno` will assign genotype probabilities at each site for each individuals.
Different options are available for this:
          
 * `-doGeno1` print major and minor alleles
 * `-doGeno2` print called genotype encoded (-1,0,1,2)
 * `-doGeno4` print called genotype directly (AA,AC,AG)
 * `-doGeno8` posterior probability of all possible genotypes
 * `-doGeno16` posterior probability of called genotype
 * `-doGeno32` posterior probability of called genotype as binary

You can combine the output by summing the numbers. For instance, `-doGeno9`(1+8) will gives you the major and minor alleles and the probability for the genotypes [(MM), (Mm), (mm)]

`-doPost` will calculate the posterior probability of genotypes defined by a choosen model. Again different options are available:

 * `-doPost1` will use frequency as a prior
 * `-doPost2` will use uniform prior

`-doMaf` will estimate  minor allele frequencies. You can estimate allele frequencies by:

 * `-doMaf1` calculate the frequency of major and minor allele, based on EM     algorithm with genotypes likelihood
 * `-doMaf2` calculate frequency with only fixed major and unknown minor allele
 * `-doMaf4` calculate the frequencies directly from genotype posterior probabilities
 * `-doMaf8` calculate frequencies based in allele counts

`-doMajorMinor` will nallow you to assign the major and minor alleles. Options are available:

 * `-doMajorMinor1` major and minor inferred from GL
 * `-doMajorMinor2` major and minor inferred from allele counts
 * `-doMajorMinor3` use major and minor specified in a file (sites- FILE)
 * `-doMajorMinor4` reference allele as major (ref- FILE)
 * `-doMajorMInor5` ancestral allele as major (anc- FILE)
       
`-GL` specify the genotype likelihood model that you are interested in. In this case, we are going to use SAMtools, with `-GL1`
    
2. Before running the analyses remember to unzip the files that you want to use (for instance,I'm going to use the whole population) with `gunzip`: 


`gunzip Data/pop.glf.gz`


and then we can calculate the genotype probabilities. In this specific case, since we need a binary file for the PCA, we must specify `-doGeno32`. So the command line will be:


`$ANGSD/angsd -glf Data/pop.glf -fai Data/ref.fasta.fai -out results/all -nInd 30 -r chrSIM:1-100 -doMajorMinor 1 -doMaf 1 -doGeno 32 -doPost 1`


If we need to filter sites that are actually variable (for instance, replacement of a base with another) in our samples, wee need to add SNP_pval, as shown below:
 

`$ANGSD/angsd -glf Data/pop.glf -fai Data/ref.fasta.fai -out SNP/SNPall -nInd 30 -r chrSIM:1-100 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-3 -doGeno 32 -doPost 1`


You can use the command `-r` if you want to filter the data. In my case I am only analysing chromosome 1-100

3. ngsCovar will allow us to estimate the covariance matrix between individuals based on genotypes probabilities. Since the data are simulated, I am not using any filter, and I can type the initial number of sites (NSITES=10000).

Before running the ngsCovar, we need to gunzip our geno.gz files with

`gunzip results/ALL.geno.gz`

and then we can calculate the covariance matrix with 

`$NGSTOOLS/ngsPopGen/ngsCovar -probfile results/ALL.geno -outfile results/matrix.covar -nind 30 -nsites 10000 -call 0 -norm 0 &> /dev/null`


If we performed a SNP calling, we need to know how many sites we are taking into consideration (581 in my case). We can see them  with these commands


`less -S SNP/SNPall.mafs.gz`
`N_SITES=`zcat SNP/SNPall.mafs.gz | tail -n+2 | wc -l``
`echo $N_SITES`


Then we can run the matrix calculation


`$NGSTOOLS/ngsPopGen/ngsCovar -probfile SNP/SNPall.geno -outfile SNP/SNPmatrix.covar -nind 30 -nsites $N_SITES -call 0 -norm 0 &> /dev/null`


Here,looking at the outup, we can see a NxN symmetric matrix with N individials 

4. Now that we have our covariance matrix, we can use the command `Rscript` to perform a PCA plot in R directly from our terminal.
With the following commands, we are going to transform our matrix into a canonical form, then create a cluster file (.clst) and plot the results.


`Rscript -e 'write.table(cbind(seq(1,30),rep(1,30),c(rep("POP1",10),rep("POP2",10),rep("POP3",10))), row.names=F, sep=" ", col.names=c("FID","IID","CLUSTER"), file="results/ALL.clst", quote=F)'`


`Rscript $NGSTOOLS/Scripts/plotPCA.R -i results/matrix.covar -c 1-2 -a results/ALL.clst -o results/ALL.pca.pdf`


`evince results/ALL.pca.pdf`


 *We can run the same commands changing the input and output files if we performed the SNP calling*

[Here](https://github.com/beatriceparodi/ngsSim-simulation/blob/master/ALL.pca.pdf) my plot without filters 

[Here](https://github.com/beatriceparodi/ngsSim-simulation/blob/master/SNPALL.pca.pdf) the plot with the SNP calling


## Admixture proportions
### How to calculate admixture proportions with NGSadmix 
Admixture occurs while populations begin interbreeding, and their offsprings represent a mixture of alleles from different ancestral populations. 
[Here](http://www.genetics.org/content/195/3/693) a useful paper.

1. We are going to use NGSadmix, so first of all we`ll set our new directories with


`$NGSADMIX=~/Software/NGSadmix`


2. NGSadmix require BEAGLE files as input format. We can create them with ANGSD
 

`$ANGSD/angsd -glf Data/pop.glf -fai Data/ref.fasta.fai -out results/admix -nInd 30 -r chrSIM:1-100 -doMajorMinor 1 -doGlf 2 -doMaf 1  -snp_pval 1e-3 &> /dev/null`


3. Type `$NGSADMIX` and make sure that you are familiar with filters and options available before running the admixture analysis. It will outputs arguments, setups, filters and other option that you may want to use in your command.

4. Now, let's assume that we want to test our dataset for 3 ancestral components. We first need to set these with


`K=3`


and then we run the analysis for our admixtur proportions


`$NGSADMIX -likes results/admix.beagle.gz -K $K -outfiles results/admixall -minMaf 0.0001 -seed 1 -minInd 3 &> /dev/null`


4. This analysis will outputs us 3 different files:
* admixall.log that summarises the analysis
* admixall.fopt.gz is a compressed file containing the allele frequency in each of the 3 ancestral populations
* admixall.qopt with a line with ancestry proportion for each individual

5. Plot the admixture proportion estimates. We can use R directly from our terminal typing


`R`


and then plotting our barplot for admixture proportion with the following commands:


`r=read.table("results/admixall.qopt")`


`barplot(t(as.matrix(r,nrow=30)),beside=F , col=rainbow(3), xlab= "Individual #" , ylab= "Ancestry")`


Finally, saving the plot in our results


`pdf("./results/barplot.pdf")`


`barplot(t(as.matrix(r,nrow=30)),beside=F, col=rainbow(3), xlab= "Individual #" , ylab= "Ancestry")`


`dev.off()`

You can see my plot in this [link](https://github.com/beatriceparodi/ngsSim-simulation/blob/master/barplot.pdf)


## Summary statistics: Site Frequency Spectrum (SFS)
### How to estimate the proportions of sites at different allele frequencies with ANGSD
Site frequency spectrum is the distribution of the alelle frequencies of a given set of loci in a population, and represent one of the most powerful method for summarising genomic data. 
More info about SFS in this [link](https://www.ncbi.nlm.nih.gov/pubmed/23770700).

With the following commands, we are going to use ANGSD and a simulated dataset with 3 populations of 10 individuals each.

1. First of all, we are going to perform our analysis for each population separately with 


`for POP in pop1 pop2 pop3`

`do`

    `echo $POP`

    `$ANGSD/angsd -glf Data/$POP.glf -fai Data/ref.fasta.fai -out results/$POP    -nInd 10 -doglf 1 -anc Data/ref.fasta -doSaf 1 &> /dev/null`

`done`


Take a look at the output (for instance at pop1) with


`$ANGSD/misc/realSFS print results/sfspop1.saf.idx | less -S`


You should see a table with different values. For example my output is 

**chrSIM  113  |  -78.105034   |  -56.263031   |  -44.215683  |  -33.021896**

These represent the allele frequency likelihoods per each site. In my case, for chromosome 113,the first value represent the likelihood of having 0 copies of the derived allele, the second value the likelihood of having 1 copy and so on..

2. Now, we are going to perform our SFS with `realSFS`. 
Type `$ANGSD/misc/realSFS` and make sure that you have understood all the options available for this program.

Then, run the analysis with


`for POP in pop1 pop2 pop3`

`do`

   `echo $POP`
    
    `$ANGSD/misc/realSFS results/sfs$POP.saf.idx 2> /dev/null > results/$POP.sfs`

`done`


Once again, look at the output that we have:


`cat results/sfspop1.sfs`


Here, for instance, what I have

**9478.656514 81.942662 45.773892 99.836446 1.958466 80.032460 11.068303 2.324720 58.955779 21.461398 0.000507 0.003621 40.847895 0.000000 28.613466 11.044443 0.000000 0.000000 27.445907 0.000143 10.033377**

These values are the expected number of sites with derived allele frequency equal to 0 (first number), 1 (second number), 2 (third number)... in my first population.

3. You can finally plot your results in R


`Rscript $NGSTOOLS/Scripts/plotSFS.R results/sfspop1.sfs results/sfspop2.sfs results/sfspop3.sfs`


`evince results/pop1-2-3.pdf`


[Here](https://github.com/beatriceparodi/ngsSim-simulation/blob/master/sfspop1_sfspop2_sfspop3.pdf) my plot

4. Generate replicates may help you to get confidence intervals when studing populations demography. Here how to perform a bootstrapped replicates of the SFS


`$ANGSD/misc/realSFS results/sfspop1.saf.idx -bootstrap 10  2> /dev/null > results/pop1.boots.sfs`


`cat Results/pop1.boots.sfs`


5. Also, you can estimate a multi-dimensional SFS, for example the joint SFS between two (2D) or three (3D) populations. This could help you while studying the divergence processes (for instance, if you want to track back two or more species to their common ancestor after a migration).

When running this analysis **make sure that you are comparing exactly the same sites** between populations. 
Let's start with 2D SFS between pop1 and all other populations, and after that we'll move to a 3D SFS comparison.


`for POP in pop2 pop3`

`do`

   `echo $POP`
    
    `$ANGSD/misc/realSFS results/sfs$POP.saf.idx results/sfspop1.saf.idx 2> /dev/null > results/$POP.pop1.sfs`

`done`


Between population 2 and 3:


`$ANGSD/misc/realSFS results/sfspop2.saf.idx results/sfspop3.saf.idx 2> /dev/null > results/pop2.pop3.sfs`


Between population 3 and 1:


`$ANGSD/misc/realSFS results/sfspop3.saf.idx results/sfspop1.saf.idx 2> /dev/null > results/pop3.pop1.sfs`


6. Now plot yopur results with R, and make sure that you have the rigth number of samples per population (for instance, I have 10 individuals for each population)


`Rscript $NGSTOOLS/Scripts/plot2DSFS.R Results/pop2.pop3.sfs 10 10`


Here a [link](https://github.com/beatriceparodi/ngsSim-simulation/blob/master/pop2.pop3.sfs.pdf) to my plot 

7. You may want to perform a 3D SFS. The command is:


`$ANGSD/misc/realSFS results/sfspop1.saf.idx results/sfspop2.saf.idx results/sfspop3.saf.idx 2> /dev/null > results/pop1.pop2.pop3.sfs`


## Fixation index (FST) and population branch statistics (PBS)
### Population genetic differentiation using ANGSD
FST is a measure of population substructure and is a useful tool when examining the overall genetic divergence among subpopulations, while PBS is a statistical method to detect selection.

More info about [FST](http://php.scripts.psu.edu/users/n/x/nxm2/1977%20Publications/1977-nei2.pdf) and [PCB](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3708501/#mst089-B50)

1. Compute per site FST index with the following command:


`$ANGSD/misc/realSFS fst index results/sfspop2.saf.idx results/sfspop3.saf.idx results/sfspop1.saf.idx -sfs results/pop2.pop3.sfs -sfs results/pop2.pop1.sfs -sfs results/pop3.pop1.sfs -fstout results/FSTpop1.pbs`


and look at your output with


`$ANGSD/misc/realSFS fst print results/FSTpop1.pbs.fst.idx | less -S`


Here what I have found:

**chrSIM  21   |   0.000000    |    0.000000    |    0.120467   |    0.348769**

In the first column you can see the chromosome and it position, (a) and the values for the three FST comparisons.
As you can see if you take a look at the whole output, FST value range from 0 to 1, where 0 means complete sharing of genetic traits and 1 means no sharing at all.

2. Now we are going to perform a [sliding window analysis](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/14/5/10.1093_bioinformatics_14.5.467/1/140467.pdf?Expires=1495818788&Signature=CpuZVvrrJz9CGkX9lC6l8y3OsVsqx7a84u1q2FKZKgW2uYGViPILt5ihk3M9d2DeKQvuQfsZLTN3u0FAFD5JjwHreKHezYVL-1I8f1kZ74Vt9YQAkjfbeClX7owfMR8~0eB419tiorEQnjbQvWK4Uj5lzVlkxU9VEqhwLVtYptdx~saKk-XDv4n8rPbYXBjchhWxfAatcv5J-DiHuCG7BANj4HiTfmqprUhfXcK-2AsCIqRhKJV~YQeVvpvH7360XU5LT~xV2qVgZOhuKoiIGvqcpCrHUsjkNw0I5NPUsUeH1vgb1IsdtdUHCUUbyX-9NGjtxB3EjoIXhTfsvH2Y8Q__&Key-Pair-Id=APKAIUCZBIA4LVPAVW3Q). With this method, data are plotted as moving averages of a particular criterion, such as the number of changing nucleotides, for a window of a certain length


`$ANGSD/misc/realSFS fst stats2 results/FSTpop1.pbs.fst.idx -win 50000 -step 10000 > results/FSTpop.pbs.txt`


And again take a look at the output 


`less -S results/FSTpop.pbs.txt`


You'll see the following headers:

**region  chr   midPos  Nsites  Fst01   Fst02   Fst12   PBS0    PBS1    PBS2**

PBS values are given assuming a target population (in my case population 1), and represent the differentiation between populations.


`$ANGSD/misc/realSFS fst index results/sfspop2.saf.idx results/sfspop3.saf.idx -sfs results/pop2.pop3.sfs -fstout results/POP2.POP3`


## Nucleotide diversity
### How to assess nucleotide diversity in a population
Nucleotide diversity is used to measure genetic variation, and more specifically the degree of polymorphism within a population.
It may provide valuable information when studying demographic history of populations,patterns of migrations or signatures of selection events.
[Here](http://www.genetics.org/content/166/3/1375.long) a useful paper about nucleotide diversity.

1. First, we have to compute the allele frequency posterior probabilities and associated statistics (-doThetas). SFS is used to calculate these as a proor information (-pest).


`for POP in pop1 pop2 pop3`

`do`

    `echo $POP`
   
    `$ANGSD/angsd -glf Data/$POP.glf -fai Data/ref.fasta.fai -out results/nd$POP   -nInd 10 -minInd 3 -doglf 1 -doSaf 1 -doThetas 1 -pest results/sfs$POP.sfs`

`done`


2. We need to index these files with the following commands, then perform a sliding-window analysis. The third line will allow you to index files and the following will perform the sliding-window analysis. 


`for POP in pop1 pop2 pop3`

`do`

    `echo $POP`
  
   `$ANGSD/misc/thetaStat print results/nd$POP.glf.gz`

  `$ANGSD/misc/thetaStat do_stat results/nd$POP.thetas.idx -win 5000 -step 5000 -   outnames results/$POP.thetas`


`done`


and have a look at the output (for instance I want too see population 3) with


`less -S results/pop3.thetas.pestPG`


The first colums contains information about region, the second and the third are reference name and centre of window.
Then we have 5 different Theta estimators (Watterson, pairwise, FuLi, fayH, L) and 5 neutrality test statistics (Tajima's D, Fu&Li, F's, Fu&li's D Fay's H, Zeng's).

3. You may be interested in calculate allele frequencies for single SNPs.
ANGSD will allow you to choose a subset of positions using `-sites` option. 
IN this case, you'll have to create a .txt file with your position, and then create an index with


`$ANGSD/angsd sites index Data/snps.txt`


Then, we are going to use`-doMajorMinor5` to polarise our alleles. 


`for POP in pop1 pop2 pop3`


`do`


   `echo $POP`

    `$ANGSD/angsd -glf Data/$POP.glf -fai Data/ref.fasta.fai -anc Data/  ref.fasta -out results/sites$POP -doglf 1 -nInd 10 -doMajorMinor 5 -doMaf 1 -sites Data/snps.txt` 


`done`


Finally, we can have a look at our output with


`zcat results/sitespop3.mafs.gz`

## ngsLD
### Linkage disequilibrium with ngsLD
Linkage disequilibrium (LD) refers to the nonrandom association of alleles in haplotypes. In genreal, you'll found that chromosomes sampled from non-related individuals are much more distantly related than those sampled from members of, for example, closed populations.

[Here](https://www.ncbi.nlm.nih.gov/pubmed/11818140) more info about linkage disequilibrium.
 
With ngsLD you can estimate pairwise LD taking into account the uncertainty of genotype's assignation by using genotype likelihoods or posterior probabilities.

1. We have to create a reference file with the site coordinates with 


`cat $SIM_DATA/testA.geno | perl -s -p -e 's/0 0/0/g; s/(\w) \1/2/g; s/\w \w/1/g; $n=s/2/2/g; tr/02/20/ if($n>$n_ind/2)' -- -n_ind=$N_IND | awk '{print "chrSIM\t"NR"\t"$0}' | gzip -cfn --best > testLD_T.geno.gz`


`zcat testLD_T.geno.gz | awk 'BEGIN{cnt=1} pos > 10000 {pos=0; cnt++} {pos+=int(rand()*1000+1); print $1"_"cnt"\t"pos}' > testLD.pos`


3. Then, use ANGSD and the option `-doGeno32` 


`$ANGSD/angsd -glf Data/pop.glf -fai Data/ref.fasta.fai -nInd 30 -doMajorMinor 1 -doPost 1 -doMaf 1 -doGeno 32 -out results/testLD_32`


4. Unzip the output with


`gunzip -f testLD_32.geno.gz`

Before running the analysis, let's have a look at the parameters available for ngsLD

Parameter | Usage
--------- | -------
`--geno FILE` |  input file with genotypes, genotype likelihoods or genotype posterior probabilities
`--probs`| is the input genotype probabilities (likelihoods or posteriors)?
`--log_scale`| if the input in log-scale
`--n_ind INT`| sample size (number of individuals)
`--n_sites INT`| total number of sites
`--pos` | input file with site coordinates
`--max_kb_dist DOUBLE`| maximum distance between SNPs (in Kb) to calculate LD. If set to 0 (zero) will perform all comparisons
`--max_snp_dist INT`| maximum distance between SNPs (in number of SNPs) to calculate LD. If set to 0 (zero) will perform all comparisons
`--min_maf DOUBLE` | minimum SNP minor allele frequency
`--call_geno` | call genotypes before running analyses
`--N_thresh DOUBLE` |  minimum threshold to consider site; missing data if otherwise (assumes -call_geno)
`--call_thresh DOUBLE` |  minimum threshold to call genotype; left as is if otherwise (assumes -call_geno)
`--rnd_sample DOUBLE` | proportion of comparisons to randomly sample
`--seed INT` | random number generator seed for random sampling (--rnd_sample)
`--out FILE` | output file name
`--n_threads INT` |  number of threads to use
`--version` |  prints program version and exits
`--verbose INT` | selects verbosity level


5. Perform the analysis with ngsLD. The option `-prob` will allow you to take genotype likelihoods or posterior probabilities into account

`$NGSLD --geno results/testLD_32.geno --n_ind 30 --n_sites 10000 --probs --pos testLD.pos --out results/LDpop`

6. Now we can have a look at our output with 

`cat LDpop`

You'll see results for all pairs of sites for which LD was calculated as: site1 label, site2 label, distance between sites (bp), r^2 from pearson correlation between expected genotypes, D from EM algorithm, D' from EM algorithm, and r^2 from EM algorithm.
