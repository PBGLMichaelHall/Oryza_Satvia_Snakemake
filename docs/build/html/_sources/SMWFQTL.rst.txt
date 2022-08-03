================================================================================
A PBGL Snakemake Workflow and Quatitative Trait Locus on Bulk Segregant Analysis 
================================================================================

=====================
:Author: Michael Hall
:Date: 07/18/2022
====================



Software Prerequisites
======================

#sra-toolkit
#Download git repository:

.. code:: shell

   #Clone pbgl DNA Proto Workflow
	git clone https://github.com/pbgl/dna-proto-workflow.git

	#Change Directory to root or project
	cd dna-proto-workflow

	#Create virtual environment
	Mamba env create --file env/all-dependencies.yml

	#Activate the environment
	conda activate dna-proto

	#Make a new directory FASTQ
	mkdir FASTQ

	#Change directory to FASTQ
	cd FASTQ

	#Download data from NCBI
	#Start with Extremely Tolerant Rice Pool ("Mutant")
	wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR834931/SRR834931

	#Use sra-toolkit to split SRA into foward and reverse FASTQ reads and compress the files
	fastq-dump --gzip --split-3

	#Rename foward and reverse read files
	ET_385_1.fq.gz
	ET_385_2.fq.gz

	#Now get Extremely Sensitive Rice Pool ("Wild-Type")
	wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR834927/SRR834927

	#Use sra-toolkit to split SRA into foward and reverse FASTQ reads and compress the files
	fastq-dump --gzip --split-3

	#Rename foward and reverse read files
	ES_430_1.fq.gz
	ES_430_2.fq.gz

	#Follow technical documentation for the DNA Prototype Workflow
	https://dna-proto-workflow-master.readthedocs.io/en/latest/

	#Download Reference Genome and Index it
	#Download annotation files (gff, gtf, protein)
	#Build a snpEFF configuration file with soft links if necessary or preferred.
	#Provide Meta Data information (Contigs of interest, sample sets and definitions)
	#Configure config.yml and snpEff.config to match file pathways and names.\

	#Test your workflow on a dry-run
	snakemake -npr



.. image:: ../images/screen1.png

.. code:: shell

	snakemake --dag -npr -j 1 | dot -Tsvg > dag.svg

.. image:: ../images/dag.svg
   :alt: Directed Acrylic Graph
   :scale: 250 %

.. _a link: https://raw.githubusercontent.com/PBGLMichaelHall/Oryza_Satvia_Snakemake/main/dag.svg

.. code:: shell

	snakemake -j 4 -kpr 

	#After PBGL Dna prototype pipeline root directory contains 634 items amounting to 271.4 Gigabytes 
	#There should be a vcf file in output variants final, lets take a look!

	bcftools view /home/michael/dna-proto-workflow/output/variants/final/freebayes~bwa~GCF_001433935.1_IRGSP-1.0~all_samples~filtered-strict.vcf.gz | less -S
   
	#Filter variants for biallelic sites only and SNPS types

	bcftools view -m2 -M2 -v snps -o freebayes~bwa~GCF_001433935.1_IRGSP-1.0~all_samples~filtered-strict.Biallelic.vcf.gz freebayes~bwa~GCF_001433935.1_IRGSP-1.0~all_samples~filtered-strict.vcf.gz

	#Follow pbgl online documentation

	https://github.com/pbgl/QTLseqr


	conda deactivate dna-proto


	#Create and activate new R Environment

	mamba env create --file env/R.yaml

	
	#Set samples
	HighBulk <- "ET_385" # SRA-run: SRR834931
	LowBulk <- "ES_430" # SRA-run: SRR834927 

	#set file name of the VCF file to load
	file <- "freebayes~bwa~GCF_001433935.1_IRGSP-1.0~all_samples~filtered-strict~snpEff.Biallelic.vcf.gz"

	#Specify which chromosomes should be included in the analysis (i.e., exclude smaller contigs)
	Chroms <- c("NC_029256.1",
        "NC_029257.1",
        "NC_029258.1",
        "NC_029259.1",
        "NC_029260.1",
        "NC_029261.1",
        "NC_029262.1",
        "NC_029263.1",
        "NC_029264.1",
        "NC_029265.1",
        "NC_029266.1",
        "NC_029267.1")
	df <- 
	importFromVCF(
	file = file,
	highBulk = HighBulk,
	lowBulk = LowBulk,
	chromList = Chroms
	)


	df_filt <-
	filterSNPs(
	SNPset = df,
	refAlleleFreq = 0.20,
	minTotalDepth = 70,
	maxTotalDepth = 200,
        minSampleDepth = 30,
        depthDifference = 100,
        #minGQ = 150,
        verbose = TRUE
        )

	df_filt <- 
        runGprimeAnalysis(
        SNPset = df_filt,
        windowSize = 1e6,
        outlierFilter = "deltaSNP",
        filterThreshold = 0.1
        )

	df_filt <- 
        runQTLseqAnalysis(
        SNPset = df_filt,
        windowSize = 1e6,
        popStruc = "F2",
        bulkSize = c(385, 430), 
        replications = 10000,
        intervals = c(95, 99)
        )

	plotQTLStats(SNPset = df_filt, var = "nSNPs", plotIntervals = TRUE)


.. figure:: ../images/nsnps.png
	
.. code:: shell

	plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)


.. figure:: ../images/delta.png
	
.. code:: shell

	plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.02)


.. figure:: ../images/Gprime.png

.. code:: shell

	df_filt$p1 <- df_filt$AD_ALT.LOW/df_filt$DP.LOW
	df_filt$p2 <- df_filt$AD_ALT.HIGH/df_filt$DP.HIGH
	df_filt3 <- df_filt %>% dplyr::filter(CHROM=="NC_029263.1")
	#Perhaps my filter is biased, but here it is. I want observed variants with pvalues below a significant threshold.
	df_filt4 <- df_filt3 %>% dplyr::filter(pvalue < .0000962)
	write.csv(df_filt3,file = "QTL.csv",sep=",")


.. image:: ../images/pvalue.png

.. csv-table:: Significant Variants
   :url: https://github.com/PBGLMichaelHall/Oryza_Satvia_Snakemake/blob/main/QTL.csv


	#Now, we open up Integrative Genomic Viewer to further analyze the most significant variant called to validate information.

.. code:: shell

	bash IGV.sh

.. image:: ../images/IGV.png
   :target: ../images/IGV.png
   :alt: Integrated Genomic Viewer
