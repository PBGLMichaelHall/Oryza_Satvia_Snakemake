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

	
