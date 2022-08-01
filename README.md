# Oryza_Satvia_Snakemake
```r

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

#Test your workflow
snakemake -npr

```

![Screenshot from 2022-08-01 13-43-14](https://user-images.githubusercontent.com/93121277/182141099-3ec98d0e-cdb1-409e-ba9d-a2b4c0352e55.png)


![dag](https://user-images.githubusercontent.com/93121277/182140684-39e6ba3a-d5cb-4a9d-8021-2a09ccce43f3.svg)

```r
#Run it if everything looks clean!
snakemake -j 4 -kpr
```


