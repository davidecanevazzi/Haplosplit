# Haplosplit
Haplotyping of long Oxford Nanopore reads

---
## Setup on the cluster:

1. Create a new conda environment with the command:  
```
conda env create -f config.json 
```
2. Then load some necessary modules:  module load PEPPER MARGIN
```
module load PEPPER MARGIN
```


---
## How to run

1. Change the configuration file data.txt, it comprehends 6  fields:
    - fastq_file --> path to the fastq file
    - reference --> path to the reference.fasta (**already indexed**)
    - name_of_the_gene --> name of the output region/gene
    - position --> coordinates of the region of interest
    - path_to_pepper_model --> path to the pepper model PEPPER_SNP_R941_ONT_V4.pkl
    - Margin_json --> path to Margin parameters

2. If you want to run all the pipeline until the haplotypes: 


```
snakemake --cores 24 out/flye_hap1/{sample}_assembly.log
snakemake --cores 24 out/flye_hap2/{sample}_assembly.log
```


3. If you want to run the pipeline step by step just run until the correct output
