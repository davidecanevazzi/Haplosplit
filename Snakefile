###############################################################################

################################# Haplosplit ##################################

###############################################################################

# import modules
import numpy as np

import pandas as pd

import os

import sys

import subprocess


def substring_after(s, delim):
    return s.partition(delim)[2]


data = open("data.txt")

fastq = data.readline()
ref = data.readline()
gene = data.readline()
pos = data.readline()
p_mod = data.readline()
m_json = data.readline()

fastq=substring_after(fastq, ",")
ref=substring_after(ref, ",")
gene=substring_after(gene, ",")
pos=substring_after(pos, ",")
p_mod=substring_after(p_mod, ",")
m_json=substring_after(m_json, ",")


chr=pos.partition(":")[0]


fastq=fastq[:-1]
ref=ref[:-1]
gene=gene[:-1]
pos=pos[:-1]
p_mod=p_mod[:-1]
m_json=m_json[:-1]

configfile:
    "config.json"

rule all:
    input:
        "data.txt"

rule minimap:
    output:
        "data/{sample}.sam"
    shell:
        "minimap2 -a -t 24 -z 600,200 -x map-ont {ref} {fastq} > {output}"



rule samtools_view:
    input:
        s="data/{sample}.sam"
    output:
        "data/{sample}.bam"
    shell:
        "samtools view -@12 -bS {input.s} > {output}"

rule samtools_sort:
    input:
        "data/{sample}.bam"
    output:
        "out/{sample}.sorted.bam"
    shell:
        "santools sort -@ 12 {input} -o {output}"

rule samtools_index:
    input:
        bam="out/{sample}.sorted.bam"
    output:
        "log/{sample}.index.txt"
    shell:
        "samtools index {input.bam}|tee {output}"

rule gene_selection:
    input:
        bam= "out/{sample}.sorted.bam",
        bai= "out/{sample}.sorted.bam.bai",
        log="log/{sample}.index.txt"
    output:
        "out/{sample}_region.bam"
    shell:
        "samtools view -bS {input.bam} {pos} > {output}"

rule reg_index:
    input:
        "out/{sample}_region.bam"
    output:
        "out/reg_{sample}.index.txt"
    shell:
        "samtools index {input}|tee {output}"

rule pepper:
    input:
        bam="out/{sample}_region.bam" ,
        ind="out/reg_{sample}.index.txt"
    output:
        "out/logs/{sample}_pepper_snp.log"
    shell:
        "singularity exec /apps/PEPPER/0.4/pepper pepper_snp  call_variant -b  {input.bam} -f  {ref} -t 4 -m {p_mod} -o out -r {chr} -s Sample -w 4 -bs 64 --ont 2>&1|tee {output}"

rule bgzip:
    input:
        "out/logs/{sample}_pepper_snp.log",
    output:
        "out/pepper_snp/{sample}_log_gzip.txt"
    shell:
        "bgzip out/PEPPER_SNP_OUTPUT.vcf|tee {output}"

rule mv:
    input:
        "out/pepper_snp/{sample}_log_gzip.txt"
    output:
        "out/{sample}_PEPPER_SNP_OUTPUT.vcf.gz"
    shell:
        "mv out/PEPPER_SNP_OUTPUT.vcf.gz {output}"

rule tabix:
    input:
        "out/{sample}_PEPPER_SNP_OUTPUT.vcf.gz"
    output:
        "out/{sample}_PEPPER_SNP_OUTPUT.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}"


rule margin:
    input:
        bam="out/{sample}_region.bam",
        vcf="out/{sample}_PEPPER_SNP_OUTPUT.vcf.gz",
        vcfi="out/{sample}_PEPPER_SNP_OUTPUT.vcf.gz.tbi",
    output:
        "out/logs/{sample}_margin_haplotag.log"
    shell:
        "time margin phase {input.bam}  {ref}  {input.vcf} {m_json} -t 4 -r {chr} -V -o out/MARGIN_PHASED.PEPPER_SNP_MARGIN 2>&1 | tee {output}"

rule samtools_index2:
    input:
        "out/logs/{sample}_margin_haplotag.log"
    output:
        "out/logs/{sample}_index2.log"
    shell:
        "samtools index out/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam| tee {output}"

rule split:
    input:
        ind="out/logs/{sample}_index2.log",
        bam="out/logs/{sample}_margin_haplotag.log"
    output:
        "out/logs/{sample}_split.log"
    shell:
        "bamtools split -in out/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam -tag HP|tee {output}"

rule to_fasta_A1:
    input:
        "out/logs/{sample}_split.log"
    output:
        "out/{sample}_region_hap1.fa"
    shell:
        "samtools fasta out/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.TAG_HP_1.bam -0 {output}"


rule to_fasta_A2:
    input:
        "out/logs/{sample}_split.log"
    output:
        "out/{sample}_region_hap2.fa"
    shell:
        "samtools fasta out/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.TAG_HP_2.bam -0 {output}"

rule flye_A_1:
    input:
        "out/{sample}_region_hap1.fa"
    output:
        "out/flye_hap1/{sample}_assembly.log"
    shell:
        "flye --nano-raw {input} --out-dir out/flye_hap1/ --threads 4 -m 1000 -i 2|tee {output}"

rule flye_A_2:
    input:
        "out/{sample}_region_hap2.fa"
    output:
        "out/flye_hap2/{sample}_assembly.log"
    shell:
        "flye --nano-raw {input} --out-dir out/flye_hap2/ --threads 4 -m 1000 -i 2|tee {output}"
