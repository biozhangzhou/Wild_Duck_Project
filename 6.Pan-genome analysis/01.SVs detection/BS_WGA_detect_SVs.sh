#!/bin/bash

###################################
### ENV			          #
#conda activate SyRI              #
### Software                      #
# `minimap2` == 2.24-r1122        #
# `paftools.js` == 2.24-r1122     #
# `vcflib` == 1.0.0_rc2           #
# `jasmine` == 1.1.4              #
# `bcftools` == 1.9               # 
###################################

:<<!
### WGA SVs detect 
# Mapping and detect SVs
# Select SVs length >= 50
for sample in BJ SX Ma SB ZW LC LW PZ HL
do       
    ##### asm5
    ### Mapping & call SVs
    minimap2  -cx asm5 -t 18 \
    --cs ../genome/1.A.3.pan_duck_genome/BJ_chr.fa \
    ../genome/1.A.3.pan_duck_genome/${sample}_chr.fa | sort -k6,6 -k8,8n | paftools.js call -f ../genome/1.A.3.pan_duck_genome/BJ_chr.fa - > BJ_chr_align_${sample}_asm5_sort.vcf;
    ### breakmultiallele, select SVs length >= 50
    vcfbreakmulti BJ_chr_align_${sample}_asm5_sort.vcf | vcflength | bcftools sort | bcftools view -i "length >= 50 | length <= -50 "|bcftools sort > BJ_chr_align_${sample}_asm5_sort_50.vcf
    
    ##### asm10
    ### Mapping & call SVs
    minimap2  -cx asm10 -t 18 \
    --cs ../genome/1.A.3.pan_duck_genome/BJ_chr.fa \
    ../genome/1.A.3.pan_duck_genome/${sample}_chr.fa | sort -k6,6 -k8,8n | paftools.js call -f ../genome/1.A.3.pan_duck_genome/BJ_chr.fa - > BJ_chr_align_${sample}_asm10_sort.vcf;
    ### breakmultiallele, select SVs length >= 50    
    vcfbreakmulti BJ_chr_align_${sample}_asm10_sort.vcf | vcflength | bcftools sort | bcftools view -i "length >= 50 | length <= -50 "|bcftools sort > BJ_chr_align_${sample}_asm10_sort_50.vcf

    ##### asm20
    ### Mapping & call SVs
    minimap2  -cx asm20 -t 18 \
    --cs ../genome/1.A.3.pan_duck_genome/BJ_chr.fa \
    ../genome/1.A.3.pan_duck_genome/${sample}_chr.fa | sort -k6,6 -k8,8n | paftools.js call -f ../genome/1.A.3.pan_duck_genome/BJ_chr.fa - > BJ_chr_align_${sample}_asm20_sort.vcf;
    ### breakmultiallele, select SVs length >= 50        
    vcfbreakmulti BJ_chr_align_${sample}_asm20_sort.vcf | vcflength | bcftools sort | bcftools view -i "length >= 50 | length <= -50 "|bcftools sort > BJ_chr_align_${sample}_asm20_sort_50.vcf
done
!

### merge population-level SVs VCF
# merge
ls BJ_chr_align_??_asm20_sort_50.vcf > sample_list.txt
jasmine \
file_list=sample_list.txt out_file=WGA_asm20_50.vcf \
max_dist=1000 \
min_dist=50 \
min_overlap=0.8 \
min_support=1 \
threads=10 \
spec_len=50 \
spec_reads=10 \
--nonlinear_dist
# Sort
bcftools sort -m 10G -Ov WGA_asm20_50.vcf -o WGA_asm20_50_sort.vcf

### rm cache
rm BJ_chr_align_??_asm5_sort.vcf
rm BJ_chr_align_??_asm10_sort.vcf
rm BJ_chr_align_??_asm20_sort.vcf

