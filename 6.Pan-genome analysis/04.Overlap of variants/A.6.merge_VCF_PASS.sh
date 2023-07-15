### Pre-processing VCF file
# Select SVs in reference genome
# Fix the VCF header
# compress and construct index
for i in `ls A.6.graph_??_A_a.vcf`
do
    ### filter and vcf2bed
    # select VCF chromosome
    bcftools view -i 'CHROM ~ "chr[0123456789][0123456789]"' -Ov -o cache1_1.vcf $i
    bcftools view -i 'CHROM ~ "chrZ"' -Ov -o cache1_2.vcf $i
    grep -v "#" cache1_2.vcf > cache1_2_grep.vcf
    cat cache1_1.vcf cache1_2_grep.vcf > $i.filtered

    # reheader & PASS
    grep -v "##contig=<ID=Chr" $i.filtered |  grep "##contig=<ID=chr" | sort > sort_contig.txt
    grep -v "contig" $i.filtered | sed '1r sort_contig.txt'| bcftools sort -m 10G -Ov | bcftools view -i 'FILTER="PASS"' | bcftools sort > $i.filtered.reheader.PASS

    # bgzip and tabix
    bgzip -f -k $i.filtered.reheader.PASS
    tabix -p vcf $i.filtered.reheader.PASS.gz

    # delete cache file
    rm $i.filtered
    rm $i.filtered.reheader.PASS
done

### norm, del 0/0, select PASS, SVLEN>=50
# merge {sample}.vcf to population.vcf
ls A.6.graph_??_A_a.vcf.filtered.reheader.PASS.gz > A.6.sample_A_a_list.txt
bcftools merge --threads 10 --force-samples -m all --file-list A.6.sample_A_a_list.txt -Ov | bcftools sort -m 10G -Ov > A.6.population_PASS_A_a_sort.vcf
bcftools norm -m -both -f ../genome/1.A.3.pan_duck_genome/BJ_chr.fa A.6.population_PASS_A_a_sort.vcf | bcftools sort > A.6.population_PASS_A_a_split_sort.vcf
bcftools view -i 'GT!="0/0"' A.6.population_PASS_A_a_sort.vcf | bcftools sort > A.6.population_PASS_A_a_v00_sort.vcf
bcftools norm -m -both -f ../genome/1.A.3.pan_duck_genome/BJ_chr.fa A.6.population_PASS_A_a_v00_sort.vcf | bcftools sort > A.6.population_PASS_A_a_v00_split_sort.vcf

### remove cache
rm A.6.graph_??_A_a.vcf.filtered.reheader.PASS.gz

