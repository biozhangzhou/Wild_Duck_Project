#!/bin/bash
#180 sample SNP-calling#
NA
#vcf QC
#maf0.01 geno0.2 and to plink #map ped to binary
for i in `cat SNP_VCF_List`
do
vcftools --vcf $i --remove-filtered-all  --maf 0.01 --max-missing 0.8 --plink --out tmp
ID=`echo $i|sed 's/.vcf//g'`
~/plink2/plink --file tmp --make-bed --out $ID\.maf0.01.geno0.2 --chr-set 40;
rm tmp.map tmp.ped
done &

#Merge#
~/plink2/plink --bfile chr01.100-110.SNP_filter.maf0.01.geno0.2 --merge-list Final-Autosome/Autosome-List --make-bed --out Final-Autosome/180 --chr-set 40 
~/plink2/plink --bfile ../0.Plink/chrZ-1-20.SNP_filter.maf0.01.geno0.2 --merge-list Z-List --make-bed --out 180-Z --chr-set 40   &

#Female Z
nohup ./get_Z_haplotype.sh &
~/plink2/plink --bfile 180-Z --remove tmp --out 180-Z-Male --make-bed --chr-set 40  &
~/plink2/plink --bfile 180-Z-Male --merge-list Female-List --maf 0.01 --geno 0.2 --chr 0 --from-bp 2000000 --to-bp 86000000 --make-bed --out 180-Hapolotype-NonPAR --chr-set 40 &



#Derived SNP
nohup ./get-derived-SNP.sh &
#

#SNP calling of 20 low sequencing depth samples.
SNP-Lists are from QC-maf0.01-geno0.2 Final-Autosome and Final-Z
1.HaplotypeCaller
#Autosom
~/gatk-4.1.8.0/gatk HaplotypeCaller -R CAU-Pek-ZJU-W.fa -I ${ID}_mkdup.bam -L Autosome-SNP-List.intervals -O ${ID}-autosome.vcf --emit-ref-confidence BP_RESOLUTION &
#Z
~/gatk-4.1.8.0/gatk HaplotypeCaller -R CAU-Pek-ZJU-W.fa -I ${ID}_mkdup.bam -L Z-SNP-List.intervals -O ${ID}-Z.vcf --emit-ref-confidence BP_RESOLUTION 
2.CombineGVCFs
3.GenotypeGVCFs
nohup ~/gatk-4.1.8.0/gatk --java-options -Xmx200G GenotypeGVCFs  -R CAU-Pek-ZJU-W.fa  --include-non-variant-sites --variant chr0${i}-Ma-autosome.vcf  -O 2chr0${i}-Ma-autosome.vcf &
cp 2chr01-Ma-autosome.vcf Autosome.vcf
for i in `seq 2 9`; do grep -v ^\# 2chr0${i}-Ma-autosome.vcf >> Autosome.vcf ; done &
for i in `seq 10 40`; do grep -v ^\# 2chr${i}-Ma-autosome.vcf >> Autosome.vcf ; done &

4.QC gatk SNP hardfilter and minDP 5
/home/zhangzhou209/1.bin/gatk-4.1.8.0/gatk --java-options -Xmx40G VariantFiltration -V 2Ma-Z.vcf --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Filter" -O 2Ma-Z.SNP_filter.vcf
vcftools --vcf Ma-Z.SNP_filter.vcf --remove-filtered-all --min-meanDP 5 --max-missing 1 --plink --out Ma-ChrZ &
vcftools --vcf Ma-Autosome.SNP_filter.vcf --remove-filtered-all --min-meanDP 5 --max-missing 1 --plink --out Ma-Autosome &

5.Combine
#Z
nohup ./get_Z_haplotype.sh &
~/plink2/plink --file ~/tmp-2022/Ma-Z/Ma-ChrZ --keep 20-Male --make-bed --out Male --chr-set 40 &
~/plink2/plink --bfile Male --merge-list Homo-list --chr 0 --from-bp 2000000 -to-bp 90000000 --make-bed --out 20-NONPAR --chr-set 40 
~/plink2/plink --bfile ../180-Hapolotype-NonPAR --bmerge 20-NONPAR.bed 20-NONPAR.bim 20-NONPAR.fam --chr 0 --from-bp 2000000 -to-bp 90000000 --make-bed --out 200-maf0.05-geno0.1 --maf 0.01 --geno 0.1 --chr-set 40  
~/plink2/plink --file 20-NONPAR --exclude 200-maf0.05-geno0.1-merge.missnp --make-bed --out 20-NONPAR_exclude

~/plink2/plink --bfile ../180-Hapolotype-NonPAR --bmerge ../20-Ma/20-NONPAR_exclude.bed ../20-Ma/20-NONPAR_exclude.bim ../20-Ma/20-NONPAR_exclude.fam --chr 0 --from-bp 2000000 -to-bp 90000000 --make-bed --out 200-maf0.05-geno0.1 --maf 0.01 --geno 0.1 --chr-set 40

#Autosome
~/plink2/plink --bfile 180-Autosome --merge 20-Ma-Autosome.ped 20-Ma-Autosome.map --maf 0.01 --geno 0.1 --make-bed --out 200-Autosome --chr-set 40
~/plink2/plink --bfile 20-Ma-Autosome --exclude 200-Autosome-merge.missnp --make-bed --out Ma-Autosome_exclude --chr-set 40
~/plink2/plink --bfile 180-Autosome --bmerge Ma-Autosome_exclude.bed Ma-Autosome_exclude.bim Ma-Autosome_exclude.fam  --maf 0.01 --geno 0.1 --make-bed --out 200-Autosome --chr-set 40 &
~/plink2/plink --bfile 200-Autosome --indep-pairwise 50 10 0.5 --out 200-Autosome-maf0.01-geno0.1-LD-50-10-0.5 --chr-set 40
~/plink2/plink --bfile 200-Autosome --extract 210-Autosome-maf0.01-geno0.1-LD-50-10-0.5.prune.in --make-bed --out 200-Autosome-maf0.01-geno0.1-LD-50-10-0.5 --chr-set 40 










