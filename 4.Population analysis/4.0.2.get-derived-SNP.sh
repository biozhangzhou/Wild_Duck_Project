#!/bin/bash

while read line
do
echo $line
ls_date=`date +m%d%H%M%S`
ID=`echo $line|sed 's/.vcf//g'`
vcftools --vcf $line --keep Goose --remove-filtered-all --maf 0 --max-missing 1 --plink --out goose-tmp
awk '{print $2}' goose-tmp.map|sed 's/:/\t/g' > snp-list
#vcftools --vcf $line --keep Duck-Out --remove-filtered-all --max-maf 0 --max-missing 1 --plink --out duck-tmp
#awk '{print $2}' duck-tmp.map|sed 's/:/\t/g' > snp-list
#vcftools --vcf $line --positions snp-list --keep Derived/Plink-Pop-File --mac 1 --remove-filtered-all --max-missing 1 --plink --out Derived/Derived-$ID 
vcftools --vcf $line --positions snp-list --mac 1 --remove-filtered-all --max-missing 1 --plink --out Derived-Goose/Derived-$ID
done<SNP_VCF_List


#此处的goose derived相当取出鹅中没有缺失的位点 而后又取出185个个体中没有缺失的位点，过于严格了。


vcftools --vcf ../192-duck-tmp.vcf --keep Goose --remove-filtered-all --max-maf 0 --max-missing 1 --plink --out Goose-tmp
awk '{print $2}' Goose-tmp.map|sed 's/:/\t/g' > snp-list

vcftools --vcf ../192-duck-tmp.vcf --positions snp-list --remove-filtered-all --plink --out Derived-Goose
awk '{print $2}' Goose-tmp.map|sed 's/:/\t/g' > snp-list

vcftools --vcf ../192-duck-tmp.vcf --positions snp-list --remove-filtered-all --plink --out Derived-Goose
/opt/software/plink2/plink --file Derived-Goose --recode vcf -out Goose-Derived --chr-set 40 &



/opt/software/plink2/plink --bfile ../Final-Autosome/180-Autosome --recode vcf -out 180-duck-tmp --chr-set 40
vcftools --vcf 180-duck-tmp.vcf --keep Goose --remove-filtered-all --max-maf 0 --max-missing 1 --plink --out Goose-tmp &
vcftools --vcf ../Derived-Goose/180-duck-tmp.vcf --keep Duck --remove-filtered-all --max-maf 0 --max-missing 1 --plink --out Duck-tmp &
awk '{print $1"\t"$4}' Duck-tmp.map > snp-list &
vcftools --vcf ../Derived-Goose/180-duck-tmp.vcf --positions snp-list --remove-filtered-all --plink --out Derived-Duck &
/opt/software/plink2/plink --file Derived-Duck --recode vcf -out Derived-Duck --chr-set 40

#Z
vcftools --vcf ../Final-Z/180-Hapolotype-NonPAR.vcf --keep Duck2 --remove-filtered-all --max-maf 0 --max-missing 1 --plink --out Duck-tmp
awk '{print $1"\t"$4}' Duck-tmp.map > snp-list &
vcftools --vcf ../Final-Z/180-Hapolotype-NonPAR.vcf --positions snp-list --remove-filtered-all --plink --out Derived-Duck-Z &
/opt/software/plink2/plink --file Derived-Duck-Z --recode vcf -out Derived-Duck-Z











zcat Derived-Goose-Z.recode.vcf.gz|grep ^\# > TMP3
paste -d '\t' <(zcat Derived-Goose-Z.recode.vcf.gz|cut -f 1-7|grep -v ^##) <(awk '{if(NR==1)print "INFO";else print $4}' freq_stat.frq) <(zcat Derived-Goose-Z.recode.vcf.gz|cut -f 9-|grep -v ^##) | awk '{FS="\t";OFS="\t"}{if($4==$8){$8=".";print $0}else if($5==$8){$5=$4;$4=$8;$8=".";for(i=10;i<=NF;i++){gsub("1","a",$i);gsub("0","1",$i);gsub("a","0",$i)};print $0}}' >> TMP3

#检查
paste freq_stat.frq <(cat Derived-Goose-Duck-Z.recode-ancestral.vcf|grep -v ^##)|awk '{print $4,$10}'|awk '{if($1!=$2)print $0}'|wc


















