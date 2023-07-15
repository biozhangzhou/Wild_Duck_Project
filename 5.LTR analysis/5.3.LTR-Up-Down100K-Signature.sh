#observed variation
for j in Pek IND Ma SB
do
	awk '{if($2>-5000 && $2<=5000)print $0}' CAU-cLTR-0-100k.bed-${j}.map.bed-Intersect4|wc|awk '{print "'"$j"'","-5-5",$1}' >> True-Variant-Number
	awk '{if($2>-10000 && $2<=10000)print $0}' CAU-cLTR-0-100k.bed-${j}.map.bed-Intersect4|wc|awk '{print "'"$j"'","-10-10",$1}' >> True-Variant-Number
	awk '{if($2>-100000 && $2<=100000)print $0}' CAU-cLTR-0-100k.bed-${j}.map.bed-Intersect4|wc|awk '{print "'"$j"'","-100-100",$1}' >> True-Variant-Number
done

#permutation
for j in Pek IND Ma SB
do
{
for i in `seq 1 1000`
do
bedtools shuffle -i <(sort -k1,1 -k2,2n ../../../CAU-cLTR.bed|grep -v chrZ|grep -v HiC) -g <(grep -v HiC ../../../../1.With-W/CAU-length|grep -v ZJU-W|grep -v Z) > ${j}-SHUFFLE-CAU-cLTR.bed
bedtools flank -i ${j}-SHUFFLE-CAU-cLTR.bed -g <(grep -v HiC ../../../../1.With-W/CAU-length|grep -v ZJU-W|grep -v Z) -b 100000 > ${j}-SHUFFLE-CAU-cLTR-0-100k.bed
bedtools flank -i ${j}-SHUFFLE-CAU-cLTR.bed -g <(grep -v HiC ../../../../1.With-W/CAU-length|grep -v ZJU-W|grep -v Z) -b 10000 > ${j}-SHUFFLE-CAU-cLTR-0-10k.bed
bedtools flank -i ${j}-SHUFFLE-CAU-cLTR.bed -g <(grep -v HiC ../../../../1.With-W/CAU-length|grep -v ZJU-W|grep -v Z) -b 5000 > ${j}-SHUFFLE-CAU-cLTR-0-5k.bed
bedtools intersect -a ${j}-SHUFFLE-CAU-cLTR-0-100k.bed -b  ~/2.analysis/10.Resequencing/2.185/Clean/QC-Maf0.01-geno0.2/FST/180-Final/Ancient/${j}.map.bed -wo|wc|awk '{print "'"$j"'","-100-100",$1}' >>  Permuation1000
bedtools intersect -a ${j}-SHUFFLE-CAU-cLTR-0-10k.bed -b  ~/2.analysis/10.Resequencing/2.185/Clean/QC-Maf0.01-geno0.2/FST/180-Final/Ancient/${j}.map.bed -wo|wc|awk '{print "'"$j"'","-10-10",$1}' >>  Permuation1000
bedtools intersect -a ${j}-SHUFFLE-CAU-cLTR-0-5k.bed -b  ~/2.analysis/10.Resequencing/2.185/Clean/QC-Maf0.01-geno0.2/FST/180-Final/Ancient/${j}.map.bed -wo|wc|awk '{print "'"$j"'","-5-5",$1}' >>  Permuation1000
done
} &
done

#observed Fst
for j in Other Pekin-SX Pekin-Sp chr
do

	bedtools flank -i <(grep $j 2022-Final-Three-Class|grep -v chrZ) -g <(grep -v HiC ../../../../../1.With-W/CAU-length|grep -v ZJU-W|grep -v Z) -b 100000 > ${j}-CAU-cLTR-0-100k.bed
	bedtools flank -i <(grep $j 2022-Final-Three-Class|grep -v chrZ) -g <(grep -v HiC ../../../../../1.With-W/CAU-length|grep -v ZJU-W|grep -v Z) -b 10000 > ${j}-CAU-cLTR-0-10k.bed
	bedtools flank -i <(grep $j 2022-Final-Three-Class|grep -v chrZ) -g <(grep -v HiC ../../../../../1.With-W/CAU-length|grep -v ZJU-W|grep -v Z) -b 5000 > ${j}-CAU-cLTR-0-5k.bed
	bedtools intersect -a ${j}-CAU-cLTR-0-100k.bed -b 2Ma-Pek_IND-Fst.weir.fst -wo|awk '{sum+=$(NF-1)}END{print "'"$j"'","Ma-Home","-100-100",sum,NR,sum/NR}' >> True-Fst
	bedtools intersect -a ${j}-CAU-cLTR-0-10k.bed -b 2Ma-Pek_IND-Fst.weir.fst -wo|awk '{sum+=$(NF-1)}END{print "'"$j"'","Ma-Home","-10-10",sum,NR,sum/NR}' >> True-Fst
	bedtools intersect -a ${j}-CAU-cLTR-0-5k.bed -b 2Ma-Pek_IND-Fst.weir.fst -wo|awk '{sum+=$(NF-1)}END{print "'"$j"'","Ma-Home","-5-5",sum,NR,sum/NR}' >> True-Fst
	bedtools intersect -a ${j}-CAU-cLTR-0-100k.bed -b 2Pek-IND-Fst.weir.fst -wo|awk '{sum+=$(NF-1)}END{print "'"$j"'","Pek-IND","-100-100",sum,NR,sum/NR}' >> True-Fst
	bedtools intersect -a ${j}-CAU-cLTR-0-10k.bed -b 2Pek-IND-Fst.weir.fst -wo|awk '{sum+=$(NF-1)}END{print "'"$j"'","Pek-IND","-10-10",sum,NR,sum/NR}' >> True-Fst
	bedtools intersect -a ${j}-CAU-cLTR-0-5k.bed -b 2Pek-IND-Fst.weir.fst -wo|awk '{sum+=$(NF-1)}END{print "'"$j"'","Pek-IND","-5-5",sum,NR,sum/NR}' >> True-Fst

	bedtools intersect -a ${j}-CAU-cLTR-0-100k.bed -b 2Ma-Pek_IND-Fst.weir.fst -wo|awk '{print "'"$j"'","Ma-Home -100-100",$0}' >> True-Fst-ALL
	bedtools intersect -a ${j}-CAU-cLTR-0-10k.bed -b 2Ma-Pek_IND-Fst.weir.fst -wo|awk '{print "'"$j"'","Ma-Home -10-10",$0}'  >> True-Fst-ALL
	bedtools intersect -a ${j}-CAU-cLTR-0-5k.bed -b 2Ma-Pek_IND-Fst.weir.fst -wo|awk '{print "'"$j"'","Ma-Home -5-5",$0}'  >> True-Fst-ALL
	bedtools intersect -a ${j}-CAU-cLTR-0-100k.bed -b 2Pek-IND-Fst.weir.fst -wo|awk '{print "'"$j"'","Pek-IND -100-100",$0}'  >> True-Fst-ALL
	bedtools intersect -a ${j}-CAU-cLTR-0-10k.bed -b 2Pek-IND-Fst.weir.fst -wo|awk '{print "'"$j"'","Pek-IND -10-10",$0}'  >> True-Fst-ALL
	bedtools intersect -a ${j}-CAU-cLTR-0-5k.bed -b 2Pek-IND-Fst.weir.fst -wo|awk '{print "'"$j"'","Pek-IND -5-5",$0}'  >> True-Fst-ALL
done

awk -v ID1=0 -v ID2=0 -v ID3=0 -v ID4=0 -v count=0 '{if($1==ID1 && $2==ID2 && $3==ID3 && $7==ID4){sum+=$(NF-1);count+=1}else{if(NR>1)print ID1,ID2,ID3,ID4,sum,count,sum/count;ID1=$1;ID2=$2;ID3=$3;ID4=$7;count=1;sum=$(NF-1)}}' True-Fst-ALL|sed 's/=/ /g' > True-Fst-ALL2

#Fst permutation
for j in Other Pekin-SX Pekin-Sp chr
do
{
for i in `seq 1 1000`
do
bedtools shuffle -i <(sort -k1,1 -k2,2n ../2022-Final-Three-Class|grep $j|grep -v chrZ|awk '{print $1"\t"$2"\t"$3}') -g <(grep -v HiC ../../../../../1.With-W/CAU-length|grep -v ZJU-W|grep -v Z) > ${i}-${j}-SHUFFLE-CAU-cLTR.bed
bedtools flank -i ${i}-${j}-SHUFFLE-CAU-cLTR.bed -g <(grep -v HiC ../../../../../1.With-W/CAU-length|grep -v ZJU-W|grep -v Z) -b 100000 > ${i}-${j}-SHUFFLE-CAU-cLTR-0-100k.bed
bedtools flank -i ${i}-${j}-SHUFFLE-CAU-cLTR.bed -g <(grep -v HiC ../../../../../1.With-W/CAU-length|grep -v ZJU-W|grep -v Z) -b 10000 > ${i}-${j}-SHUFFLE-CAU-cLTR-0-10k.bed
bedtools flank -i ${i}-${j}-SHUFFLE-CAU-cLTR.bed -g <(grep -v HiC ../../../../../1.With-W/CAU-length|grep -v ZJU-W|grep -v Z) -b 5000 > ${i}-${j}-SHUFFLE-CAU-cLTR-0-5k.bed
bedtools intersect -a ${i}-${j}-SHUFFLE-CAU-cLTR-0-100k.bed -b  2Pek-IND-Fst.weir.fst -wo|awk '{sum+=$7}END{if(NR>0)print "'"$j"'","-100-100",sum,NR,sum/NR;else print "'"$j"'","-100-100",sum,NR,"0"}' >>  Pek-IND-Permuation1000
bedtools intersect -a ${i}-${j}-SHUFFLE-CAU-cLTR-0-10k.bed -b  2Pek-IND-Fst.weir.fst -wo|awk '{sum+=$7}END{if(NR>0)print "'"$j"'","-10-10",sum,NR,sum/NR;else print "'"$j"'","-10-10",sum,NR,"0"}' >>  Pek-IND-Permuation1000
bedtools intersect -a ${i}-${j}-SHUFFLE-CAU-cLTR-0-5k.bed -b  2Pek-IND-Fst.weir.fst -wo|awk '{sum+=$7}END{if(NR>0)print "'"$j"'","-5-5",sum,NR,sum/NR;else print "'"$j"'","-5-5",sum,NR,"0"}' >>  Pek-IND-Permuation1000
bedtools intersect -a ${i}-${j}-SHUFFLE-CAU-cLTR-0-100k.bed -b  2Ma-Pek_IND-Fst.weir.fst -wo|awk '{sum+=$7}END{if(NR>0)print "'"$j"'","-100-100",sum,NR,sum/NR;else print "'"$j"'","-100-100",sum,NR,"0"}' >>  Ma-home-Permuation1000
bedtools intersect -a ${i}-${j}-SHUFFLE-CAU-cLTR-0-10k.bed -b  2Ma-Pek_IND-Fst.weir.fst -wo|awk '{sum+=$7}END{if(NR>0)print "'"$j"'","-10-10",sum,NR,sum/NR;else print "'"$j"'","-10-10",sum,NR,"0"}' >>  Ma-home-Permuation1000
bedtools intersect -a ${i}-${j}-SHUFFLE-CAU-cLTR-0-5k.bed -b  2Ma-Pek_IND-Fst.weir.fst -wo|awk '{sum+=$7}END{if(NR>0)print "'"$j"'","-5-5",sum,NR,sum/NR;else print "'"$j"'","-5-5",sum,NR,"0"}' >>  Ma-home-Permuation1000 
rm ${i}-${j}-SHUFFLE-CAU-cLTR*bed
done
} &
done



#variation in differ populations
for i in `ls *ALl*Populaiton*ID`; do
cat CAU-cLTR-0-100k.bed-${i}.map.bed-Intersect)｜awk '{if($5!=-1)print $0}'｜awk -v CHR=0 -v ID=0 -v Direct=1 '{if(ID!=$2 || CHR!=$1){Direct=Direct*(-1);if(sum>0)print CHR,ID,sum/count,count;CHR=$1;ID=$2;sum=0;count=1};count+=1;if(Direct==-1)sum+=$3-$5;else sum+=$5-$2}'  > CAU-cLTR-0-100k.bed-${i}.map.bed-Intersect3
cat CAU-cLTR-0-100k.bed-${i}.map.bed-Intersect)｜awk '{if($5!=-1)print $0}'｜awk -v CHR=0 -v ID=0 -v Direct=1 '{if(ID!=$2 || CHR!=$1){Direct=Direct*(-1);CHR=$1;ID=$2;if(sum>0)print sum/count,count;sum=0;count=1};count+=1;if(Direct==-1)sum+=$3-$5;else sum+=$5-$2}'|awk '{sum+=$1}END{print sum,sum/NR}' > CAU-cLTR-0-100k.bed-${i}.map.bed-Intersect2
done

while read line; do Count=`expr $Count + 1`; UpChr=`echo $line|awk '{print $1}'`; UpStart=`echo $line|awk '{print $2}'`; DownChr=`echo $line|awk '{print $4}'`; DownStart=`echo $line|awk '{print $5}'`; awk '{if($1=="'"$UpChr"'" && $2=="'"$UpStart"'")print $1"-'"$Count"'",$6-$3;else if($1=="'"$DownChr"'" && $2=="'"$DownStart"'")print $1"-'"$Count"'",$6-$2}' CAU-cLTR-0-100k.bed-Ma.map.bed-Intersect >> CAU-cLTR-0-100k.bed-Ma.map.bed-Intersect4; done <CAU-cLTR-0-100k.bed2 &

#Fst calculate
~/plink2/plink --bfile ~/2.analysis/10.Resequencing/2.185/Clean/QC-Maf0.01-geno0.2/Final-Autosome/180-Autosome --keep Ma-IND-Pekin-ID --recode --maf 0.01 --out Ma-IND-Pekin --chr-set 40
~/plink2/plink --bfile ~/2.analysis/10.Resequencing/2.185/Clean/QC-Maf0.01-geno0.2/Final-Autosome/180-Autosome --keep IND-Pekin-ID --recode --maf 0.01 --out IND-Pekin --chr-set 40 &
~/plink2/plink --file Ma-IND-Pekin --recode vcf --out Ma-IND-Pekin --chr-set 40
~/plink2/plink --file IND-Pekin --recode vcf --out IND-Pekin --chr-set 40 
vcftools --vcf Ma-IND-Pekin.vcf --weir-fst-pop Pek-IND-ID2 --weir-fst-pop Ma-ID2 --out Ma-Pek_IND-Fst &
vcftools --vcf IND-Pekin.vcf --weir-fst-pop Pek-ID2 --weir-fst-pop IND-ID2 --out Pek-IND-Fst &
#
awk 'NR>1{if($3>0){if($1<10)print "chr0"$1"\t"$2-1"\t"$2"\t"$3;else{print "chr"$1"\t"$2-1"\t"$2"\t"$3}}}' Pek-IND-Fst.weir.fst > 2Pek-IND-Fst.weir.fst 

