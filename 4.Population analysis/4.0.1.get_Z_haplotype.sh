#!/bin/bash

while read line
do
echo $line
ls_date=`date +m%d%H%M%S`
ID=`echo $line|awk '{print $1}'`
echo $line|awk '{print $1,$1}' > tmp
/opt/software/plink2/plink --bfile CHRZ --keep tmp --out tmp --make-bed --chr-set 40
/opt/software/plink2/plink --bfile tmp --maf 0.5 --out tmp-maf0.5 --recode --chr-set 40
awk '{print $2}' tmp-maf0.5.map > keep
/opt/software/plink2/plink --bfile tmp --exclude keep --make-bed --out tmp-F-Homo --chr-set 40
mv tmp-F-Homo.bed ${ID}-F-Homo.bed
mv tmp-F-Homo.fam ${ID}-F-Homo.fam
mv tmp-F-Homo.bim ${ID}-F-Homo.bim
mv tmp-F-Homo.log ${ID}-F-Homo.log
mv tmp-F-Homo.nosex ${ID}-F-Homo.nosex
done<180-ID-Sex-Female
~                       