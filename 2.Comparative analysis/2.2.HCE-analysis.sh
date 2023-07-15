#spliting maf file
for i in `ls CAU-Chr*-25-Bird.maf`; do { sed '2d' $i > ../2.PhyloFit/split-maf/tmp-$i; }& done
for i in `seq 1 40` Z; do mkdir $i; mv tmp-CAU-Chr$i-25-Bird.maf $i; done

~/1.bin/phast/bin/phyloFit --tree 2ALL-except-WZ-Filter-CAU-Chr-25-Bird.4d.fa.raxml.bestTree -i FASTA 2ALL-except-WZ-Filter-CAU-Chr-25-Bird.4d.fa --out-root nonconserved

for i in `seq 1 40` Z;do cd $i; { for j in `ls *`; do ~/1.bin/phast/bin/msa_split ${j} --in-format MAF  --windows 1000000,0 --out-root split --out-format SS --min-informative 1000 --between-blocks 5000; done; }& cd ..; done

for i in `seq 1 115`; do { for j in `seq 1 10`; do id=`expr $i \* 10 + $j`; echo $id; ss=`sed -n ''"$id"','"$id"'p' ss.list `; ~/1.bin/phast/bin/phastCons --estimate-rho $ss --no-post-probs $ss ../nonconserved-4d-site.mod; done; }& done

ls */*.noncons.mod > all.noncons.mod.list
ls */*.cons.mod > all.cons.mod.list

~/1.bin/phast/bin/phyloBoot --read-mods '*all.cons.mod.list' --output-average ave.cons.mod
~/1.bin/phast/bin/phyloBoot --read-mods '*all.noncons.mod.list' --output-average ave.noncons.mod

#PhastCons
 for i in `seq 1 41`
 do 
 { id=`sed -n ''"$i"','"$i"'p' maf-list`
~/1.bin/phast/bin/phastCons --most-conserved Birds.${i}\.bed --score $id ave.cons.mod,ave.noncons.mod > ${i}\.wig
}& done

#Total conserved region
cat Birds.*.bed|awk '{sum+=$3-$2+1}END{print sum}'
212778139
#
for i in `seq 1 40` Z; do ~/1.bin/phast/bin/phyloP --subtree duck --msa-format MAF --features ../Birds.${i}\.bed --method LRT --mode CONACC ../../rooted-nonconserved-4d-site.mod ../$i/tmp-CAU-Chr$i-25-Bird.maf > Birds.${i}\.bed-duck-element-scores.txt & done
for i in `seq 1 40` Z; do ~/1.bin/phast/bin/phyloP --subtree home --msa-format MAF --features ../Birds.${i}\.bed --method LRT --mode CONACC ../../rooted-nonconserved-4d-site.mod ../$i/tmp-CAU-Chr$i-25-Bird.maf > Birds.${i}\.bed-home-element-scores.txt & done

for i in `seq 1 40` Z;do awk '{if($NF>=-0.05 && $NF<=0)print $0}' Birds.$i\.bed-element-scores.txt > Acc.$i\.bed; awk '{if($NF<=0.05 && $NF>=0)print $0}' Birds.$i\.bed-element-scores.txt > Con.$i\.bed; done
for i in `seq 1 40` Z; do awk '{if($NF>=-0.05 && $NF<=0)print $0}' Birds.$i\.bed-duck-element-scores.txt > Acc.$i\.bed; awk '{if($NF<=0.05 && $NF>=0)print $0}' Birds.$i\.bed-duck-element-scores.txt > Con.$i\.bed; done

for i in `seq 1 40` Z; do perl maf_bed_filter_Type1.pl Birds.$i\.bed ../../../0.Datas-Orthgene-4dsites/0.Maf/CAU-Chr$i\-25-Bird.maf.lst ./workout-$i &  done
for i in `seq 1 40` Z; do perl maf_bed_filter_Type2.pl 2.Type1/Duck/Acc.$i\.bed ../../../0.Datas-Orthgene-4dsites/0.Maf/CAU-Chr$i-25-Bird.maf.lst Acc-workout-$i/ & done
for i in `seq 1 40` Z; do perl maf_bed_filter_Type2.pl 2.Type1/Home/Acc.$i\.bed ../../../0.Datas-Orthgene-4dsites/0.Maf/CAU-Chr$i-25-Bird.maf.lst Acc-workout-$i/ & done
for i in `seq 1 40` Z; do perl maf_bed_filter_Type2.pl 2.Type1/Duck/Con.$i\.bed ../../../0.Datas-Orthgene-4dsites/0.Maf/CAU-Chr$i-25-Bird.maf.lst workout-$i/ & done
for i in `seq 1 40` Z; do perl maf_bed_filter_Type2.pl 2.Type1/Home/Con.$i\.bed ../../../0.Datas-Orthgene-4dsites/0.Maf/CAU-Chr$i-25-Bird.maf.lst workout-$i/ & done

cat workout-*/*bed > Con.bed 
cat workout-*/*bed > Acc.bed

awk '{if($3-$2 >=19)sum+=$3-$2+1}END{print sum}' Acc.bed;awk '{if($3-$2 >=19)print $0}' Acc.bed
awk '{if($3-$2 >=99)sum+=$3-$2+1}END{print sum}' Acc.bed;awk '{if($3-$2 >=99)print $0}' Acc.bed


#distribution of HCE
for i in `ls *.bed`
do
for j in `ls CAU-bed/CAU-**bed|grep -v TE|sed 's/CAU-bed\///g'`
do
	sed 's/chr0/Chr/g' CAU-bed/${j}|sed 's/chr/Chr/g' > Tmp.bed
	bedtools intersect -a Tmp.bed -b $i -wo|awk -v DES=0 -v OLD=NULL -v sum=0 '{des=$1"-"$2"-"$3;if(des==DES)sum+=$NF;else{print OLD"\t"sum;OLD=$1"\t"$2"\t"$3;DES=$1"-"$2"-"$3;sum=$NF}}'|sed '1d' > $j\-${i}-len
done
done
#
for i in `ls *.bed`
do
awk '{sum+=$3-$2+1}END{print "'"$i"'\tsum\t"sum}' ${i} > ${i}-Distribution
awk '{sum+=$4}END{print "'"$i"'\tgene\t"sum}' HCE-Distribution/CAU-gene.bed-${i}-len >> ${i}-Distribution
awk '{sum+=$4}END{print "'"$i"'\tCDS\t"sum}' HCE-Distribution/CAU-CDS.bed-${i}-len >> ${i}-Distribution
awk '{sum+=$4}END{print "'"$i"'\texon\t"sum}' HCE-Distribution/CAU-exon.bed-${i}-len >> ${i}-Distribution
awk '{sum+=$4}END{print "'"$i"'\tintron\t" sum}' HCE-Distribution/CAU-intro.bed-${i}-len >> ${i}-Distribution
awk '{sum+=$4}END{print "'"$i"'\t0-10k\t"sum}' HCE-Distribution/CAU-gene10k.bed-${i}-len >> ${i}-Distribution
awk '{sum+=$4}END{print "'"$i"'\t10-20k\t"sum}' HCE-Distribution/CAU-gene20k.bed-${i}-len >> ${i}-Distribution
awk '{sum+=$4}END{print "'"$i"'\t20-30k\t"sum}' HCE-Distribution/CAU-gene30k.bed-${i}-len >> ${i}-Distribution
awk '{sum+=$4}END{print "'"$i"'\t30-40k\t"sum}' HCE-Distribution/CAU-gene40k.bed-${i}-len >> ${i}-Distribution
awk '{sum+=$4}END{print "'"$i"'\t40-50k\t"sum}' HCE-Distribution/CAU-gene50k.bed-${i}-len >> ${i}-Distribution
done

awk '{FS="\t";OFS="\t"}{$4=$4-10000;$5=$5+10000;if($4<0){$4=0}print $0}' CAU-GFF-Gene.gff > CAU-GFF-Gene-up10k-down10k.gff 
for i in `ls Type*bed`; do bedtools intersect -a ../CAU-GFF-Gene.gff -b $i -wa|uniq|awk -F\; '{print $2}'|sort|sed 's/gene=//g' > $i-Gene; done
for i in `ls Type*bed`; do bedtools intersect -a ../CAU-GFF-Gene-up10k-down10k.gff -b $i -wa|uniq|awk -F\; '{print $2}'|sort|sed 's/gene=//g' > $i-Gene-Up10-Down10k; done




