#splitting
#1M 500K 100K
for i in `seq 1 30` Z
do
cd Window-Tree/100K
mkdir $i && cd $i
~/phast/bin/msa_split ./${i}/*maf --in-format MAF  --windows 1000000,0 --out-root 1M/split --out-format FASTA --min-informative 1000 --between-blocks 5000 &
~/phast/bin/msa_split ./${i}/*maf --in-format MAF  --windows 500000,0 --out-root 500K/split --out-format FASTA --min-informative 1000 --between-blocks 5000 &
~/phast/bin/msa_split ./${i}/*maf --in-format MAF  --windows 100000,0 --out-root 100K/split --out-format FASTA --min-informative 1000 --between-blocks 5000 &
done

#tree building
for j in 1M-list 500k-list 100k-list
for i in `cat $j`
do
ID=`echo $i|awk -F\/ '{print $(NF-1)"-"$NF}'`;
seqkit grep -v -r -p 'Ana|ZJU' $i > $ID
~/1.bin/Gblocks_0.91b/Gblocks $ID -b1
seqkit seq -w 0 $ID-gb|sed 's/ //g'|sed 's/\*/-/g'> $ID-gb.fa 
~/1.bin/raxml-ng/raxml-ng --msa $ID\-gb.fa --model GTR+G --threads 40 
rm $ID $ID-gb $ID-gb.htm
done
done

#10k 1k
for i in `cat 100k-list`
do
ID=`echo $i|awk -F\/ '{print $(NF-1)"-"$NF}'`;
seqkit grep -v -r -p 'Ana|ZJU' $i|seqkit seq -w 0 > $ID
~/1.bin/phast/bin/msa_split ${ID} --in-format FASTA  --windows 1000,0 --out-root 1K/tmp-split --out-format FASTA --min-informative 30 --between-block 0
for i in `ls 1K/tmp-split* `
do
#keep or filter
#10K window 9k nogap #1K window 900bp nogap
#gapratio <=10%
KOF=`seqkit seq -w 0 $i|sed 's/*//g'|sed 's/[-N]//g'|grep -v \>|awk '{if(length($0)<900)filter="YES"}END{if(filter=="YES")print "filter";else print "keep"}'` 

if [ $KOF == "filter" ];then 
rm $i
else sed -i 's/*/-/g' $i
fi
done
for i in `ls tmp-1K-split*`
do
#consensus qc
POLY=`seqkit seq -w 0 $i|grep -v \>|awk '{if(NR==1){split($0,Seq,"");for(i=1;i<=length(Seq);i++) a[i,1]=Seq[i]}else{split($0,Seq,"");for(i=1;i<=length(Seq);i++){if(a[i,1]!=Seq[i]){a[i,1]="polymorphism"}}}}END{for(i=1;i<=length(Seq);i++){print a[i,1]}}'|grep polymorphism|wc|awk '{print $1}'`
iqtree2 -s $i -m MFP -T 4
#Failed or Pass
FOP=`grep 'sequences failed composition' ${i}.log|awk '{if($4>0)print "Failed";else print "PASS"}'`
if [ $FOP == "PASS" ];then
cat ${i}.treefile|awk '{print "'"$j"'-'"$i"' '"$POLY"'",$0}' >> 10
fi
rm ${i}*
done
rm ${ID}
done


#statistic
nw_reroot <(cat *ALL_Tree_In_Specific_Window*) CAU|nw_reroot - HL PZ|nw_topology -|nw_order -|sort|uniq -c|sort -k1,1nr > TreeFile-stat
