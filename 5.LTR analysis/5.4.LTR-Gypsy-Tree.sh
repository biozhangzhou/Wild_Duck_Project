#LTR Gypsy
for i in `ls *All*genome*ID`
do
awk '{if($3=="long_terminal_repeat")print $0}' ./${i}*EDTA*raw/LTR/*.mod.LTR.intact.gff3 > $i\-LTR; done
awk '{if($1~/Chr/ && $7=="+" && $9~/Gypsy/)print $1"\t"$4-1"\t"$5-1"\t"$7,$9}' $i\-LTR > $i\-LTR-F
awk '{if($1~/Chr/ && $7=="-" && $9~/Gypsy/)print $1"\t"$4-1"\t"$5-1"\t"$7,$9}' $i\-LTR > $i\-LTR-R
done
#extract fasta
for i in `ls *All*genome*ID`; do seqkit subseq --bed $i\-LTR-F ${i}.fa|seqkit seq -w0 > $i\-LTR-F.fa &  seqkit subseq --bed $i\-LTR-R ${i}.fa|seqkit seq -w0 -r -p  > $i\-LTR-R.fa & done

#Combine fasta
for i in `ls *All*genome*ID`
do
awk '{if($0~/>/){if($0~/lLTR/)print $1"-lLTR";else print $1"-rLTR"}else print $0}' $i\-LTR-F.fa|sed 's/>/>'"$i"'-/g'|sed 's/\:\.//g' >> ALL-LTR.fa
awk '{if($0~/>/){if($0~/lLTR/)print $1"-lLTR";else print $1"-rLTR"}else print $0}' $i\-LTR-R.fa|sed 's/>/>'"$i"'-/g'|sed 's/\:\.//g' >> ALL-LTR.fa
done
awk '{if($1~/>/)printf $1"\t";else print $0}' ALL-LTR.fa |sort|uniq|sed 's/\t/\n/g' > 2ALL-LTR.fa 
mafft --thread 72 2ALL-LTR.fa > ALL-LTR.msa
iqtree2 -s ALL-LTR.msa -m K80 -T 100
