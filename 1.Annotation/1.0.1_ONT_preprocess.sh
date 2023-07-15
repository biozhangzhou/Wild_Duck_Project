##
#1 get full length seqs
~/cdna_classifier.py -r report.pdf -u unclassified.fq -w rescued.fq ${fq} full_length_output.fq
#2 ngs polish
./fmlrc [options] <comp_msbwt.npy> <long_reads.fa> <corrected_reads.fa>
#3 remove redundent transcripts
cd-hit-est -i correct_full_length_output.fq -o nr_correct_full_length_output.fa -M 20000 -c 0.95 -T $thread -G 0 -aL 0.00 -aS 0.99
#4 define coding or non-coding transcripts
#4.1 PLEK
#seqkit split -p 10 nr_correct_full_length_output.fa -j $thread
python2.7 ~/PLEK.1.2/PLEK.py -minlength 200 -fasta nr_correct_full_length_output.part_010.fa -out Sample010_predicted -thread 10 &
#4.2 CNCI
python2.7 ~/CNCI-master/CNCI.py -f $nr_INFile -o cnci_predict -p $thread -m ve
#4.3 CPC
~/py27/bin/CPC2.py -i $nr_INFile -o CPC.txt

cat PLEK_predicted|grep ^Non-coding|awk '{print $3}'|sed 's/^>//g'  >> PLEK_Noncoding_ID
cat PLEK_predicted|grep ^Coding|awk '{print $3}'|sed 's/^>//g'  >> PLEK_Coding_ID

cat cnci_predict/CNCI.index|grep -w coding |awk '{print $1}' > CNCI_Coding_ID
cat cnci_predict/CNCI.index|grep -w noncoding |awk '{print $1}' > CNCI_Noncoding_ID

cat CPC.txt |grep -w coding |awk '{print $1}' >> CPC_Coding_ID
cat CPC.txt |grep -w noncoding |awk '{print $1}' >> CPC_Noncoding_ID

comm -12 <(sort CNCI_Coding_ID) <(sort PLEK_Coding_ID) > CNCI_PLEK_Coding_ID
comm -12 <(sort CPC_Coding_ID) <(sort CNCI_PLEK_Coding_ID) > ../CNCI_PLEK_CPC_Coding_ID

comm -12 <(sort CNCI_Noncoding_ID) <(sort PLEK_Noncoding_ID) > CNCI_PLEK_Noncoding_ID
comm -12 <(sort CPC_Noncoding_ID) <(sort CNCI_PLEK_Noncoding_ID) > ../CNCI_PLEK_CPC_Noncoding_ID

#5 get Coding seqs
seqkit grep -f CNCI_PLEK_CPC_Coding_ID -w 0 $nr_INFile -o Coding.fa
#6 minimap mapping
minimap2 -ax splice -t 10 Sample_polish.fasta FINAL_Coding.fa
~/bin/samtools view -Su ${Output}/${ID}_2.sam  -@ 10 |/opt/software/anaconda2/bin/samtools sort - -o ${Output}/${ID}_2.sorted.bam -@ 10
#7 stringtie to get gtf and merge
/opt/software/stringtie2/stringtie Sample_Filter_LncRNA_XX.sorted.bam -o Sample_Filter_LncRNA_XX.GTF -p 20 -L
#8 merge
~/gffcompare-0.11.6/gffcompare -i $i\-GFT-List -o $i\-ALL
#9 convert to MAKER format
~/cufflinks2gff3 2${i}-ALL.combined.gtf > 2${i}\-ALL.combined.maker.gtf




