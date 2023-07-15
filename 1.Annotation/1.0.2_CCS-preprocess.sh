##
for outputfile in 1002_bu 1002
do
lima ${outputfile}.ccs.bam  IsoSeqPrimers.fasta ${outputfile}.fl.bam --isoseq --peek-guess --num-threads 40
isoseq3 refine ${outputfile}.fl.primer_5p--primer_3p.bam IsoSeqPrimers.fasta ${outputfile}.flnc.bam --require-polya
isoseq3 cluster ${outputfile}.flnc.bam ${outputfile}_clustered.bam --verbose --use-qvs
done
#polish
isoseq3 polish ${outputfile}_clustered.bam ${outputfile}.subreads.bam  ${outputfile}.polished.bam

#Get_Coding_Seqs
~/1.Fulllength_processing.sh -a 1002.polished.hq.fasta -o 1002 -t 10

############qsub#Cupcake ToFU#############################################################
#!/bin/bash
#PBS -N corre_collapsed
#PBS -l nodes=1:ppn=15,mem=15gb
#PBS -e /home/qsub_log
#PBS -o /home/qsub_log
#PBS -q cu
#PBS -t 0-1


# Kill script if any commands fail
set -e
source /opt/software/anaconda2/bin/activate

collapse=/opt/software/cDNA_Cupcake/cupcake/tofu
fa2fq=/opt/software/cDNA_Cupcake/sequence

#Setup tempdisk for output
uid="zhangzhou" #user id
ls_date=`date +m%d%H%M%S` #Set Random date and time
Work_dir=${uid}_${ls_date}_${PBS_ARRAYID}  #temp directory
cd /tmpdisk/ 
mkdir ${Work_dir} && cd ${Work_dir}

DataID=(`ls ~/*1006*fasta.gz|awk -F\/ '{print $9}'|sed 's/.polished.hq.fasta.gz//g'`)

cp ~/${DataID[$PBS_ARRAYID]}/Coding.fa ./ccs.fa

inputfile=ccs.fa
outputfile=ccs.gmap

gmap -t 15 -D ~/Genome/ -d Ma03  -f samse -n 1 ${inputfile}  > ${outputfile}.sam
sort -k 3,3 -k 4,4n ${outputfile}.sam > ${outputfile}.sorted.sam

python ${collapse}/collapse_isoforms_by_sam.py --input ${inputfile} -s ${outputfile}.sorted.sam --dun-merge-5-shorter -o ${outputfile}

python ${collapse}/get_abundance_post_collapse.py ${outputfile}.collapsed ~/0.Isoseq-out/${DataID}.polished.cluster_report.csv

python ${fa2fq}/fa2fq.py ${outputfile}.collapsed.rep.fa

mv ${outputfile}.collapsed.rep.fastq ${outputfile}.collapsed.rep.fq

python ${collapse}/filter_away_subset.py ${outputfile}.collapsed

python ${fa2fq}/fq2fa.py ${outputfile}.collapsed.filtered.rep.fq
###################END#############END######################################################

#Merge
~/1.bin/gffcompare-0.11.6/gffcompare -i Ma_all_isoforms.txt -o Ma_all_isoforms
./cufflinks2gff3 BT-1002.ccs.gmap.collapsed.filtered.gff > BT-1002.ccs.gmap.collapsed.filtered.maker.gff &






