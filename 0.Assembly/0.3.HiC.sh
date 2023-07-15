#Hic scaffolding #
#3.1 QC
fastp -f 10 -t 30 -F 10 -T 30 -i file_R1 -o file_clean_R1 -I file_R2 -O file_clean_R2
#3.2 File prepare
SampleName=Sample
mkdir fastq && cd fastq
ln -s ${file_clean_R1} ./hic_R1.fastq
ln -s ${file_clean_R2} ./hic_R2.fastq
cd ..
mkdir reference && cd reference
ln -s ${Bionano_out.fasta} ./${SampleName}.fasta
bwa index ${SampleName}.fasta
python2.7 /path/to/juicer/misc/generate_site_positions.py DpnII ${SampleName} ${SampleName}.fasta
awk 'BEGIN{OFS="\t"}{print $1, $NF}' ${SampleName}_DpnII.txt > ${SampleName}.chrom.size
cd ..
#3.3 Run juicer
ln -s ~/1.bin/juicer-master/CPU/ ./scripts

/path/to/juicer/scripts/juicer.sh \
-g genome \
-s DpnII \
-z reference/${SampleName}.fa \
-y reference/${SampleName}_DpnII.txt \
-p reference/${SampleName}.chrom.size \
-D /path/to/scripts/ \
-t 72 &> juicer.log &
#3.4 Run 3d-dna
#3.4.1
mkdir 3d-dna && cd 3d-dna
/path/to/3d-dna_modify/run-asm-pipeline.sh \
-m haploid \
-i 10000 \
-r 0 \
../references/${SampleName}.fasta \
/path/to/juicer_output/aligned/merged_nodups.txt
#3.4.2
Visualization by juicerBox and Revision
#3.4.3
mkdir post-review && post-review
/path/to/3d-dna_modify/run-asm-pipeline-post-review_zz.sh \
-r ${SampleName}.rawchrom.review.assembly \
../../references/${SampleName}.fasta \
/path/to/juicer_output/aligned/merged_nodups.txt