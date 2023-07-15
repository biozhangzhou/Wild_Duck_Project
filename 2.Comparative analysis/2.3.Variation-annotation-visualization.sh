#construct database
~/gffread-0.11.8/gffread Final.gff -T -o Final.gtf
209:gtfToGenePred -genePredExt Final.gtf Final.txt
perl retrieve_seq_from_fasta.pl -format refGene -seqfile Genoma.fasta Final.txt --outfile Final.mRNA.fa

#variants filteration #--maf 0.5
vcftools --remove-filtered-all --maf 0.5 --gzvcf Sample_filter_variants.vcf.gz --recode --stdout | gzip -c > Sample_Pass.vcf.gz &

#format transform
perl ~/1.bin/annovar/convert2annovar.pl -format vcf4 Sample_Pass.vcf.gz -out Sample-annovar

#variants annotation
i=Sample
perl ~/1.bin/annovar/table_annovar.pl $i\-annovar Anno_db/$i/ -buildver $i -out $i --otherinfo -remove -protocol refGene -operation g -nastring NA

#Merge
for i in `ls *All*Sample*`; do awk '{OFS="\t";FS="\t"}{if(NR>1)print $0,"'"$i"'"}' ${i}\_multianno.txt >> TTT; done

#visualization

library(maftools)
var.annovar.maf = annovarToMaf(annovar="TTT")
laml = read.maf(maf = var.annovar.maf)
#plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE,textSize = 0.5,showBarcodes=TRUE)
dev.off()
