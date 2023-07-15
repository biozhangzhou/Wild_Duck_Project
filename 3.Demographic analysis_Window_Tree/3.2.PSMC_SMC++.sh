#PSMC
for i in `ls *All*Sample`
do
  bcftools mpileup -C50 -f $i\_genome.fa $i\_mkdup.bam | bcftools call -c - | vcfutils.pl vcf2fq -d 50 -D 200 | gzip > $i\.diploid.fq.gz &
  ~/psmc-0.6.5/utils/fq2psmcfa -q20 $i\.diploid.fq.gz > $i\.diploid.psmcfa
  ~/psmc-0.6.5/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o $i\.diploid.psmc $i\.diploid.psmcfa
done

perl ~/1.bin/psmc-0.6.5/utils/psmc_plot.pl -g 1 -u 1.91e-9 -M "PZ,ZW,LC,LW,HL,Ma,Ma2,Bs,Bs2,Sx" -p ALL PZ.diploid.psmc ZW.diploid.psmc LC.diploid.psmc LW.diploid.psmc HL.diploid.psmc Ma01.diploid.psmc Ma03.diploid.psmc BT02.diploid.psmc BS02.diploid.psmc SX31.diploid.psmc 



#SMC++
for i in `seq 1 40`; do smc++ vcf2smc chr$i\-Derived.SNP_filter.vcf.gz $i\.smc.txt $i SB:BS13_BS13,BS14_BS14,BS15_BS15,BS16_BS16,BS17_BS17,BS18_BS18,SB01_SB01,SB02_SB02,SB08_SB08,SB09_SB09; done

smc++ estimate -o SB-out/ 1.91e-9 smc-SB/*smc.txt --cores 120 --timepoints 1 1000000 -c 50000

smc++ plot plot-Ma-SB-Pek-IND.png Ma-out/model.final.json SB-out/model.final.json IND-out/model.final.json Pek-out/model.final.json

#Z
smc++ vcf2smc chrZ-Derived.SNP_filter.vcf.gz IND-Female-Z.smc.txt 0 IND:SRR7091449_SRR7091449,SRR7091452_SRR7091452,SRR7091456_SRR7091456,SX31_SX31,SRR7091511_SRR7091511
smc++ vcf2smc chrZ-Derived.SNP_filter.vcf.gz Pek-Female-Z.smc.txt 0 Pek:CAU_CAU,SRR7091411_SRR7091411,SRR7091412_SRR7091412,SRR7091414_SRR7091414,SRR7091416_SRR7091416
smc++ vcf2smc chrZ-Derived.SNP_filter.vcf.gz Ma-Female-Z.smc.txt 0 Ma:Ma03_Ma03,Ma04_Ma04,Ma05_Ma05,Ma06_Ma06,Ma07_Ma07
smc++ vcf2smc chrZ-Derived.SNP_filter.vcf.gz SB-Female-Z.smc.txt 0 SB:SB02_SB02,SB08_SB08,SB09_SB09,BS13_BS13,BS18_BS18


for i in `seq 1 39`;
do
    smc++ vcf2smc chr${i}-Derived.SNP_filter.vcf.gz Ma-SB-split/${i}.smc.txt $i Ma:Ma01_Ma01,Ma02_Ma02,Ma03_Ma03,Ma04_Ma04,Ma05_Ma05,Ma06_Ma06,Ma07_Ma07,Ma14_Ma14,Ma15_Ma15,Ma16_Ma16 SB:BS13_BS13,BS14_BS14,BS15_BS15,BS16_BS16,BS17_BS17,BS18_BS18,SB01_SB01,SB02_SB02,SB08_SB08,SB09_SB09;
done
for i in `seq 1 39`;
do
    smc++ vcf2smc chr${i}-Derived.SNP_filter.vcf.gz Ma-SB-split/${i}.smc.txt $i Ma:Ma01_Ma01,Ma02_Ma02,Ma03_Ma03,Ma04_Ma04,Ma05_Ma05,Ma06_Ma06,Ma07_Ma07,Ma14_Ma14,Ma15_Ma15,Ma16_Ma16 SB:BS13_BS13,BS14_BS14,BS15_BS15,BS16_BS16,BS17_BS17,BS18_BS18,SB01_SB01,SB02_SB02,SB08_SB08,SB09_SB09;
done
smc++ split -o Ma-SB-split/ Ma/model.final.json SB/model.final.json Ma-SB-split/*smc.txt

for i in `seq 1 39`; do smc++ vcf2smc chr$i\-Derived.SNP_filter.vcf.gz Ma-Pek-split/$i\.smc.txt $i Pek:CAU_CAU,SRR7091411_SRR7091411,SRR7091412_SRR7091412,SRR7091413_SRR7091413,SRR7091414_SRR7091414,SRR7091415_SRR7091415,SRR7091416_SRR7091416,SRR7091417_SRR7091417,SRR7091419_SRR7091419,SRR7091431_SRR7091431 Ma:Ma01_Ma01,Ma02_Ma02,Ma03_Ma03,Ma04_Ma04,Ma05_Ma05,Ma06_Ma06,Ma07_Ma07,Ma14_Ma14,Ma15_Ma15,Ma16_Ma16; done &
for i in `seq 1 39`; do smc++ vcf2smc chr${i}-Derived.SNP_filter.vcf.gz Ma-IND-split/${i}.smc.txt $i IND:SRR7091447_SRR7091447,SRR7091448_SRR7091448,SRR7091449_SRR7091449,SRR7091450_SRR7091450,SRR7091451_SRR7091451,SRR7091452_SRR7091452,SRR7091456_SRR7091456,SRR7091511_SRR7091511,SX04_SX04,SX31_SX31 Ma:Ma01_Ma01,Ma02_Ma02,Ma03_Ma03,Ma04_Ma04,Ma05_Ma05,Ma06_Ma06,Ma07_Ma07,Ma14_Ma14,Ma15_Ma15,Ma16_Ma16;  done &
for i in `seq 1 39`; do smc++ vcf2smc chr$i\-Derived.SNP_filter.vcf.gz Pek-IND-split/$i\.smc.txt $i Pek:CAU_CAU,SRR7091411_SRR7091411,SRR7091412_SRR7091412,SRR7091413_SRR7091413,SRR7091414_SRR7091414,SRR7091415_SRR7091415,SRR7091416_SRR7091416,SRR7091417_SRR7091417,SRR7091419_SRR7091419,SRR7091431_SRR7091431 IND:SRR7091447_SRR7091447,SRR7091448_SRR7091448,SRR7091449_SRR7091449,SRR7091450_SRR7091450,SRR7091451_SRR7091451,SRR7091452_SRR7091452,SRR7091456_SRR7091456,SRR7091511_SRR7091511,SX04_SX04,SX31_SX31;done &

#Z
smc++ split -o Ma-Pek-Z-split/ Ma-Z-out/model.final.json Pek-Z-out/model.final.json Ma-Pek-Z-split/*smc.txt
smc++ split -o Ma-IND-Z-split/ Ma-Z-out/model.final.json Pek-Z-out/model.final.json Ma-IND-Z-split/*smc.txt
smc++ split -o Pek-IND-Z-split/ Pek-Z-out/model.final.json IND-Z-out/model.final.json Pek-IND-Z-split/*smc.txt
smc++ split -o Ma-SB-Z-split/ Ma-Z-out/model.final.json SB-Z-out/model.final.json Ma-SB-Z-split/*smc.txt
smc++ plot Pek-IND-Z-split.png Pek-IND-Z-split/model.final.json
smc++ plot Ma-Pek-Z-split.png Ma-Pek-Z-split/model.final.json
smc++ plot Ma-SB-Z-split.png Ma-SB-IND-Z-split/model.final.json
smc++ plot Ma-IND-Z-split.png Ma-IND-Z-split/model.final.json

