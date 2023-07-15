#!/bin/bash

#NGS Reads
scallop2 -i MA3-$ID\_sort.bam -o MA3-$ID\.gtf 
stringtie --merge Ma03-BA.gtf Ma03-BM.gtf Ma03-LI.gtf Ma03-LU.gtf Ma03-SP.gtf -p 20 -o Ma_NGS_scallop_merged.gtf

#######
#get trans directly from long reads.
~/gffread-0.11.8/gffread Ma_all_isoforms.combined.gtf -g Ma_Ref.fa -w Ma_all_isoforms.combined.gff3.exon.fa
#ORFfinder
~/1.bin/maker/exe/ORFfinder -in Ma_all_isoforms.combined.gff3.exon.fa -s 0 -out Ma_all_isoforms.combined.gff3.exon.fa-ORF  &
#length filter
seqkit seq -w0 Ma_all_isoforms.combined.gff3.exon.fa-ORF|awk '{if($0~/>/)printf $0"\t";else print $0}'|awk '{if(length($NF)>=100)print $0}' |sed 's/\t/\n/g' > Ma_all_isoforms.combined.gff3.exon.fa-ORF_100aa
#blastp
blastp -query Ma_NGS_scallop_merged.gtf.exon.fa-ORF_100aa \
	-db ~/1.bin/Database/uniprot_sprot.fa.blastdb \
	-max_target_seqs 5 \
	-outfmt "6 qacc qlen sacc slen pident mismatch evalue score qstart qend sstart send" \
	-evalue 1e-10 -num_threads 40 \
	> Ma_NGS_scallop_merged.gtf.exon.fa-ORF_100aa_blastp.outfmt6 &
#blastp uniq
awk -F '[_:]' '{print $0,$2"_"$3}' Ma_all_isoforms.combined.gff3.exon.fa-ORF_100aa_blastp.outfmt6|sort -k13,13 -k12,12nr|uniq -f12 > Ma_all_isoforms.combined.gff3.exon.fa-ORF_100aa_blastp.outfmt6-uniq &
awk '{print $1}' Ma_all_isoforms.combined.gff3.exon.fa-ORF_100aa_blastp.outfmt6-uniq > Ma_all_isoforms.combined.gff3.exon.fa-ORF_100aa_blastp.outfmt6-uniq-ID
awk -F '[:_]' '{print $2"_"$3}' Ma_all_isoforms.combined.gff3.exon.fa-ORF_100aa_blastp.outfmt6-uniq-ID > Ma_all_isoforms.combined.gff3.exon.fa-ORF_100aa_blastp.outfmt6-uniq-ID2
seqkit grep -f Ma_all_isoforms.combined.gff3.exon.fa-ORF_100aa_blastp.outfmt6-uniq-ID Ma_all_isoforms.combined.gff3.exon.fa-ORF_100aa  > Ma_all_isoforms.combined.gff3.exon.fa-ORF_100aa_blastp.outfmt6-uniq.faa
seqkit grep -f Ma_all_isoforms.combined.gff3.exon.fa-ORF_100aa_blastp.outfmt6-uniq-ID2 Ma_all_isoforms.combined.gff3.exon.fa > Ma_all_isoforms.combined.gff3.exon.fa-ORF_100aa-uniq.fa
#Get fake Transdecoder gff
python ORFfinder2transdecoder_gff3.py -cds_fa Ma_all_isoforms.combined.gff3.exon.fa-ORF_100aa-uniq.fa -ORFfinder_out Ma_all_isoforms.combined.gff3.exon.fa-ORF_100aa_blastp.outfmt6-uniq.faa -out_gff3 Ma_all_isoforms.combined.gff3.exon.fa-ORF.gff3
~/TransDecoder-TransDecoder-v5.5.0/util/gtf_to_alignment_gff3.pl Ma_all_isoforms.combined.gtf > Ma_all_isoforms.combined.gff3
~/TransDecoder-TransDecoder-v5.5.0/util/cdna_alignment_orf_to_genome_orf.pl Ma_all_isoforms.combined.gff3.exon.fa-ORF.gff3 Ma_all_isoforms.combined.gff3 Ma_all_isoforms.combined.gff3.exon.fa-ORF_100aa-uniq.fa > Ma_all_isoforms.combined.gff3.exon.fa-ORF-genome.gff3

##maker
nohup /opt/software/mpi/mpich/bin/mpiexec -n 50 ../bin/maker -base output_file_name -fix_nucleotides &

##Liftoff v1.6.1
liftoff -m minimap2 -p 50 -u Ma_txt -db CAU.gff_db -o Ma_Liftoff Ma.genome.fa CAU.fasta -polish
##EDTA
perl ~/1.bin/EDTA-1.9.6/EDTA.pl --anno 1 --genome genome.fa --species others --step all -t 20

#Repeat annotation
#cat EDTA Lib and Repeatmasker lib;
cat ~/1.bin/miniconda3/envs/EDTA/share/RepeatMasker/Libraries/RepeatMasker.lib >> TELib.fa
RepeatMasker -e ncbi -pa 40 -q -no_is -norna -nolow -div 40 -lib TElib.fa genome.fa



#Function annotation
#eggmapper
~/1.bin/miniconda2/envs/py37/bin/python ~/1.bin/eggnog-mapper-2.1.2/emapper.py -i Merge-maker-lift-ccs-ordered.gff.faa --output Merge-maker-lift-ccs -m diamond --cpu 80  &
awk -F '[\t]' '{print $1,$9}' Merge-maker-lift-ccs.emapper.annotations|grep -v ^\# > emapperID-mRNA

#Uniport
awk '{print $0,$1}' Merge-maker-lift-ccs-ordered.gff.faa-Uniport_protein-outfmt6 |uniq -f12 > Merge-maker-lift-ccs-ordered.gff.faa-Uniport_protein-outfmt6-Uniq
awk '{print $1,$3}' Merge-maker-lift-ccs-ordered.gff.faa.outfmt6-Uniq > UniportID-mRNA
sed 's/|/\t/g' UniportID-mRNA|awk '{print $3}'|sort|uniq  > UniportID-mRNA2
sed 's/|/\t/g' UniportID-mRNA > UniportID-mRNA3
#UniportID-GeneID.txt from are uniport webpage
sed '1d' UniportID-GeneID.txt|awk '{print $1,$2}'|sort|uniq > UniportID-GeneID.txt2
#join -1 1 -2 3 <(sort -k1,1 Uniport-mapping-list.txt2) <(sort -k3,3 UniportID-mRNA3)|awk '{print $3,$2}'  > Uniport-Anno
join -1 1 -2 3 <(awk '{print $0,$1}' BT-Uniport-mapping-list.txt2|sort -k1,1|uniq -f2|awk '{print $1,$2}') <(sort -k3,3 UniportID-mRNA3)|awk '{print $3,$2}'  > Uniport-Anno2



#CAU ID relationship
awk '{if($3=="gene")print $0}' ~/2.analysis/1.Annotation/0.Genome/CAU.gff > CAU-Gene-mRNA
awk '{print $9}' CAU-Gene-mRNA|awk -F '[=;]' '{print $2,$4}' > CAU-Gene2symble

awk '{print $0,$1}'  Merge-maker-lift-ccs-ordered.gff.faa-CAU_protein-outfmt6|uniq -f12 > Merge-maker-lift-ccs-ordered.gff.faa-CAU_protein-outfmt6-Uniq &
awk '{print $1,$3}' Merge-maker-lift-ccs-ordered.gff.faa-CAU_protein-outfmt6-Uniq > CAUID-mRNA
awk -F\. '{print $1"."$2}' CAUID-mRNA  > CAUID-mRNA2
join -1 1 -2 2 <(sort -k1,1 CAU-Gene2symble) <(sort -k2,2 CAUID-mRNA2)|awk '{print $3,$2}' > CAU-Anno 

#ZJU ID relationship
awk -F '[\t]' '{if($3=="CDS")print $9}' ~/0.data/0.Reference/0.Genome/Pek/ZJU1.0/GCF_015476345.1_ZJU1.0_genomic.gff |uniq > ZJU_CDS_Description
awk -F\; '{for(i=1;i<=NF;i++){if($i~/^ID/ || $i~/gene=/)printf $i"\t"};print ""}' ZJU_CDS_Description |sed 's/ID=cds-//g'|sed 's/gene=//g' > ZJU-proteion2symble 

awk '{print $0,$1}'  Merge-maker-lift-ccs-ordered.gff.faa-ZJU_protein-outfmt6|uniq -f12 > Merge-maker-lift-ccs-ordered.gff.faa-ZJU_protein-outfmt6-Uniq &
awk '{print $1,$3}' Merge-maker-lift-ccs-ordered.gff.faa-ZJU_protein-outfmt6-Uniq > ZJUID-mRNA
join -1 1 -2 2 <(sort -k1,1 ZJU-proteion2symble) <(sort -k2,2 ZJUID-mRNA)|awk '{print $3,$2}'|sort|uniq > ZJU-Anno

join -1 1 -2 1 -a1 <(sort All_ID) <(sort -k1,1 CAU-Anno)|awk '{if(NF<2){print $0,"missing"}else print $0}' > All_ID_CAU
join -1 1 -2 1 -a1 <(sort -k1,1 All_ID_CAU) <(sort -k1,1 ZJU-Anno)|awk '{if(NF<3){print $0,"missing"}else print $0}' > All_ID_CAU_ZJU
join -1 1 -2 1 -a1 <(sort -k1,1 All_ID_CAU_ZJU) <(sort -k1,1 emapperID-mRNA)|awk '{if(NF<4){print $0,"missing"}else print $0}' > All_ID_CAU_ZJU_emmaper
join -1 1 -2 1 -a1 <(sort -k1,1 All_ID_CAU_ZJU_emmaper) <(sort -k1,1 Uniport-Anno|awk '{if(NF<5){print $0,"missing"}else print $0}' > All_ID_CAU_ZJU_emmaper_Uniport

awk '{if(toupper($2)==toupper($3) && $2!="missing" && $2!~/^LOC/)print $1,$2}' All_ID_CAU_ZJU_emmaper_Uniport > Final-Ma-Anno-class1


awk '{if(toupper($2)==toupper($3) && $2=="missing" && $3=="missing") \
			{if($4!="-" && $4!="missing")print $1,$4;else if($5!="missing")print $1,$5;else print $1,"missing"}}' All_ID_CAU_ZJU_emmaper_Uniport > Final-Ma-Anno-class4

awk '{if(toupper($2)==toupper($3) && $2~/^LOC/) \
{if(toupper($4)==toupper($5) && $4!="missing")print $1,$4,$2,"emapper-uniport"; \
else if($4!="-" && $4!="missing")print $1,$4,$2,"emapper"; \
else if($5!="missing")print $1,$5,$2,"uniport";else print $1,$2}
}' All_ID_CAU_ZJU_emmaper_Uniport > Final-Ma-Anno-class2

awk '{if(toupper($2)!=toupper($3)) \
{if($2~/^LOC/ && $3~/^LOC/) \
{if($4!="-" && $4!="missing")print $1,$4;else if($5!="missing")print $1,$5;else print $1,$2;} \
else if((toupper($2)==toupper($4)||toupper($2)==toupper($5)) && $2!="missing")print $1,$2; \
else if((toupper($3)==toupper($4)||toupper($3)==toupper($5)) && $3!="missing")print $1,$3; \
else if(toupper($4)==toupper($5) && $4!="missing")print $1,$4; \
else if($2!="missing" && $2!~/^LOC/){print $1,$2} \
else if($3!="missing" && $3!~/^LOC/){print $1,$3} \
else if($3=="missing" && $2~/^LOC/){print $1,$2} \
else if($2=="missing" && $3~/^LOC/){print $1,$3} \
}}' \
 All_ID_CAU_ZJU_emmaper_Uniport > Final-Ma-Anno-class3

awk '{if(NR==FNR){if($1~/gene/ && $2!~/missing/ && $2!~/^LOC/){split($1,mRNAID,"-mRNA");geneID=mRNAID[1];Gene[geneID]=$2;if($3~/^LOC/)Gene2[geneID]=$3} \
else if($1~/mrna/ && $2!~/missing/ && $2!~/^LOC/){split($1,mRNAID,".mrna");geneID=mRNAID[1];Gene[geneID]=$2;if($3~/^LOC/)Gene2[geneID]=$3} \
else if($1~/TCONS/ && $2!~/missing/ && $2!~/^LOC/){split($1,mRNAID,".p");geneID=mRNAID[1];Gene[geneID]=$2;if($3~/^LOC/)Gene2[geneID]=$3}} \
else if(NR>FNR) \
	{if($1~/gene/ && ($2~/missing/ || $2~/^LOC/)){split($1,mRNAID,"-mRNA");geneID=mRNAID[1];if(geneID in GeneID && Gene2[geneID]==$2){print $1,GeneID[geneID]}else if(geneID in GeneID && $2~/missing/){print $1,GeneID[geneID]}else{print $1,$2}} \
	else if($1~/mrna/ && ($2~/missing/ || $2~/^LOC/)){split($1,mRNAID,".mrna");geneID=mRNAID[1];if(geneID in GeneID && Gene2[geneID]==$2){print $1,GeneID[geneID]}else if(geneID in GeneID && $2~/missing/){print $1,GeneID[geneID]}else{print $1,$2}}  \
	else if($1~/TCONS/ && ($2~/missing/ || $2~/^LOC/)){split($1,mRNAID,".p");geneID=mRNAID[1];if(geneID in GeneID && Gene2[geneID]==$2){print $1,GeneID[geneID]}else if(geneID in GeneID && $2~/missing/){print $1,GeneID[geneID]}else{print $1,$2}} \
	else{print $1,$2}}}' Final-Ma-Anno-class Final-Ma-Anno-class|uniq > 2Final-Ma-Anno-class

#Then use the mRNA to Gene scripts. all source can change to "Merge"
awk -v chr=0 -v source=0 -v pos1=0 -v pos2=0 -v len=0 '
{if(NR==FNR)Gene[$1]=$2;else{split($9,DES,";");gsub("ID=","",DES[1]);ID=DES[1];if(FNR==1){pos1=$4;pos2=$5;len=pos2-pos1+1;count=1;chr=$1;strand=$7;source=$2;mRNA[count]=$0;oldID=ID} \
else if(FNR!=1&&chr!=$1){ \
        if(1 in mRNAAAAA){printf chr"\t"sourceeeee"\tgene\t"posssss1"\t"posssss2"\t.\t"stranddddd"\t.\tID="Gene[oldID5]"\n";for(i=1;i<=counttttt+1;i++){print mRNAAAAA[i]};delete mRNAAAAA}; \
        if(1 in mRNAAAA){printf chr"\t"sourceeee"\tgene\t"possss1"\t"possss2"\t.\t"strandddd"\t.\tID="Gene[oldID4]"\n";for(i=1;i<=countttt+1;i++){print mRNAAAA[i]};delete mRNAAAA}; \
        if(1 in mRNAAA){printf chr"\t"sourceee"\tgene\t"posss1"\t"posss2"\t.\t"stranddd"\t.\tID="Gene[oldID3]"\n";for(i=1;i<=counttt+1;i++){print mRNAAA[i]};delete mRNAAA}; \
        if(1 in mRNAA){printf chr"\t"sourcee"\tgene\t"poss1"\t"poss2"\t.\t"strandd"\t.\tID="Gene[oldID2]"\n";for(i=1;i<=countt+1;i++){print mRNAA[i]};delete mRNAA};\
        printf chr"\t"source"\tgene\t"pos1"\t"pos2"\t.\t"strand"\t.\tID="Gene[oldID]"\n";for(i=1;i<=count+1;i++){print mRNA[i]};\
        pos1=$4;pos2=$5;len=pos2-pos1+1;count=1;chr=$1;strand=$7;source=$2;delete mRNA;mRNA[count]=$0;oldID=ID} \
else{ \
        if($4<=(pos2-len*0.1) && $5>=(pos1+len*0.1)) \
                {if((ID in Gene && oldID in Gene && Gene[ID]!=Gene[oldID])||strand!=$7){ \
                        if($4<=(poss2-lenn*0.1) && $5>=(poss1+lenn*0.1) && ID in Gene && oldID2 in Gene && Gene[ID]==Gene[oldID2] && strandd==$7){countt+=1;mRNAA[countt]=$0;if($4<=poss1)poss1=$4;if($5>=poss2)poss2=$5;lenn=poss2-poss1+1;chr=$1;source=$2} \
                        else if($4<=(posss2-lennn*0.1) && $5>=(posss1+lennn*0.1) && ID in Gene && oldID3 in Gene && Gene[ID]==Gene[oldID3] && stranddd==$7){counttt+=1;mRNAAA[counttt]=$0;if($4<=posss1)posss1=$4;if($5>=posss2)posss2=$5;lennn=posss2-posss1+1;chr=$1;source=$2} \
                        else if($4<=(possss2-lennnn*0.1) && $5>=(possss1+lennnn*0.1) && ID in Gene && oldID4 in Gene && Gene[ID]==Gene[oldID4] && strandddd==$7){countttt+=1;mRNAAAA[countttt]=$0;if($4<=possss1)possss1=$4;if($5>=possss2)possss2=$5;lennnn=possss2-possss1+1;chr=$1;source=$2} \
                        else if($4<=(posssss2-lennnnn*0.1) && $5>=(posssss1+lennnnn*0.1) && ID in Gene && oldID5 in Gene && Gene[ID]==Gene[oldID5] && stranddddd==$7){counttttt+=1;mRNAAAAA[counttttt]=$0;if($4<=posssss1)posssss1=$4;if($5>=posssss2)posssss2=$5;lennnnn=posssss2-posssss1+1;chr=$1;source=$2} \
                        else{if(1 in mRNAAAAA){printf chr"\t"source"\tgene\t"posssss1"\t"posssss2"\t.\t"stranddddd"\t.\tID="Gene[oldID5]"\n"; \
                                for(i=1;i<=counttttt+1;i++){print mRNAAAAA[i]};delete mRNAAAAA};\
                                for(i in mRNAAAA)mRNAAAAA[i]=mRNAAAA[i];oldID5=oldID4;posssss1=possss1;posssss2=possss2;lennnnn=lennnn;counttttt=countttt;stranddddd=strandddd;sourceeeee=sourceeee;delete mRNAAAA;\
                                for(i in mRNAAA)mRNAAAA[i]=mRNAAA[i];oldID4=oldID3;possss1=posss1;possss2=posss2;lennnn=lennn;countttt=counttt;strandddd=stranddd;sourceeee=sourceee;delete mRNAAA;\
                                for(i in mRNAA)mRNAAA[i]=mRNAA[i];oldID3=oldID2;posss1=poss1;posss2=poss2;lennn=lenn;counttt=countt;stranddd=strandd;sourceee=sourcee;delete mRNAA;\
                                for(i in mRNA)mRNAA[i]=mRNA[i];oldID2=oldID;poss1=pos1;poss2=pos2;lenn=len;countt=count;strandd=strand;sourcee=source;delete mRNA;\
                                pos1=$4;pos2=$5;len=$5-$4;count=1;chr=$1;strand=$7;source=$2;delete mRNA;mRNA[count]=$0;oldID=ID}}\
                else{if(ID in Gene)oldID=ID;count+=1;mRNA[count]=$0;if($4<=pos1)pos1=$4;if($5>=pos2)pos2=$5;len=pos2-pos1+1;chr=$1;source=$2}}\
        else if($4<=(poss2-lenn*0.1) && $5>=(poss1+lenn*0.1) && ID in Gene && oldID2 in Gene && Gene[ID]==Gene[oldID2] && strandd==$7){countt+=1;mRNAA[countt]=$0;if($4<=poss1)poss1=$4;if($5>=poss2)poss2=$5;lenn=poss2-poss1+1;chr=$1;sourcee=$2}\
        else if($4<=(posss2-lennn*0.1) && $5>=(posss1+lennn*0.1) && ID in Gene && oldID3 in Gene && Gene[ID]==Gene[oldID3] && stranddd==$7){counttt+=1;mRNAAA[counttt]=$0;if($4<=posss1)posss1=$4;if($5>=posss2)posss2=$5;lennn=posss2-posss1+1;chr=$1;sourceee=$2}\
        else if($4<=(possss2-lennnn*0.1) && $5>=(possss1+lennnn*0.1) && ID in Gene && oldID4 in Gene && Gene[ID]==Gene[oldID4] && strandddd==$7){countttt+=1;mRNAAAA[countttt]=$0;if($4<=possss1)possss1=$4;if($5>=possss2)possss2=$5;lessnn=possss2-possss1+1;chr=$1;sourceeee=$2}\
        else if($4<=(posssss2-lennnnn*0.1) && $5>=(posssss1+lennnnn*0.1) && ID in Gene && oldID5 in Gene && Gene[ID]==Gene[oldID5] && stranddddd==$7){counttttt+=1;mRNAAAAA[counttttt]=$0;if($4<=posssss1)posssss1=$4;if($5>=posssss2)posssss2=$5;lennnnn=posssss2-posssss1+1;chr=$1;sourceeeee=$2}\
        else if(!($4<=(pos2-len*0.1) && $5>=(pos1+len*0.1))){ \
                if(ID in Gene && oldID in Gene && Gene[ID]==Gene[oldID] && strand==$7){count+=1;mRNA[count]=$0;if($4<=pos1)pos1=$4;if($5>=pos2)pos2=$5;len=pos2-pos1+1;chr=$1;source=$2} \
                else{if(1 in mRNAAAAA){printf chr"\t"sourceeeee"\tgene\t"posssss1"\t"posssss2"\t.\t"stranddddd"\t.\tID="Gene[oldID5]"\n";for(i=1;i<=counttttt+1;i++){print mRNAAAAA[i]};delete mRNAAAAA};\
                        for(i in mRNAAAA)mRNAAAAA[i]=mRNAAAA[i];oldID5=oldID4;posssss1=possss1;posssss2=possss2;lennnnn=lennnn;counttttt=countttt;stranddddd=strandddd;sourceeeee=sourceeee;delete mRNAAAA;\
                        for(i in mRNAAA)mRNAAAA[i]=mRNAAA[i];oldID4=oldID3;possss1=posss1;possss2=posss2;lennnn=lennn;countttt=counttt;strandddd=stranddd;sourceeee=sourceee;delete mRNAAA;\
                        for(i in mRNAA)mRNAAA[i]=mRNAA[i];oldID3=oldID2;posss1=poss1;posss2=poss2;lennn=lenn;counttt=countt;stranddd=strandd;sourceee=sourcee;delete mRNAA;\
                        for(i in mRNA)mRNAA[i]=mRNA[i];oldID2=oldID;poss1=pos1;poss2=pos2;lenn=len;countt=count;strandd=strand;sourcee=source;delete mRNA;\
                        pos1=$4;pos2=$5;len=$5-$4;count=1;chr=$1;strand=$7;source=$2;delete mRNA;mRNA[count]=$0;oldID=ID}}}}} \
END{ \
        if(1 in mRNAAAAA){printf chr"\t"sourceeeee"\tgene\t"posssss1"\t"posssss2"\t.\t"stranddddd"\t.\tID="Gene[oldID5]"\n";for(i=1;i<=counttttt+1;i++){print mRNAAAAA[i]};delete mRNAAAAA}; \
        if(1 in mRNAAAA){printf chr"\t"sourceeee"\tgene\t"possss1"\t"possss2"\t.\t"strandddd"\t.\tID="Gene[oldID4]"\n";for(i=1;i<=countttt+1;i++){print mRNAAAA[i]};delete mRNAAAA}; \
        if(1 in mRNAAA){printf chr"\t"sourceee"\tgene\t"posss1"\t"posss2"\t.\t"stranddd"\t.\tID="Gene[oldID3]"\n";for(i=1;i<=counttt+1;i++){print mRNAAA[i]};delete mRNAAA}; \
        if(1 in mRNAA){printf chr"\t"sourcee"\tgene\t"poss1"\t"poss2"\t.\t"strandd"\t.\tID="Gene[oldID2]"\n";for(i=1;i<=countt+1;i++){print mRNAA[i]};delete mRNAA};\
        printf chr"\t"source"\tgene\t"pos1"\t"pos2"\t.\t"strand"\t.\tID="Gene[oldID]"\n";for(i=1;i<=count+1;i++){print mRNA[i]};
        }' BlastP/2Final-Ma-Anno-class Merge-maker-lift-ccs.gff-sort-mRNA > mRNA2gene.out

#Then added exon and cds info.
#no change ID
awk -v old=0 -v oldID=0 '{if(NR==FNR){ \
				if($3=="exon"){split($9,DES,";");gsub("Parent=","",DES[1]);Parent=DES[1];if(Parent==oldID && old=="exon"){count+=1}else count=1;GFF["exon"][Parent][count]=$0;oldID=Parent;old=exon} \
				if($3=="CDS"){split($9,DES,";");gsub("Parent=","",DES[1]);Parent=DES[1];if(Parent==oldID && old=="CDS"){count+=1}else count=1;GFF["CDS"][Parent][count]=$0;oldID=Parent;old=CDS} \
						} \
				else{ \
					if($3=="gene")print $0; \
					else if($3=="mRNA"){split($9,DES,";");gsub("ID=","",DES[1]);Parent=DES[1];print $0; \
								if(Parent in GFF["exon"]){for(i=1;i<=length(GFF["exon"][Parent]);i++)print GFF["exon"][Parent][i]}; \
								if(Parent in GFF["CDS"]){for(i=1;i<=length(GFF["CDS"][Parent]);i++)print GFF["CDS"][Parent][i]}} \
					else{print $0}}}' ../Merge-maker-lift-ccs.gff Ma-Gene-pos-ID > tmp




#added exon and cds and rename
awk -v oldID=0 -v genecount=0 '{if(NR==FNR){ \
					if($3=="exon"){split($9,DES,";");gsub("Parent=","",DES[1]);Parent=DES[1];if(Parent==oldID &&  old=="exon"){count+=1}else count=1;GFF["exon"][Parent][count]=$0;oldID=Parent;old="exon"} \
					if($3=="CDS"){split($9,DES,";");gsub("Parent=","",DES[1]);Parent=DES[1];if(Parent==oldID &&  old=="CDS"){count+=1}else count=1;GFF["CDS"][Parent][count]=$0;oldID=Parent;old="CDS"} \
							} \
				else{ \
					if($3=="gene"){genecount+=1;tmp="0000";tmpp=tmp"0"genecount;genecountt=substr(tmpp,length(tmpp)-4);split($NF,DES,";");gsub("ID=","",DES[1]);print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tID=Mallard-gene"genecountt";gene="DES[1]";gene_biotype=protein_coding";mRNAID=0} \
					else if($3=="mRNA"){mRNAID+=1;split($9,DES,";");gsub("ID=","",DES[1]);Parent=DES[1];print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tID=Mallard-gene"genecountt".mrna"mRNAID";Parent=Mallard-gene"genecountt; \
										if(Parent in GFF["exon"]){for(i=1;i<=length(GFF["exon"][Parent]);i++){exon=GFF["exon"][Parent][i];split(exon,exons,"\t");print exons[1]"\t"exons[2]"\t"exons[3]"\t"exons[4]"\t"exons[5]"\t"exons[6]"\t"exons[7]"\t"exons[8]"\tID=exon-"i";Parent=Mallard-gene"genecountt".mrna"mRNAID}}; \
										if(Parent in GFF["CDS"]){for(i=1;i<=length(GFF["CDS"][Parent]);i++){CDS=GFF["CDS"][Parent][i];split(CDS,CDSs,"\t");print CDSs[1]"\t"CDSs[2]"\t"CDSs[3]"\t"CDSs[4]"\t"CDSs[5]"\t"CDSs[6]"\t"CDSs[7]"\t"CDSs[8]"\tID=CDS-"i";Parent=Mallard-gene"genecountt".mrna"mRNAID}}} \
					else{print $0}}}' Merge-maker-lift-ccs.gff mRNA2gene.out > add-exon-cds.out


#added UTR
awk -v pos1=0 -v pos2=0 '{if(NR==FNR && $3=="CDS") \
							{split($9,DES,";");gsub("Parent=","",DES[2]);mRNAID=DES[2]; \
								if(mRNAID==mRNAIDold){CDS[mRNAID]["pos2"]=$5}else{CDS[mRNAID]["pos1"]=$4;CDS[mRNAID]["pos2"]=$5};mRNAIDold=mRNAID} \
						else{if(NR>FNR){if($3=="exon"){print $0;split($0,DES,";");gsub("Parent=","",DES[2]);mRNAID=DES[2]; \
											if($5<CDS[mRNAID]["pos1"]){mRNAUTRS=mRNAUTRS"\n"$1"\t"$2"\tUTR\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tID=UTR-"count";Parent="mRNAID;count+=1} \
											else if($4<CDS[mRNAID]["pos1"]&&$5>=CDS[mRNAID]["pos1"]){mRNAUTRS=mRNAUTRS"\n"$1"\t"$2"\tUTR\t"$4"\t"CDS[mRNAID]["pos1"]-1"\t"$6"\t"$7"\t"$8"\tID=UTR-"count";Parent="mRNAID;count+=1} \
											else if($4<=CDS[mRNAID]["pos2"]&&$5>CDS[mRNAID]["pos2"]){mRNAUTRL=mRNAUTRL"\n"$1"\t"$2"\tUTR\t"CDS[mRNAID]["pos2"]+1"\t"$5"\t"$6"\t"$7"\t"$8"\tID=UTR-"count";Parent="mRNAID;count+=1} \
											else if($4>CDS[mRNAID]["pos2"]){mRNAUTRL=mRNAUTRL"\n"$1"\t"$2"\tUTR\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tID=UTR-"count";Parent="mRNAID;count+=1};} \
											else if($3=="mRNA" && oldtype!="gene"){if(strand=="+"){gsub("UTR","3-UTR",mRNAUTRL)}else{gsub("UTR","5-UTR",mRNAUTRL)};print mRNAUTRL;oldtype="mRNA";count=1};mRNAUTRL=NULL;print $0} \
											else if($3=="mRNA" && oldtype=="gene"){count=1;print $0} \
											else if($3=="CDS" && oldtype!="CDS"){if(strand=="+"){gsub("UTR","5-UTR",mRNAUTRS)}else{gsub("UTR","3-UTR",mRNAUTRS)};print mRNAUTRS;mRNAUTRS=NULL;print $0;oldtype="CDS"} \
											else if($3=="CDS" && oldtype=="CDS"){print $0} \
											else if($3=="gene"){if(strand=="+"){gsub("UTR","3-UTR",mRNAUTRL)}else{gsub("UTR","5-UTR",mRNAUTRL);print mRNAUTRL;mRNAUTRL=NULL;print $0;oldtype="gene";strand=$7}}}}}' add-exon-cds.out add-exon-cds.out|sed '/^$/d' > add-UTR.out


#then add EDTA Repeats and repeatrunner results
grep -v match_part Maker-repeat.gff > Maker-repeat.gff2
ln -s ~/2.analysis/1.Annotation/2.EDTA-RepeatMasker/2.TELib_SoftMask/Ma/Ma03-Final_PLUS_Unplaced.fa.gff ./EDTA-repeat.gff

bedtools intersect -a Maker-repeat.gff2 -b EDTA-repeat.gff -wo > tmp-Overlap.gff
awk -F '[\t]' '{print $5-$4+1,$14-$13+1,$NF}' tmp-Overlap.gff |awk '{print $3/$1,"type1";print $3/$2,"type2"}' > overlap_distribution

###############################
bedtools intersect -a Maker-repeat.gff2 -b EDTA-repeat.gff -F 0.5 -f 0.5 -v > Maker-noredundent.gff
grep runner Maker-noredundent.gff > Maker-noredundent-repeatrunner.gff
cat EDTA-repeat.gff Maker-noredundent-repeatrunner.gff > Final-Repeat.gff

#mark the repeat related gene
bedtools intersect -a add-UTR.out-gene.gff -b Final-Repeat.gff -wo > tmp-overlap-gene.gff
awk -F '[\t]' '{print $5-$4+1,$14-$13+1,$NF}' tmp-overlap-gene.gff|awk '{print $3/$1,"type1";print $3/$2,"type2"}' > overlap-gene-distribution
Rscript overlap-distribution.R
mv Rplots.pdf overlap-gene-distribution.pdf

bedtools intersect -a add-UTR.out-gene.gff -b Final-Repeat.gff -f 0.5 -wo > tmp-overlap-gene0.5.gff
awk -F '[\t]' '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' tmp-overlap-gene0.5.gff|sed 's/Merge/Merge-Repeat-Relate/g'|uniq > Repeat-Relate-gene.gff
awk '{if(NR==FNR){ID=$1"-"$4"-"$5;Gene[ID]=$0}else{if($3=="gene"){ID=$1"-"$4"-"$5;if(ID in Gene){gsub("Merge","Merge-Repeat-Relate",$2);print $0;}else print $0;}else print $0}}' Repeat-Relate-gene.gff add-UTR.out > new-add-UTR.out
cat new-add-UTR.out Final-Repeat.gff >> add-UTR-Repeat.gff

awk -v note=0 '{if($3~/lnc/||note=="lnc"){if($3!~/gene/){print $0;note="lnc"}else if($3~/gene/)note=0}}' ~/Ma_Liftoff > Ma_ZJU_lnc_Liftoff
awk '{if($3~/lnc/)print $0}' Ma_ZJU_lnc_Liftoff > Ma_ZJU_lnc_Liftoff-lnc.gff
#mark the lnc related gene
bedtools intersect -a add-UTR.out-gene.gff -b Ma_ZJU_lnc_Liftoff-lnc.gff -wo > tmp-overlap-lnc-gene.gff
awk -F '[\t]' '{print $5-$4+1,$14-$13+1,$NF}' tmp-overlap-lnc-gene.gff|awk '{print $3/$1,"type1";print $3/$2,"type2"}' > overlap-lnc-gene-distribution
Rscript overlap-distribution.R
mv Rplots.pdf overlap-lnc-gene-distribution.pdf

bedtools intersect -a add-UTR.out-gene.gff -b Ma_ZJU_lnc_Liftoff-lnc.gff -f 0.5 -wo > tmp-overlap-lnc-gene0.5.gff
awk -F '[\t]' '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' tmp-overlap-lnc-gene0.5.gff|sed 's/Merge/Merge-Lnc-Relate/g'|uniq > Lnc-Relate-gene.gff
awk '{OFS="\t"}{if(NR==FNR){ID=$1"-"$4"-"$5;Gene[ID]=$0}else{if($3=="gene"){ID=$1"-"$4"-"$5;if(ID in Gene){gsub("Merge","Merge-Lnc-Relate",$2);print $0;}else print $0}else print $0}}' Lnc-Relate-gene.gff new-add-UTR.out > new-new-add-UTR.out

sort -suV Final-Repeat.gff|grep -v ^\# > 2Final-Repeat.gff 

#merge
cat new-new-add-UTR.out 2Final-Repeat.gff ../3.Lift-Lnc/Ma_ZJU_lnc_Liftoff > Ma-Final.gff  &





