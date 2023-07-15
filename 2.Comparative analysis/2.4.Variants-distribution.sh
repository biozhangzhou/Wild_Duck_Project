
#cLTR ALL-LTR ALL-Repeat 

for i in `ls *ALL*Genome*ID`
do 
bedtools closest -a <(sort -k1,1 -k2,2n ${i}-cLTR.bed ) -b ${i}-gene.bed -D b > ${i}-gene-cLTR-Closet
bedtools closest -a <(sort -k1,1 -k2,2n ${i}-ALL-LTR.bed ) -b ${i}-gene.bed -D b > ${i}-gene-ALL-LTR-Closet
bedtools closest -a <(sort -k1,1 -k2,2n ${i}-ALL-Repeat.bed ) -b ${i}-gene.bed -D b > ${i}-gene-ALL-RepeatCloset
done


#Population
for i in `ls *ALL*Genome*ID`
do 
bedtools closest -a <(sort -k1,1 -k2,2n ../8.Resequencing/${i}.map.bed) -b CAU-gene.bed -D b > ${i}-Resequence-Closet &
done


for i in `ls *ALL*Genome*ID`
do 
bedtools closest -a <(sort -k1,1 -k2,2n ${i}-INDEL.map.bed|grep -v chrZ) -b CAU-gene.bed -D b > ${i}-INDEL-Resequence-Closet &
done

#NGS
for i in `ls *ALL*Genome*ID`
do 
{
bedtools closest -a  <(sort -k1,1 -k2,2n 5.NGS/${i}*SNP.vcf|grep -v chrZ) -b ${i}-gene.bed -D b > ${i}-SNP-Closet
bedtools closest -a  <(sort -k1,1 -k2,2n 5.NGS/${i}*INDEL.vcf|grep -v chrZ) -b ${i}-gene.bed -D b > ${i}-INDEL-Closet
} &
done

#TGS
for i in `ls *ALL*Genome*ID`
do 
cp 6.TGS/${i}-TGS2.VCF tmp.bed
bedtools closest -a <(sort -k1,1 -k2,2n tmp.bed|grep -v chrZ) -b ${i}-gene.bed -D b > ${i}-TGS-Closet
done

#PAN
for j in snp ins del
do
for i in `ls *ALL*Genome*ID`
do
        Info=7.PAN/0.New-ID/${i}-ID
        grep -v ^# 7.PAN/${i}-${j}.bed*|sed 's/^'$i'-//g'|sed 's/ /\t/g' > ${i}-${j}-tmp.bed
        BED=${i}-${j}-tmp.bed
        bedtools closest -a <(sort -k1,1 -k2,2n 2${i}-${j}-tmp.bed) -b ${i}-gene.bed -D b > ${j}-${i}-PAN-Closet
} &
done
done

#STAT 未包含0
#Pop
for i in `ls *ALL*Genome*ID`; do awk '{if($NF!=0 && $(NF-1)!=".")print "'"$i"' Population-SNP",$NF}' ${i}-Resequence-Closet >> ALL-stat; done 
for i in `ls *ALL*Genome*ID`; do awk '{if($NF!=0 && $(NF-1)!=".")print "'"$i"' Population-INDEL",$NF}' ${i}-INDEL-Resequence-Closet >> ALL-stat; done 
#PAN
for j in snp ins del; do for i in `ls *ALL*Genome*ID`; do awk '{if($NF!=0 && $(NF-1)!=".")print "'"$i"' PAN-'"$j"'",$NF}' ${j}-${i}-PAN-Closet >> ALL-stat; done; done

#TGS
for i in `ls *ALL*Genome*ID`
do
	awk '{if($NF!=0 && $(NF-1)!=".")print "'"$i"' TGS-SV",$NF}' ${i}-TGS-Closet >> ALL-stat
done
#NGS
for j in SNP INDEL; do for i in `ls *ALL*Genome*ID`; do awk '{if($NF!=0 && $(NF-1)!=".")print "'"$i"' NGS-'"$j"'",$NF}' ${i}-${j}-Closet >> ALL-stat; done; done
#Repeat
for i in ALL-Repeat ALL-LTR cLTR
do
for i in `ls *ALL*Genome*ID`
do
	awk '{if($NF!=0 && $(NF-1)!=".")print "'"$j"' '"$i"'",$NF}' Merge-Gene/${j}-gene-${i}*Closet >> ALL-stat
done
done
#HCE
for i in Duck Home; do for j in Con Acc; do for f in 20 100; do awk '{if($NF!=0 && $(NF-1)!=".")print "'"$i"' Type2 '"$j"' '"$f"'bp",$NF}' CAU-gene-Type2-${i}-${j}-${f}.bed-Closet >> HCE-stat; done; done; done
for i in Duck Home; do for f in 20 100; do awk '{if($NF!=0 && $(NF-1)!=".")print "'"$i"' Type1 Con '"$f"'bp",$NF}' CAU-gene-Type1-${i}-Con-${f}bp.bed-Closet >> HCE-stat; done; done

##Plotting
library(ggplot2)
data <- read.table("ALL-stat")
colnames(data) <- c("Species","Type","Distance")
data$Species <- gsub("CAU","Pek",data$Species)
data$Species <- gsub("SX","IND",data$Species)
data$Species <- factor(data$Species,ordered = TRUE,level=c("Pek","IND","Ma","SB","ZW","LC","LW","CJ","PZ","HL","G"))
data$Distance <- data$Distance/1000
data$Type <- factor(data$Type,ordered = TRUE,level=c("ALL-Repeat","ALL-LTR","cLTR","NGS-SNP","NGS-INDEL","TGS-SV","PAN-snp","PAN-ins","PAN-del","Population-SNP","Population-INDEL"))
#整体
pdf("ALL.pdf",height=7,width=7)
ggplot(data,aes(x=Distance,color=Species))+geom_density()+facet_wrap(~Type,ncol=3)+theme_bw()
dev.off()
pdf("ALL2.pdf",height=7,width=7)
ggplot(data,aes(x=Distance,color=Species))+geom_density()+facet_wrap(~Type,ncol=3)+theme_bw()+scale_color_brewer(palette = 'Set1')
dev.off()
pdf("ALL-100-100.pdf",height=7,width=7)
ggplot(data,aes(x=Distance,color=Species))+geom_density()+facet_wrap(~Type,ncol=3)+theme_bw()+xlim(-100,100)
dev.off()
pdf("ALL2-100-100.pdf",height=7,width=7)
ggplot(data,aes(x=Distance,color=Species))+geom_density()+facet_wrap(~Type,ncol=3)+theme_bw()+scale_color_brewer(palette = 'Set1')+xlim(-100,100)
dev.off()
pdf("ALL-0-100.pdf",height=7,width=7)
ggplot(data,aes(x=Distance,color=Species))+geom_density()+facet_wrap(~Type,ncol=3)+theme_bw()+xlim(0,100)
dev.off()
pdf("ALL2-0-100.pdf",height=7,width=7)
ggplot(data,aes(x=Distance,color=Species))+geom_density()+facet_wrap(~Type,ncol=3)+theme_bw()+scale_color_brewer(palette = 'Set1')+xlim(0,100)
dev.off()
pdf("ALL-100-0.pdf",height=7,width=7)
ggplot(data,aes(x=Distance,color=Species))+geom_density()+facet_wrap(~Type,ncol=3)+theme_bw()+xlim(-100,0)
dev.off()
pdf("ALL2-100-0.pdf",height=7,width=7)
ggplot(data,aes(x=Distance,color=Species))+geom_density()+facet_wrap(~Type,ncol=3)+theme_bw()+scale_color_brewer(palette = 'Set1')+xlim(-100,0)
dev.off()
#部分
pdf("Repeat-100-100.pdf",height=7,width=7)
ggplot(data[grep("Repeat|LTR",data$Type),],aes(x=Distance,color=Species))+geom_density()+facet_wrap(~Type,ncol=3)+theme_bw()+scale_color_brewer(palette = 'Set1')+xlim(-100,100)
dev.off()
pdf("Repeat2-100-100.pdf",height=7,width=7)
ggplot(data[grep("Repeat|LTR",data$Type),],aes(x=Distance,color=Species))+geom_density()+facet_wrap(~Type,ncol=3)+theme_bw()+xlim(-100,100)
dev.off()
pdf("PAN-100-100.pdf",height=7,width=7)
ggplot(data[grep("PAN",data$Type),],aes(x=Distance,color=Species))+geom_density()+facet_wrap(~Type,ncol=3)+theme_bw()+scale_color_brewer(palette = 'Set1')+xlim(-100,100)
dev.off()
pdf("PAN2-100-100.pdf",height=7,width=7)
ggplot(data[grep("PAN",data$Type),],aes(x=Distance,color=Species))+geom_density()+facet_wrap(~Type,ncol=3)+theme_bw()+xlim(-100,100)
dev.off()
pdf("NGS-TGS-100-100.pdf",height=7,width=7)
ggplot(data[grep("NGS|TGS",data$Type),],aes(x=Distance,color=Species))+geom_density()+facet_wrap(~Type,ncol=3)+theme_bw()+scale_color_brewer(palette = 'Set1')+xlim(-100,100)
dev.off()
pdf("NGS-TGS2-100-100.pdf",height=7,width=7)
ggplot(data[grep("NGS|TGS",data$Type),],aes(x=Distance,color=Species))+geom_density()+facet_wrap(~Type,ncol=3)+theme_bw()+xlim(-100,100)
dev.off()
pdf("Pop-100-100.pdf",height=7,width=7)
ggplot(data[grep("Pop",data$Type),],aes(x=Distance,color=Species))+geom_density()+facet_wrap(~Type,ncol=3)+theme_bw()+scale_color_brewer(palette = 'Set1')+xlim(-100,100)
dev.off()
pdf("Pop2-100-100.pdf",height=7,width=7)
ggplot(data[grep("Pop",data$Type),],aes(x=Distance,color=Species))+geom_density()+facet_wrap(~Type,ncol=3)+theme_bw()+xlim(-100,100)
dev.off()



data2 <- read.table("HCE-stat")
colnames(data2) <- c("Node","Type2","Type","Length","Distance")
data2$Distance <- data2$Distance/1000
pdf("HCE.pdf",height=7,width=7)
ggplot(data2,aes(x=Distance,color=Length))+geom_density()+facet_grid(Node~Type2+Type)+theme_bw()+scale_color_brewer(palette = 'Set1')
dev.off()
pdf("HCE2.pdf",height=7,width=7)
ggplot(data2,aes(x=Distance,color=Length))+geom_density()+facet_grid(Node~Type2+Type)+theme_bw()
dev.off()
pdf("HCE-100-100.pdf",height=7,width=7)
ggplot(data2,aes(x=Distance,color=Length))+geom_density()+facet_grid(Node~Type2+Type)+theme_bw()+scale_color_brewer(palette = 'Set1')+xlim(-100,100)
dev.off()
pdf("HCE2-100-100.pdf",height=7,width=7)
ggplot(data2,aes(x=Distance,color=Length))+geom_density()+facet_grid(Node~Type2+Type)+theme_bw()+xlim(-100,100)
dev.off()
pdf("HCE-0-100.pdf",height=7,width=7)
ggplot(data2,aes(x=Distance,color=Length))+geom_density()+facet_grid(Node~Type2+Type)+theme_bw()+scale_color_brewer(palette = 'Set1')+xlim(0,100)
dev.off()
pdf("HCE2-0-100.pdf",height=7,width=7)
ggplot(data2,aes(x=Distance,color=Length))+geom_density()+facet_grid(Node~Type2+Type)+theme_bw()+xlim(0,100)
dev.off()
pdf("HCE-100-0.pdf",height=7,width=7)
ggplot(data2,aes(x=Distance,color=Length))+geom_density()+facet_grid(Node~Type2+Type)+theme_bw()+scale_color_brewer(palette = 'Set1')+xlim(-100,0)
dev.off()
pdf("HCE2-100-0.pdf",height=7,width=7)
ggplot(data2,aes(x=Distance,color=Length))+geom_density()+facet_grid(Node~Type2+Type)+theme_bw()+xlim(-100,0)
dev.off()







