#LTR sharing
#1# get seq & mappinng
for i in `ls *All*genome*ID`; do awk '{OFS="\t"}{if($1~/Chr/ && $3=="repeat_region"){print $1,$4-10000,$5+10000 > "'"$i"'-LTR-Up-Down10k.bed"}}' ../../${i}*raw/LTR/*mod.LTR.intact.gff3 ; done
for i in `ls *All*genome*ID`; do seqkit subseq --bed $i\-LTR-Up-Down10k.bed ../../../${i}*.fa > ${i}-LTR-Up-Down10k.fa & done
for DB in `ls *All*genome*ID`
do
{
    mkdir $DB-genome2 && cd $DB-genome2
    for Query in `ls *All*genome*ID`
    do 
        if [ $Query != $DB ];then blastn -query ../${Query}-LTR-Up-Down10k.fa -db ../${DB}-genome -evalue 1e-10 -outfmt "6 qacc qlen sacc slen pident mismatch evalue score qstart qend sstart send" -out ${Query}-Up-map-to-${DB} -num_threads 20
        fi
    done
cd /home/zhangzhou209/1.bin/chromeister-1.5.a/New/Final-ALL/1.EDTA/4.Intact-LTR/2.LTR-Overlap_Up-Down-10K/
} &
done

#2# results process
for i in `ls *`
do
    echo $i
    awk -v ID=0 '{if($1!=ID){print $0;ID=$1}}' $i > ${i}2
    awk '{if($(NF-4)>11000 && $5>95)print $0}' ${i}2 |wc
done
#3
#one-one share
for i in `ls *All*genome*ID`
do
    for j in `ls *All*genome*ID`
    do
        cat ${i}*genome2/${j}*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort|uniq|grep -v HiC|wc
    done
done
#nonSpecific

for i in `ls *All*genome*ID`
do
cat *genome2/${i}*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort|uniq|grep -v HiC|wc
done

#specific shared
for i in `ls *All*genome*ID`
do
cat *genome2/${i}*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort|uniq -c|awk '{if($1==1)print $2}' > ${i}-uniq-tmp
done

for i in `ls *All*genome*ID`
do
    cat ${i}*genome2/CAU*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort > tmp-uniq2
    comm -12 tmp-uniq2 <(grep -v HiC CAU-uniq-tmp)|wc
    cat ${i}*genome2/SX*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort > tmp-uniq2
    comm -12 tmp-uniq2 SX-uniq-tmp|wc
    cat ${i}*genome2/Ma*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort > tmp-uniq2
    comm -12 tmp-uniq2 Ma-uniq-tmp|wc
    cat ${i}*genome2/SB*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort > tmp-uniq2
    comm -12 tmp-uniq2 SB-uniq-tmp|wc
    cat ${i}*genome2/ZW*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort > tmp-uniq2
    comm -12 tmp-uniq2 ZW-uniq-tmp|wc
    cat ${i}*genome2/LC*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort > tmp-uniq2
    comm -12 tmp-uniq2 LC-uniq-tmp|wc
    cat ${i}*genome2/LW*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort > tmp-uniq2
    comm -12 tmp-uniq2 LW-uniq-tmp|wc
    cat ${i}*genome2/PZ*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort > tmp-uniq2
    comm -12 tmp-uniq2 PZ-uniq-tmp|wc
    cat ${i}*genome2/HL*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort > tmp-uniq2
    comm -12 tmp-uniq2 HL-uniq-tmp|wc
done