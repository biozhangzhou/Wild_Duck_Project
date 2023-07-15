#classify 478 CAU-cLTR into three types
bedtools subtract -f 1 -a <(awk '{print "chr"$2"\t"$3"\t"$4"\t"$1"-"$5}' CAU-LTR-Sim-ID|sed 's/chrchr/chr/g'|sort -k1,1 -k2,2n) -b <(cat *genome2/CAU*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort|uniq|sed 's/-/\t/g'|sed 's/_/\t/g'|sed 's/:.//g'|awk '{if($0~/HiC/)print $1"_"$2"_"$3"\t"$4"\t"$5;else print $0}'|sort -k1,1 -k2,2n)| \
awk '{print $0,"Pekin-Specific"}' > Diff-Type-LTR/2022-Final-Three-Class 

bedtools intersect -f 0.9 -a <(awk '{print "chr"$2"\t"$3"\t"$4"\t"$1"-"$5}' CAU-LTR-Sim-ID|sed 's/chrchr/chr/g'|sort -k1,1 -k2,2n) -b <(comm -12 <(cat SX*genome2/CAU*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort) CAU-uniq-tmp|grep -v HiC|sed 's/-/\t/g'|sed 's/_/\t/g'|sed 's/:.//g'|awk '{if($0~/HiC/)print $1"_"$2"_"$3"\t"$4+10000"\t"$5-10000;else print $1"\t"$2+10000"\t"$3-10000}'|sort -k1,1 -k2,2n) -wa|sort|uniq| \
awk '{print $0,"Pekin-SX-Specific"}' >> Diff-Type-LTR/2022-Final-Three-Class

bedtools subtract -f 1 -a <(awk '{print "chr"$2"\t"$3"\t"$4"\t"$1"-"$5}' CAU-LTR-Sim-ID|sed 's/chrchr/chr/g'|sort -k1,1 -k2,2n) -b Diff-Type-LTR/2022-Final-Three-Class| \
awk '{print $0,"Other"}' > FFF
cat FFF >> Diff-Type-LTR/2022-Final-Three-Class && rm FFF
awk -F\= '{print $NF}' 2022-Final-Three-Class > 2022-Final-Three-Class2
readlink -f 2022-Final-Three-Class 

#classification of cLTR from other genomes
bedtools subtract -f 1 -a <(awk '{print $2"\t"$3"\t"$4"\t"$1"-"$5}' SX-LTR-Sim-ID|sort -k1,1 -k2,2n) -b <(cat *genome2/SX*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort|uniq|sed 's/\_\-[0-9]*/\t0/g'|sed 's/-/\t/g'|sed 's/_/\t/g'|sed 's/:.//g'|sort -k1,1 -k2,2n|sed 's/\t\t/\t/g')| \
awk '{print $0,"SX-Specific"}' > Diff-Type-LTR/2022-Final-Other-Class 
bedtools subtract -f 1 -a <(awk '{print $2"\t"$3"\t"$4"\t"$1"-"$5}' Ma-LTR-Sim-ID|sort -k1,1 -k2,2n) -b <(cat *genome2/Ma*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort|uniq|sed 's/\_\-[0-9]*/\t0/g'|sed 's/-/\t/g'|sed 's/_/\t/g'|sed 's/:.//g'|sort -k1,1 -k2,2n|sed 's/\t\t/\t/g')| \
awk '{print $0,"Ma-Specific"}' >> Diff-Type-LTR/2022-Final-Other-Class
bedtools subtract -f 1 -a <(awk '{print "chr"$2"\t"$3"\t"$4"\t"$1"-"$5}' CAU-LTR-Sim-ID|sed 's/chrchr/chr/g'|sort -k1,1 -k2,2n) -b <(cat *genome2/CAU*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort|uniq|sed 's/-/\t/g'|sed 's/_/\t/g'|sed 's/:.//g'|awk '{if($0~/HiC/)print $1"_"$2"_"$3"\t"$4"\t"$5;else print $0}'|sort -k1,1 -k2,2n)| \
awk '{print $0,"Pekin-Specific"}' >> Diff-Type-LTR/2022-Final-Other-Class 
bedtools intersect -f 1 -a <(awk '{print $2"\t"$3"\t"$4"\t"$1"-"$5}' SX-LTR-Sim-ID|sort -k1,1 -k2,2n) -b <(cat *genome2/SX*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort|uniq|sed 's/\_\-[0-9]*/\t0/g'|sed 's/-/\t/g'|sed 's/_/\t/g'|sed 's/:.//g'|sort -k1,1 -k2,2n|sed 's/\t\t/\t/g') -wa| \
awk '{print $0,"SX-nonSpecific"}'|sort|uniq >> Diff-Type-LTR/2022-Final-Other-Class 
bedtools intersect -f 1 -a <(awk '{print $2"\t"$3"\t"$4"\t"$1"-"$5}' Ma-LTR-Sim-ID|sort -k1,1 -k2,2n) -b <(cat *genome2/Ma*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort|uniq|sed 's/\_\-[0-9]*/\t0/g'|sed 's/-/\t/g'|sed 's/_/\t/g'|sed 's/:.//g'|sort -k1,1 -k2,2n|sed 's/\t\t/\t/g') -wa| \
awk '{print $0,"Ma-nonSpecific"}'|sort|uniq >> Diff-Type-LTR/2022-Final-Other-Class
bedtools intersect -f 1 -a <(awk '{print "chr"$2"\t"$3"\t"$4"\t"$1"-"$5}' CAU-LTR-Sim-ID|sed 's/chrchr/chr/g'|sort -k1,1 -k2,2n) -b <(cat *genome2/CAU*2|awk '{if($(NF-4)>11000 && $5>95)print $1}'|sort|uniq|sed 's/-/\t/g'|sed 's/_/\t/g'|sed 's/:.//g'|awk '{if($0~/HiC/)print $1"_"$2"_"$3"\t"$4"\t"$5;else print $0}'|sort -k1,1 -k2,2n) -wa| \
awk '{print $0,"Pekin-nonSpecific"}'|sort|uniq >> Diff-Type-LTR/2022-Final-Other-Class 
awk -F\= '{print $NF}' 2022-Final-Other-Class > 2022-Final-Other-Class2



