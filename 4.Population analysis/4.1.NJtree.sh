#200
#autosome
~/plink2/plink --bfile ../admixture-200/200-maf0.01-geno0.1-LD-50-10-0.5 --distance square 1-ibs flat-missing --out 200-maf0.01-geno0.1-LD-50-10-0.5 --chr-set 40
awk '{print $2"          "}' 200-Pop.txt| cut -c1-10 > IDforPhylip.txt
paste -d' ' IDforPhylip.txt 200-maf0.01-geno0.1-LD-50-10-0.5.mdist| sed 's/ //; 1i 200' > LD-phylip.mdist
#Z
~/plink2/plink --bfile ../200-Z-LD-50-10-0.5 --recode --out 200-Hapolotype-NonPAR-LD-50-10-0.5 --chr-set 40 &
sed -i 's/41/1/g' 200-Hapolotype-NonPAR-LD-50-10-0.5.map 
~/plink2/plink --file 200-Hapolotype-NonPAR-LD-50-10-0.5 --distance square 1-ibs flat-missing --out 200-Hapolotype-NonPAR-LD-50-10-0.5 --chr-set 40
awk '{print $1"          "}' 200-Pop.txt| cut -c1-10 > IDforPhylip.txt
paste -d' ' IDforPhylip.txt 200-Hapolotype-NonPAR-LD-50-10-0.5.mdist| sed 's/ //; 1i 200' > LD-phylip.mdist  


#180
#autosome
~/plink2/plink --bfile ../180-maf0.01-geno0.1-LD-50-10-0.5 --distance square 1-ibs flat-missing --out 180-maf0.01-geno0.1-LD-50-10-0.5 --chr-set 40
awk '{print $1"          "}' 180-Pop.txt| cut -c1-10 > IDforPhylip.txt
paste -d' ' IDforPhylip.txt 180-maf0.01-geno0.1-LD-50-10-0.5.mdist| sed 's/ //; 1i 180' > LD-phylip.mdist  

#Z
~/plink2/plink --bfile ../180-Hapolotype-NonPAR-LD-50-10-0.5 --recode --out 180-Hapolotype-NonPAR-LD-50-10-0.5 --chr-set 40 &
sed -i 's/41/1/g' 180-Hapolotype-NonPAR-LD-50-10-0.5.map 
~/plink2/plink --file 180-Hapolotype-NonPAR-LD-50-10-0.5 --distance square 1-ibs flat-missing --out 180-Hapolotype-NonPAR-LD-50-10-0.5 --chr-set 40
awk '{print $1"          "}' 180-Pop.txt| cut -c1-10 > IDforPhylip.txt
paste -d' ' IDforPhylip.txt 180-Hapolotype-NonPAR-LD-50-10-0.5.mdist| sed 's/ //; 1i 180' > LD-phylip.mdist  
