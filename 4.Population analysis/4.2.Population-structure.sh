#Z
#change chrZ to chrX
for i in `seq 2 7`
do
~/admixture_linux-1.3.0/admixture --cv 200-maf0.01-geno0.1-LD-50-10-0.5.bed $i --haploid="male:X"|tee log${i}.out &
done
#Autosome
for i in `seq 2 7`
do 
~/admixture_linux-1.3.0/admixture --cv 200-Autosome-maf0.01-geno0.1-LD-50-10-0.5.bed $i|tee log${i}.out &
done
