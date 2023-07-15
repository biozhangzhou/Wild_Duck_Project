nohup cactus /tmp/ evolver25Anas.txt 25-Bird.hal &


for i in `seq 1 31` Z W `seq 34 41`
do 
  ~/1.bin/cactus-bin-v1.2.0/bin/hal2maf ../../../25-Bird.hal Chr${i}\-25-Bird.maf \
  --refGenome CAU \
  --refSequence Chr${i} \
  --onlyOrthologs \
  --noDupe & 
done


for i in `ls Chr*maf`
do 
  {
  mafFilter -minCol=100 $i > Filter-${i}
  perl ./01.convertMaf2List.pl Filter-${i} CAU
  }&
done 
