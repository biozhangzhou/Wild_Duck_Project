#treeshrink
nohup python ~/1.bin/TreeShrink/run_treeshrink.py -t 2-Chr1-30-1kbp.treeFile -O treeshrink &
nohup python ~/1.bin/TreeShrink/run_treeshrink.py -t 2-ChrZ-1kbp.treeFile -O treeshrink &

#Astral and discovista
awk '{if($0~/CAU/ && $0~/Ma/ && $0~/SB/ && $0~/SX/ && $0~/ZW/ && $0~/LC/ && $0~/LW/ && $0~/PZ/ && $0~/HL/)print $0}' ../../2-ChrZ-1K_treeshrink/treeshrink.treeFile > estimated_gene_trees.tree
java -jar ~/1.bin/Astral/astral.5.7.8.jar -i estimated_gene_trees.tree -o estimated_species_tree.tree 2>out.log 
python2 ~/1.bin/DiscoVista/src/utils/discoVista.py -a annotation-1.txt -m 5 -p Z/ -o Z/results/anno1 -g Base


for i in 1K 10K 100K
do
  cd Autosome
  ln -s ../${i}/2.Astral-DiscoVista/1-30/estimated_gene_trees.tree ./
  cd ../
  python2 ~/1.bin/DiscoVista/src/utils/discoVista.py -a ./annotation-1.txt -m 5 -p Autosome/ -o Autosome/result-${i} -g Base
  cd Z
  ln -s ../${i}/2.Astral-DiscoVista/Z/estimated_gene_trees.tree ./
  cd ../
  python2 ~/1.bin/DiscoVista/src/utils/discoVista.py -a ./annotation-1.txt -m 5 -p Z -o Z/result-${i} -g Base
  rm Autosome/estimated_gene_trees.tree Z/estimated_gene_trees.tree
done

#plot

for j in Autosome Z
do
for i in 1K 10K 100K
do
  awk -F\t '{print "'"$j"'\t'"$i"'\t"$0}' ${j}/result-${i}/freqQuadCorrected.csv >> ALL.freqQuadCorrected.csv
done
done

red='#d53e4f';orange='#1d91c0';blue='#41b6c4';colormap = c(red,orange,blue)
freq <- read.csv("ALL.freqQuadCorrected.csv",sep="\t",header=F)
freq$value = freq$V7/freq$V8
#freq$V10<-reorder(freq$V10,-freq$value)
freq$V2 <- factor(freq$V2,ordered = TRUE,level=c("100K","10K","1K"))
a<-length(levels(as.factor(freq$V9)))*3.7; b<-4; sizes <- c(a,b);
#ggplot(data=freq[freq$V9==1,])+aes(x=V10,y=value,fill=V11)+geom_bar(stat='identity',color=1,width=0.8,position='dodge')+theme_bw()+theme(axis.text.x=element_text(angle=90))+scale_fill_manual(values=colormap,name='Topology')+geom_hline(yintercept=1/3,size=0.4,linetype=2)+ylab('relative freq.')+facet_wrap(V1~V2,ncol=4,scales='free_x')+xlab('')
ggplot(freq)+aes(x=V11,y=value,fill=V4)+geom_bar(stat='identity',color=1,width=0.8,position='dodge')+theme_bw()+theme(axis.text.x=element_text(angle=90))+scale_fill_manual(values=colormap,name='Topology')+geom_hline(yintercept=1/3,size=0.4,linetype=2)+ylab('relative freq.')+facet_grid(V3~V1+V2,scales='free_x')+xlab('')
