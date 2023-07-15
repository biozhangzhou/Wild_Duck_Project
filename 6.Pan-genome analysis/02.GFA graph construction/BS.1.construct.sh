### Running minigraph 
minigraph -cxggs -t 16 \
../genome/1.A.3.pan_duck_genome/BJ_chr.fa \
../genome/1.A.3.pan_duck_genome/SX_chr.fa \
../genome/1.A.3.pan_duck_genome/MA_chr.fa \
../genome/1.A.3.pan_duck_genome/BT_chr.fa \
../genome/1.A.3.pan_duck_genome/ZW_chr.fa \
../genome/1.A.3.pan_duck_genome/LC_chr.fa \
../genome/1.A.3.pan_duck_genome/LW_chr.fa \
../genome/1.A.3.pan_duck_genome/PZ_chr.fa \
../genome/1.A.3.pan_duck_genome/HL_chr.fa > A.6.pan_duck_graph_chr_9.gfa 

### Extract bubble from GFA graph 
gfatools bubble A.6.pan_duck_graph_chr_9.gfa > A.6.pan_duck_graph_chr_9.bed
