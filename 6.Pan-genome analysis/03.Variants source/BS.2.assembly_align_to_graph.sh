##### Pipeline #####
### Mapping assembly to graph
for i in BJ SX MA BT ZW LC LW PZ HL
do
minigraph -xasm --cov -t 20 A.6.pan_duck_graph_chr_9.gfa ../genome/1.A.3.pan_duck_genome/${i}_chr.fa  > A.6.graph_${i}_chr.gfa
done

### Extract info from <graph_sample.gfa>
# Python script modified from cattle pangenome  
# https://github.com/AnimalGenomicsETH/bovine-graphs/blob/main/scripts/comb_coverage.py
python3 scripts/comb_coverage.py -g A.6.graph -a BJ_chr SX_chr MA_chr BT_chr ZW_chr LC_chr LW_chr PZ_chr HL_chr

