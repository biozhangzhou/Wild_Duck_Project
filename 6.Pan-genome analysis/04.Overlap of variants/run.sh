# Tuvari == 3.5.0
#-C, --chunksize: Max reference distance to compare calls (1000)
#-r, --refdist: Max reference location distance (500)
#-p, --pctsim: Min percent allele sequence similarity. Set to 0 to ignore.(0.7)
#-P, --pctsize: Min pct allele size similarity (minvarsize/maxvarsize)
#-O, --pctovl: Min pct reciprocal overlap (0.0)
#-t, --typeignore: Variant types don't need to match to compare (False)

base=$1
comp=$2
dir=$3

truvari bench \
-b $base \
-c $comp \
-f /home2/nizijia/projects/1.A.3.duck_graph_genome/genome/1.A.3.pan_duck_genome/BJ_chr.fa \
-o $dir \
-C 2000 -r 2000 -p 0.8 -P 0.7 -O 0.5 --sizemax 200000 -t --multimatch
