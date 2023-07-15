### Env
# conda activate R

### Software
# vg == v1.43.0

########## HiFi reads #####################

for i in Ma SB
do
### Alignment
vg giraffe \
-t 20 --align-from-chains \
-Z A.6.pan_duck_graph_chr_9_giraffe.giraffe.gbz \
-m A.6.pan_duck_graph_chr_9_giraffe.min \
-d A.6.pan_duck_graph_chr_9_giraffe.dist \
-x A.6.pan_duck_graph_chr_9_giraffe.xg \
-f ../data/1.A.3.pan_duck_3rd/hifi_${i}.fastq.gz \
-o gaf > A.6.graph_giraffe_${i}.gaf

### Split file
split -a 2 -l 100000 -d A.6.graph_giraffe_${i}.gaf cache/A.6.graph_giraffe_${i}_split_gaf

### Filter
	for split_num in `ls cache/A.6.graph_giraffe_${i}_split_gaf??`    
	do
	Rscript S1.filter.gaf.R \
	-i ${split_num} \
	-o ${split_num}.filter
	done

### Merge
cat cache/A.6.graph_giraffe_${i}_split_gaf*.filter > A.6.graph_giraffe_${i}_filter.gaf

### vg pack
#vg pack -t 20 \
#-x A.6.pan_duck_graph_chr_9_giraffe.xg \
#-a A.6.graph_giraffe_${i}_filter.gaf \
#-o A.6.graph_giraffe_${i}_filter.pack

### Circulate depth
#vg pack -D -t 20 \
#-i A.6.graph_giraffe_${i}_filter.pack \
#-x A.6.pan_duck_graph_chr_9_giraffe.xg > A.6.graph_giraffe_${i}_depth.txt
done



########## ONT reads #####################
for i in BJ SX PZ ZW LC LW HL
do
### Alignment
vg giraffe \
-t 20 --align-from-chains \
-Z A.6.pan_duck_graph_chr_9_giraffe.giraffe.gbz \
-m A.6.pan_duck_graph_chr_9_giraffe.min \
-d A.6.pan_duck_graph_chr_9_giraffe.dist \
-x A.6.pan_duck_graph_chr_9_giraffe.xg \
-f ../data/1.A.3.pan_duck_3rd/ont_${i}.fastq.gz \
-o gaf > A.6.graph_giraffe_${i}.gaf

### Split file
split -a 2 -l 100000 -d A.6.graph_giraffe_${i}.gaf cache/A.6.graph_giraffe_${i}_split_gaf

### Filter
        for split_num in `ls cache/A.6.graph_giraffe_${i}_split_gaf??`
        do
        Rscript S1.filter.gaf.R \
        -i ${split_num} \
        -o ${split_num}.filter \
	-q 50
        done

### Merge
cat cache/A.6.graph_giraffe_${i}_split_gaf*.filter > A.6.graph_giraffe_${i}_filter.gaf

### vg pack
#vg pack -t 20 \
#-x A.6.pan_duck_graph_chr_9_giraffe.xg \
#-a A.6.graph_giraffe_${i}_filter.gaf \
#-o A.6.graph_giraffe_${i}_filter.pack

### Circulate depth
#vg pack -D -t 20 \
#-i A.6.graph_giraffe_${i}_filter.pack \
#-x A.6.pan_duck_graph_chr_9_giraffe.xg > A.6.graph_giraffe_${i}_depth.txt

done









