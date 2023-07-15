### Create chromosome level graph
# chr1-9
for i in {1..9..1}
do
vg find -x A.6.pan_duck_graph_chr_9.xg -p chr0${i} -c 10000 > A.6.pan_duck_graph_chr_9_chr0${i}.xg
done
# chr10-39
for i in {10..39..1}
do
vg find -x A.6.pan_duck_graph_chr_9.xg -p chr${i} -c 10000 > A.6.pan_duck_graph_chr_9_chr${i}.xg
done
# chrZ
vg find -x A.6.pan_duck_graph_chr_9.xg -p chrZ -c 10000 > A.6.pan_duck_graph_chr_9_chrZ.xg

### Statistic
# rm stats
rm Stats3_A.6.graph.txt
# rm cache
rm `ls cache*`

for i in `ls A.6.pan_duck_graph_chr_9_chr*.xg`
do
echo "${i}" > cache1
vg stats -N $i > cache2 
vg stats -E $i > cache3
vg stats -l $i > cache4
paste `ls cache*` >> Stats3_A.6.graph.txt

# rm cache
rm `ls cache*`
done


