#!/bin/bash

##################################
### ENV                          #
#conda activate base             #
### Software                     #
#`vg == v1.43.0` binary version  #
#`bcftools` == 1.15.1            #
##################################


###### Construct graph index ###### 
### VG autoindex
#vg autoindex \
#-t 30 \
#--workflow giraffe \
#--prefix A.6.pan_duck_graph_chr_9_giraffe \
#--gfa ../report_minigraph/A.6.pan_duck_graph_chr_9.gfa \
#-R XG

###### Reads mapping & Calling SVs ###### 
for i in BJ SX Ma SB ZW LC LW PZ HL
do
    ### VG giraffe map
    # Pair-end
    vg giraffe \
    -t 30 \
    -Z A.6.pan_duck_graph_chr_9_giraffe.giraffe.gbz \
    -m A.6.pan_duck_graph_chr_9_giraffe.min \
    -d A.6.pan_duck_graph_chr_9_giraffe.dist \
    -x A.6.pan_duck_graph_chr_9_giraffe.xg \
    -f ../data/ngs_${i}_1.fq.gz -f ../data/ngs_${i}_2.fq.gz > A.6.graph_${i}.gam

    ### Using XG format graph
    vg pack -t 26 \
    -x A.6.pan_duck_graph_chr_9_giraffe.xg \
    -g A.6.graph_${i}.gam \
    -Q 5 -s 5 \
    -o A.6.graph_${i}.pack

    ### VG call
    # generate vcf file from graph genome
    # min-length INT，only call variants > INT
    vg call \
    A.6.pan_duck_graph_chr_9_giraffe.xg \
    -k A.6.graph_${i}.pack \
    -t 30 -A -a > A.6.graph_${i}_A_a.vcf
done


