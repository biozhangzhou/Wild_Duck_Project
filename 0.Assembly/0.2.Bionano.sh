#Bionano Denovo Assembly#
#Relate scripts#
fa2map=/opt/software/Solve3.3_10252018/HybridScaffold/10252018/scripts/fa2cmap_multi_color.pl
hybridScaffoldDir=/opt/software/Solve3.3_10252018/HybridScaffold/10252018
pipelineDir=/opt/software/Solve3.3_10252018/Pipeline/10252018
RefAlignerDir=/opt/software/Solve3.3_10252018/RefAligner/7915.7989rel
#2.1 BNX Filter
perl ${pipelineDir}/filter_SNR_dynamic.pl -i ${input_bnx} -o ${output} -P diag_hist.pdf 2>snr.log
#2.2 BNX filter minlength(kb)
${RefAlignerDir}/RefAligner -i ${output}.bnx -minlen 120 -merge -o ${output}.filter -bnx 2>run.log
#2.3 Merge (Optional)
${RefAlignerDir}/RefAligner -i Bionano_clean1.filter.bnx  -i Bionano_clean2.filter.bnx -merge -o Bionano_clean_Total.filter -bnx 
#2.4 Denovo assemble the BNX file without Reference (based on Solve)
/opt/software/anaconda2/bin/python ${pipelineDir}/pipelineCL.py \
-T 72 -j 72 -N 4 -f -R \
-i 5 \
-b ${WorkDir}/Input_clean_Total.filter.bnx \
-l ${WorkDir}/bnxAssembly \
-t /opt/software/Solve3.3_10252018/RefAligner/7915.7989rel/ \
-a /opt/software/Solve3.3_10252018/RefAligner/7915.7989rel/optArguments_nonhaplotype_noES_noCut_DLE1_saphyr.xml &>denovoBNX1.log
#2.5 Bionano Hybrid Assembly#
/usr/bin/perl $hybridScaffoldDir/hybridScaffold.pl -n ${input_fasta} \
-b ${output}.cmap \
-c $hybridScaffoldDir/hybridScaffold_DLE1_config.xml \
-r /opt/software/Solve3.3_10252018/RefAligner/7915.7989rel/RefAligner \
-o ${output} \
-B 2 -N 2 -f

combine SCAFFOLD and NOT SCAFFOLD .fasta

#2.5 Optional（Bionano Manual adjustment） 