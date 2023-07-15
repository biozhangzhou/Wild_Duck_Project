#1 Reads prepare
combine all ONT pass reads;
change .fastq to .fasta
ls ONT_pass.fasta > input.fofn
#2 assembly
nextDenovo2.2 run.cfg
######about run.cfg ########
[General]
job_type = local
job_prefix = nextDenovo
task = all # 'all', 'correct', 'assemble'
rewrite = yes # yes/no
deltmp = yes
rerun = 3
parallel_jobs = 20 #could be imporve,be care about para "-t 6" * 20 =120
input_type = raw
input_fofn = ./input.fofn
workdir = ./01_rundir

[correct_option]
read_cutoff = 1k
seed_cutoff = 20000
blocksize = 1g
pa_correction = 6
seed_cutfiles = 2
sort_options = -m 1g -t 2 -k 50
minimap2_options_raw = -x ava-ont -t 6
correction_options = -p 15 #should be match with pa_correction

[assemble_option]
random_round = 20
minimap2_options_cns = -x ava-ont -t 8 -k17 -w17
nextgraph_options = -a 1
#############################
#3 nextpolish
#3.1 NGS reads QC
fastp -i file_R1 -o file_clean_R1 -I file_R2 -O file_clean_R2
ls NGS_clean_reads_R1.fq NGS_clean_reads_R2.fq > sgs.fofn
ls ONT_reads.fa > lgs.fofn 
#3.2 polish
nextpolish1.1.0 run.cfg
######about run.cfg ########
[General]
job_type = local
job_prefix = nextPolish
task = default
rewrite = yes
rerun = 3
parallel_jobs = 30
multithread_jobs = 3
genome = round05.fasta #Choose from the assembly results by N50 and number of contigs 
genome_size = auto
workdir = ./01_rundir
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = sgs.fofn
sgs_options = -max_depth 100

[lgs_option]
lgs_fofn = lgs.fofn
lgs_options = -min_read_len 10k -max_read_len 150k -max_depth 60
lgs_minimap2_options = -x map-ont
##############################