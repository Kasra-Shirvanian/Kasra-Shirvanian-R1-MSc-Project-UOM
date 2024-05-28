#!/bin/bash --login
#$ -cwd
#$ -l mem256
#$ -pe smp.pe 4
#$ -t 4-4
#$ -N BRASS_SP135194

prefix="Sample_Code"
INDEX=$((SGE_TASK_ID-1))
mapfile -t valuelist < "valuelist.txt"

name_bit=$(echo ${valuelist[$INDEX]})
name="${prefix}_T_${name_bit}.bam"
normalbam="${prefix}_N_30.bam"

new_dir="/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D/BRASS_M/SP135194_T_${name_bit}/"

cd ${new_dir}
# Run the brass.pl script using Singularity
singularity exec \
--cleanenv \
-i --bind /mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D:/mnt/in,${new_dir}:/mnt/out,/mnt/bmh01-rds/UoOxford_David_W/s99384ml/ML_software/cgpwgs_GRCh37d5reference/core_ref_GRCh37d5:/mnt/reference,/mnt/bmh01-rds/UoOxford_David_W/s99384ml/ML_software/to_run_ascat_brass:/mnt/reference2,/mnt/bmh01-rds/UoOxford_David_W/s99384ml/ML_software/to_run_ascat_brass/brass:/mnt/reference3,/mnt/bmh01-rds/UoOxford_David_W/s99384ml/ML_software/cgpwgs_GRCh37d5reference/VAGrENT_ref_GRCh37d5_ensembl_75/vagrent:/mnt/reference4  \
"/mnt/bmh01-rds/UoOxford_David_W/shared/pipelines/Sanger_hg37/dockstore-cgpwgs_2.1.0.sif" brass.pl \
-outdir /mnt/out  \
-tumour "/mnt/in/${name}" \
-normal "/mnt/in/${normalbam}" \
-genome /mnt/reference/genome.fa -protocol WGS -g_cache /mnt/reference4/vagrent.cache.gz -viral /mnt/reference3/viral.genomic.fa.2bit -microbe /mnt/reference3/all_ncbi_bacteria -gcbins /mnt/reference3/500bp_windows.gc.bed.gz -cytoband /mnt/reference3/cytoband.txt -centtel /mnt/reference3/CentTelo.tsv -species Human -assembly GRCh37d5 -cpus 4 -sampstat /mnt/out/7b482f46-adcc-48ae-a4d6-17f4edd998dc.samplestatistics.txt -depth /mnt/reference3/HiDepth.bed.gz
# Execute the brass.pl script within the Singularity container with the specified parameters