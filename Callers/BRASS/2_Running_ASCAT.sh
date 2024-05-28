#!/bin/bash --login
#$ -cwd
#$ -l mem256
#$ -pe smp.pe 4
#$ -t 4-4
#$ -N ascat_SP135194

prefix="Sample_Code"

# Calculate the index for the current task (SGE_TASK_ID starts from 1)
INDEX=$((SGE_TASK_ID-1))
# Read the contents of valuelist.txt into an array named valuelist
mapfile -t valuelist < "valuelist.txt"

name_bit=$(echo ${valuelist[$INDEX]})
name="${prefix}_T_${name_bit}.bam"
normalbam="${prefix}_N_30.bam"

new_dir="/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D/BRASS_M/SP135194_T_${name_bit}/"

mkdir -p ${new_dir}
# Run the ascat.pl script using Singularity
singularity exec \
--cleanenv \
-i --bind /mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D:/mnt/in,${new_dir}:/mnt/out,/mnt/bmh01-rds/UoOxford_David_W/s99384ml/ML_software/cgpwgs_GRCh37d5reference:/mnt/reference,/mnt/bmh01-rds/UoOxford_David_W/s99384ml/ML_software/to_run_ascat_brass:/mnt/reference2 \
--home ${PWD}:/home  \
--workdir ${PWD} \
"/mnt/bmh01-rds/UoOxford_David_W/shared/pipelines/Sanger_hg37/dockstore-cgpwgs_2.1.0.sif" \
 ascat.pl \
  -outdir /mnt/out  \
  -tumour "/mnt/in/${name}" \
  -normal "/mnt/in/${normalbam}" \
  -reference "/mnt/reference/core_ref_GRCh37d5.tar.gz" \
  -snp_gc /mnt/reference2/SnpGcCorrections.tsv \
  -gender XX \
  -genderChr Y \
  -protocol WGS \
  -platform ILLUMINA \
  -species Human \
  -assembly GRCh37d5 \
  -cpus 4
# Execute the ascat.pl script within the Singularity container with the specified parameters