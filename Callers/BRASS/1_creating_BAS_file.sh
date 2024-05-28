#!/bin/bash --login
#$ -cwd
#$ -pe smp.pe 6
#$ -t 5-5
#$ -o ./eofiles_nBase
#$ -e ./eofiles_nBaso
#$ -N normalBas_create
INDEX=$((SGE_TASK_ID-1))
mapfile -t bamlist < "bamlist.txt"

# Get the BAM file name corresponding to the current index
name=$(echo ${bamlist[$INDEX]})


# Run a Singularity container with the specified bindings and command
singularity run \
--cleanenv \
--bind \
/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D:/mnt/in \
--bind /mnt/bmh01-rds/UoOxford_David_W/s99384ml/ML_software/cgpwgs_GRCh37d5reference/core_ref_GRCh37d5:/mnt/in2 \
--home ${PWD}:/home  \
--workdir ${PWD} \
"/mnt/bmh01-rds/UoOxford_David_W/s99384ml/ML_software/cgpmap_3.3.0/dockstore-cgpmap_3.3.0.sif" \
bam_stats -i /mnt/in/${name}  -o ${name}.bas -r /mnt/in2/genome.fa.fai
# Run the bam_stats command within the Singularity container, processing the specified BAM file and outputting the results