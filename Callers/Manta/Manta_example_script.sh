#!/bin/bash --login
#$ -cwd
#$ -pe smp.pe 4
#$ -o ./eofiles_mantaop
#$ -e  ./erfiles_mantaep
#$ -N MANTA_SP135194
 
module load apps/manta/1.4.0
module load apps/bioinf
module load  apps/python3/3.6.4/gcc-4.8.5
module load apps/gcc/bcftools/1.11
 
 
# make a directory to run it in for all samples no evidence
 
prefix1="Sample_Code"
prefix2="Bamfile_Coverage"
parent_dir="/mnt/bmh01-rds/UoOxford_David_W/s99384ml/GACA_CN_Projects/ISR684698/data/decrypted_files/"
basefolder="/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D/Callers/Manta/"
name_bit="MANTA_${prefix1}_${prefix2}"
new_dir=${basefolder}Manta_${prefix1}_${prefix2}
mkdir -p ${new_dir}
# Define the paths to the tumour and normal BAM files
tumour="${parent_dir}EGAZ00001224996/EGAZ00001224996_814a3122-bdc1-459c-ab8d-176d1fcc3694_cae75e168fdab75b63353fd276c62492.bam"
normal="${parent_dir}EGAZ00001225027/EGAZ00001225027_cda26432-0d12-4b44-82fd-e6a35a708b1c_f94b3e298940bd3a123b54d3276cecfd.bam"

# Define the output VCF file name

output="MANTAPASS_${prefix1}_${prefix2}.vcf"

# Configure Manta for the analysis 
configManta.py \
--normalBam ${normal} \
--tumorBam  ${tumour} \
--referenceFasta /mnt/bmh01-rds/UoOxford_David_W/s99384ml/ML_software/cgpwgs_GRCh37d5reference/core_ref_GRCh37d5/genome.fa \
--runDir ${new_dir}
 
wait
cd ${new_dir}
# Run the Manta workflow
./runWorkflow.py -m local -j 8


# Filter the resulting VCF file for high-quality structural variants
 
wait
# SV quality filtering
            bcftools filter \
                -O v  \
                -o  ${output} \
                -i "FILTER == 'PASS'" \
               "${new_dir}/results/variants/somaticSV.vcf.gz"

