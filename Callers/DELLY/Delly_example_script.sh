#!/bin/bash --login
#$ -cwd
#$ -pe smp.pe 4
#$ -o ./eofiles_mantao_og
#$ -e  ./erfiles_mantae_og
#$ -N Delly_SP135194_og

prefix1="Sample_Code"
prefix2="Coverage"
module load apps/gcc/bcftools/1.11
output_folder="/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D/Callers/Delly/Delly_${prefix1}_${prefix2}"

mkdir -p ${output_folder}
# Change directory to Delly software location
cd /mnt/bmh01-rds/UoOxford_David_W/s99384ml/ML_software/delly_v1.2.6/
# Run Delly using Singularity
singularity exec \
        --bind /mnt/bmh01-rds/UoOxford_David_W/s99384ml/ML_software/cgpwgs_GRCh37d5reference/core_ref_GRCh37d5:/ref1:ro  \
         --bind /mnt/bmh01-rds/UoOxford_David_W/shared/pipelines/EO_hg37/SV_pipeline/EO_SV_ref_hg37/reference-delly:/ref2:ro \
                 --bind  /mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D:/data1:ro  \
--bind  /mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D:/data2:ro  \
        --bind ${output_folder}:/output:rw  \
        --workdir ${PWD} \
 /mnt/bmh01-rds/UoOxford_David_W/s99384ml/ML_software/delly_v1.2.6/delly_v1.2.6.sif delly call -g /ref1/genome.fa -o /output/delly_${prefix1}.bcf -x /ref2/human.hg19.excl.tsv /data2/EGAZ00001224996_814a3122-bdc1-459c-ab8d-176d1fcc3694_cae75e168fdab75b63353fd276c62492.bam  /data1/EGAZ00001225027_cda26432-0d12-4b44-82fd-e6a35a708b1c_f94b3e298940bd3a123b54d3276cecfd.bam

wait


# file has both germline and somatic and must be filtered

# first need to make a table with the headers from the vcf file for the tumor and control
bcftools view ${output_folder}/delly_${prefix1}.bcf|grep '^#CHROM' | awk -F '\t' '{print $10 "\ttumor\n"$11"\tcontrol"}'>${output_folder}/samples.txt 
wait
# now can run Delly's optimised filtering that gets rid of any with possible evidence in controls
singularity exec \
        --bind ${output_folder}:/output:rw  \
        --workdir ${PWD} \
 /mnt/bmh01-rds/UoOxford_David_W/s99384ml/ML_software/delly_v1.2.6/delly_v1.2.6.sif delly filter -f somatic -s /output/samples.txt  -o /output/bdelly_${prefix1}.bcf /output/delly_${prefix1}.bcf

wait

# Filter the resulting BCF file for high-quality structural variants
            bcftools filter \
                -O v  \
                -o   ${output_folder}/dellyPASS_${prefix1}_new.vcf  \
                -i "FILTER == 'PASS'" \
               ${output_folder}/bdelly_${prefix1}.bcf
