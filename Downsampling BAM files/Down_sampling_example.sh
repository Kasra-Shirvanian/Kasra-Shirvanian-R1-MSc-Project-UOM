#!/bin/bash --login
#$ -cwd 
#$ -pe  smp.pe 4
#$ -o  ./eo_samtool_SP135194D_ds20x
#$ -e  ./errorfile_samtool_SP135194D_ds20x
#$ -N  samtool_SP135194D_ds20x

module load apps/bioinf
module load apps/gcc/samtools/1.13

prefix1="sample_code"
prefix2="final_coverage"
prefix3="initial_coverage"
bam_stored="/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D/"

input_bam="${bam_stored}${prefix1}_T_${prefix3}

tumour="${prefix1}_T_${prefix2}"
normal="${prefix1}_N_30"


#input_bam="${bam_stored}${normal}.bam"

downsampled_bam="${bam_stored}${tumour}.bam"
# example flags and random seed for downsampling from 30 times to 20 
samtools  view --subsample 0.666 --subsample-seed 225  -b $input_bam > $downsampled_bam

wait
# calculating the coverage
echo "average depth not including zero coverage" >  ${tumour}_coverageds.txt

samtools depth  $downsampled_bam | awk '{sum+=$3} END {print sum/NR}' >> ${tumour}_coverageds.txt

echo "average depth  including zero coverage" >>  ${tumour}_coverageds.txt

samtools  depth -a  $downsampled_bam | awk '{sum+=$3} END {print sum/NR}' >> ${tumour}_coverageds.txt

wait
#createing the index file
samtools index $downsampled_bam
