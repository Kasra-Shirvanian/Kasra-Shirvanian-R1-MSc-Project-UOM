#!/bin/bash --login
#$ -cwd
#$ -pe smp.pe 4
#$ -o ./eofiles_mantao_X4
#$ -e  ./erfiles_mantae_X4
#$ -N LUMPY_SP135194_X4

module load apps/bioinf
module load apps/gcc/samtools/1.13
module load apps/binapps/anaconda3/2021.11
source activate lumpy_sv_v0.3.1KS

prefix1="Sample_code"
prefix2="Coverage"

bam_stored=/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D/
 
basefolder=/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D/Callers/Lumpy/

name_bit="Lumpy_${prefix1}_${prefix2}"
new_dir=${basefolder}${name_bit}
mkdir -p ${new_dir}

# Define the paths to the input BAM files
inputTumorBAM=${bam_stored}${prefix1}_T_${prefix2}.bam
inputControlBAM=${bam_stored}${prefix1}_N_30.bam
tumourid="${prefix1}_T_${prefix2}"
controlid="${prefix1}_N_30"

# Extract the discordant paired-end alignments.
samtools view -b -F 1294 ${inputTumorBAM} > ${tumourid}".discordants.unsorted.bam"
wait
# Extract the split-read alignments
samtools view -h  ${inputTumorBAM} \
    | ~/.conda/envs/lumpy_sv_v0.3.1KS/bin/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
   > ${tumourid}".splitters.unsorted.bam"

wait
# Sort both alignments

samtools sort -@ ${NSLOTS}  -o ${tumourid}".discordants.bam" -O bam ${tumourid}".discordants.unsorted.bam"
samtools sort -@ ${NSLOTS}  -o ${tumourid}".splitters.bam" -O bam ${tumourid}".splitters.unsorted.bam"

wait
rm ${tumourid}".discordants.unsorted.bam" ${tumourid}".splitters.unsorted.bam"

wait 
# Extract the discordant paired-end alignments.
samtools view -b -F 1294 ${inputControlBAM} > ${controlid}".discordants.unsorted.bam"
wait
# Extract the split-read alignments
samtools view -h  ${inputControlBAM} \
  |  ~/.conda/envs/lumpy_sv_v0.3.1KS/bin/extractSplitReads_BwaMem -i stdin \
  | samtools view -Sb - \
  > ${controlid}".splitters.unsorted.bam"

wait
# Sort both alignments

samtools sort -@ ${NSLOTS}  -o ${controlid}".discordants.bam" -O bam ${controlid}".discordants.unsorted.bam"
samtools sort -@ ${NSLOTS}  -o ${controlid}".splitters.bam" -O bam ${controlid}".splitters.unsorted.bam"


rm ${controlid}".discordants.unsorted.bam" ${controlid}".splitters.unsorted.bam"

wait 

# Run Lumpy for structural variant detection
lumpyexpress \
   -B ${inputTumorBAM},${inputControlBAM}  \
  -S  ${basefolder}${tumourid}".splitters.bam",${basefolder}${controlid}".splitters.bam" \
 -D ${basefolder}${tumourid}".discordants.bam",${basefolder}${controlid}".discordants.bam" \
  -o ${basefolder}${tumourid}"_lumpyunfiltered.vcf"

wait

# sorting the vcf so that any with evidence in the control are eliminated and there are at least 2 reads that validate the SV in the tumour
cat ${tumourid}"_lumpyunfiltered.vcf"| grep "#">${tumourid}"_lumpyfiltered.vcf"
awk -F '\t' '!/^#/ && split($11, arr, ":") && arr[2] == 0' ${tumourid}"_lumpyunfiltered.vcf"|awk -F '\t' '!/^#/ && split($10, arr, ":") && arr[2] > 2' >>${tumourid}"_lumpyfiltered.vcf"
