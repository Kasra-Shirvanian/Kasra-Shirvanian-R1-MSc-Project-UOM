#!/bin/bash --login
#$ -cwd
#$ -l mem256
#$ -o ./eofiles_survivoroSP135194M
#$ -e  ./eofiles_survivoreSP135194M
#$ -N survivor_SP135194M

module load apps/binapps/anaconda3/2021.11
source activate survivor1.0.7KS

# inputs
prefix1="Sample_Code"
prefix2="Coverage"
manta_stored="/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D/Callers/Manta/"
pcawg_vcf="/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/ground_truths/${prefix1}bedpe.vcf"
MANTA_vcf="${manta_stored}MANTA_${prefix1}_${prefix2}/MANTAPASS_${prefix1}_${prefix2}.vcf"
#generate a pseudo random number so we can delete only these files
random_number=$$

new_dir=/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D/Survivor/Manta_${prefix1}_${prefix2}/
# we make this folder and all folders to it if they don't already exist
mkdir -p "${new_dir}"

file_name=${new_dir}"temp_${random_number}.txt"

#we need to echo our vcfs into our file 
echo -e ${MANTA_vcf}>${file_name}
echo -e ${pcawg_vcf}>>${file_name}

#our prefix - change to your naming scheme
prefix=SP135436

# will hold the actual results from survivor intersection
orimerged_vcf=${new_dir}survivor_merged_${prefix1}_M.vcf

# will hold the results of the SVs In the original PCAWG vcf that were identified
merged_vcf=${new_dir}survivor_merged_${prefix1}_M_tidy.vcf

# a table for sensitivity and precision and F1 score
output_table=${new_dir}survivor_merged_${prefix1}_M_results.txt

#get ones in common - using suggested default filters - runs the survivor command.
#parameters can be adjusted
SURVIVOR merge ${file_name}  1000 2 1 1 0 30 ${orimerged_vcf}

wait
# clean up
#rm ${file_name}

# As want to have everything the same so can compare Id's we will grab the ones that are as in the original file
cat ${orimerged_vcf}|grep -v "#"|awk -F ":" '{print $(NF-3)}'  >  merged_${random_number}_ids.txt
cat ${pcawg_vcf}|grep "#"> ${merged_vcf}
cat ${pcawg_vcf}| grep -wFf merged_${random_number}_ids.txt|grep -v "#">>${merged_vcf}
wait

#tidy
#rm ${file_name} merged_${random_number}_ids.txt

###################################################################
# Now need to extract SV types and compare for sensitivity precision and F1
###################################################################

#from original filtered MANTA file  - did all weirdly in case insertions 
MANTAdel=$(cat ${MANTA_vcf}| grep -v "#"|awk -F 'SVTYPE=' '{print $2}'|awk -F ';' '{print $1}'|grep DEL| wc -l)
MANTAdup=$(cat ${MANTA_vcf}| grep -v "#"|awk -F 'SVTYPE=' '{print $2}'|awk -F ';' '{print $1}'|grep DUP| wc -l)
MANTAinv=$(cat ${MANTA_vcf}| grep -v "#"|awk -F 'SVTYPE=' '{print $2}'|awk -F ';' '{print $1}'|grep INV| wc -l)
MANTAtra=$(cat ${MANTA_vcf}| grep -v "#"|awk -F 'SVTYPE=' '{print $2}'|awk -F ';' '{print $1}'|grep BND| wc -l)
MANTAall=$((MANTAdel + MANTAdup + MANTAinv + MANTAtra))

# from the survivor tidied output 
mergeddel=$(cat ${merged_vcf}| grep -v "#"|awk -F '\t' '{print $5}'|grep DEL| wc -l)
mergeddup=$(cat ${merged_vcf}| grep -v "#"|awk -F '\t' '{print $5}'|grep DUP| wc -l)
mergedinv=$(cat ${merged_vcf}| grep -v "#"|awk -F '\t' '{print $5}'|grep INV| wc -l)
mergedtra=$(cat ${merged_vcf}| grep -v "#"|awk -F '\t' '{print $5}'|grep TRA| wc -l)
mergedall=$((mergeddel + mergeddup + mergedinv + mergedtra))

# numbers from the pcawg_vcf
pcawgdel=$(cat ${pcawg_vcf}| grep -v "#"|awk -F '\t' '{print $5}'|grep DEL| wc -l)
pcawgdup=$(cat ${pcawg_vcf}| grep -v "#"|awk -F '\t' '{print $5}'|grep DUP| wc -l)
pcawginv=$(cat ${pcawg_vcf}| grep -v "#"|awk -F '\t' '{print $5}'|grep INV| wc -l)
pcawgtra=$(cat ${pcawg_vcf}| grep -v "#"|awk -F '\t' '{print $5}'|grep TRA| wc -l)
pcawgall=$((pcawgdel + pcawgdup + pcawginv + pcawgtra))

# Precision tp/(tp+fp) - math was going weird so using Python
allPr=$(python -c "mergedall=$mergedall; MANTAall=$MANTAall; print('{:.3f}'.format(mergedall / MANTAall))")
delPr=$(python -c "mergeddel=$mergeddel; MANTAdel=$MANTAdel; print('{:.3f}'.format(mergeddel / MANTAdel))")
dupPr=$(python -c "mergeddup=$mergeddup; MANTAdup=$MANTAdup; print('{:.3f}'.format(mergeddup / MANTAdup))")
invPr=$(python -c "mergedinv=$mergedinv; MANTAinv=$MANTAinv; print('{:.3f}'.format(mergedinv / MANTAinv))")
traPr=$(python -c "mergedtra=$mergedtra; MANTAtra=$MANTAtra; print('{:.3f}'.format(mergedtra / MANTAtra))")

# Sensitivity/recall tp/(tp+fn) - math was going weird so using Python
allS=$(python -c "mergedall=$mergedall; pcawgall=$pcawgall; print('{:.3f}'.format(mergedall / pcawgall))")
delS=$(python -c "mergeddel=$mergeddel; pcawgdel=$pcawgdel; print('{:.3f}'.format(mergeddel / pcawgdel))")
dupS=$(python -c "mergeddup=$mergeddup; pcawgdup=$pcawgdup; print('{:.3f}'.format(mergeddup / pcawgdup))")
invS=$(python -c "mergedinv=$mergedinv; pcawginv=$pcawginv; print('{:.3f}'.format(mergedinv / pcawginv))")
traS=$(python -c "mergedtra=$mergedtra; pcawgtra=$pcawgtra; print('{:.3f}'.format(mergedtra / pcawgtra))")

#F1 score = 2* precision* recall/(precision +recall)
# Perform the calculations with Python
allF1=$(python -c "print('{:.3f}'.format(2 * $allPr * $allS / ($allPr + $allS)))")
delF1=$(python -c "print('{:.3f}'.format(2 * $delPr * $delS / ($delPr + $delS)))")
dupF1=$(python -c "print('{:.3f}'.format(2 * $dupPr * $dupS / ($dupPr + $dupS)))")
invF1=$(python -c "print('{:.3f}'.format(2 * $invPr * $invS / ($invPr + $invS)))")
traF1=$(python -c "print('{:.3f}'.format(2 * $traPr * $traS / ($traPr + $traS)))")

# now make a tsv
echo -e "measure\tALL\tDEL\tDUP\tINV\tTRA" > "${output_table}"
echo -e "pcawgcount\t${pcawgall}\t${pcawgdel}\t${pcawgdup}\t${pcawginv}\t${pcawgtra}" >> "${output_table}"
echo -e "MANTAcount\t${MANTAall}\t${MANTAdel}\t${MANTAdup}\t${MANTAinv}\t${MANTAtra}" >> "${output_table}"
echo -e "intersect\t${mergedall}\t${mergeddel}\t${mergeddup}\t${mergedinv}\t${mergedtra}" >> "${output_table}"
echo -e "precision\t${allPr}\t${delPr}\t${dupPr}\t${invPr}\t${traPr}" >> "${output_table}"
echo -e "sensitivity\t${allS}\t${delS}\t${dupS}\t${invS}\t${traS}" >> "${output_table}"
echo -e "F1_score\t${allF1}\t${delF1}\t${dupF1}\t${invF1}\t${traF1}" >> "${output_table}"
