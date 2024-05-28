#!/bin/bash --login
#$ -cwd
#$ -l mem256
#$ -o ./eofiles_survivoroSP135194D
#$ -e  ./eofiles_survivoreSP135194D
#$ -N survivor_SP135194D
 
module load apps/binapps/anaconda3/2021.11
source activate survivor1.0.7KS
 
# inputs
prefix1="Sample_Code"
prefix2="Coverage"
delly_stored="/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D/Callers/Delly/"
pcawg_vcf="/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/ground_truths/{prefix1}bedpe.vcf"
#DELLY_vcf="${delly_stored}Delly_${prefix1}_${prefix2}/dellyPASS_${prefix1}_new.vcf"
DELLY_vcf="/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D/Callers/Delly/Delly_${prefix1}_${prefix2}/dellyPASS_${prefix1}_new.vcf"
#generate a pseudo random number so we can delete only these files
random_number="123"
 
new_dir=/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D/Survivor/test/Delly_${prefix1}_${prefix2}/
 
# we make this folder and all folders to it if they don't already exist
mkdir -p "${new_dir}"
 
file_name=${new_dir}"temp_${random_number}.txt"
 
#we need to echo our vcfs into our file. Always with the pcawg_vcf last
echo -e ${DELLY_vcf}>${file_name}
echo -e ${pcawg_vcf}>>${file_name}
 
#our prefix - change to your naming scheme
 
# will hold the actual results from survivor intersection
orimerged_vcf=${new_dir}survivor_merged_${prefix1}_D.vcf
 
# will hold the results of the SVs In the original PCAWG vcf that were identified
merged_vcf=${new_dir}survivor_merged_${prefix1}_D_tidy.vcf
 
# a table for sensitivity anfd precision and F1 score
output_table=${new_dir}survivor_merged_${prefix1}_D_results.txt
 
#get ones in common - using suggested default filters - runs the survivor command.
#parameters can be adjusted
SURVIVOR merge ${file_name}  1000 2 1 1 0 30 ${orimerged_vcf}
 
wait
# clean up
#rm ${file_name}
 
# As want to have everything the same so can compare Id's we will grab the ones that are as in the original file
# change this to fish as we did for MANTA or we won't get them all
cat ${orimerged_vcf}|grep -v "#"|awk -F ":" '{print $(NF-3)}' > merged_${random_number}_ids.txt
cat ${pcawg_vcf}|grep "#"> ${merged_vcf}
cat ${pcawg_vcf}| grep -wFf merged_${random_number}_ids.txt|grep -v "#">>${merged_vcf}
wait
 
#tidy
#rm ${file_name} merged_${random_number}_ids.txt
 
###################################################################
# Now need to extract SV types and compare for sensitivity precision and F1
###################################################################
 
#from original filtered DELLY file  - did all weirdly in case insertions
DELLYdel=$(cat ${DELLY_vcf}| grep -v "#"|awk -F 'SVTYPE=' '{print $2}'|awk -F ';' '{print $1}'|grep DEL| wc -l)
DELLYdup=$(cat ${DELLY_vcf}| grep -v "#"|awk -F 'SVTYPE=' '{print $2}'|awk -F ';' '{print $1}'|grep DUP| wc -l)
DELLYinv=$(cat ${DELLY_vcf}| grep -v "#"|awk -F 'SVTYPE=' '{print $2}'|awk -F ';' '{print $1}'|grep INV| wc -l)
DELLYtra=$(cat ${DELLY_vcf}| grep -v "#"|awk -F 'SVTYPE=' '{print $2}'|awk -F ';' '{print $1}'|grep BND| wc -l)
DELLYall=$((DELLYdel + DELLYdup + DELLYinv + DELLYtra))
 
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
allPr=$(python -c "mergedall=$mergedall; DELLYall=$DELLYall; print('{:.3f}'.format(mergedall / DELLYall))")
delPr=$(python -c "mergeddel=$mergeddel; DELLYdel=$DELLYdel; print('{:.3f}'.format(mergeddel / DELLYdel))")
dupPr=$(python -c "mergeddup=$mergeddup; DELLYdup=$DELLYdup; print('{:.3f}'.format(mergeddup / DELLYdup))")
invPr=$(python -c "mergedinv=$mergedinv; DELLYinv=$DELLYinv; print('{:.3f}'.format(mergedinv / DELLYinv))")
traPr=$(python -c "mergedtra=$mergedtra; DELLYtra=$DELLYtra; print('{:.3f}'.format(mergedtra / DELLYtra))")
 
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
echo -e "DELLYcount\t${DELLYall}\t${DELLYdel}\t${DELLYdup}\t${DELLYinv}\t${DELLYtra}" >> "${output_table}"
echo -e "intersect\t${mergedall}\t${mergeddel}\t${mergeddup}\t${mergedinv}\t${mergedtra}" >> "${output_table}"
echo -e "precision\t${allPr}\t${delPr}\t${dupPr}\t${invPr}\t${traPr}" >> "${output_table}"
echo -e "sensitivity\t${allS}\t${delS}\t${dupS}\t${invS}\t${traS}" >> "${output_table}"
echo -e "F1_score\t${allF1}\t${delF1}\t${dupF1}\t${invF1}\t${traF1}" >> "${output_table}"