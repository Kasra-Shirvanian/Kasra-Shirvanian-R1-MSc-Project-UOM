#!/bin/bash --login
#$ -cwd
#$ -l mem256
#$ -o ./eo_SP135194B
#$ -e ./er_SP135194B
#$ -N SP135194B_survivor
 
module load apps/binapps/anaconda3/2021.11
source activate survivor1.0.7KS
 
# inputs
prefix1="Sample_Code"
prefix2="Coverage"

pcawg_vcf=/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/ground_truths/${prefix1}bedpe.vcf

brass_vcf=/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D/BRASS_M/${prefix1}_T_${prefix2}/Brass_${prefix1}_${prefix2}.vcf
#generate a pseudo random number so we can delete only these files
random_number=123
 
new_dir=/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D/Survivor/Brass_${prefix1}_${prefix2}_NS/
 
# we make this folder and all folders to it if they don't already exist
mkdir -p "${new_dir}"
 
file_name=${new_dir}"temp_${random_number}.txt"
 
#we need to echo our vcfs into our file
echo -e ${brass_vcf}>${file_name}
echo -e ${pcawg_vcf}>>${file_name}
 
#our prefix - change to your naming scheme
prefix=SP135436
 
# will hold the actual results from survivor intersection
orimerged_vcf=${new_dir}survivor_merged_${prefix1}_B.vcf
 
# will hold the results of the SVs In the original PCAWG vcf that were identified
merged_vcf=${new_dir}survivor_merged_${prefix1}_B_tidy.vcf
 
# a table for sensitivity anfd precision and F1 score
output_table=${new_dir}survivor_merged_${prefix1}_B_results.txt
 
#get ones in common - using suggested default filters - runs the survivor command.
#parameters can be adjusted
SURVIVOR merge ${file_name}  1000 2 1 0 0 30 ${orimerged_vcf}
 
wait
# clean up
#rm ${file_name}
# change this to fish as we did for MANTA or we won't get them all
# As want to have everything the same so can compare Id's we will grab the ones that are as in the original file
cat ${orimerged_vcf}|grep -v "#"|awk -F ":" '{print $(NF-3)}'  >merged_${random_number}_ids.txt
cat ${pcawg_vcf}|grep "#"> ${merged_vcf}
cat ${pcawg_vcf}| grep -wFf merged_${random_number}_ids.txt|grep -v "#">>${merged_vcf}
wait
 
#tidy
#rm ${file_name} merged_${random_number}_ids.txt
 
###################################################################
# Now need to extract SV types and compare for sensitivity precision and F1
###################################################################
 
#from original filtered brass file  - did all weirdly in case insertions
#from original filtered brass file  - did all weirdly in case insertions
brassdel=$(cat ${brass_vcf}| awk -F '\t' '{print $5}'|grep DEL| wc -l)
brassdup=$(cat ${brass_vcf}| awk -F '\t' '{print $5}'|grep DUP| wc -l)
brassinv=$(cat ${brass_vcf}| awk -F '\t' '{print $5}'|grep INV| wc -l)
brasstra=$(cat ${brass_vcf}| awk -F '\t' '{print $5}'|grep TRA| wc -l)
brassall=$((brassdel + brassdup + brassinv + brasstra))
 
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
allPr=$(python -c "mergedall=$mergedall; brassall=$brassall; print('{:.3f}'.format(mergedall / brassall))")
delPr=$(python -c "mergeddel=$mergeddel; brassdel=$brassdel; print('{:.3f}'.format(mergeddel / brassdel))")
dupPr=$(python -c "mergeddup=$mergeddup; brassdup=$brassdup; print('{:.3f}'.format(mergeddup / brassdup))")
invPr=$(python -c "mergedinv=$mergedinv; brassinv=$brassinv; print('{:.3f}'.format(mergedinv / brassinv))")
traPr=$(python -c "mergedtra=$mergedtra; brasstra=$brasstra; print('{:.3f}'.format(mergedtra / brasstra))")
 
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
echo -e "brasscount\t${brassall}\t${brassdel}\t${brassdup}\t${brassinv}\t${brasstra}" >> "${output_table}"
echo -e "intersect\t${mergedall}\t${mergeddel}\t${mergeddup}\t${mergedinv}\t${mergedtra}" >> "${output_table}"
echo -e "precision\t${allPr}\t${delPr}\t${dupPr}\t${invPr}\t${traPr}" >> "${output_table}"
echo -e "sensitivity\t${allS}\t${delS}\t${dupS}\t${invS}\t${traS}" >> "${output_table}"
echo -e "F1_score\t${allF1}\t${delF1}\t${dupF1}\t${invF1}\t${traF1}" >> "${output_table}"
