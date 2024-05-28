#!/bin/bash --login
#$ -cwd
#$ -l mem256
#$ -o ./eo_SP135194L
#$ -e ./er_SP135194L
#$ -N survivor_SP135194L
 
module load apps/binapps/anaconda3/2021.11
source activate survivor1.0.7KS
 
# inputs
prefix1="Sample_Code"
prefix2="Coverage"
Lumpy_stored="/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D/Callers/Lumpy/"
pcawg_vcf=/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/ground_truths/{prefix1}bedpe.vcf
lumpy_vcf="/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D/Callers/Lumpy/Lumpy_${prefix1}_${prefix2}/lumpy_playlumpynew_${prefix1}_union.vcf"
 
 
#generate a pseudo random number so we can delete only these files
random_number=123
 
#new_dir=MAKE a new directory for your results
new_dir=/mnt/bmh01-rds/UoOxford_David_W/s52440ks/Samtool/SP135194D/Survivor/Lumpy_${prefix1}_${prefix2}/
 
# we make this folder and all folders to it if they don't already exist
mkdir -p "${new_dir}"
 
file_name=${new_dir}"temp_${random_number}.txt"
 
#we need to echo our vcfs into our file
echo -e ${lumpy_vcf}>${file_name}
echo -e ${pcawg_vcf}>>${file_name}
 
#our prefix - change to your naming scheme
prefix=SP135436
 
# will hold the actual results from survivor intersection
orimerged_vcf=${new_dir}survivor_merged_${prefix1}_L.vcf
 
# will hold the results of the SVs In the original PCAWG vcf that were identified
merged_vcf=${new_dir}survivor_merged_${prefix1}_L_tidy.vcf
 
# a table for sensitivity anfd precision and F1 score
output_table=${new_dir}survivor_merged_${prefix1}_L_results.txt
 
#get ones in common - using suggested default filters - runs the survivor command.
#parameters can be adjusted
SURVIVOR merge ${file_name}  1000 2 1 1 0 30 ${orimerged_vcf}
 
wait
# clean up
#rm ${file_name}
 
# As want to have everything the same so can compare Id's we will grab the ones that are as in the original file
# change this to fish as we did for MANTA or we won't get them all

cat ${orimerged_vcf}|grep -v "#"|awk -F ":" '{print $(NF-3)}'  >merged_${random_number}_ids.txt
cat ${pcawg_vcf}|grep "#"> ${merged_vcf}
cat ${pcawg_vcf}| grep -wFf merged_${random_number}_ids.txt|grep -v "#">>${merged_vcf}
wait
 
#tidy
#rm ${file_name} merged_${random_number}_ids.txt
 
###################################################################
# Now need to extract SV types and compare for sensitivity precision and F1
###################################################################
 
 
#from original filtered lumpy file  - did all weirdly in case insertions
lumpydel=$(cat ${lumpy_vcf}| grep -v "#"| grep "SVTYPE=DEL"| wc -l)
lumpydup=$(cat ${lumpy_vcf}| grep -v "#"|grep "SVTYPE=DUP"| wc -l)
lumpyinv=$(cat ${lumpy_vcf}| grep -v "#"|grep "SVTYPE=INV"| wc -l)
lumpytra=$(cat ${lumpy_vcf}| grep -v "#"|grep "SVTYPE=TRA"| wc -l)
lumpyall=$((lumpydel + lumpydup + lumpyinv + lumpytra))
 
 
# from the survivor tidied output
mergeddel=$(cat ${merged_vcf}| grep -v "#"|grep "SVTYPE=DEL"| wc -l)
mergeddup=$(cat ${merged_vcf}| grep -v "#"|grep "SVTYPE=DUP"| wc -l)
mergedinv=$(cat ${merged_vcf}| grep -v "#"|grep "SVTYPE=INV"| wc -l)
mergedtra=$(cat ${merged_vcf}| grep -v "#"|grep "SVTYPE=TRA"| wc -l)
mergedall=$((mergeddel + mergeddup + mergedinv + mergedtra))
 
 
# numbers from the pcawg_vcf
pcawgdel=$(cat ${pcawg_vcf}| grep -v "#"|grep "SVTYPE=DEL"| wc -l)
pcawgdup=$(cat ${pcawg_vcf}| grep -v "#"|grep "SVTYPE=DUP"| wc -l)
pcawginv=$(cat ${pcawg_vcf}| grep -v "#"|grep "SVTYPE=INV"| wc -l)
pcawgtra=$(cat ${pcawg_vcf}| grep -v "#"|grep "SVTYPE=TRA"| wc -l)
pcawgall=$((pcawgdel + pcawgdup + pcawginv + pcawgtra))
 
 
# Precision tp/(tp+fp) - math was going weird so using Python
allPr=$(python -c "mergedall=$mergedall; lumpyall=$lumpyall; print('{:.3f}'.format(mergedall / lumpyall))")
delPr=$(python -c "mergeddel=$mergeddel; lumpydel=$lumpydel; print('{:.3f}'.format(mergeddel / lumpydel))")
dupPr=$(python -c "mergeddup=$mergeddup; lumpydup=$lumpydup; print('{:.3f}'.format(mergeddup / lumpydup))")
invPr=$(python -c "mergedinv=$mergedinv; lumpyinv=$lumpyinv; print('{:.3f}'.format(mergedinv / lumpyinv))")
traPr=$(python -c "mergedtra=$mergedtra; lumpytra=$lumpytra; print('{:.3f}'.format(mergedtra / lumpytra))")
 
 
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
echo -e "lumpycount\t${lumpyall}\t${lumpydel}\t${lumpydup}\t${lumpyinv}\t${lumpytra}" >> "${output_table}"
echo -e "intersect\t${mergedall}\t${mergeddel}\t${mergeddup}\t${mergedinv}\t${mergedtra}" >> "${output_table}"
echo -e "precision\t${allPr}\t${delPr}\t${dupPr}\t${invPr}\t${traPr}" >> "${output_table}"
echo -e "sensitivity\t${allS}\t${delS}\t${dupS}\t${invS}\t${traS}" >> "${output_table}"
echo -e "F1_score\t${allF1}\t${delF1}\t${dupF1}\t${invF1}\t${traF1}" >> "${output_table}"
