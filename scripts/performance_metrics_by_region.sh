#!/bin/bash
#SBATCH --job-name=performance_metrics_by_region # Job name
#SBATCH -o slurm.performance_metrics_by_region.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.performance_metrics_by_region.err                # STDERR (%j = JobId)
#SBATCH -n 1
#SBATCH -t 7-00:00:00
#SBATCH -q tempboost
#SBATCH --mail-user=amtarave@asu.edu # send-to address
#SBATCH --mail-type=ALL
#cd /scratch/amtarave/variant_calling_simulations_project/performance_metrics_by_region

# For this I want to subset FPs, TPs, FNs by:
# 1) PARs, nonPARs (minus XTR), XTR
# 2) Complex PARs, not complex PARs, complex nonPARs (minus XTR),
#    not complex nonPARs (minus XTR), XTR
# Do this for SCC and default alignment for males and females



################################################################################
#                                                                              #
# SCC aligned data #
#                                                                              #
################################################################################
#------------------------------------------------------------------------------#
# MALES - Chr 8 #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chr8_autos"
ploidy="diploid"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${pop}/${sex}/autos"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

cd ${inputdatapath}

for i in $samplelist
do
   for t in TP FP FN
   do
      # First, convert results files to bed format
      awk '{ print $1"\t"$2-1"\t"$2 }' ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.txt > ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.bed

      # Then subset on whether these calls are in complex, XTR or non complex regions
      # complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chr8.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chr8.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SimpleRepeats.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex.bed

      # Not complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chr8.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chr8.human.GRCh38.Repeats.RepeatMasker.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_not_complex.bed

   done

done

# Then make results table but counting the total number of lines in each file.
cd ${inputdatapath}/by_region

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt

for i in $samplelist
do
   # get counts of TP, FP, and FN for each category
   tpcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_complex.bed | wc -l`
   fpcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_complex.bed | wc -l`
   fncomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_complex.bed | wc -l`

   tpnotcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_not_complex.bed | wc -l`
   fpnotcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_not_complex.bed | wc -l`
   fnnotcomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_not_complex.bed | wc -l`

   # add to results file
   echo -e "${i}\t${tpcomplex}\t${fpcomplex}\t${fncomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
   echo -e "${i}\t${tpnotcomplex}\t${fpnotcomplex}\t${fnnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
done



#------------------------------------------------------------------------------#
# MALES - Chr X PARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chrX_PARs"
ploidy="diploid"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${pop}/${sex}/${chr}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

cd ${inputdatapath}

for i in $samplelist
do
   for t in TP FP FN
   do
      # First, convert results files to bed format
      awk '{ print $1"\t"$2-1"\t"$2 }' ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.txt > ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.bed

      # Then subset on whether these calls are in complex, XTR or non complex regions
      # complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex.bed

      # Not complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_not_complex.bed

   done

done

# Then make results table but counting the total number of lines in each file.
cd ${inputdatapath}/by_region

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt

for i in $samplelist
do
   # get counts of TP, FP, and FN for each category
   tpcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_complex.bed | wc -l`
   fpcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_complex.bed | wc -l`
   fncomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_complex.bed | wc -l`

   tpnotcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_not_complex.bed | wc -l`
   fpnotcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_not_complex.bed | wc -l`
   fnnotcomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_not_complex.bed | wc -l`

   # add to results file
   echo -e "${i}\t${tpcomplex}\t${fpcomplex}\t${fncomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
   echo -e "${i}\t${tpnotcomplex}\t${fpnotcomplex}\t${fnnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
done


#------------------------------------------------------------------------------#
# MALES - Chr X non-PARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chrX_nonPARs"
ploidy="haploid"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${pop}/${sex}/${chr}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

cd ${inputdatapath}

for i in $samplelist
do
   for t in TP FP FN
   do
      # First, convert results files to bed format
      awk '{ print $1"\t"$2-1"\t"$2 }' ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.txt > ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.bed

      # Then subset on whether these calls are in complex, XTR or non complex regions
      # complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX_amplicons_GRCh38.bed  -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed

      # make sure XTR is not in complex file
      bedtools intersect -a by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed -b ${regionsbed}/chrX_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex.bed

      # XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX_XTR_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_XTR.bed

      # Not complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX_amplicons_GRCh38.bed ${regionsbed}/chrX_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_not_complex.bed

      # non PARs minus XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed

      # amplicons only
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX_amplicons_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_Amplicons.bed

      # non PARs minus XTR and amplicons
      bedtools intersect -a ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed -b ${regionsbed}/chrX_amplicons_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTRminusAmplicons.bed

   done

done

# Then make results table but counting the total number of lines in each file.
cd ${inputdatapath}/by_region

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt

for i in $samplelist
do
   # get counts of TP, FP, and FN for each category
   tpcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_complex.bed | wc -l`
   fpcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_complex.bed | wc -l`
   fncomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_complex.bed | wc -l`

   tpXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_XTR.bed | wc -l`
   fpXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_XTR.bed | wc -l`
   fnXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_XTR.bed | wc -l`

   tpnotcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_not_complex.bed | wc -l`
   fpnotcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_not_complex.bed | wc -l`
   fnnotcomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_not_complex.bed | wc -l`

   tpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTR.bed | wc -l`
   fpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTR.bed | wc -l`
   fnnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTR.bed | wc -l`

   tpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fnnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`

   tpAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_Amplicons.bed | wc -l`
   fpAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_Amplicons.bed | wc -l`
   fnAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_Amplicons.bed | wc -l`

   # add to results file
   echo -e "${i}\t${tpcomplex}\t${fpcomplex}\t${fncomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
   echo -e "${i}\t${tpXTR}\t${fpXTR}\t${fnXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
   echo -e "${i}\t${tpnotcomplex}\t${fpnotcomplex}\t${fnnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
   echo -e "${i}\t${tpnonPARminusXTR}\t${fpnonPARminusXTR}\t${fnnonPARminusXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt

   echo -e "${i}\t${tpnonPARminusXTRminusAmplicons}\t${fpnonPARminusXTRminusAmplicons}\t${fnnonPARminusXTRminusAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt
   echo -e "${i}\t${tpAmplicons}\t${fpAmplicons}\t${fnAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt
done



#------------------------------------------------------------------------------#
# MALES - Chr X non-PARs #
# DIPLOID #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chrX_nonPARs"
ploidy="diploid"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${pop}/${sex}/${chr}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

cd ${inputdatapath}

for i in $samplelist
do
   for t in TP FP FN
   do
      # First, convert results files to bed format
      awk '{ print $1"\t"$2-1"\t"$2 }' ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.txt > ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.bed

      # Then subset on whether these calls are in complex, XTR or non complex regions
      # complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX_amplicons_GRCh38.bed  -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed

      # make sure XTR is not in complex file
      bedtools intersect -a by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed -b ${regionsbed}/chrX_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex.bed

      # XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX_XTR_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_XTR.bed

      # Not complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX_amplicons_GRCh38.bed ${regionsbed}/chrX_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_not_complex.bed

      # non PARs minus XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed

      # amplicons only
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX_amplicons_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_Amplicons.bed

      # non PARs minus XTR and amplicons
      bedtools intersect -a ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed -b ${regionsbed}/chrX_amplicons_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTRminusAmplicons.bed

   done

done

# Then make results table but counting the total number of lines in each file.
cd ${inputdatapath}/by_region

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt


for i in $samplelist
do
   # get counts of TP, FP, and FN for each category
   tpcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_complex.bed | wc -l`
   fpcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_complex.bed | wc -l`
   fncomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_complex.bed | wc -l`

   tpXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_XTR.bed | wc -l`
   fpXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_XTR.bed | wc -l`
   fnXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_XTR.bed | wc -l`

   tpnotcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_not_complex.bed | wc -l`
   fpnotcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_not_complex.bed | wc -l`
   fnnotcomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_not_complex.bed | wc -l`

   tpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTR.bed | wc -l`
   fpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTR.bed | wc -l`
   fnnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTR.bed | wc -l`

   tpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fnnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`

   tpAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_Amplicons.bed | wc -l`
   fpAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_Amplicons.bed | wc -l`
   fnAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_Amplicons.bed | wc -l`

   # add to results file
   echo -e "${i}\t${tpcomplex}\t${fpcomplex}\t${fncomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
   echo -e "${i}\t${tpXTR}\t${fpXTR}\t${fnXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
   echo -e "${i}\t${tpnotcomplex}\t${fpnotcomplex}\t${fnnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
   echo -e "${i}\t${tpnonPARminusXTR}\t${fpnonPARminusXTR}\t${fnnonPARminusXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt
   echo -e "${i}\t${tpnonPARminusXTRminusAmplicons}\t${fpnonPARminusXTRminusAmplicons}\t${fnnonPARminusXTRminusAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt
   echo -e "${i}\t${tpAmplicons}\t${fpAmplicons}\t${fnAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt
done



#------------------------------------------------------------------------------#
# MALES - Chr Y non-PARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chrY"
ploidy="haploid"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${pop}/${sex}/${chr}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"


cd ${inputdatapath}

for i in $samplelist
do
   for t in TP FP FN
   do
      # First, convert results files to bed format
      awk '{ print $1"\t"$2-1"\t"$2 }' ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.txt > ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.bed

      # Then subset on whether these calls are in complex, XTR or non complex regions
      # complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrY.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrY_amplicons_GRCh38.bed  -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed

      # make sure XTR is not in complex file
      bedtools intersect -a by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed -b ${regionsbed}/chrY_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex.bed

      # XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY_XTR_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_XTR.bed

      # Not complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrY.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrY_amplicons_GRCh38.bed ${regionsbed}/chrY_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_not_complex.bed

      # non PARs minus XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed

      # amplicons only
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY_amplicons_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_Amplicons.bed

      # non PARs minus XTR and amplicons
      bedtools intersect -a ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed -b ${regionsbed}/chrY_amplicons_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTRminusAmplicons.bed

   done

done

# Then make results table but counting the total number of lines in each file.
cd ${inputdatapath}/by_region

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt

for i in $samplelist
do
   # get counts of TP, FP, and FN for each category
   tpcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_complex.bed | wc -l`
   fpcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_complex.bed | wc -l`
   fncomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_complex.bed | wc -l`

   tpXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_XTR.bed | wc -l`
   fpXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_XTR.bed | wc -l`
   fnXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_XTR.bed | wc -l`

   tpnotcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_not_complex.bed | wc -l`
   fpnotcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_not_complex.bed | wc -l`
   fnnotcomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_not_complex.bed | wc -l`

   tpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTR.bed | wc -l`
   fpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTR.bed | wc -l`
   fnnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTR.bed | wc -l`

   tpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fnnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`

   tpAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_Amplicons.bed | wc -l`
   fpAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_Amplicons.bed | wc -l`
   fnAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_Amplicons.bed | wc -l`

   # add to results file
   echo -e "${i}\t${tpcomplex}\t${fpcomplex}\t${fncomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
   echo -e "${i}\t${tpXTR}\t${fpXTR}\t${fnXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
   echo -e "${i}\t${tpnotcomplex}\t${fpnotcomplex}\t${fnnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
   echo -e "${i}\t${tpnonPARminusXTR}\t${fpnonPARminusXTR}\t${fnnonPARminusXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt
   echo -e "${i}\t${tpnonPARminusXTRminusAmplicons}\t${fpnonPARminusXTRminusAmplicons}\t${fnnonPARminusXTRminusAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt
   echo -e "${i}\t${tpAmplicons}\t${fpAmplicons}\t${fnAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt

done


#------------------------------------------------------------------------------#
# MALES - Chr Y non-PARs #
# DIPLOID #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chrY"
ploidy="diploid"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${pop}/${sex}/${chr}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"


cd ${inputdatapath}

for i in $samplelist
do
   for t in TP FP FN
   do
      # First, convert results files to bed format
      awk '{ print $1"\t"$2-1"\t"$2 }' ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.txt > ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.bed

      # Then subset on whether these calls are in complex, XTR or non complex regions
      # complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrY.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrY_amplicons_GRCh38.bed  -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed

      # make sure XTR is not in complex file
      bedtools intersect -a by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed -b ${regionsbed}/chrY_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex.bed

      # XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY_XTR_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_XTR.bed

      # Not complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrY.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrY_amplicons_GRCh38.bed ${regionsbed}/chrY_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_not_complex.bed

      # non PARs minus XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed

      # amplicons only
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY_amplicons_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_Amplicons.bed

      # non PARs minus XTR and amplicons
      bedtools intersect -a ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed -b ${regionsbed}/chrY_amplicons_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTRminusAmplicons.bed

   done

done

# Then make results table but counting the total number of lines in each file.
cd ${inputdatapath}/by_region

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt

for i in $samplelist
do
   # get counts of TP, FP, and FN for each category
   tpcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_complex.bed | wc -l`
   fpcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_complex.bed | wc -l`
   fncomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_complex.bed | wc -l`

   tpXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_XTR.bed | wc -l`
   fpXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_XTR.bed | wc -l`
   fnXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_XTR.bed | wc -l`

   tpnotcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_not_complex.bed | wc -l`
   fpnotcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_not_complex.bed | wc -l`
   fnnotcomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_not_complex.bed | wc -l`

   tpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTR.bed | wc -l`
   fpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTR.bed | wc -l`
   fnnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTR.bed | wc -l`

   tpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fnnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`

   tpAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_Amplicons.bed | wc -l`
   fpAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_Amplicons.bed | wc -l`
   fnAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_Amplicons.bed | wc -l`

   # add to results file
   echo -e "${i}\t${tpcomplex}\t${fpcomplex}\t${fncomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
   echo -e "${i}\t${tpXTR}\t${fpXTR}\t${fnXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
   echo -e "${i}\t${tpnotcomplex}\t${fpnotcomplex}\t${fnnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
   echo -e "${i}\t${tpnonPARminusXTR}\t${fpnonPARminusXTR}\t${fnnonPARminusXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt
   echo -e "${i}\t${tpnonPARminusXTRminusAmplicons}\t${fpnonPARminusXTRminusAmplicons}\t${fnnonPARminusXTRminusAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt
   echo -e "${i}\t${tpAmplicons}\t${fpAmplicons}\t${fnAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt

done



#------------------------------------------------------------------------------#
# FEMALES - Chr 8 #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="females"
chr="chr8_autos"
ploidy="diploid"
samplelist="NA06985 NA06989 NA07000 NA07037 NA07056 NA10847 NA11830 NA11832 NA11840 NA11892"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${pop}/${sex}/autos"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

cd ${inputdatapath}

for i in $samplelist
do
   for t in TP FP FN
   do
      # First, convert results files to bed format
      awk '{ print $1"\t"$2-1"\t"$2 }' ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.txt > ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.bed

      # Then subset on whether these calls are in complex, XTR or non complex regions
      # complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chr8.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chr8.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SimpleRepeats.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex.bed

      # Not complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chr8.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chr8.human.GRCh38.Repeats.RepeatMasker.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_not_complex.bed

   done

done

# Then make results table but counting the total number of lines in each file.
cd ${inputdatapath}/by_region

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt

for i in $samplelist
do
   # get counts of TP, FP, and FN for each category
   tpcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_complex.bed | wc -l`
   fpcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_complex.bed | wc -l`
   fncomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_complex.bed | wc -l`

   tpnotcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_not_complex.bed | wc -l`
   fpnotcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_not_complex.bed | wc -l`
   fnnotcomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_not_complex.bed | wc -l`

   # add to results file
   echo -e "${i}\t${tpcomplex}\t${fpcomplex}\t${fncomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
   echo -e "${i}\t${tpnotcomplex}\t${fpnotcomplex}\t${fnnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
done


#------------------------------------------------------------------------------#
# FEMALES - Chr X PARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="females"
chr="chrX_PARs"
ploidy="diploid"
samplelist="NA06985 NA06989 NA07000 NA07037 NA07056 NA10847 NA11830 NA11832 NA11840 NA11892"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${pop}/${sex}/${chr}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

cd ${inputdatapath}

for i in $samplelist
do
   for t in TP FP FN
   do
      # First, convert results files to bed format
      awk '{ print $1"\t"$2-1"\t"$2 }' ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.txt > ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.bed

      # Then subset on whether these calls are in complex, XTR or non complex regions
      # complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex.bed

      # Not complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_not_complex.bed

   done

done

# Then make results table but counting the total number of lines in each file.
cd ${inputdatapath}/by_region

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt

for i in $samplelist
do
   # get counts of TP, FP, and FN for each category
   tpcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_complex.bed | wc -l`
   fpcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_complex.bed | wc -l`
   fncomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_complex.bed | wc -l`

   tpnotcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_not_complex.bed | wc -l`
   fpnotcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_not_complex.bed | wc -l`
   fnnotcomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_not_complex.bed | wc -l`

   # add to results file
   echo -e "${i}\t${tpcomplex}\t${fpcomplex}\t${fncomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
   echo -e "${i}\t${tpnotcomplex}\t${fpnotcomplex}\t${fnnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
done


#------------------------------------------------------------------------------#
# FEMALES - Chr X non-PARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="females"
chr="chrX_nonPARs"
ploidy="diploid"
samplelist="NA06985 NA06989 NA07000 NA07037 NA07056 NA10847 NA11830 NA11832 NA11840 NA11892"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${pop}/${sex}/${chr}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

cd ${inputdatapath}

for i in $samplelist
do
   for t in TP FP FN
   do
      # First, convert results files to bed format
      awk '{ print $1"\t"$2-1"\t"$2 }' ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.txt > ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.bed

      # Then subset on whether these calls are in complex, XTR or non complex regions
      # complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX_amplicons_GRCh38.bed  -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed

      # make sure XTR is not in complex file
      bedtools intersect -a by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed -b ${regionsbed}/chrX_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex.bed

      # XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX_XTR_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_XTR.bed

      # Not complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX_amplicons_GRCh38.bed ${regionsbed}/chrX_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_not_complex.bed

      # non PARs minus XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed

      # amplicons only
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX_amplicons_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_Amplicons.bed

      # non PARs minus XTR and amplicons
      bedtools intersect -a ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed -b ${regionsbed}/chrX_amplicons_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTRminusAmplicons.bed

   done

done

# Then make results table but counting the total number of lines in each file.
cd ${inputdatapath}/by_region

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt

for i in $samplelist
do
   # get counts of TP, FP, and FN for each category
   tpcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_complex.bed | wc -l`
   fpcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_complex.bed | wc -l`
   fncomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_complex.bed | wc -l`

   tpXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_XTR.bed | wc -l`
   fpXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_XTR.bed | wc -l`
   fnXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_XTR.bed | wc -l`

   tpnotcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_not_complex.bed | wc -l`
   fpnotcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_not_complex.bed | wc -l`
   fnnotcomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_not_complex.bed | wc -l`

   tpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTR.bed | wc -l`
   fpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTR.bed | wc -l`
   fnnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTR.bed | wc -l`

   tpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fnnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`

   tpAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_Amplicons.bed | wc -l`
   fpAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_Amplicons.bed | wc -l`
   fnAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_Amplicons.bed | wc -l`

   # add to results file
   echo -e "${i}\t${tpcomplex}\t${fpcomplex}\t${fncomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
   echo -e "${i}\t${tpXTR}\t${fpXTR}\t${fnXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
   echo -e "${i}\t${tpnotcomplex}\t${fpnotcomplex}\t${fnnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
   echo -e "${i}\t${tpnonPARminusXTR}\t${fpnonPARminusXTR}\t${fnnonPARminusXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt
   echo -e "${i}\t${tpnonPARminusXTRminusAmplicons}\t${fpnonPARminusXTRminusAmplicons}\t${fnnonPARminusXTRminusAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt
   echo -e "${i}\t${tpAmplicons}\t${fpAmplicons}\t${fnAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt
done



################################################################################
#                                                                              #
# Default aligned data #
#
################################################################################
#------------------------------------------------------------------------------#
# MALES - Chr 8 #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chr8_autos"
ploidy="diploid"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${pop}/${sex}/default/autos"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

cd ${inputdatapath}

for i in $samplelist
do
   for t in TP FP FN
   do
      # First, convert results files to bed format
      awk '{ print $1"\t"$2-1"\t"$2 }' ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.txt > ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.bed

      # Then subset on whether these calls are in complex, XTR or non complex regions
      # complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chr8.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chr8.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SimpleRepeats.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex.bed

      # Not complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chr8.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chr8.human.GRCh38.Repeats.RepeatMasker.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_not_complex.bed

   done

done

# Then make results table but counting the total number of lines in each file.
cd ${inputdatapath}/by_region

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt

for i in $samplelist
do
   # get counts of TP, FP, and FN for each category
   tpcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_complex.bed | wc -l`
   fpcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_complex.bed | wc -l`
   fncomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_complex.bed | wc -l`

   tpnotcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_not_complex.bed | wc -l`
   fpnotcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_not_complex.bed | wc -l`
   fnnotcomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_not_complex.bed | wc -l`

   # add to results file
   echo -e "${i}\t${tpcomplex}\t${fpcomplex}\t${fncomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
   echo -e "${i}\t${tpnotcomplex}\t${fpnotcomplex}\t${fnnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
done


#------------------------------------------------------------------------------#
# MALES - Chr X PARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chrX_PARs"
ploidy="diploid"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${pop}/${sex}/default/${chr}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

cd ${inputdatapath}

for i in $samplelist
do
   for t in TP FP FN
   do
      # First, convert results files to bed format
      awk '{ print $1"\t"$2-1"\t"$2 }' ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.txt > ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.bed

      # Then subset on whether these calls are in complex, XTR or non complex regions
      # complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex.bed

      # Not complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_not_complex.bed

   done

done

# Then make results table but counting the total number of lines in each file.
cd ${inputdatapath}/by_region

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt

for i in $samplelist
do
   # get counts of TP, FP, and FN for each category
   tpcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_complex.bed | wc -l`
   fpcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_complex.bed | wc -l`
   fncomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_complex.bed | wc -l`

   tpnotcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_not_complex.bed | wc -l`
   fpnotcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_not_complex.bed | wc -l`
   fnnotcomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_not_complex.bed | wc -l`

   # add to results file
   echo -e "${i}\t${tpcomplex}\t${fpcomplex}\t${fncomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
   echo -e "${i}\t${tpnotcomplex}\t${fpnotcomplex}\t${fnnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
done


#------------------------------------------------------------------------------#
# MALES - Chr X non-PARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chrX_nonPARs"
ploidy="haploid"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${pop}/${sex}/default/${chr}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

cd ${inputdatapath}

for i in $samplelist
do
   for t in TP FP FN
   do
      # First, convert results files to bed format
      awk '{ print $1"\t"$2-1"\t"$2 }' ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.txt > ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.bed

      # Then subset on whether these calls are in complex, XTR or non complex regions
      # complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX_amplicons_GRCh38.bed  -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed

      # make sure XTR is not in complex file
      bedtools intersect -a by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed -b ${regionsbed}/chrX_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex.bed

      # XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX_XTR_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_XTR.bed

      # Not complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX_amplicons_GRCh38.bed ${regionsbed}/chrX_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_not_complex.bed

      # non PARs minus XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed

      # amplicons only
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX_amplicons_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_Amplicons.bed

      # non PARs minus XTR and amplicons
      bedtools intersect -a ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed -b ${regionsbed}/chrX_amplicons_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTRminusAmplicons.bed
   done

done

# Then make results table but counting the total number of lines in each file.
cd ${inputdatapath}/by_region

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt

for i in $samplelist
do
   # get counts of TP, FP, and FN for each category
   tpcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_complex.bed | wc -l`
   fpcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_complex.bed | wc -l`
   fncomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_complex.bed | wc -l`

   tpXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_XTR.bed | wc -l`
   fpXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_XTR.bed | wc -l`
   fnXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_XTR.bed | wc -l`

   tpnotcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_not_complex.bed | wc -l`
   fpnotcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_not_complex.bed | wc -l`
   fnnotcomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_not_complex.bed | wc -l`

   tpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTR.bed | wc -l`
   fpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTR.bed | wc -l`
   fnnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTR.bed | wc -l`

   tpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fnnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`

   tpAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_Amplicons.bed | wc -l`
   fpAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_Amplicons.bed | wc -l`
   fnAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_Amplicons.bed | wc -l`

   # add to results file
   echo -e "${i}\t${tpcomplex}\t${fpcomplex}\t${fncomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
   echo -e "${i}\t${tpXTR}\t${fpXTR}\t${fnXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
   echo -e "${i}\t${tpnotcomplex}\t${fpnotcomplex}\t${fnnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
   echo -e "${i}\t${tpnonPARminusXTR}\t${fpnonPARminusXTR}\t${fnnonPARminusXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt
   echo -e "${i}\t${tpnonPARminusXTRminusAmplicons}\t${fpnonPARminusXTRminusAmplicons}\t${fnnonPARminusXTRminusAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt
   echo -e "${i}\t${tpAmplicons}\t${fpAmplicons}\t${fnAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt
done



#------------------------------------------------------------------------------#
# MALES - Chr X non-PARs #
# DIPLOID
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chrX_nonPARs"
ploidy="diploid"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${pop}/${sex}/default/${chr}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

cd ${inputdatapath}

for i in $samplelist
do
   for t in TP FP FN
   do
      # First, convert results files to bed format
      awk '{ print $1"\t"$2-1"\t"$2 }' ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.txt > ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.bed

      # Then subset on whether these calls are in complex, XTR or non complex regions
      # complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX_amplicons_GRCh38.bed  -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed

      # make sure XTR is not in complex file
      bedtools intersect -a by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed -b ${regionsbed}/chrX_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex.bed

      # XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX_XTR_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_XTR.bed

      # Not complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX_amplicons_GRCh38.bed ${regionsbed}/chrX_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_not_complex.bed

      # non PARs minus XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed

      # amplicons only
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX_amplicons_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_Amplicons.bed

      # non PARs minus XTR and amplicons
      bedtools intersect -a ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed -b ${regionsbed}/chrX_amplicons_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTRminusAmplicons.bed
   done

done

# Then make results table but counting the total number of lines in each file.
cd ${inputdatapath}/by_region

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt

for i in $samplelist
do
   # get counts of TP, FP, and FN for each category
   tpcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_complex.bed | wc -l`
   fpcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_complex.bed | wc -l`
   fncomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_complex.bed | wc -l`

   tpXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_XTR.bed | wc -l`
   fpXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_XTR.bed | wc -l`
   fnXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_XTR.bed | wc -l`

   tpnotcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_not_complex.bed | wc -l`
   fpnotcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_not_complex.bed | wc -l`
   fnnotcomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_not_complex.bed | wc -l`

   tpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTR.bed | wc -l`
   fpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTR.bed | wc -l`
   fnnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTR.bed | wc -l`

   tpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fnnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`

   tpAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_Amplicons.bed | wc -l`
   fpAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_Amplicons.bed | wc -l`
   fnAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_Amplicons.bed | wc -l`

   # add to results file
   echo -e "${i}\t${tpcomplex}\t${fpcomplex}\t${fncomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
   echo -e "${i}\t${tpXTR}\t${fpXTR}\t${fnXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
   echo -e "${i}\t${tpnotcomplex}\t${fpnotcomplex}\t${fnnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
   echo -e "${i}\t${tpnonPARminusXTR}\t${fpnonPARminusXTR}\t${fnnonPARminusXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt
   echo -e "${i}\t${tpnonPARminusXTRminusAmplicons}\t${fpnonPARminusXTRminusAmplicons}\t${fnnonPARminusXTRminusAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt
   echo -e "${i}\t${tpAmplicons}\t${fpAmplicons}\t${fnAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt
done


#------------------------------------------------------------------------------#
# MALES - Chr Y non-PARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chrY"
ploidy="haploid"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${pop}/${sex}/default/${chr}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"


cd ${inputdatapath}

for i in $samplelist
do
   for t in TP FP FN
   do
      # First, convert results files to bed format
      awk '{ print $1"\t"$2-1"\t"$2 }' ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.txt > ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.bed

      # Then subset on whether these calls are in complex, XTR or non complex regions
      # complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrY.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrY_amplicons_GRCh38.bed  -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed

      # make sure XTR is not in complex file
      bedtools intersect -a by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed -b ${regionsbed}/chrY_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex.bed

      # XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY_XTR_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_XTR.bed

      # Not complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrY.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrY_amplicons_GRCh38.bed ${regionsbed}/chrY_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_not_complex.bed

      # non PARs minus XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed

      # amplicons only
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY_amplicons_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_Amplicons.bed

      # non PARs minus XTR and amplicons
      bedtools intersect -a ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed -b ${regionsbed}/chrY_amplicons_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTRminusAmplicons.bed
   done

done

# Then make results table but counting the total number of lines in each file.
cd ${inputdatapath}/by_region

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt

for i in $samplelist
do
   # get counts of TP, FP, and FN for each category
   tpcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_complex.bed | wc -l`
   fpcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_complex.bed | wc -l`
   fncomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_complex.bed | wc -l`

   tpXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_XTR.bed | wc -l`
   fpXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_XTR.bed | wc -l`
   fnXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_XTR.bed | wc -l`

   tpnotcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_not_complex.bed | wc -l`
   fpnotcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_not_complex.bed | wc -l`
   fnnotcomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_not_complex.bed | wc -l`

   tpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTR.bed | wc -l`
   fpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTR.bed | wc -l`
   fnnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTR.bed | wc -l`

   tpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fnnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`

   tpAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_Amplicons.bed | wc -l`
   fpAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_Amplicons.bed | wc -l`
   fnAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_Amplicons.bed | wc -l`

   # add to results file
   echo -e "${i}\t${tpcomplex}\t${fpcomplex}\t${fncomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
   echo -e "${i}\t${tpXTR}\t${fpXTR}\t${fnXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
   echo -e "${i}\t${tpnotcomplex}\t${fpnotcomplex}\t${fnnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
   echo -e "${i}\t${tpnonPARminusXTR}\t${fpnonPARminusXTR}\t${fnnonPARminusXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt
   echo -e "${i}\t${tpnonPARminusXTRminusAmplicons}\t${fpnonPARminusXTRminusAmplicons}\t${fnnonPARminusXTRminusAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt
   echo -e "${i}\t${tpAmplicons}\t${fpAmplicons}\t${fnAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt

done


#------------------------------------------------------------------------------#
# MALES - Chr Y non-PARs #
# DIPLOID
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chrY"
ploidy="diploid"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${pop}/${sex}/default/${chr}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"


cd ${inputdatapath}

for i in $samplelist
do
   for t in TP FP FN
   do
      # First, convert results files to bed format
      awk '{ print $1"\t"$2-1"\t"$2 }' ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.txt > ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.bed

      # Then subset on whether these calls are in complex, XTR or non complex regions
      # complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrY.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrY_amplicons_GRCh38.bed  -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed

      # make sure XTR is not in complex file
      bedtools intersect -a by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed -b ${regionsbed}/chrY_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex.bed

      # XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY_XTR_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_XTR.bed

      # Not complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrY.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrY_amplicons_GRCh38.bed ${regionsbed}/chrY_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_not_complex.bed

      # non PARs minus XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed

      # amplicons only
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrY_amplicons_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_Amplicons.bed

      # non PARs minus XTR and amplicons
      bedtools intersect -a ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed -b ${regionsbed}/chrY_amplicons_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTRminusAmplicons.bed

   done

done

# Then make results table but counting the total number of lines in each file.
cd ${inputdatapath}/by_region

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt

for i in $samplelist
do
   # get counts of TP, FP, and FN for each category
   tpcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_complex.bed | wc -l`
   fpcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_complex.bed | wc -l`
   fncomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_complex.bed | wc -l`

   tpXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_XTR.bed | wc -l`
   fpXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_XTR.bed | wc -l`
   fnXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_XTR.bed | wc -l`

   tpnotcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_not_complex.bed | wc -l`
   fpnotcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_not_complex.bed | wc -l`
   fnnotcomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_not_complex.bed | wc -l`

   tpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTR.bed | wc -l`
   fpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTR.bed | wc -l`
   fnnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTR.bed | wc -l`

   tpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fnnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`

   tpAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_Amplicons.bed | wc -l`
   fpAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_Amplicons.bed | wc -l`
   fnAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_Amplicons.bed | wc -l`

   # add to results file
   echo -e "${i}\t${tpcomplex}\t${fpcomplex}\t${fncomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
   echo -e "${i}\t${tpXTR}\t${fpXTR}\t${fnXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
   echo -e "${i}\t${tpnotcomplex}\t${fpnotcomplex}\t${fnnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
   echo -e "${i}\t${tpnonPARminusXTR}\t${fpnonPARminusXTR}\t${fnnonPARminusXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt
   echo -e "${i}\t${tpnonPARminusXTRminusAmplicons}\t${fpnonPARminusXTRminusAmplicons}\t${fnnonPARminusXTRminusAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt
   echo -e "${i}\t${tpAmplicons}\t${fpAmplicons}\t${fnAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt

done



#------------------------------------------------------------------------------#
# FEMALES - Chr 8 #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="females"
chr="chr8_autos"
ploidy="diploid"
samplelist="NA06985 NA06989 NA07000 NA07037 NA07056 NA10847 NA11830 NA11832 NA11840 NA11892"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${pop}/${sex}/default/autos"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

cd ${inputdatapath}

for i in $samplelist
do
   for t in TP FP FN
   do
      # First, convert results files to bed format
      awk '{ print $1"\t"$2-1"\t"$2 }' ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.txt > ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.bed

      # Then subset on whether these calls are in complex, XTR or non complex regions
      # complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chr8.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chr8.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SimpleRepeats.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex.bed

      # Not complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chr8.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chr8.human.GRCh38.Repeats.RepeatMasker.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_not_complex.bed

   done

done

# Then make results table but counting the total number of lines in each file.
cd ${inputdatapath}/by_region

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt

for i in $samplelist
do
   # get counts of TP, FP, and FN for each category
   tpcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_complex.bed | wc -l`
   fpcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_complex.bed | wc -l`
   fncomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_complex.bed | wc -l`

   tpnotcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_not_complex.bed | wc -l`
   fpnotcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_not_complex.bed | wc -l`
   fnnotcomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_not_complex.bed | wc -l`

   # add to results file
   echo -e "${i}\t${tpcomplex}\t${fpcomplex}\t${fncomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
   echo -e "${i}\t${tpnotcomplex}\t${fpnotcomplex}\t${fnnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
done



#------------------------------------------------------------------------------#
# FEMALES - Chr X PARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="females"
chr="chrX_PARs"
ploidy="diploid"
samplelist="NA06985 NA06989 NA07000 NA07037 NA07056 NA10847 NA11830 NA11832 NA11840 NA11892"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${pop}/${sex}/default/${chr}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

cd ${inputdatapath}

for i in $samplelist
do
   for t in TP FP FN
   do
      # First, convert results files to bed format
      awk '{ print $1"\t"$2-1"\t"$2 }' ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.txt > ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.bed

      # Then subset on whether these calls are in complex, XTR or non complex regions
      # complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex.bed

      # Not complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_not_complex.bed

   done

done

# Then make results table but counting the total number of lines in each file.
cd ${inputdatapath}/by_region

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt

for i in $samplelist
do
   # get counts of TP, FP, and FN for each category
   tpcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_complex.bed | wc -l`
   fpcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_complex.bed | wc -l`
   fncomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_complex.bed | wc -l`

   tpnotcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_not_complex.bed | wc -l`
   fpnotcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_not_complex.bed | wc -l`
   fnnotcomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_not_complex.bed | wc -l`

   # add to results file
   echo -e "${i}\t${tpcomplex}\t${fpcomplex}\t${fncomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
   echo -e "${i}\t${tpnotcomplex}\t${fpnotcomplex}\t${fnnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
done


#------------------------------------------------------------------------------#
# FEMALES - Chr X non-PARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="females"
chr="chrX_nonPARs"
ploidy="diploid"
samplelist="NA06985 NA06989 NA07000 NA07037 NA07056 NA10847 NA11830 NA11832 NA11840 NA11892"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${pop}/${sex}/default/${chr}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

cd ${inputdatapath}

for i in $samplelist
do
   for t in TP FP FN
   do
      # First, convert results files to bed format
      awk '{ print $1"\t"$2-1"\t"$2 }' ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.txt > ${inputdatapath}/${i}_${chr}_golden_vs_called_${t}_pos.bed

      # Then subset on whether these calls are in complex, XTR or non complex regions
      # complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX_amplicons_GRCh38.bed  -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed

      # make sure XTR is not in complex file
      bedtools intersect -a by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.bed -b ${regionsbed}/chrX_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_complex.bed

      # XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX_XTR_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_XTR.bed

      # Not complex
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX_amplicons_GRCh38.bed ${regionsbed}/chrX_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_not_complex.bed

      # non PARs minus XTR
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX_XTR_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed

      # amplicons only
      bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos.bed -b ${regionsbed}/chrX_amplicons_GRCh38.bed -u > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_Amplicons.bed

      # non PARs minus XTR and amplicons
      bedtools intersect -a ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTR.bed -b ${regionsbed}/chrX_amplicons_GRCh38.bed -v > ${inputdatapath}/by_region/${i}_${chr}_golden_vs_called_${t}_pos_nonPARsminusXTRminusAmplicons.bed
   done

done

# Then make results table but counting the total number of lines in each file.
cd ${inputdatapath}/by_region

echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt
echo -e "Sample\tTP\tFP\tFN" > ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt

for i in $samplelist
do
   # get counts of TP, FP, and FN for each category
   tpcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_complex.bed | wc -l`
   fpcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_complex.bed | wc -l`
   fncomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_complex.bed | wc -l`

   tpXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_XTR.bed | wc -l`
   fpXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_XTR.bed | wc -l`
   fnXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_XTR.bed | wc -l`

   tpnotcomplex=`cat ${i}_${chr}_golden_vs_called_TP_pos_not_complex.bed | wc -l`
   fpnotcomplex=`cat ${i}_${chr}_golden_vs_called_FP_pos_not_complex.bed | wc -l`
   fnnotcomplex=`cat ${i}_${chr}_golden_vs_called_FN_pos_not_complex.bed | wc -l`

   tpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTR.bed | wc -l`
   fpnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTR.bed | wc -l`
   fnnonPARminusXTR=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTR.bed | wc -l`

   tpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fpnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`
   fnnonPARminusXTRminusAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_nonPARsminusXTRminusAmplicons.bed | wc -l`

   tpAmplicons=`cat ${i}_${chr}_golden_vs_called_TP_pos_Amplicons.bed | wc -l`
   fpAmplicons=`cat ${i}_${chr}_golden_vs_called_FP_pos_Amplicons.bed | wc -l`
   fnAmplicons=`cat ${i}_${chr}_golden_vs_called_FN_pos_Amplicons.bed | wc -l`

   # add to results file
   echo -e "${i}\t${tpcomplex}\t${fpcomplex}\t${fncomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_complex.txt
   echo -e "${i}\t${tpXTR}\t${fpXTR}\t${fnXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_XTR.txt
   echo -e "${i}\t${tpnotcomplex}\t${fpnotcomplex}\t${fnnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_not_complex.txt
   echo -e "${i}\t${tpnonPARminusXTR}\t${fpnonPARminusXTR}\t${fnnonPARminusXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTR.txt
   echo -e "${i}\t${tpnonPARminusXTRminusAmplicons}\t${fpnonPARminusXTRminusAmplicons}\t${fnnonPARminusXTRminusAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt
   echo -e "${i}\t${tpAmplicons}\t${fpAmplicons}\t${fnAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_golden_vs_called_performance_metrics_Amplicons.txt
done
