#!/bin/bash
#SBATCH --job-name=simulated_variants_by_region # Job name
#SBATCH -o slurm.simulated_variants_by_region.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.simulated_variants_by_region.err                # STDERR (%j = JobId)
#SBATCH -n 1
#SBATCH -t 7-00:00:00
#SBATCH -q tempboost
#SBATCH --mail-user=amtarave@asu.edu # send-to address
#SBATCH --mail-type=ALL

# For this I want to get the number of simulated variants per individual in the
# following regions:
# 1) PARs, nonPARs (minus XTR), XTR
# 2) Complex PARs, not complex PARs, complex nonPARs (minus XTR),
#    not complex nonPARs (minus XTR), XTR
# And for chr8 and chrM

cd /scratch/amtarave/variant_calling_simulations_project/performance_metrics_by_region/simulated_variants_by_region


#------------------------------------------------------------------------------#
# MALES - Chr X PARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chrX_PARs"
ploidy="diploid"
coverage="20"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/simulations/${pop}/${sex}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"


echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt

for i in $samplelist
do
    # all of PARs
    sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf | wc -l`

    # Then subset on whether these calls are in complex, XTR or non complex regions
    # complex
    simcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed -u | wc -l`

    # Not complex
    simnotcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed -v | wc -l`

    # add to results file
    echo -e "${i}\t${sim}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
    echo -e "${i}\t${simcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
    echo -e "${i}\t${simnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt

done


#------------------------------------------------------------------------------#
# MALES - Chr X nonPARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chrX_nonPARs"
ploidy="haploid"
coverage="10"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/simulations/${pop}/${sex}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"


#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_XTR.txt
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_nonPARsminusXTR.txt

echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_Amplicons.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_nonPARsminusXTRminusAmplicons.txt

for i in $samplelist
do
    # all of non PARs
    #sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf | wc -l`

    # Then subset on whether these calls are in complex, XTR or non complex regions
    # complex
    #bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX_amplicons_GRCh38.bed -u -header > ${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.vcf

    # make sure XTR is not in complex file
    #simcomplex=`bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.vcf -b ${regionsbed}/chrX_XTR_GRCh38.bed -v | wc -l`

    # remove intermediate file
    #rm ${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.vcf

    # XTR
    #simXTR=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX_XTR_GRCh38.bed -u | wc -l`

    # Not complex
    #simnotcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX_amplicons_GRCh38.bed ${regionsbed}/chrX_XTR_GRCh38.bed -v | wc -l`

    # non PARs minus XTR
    #simnonPARsminusXTR=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX_XTR_GRCh38.bed -v | wc -l`

    # Amplicons
    simAmplicons=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX_amplicons_GRCh38.bed -u | wc -l`

    # non PARs minus XTR minus amplicons
    simnonPARsminusXTRminusAmplicons=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX_XTR_GRCh38.bed ${regionsbed}/chrX_amplicons_GRCh38.bed -v | wc -l`

    # add to results file
    #echo -e "${i}\t${sim}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
    #echo -e "${i}\t${simcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
    #echo -e "${i}\t${simnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
    #echo -e "${i}\t${simXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_XTR.txt
    #echo -e "${i}\t${simnonPARsminusXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_nonPARsminusXTR.txt

    echo -e "${i}\t${simAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_Amplicons.txt
    echo -e "${i}\t${simnonPARsminusXTRminusAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_nonPARsminusXTRminusAmplicons.txt

done


#------------------------------------------------------------------------------#
# MALES - Chr Y non-PARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chrY"
ploidy="haploid"
coverage="10"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/simulations/${pop}/${sex}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"


#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_XTR.txt
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_nonPARsminusXTR.txt

echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_Amplicons.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_nonPARsminusXTRminusAmplicons.txt


for i in $samplelist
do
    # all of non PARs
    #sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf | wc -l`

    # Then subset on whether these calls are in complex, XTR or non complex regions
    # complex
    #bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrY.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrY.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrY_amplicons_GRCh38.bed -u -header > ${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.vcf

    # make sure XTR is not in complex file
    #simcomplex=`bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.vcf -b ${regionsbed}/chrY_XTR_GRCh38.bed -v | wc -l`

    # remove intermediate file
    #rm ${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.vcf

    # XTR
    #simXTR=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrY_XTR_GRCh38.bed -u | wc -l`

    # Not complex
    #simnotcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrY.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrY.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrY_amplicons_GRCh38.bed ${regionsbed}/chrY_XTR_GRCh38.bed -v | wc -l`

    # non PARs minus XTR
    #simnonPARsminusXTR=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrY_XTR_GRCh38.bed -v | wc -l`

    # Amplicons
    simAmplicons=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrY_amplicons_GRCh38.bed -u | wc -l`

    # non PARs minus XTR minus amplicons
    simnonPARsminusXTRminusAmplicons=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrY_XTR_GRCh38.bed ${regionsbed}/chrY_amplicons_GRCh38.bed -v | wc -l`

    # add to results file
    #echo -e "${i}\t${sim}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
    #echo -e "${i}\t${simcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
    #echo -e "${i}\t${simnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
    #echo -e "${i}\t${simXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_XTR.txt
    #echo -e "${i}\t${simnonPARsminusXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_nonPARsminusXTR.txt

    echo -e "${i}\t${simAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_Amplicons.txt
    echo -e "${i}\t${simnonPARsminusXTRminusAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_nonPARsminusXTRminusAmplicons.txt

done

# TO DO - Add for loop to the chr8 ones for all autosomes..no need to do complex
# and not complex
#------------------------------------------------------------------------------#
# MALES - Chr 8 #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chr8"
ploidy="diploid"
coverage="20"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/simulations/${pop}/${sex}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

# Havent separated out chr 8 complex and not complex but maybe will do this?
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt

for i in $samplelist
do
    # all of PARs
    sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf | wc -l`

    # Then subset on whether these calls are in complex, XTR or non complex regions
    # complex
    simcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chr8.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chr8.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SimpleRepeats.bed -u | wc -l`

    # Not complex
    simnotcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chr8.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chr8.human.GRCh38.Repeats.RepeatMasker.bed -v | wc -l`

    # add to results file
    #echo -e "${i}\t${sim}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
    echo -e "${i}\t${simcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
    echo -e "${i}\t${simnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt

done

#---------------#
# ALL AUTOSOMES #
#---------------#
for t in {1..22}
do
    echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_chr${t}_${ploidy}_simulated_counts.txt
    for i in $samplelist
    do
          sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_chr${t}_golden.vcf | wc -l`
          echo -e "${i}\t${sim}" >> ${pop}_${sex}_chr${t}_${ploidy}_simulated_counts.txt
    done

  echo "chr${t}"
done


#------------------------------------------------------------------------------#
# MALES - Chr M #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chrM"
ploidy="haploid"
coverage="20"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/simulations/${pop}/${sex}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

# Havent separated out chr 8 complex and not complex but maybe will do this?
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt

for i in $samplelist
do
    # all of PARs
    sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf | wc -l`

    # Then subset on whether these calls are in complex, XTR or non complex regions
    # complex
    #simcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed -u | wc -l`

    # Not complex
    #simnotcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed -v | wc -l`

    # add to results file
    echo -e "${i}\t${sim}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
    #echo -e "${i}\t${simcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
    #echo -e "${i}\t${simnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt

done



#------------------------------------------------------------------------------#
# FEMALES - Chr X PARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="females"
chr="chrX_PARs"
ploidy="diploid"
coverage="20"
samplelist="NA06985 NA06989 NA07000 NA07037 NA07056 NA10847 NA11830 NA11832 NA11840 NA11892"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/simulations/${pop}/${sex}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"


echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt

for i in $samplelist
do
    # all of PARs
    sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf | wc -l`

    # Then subset on whether these calls are in complex, XTR or non complex regions
    # complex
    simcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed -u | wc -l`

    # Not complex
    simnotcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed -v | wc -l`

    # add to results file
    echo -e "${i}\t${sim}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
    echo -e "${i}\t${simcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
    echo -e "${i}\t${simnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt

done


#------------------------------------------------------------------------------#
# FEMALES - Chr X nonPARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="females"
chr="chrX_nonPARs"
ploidy="diploid"
coverage="20"
samplelist="NA06985 NA06989 NA07000 NA07037 NA07056 NA10847 NA11830 NA11832 NA11840 NA11892"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/simulations/${pop}/${sex}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"


#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_XTR.txt
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_nonPARsminusXTR.txt

echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_Amplicons.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_nonPARsminusXTRminusAmplicons.txt

for i in $samplelist
do
    # all of non PARs
    #sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf | wc -l`

    # Then subset on whether these calls are in complex, XTR or non complex regions
    # complex
    #bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX_amplicons_GRCh38.bed -u -header > ${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.vcf

    # make sure XTR is not in complex file
    #simcomplex=`bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.vcf -b ${regionsbed}/chrX_XTR_GRCh38.bed -v | wc -l`

    # remove intermediate file
    #rm ${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.vcf

    # XTR
    #simXTR=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX_XTR_GRCh38.bed -u | wc -l`

    # Not complex
    #simnotcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX_amplicons_GRCh38.bed ${regionsbed}/chrX_XTR_GRCh38.bed -v | wc -l`

    # non PARs minus XTR
    #simnonPARsminusXTR=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX_XTR_GRCh38.bed -v | wc -l`

    # Amplicons
    simAmplicons=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX_amplicons_GRCh38.bed -u | wc -l`

    # non PARs minus XTR minus amplicons
    simnonPARsminusXTRminusAmplicons=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX_XTR_GRCh38.bed ${regionsbed}/chrX_amplicons_GRCh38.bed -v | wc -l`

    # add to results file
    #echo -e "${i}\t${sim}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
    #echo -e "${i}\t${simcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
    #echo -e "${i}\t${simnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
    #echo -e "${i}\t${simXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_XTR.txt
    #echo -e "${i}\t${simnonPARsminusXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_nonPARsminusXTR.txt

    echo -e "${i}\t${simAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_Amplicons.txt
    echo -e "${i}\t${simnonPARsminusXTRminusAmplicons}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_nonPARsminusXTRminusAmplicons.txt

done

#------------------------------------------------------------------------------#
# FEMALES - Chr 8 #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="females"
chr="chr8"
ploidy="diploid"
coverage="20"
samplelist="NA06985 NA06989 NA07000 NA07037 NA07056 NA10847 NA11830 NA11832 NA11840 NA11892"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/simulations/${pop}/${sex}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

# Havent separated out chr 8 complex and not complex but maybe will do this?
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt

#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt

for i in $samplelist
do
    # all of PARs
    sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf | wc -l`

    # Then subset on whether these calls are in complex, XTR or non complex regions
    # complex
    simcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chr8.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chr8.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SimpleRepeats.bed -u | wc -l`

    # Not complex
    simnotcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chr8.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chr8.human.GRCh38.Repeats.RepeatMasker.bed -v | wc -l`

    # add to results file
    #echo -e "${i}\t${sim}"
    #echo -e "${i}\t${sim}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt

    echo -e "${i}\t${simcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
    echo -e "${i}\t${simnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt

done


#---------------#
# ALL AUTOSOMES #
#---------------#
for t in {1..22}
do
    echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_chr${t}_${ploidy}_simulated_counts.txt
    for i in $samplelist
    do
          sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_chr${t}_golden.vcf | wc -l`
          echo -e "${i}\t${sim}" >> ${pop}_${sex}_chr${t}_${ploidy}_simulated_counts.txt
    done

  echo "chr${t}"
done


#------------------------------------------------------------------------------#
# FEMALES - Chr M #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="females"
chr="chrM"
ploidy="haploid"
coverage="20"
samplelist="NA06985 NA06989 NA07000 NA07037 NA07056 NA10847 NA11830 NA11832 NA11840 NA11892"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/simulations/${pop}/${sex}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

# Havent separated out chr 8 complex and not complex but maybe will do this?
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt

for i in $samplelist
do
    # all of PARs
    sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf | wc -l`

    # Then subset on whether these calls are in complex, XTR or non complex regions
    # complex
    #simcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed -u | wc -l`

    # Not complex
    #simnotcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed -v | wc -l`

    # add to results file
    echo -e "${i}\t${sim}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
    #echo -e "${i}\t${simcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
    #echo -e "${i}\t${simnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt

done






# HAVENT UPDATED WITH AMPLICON AND NONPARS MINUS XTR AND AMPLICONS #
################################################################################
################################################################################
# 20 Samples #
################################################################################
################################################################################

cd /scratch/amtarave/variant_calling_simulations_project/performance_metrics_by_region/simulated_variants_by_region/20_samples


#------------------------------------------------------------------------------#
# MALES - Chr X PARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chrX_PARs"
ploidy="diploid"
coverage="20"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831 NA11843 NA11881 NA11893 NA11919 NA11930 NA11932 NA11992 NA11994 NA12003 NA12005"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/simulations/${pop}/${sex}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"


echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt

for i in $samplelist
do
    # all of PARs
    sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf | wc -l`

    # Then subset on whether these calls are in complex, XTR or non complex regions
    # complex
    simcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed -u | wc -l`

    # Not complex
    simnotcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed -v | wc -l`

    # add to results file
    echo -e "${i}\t${sim}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
    echo -e "${i}\t${simcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
    echo -e "${i}\t${simnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt

done


#------------------------------------------------------------------------------#
# MALES - Chr X nonPARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chrX_nonPARs"
ploidy="haploid"
coverage="10"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831 NA11843 NA11881 NA11893 NA11919 NA11930 NA11932 NA11992 NA11994 NA12003 NA12005"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/simulations/${pop}/${sex}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"


echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_XTR.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_nonPARsminusXTR.txt

for i in $samplelist
do
    # all of non PARs
    sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf | wc -l`

    # Then subset on whether these calls are in complex, XTR or non complex regions
    # complex
    bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX_amplicons_GRCh38.bed -u -header > ${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.vcf

    # make sure XTR is not in complex file
    simcomplex=`bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.vcf -b ${regionsbed}/chrX_XTR_GRCh38.bed -v | wc -l`

    # remove intermediate file
    rm ${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.vcf

    # XTR
    simXTR=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX_XTR_GRCh38.bed -u | wc -l`

    # Not complex
    simnotcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX_amplicons_GRCh38.bed ${regionsbed}/chrX_XTR_GRCh38.bed -v | wc -l`

    # non PARs minus XTR
    simnonPARsminusXTR=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX_XTR_GRCh38.bed -v | wc -l`

    # add to results file
    echo -e "${i}\t${sim}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
    echo -e "${i}\t${simcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
    echo -e "${i}\t${simnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
    echo -e "${i}\t${simXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_XTR.txt
    echo -e "${i}\t${simnonPARsminusXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_nonPARsminusXTR.txt

done


#------------------------------------------------------------------------------#
# MALES - Chr Y non-PARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chrY"
ploidy="haploid"
coverage="10"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831 NA11843 NA11881 NA11893 NA11919 NA11930 NA11932 NA11992 NA11994 NA12003 NA12005"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/simulations/${pop}/${sex}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"


echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_XTR.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_nonPARsminusXTR.txt

for i in $samplelist
do
    # all of non PARs
    sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf | wc -l`

    # Then subset on whether these calls are in complex, XTR or non complex regions
    # complex
    bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrY.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrY.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrY_amplicons_GRCh38.bed -u -header > ${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.vcf

    # make sure XTR is not in complex file
    simcomplex=`bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.vcf -b ${regionsbed}/chrY_XTR_GRCh38.bed -v | wc -l`

    # remove intermediate file
    rm ${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.vcf

    # XTR
    simXTR=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrY_XTR_GRCh38.bed -u | wc -l`

    # Not complex
    simnotcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrY.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrY.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrY.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrY_amplicons_GRCh38.bed ${regionsbed}/chrY_XTR_GRCh38.bed -v | wc -l`

    # non PARs minus XTR
    simnonPARsminusXTR=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrY_XTR_GRCh38.bed -v | wc -l`

    # add to results file
    echo -e "${i}\t${sim}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
    echo -e "${i}\t${simcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
    echo -e "${i}\t${simnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
    echo -e "${i}\t${simXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_XTR.txt
    echo -e "${i}\t${simnonPARsminusXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_nonPARsminusXTR.txt

done

# TO DO - Add for loop to the chr8 ones for all autosomes..no need to do complex
# and not complex
#------------------------------------------------------------------------------#
# MALES - Chr 8 #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chr8"
ploidy="diploid"
coverage="20"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831 NA11843 NA11881 NA11893 NA11919 NA11930 NA11932 NA11992 NA11994 NA12003 NA12005"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/simulations/${pop}/${sex}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

# Havent separated out chr 8 complex and not complex but maybe will do this?
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt

for i in $samplelist
do
    # all of PARs
    sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf | wc -l`

    # Then subset on whether these calls are in complex, XTR or non complex regions
    # complex
    simcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chr8.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chr8.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SimpleRepeats.bed -u | wc -l`

    # Not complex
    simnotcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chr8.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chr8.human.GRCh38.Repeats.RepeatMasker.bed -v | wc -l`

    # add to results file
    #echo -e "${i}\t${sim}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
    echo -e "${i}\t${simcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
    echo -e "${i}\t${simnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt

done

#---------------#
# ALL AUTOSOMES #
#---------------#
# DONT RUN, HAVENT RAN ALL AUTOSOMES FOR THIS
#for t in {1..22}
for t in {8..8}
do
    echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_chr${t}_${ploidy}_simulated_counts.txt
    for i in $samplelist
    do
          sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_chr${t}_golden.vcf | wc -l`
          echo -e "${i}\t${sim}" >> ${pop}_${sex}_chr${t}_${ploidy}_simulated_counts.txt
    done

  echo "chr${t}"
done


#------------------------------------------------------------------------------#
# MALES - Chr M #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="males"
chr="chrM"
ploidy="haploid"
coverage="20"
samplelist="NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831 NA11843 NA11881 NA11893 NA11919 NA11930 NA11932 NA11992 NA11994 NA12003 NA12005"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/simulations/${pop}/${sex}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

# Havent separated out chr 8 complex and not complex but maybe will do this?
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt

for i in $samplelist
do
    # all of PARs
    sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf | wc -l`

    # Then subset on whether these calls are in complex, XTR or non complex regions
    # complex
    #simcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed -u | wc -l`

    # Not complex
    #simnotcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed -v | wc -l`

    # add to results file
    echo -e "${i}\t${sim}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
    #echo -e "${i}\t${simcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
    #echo -e "${i}\t${simnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt

done



#------------------------------------------------------------------------------#
# FEMALES - Chr X PARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="females"
chr="chrX_PARs"
ploidy="diploid"
coverage="20"
samplelist="NA06985 NA06989 NA07000 NA07037 NA07056 NA10847 NA11830 NA11832 NA11840 NA11892 NA11894 NA11918 NA11920 NA11931 NA11933 NA11995 NA12004 NA12006 NA12044 NA12046"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/simulations/${pop}/${sex}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"


echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt

for i in $samplelist
do
    # all of PARs
    sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf | wc -l`

    # Then subset on whether these calls are in complex, XTR or non complex regions
    # complex
    simcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed -u | wc -l`

    # Not complex
    simnotcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed -v | wc -l`

    # add to results file
    echo -e "${i}\t${sim}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
    echo -e "${i}\t${simcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
    echo -e "${i}\t${simnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt

done


#------------------------------------------------------------------------------#
# FEMALES - Chr X nonPARs #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="females"
chr="chrX_nonPARs"
ploidy="diploid"
coverage="20"
samplelist="NA06985 NA06989 NA07000 NA07037 NA07056 NA10847 NA11830 NA11832 NA11840 NA11892 NA11894 NA11918 NA11920 NA11931 NA11933 NA11995 NA12004 NA12006 NA12044 NA12046"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/simulations/${pop}/${sex}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"


echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_XTR.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_nonPARsminusXTR.txt

for i in $samplelist
do
    # all of non PARs
    sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf | wc -l`

    # Then subset on whether these calls are in complex, XTR or non complex regions
    # complex
    bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX_amplicons_GRCh38.bed -u -header > ${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.vcf

    # make sure XTR is not in complex file
    simcomplex=`bedtools intersect -a ${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.vcf -b ${regionsbed}/chrX_XTR_GRCh38.bed -v | wc -l`

    # remove intermediate file
    rm ${i}_${chr}_golden_vs_called_${t}_pos_complex_intermediate.vcf

    # XTR
    simXTR=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX_XTR_GRCh38.bed -u | wc -l`

    # Not complex
    simnotcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX_amplicons_GRCh38.bed ${regionsbed}/chrX_XTR_GRCh38.bed -v | wc -l`

    # non PARs minus XTR
    simnonPARsminusXTR=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX_XTR_GRCh38.bed -v | wc -l`

    # add to results file
    echo -e "${i}\t${sim}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
    echo -e "${i}\t${simcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
    echo -e "${i}\t${simnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
    echo -e "${i}\t${simXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_XTR.txt
    echo -e "${i}\t${simnonPARsminusXTR}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_nonPARsminusXTR.txt

done

#------------------------------------------------------------------------------#
# FEMALES - Chr 8 #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="females"
chr="chr8"
ploidy="diploid"
coverage="20"
samplelist="NA06985 NA06989 NA07000 NA07037 NA07056 NA10847 NA11830 NA11832 NA11840 NA11892 NA11894 NA11918 NA11920 NA11931 NA11933 NA11995 NA12004 NA12006 NA12044 NA12046"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/simulations/${pop}/${sex}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

# Havent separated out chr 8 complex and not complex but maybe will do this?
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt

#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt

for i in $samplelist
do
    # all of PARs
    sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf | wc -l`

    # Then subset on whether these calls are in complex, XTR or non complex regions
    # complex
    simcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chr8.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chr8.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SimpleRepeats.bed -u | wc -l`

    # Not complex
    simnotcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chr8.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chr8.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chr8.human.GRCh38.Repeats.RepeatMasker.bed -v | wc -l`

    # add to results file
    #echo -e "${i}\t${sim}"
    #echo -e "${i}\t${sim}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt

    echo -e "${i}\t${simcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
    echo -e "${i}\t${simnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt

done


#---------------#
# ALL AUTOSOMES #
#---------------#
#for t in {1..22}
for t in {8..8}
do
    echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_chr${t}_${ploidy}_simulated_counts.txt
    for i in $samplelist
    do
          sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_chr${t}_golden.vcf | wc -l`
          echo -e "${i}\t${sim}" >> ${pop}_${sex}_chr${t}_${ploidy}_simulated_counts.txt
    done

  echo "chr${t}"
done


#------------------------------------------------------------------------------#
# FEMALES - Chr M #
#------------------------------------------------------------------------------#
# Define variables
pop="EUR"
sex="females"
chr="chrM"
ploidy="haploid"
coverage="20"
samplelist="NA06985 NA06989 NA07000 NA07037 NA07056 NA10847 NA11830 NA11832 NA11840 NA11892 NA11894 NA11918 NA11920 NA11931 NA11933 NA11995 NA12004 NA12006 NA12044 NA12046"
inputdatapath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/simulations/${pop}/${sex}"
regionsbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/complex_region_beds"

# Havent separated out chr 8 complex and not complex but maybe will do this?
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
#echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt
echo -e "Sample\tNum_sim_variants" > ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt

for i in $samplelist
do
    # all of PARs
    sim=`grep -v "#" ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf | wc -l`

    # Then subset on whether these calls are in complex, XTR or non complex regions
    # complex
    #simcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed -u | wc -l`

    # Not complex
    #simnotcomplex=`bedtools intersect -a ${inputdatapath}/${i}/${i}_NEAT_simulated_${coverage}x_${ploidy}_len150_PE_${chr}_golden.vcf -b ${regionsbed}/chrX.human.GRCh38.Repeats.SegmentalDups.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SelfChain.bed ${regionsbed}/chrX.human.GRCh38.Repeats.SimpleRepeats.bed ${regionsbed}/chrX.human.GRCh38.Repeats.RepeatMasker.bed -v | wc -l`

    # add to results file
    echo -e "${i}\t${sim}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts.txt
    #echo -e "${i}\t${simcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_complex.txt
    #echo -e "${i}\t${simnotcomplex}" >> ${pop}_${sex}_${chr}_${ploidy}_simulated_counts_not_complex.txt

done
