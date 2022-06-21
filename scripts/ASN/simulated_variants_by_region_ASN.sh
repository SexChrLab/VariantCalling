#!/bin/bash
#SBATCH --job-name=simulated_variants_by_region_ASN # Job name
#SBATCH -o slurm.simulated_variants_by_region_ASN.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.simulated_variants_by_region_ASN.err                # STDERR (%j = JobId)
#SBATCH -n 1
#SBATCH -t 7-00:00:00
#SBATCH -q tempboost
#SBATCH --mail-user=amtarave@asu.edu # send-to address
#SBATCH --mail-type=ALL


cd /scratch/amtarave/variant_calling_simulations_project/performance_metrics_by_region/simulated_variants_by_region/20_samples


#------------------------------------------------------------------------------#
# MALES - Chr X PARs #
#------------------------------------------------------------------------------#
# Define variables
pop="ASN"
sex="males"
chr="chrX_PARs"
ploidy="diploid"
coverage="20"
samplelist="NA18530 NA18534 NA18536 NA18543 NA18544 NA18546 NA18548 NA18549 NA18557 NA18558 NA18559 NA18561 NA18562 NA18563 NA18572 NA18603 NA18605 NA18606 NA18608 NA18609"
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
pop="ASN"
sex="males"
chr="chrX_nonPARs"
ploidy="haploid"
coverage="10"
samplelist="NA18530 NA18534 NA18536 NA18543 NA18544 NA18546 NA18548 NA18549 NA18557 NA18558 NA18559 NA18561 NA18562 NA18563 NA18572 NA18603 NA18605 NA18606 NA18608 NA18609"
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
pop="ASN"
sex="males"
chr="chrY"
ploidy="haploid"
coverage="10"
samplelist="NA18530 NA18534 NA18536 NA18543 NA18544 NA18546 NA18548 NA18549 NA18557 NA18558 NA18559 NA18561 NA18562 NA18563 NA18572 NA18603 NA18605 NA18606 NA18608 NA18609"
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
pop="ASN"
sex="males"
chr="chr8"
ploidy="diploid"
coverage="20"
samplelist="NA18530 NA18534 NA18536 NA18543 NA18544 NA18546 NA18548 NA18549 NA18557 NA18558 NA18559 NA18561 NA18562 NA18563 NA18572 NA18603 NA18605 NA18606 NA18608 NA18609"
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
pop="ASN"
sex="males"
chr="chrM"
ploidy="haploid"
coverage="20"
samplelist="NA18530 NA18534 NA18536 NA18543 NA18544 NA18546 NA18548 NA18549 NA18557 NA18558 NA18559 NA18561 NA18562 NA18563 NA18572 NA18603 NA18605 NA18606 NA18608 NA18609"
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
pop="ASN"
sex="females"
chr="chrX_PARs"
ploidy="diploid"
coverage="20"
samplelist="NA18525 NA18526 NA18528 NA18531 NA18532 NA18533 NA18535 NA18537 NA18538 NA18539 NA18541 NA18542 NA18545 NA18547 NA18550 NA18552 NA18553 NA18555 NA18560 NA18564"
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
pop="ASN"
sex="females"
chr="chrX_nonPARs"
ploidy="diploid"
coverage="20"
samplelist="NA18525 NA18526 NA18528 NA18531 NA18532 NA18533 NA18535 NA18537 NA18538 NA18539 NA18541 NA18542 NA18545 NA18547 NA18550 NA18552 NA18553 NA18555 NA18560 NA18564"
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
pop="ASN"
sex="females"
chr="chr8"
ploidy="diploid"
coverage="20"
samplelist="NA18525 NA18526 NA18528 NA18531 NA18532 NA18533 NA18535 NA18537 NA18538 NA18539 NA18541 NA18542 NA18545 NA18547 NA18550 NA18552 NA18553 NA18555 NA18560 NA18564"
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
pop="ASN"
sex="females"
chr="chrM"
ploidy="haploid"
coverage="20"
samplelist="NA18525 NA18526 NA18528 NA18531 NA18532 NA18533 NA18535 NA18537 NA18538 NA18539 NA18541 NA18542 NA18545 NA18547 NA18550 NA18552 NA18553 NA18555 NA18560 NA18564"
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
