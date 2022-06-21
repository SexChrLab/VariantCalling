#!/usr/bin/env python3
import csv
from optparse import  OptionParser

"""
# Background:
# For this script, I want to take as input 1) golden simulated VCF for a sample.
# this VCF must be produced by the software NEAT, and 2) called VCF. This is a
# VCF for an individual that was called using simulated reads from NEAT.

# With this input, I want to output:
# 1) Table of performance metrics for a given individual with the following cols:
# Sample
# TP
# FP
# FN
# Validated Genotypes (number of SNPs where the genotypes exactly match between
# golden and called VCFs)

# The golden VCF entries look like this (and will need to be reformatted)
grep -v "#" NA06984_NEAT_simulated_20x_haploid_len101_PE_chrY_nonPARs_golden.vcf | head
chrY    2793653 .       A       G       .       PASS    WP=1
chrY    2868868 .       A       T       .       PASS    WP=1
chrY    2968390 .       A       C       .       PASS    WP=1
chrY    3028467 .       C       T       .       PASS    WP=1
chrY    3491379 .       C       G       .       PASS    WP=1
chrY    3491409 .       A       T       .       PASS    WP=1
chrY    3826492 .       T       G       .       PASS    WP=1
chrY    3846274 .       T       C       .       PASS    WP=1
chrY    4001704 .       G       A       .       PASS    WP=1
chrY    4472585 .       A       G       .       PASS    WP=1

# NOTE TO SELF: For X and Y non PARs the golden VCF is haploid but the called
# VCF is diploid so will need to account for this (make WP diploid so 1/1). Can
# add an if statement in this

# The called VCF
grep -v "#" NA06984_chrY_filtered.vcf | head
chrY    2791306 .       A       C       1604.92 PASS    .       GT:AD:DP:GQ:PL  0/0:11,0:11:33:0,33,441
chrY    2793653 .       A       G       5683.16 PASS    .       GT:AD:DP:GQ:PL  1/1:0,12:12:36:449,36,0
chrY    2795671 .       T       A       1118.8  PASS    .       GT:AD:DP:GQ:PL  0/0:15,0:15:33:0,33,495
chrY    2800415 .       T       C       1860.97 PASS    .       GT:AD:DP:GQ:PL  0/0:13,0:13:39:0,39,501
chrY    2806199 .       G       A       1384.96 PASS    .       GT:AD:DP:GQ:PL  0/0:15,0:15:45:0,45,611
chrY    2815201 .       C       T       1568.11 PASS    .       GT:AD:DP:GQ:PL  0/0:12,0:12:36:0,36,463
chrY    2818514 .       A       G       1590.95 PASS    .       GT:AD:DP:GQ:PL  0/0:18,0:18:42:0,42,630
chrY    2818686 .       A       G       1280.07 PASS    .       GT:AD:DP:GQ:PL  0/0:18,0:18:48:0,48,720
chrY    2820401 .       T       A       1362.91 PASS    .       GT:AD:DP:GQ:PL  0/0:13,0:13:39:0,39,498
chrY    2839031 .       C       T       1574.93 PASS    .       GT:AD:DP:GQ:PL  0/0:10,0:10:30:0,30,400

# NOTE TO SELF: For called VCFs there are 0/0 calls because of joint calling so
# need to ignore 0/0 calls

# TO THINK ABOUT: Note sure how to call true negatives (TN)


else: # called variant is not in golden VCF # need a TRUE FALSE variable
    # so it is a false positive
    print("variant called but not simulated: ", calledVCFLine)
    fpCount += 1

"""

################################################################################
USAGE = """
python3 compare_VCFs_fix.py	--golden <path and name to golden vcf >
                                --called <path and name to called vcf >
                                --gploidy <ploidy of golden VCF>
								--out <path and name for output file>

golden == path and name to golden vcf
called == path and name to called vcf
gploidy == ploidy of golden VCF
out == path and name for output file
"""

parser = OptionParser(USAGE)
parser.add_option('--golden',dest='golden', help = 'path and name to golden vcf')
parser.add_option('--called',dest='called', help = 'path and name to called vcf')
parser.add_option('--gploidy',dest='gploidy', help = 'ploidy of golden VCF')
parser.add_option('--out',dest='out', help = 'path and name to output file')

(options, args) = parser.parse_args()

parser = OptionParser(USAGE)
if options.golden is None:
	parser.error('path and name to golden vcf not given')
if options.called is None:
	parser.error('path and name to called vcf not give')
if options.gploidy is None:
	parser.error('ploidy of golden VCF not give. Please provide number (only excepts the following: haploid, diploid)')
if options.out is None:
	parser.error('path and name for output file not give')
################################################################################
# Set up output file
# 1. performance metrics
outfn1 = (options.out + '_performance_metrics.txt')
outfile1 = open(outfn1, 'w')
header = 'Sample\tTP\tFP\tFN\tValidated_Genotypes\n'
outfile1.write(header)

# 2. False positives from call VCF
outfn2 = (options.out + '_FPs_call_VCF.txt')
outfile2 = open(outfn2, 'w')

# 3. TPs that have called missmatches
outfn3 = (options.out + '_called_genotype_missmatches.txt')
outfile3 = open(outfn3, 'w')

# 4. Alleles that are missmatches (REF and ALT dont match between sim and called)
outfn4 = (options.out + '_called_REF_ALT_missmatches.txt')
outfile4 = open(outfn4, 'w')

# 5. TP positions
outfn5 = (options.out + '_TP_pos.txt')
outfile5 = open(outfn5, 'w')

# 6. FN positions
outfn6 = (options.out + '_FN_pos.txt')
outfile6 = open(outfn6, 'w')

# 7. FP positions (is in 2.)
outfn7 = (options.out + '_FP_pos.txt')
outfile7 = open(outfn7, 'w')


calledVCFentries = []
goldenVCFentries = []

print("Begun reading called VCF")
with open(options.called, 'r') as calledVCF:
    for line in calledVCF:
        if '#CHROM' in line:
            smpl = line.strip().split()[9]
            #print(smpl)
        if '#' in line: # skip VCF header
            continue
        else:
            calledVCFLine = line.strip().split() # process called VCF line
            #print(calledVCFLine)
            # get GT and see if REF/REF
            gt = calledVCFLine[9].split(":")[0]
            #print(gt)
            if gt == '0/0' or gt == './.' or gt == '0' or gt == '.': # don't analyze HOM REF and missing calls
                #print("HOM REF call: ", gt)
                #print(calledVCFLine)
                continue
            elif gt == '0|1':
                #print(calledVCFLine)
                calledVCFLine[9] = '0/1'
                #print(calledVCFLine)
            elif gt == '1|0':
                #print(calledVCFLine)
                calledVCFLine[9] = '0/1'
            elif gt == '1|1' or gt == '1':
                calledVCFLine[9] = '1/1'
            else:
                calledVCFLine[9] = gt
            calledVCFentries.append(calledVCFLine)
#print(calledVCFentries)
print("Finished reading called VCF")


print("Begun reading simulated VCF")
with open(options.golden, 'r') as gldVCF:
    for g in gldVCF:
        if '#' in g: # skip VCF header
            continue
        else:
            gldVCFLine = g.strip().split() # process called VCF line
            #print(gldVCFLine)
            goldengt = gldVCFLine[7].split('=')[1]

            if options.gploidy == 'haploid':
                goldengt = gldVCFLine[7].split('=')[1]
                #print(goldengt)
                goldengtReform1 = goldengt + '/' + goldengt
                gldVCFLine[7] = goldengtReform1
            else:
                gldVCFLine[7] = goldengt
            goldenVCFentries.append(gldVCFLine)
#print(goldenVCFentries)
print("Finished reading simulated VCF")


# TP is defined as variant that is in both call file and sim file (Ref and ALT
# must match between sim and called vcf)
tpCount = 0 # true positives
vgCount = 0 # validated genotypes
fpCount = 0 # false positives
fnCount = 0 # false negatives

print("Getting TPs, FPs, and Valid Genotypes")
for i in calledVCFentries:
    #print(i)
    goldenVCFentriesSubset = [ln for ln in goldenVCFentries if ln[0] == i[0] and ln[1] == i[1]]
    #print(goldenVCFentriesSubset)
    if len(goldenVCFentriesSubset) == 0:
        # this means that call var not in sim vcf, this is a FP
        fpCount += 1
        outfile2.write('\t'.join(i))
        #outfile2.write(str(i))
        outfile2.write('\n')

        fp_pos = i[0] + '\t' + i[1] + '\n'
        outfile7.write(fp_pos)

    else:
        # this means that the called var is in sim vcf so it may be a TP
        # TO DO: Need to check to make sure REF and ALT alleles match
        # cols 3 and 4 need to match
        # If they match then its a TP if not then its what, FP?
        calledREF = i[3]
        calledALT = i[4]
        goldenREF = goldenVCFentriesSubset[0][3]
        goldenALT = goldenVCFentriesSubset[0][4]
        if calledREF == goldenREF and calledALT == goldenALT:
            tpCount += 1

            tp_pos = i[0] + '\t' + i[1] + '\n'
            outfile5.write(tp_pos)
            #print(i)
            #print(goldenVCFentriesSubset[0])
            # now test to see if var is correct (validated genotype)
            calledGT = i[9]
            goldenGT = goldenVCFentriesSubset[0][7]
            if calledGT == goldenGT:
                vgCount += 1
            else:
                #print('\t'.join(goldenVCFentriesSubset[0]))
                #print('\t'.join(i))
                outfile3.write('\t'.join(i))
                outfile3.write('\t')
                outfile3.write('\t'.join(goldenVCFentriesSubset[0]))
                #outfile3.write(str(goldenVCFentriesSubset))
                outfile3.write('\n')
        else:
            print("Mismatch in alleles")
            #print(i)
            #print(goldenVCFentriesSubset[0])
            outfile4.write('\t'.join(i))
            outfile4.write('\t')
            outfile4.write('\t'.join(goldenVCFentriesSubset[0]))
            #outfile3.write(str(goldenVCFentriesSubset))
            outfile4.write('\n')

print('TP: ', tpCount, 'FP: ', fpCount, 'Valid Genotypes: ', vgCount)

# FN is a site that is in golden VCF but not called VCF
for i in goldenVCFentries:
    #print("i test", i)
    calledVCFentriesSubset = [ln for ln in calledVCFentries if ln[0] == i[0] and ln[1] == i[1]]
    if len(calledVCFentriesSubset) == 0:
        #print("i not in called vcf", i)
        # this means that sim var is not in called VCF, this is a FN
        fnCount += 1

        fn_pos = i[0] + '\t' + i[1] + '\n'
        outfile6.write(fn_pos)


print('TP: ', tpCount, 'FP: ', fpCount, 'FN: ', fnCount, 'Valid Genotypes: ', vgCount)
lntowrite = str(smpl) + '\t' + str(tpCount) + '\t' + str(fpCount) + '\t' + str(fnCount) + '\t' + str(vgCount) + '\n'
outfile1.write(lntowrite)



# Close output file when done
outfile1.close()
outfile2.close()
outfile3.close()
outfile4.close()
outfile5.close()
outfile6.close()
outfile7.close()
