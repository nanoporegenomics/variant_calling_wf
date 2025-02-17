import argparse
import pandas as pd
import pysam
import os
import sys
import time
from multiprocessing import Pool

"""
	convert vcf to samples as rows and metadata and variant positions as columns
	input vcf must be bgzipped and indexed with tabix 
	ex: tabix multisample.vcf.gz
	
	Usage: 
	python vcf_to_meta_csv.parallel.py -i multisample.vcf.gz -c HBCC_covs.csv -o combined.meta.tsv

	output: 
	A tsv ( or csv if desired ) for each chromosome with metadata and variant information for each sample. 

	These can be used individually or combind with a series of cut,paste commands.

	Melissa Meredith
	02/2025
"""

def vcf_to_csv(chrom, vcf_file, covariant_file, csv_file):

	start_time = time.time()  
	print(f"Process {chrom} at: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}")

	# read in metadata and reformat ids to match vcf
	covs = pd.read_csv(covariant_file)
	covs['SampleId'] = covs['SampleId']+"_FTX"
	covs.set_index("SampleId", inplace=True)

	cov_samples = list(covs.index)
	# print(cov_samples)

	

	# open the vcf with pysam
	vcf = pysam.VariantFile(vcf_file)

	# set up a dictionary to store the samples as rows
	# sample_names = list(vcf.header.samples)
	# only use samples in the covariate file instead of all the ones in the vcf
	sample_names = cov_samples
	sample_records = {}
	for s in sample_names:
		for hap in ['H0','H1']:
			sample_records[s+"_"+hap] = {}
			sample_records[s+"_"+hap] = {"hap":hap[-1], 
										"age":covs.loc[s]['AgeDeath'], 
										'sex':covs.loc[s]['Sex at Birth'],
										'PMI':covs.loc[s]['PMI'],
										'ancestry':covs.loc[s]['Race'],
										'brainRegion':covs.loc[s]['Region']}
	# reformat to make samples's rows
	sdf = pd.DataFrame(sample_records)
	sdf = sdf.transpose()
	# print(sdf)

	# make a variant dictionary 
	variants = {}

	for record in vcf.fetch(chrom):
		# just get the first alt record... I guess we need to check if there are multi-allelic positions. there souldn't be
		variant_id = f"{record.chrom}_{record.pos}_{record.ref}_{record.alts[0]}" 
		variants[variant_id] = {}

        # for each sample in the vcf
		for sample in record.samples:

			# if this sample is in the cov file store it's genotypes
			if sample in cov_samples:

				# get the genotypes 
				genotype = record.samples[sample]["GT"]
				# print(sample, genotype[0], genotype[1], variant_id)
				
				# store each genotype under it's variant id and sample haplotype id
				if genotype[0] is not None:
					variants[variant_id][sample+"_H0"]=str(genotype[0])
					variants[variant_id][sample+"_H1"]=str(genotype[1])
				else:
					# could change to None if that is better
					variants[variant_id][sample+"_H0"]="."
					variants[variant_id][sample+"_H1"]="."


	# convert the variants to  a DataFrame
	vdf = pd.DataFrame(variants)


	# combine the meta data and the variants 
	meta_variant_sdf = pd.concat([sdf, vdf], axis=1)

	# # Save to tsv/csv file 
	csv_file_chr = f"{chrom}_{csv_file}"
	meta_variant_sdf.to_csv(csv_file_chr, index=True, sep="\t")

	print(f"{chrom} TSV file saved as {csv_file_chr}")

	vcf.close()

	end_time = time.time()  # Record end time
	print(f"{chrom} completed at: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}")
	print(f"Total execution time for {chrom}: {end_time - start_time:.2f} seconds\n")


# process each chromosome in parallel
def process_chromosomes_in_parallel(vcf_file, cov_csv, output_csv, chromosomes):
    # create a pool of worker processes
    with Pool() as pool:
        pool.starmap(vcf_to_csv, [(chrom, vcf_file, cov_csv, output_csv) for chrom in chromosomes])



if __name__ == "__main__":

	# Make an argument parser
	parser = argparse.ArgumentParser(description="Process a vcf file using pysam and replace genotypes.")

	# add argument for the input vcf file
	parser.add_argument(
		"-i","--in_vcf",
		required=True,
		help="Path to the input multisample vcf file to be analyzed. Can be bgzipped, having an index would increase processing speed."
	)

	# add argument for the input vcf file
	parser.add_argument(
		"-c","--cov_csv",
		required=True,
		help="Path to the cohort covariant data."
	)

	# add argument for the output file
	parser.add_argument(
		"-o","--out_csv",
		type=str,
		required=True,
		help="Path to the output csv file."
	)

	# add optional argument for the chromosomes to run on
	parser.add_argument(
		'--chromosomes', 
		nargs='*', 
		default=None, 
		help="List of chromosomes to process (optional). If not provided, all chromosomes will be used.")


	if len(sys.argv) == 0:
		parser.print_help(sys.stderr)
		sys.exit(1)

	# Parse arguments
	args = parser.parse_args()

	# process whole file at once
	# vcf_to_csv(args.in_vcf, args.cov_csv, args.out_csv)

	# if no chromosomes are provided, use all chromosomes from the VCF header
	if args.chromosomes is None:

		# make a list of chromosomes in the vcf header 
		vcfh = pysam.VariantFile(args.in_vcf)
		args.chromosomes = list(vcfh.header.contigs.keys())
		# print(chromosomes)
		vcfh.close()

	# run each chromosome in a different thread 
	process_chromosomes_in_parallel(args.in_vcf, args.cov_csv, args.out_csv, args.chromosomes)




