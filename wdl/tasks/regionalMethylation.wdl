version 1.0

workflow regionalMethylation {
	meta {
	    author: "Shloka Negi"
	    email: "shnegi@ucsc.edu"
	    description: "Calculates haplotype-specific average percent modification per region, using bedMethyl files generated from modKit"
	}	
	
	parameter_meta {
	    BEDS: "List of 1 diploid BED or 2 haplotype-specific BED files with modified base information (Should be gzipped)"
	    SAMPLE: "Sample Name. Will be used in output file"
	    REGIONS_BED: "BED file with regions of interest"
	    REF_NAME: "Reference name. Will be used in output bed file."
	    REGION_TYPE: "Descriptor of the regions of interest in the bed. eg: promoters, CpG_Islands, CCREs, etc.."
	}

    input {
		Array[File] BEDS
		String SAMPLE
		String REF_NAME
		String REGION_TYPE
		File REGIONS_BED
	}

    if (length(BEDS) == 1){
        call calcMeth {
            input:
            bedMethyls = BEDS,
            sample    = SAMPLE,
            ref_name  = REF_NAME,
            region_type = REGION_TYPE,
            regions_bed = REGIONS_BED
        }

    }

    if (length(BEDS) == 2){
        call calcMeth_phased {
            input:
            bedMethyls = BEDS,
            sample     = SAMPLE,
            ref_name   = REF_NAME,
            region_type = REGION_TYPE,
            regions_bed = REGIONS_BED
		}
	}

	output {
		File? output_bed = calcMeth.comb_output_bed
		File? output_qc_tsv = calcMeth.qc_tsv
		File? output_pdf = calcMeth.pdf
		File? hp1_output_bed = calcMeth_phased.hp1_output_bed
		File? hp2_output_bed = calcMeth_phased.hp2_output_bed
		File? hp1_output_qc_tsv = calcMeth_phased.hp1_qc_tsv
		File? hp2_output_qc_tsv = calcMeth_phased.hp2_qc_tsv
		File? hp1_output_pdf = calcMeth_phased.hp1_pdf
		File? hp2_output_pdf = calcMeth_phased.hp2_pdf
	}
}

task calcMeth {
    parameter_meta {
        cov: "(OPTIONAL) Valid coverage threshold for filtering CpGs. Default: 5"
        min_cpg: "(OPTIONAL) Minimum number of CpGs required for regional methylation calculation. Default: 10"
    }

    input {
		Array[File] bedMethyls
		String sample
		String ref_name
		String region_type
		File regions_bed
		Int cov = 5
		Int min_cpg = 10
		Int memSizeGB = 40
		Int diskSizeGB = 5*round(size(bedMethyls, "GB")) + 40
    }

    File bedMethyl = select_first(bedMethyls)
    String suffix = "combined."+"~{ref_name}"+"."+"~{region_type}"

    command <<<

        set -eux -o pipefail

        zcat ~{bedMethyl} | sort -k1,1 -k2,2n -k3,3n > samp.comb.UF.bed
        # Remove mods low coverage and sort file to match sorting criteria with regions_bed
        awk -v c=~{cov} '$5>=c' samp.comb.UF.bed > samp.comb.bed
        head samp.comb.bed
        wc -l samp.comb.bed

        # Sort regions file
        sort -k1,1 -k2,2n -k3,3n ~{regions_bed} > regions.sort.bed

        # Count number of "total" measured CpGs per target region
        bedtools map -a regions.sort.bed -b samp.comb.UF.bed -c 11 -o count | awk '{print $NF}' > n_total_cpgs
        paste regions.sort.bed n_total_cpgs | awk '{OFS="\t"}{print $0}' > regions.nt.bed
        rm regions.sort.bed

        # Generate other stats only using well-covered CpGs
        bedtools map -a regions.nt.bed -b samp.comb.bed -c 11,11,5,11,5,2 -o count,mean,mean,collapse,collapse,collapse > bedmap.tsv
        rm samp.comb.UF.bed

        # Process files in R
        Rscript /getMethylStats.R -i bedmap.tsv -m ~{min_cpg} -o ~{sample}.~{suffix}
        
    >>>
    
    output {
        File comb_output_bed = "~{sample}.~{suffix}.methyl_stats.bed"
        File qc_tsv = "qc_stats.tsv"
        File pdf = "methyl_plots.pdf"
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/shnegi/regionalmethylation:latest"
        preemptible: 1
    }
}

task calcMeth_phased {

	parameter_meta {
		cov: "(OPTIONAL) Valid coverage threshold for filtering CpGs. Default: 3"
		min_cpg: "(OPTIONAL) Minimum number of CpGs required for regional methylation calculation. Default: 10"
	}

	input {
		Array[File] bedMethyls
		String sample
		String ref_name
		String region_type
		File regions_bed
		Int cov = 3
		Int min_cpg = 10
		Int memSizeGB = 40
		Int diskSizeGB = 5*round(size(bedMethyls, "GB")) + 40
	}

	String suffix = "~{ref_name}"+"."+"~{region_type}"

	command <<<

		set -eux -o pipefail

		# Sort regions file
		sort -k1,1 -k2,2n -k3,3n ~{regions_bed} > regions.sort.bed
		# Sort BED files
		FID=1
		for FF in ~{sep=" " bedMethyls}
		do
			zcat $FF | sort -k1,1 -k2,2n -k3,3n > samp.hp${FID}.UF.bed
			# Remove mods with low coverage and sort file to match sorting criteria with regions_bed
			awk -v c=~{cov} '$5>=c' samp.hp${FID}.UF.bed > samp.hp${FID}.bed
			# Count number of "total" measured CpGs per target region
			bedtools map -a regions.sort.bed -b samp.hp${FID}.UF.bed -c 11 -o count | awk '{print $NF}' > n_total_cpgs
			paste regions.sort.bed n_total_cpgs | awk '{OFS="\t"}{print $0}' > regions.nt.bed
			# Generate other stats only using well-covered CpGs
			bedtools map -a regions.nt.bed -b samp.hp${FID}.bed -c 11,11,5,11,5,2 -o count,mean,mean,collapse,collapse,collapse > bedmap.tsv
			rm samp.hp${FID}.UF.bed
			# Process files in R
			Rscript /getMethylStats.R -i bedmap.tsv -m ~{min_cpg} -o ~{sample}.hp${FID}.~{suffix}
			FID=$((FID+1))
		done

	>>>

	output {
		File hp1_output_bed = "~{sample}.hp1.~{suffix}.methyl_stats.bed"
		File hp2_output_bed = "~{sample}.hp2.~{suffix}.methyl_stats.bed"
		File hp1_qc_tsv = "~{sample}.hp1.~{suffix}.qc_stats.tsv"
		File hp2_qc_tsv = "~{sample}.hp2.~{suffix}.qc_stats.tsv"
		File hp1_pdf = "~{sample}.hp1.~{suffix}.methyl_plots.pdf"
		File hp2_pdf = "~{sample}.hp2.~{suffix}.methyl_plots.pdf"
	}

	runtime {
		memory: memSizeGB + " GB"
		disks: "local-disk " + diskSizeGB + " SSD"
		docker: "quay.io/shnegi/regionalmethylation:latest"
		preemptible: 1
	}
}
