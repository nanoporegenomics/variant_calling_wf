version 1.0

workflow bedmethylCalcMeth {
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
		File? hp1_output_bed = calcMeth_phased.hp1_output_bed
        File? hp2_output_bed = calcMeth_phased.hp2_output_bed
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
        Int threadCount = 2
        Int diskSizeGB = 5*round(size(bedMethyls, "GB")) + 40
    }

    File bedMethyl = select_first(bedMethyls)
    String suffix = "combined."+"~{ref_name}"+"."+"~{region_type}"+".bed"

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
        bedtools map -a regions.sort.bed -b samp.comb.UF.bed -c 11 -o count > ~{sample}.comb.unfilt.regionCount.bed
        # Count number of "well-covered" CpGs per target region
        bedtools map -a regions.sort.bed -b samp.comb.bed -c 11 -o count > ~{sample}.comb.filt.regionCount.bed

        # For regions with total number of CpGs < ~{min_cpg}, output "-1"
        awk -v m=~{min_cpg} 'FNR==NR { if ($NF != 0) a[FNR] = $NF; next } { if (a[FNR] >= m) print $NF/a[FNR]; else print "-1" }' ~{sample}.comb.unfilt.regionCount.bed ~{sample}.comb.filt.regionCount.bed > ratios.txt

        rm -rf *regionCount.bed samp.comb.UF.bed
        # Add frac_wellcov_CpGs column to regionalMethyl BED
        paste regions.sort.bed ratios.txt | awk '{OFS="\t"}{print $0}' > regions.ratios.bed

        # Run bedtools map
        bedtools map -a regions.ratios.bed -b samp.comb.bed -c 11 -o mean > "~{sample}.~{suffix}"
        

    >>>
    
    output {
        File comb_output_bed = "~{sample}.~{suffix}"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/bedtools@sha256:1b3f971e340dcd148d46e8e6c128f969850d1a1a45aa40f83060f6937e5f1a60"
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
        Int threadCount = 16
        Int diskSizeGB = 5*round(size(bedMethyls, "GB")) + 40
    }

    String suffix = "~{ref_name}"+"."+"~{region_type}"+".bed"

    command <<<

        set -eux -o pipefail

        # Sort BED files
        FID=1
        for FF in ~{sep=" " bedMethyls}
        do
            zcat $FF | sort -k1,1 -k2,2n -k3,3n > samp.hp${FID}.UF.bed
            # Remove mods with coverage<4 and sort file to match sorting criteria with regions_bed
            awk -v c=~{cov} '$5>=c' samp.hp${FID}.UF.bed > samp.hp${FID}.bed

            # Sort regions file
            sort -k1,1 -k2,2n -k3,3n ~{regions_bed} > regions.sort.bed

            # Count number of "total" measured CpGs per target region
            bedtools map -a regions.sort.bed -b samp.hp${FID}.UF.bed -c 11 -o count > ~{sample}.hp${FID}.unfilt.regionCount.bed
            # Count number of "well-covered" CpGs per target region
            bedtools map -a regions.sort.bed -b samp.hp${FID}.bed -c 11 -o count > ~{sample}.hp${FID}.filt.regionCount.bed

            # For regions with total number of CpGs < ~{min_cpg}, output "-1"
            awk -v m=~{min_cpg} 'FNR==NR { if ($NF != 0) a[FNR] = $NF; next } { if (a[FNR] >= m) print $NF/a[FNR]; else print "-1" }' ~{sample}.hp${FID}.unfilt.regionCount.bed ~{sample}.hp${FID}.filt.regionCount.bed > ratios.txt

            rm -rf *regionCount.bed samp.hp${FID}.UF.bed
            # Add frac_wellcov_CpGs column to regionalMethyl BED
            paste regions.sort.bed ratios.txt | awk '{OFS="\t"}{print $0}' > regions.ratios.bed

            # Calculate mean methylation per region
            bedtools map -a regions.ratios.bed -b samp.hp${FID}.bed -c 11 -o mean > "~{sample}.hp${FID}.~{suffix}"
            FID=$((FID+1))
        done        

    >>>
    
    output {
        File hp1_output_bed = "~{sample}.hp1.~{suffix}"
        File hp2_output_bed = "~{sample}.hp2.~{suffix}"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/bedtools@sha256:1b3f971e340dcd148d46e8e6c128f969850d1a1a45aa40f83060f6937e5f1a60"
        preemptible: 1
    }
}