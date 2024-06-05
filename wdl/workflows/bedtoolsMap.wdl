version 1.0

workflow bedmethylCalcMeth {
    meta {
        author: "Shloka Negi"
        email: "shnegi@ucsc.edu"
        description: "Calculates combined average regional methylation, using bedMethyl generated from modkit"
    }

    parameter_meta {
        BED: "combined BED file with modified base information (Should be gzipped)"
        SAMPLE: "Sample Name. Will be used in output file"
        REGIONS_BED: "BED file with regions of interest"
    }

    input {
        File BED
        String SAMPLE
        File REGIONS_BED
    }

    call calcMeth {
        input:
        bedMethyl=BED,
        sample=SAMPLE,
        regions_bed=REGIONS_BED
    }
    
    output {
        File comb_output_bed = calcMeth.comb_output_bed
    }
}

task calcMeth {
    parameter_meta {
        cov: "(OPTIONAL) Valid coverage threshold for filtering CpGs. Default: 5"
        min_cpg: "(OPTIONAL) Minimum number of CpGs required for regional methylation calculation. Default: 10"
    }

    input {
        File bedMethyl
        String sample
        File regions_bed
        Int cov = 5
        Int min_cpg = 10
        Int memSizeGB = 50
        Int threadCount = 8
        Int diskSizeGB = 5*round(size(bedMethyl, "GB")) + 40
    }

    command <<<

        set -eux -o pipefail

        zcat ~{bedMethyl} | sort -k1,1 -k2,2n -k3,3n > samp.comb.UF.bed
        # Remove mods with coverage<4 and sort file to match sorting criteria with regions_bed
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

        # Calculate mean methylation per region
        bedtools map -a regions.ratios.bed -b samp.comb.bed -c 11 -o mean > ~{sample}.comb.regionMethyl.bed

    >>>
    
    output {
        File comb_output_bed = "~{sample}.comb.regionMethyl.bed"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/bedtools@sha256:1b3f971e340dcd148d46e8e6c128f969850d1a1a45aa40f83060f6937e5f1a60"
        preemptible: 1
    }
}
