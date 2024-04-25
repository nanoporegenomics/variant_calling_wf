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
        sample_name: "Sample name. Will be used in output bed file."
        REF_NAME: "Reference name. Will be used in output bed file."
        REGION_TYPE: "Descriptor of the regions of interest in the bed. eg: promoters, CpG_Islands"
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
            bedMethyl = select_first(BEDS),
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
        File? output_bed = calcMeth.output_bed
		File? hp1_output_bed = calcMeth_phased.hp1_output_bed
        File? hp2_output_bed = calcMeth_phased.hp2_output_bed
    }
}

task calcMeth {
    input {
        File bedMethyl
        String sample
        String ref_name
        String region_type
        Int column = 11
        String function = "mean"
        File regions_bed
        Int memSizeGB = 100
        Int threadCount = 10
        Int diskSizeGB = 5*round(size(bedMethyl, "GB")) + 40
    }

    String out_sample_name_ref = "calcMeth_"+"~{sample}"+"_"+"~{ref_name}"+"_"+"~{region_type}"+".bed"

    command <<<

        set -eux -o pipefail

        # Sort BED file
        # Remove mods with 0 coverage/score and sort file to match sorting criteria with regions_bed
        zcat ~{bedMethyl} | \
         awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' | \
         awk '$5!=0' | sort -k1,1 -k2,2n -k3,3n > sorted.bed
        # Run bedtools map
        bedtools map -a ~{regions_bed} -b sorted.bed -c ~{column} -o ~{function} > ~{out_sample_name_ref}
        

    >>>
    
    output {
        File output_bed = "~{out_sample_name_ref}"
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
    input {
        Array[File] bedMethyls
        String sample
        String ref_name
        String region_type
        File regions_bed
        Int memSizeGB = 100
        Int threadCount = 10
        Int diskSizeGB = 5*round(size(bedMethyls, "GB")) + 40
    }

    String out_sample_name_ref = "_"+"~{ref_name}"+"_"+"~{region_type}"+".bed"

    command <<<

        set -eux -o pipefail

        # Sort BED files
        FID=1
        for FF in ~{sep=" " bedMethyls}
        do
            # Remove mods with 0 coverage/score and sort file to match sorting criteria with regions_bed
            zcat $FF | awk '$5!=0' | sort -k1,1 -k2,2n -k3,3n > samp.hp${FID}.bed
            # Run bedtools map
            bedtools map -a ~{regions_bed} -b samp.hp${FID}.bed -c 11 -o mean > ~{sample}.hp${FID}~{out_sample_name_ref}
            FID=$((FID+1))
        done
        

    >>>
    
    output {
        File hp1_output_bed = "~{sample}.hp1~{out_sample_name_ref}"
        File hp2_output_bed = "~{sample}.hp2~{out_sample_name_ref}"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/bedtools@sha256:1b3f971e340dcd148d46e8e6c128f969850d1a1a45aa40f83060f6937e5f1a60"
        preemptible: 1
    }
}