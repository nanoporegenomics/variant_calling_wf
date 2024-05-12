version 1.0

workflow run_glnexus{
    input {
        Array[File] vcfFiles = []
        File? regional_bed
        String out_name
        String configuration = "DeepVariantWGS"
        String? extraArgs = ""
        String dockerImage = "quay.io/mlin/glnexus:v1.2.7"

    }

    # Scatter the filterVCF task over the input VCF files to remove NoCall lines
    scatter (input_vcf in vcfFiles) {
        call filterVCF as filter_vcf {
            input:
                input_vcf = input_vcf
        }
    }

    call glnexus as glnexux_merge{
        input:
					vcfFiles = filter_vcf.filtered_vcf,
					out_name = out_name,
					regional_bed = regional_bed,
                    dockerImage = dockerImage
    }

    output{
        File? merged_gvcf = glnexux_merge.merged_gvcf
    }
}

task glnexus {
    input {
        Array[File] vcfFiles = []
        File? regional_bed
        String out_name
        String configuration = "DeepVariantWGS"
        String? extraArgs = ""
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 5 * round(size(vcfFiles, 'G')) + 300
        String dockerImage = "quay.io/mlin/glnexus:v1.2.7"

    }

    Boolean regionaly =  defined(regional_bed)

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        # note output format: https://nanoporetech.github.io/modkit/intro_bedmethyl.html
        # filtering threshold default 10-th percentile of calls https://github.com/nanoporetech/modkit/blob/master/filtering.md
        if [ ~{regionaly} == true ]
        then
          glnexus_cli \
          --config ~{configuration} \
          --bed ~{regional_bed} \
          ~{sep=" " select_all(vcfFiles)} | bcftools view \
         | bgzip -@ ~{threadCount} > ~{out_name}.deepvariant.cohort.vcf.gz

        else
          glnexus_cli \
          --config ~{configuration} \
          ~{sep=" " vcfFiles} | bcftools view \
         | bgzip -@ ~{threadCount} > ~{out_name}.deepvariant.cohort.vcf.gz
        fi


    >>>

        output {
            File? merged_gvcf = "~{out_name}.deepvariant.cohort.vcf.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}

task filterVCF {
    input {
        File input_vcf
        Int memSizeGB = 128
        Int threadCount = 4
        Int diskSizeGB = 2 * round(size(input_vcf, 'G')) + 30
        String dockerImage = "quay.io/mlin/glnexus:v1.2.7"
    }
    
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace
        # Define a temp file for filtering
        tmp_file=~{input_vcf}.tmp
        
        # Remove 'NoCall' lines and write to the temporary file
        bgzip -dc ~{input_vcf} | grep -v NoCall | bgzip > $tmp_file
        
        # Rename the temporary file to the input VCF name to save space
        mv $tmp_file ~{input_vcf}
    >>>

    output {
        File filtered_vcf = "~{input_vcf}"
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}
