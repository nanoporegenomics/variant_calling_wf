version 1.0

task modkit {
    input {
        File haplotaggedBam
        File haplotaggedBamBai
        File ref
        File? regional_bed
        String sample_name
        String ref_name = "ref"
        String out_type_filter = "--cpg"
#        String partitionTag = "--partition-tag HP"
        Boolean partitionTag
        String? extraArgs = ""
        Int memSizeGB = 64
        Int threadCount = 64
        Int diskSizeGB = 4 * round(size(haplotaggedBam, 'G')) + round(size(ref, 'G')) + 100
        String dockerImage = "meredith705/modkit:latest"
        File? resourceLogScript
    }

    parameter_meta {
        haplotaggedBam: "Guppy with Remora reads aligned to assembly. in BAM format."
        ref: "Assembly (reference) to that reads are aligned to."
        sample_name: "Sample name. Will be used in output bed file."
        ref_name: "Reference name. Will be used in output bed file."
        #modType: "Modified base of interest, one of: 5mC (default), 5hmC, 5fC, 5caC, 5hmU, 5fU, 5caU, 6mA, 5oxoG, Xao."
        out_type_filter: " Output records filtered to these sites: cpg (default), chg, chh"
    }

    String sample_name_ref = "~{sample_name}"+"_"+"~{ref_name}"
    Boolean regionaly =  defined(regional_bed)
    String regional_arg = "--include-bed "+"~{regional_bed}"
    String unphased_outFile = "modkit_out/" +"~{sample_name_ref}"+".bed"

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        ## run a recurrent "top" in the background to monitor resource usage
        if [ ~{resourceLogScript} != "" ]
        then
            bash ~{resourceLogScript} 20 top.log &
        fi

        ln -s ~{haplotaggedBam} reads.bam
        ln -s ~{haplotaggedBamBai} reads.bam.bai


        mkdir modkit_out

        # note output format: https://nanoporetech.github.io/modkit/intro_bedmethyl.html
        # filtering threshold default 10-th percentile of calls https://github.com/nanoporetech/modkit/blob/master/filtering.md
        if [ ~{partitionTag} == true ]
        then
            modkit pileup ~{out_type_filter} --partition-tag HP --prefix ~{sample_name_ref} --ref ~{ref} \
                   --threads ~{threadCount} --only-tabs ~{extraArgs} reads.bam modkit_out
        fi

        if [ ~{partitionTag} == false ]
        then
            # modkit command with reference on input reads
            modkit pileup ~{out_type_filter} --prefix ~{sample_name_ref} --ref ~{ref} \
                    --threads ~{threadCount} --only-tabs ~{extraArgs} reads.bam ~{unphased_outFile}
        fi

        if [ ~{partitionTag} == true ]
        then
            bgzip -@ ~{threadCount} modkit_out/~{sample_name_ref}_1.bed
            bgzip -@ ~{threadCount} modkit_out/~{sample_name_ref}_2.bed
            bgzip -@ ~{threadCount} modkit_out/~{sample_name_ref}_ungrouped.bed
        else
            bgzip -@ ~{threadCount} modkit_out/~{sample_name_ref}.bed
        fi

    >>>

        output {
            File? hap1bedOut      = "modkit_out/~{sample_name_ref}_1.bed.gz"
            File? hap2bedOut      = "modkit_out/~{sample_name_ref}_2.bed.gz"
            File? ungroupedBedOut      = "modkit_out/~{sample_name_ref}_ungrouped.bed.gz"
            File? wholeGenomeOut  = "modkit_out/~{sample_name_ref}.bed.gz"
            File? toplog = "top.log"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 2
    }
}

task regionalMethylation {
    input {
        File haplotaggedBam
        File haplotaggedBamBai
        File ref
        File regional_bed
        String region_type
        String sample_name
        String ref_name = "ref"
        Int min_ml_prob = 200
        Int make_heatmap = 0
        #String hap_sparated_stats = "--hap"
        String? extraArgs = ""
        Int memSizeGB = 64
        Int threadCount = 16
        Int diskSizeGB = 4 * round(size(haplotaggedBam, 'G')) + round(size(ref, 'G')) + 100
        String dockerImage = "meredith705/modkit:latest"
        File? resourceLogScript
    }

    parameter_meta {
        haplotaggedBam: "Guppy with Remora reads aligned to assembly. in BAM format."
        ref: "Assembly (reference) to that reads are aligned to."
        sample_name: "Sample name. Will be used in output bed file."
        ref_name: "Reference name. Will be used in output bed file."
        regional_bed: "Bedfile with regions of interst"
        region_type: "Descriptor of the regions of interest in the bed. eg: promoters, CpG_Islands"

    }

    String out_sample_name_ref = "~{sample_name}"+"_"+"~{ref_name}"+"_"+"~{region_type}"+".bed"


    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## run a recurrent "top" in the background to monitor resource usage
        if [ ~{resourceLogScript} != "" ]
        then
            bash ~{resourceLogScript} 20 top.log &
        fi

        python3 /home/apps/methyl_read_grouping/methlotype/regionalMethylation.py -b ~{haplotaggedBam} -s ~{sample_name} -d ~{regional_bed} -e ~{ref_name}_~{region_type} -p ~{min_ml_prob} --heatmap ~{make_heatmap}

    >>>
      output {
        File regionally_aggregated_bed = "~{sample_name}"+"_"+"~{ref_name}"+"_"+"~{region_type}"+".bed"
        File? toplog = "top.log"
      }

      runtime {
        preemptible: 2
        docker: "meredith705/modkit:latest"
        cpu: threadCount
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
      }
}




task plot_ML_hist {
    input {
        File modbam
        File modbam_index
        String sample
        Int buckets = 128
        Int memSizeGB = 50
        Int threadCount = 64
        Int diskSizeGB = 5*round(size(modbam, "GB")) + 20
    }

    command <<<
        set -eux -o pipefail
        set -o xtrace

        # link the modbam to make sure it's index can be found
        ln -s ~{modbam} input.bam
        ln -s ~{modbam_index} input.bam.bai

        # Generate ML distribution (Default: 128 buckets)
        mkdir output_dir/
        modkit sample-probs input.bam -t ~{threadCount} \
            --only-mapped --hist ~{"--buckets " + buckets} \
            --force --prefix ~{sample} -o output_dir/

        # Plot
#        Rscript /opt/scripts/plot-dist.R output_dir/~{sample}_probabilities.tsv output_dir/~{sample}.ML_hist.png
    >>>

    output {
#        File out = "output_dir/~{sample}.ML_hist.png"
        File out_ML_tsv = "output_dir/~{sample}_probabilities.tsv"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "meredith705/modkit:latest"
        preemptible: 1
    }
}

