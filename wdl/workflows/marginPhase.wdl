version 1.0

workflow runMarginPhase {
    input {
        File smallVariantsFile
        File structuralVariantsFile
        File refFile
        File bamFile
        String sampleName
        String dockerImage = "mkolmogo/card_harmonize_vcf:0.1"
        File? resourceLogScript
    }

    call combineVcfs {
        input:
            smallVariantsFile = smallVariantsFile,
            structuralVariantsFile = structuralVariantsFile,
            sampleName = sampleName,
            dockerImage = dockerImage
    }

    call marginPhase {
        input:
        combinedVcfFile = combineVcfs.outVcf,
        refFile = refFile,
        bamFile = bamFile,
        sampleName = sampleName,
        dockerImage = dockerImage,
        resourceLogScript = resourceLogScript
    }

    output {
        File out_margin_phase_svs = marginPhase.phasedVcf
        File out_margin_phase_bam = marginPhase.haplotaggedBam
        File out_margin_phase_bam_bai = marginPhase.haplotaggedBamIdx
    }
}

task combineVcfs {
    input {
        File smallVariantsFile
        File structuralVariantsFile
        String sampleName
        String dockerImage
        Int svLengthCutoff = 25
        Int threads = 32
        Int memSizeGb = 128
        Int diskSizeGb = 256
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # check if input is bgzipped or not
        SV_FILTERED=~{structuralVariantsFile}_size_filtered.vcf
        SMALL_FILTERED=~{smallVariantsFile}_size_filtered.vcf

        echo ~{sampleName} > samplename.txt
        
        #-f option supports unzgipped input
        zcat -f ~{structuralVariantsFile} | python3 /opt/vcf_filter_size.py greater ~{svLengthCutoff} | bcftools reheader -s samplename.txt | bgzip > $SV_FILTERED
        tabix -p vcf $SV_FILTERED
        zcat -f ~{smallVariantsFile} | python3 /opt/vcf_filter_size.py less ~{svLengthCutoff} | bcftools reheader -s samplename.txt | bgzip > $SMALL_FILTERED
        tabix -p vcf $SMALL_FILTERED

        bcftools concat -a $SMALL_FILTERED $SV_FILTERED -o ~{sampleName}.merged_small_svs.vcf
    >>>
    output {
        File outVcf = "~{sampleName}.merged_small_svs.vcf"
    }
    runtime {
        preemptible: 2
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        docker: dockerImage
    }
}

task marginPhase {
    input {
        File combinedVcfFile
        File refFile
        File bamFile
        String sampleName
        String dockerImage
        String marginOtherArgs = ""
        Int threads = 32
        Int memSizeGb = 64
        Int diskSizeGb = 256
        File? resourceLogScript
    }
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
        
        samtools index -@ ~{threads} ~{bamFile}
        samtools faidx ~{refFile}
        mkdir output/
        margin phase ~{bamFile} ~{refFile} ~{combinedVcfFile} /opt/margin/params/phase/allParams.phase_vcf.ont.sv.json -t ~{threads} ~{marginOtherArgs} -o output/~{sampleName}_hvcf 
        bgzip output/~{sampleName}_hvcf.phased.vcf
        samtools index -@ ~{threads} output/~{sampleName}_hvcf.haplotagged.bam

    >>>
    output {
        File phasedVcf = "output/~{sampleName}_hvcf.phased.vcf.gz"
        File haplotaggedBam = "output/~{sampleName}_hvcf.haplotagged.bam"
        File haplotaggedBamIdx = "output/~{sampleName}_hvcf.haplotagged.bam.bai"
        File? toplog = "top.log"
    }

    runtime {
        preemptible: 2
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        docker: dockerImage
    }
}
