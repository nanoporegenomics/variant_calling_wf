version 1.0

task dv_t {
  input {
    Int threads
    File reference
    File? bamAlignment
    File bamAlignmentIndex
    String sampleName
    String extraArguments = ""
    String dvModel = "ONT_R104"
    Boolean oneChr = false
    Int memSizeGb = 128
    Int diskSizeGb = 1024
    Int preemptible = 0
    File? resourceLogScript
  }  

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    ## run a recurrent "top" in the background to monitor resource usage
    if [[ "~{resourceLogScript}" != "" ]]
    then
        bash ~{resourceLogScript} 20 top.log &
    fi

    ln -s ~{reference} ref.fa
    samtools faidx ref.fa

    ln -s ~{bamAlignment} reads.bam
    ln -s ~{bamAlignmentIndex} reads.bam.bai

    ## if BAM has reads only for one chromosome
    ## figure out which one and add argument
    REGION_ARG=""
    if [ ~{oneChr} == true ]
    then
        CONTIG_ID=`head -1 < <(samtools view ~{bamAlignment}) | cut -f3`
        REGION_ARG="--regions $CONTIG_ID"
    fi

    ## run DeepVariant
    /opt/deepvariant/bin/run_deepvariant \
        --model_type ~{dvModel} \
        --ref ref.fa \
        --reads reads.bam \
        --sample_name ~{sampleName} \
        --output_vcf dv.vcf.gz $REGION_ARG \
        --output_gvcf dv.g.vcf.gz \
        --num_shards ~{threads} ~{extraArguments}
  >>>

  output {
    File dvVcf = "dv.vcf.gz"
    File dvgVcf = "dv.g.vcf.gz"
    File? toplog = "top.log"
  }

  runtime {
    preemptible: preemptible
    docker: "google/deepvariant:cl508467184"
    cpu: threads
    memory: memSizeGb + " GB"
    disks: "local-disk " + diskSizeGb + " SSD"
  }
}

task margin_t {
  input {
    File reference
    File vcfFile
    File gvcfFile
    File bamAlignment
    File bamAlignmentIndex
    String sampleName
    Int threads
    String marginOtherArgs = ""
    Int memSizeGb = 64
    Int diskSizeGb = 1024
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

    ln -s ~{reference} ref.fa
    samtools faidx ref.fa

    ln -s ~{bamAlignment} reads.bam
    ln -s ~{bamAlignmentIndex} reads.bam.bai
    
    mkdir output/
    margin phase reads.bam ref.fa ~{vcfFile} /opt/margin/params/phase/allParams.haplotag.ont-r104q20.json -t ~{threads} ~{marginOtherArgs} -o output/~{sampleName}

    bgzip output/~{sampleName}.phased.vcf

    # Don't output a bam (-M) for gVCF phasing
    margin phase reads.bam ref.fa ~{gvcfFile} /opt/margin/params/phase/allParams.haplotag.ont-r104q20.json -t ~{threads} ~{marginOtherArgs} -o output/~{sampleName}.g -M

    bgzip output/~{sampleName}.g.phased.vcf

    samtools index -@ ~{threads} output/~{sampleName}.haplotagged.bam
  >>>

  output {
      File phasedVcf = "output/~{sampleName}.phased.vcf.gz"
      File phasedgVcf = "output/~{sampleName}.g.phased.vcf.gz"
      File haplotaggedBam = "output/~{sampleName}.haplotagged.bam"
      File haplotaggedBamIdx = "output/~{sampleName}.haplotagged.bam.bai"
      File? toplog = "top.log"
  }

  runtime {
    preemptible: 2
    docker: "mkolmogo/card_harmonize_vcf:0.1"
    cpu: threads
    memory: memSizeGb + " GB"
    disks: "local-disk " + diskSizeGb + " SSD"
  }
}

task mergeVCFs {
  input {
    Array[File] vcfFiles
    Array[File] gvcfFiles
    String outname = "merged"
    Int memSizeGb = 64   # this can probably be reduced again..
    Int diskSizeGb = 5 * round(size(vcfFiles, 'G')) + 5 * round(size(gvcfFiles, 'G')) + 500
  }

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    mkdir bcftools.tmp
    # vcf merging
    bcftools concat -n ~{sep=" " vcfFiles} | bcftools sort -T bcftools.tmp -O z -o ~{outname}.vcf.gz -
    bcftools index -t -o ~{outname}.vcf.gz.tbi ~{outname}.vcf.gz

    # gvcf merging
    mkdir bcftools.tmp
    bcftools concat -n ~{sep=" " gvcfFiles} | bcftools sort -T bcftools.tmp -O z -o ~{outname}.g.vcf.gz -
    bcftools index -t -o ~{outname}.g.vcf.gz.tbi ~{outname}.g.vcf.gz
  >>>

  output {
      File vcf = "~{outname}.vcf.gz"
      File vcfIndex = "~{outname}.vcf.gz.tbi"
      File gvcf = "~{outname}.g.vcf.gz"
      File gvcfIndex = "~{outname}.g.vcf.gz.tbi"
  }

  runtime {
    preemptible: 2
    docker: "quay.io/biocontainers/bcftools@sha256:95c212df20552fc74670d8f16d20099d9e76245eda6a1a6cfff4bd39e57be01b"
    cpu: 1
    memory: memSizeGb + " GB"
    disks: "local-disk " + diskSizeGb + " SSD"
  }
}
