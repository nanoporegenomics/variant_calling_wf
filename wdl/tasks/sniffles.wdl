version 1.0

task sniffles_t {
  input {
    Int threads = 22
    File bamAlignment
    File bamAlignmentIndex
    File? vntrAnnotations
    String sample = "sniffles"
    Boolean phaseVariants = true
    String extraArgs = ""
    Int minSvLen = 25
    Int memSizeGb = 32
    Int diskSizeGb = 256
    File? resourceLogScript
  }

  String trfString = if defined(vntrAnnotations) then "--tandem-repeats " else ""
  String phaseArg = if phaseVariants then "--phase " else ""

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

    ln -s ~{bamAlignment} reads.bam
    ln -s ~{bamAlignmentIndex} reads.bam.bai
    
    sniffles -i reads.bam -v ~{sample}.sniffles.vcf --snf ~{sample}.snf -t ~{threads} ~{trfString}~{vntrAnnotations} \
      ~{phaseArg} --minsvlen ~{minSvLen} ~{extraArgs} 2>&1 | tee ~{sample}_sniffles.log

  >>>

  output {
    File snifflesVcf = "~{sample}.sniffles.vcf"
    File? snifflesLog = "~{sample}_sniffles.log"
    File snifflesSnf = "~{sample}.snf"
    File? toplog = "top.log"
  }

  runtime {
    preemptible: 2
    docker: "meredith705/card_sniffles:2.2"
    cpu: threads
    memory: memSizeGb + " GB"
    disks: "local-disk " + diskSizeGb + " SSD"
  }
}
