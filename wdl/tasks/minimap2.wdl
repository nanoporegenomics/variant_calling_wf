version 1.0

task minimap2_t {
  input {
    Int threads
    File reference
    File reads
    String sampleId = ""
    String mapMode = "map-ont"
    Boolean useMd = false
    Boolean useEqx = true
    Int memSizeGb = 128
    Int diskSizeGb = 1024
    Int kmerSize = 17
    String minibatchSize = "5G"
    Boolean customAln = false
    String customArgs = ""
    Int sortMemgb = "4"
    Int preemptible = 0
  }

  String mdString = if useMd then "--MD" else ""
  String eqxString = if useEqx then "--eqx" else ""

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    MM_INPUT=~{reads}
    if [ "${MM_INPUT: -3}" == "bam" ]
    then
      samtools fastq -TMm,Ml,MM,ML ~{reads} | \
        minimap2 -ax ~{mapMode} ~{reference} - -k ~{kmerSize} -y -K ~{minibatchSize} -t ~{threads} ~{mdString} ~{eqxString} | samtools sort -@4 -m ~{sortMemgb}G > ~{sampleId}minimap2.bam
    elif [ ~{customAln} == true ]
    then
      minimap2 -a ~{customArgs} -t ~{threads} ~{reference} ~{reads} | samtools sort -@4 -m ~{sortMemgb}G > ~{sampleId}minimap2.bam
    else
      minimap2 -ax ~{mapMode} ~{reference} ~{reads} -k ~{kmerSize} -y -K ~{minibatchSize} -t ~{threads} ~{mdString} ~{eqxString} | samtools sort -@4 -m ~{sortMemgb}G > ~{sampleId}minimap2.bam
    fi

    samtools index -@ ~{threads} ~{sampleId}minimap2.bam
  >>>

  output {
        File bam = "~{sampleId}minimap2.bam"
        File bamIndex = "~{sampleId}minimap2.bam.bai"
  }

  runtime {
    docker: "mkolmogo/card_minimap2:2.23"
    cpu: threads
	memory: memSizeGb + " GB"
	disks: "local-disk " + diskSizeGb + " SSD"
    preemptible: preemptible
  }
}

task splitReads {
    input {
        File reads
        Int readsPerChunk
        Int threads = 8
        Int diskGb = 5 * round(size(reads, "G")) + 20
        Int preemptible = 2
    }

    Int gzThreads = if threads > 1 then threads - 1 else 1
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        CHUNK_LINES=$(( ~{readsPerChunk} * 4 ))

        MM_INPUT=~{reads}
        if [ "${MM_INPUT: -3}" == "bam" ]
        then
            samtools fastq -TMm,Ml,MM,ML ~{reads} | split -l $CHUNK_LINES --filter='pigz -p ~{gzThreads} > ${FILE}.fq.gz' - "fq_chunk.part."
        else
            gzip -cd ~{reads} | split -l $CHUNK_LINES --filter='pigz -p ~{gzThreads} > ${FILE}.fq.gz' - "fq_chunk.part."
        fi
    >>>
    output {
        Array[File] readChunks = glob("fq_chunk.part.*")
    }
    runtime {
        preemptible: preemptible
        time: 120
        cpu: threads
        memory: "4 GB"
        disks: "local-disk " + diskGb + " SSD"
        docker: "quay.io/jmonlong/minimap2_samtools:v2.24_v1.16.1_pigz"
    }
}

task mergeBAM {
    input {
        Array[File] bams
        String outname = "merged"
        Array[String] chrs = []
        Int threads = 8
        Int diskGb = round(5 * size(bams, 'G')) + 20
        Int memGb = 8
        Int preemptible = 2
    }

    Boolean anyChrs = length(chrs) > 0

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        samtools merge -f -p -c --threads ~{threads} ~{outname}.bam ~{sep=" " bams}
        samtools index -@ ~{threads} ~{outname}.bam

        ## split by chromosome, if any chrs specified
        if [ ~{anyChrs} == true ]
        then
            mkdir bamPerChrs
            while read -r chrn
            do
                samtools view -@ ~{threads} -h -O BAM ~{outname}.bam ${chrn} -o bamPerChrs/~{outname}.${chrn}.bam
                samtools index -@ ~{threads} bamPerChrs/~{outname}.${chrn}.bam
            done < ~{write_lines(chrs)}
        fi
    >>>
    output {
        File bam = "~{outname}.bam"
        File bamIndex = "~{outname}.bam.bai"
        Array[File]? bamPerChrs = glob("bamPerChrs/*.bam")
        Array[File]? bamPerChrsIndex = glob("bamPerChrs/*.bam.bai")
    }
    runtime {
        preemptible: preemptible
        time: 240
        memory: memGb + " GB"
        cpu: threads
        disks: "local-disk " + diskGb + " SSD"
        docker: "quay.io/jmonlong/minimap2_samtools:v2.24_v1.16.1_pigz"
    }
}

task indexBAM {
    input {
        File bam
        Array[String] chrs = []
        Int threads = 8
        Int diskGb = round(5 * size(bam, 'G')) + 20
        Int memGb = 8
    }
    Boolean anyChrs = length(chrs) > 0
    String outname = basename(bam, ".bam")
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ln -s ~{bam} reads.bam
        samtools index -@ ~{threads} reads.bam

        ## split by chromosome, if any chrs specified
        if [ ~{anyChrs} == true ]
        then
            mkdir bamPerChrs
            while read -r chrn
            do
                samtools view -@ ~{threads} -h -O BAM reads.bam ${chrn} -o bamPerChrs/~{outname}.${chrn}.bam
                samtools index -@ ~{threads} bamPerChrs/~{outname}.${chrn}.bam
            done < ~{write_lines(chrs)}
        fi
    >>>
    output {
        File bamIndex = "reads.bam.bai"
        Array[File]? bamPerChrs = glob("bamPerChrs/*.bam")
        Array[File]? bamPerChrsIndex = glob("bamPerChrs/*.bam.bai")
    }
    runtime {
        preemptible: 2
        time: 240
        memory: memGb + " GB"
        cpu: threads
        disks: "local-disk " + diskGb + " SSD"
        docker: "quay.io/jmonlong/minimap2_samtools:v2.24_v1.16.1_pigz"
    }
}

task mergeFASTQ {
    input {
        Array[File] reads
        String outname = "merged"
        Int threads = 8
        Int diskGb = round(3 * size(reads, 'G')) + 20
        Int memGb = 8
        Int preemptible = 2

    }

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        cat ~{sep=" " reads} > ~{outname}.fastq.gz

    >>>
    output {
        File fq = "~{outname}.fastq.gz"
    }
    runtime {
        preemptible: preemptible
        time: 240
        memory: memGb + " GB"
        cpu: threads
        disks: "local-disk " + diskGb + " SSD"
        docker: "quay.io/jmonlong/minimap2_samtools:v2.24_v1.16.1_pigz"
    }
}
