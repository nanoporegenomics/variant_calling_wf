version 1.0

import "../tasks/minimap2.wdl" as minimap2_t
import "../tasks/shasta.wdl" as shasta_t
import "../tasks/hapdup.wdl" as hapdup_t

workflow structuralVariantsDenovoAssembly {

    input {
        File readsFile
        File? shastaFasta
        File? shastaGFA
        File? shastaHTML
        File? shastaLOG
        File? readAlign
        Boolean shastaInMem = false
        String extraShastaArgs = ""
        Array[File] chunkedReadsFiles = []
        Int threads
        Int shastaDiskSizeGB = 1024
        Int hapdupDiskSizeGB = 1024
        Int preemptible = 2
    }

    
    ### Shasta assembly ###
    ## skip if a shastaFasta is 
    Array[File] readArray = [readsFile]
    if(!defined(shastaFasta)){
        if ((basename(readsFile, ".fasta") == basename(readsFile)) && (basename(readsFile, ".fa") == basename(readsFile))){
            call shasta_t.convertToFasta {
                input:
                readfiles=readArray,
                preemptible=preemptible
            }
        }
        File readsFasta = select_first([convertToFasta.fasta, readsFile])
        if(shastaInMem){
            call shasta_t.shasta_inmem_t as shasta_t {
                input:
                reads=readsFasta,
                shastaArgs=extraShastaArgs,
                diskSizeGb=shastaDiskSizeGB
            }
        }
        if(!shastaInMem){
            call shasta_t.shasta_t as shasta_inmem_t {
                input:
                reads=readsFasta,
                shastaArgs=extraShastaArgs,
                diskSizeGb=shastaDiskSizeGB
            }
        }
    }
    File ambFasta = select_first([shasta_t.shastaFasta, shasta_inmem_t.shastaFasta, shastaFasta])
    File? shastaGfa_t = select_first([shasta_t.shastaGfa, shasta_inmem_t.shastaGfa, shastaGFA])
    File? shastaLog_t = select_first([shasta_t.shastaLog, shasta_inmem_t.shastaLog, shastaLOG])
    File? shastaHtml_t = select_first([shasta_t.shastaHtml, shasta_inmem_t.shastaHtml, shastaHTML])

    if(!defined(readAlign)){
        ### minimap2 alignment ###
        if(length(chunkedReadsFiles) == 0){
            call minimap2_t.minimap2_t as minimap2 {
                input:
                reads = readsFile,
                reference=ambFasta,
                useEqx=false,
                threads = threads
            }
        }
        if(length(chunkedReadsFiles) > 0){
            scatter (readChunk in chunkedReadsFiles){
                call minimap2_t.minimap2_t as minimap2_chunk {
                    input:
                    reads = readChunk,
                    reference=ambFasta,
                    useEqx=false,
                    preemptible=preemptible,
                    threads = threads
                }
            }

            call minimap2_t.mergeBAM as mergeBAMhapdup {
                input:
                bams = minimap2_chunk.bam
            }
        }


    }
	

    File bamFile = select_first([minimap2.bam, mergeBAMhapdup.bam, readAlign])

	### hapdup
	call hapdup_t.hapdup_t as hapdup_t {
		input:
			threads=threads,
			alignedBam=bamFile,
			contigs=ambFasta,
			diskSizeGb=hapdupDiskSizeGB
    }

	output {
        File asmDual1 = hapdup_t.hapdupDual1
        File asmDual2 = hapdup_t.hapdupDual2
        File asmPhased1 = hapdup_t.hapdupPhased1
        File asmPhased2 = hapdup_t.hapdupPhased2
        File phaseBed1 = hapdup_t.hapdupPhaseBed1
        File phaseBed2 = hapdup_t.hapdupPhaseBed2 
		File shastaHaploid = ambFasta
        File? shastaGfa = shastaGfa_t
		File? shastaLog = shastaLog_t
        File? shastaHtml = shastaHtml_t
	}
}
