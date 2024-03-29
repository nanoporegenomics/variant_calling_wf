version 1.0

import "../tasks/sniffles.wdl" as sniffles_t

workflow snifflesWf {

    input {
        File bamAlignment
        File bamAlignmentIndex
        String sample
        File vntrAnnotations
        Boolean phaseVariants
        Int threads
    }

	### Sniffles
    call sniffles_t.sniffles_t as sniffles_t {
        input:
            threads=threads,
			bamAlignment=bamAlignment,
            bamAlignmentIndex=bamAlignmentIndex,
            sample=sample,
			vntrAnnotations=vntrAnnotations,
            phaseVariants=phaseVariants
    }

	output {
        File structuralVariantsSniffles = sniffles_t.snifflesVcf
        File snifflesSnf = sniffles_t.snifflesSnf
	}
}
