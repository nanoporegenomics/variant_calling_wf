version 1.0

import "../tasks/sniffles.wdl" as sniffles_t

workflow snifflesWf {

    input {
        File bamAlignment
        File vntrAnnotations
        Int threads
        String sample = "sniffles"
    }

	### Sniffles
    call sniffles_t.sniffles_t as sniffles_t {
        input:
            threads=threads,
            sample=sample,
			bamAlignment=bamAlignment,
			vntrAnnotations=vntrAnnotations
    }

	output {
        File structuralVariantsSniffles = sniffles_t.snifflesVcf
        File snifflesSnf = sniffles_t.snifflesSnf
	}
}
