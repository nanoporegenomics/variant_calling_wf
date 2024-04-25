version 1.0

import "../tasks/modkit.wdl" as modkit_t

workflow modkit_wf
{
	input {
		File  haplotaggedBam	
		File  haplotaggedBamIdx
		File  referenceFasta
		File? regional_bed
		Boolean run_pileup = true
		Boolean run_probs = false
		String  sampleName = "sample"
		String  referenceName = "ref"
	}

	if (run_pileup){
			call modkit_t.modkit as modkit {
				input:
					haplotaggedBam = haplotaggedBam,
					haplotaggedBamBai = haplotaggedBamIdx,
					ref = referenceFasta,
					sample_name = sampleName,
					ref_name = referenceName,
					regional_bed = regional_bed
			}
    }

	if (run_probs){
		call modkit_t.plot_ML_hist {
            input:
				modbam=haplotaggedBam,
				modbam_index=haplotaggedBamIdx,
				sample=sampleName
        }
	}



	output {
		File? methylationBed1 = modkit.hap1bedOut
		File? methylationBed2 = modkit.hap2bedOut
		File? methylationBedUngrouped = modkit.ungroupedBedOut
		File? methylationBed = modkit.wholeGenomeOut
#		File? ML_histogram = plot_ML_hist.out
		File? ML_tsv = plot_ML_hist.out_ML_tsv
	}
}
