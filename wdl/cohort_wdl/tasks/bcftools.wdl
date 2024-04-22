version 1.0

workflow run_bcftools{
	input {
		File gvcfFile
		String out_name
		String? extraArgs = ""
		String dockerImage = "meredith705/truvari"

	}

	call filterVCF {
		input:
			gvcfFile = gvcfFile,
			out_name = out_name,
			dockerImage = dockerImage
	}

	output{
		File filtered_gvcf = filterVCF.filtered_gvcf
	}
}

task filterVCF {
	input {
		File gvcfFile
		String out_name
		Int minDepth = 1
		Int minGenoTQual = 1
		String? extraArgs = ""
		Int memSizeGB = 128
		Int threadCount = 64
		Int diskSizeGB = round(size(gvcfFile, 'G')) + 300
		String dockerImage = "meredith705/truvari"

	}

	command <<<

	bcftools view -i 'MIN(FMT/DP)>~{minDepth} & MIN(FMT/GQ)>~{minGenoTQual}' \
		--threads ~{threadCount} ~{extraArgs} ~{gvcfFile} \
		| bgzip -@ ~{threadCount} > ~{out_name}.filt.g.vcf.gz

	>>>

	output {
		File filtered_gvcf = "~{out_name}.filt.g.vcf.gz"
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + diskSizeGB + " SSD"
		docker: dockerImage
	}
}