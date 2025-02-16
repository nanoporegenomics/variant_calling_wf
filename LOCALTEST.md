## Test locally

Uses the tiny test dataset in [test](test).

Test locally, for example using [miniwdl](https://github.com/chanzuckerberg/miniwdl):

```
## End-to-end workflow (no methylation)
miniwdl run --as-me  --copy-input-files wdl/workflows/cardEndToEndVcf.wdl -i test/test.input.endtoendfastq.json

## End-to-end workflow (no methylation) with multiple input FASTQs and cheaper by chunking reads (and preempting)
miniwdl run --as-me  --copy-input-files wdl/workflows/cardEndToEndVcf.wdl -i test/test.input.endtoendfastq.chunks.json

## End-to-end workflow (no methylation) starting with an available assembly (e.g. from shasta workflow below)
miniwdl run --as-me  --copy-input-files wdl/workflows/cardEndToEndVcf.wdl -i test/test.input.endtoendfastq.amb.json

## End-to-end workflow with uBAM input (for methylation)
miniwdl run --as-me  --copy-input-files wdl/workflows/cardEndToEndVcf.wdl -i test/test.input.endtoendfastq.ubam.json

## End-to-end workflow with multiple mapped BAM input and chunking
miniwdl run --as-me  --copy-input-files wdl/workflows/cardEndToEndVcf.wdl -i test/test.input.endtoendfastq.bam.json

## DV+margin workflow
miniwdl run --as-me  --copy-input-files wdl/workflows/dvMargin.wdl -i test/test.input.dvmargin.json

## shasta workflow
miniwdl run --as-me  --copy-input-files wdl/tasks/shasta.wdl -i test/test.input.shasta.json
miniwdl run --as-me  --copy-input-files wdl/tasks/shasta.wdl -i test/test.input.shasta.inmem.json

## regional methylation workflow
miniwdl run --as-me  -i wdl/tasks/regionalMethylation.wdl -i test/test/test.input.regionalMethylation.phased.json

```

### Next

- Simulate unmapped BAMs with methylation information
