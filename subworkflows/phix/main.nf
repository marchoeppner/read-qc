include { BWAMEM2_MEM }     from "./../../modules/bwamem2/mem"
include { SAMTOOLS_MERGE }  from './../../modules/samtools/merge'
include { SAMTOOLS_INDEX }  from './../../modules/samtools/index'
include { SAMTOOLS_STATS }  from './../../modules/samtools/stats'

workflow PHIX {

    take:
    reads
    idx

    main:
    ch_multiqc = Channel.from([])
    ch_versions = Channel.from([])

    reads_with_index = reads.combine(idx)

    // align each set of reads to the phix genome
    BWAMEM2_MEM(
        reads_with_index
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)

    // merge all alignments
    SAMTOOLS_MERGE(
        BWAMEM2_MEM.out.bam.map { m,b ->
            def meta = [:]
            meta.sample_id = file(params.input).getBaseName()
            tuple(meta,b)
        }.groupTuple()
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    // index merged alignments
    SAMTOOLS_INDEX(
        SAMTOOLS_MERGE.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    // compute stats for the merged alignment (contains % aligned and number of mismatches)
    SAMTOOLS_STATS(
        SAMTOOLS_INDEX.out.bam
    )
    ch_multiqc = ch_multiqc.mix(SAMTOOLS_STATS.out.stats)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    emit:

    qc = ch_multiqc
    versions = ch_versions
}