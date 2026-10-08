include { BCL2FASTQ }   from './../../modules/bcl2fastq'
include { BCLCONVERT }  from './../../modules/bclconvert'

workflow BCL_DEMULTIPLEX {

    take:
    illumina_run

    main:
    ch_versions = channel.from([])
    ch_fastq    = channel.from([])
    ch_interop  = channel.empty()
    ch_stats    = channel.empty()
    ch_qc       = channel.from([])
    ch_demuxed  = channel.empty()

    if (params.bclconvert) {
        BCLCONVERT(
            illumina_run
        )
        ch_fastq = BCLCONVERT.out.fastq
        ch_interop = BCLCONVERT.out.interop
        ch_versions = ch_versions.mix(BCLCONVERT.out.versions)
        ch_demuxed = BCLCONVERT.out.demuxed
    } else if (params.bcl2fastq) {
        BCL2FASTQ(
            illumina_run
        )
        ch_fastq = BCL2FASTQ.out.fastq
        ch_interop = BCL2FASTQ.out.interop
        ch_versions = ch_versions.mix(BCL2FASTQ.out.versions)
        ch_stats = BCL2FASTQ.out.stats
        ch_demuxed = BCL2FASTQ.out.demuxed
    }



    emit:
    qc = ch_qc
    reads = ch_fastq
    stats = ch_stats
    interop = ch_interop
    demuxed = ch_demuxed
}