include { FASTP }                        from './../../modules/fastp'
include { BIOBLOOMTOOLS_CATEGORIZER }    from './../../modules/biobloomtools/categorizer'


workflow CONTAMINATION {

    take:
    reads
    bloomfilters

    main:

    ch_versions = channel.from([])
    
    FASTP(
        reads
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    BIOBLOOMTOOLS_CATEGORIZER(
        FASTP.out.reads,
        bloomfilters
    )
    ch_versions = ch_versions.mix(BIOBLOOMTOOLS_CATEGORIZER.out.versions)

    emit:
    versions = ch_versions
    qc = BIOBLOOMTOOLS_CATEGORIZER.out.results
}