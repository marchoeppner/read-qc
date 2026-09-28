include { FASTP }                        from './../../modules/fastp'
include { BIOBLOOMTOOLS_CATEGORIZER }    from './../../modules/biobloomtools/categorizer'
include { BIOBLOOM_SUMMARY }             from './../../modules/helper/biobloom_summary'


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

    BIOBLOOMTOOLS_CATEGORIZER.out.results.map { t ->
        [
            [ id: params.run_name ], t
        ]
    }.groupTuple()
    .set { contamination_jsons }

    BIOBLOOM_SUMMARY(
        contamination_jsons
    )

    emit:
    versions = ch_versions
    fastp_json = FASTP.out.json
    biobloom_json = BIOBLOOM_SUMMARY.out.json
    qc = BIOBLOOMTOOLS_CATEGORIZER.out.results
}