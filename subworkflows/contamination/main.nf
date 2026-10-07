include { FASTP }               from './../../modules/fastp'
include { KRAKEN2_KRAKEN2 }     from './../../modules/kraken2/kraken2'            

workflow CONTAMINATION {

    take:
    reads
    kraken_db

    main:

    ch_qc = channel.from([])
    ch_versions = channel.from([])
    
    FASTP(
        reads
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)
    ch_qc = ch_qc.mix(FASTP.out.json)

    KRAKEN2_KRAKEN2(
        reads,
        kraken_db,
        false,
        false
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)
    qc_qc = ch_qc.mix(KRAKEN2_KRAKEN2.out.report)
   
    emit:
    versions = ch_versions

    qc = ch_qc
}