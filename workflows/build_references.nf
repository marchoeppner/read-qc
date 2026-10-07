include { UNTAR } from "./../modules/untar"

workflow BUILD_REFERENCES {

    main:

    ch_versions = channel.from([])

    ch_kraken_db = channel.fromPath(params.references.kraken2.url)

    UNTAR(
        ch_kraken_db.map { f ->
            [
                [ id: f.getBaseName() ],
                f
            ]
        }
    )
    ch_versions = ch_versions.mix(UNTAR.out.versions)

    emit:
    versions = ch_versions

}