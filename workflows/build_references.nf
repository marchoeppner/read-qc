
include { BIOBLOOMTOOLS_MAKER }         from "./../modules/biobloomtools/maker"
include { GUNZIP  as GUNZIP_GENOME }    from "./../modules/gunzip"

workflow BUILD_REFERENCES {

    main:
    
    hosts = params.bloomfilter.keySet()

    host_files = []

    hosts.each { h ->
        if (params.bloomfilter[h].name) {
            host_files << [
                [sample_id: params.bloomfilter[h].name],
                file(params.bloomfilter[h].url, checkIfExists: true)
            ]
        }
    }

    ch_hosts = Channel.from(host_files)

    ch_hosts.branch { m, f ->
        gzipped: f.toString().contains(".gz")
        uncompressed: !f.toString().contains(".gz")
    }.set { genomes_by_compression }
    
    GUNZIP_GENOME(
        genomes_by_compression.gzipped
    )

    ch_genomes = GUNZIP_GENOME.out.gunzip.mix(genomes_by_compression.uncompressed)

    BIOBLOOMTOOLS_MAKER(
        ch_genomes
    )

}