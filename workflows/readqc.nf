/*
Import modules
*/
include { FASTQC }                                  from './../modules/fastqc'
include { MULTIQC }                                 from './../modules/multiqc'
include { MD5SUM }                                  from './../modules/md5sum'
include { CUSTOM_DUMPSOFTWAREVERSIONS }             from './../modules/custom/dumpsoftwareversions'
include { INTEROP_SUMMARY }                         from './../modules/interop/summary'
/*
Import sub workflows
*/
include { CONTAMINATION }                           from './../subworkflows/contamination'
include { BCL_DEMULTIPLEX }                         from './../subworkflows/bcl_demultiplex'

workflow READQC {

    main:

    ch_versions     = channel.from([])
    multiqc_files   = channel.from([])
    ch_reports      = channel.from([])

    ch_kraken_db = params.kraken_db ? channel.fromPath(params.kraken_db, checkIfExists: true) : channel.fromPath(params.references.kraken2.db, checkIfExists: true)

    ch_multiqc_config = params.multiqc_config   ? channel.fromPath(params.multiqc_config, checkIfExists: true).collect() : channel.value([])
    ch_multiqc_logo   = params.multiqc_logo     ? channel.fromPath(params.multiqc_logo, checkIfExists: true).collect() : channel.value([])

    ch_demux_reads  = channel.empty()
    ch_demux_dir    = channel.empty()
    ch_interop      = channel.empty()
    ch_stats        = channel.empty()

    // Turn input into a channel
    channel.fromPath(params.input, checkIfExists: true).map { f -> 
        [
            [id: file(params.input).getBaseName()],
            file(f)
        ]
    }.set { illumina_folder }
        
    samplesheet = params.samplesheet ? channel.fromPath(params.samplesheet, checkIfExists: true) : channel.fromPath("${params.input}/SampleSheet.csv", checkIfExists: true)

    /*
    Check if we should demultiplex or use existing reads
    */
    if (params.bclconvert || params.bcl2fastq) {

        BCL_DEMULTIPLEX(
            illumina_folder.merge(samplesheet)
        )
        ch_interop      = BCL_DEMULTIPLEX.out.interop
        ch_demux_dir    = BCL_DEMULTIPLEX.out.demuxed
        ch_stats        = BCL_DEMULTIPLEX.out.stats

        BCL_DEMULTIPLEX.out.reads.flatMap { m, fastqs ->
            fastqs.collect { fastq ->
                def meta = [:]
                meta.id = m.id
                meta.sample_id = fastq.getSimpleName()
                tuple(meta,fastq)
            }
        }.set { ch_demux_reads }
    
    } else {
        // Find the demuxed reads if any
        channel.fromPath(params.input + "/**/*.fastq.gz").ifEmpty(false).map { f ->
            [
                [sample_id: file(params.input).getBaseName()],
                f
            ]
        }.set { ch_demux_reads }

        ch_interop      = channel.fromPath("${params.input}/InterOp/*.bin")
        ch_demux_dir    = illumina_folder
    }

    /*
    Run Illumina interop summary
    */
    INTEROP_SUMMARY(
        ch_demux_dir
    )
    ch_versions = ch_versions.mix(INTEROP_SUMMARY.out.versions)
    multiqc_files = multiqc_files.mix(INTEROP_SUMMARY.out.csv.map { _m,t -> t})
    ch_reports = ch_reports.mix(INTEROP_SUMMARY.out.csv)

    /* 
    Perform contamination check on reads
    group by library ID and Lane
    */
    ch_demux_reads.map { _m,f ->
            def sample = sample_from_library(f.getBaseName())
            tuple(sample,f)
        }.groupTuple()
        .map { k,files ->
            def meta = [:]
            meta.sample_id = k
            meta.single_end = false
            tuple(meta,files)
        }.set { ch_grouped_reads }
        
    CONTAMINATION(
        ch_grouped_reads,
        ch_kraken_db.collect()
    )
    ch_versions     = ch_versions.mix(CONTAMINATION.out.versions)
    multiqc_files   = multiqc_files.mix(CONTAMINATION.out.qc)
    ch_reports      = ch_reports.mix(CONTAMINATION.out.qc)
    
    /*
    Perform basic read qc
    */
    FASTQC(
        ch_grouped_reads
    )
    ch_versions     = ch_versions.mix(FASTQC.out.versions)
    multiqc_files   = multiqc_files.mix(FASTQC.out.zip.map {_m,z -> z})

    // Compute md5sum and store reads by sequencing project
    MD5SUM(
        ch_demux_reads.map { m,reads -> 
            def meta = [:]
            def project = reads.getParent().getName()
            meta.sample_id = m.sample_id
            meta.project = project
            tuple(meta, reads)
        },
        true
    )

    ch_versions = ch_versions.mix(MD5SUM.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    ch_reports = ch_reports.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.yml.map { y -> [[id: params.run_name], y] })

    multiqc_files = multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml)
    /*
    Combine QC results
    */    
    MULTIQC(
        multiqc_files.collect(),
        ch_multiqc_config,
        ch_multiqc_logo
    )

    emit:
    qc = MULTIQC.out.html
}

def sample_from_library(lib)  {
    def sample = lib.split("_L00")[0]
    return sample
}