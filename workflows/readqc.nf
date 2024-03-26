/*
Import modules
*/
include { FASTQC }                      from './../modules/fastqc/main'
include { MULTIQC }                     from './../modules/multiqc/main'
include { MULTIQC as MULTIQC_RUN }      from './../modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './../modules/custom/dumpsoftwareversions'

/*
Import sub workflows
*/
include { CONTAMINATION }               from './../subworkflows/contamination'

/*
Read the sub folders and any fastQ files therein
*/
reads = Channel.fromFilePairs("${params.input}/Alignment*/[0-9]*/Fastq/*_R{1,2}*.fastq.gz").filter { k,f -> !f[0].toString().contains('Undetermined') }

/*
Read the fastQ files and get the
parent folder name as project name
*/
reads.map { key, reads ->
    def meta = [:]
    meta.id = "MiSeq" //f.getParent().toString().split('/')[-1]
    meta.sample_id = key
    [
        meta, reads
    ]
}.set { ch_reads_by_project }

ch_multiqc_config = params.multiqc_config   ? Channel.fromPath(params.multiqc_config, checkIfExists: true).collect() : Channel.value([])
ch_multiqc_logo   = params.multiqc_logo     ? Channel.fromPath(params.multiqc_logo, checkIfExists: true).collect() : Channel.value([])

/*
Read the sequencer stats for this run
*/

ch_versions = Channel.from([])
multiqc_files = Channel.from([])

workflow READQC {
    main:

    /*
    Check for contaminations
    */
    CONTAMINATION(
        ch_reads_by_project
    )
    ch_versions = ch_versions.mix(CONTAMINATION.out.versions)
    multiqc_files = multiqc_files.mix(CONTAMINATION.out.qc)

    /*
    Perform basic read qc
    */
    FASTQC(
        ch_reads_by_project
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    multiqc_files = multiqc_files.mix(FASTQC.out.zip.map {m,z -> z})

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

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
