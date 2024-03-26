// Modules
include { FASTQC }                      from '../modules/fastqc/main'
include { MULTIQC }                     from './../modules/multiqc/main'
include { MULTIQC as MULTIQC_RUN }      from './../modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './../modules/custom/dumpsoftwareversions'

/*
Read the sub folders and any fastQ files therein
*/
reads = Channel.fromPath("${params.input}/[A-Za-z0-9_-]*/*_R{1,2}*.fastq.gz").filter { f -> f.toString().contains('Undetermined') }

/*
Read the fastQ files and get the
parent folder name as project name
*/
reads.map { f ->
    def meta = [:]
    meta.id = f.getParent().toString().split('/')[-1]
    [
        meta ,file(f)
    ]
}.set { ch_reads_by_project }

ch_multiqc_config = params.multiqc_config   ? Channel.fromPath(params.multiqc_config, checkIfExists: true).collect() : Channel.value([])
ch_multiqc_logo   = params.multiqc_logo     ? Channel.fromPath(params.multiqc_logo, checkIfExists: true).collect() : Channel.value([])

/*
Read the sequencer stats for this run
*/
stats = file("${params.input}/Stats/Stats.json")
ch_stats = stats.exists() ? Channel.fromPath(stats).collect() : []

ch_versions = Channel.from([])
multiqc_files = Channel.from([])

workflow READQC {
    main:

    /*
    Perform basic read qc
    */
    FASTQC(
        ch_reads_by_project
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    multiqc_files = multiqc_files.mix(FASTQC.out.zip)

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    /*
    Group FastQC by meta hash and add the software versions
    */
    ch_multiqc_input = multiqc_files.groupTuple().cross(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml)

    /*
    Combine QC results
    */
    MULTIQC(
        ch_multiqc_input,
        ch_multiqc_config,
        ch_multiqc_logo
    )

    MULTIQC_RUN(
        ch_stats.cross(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml),
        ch_multiqc_config,
        ch_multiqc_logo
    )

    emit:
    qc = MULTIQC.out.html
}
