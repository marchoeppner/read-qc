/*
Import modules
*/
include { FASTQC }                                  from './../modules/fastqc'
include { MULTIQC }                                 from './../modules/multiqc'
include { CHECKQC }                                 from './../modules/checkqc'
include { BCL2FASTQ }                               from './../modules/bcl2fastq'
include { MULTIQC as MULTIQC_RUN }                  from './../modules/multiqc'
include { MD5SUM }                                  from './../modules/md5sum'
include { CUSTOM_DUMPSOFTWAREVERSIONS }             from './../modules/custom/dumpsoftwareversions'
include { BIOBLOOMTOOLS_CATEGORIZER as CHECK_PHIX } from './../modules/biobloomtools/categorizer'

/*
Import sub workflows
*/
include { CONTAMINATION }               from './../subworkflows/contamination'
include { PHIX }                        from './../subworkflows/phix'

workflow READQC {
    main:

    ch_versions = Channel.from([])
    multiqc_files = Channel.from([])
    
    ch_multiqc_config = params.multiqc_config   ? Channel.fromPath(params.multiqc_config, checkIfExists: true).collect() : Channel.value([])
    ch_multiqc_logo   = params.multiqc_logo     ? Channel.fromPath(params.multiqc_logo, checkIfExists: true).collect() : Channel.value([])

    illumina_folder = Channel.fromPath(params.input, checkIfExists: true).map { f -> 
        [
            [sample_id: file(params.input).getBaseName()],
            f
        ]
    }
    samplesheet = params.samplesheet ? Channel.fromPath(params.samplesheet, checkIfExists: true) : Channel.fromPath("${params.input}/SampleSheet.csv", checkIfExists: true)

    // bloom filter - built-in plus installed
    bloomfilters = [ "$baseDir/assets/bloomfilters/phix.bf", "$baseDir/assets/bloomfilters/univec_core.bf"]
    params.bloomfilter.keySet().each { b ->
        bloomfilters << params.bloomfilter[b].bf
    }
    ch_bloomfilters = bloomfilters.join(" ")

    // The PhiX BWA index
    bwa_index = Channel.fromPath("${baseDir}/assets/bloomfilters/phix.fasta*").collect().toList()

    /*
    Read the sub folders and any fastQ files therein
    */

    // Demux reads from scratch to obtain relevant metrics
    BCL2FASTQ(
        illumina_folder.combine(samplesheet)
    )
    ch_versions = ch_versions.mix(BCL2FASTQ.out.versions)
    multiqc_files = multiqc_files.mix(BCL2FASTQ.out.stats).map {m,s -> s }

    BCL2FASTQ.out.fastq.flatMap { m, fastqs ->
        fastqs.collect { fastq ->
            def meta = [:]
            meta.id = m.sample_id
            meta.sample_id = fastq.getSimpleName()
            tuple(meta,fastq)
        }
    }.set { ch_reads}

    ch_undetermined = BCL2FASTQ.out.undetermined.map {m,f -> f}.flatten()
    ch_fastqs = BCL2FASTQ.out.fastq.map {m,f -> f}.flatten()
    ch_undetermined.concat(ch_fastqs).map { f ->
        def meta = [:]
        meta.sample_id = f.getBaseName().split("_L00")[0]
        tuple(meta,f)
    }.groupTuple().set { ch_all_reads }

    // Align reads against phix and compute fraction and error rate
    PHIX(
        ch_all_reads,
        bwa_index
    )
    ch_versions = ch_versions.mix(PHIX.out.versions)
    //multiqc_files = multiqc_files.mix(PHIX.out.qc.map{ m,s -> s}) // not useful to have in MultiQC as is

    // forward the illumina folder after demuxing; we use one of the outputs of bcl2fastq to trigger this
    illumina_folder.combine(BCL2FASTQ.out.versions).map { m,f,v ->
        tuple(m,f)
    }.set { demux_folder }

    // Perform basic check on demuxed data
    CHECKQC(
        demux_folder
    )
    multiqc_files = multiqc_files.mix(CHECKQC.out.json)
    ch_versions = ch_versions.mix(CHECKQC.out.versions)

    // Perform contamination check on reads
    // group by library ID and Lane
    CONTAMINATION(
        ch_reads.map { m,f ->
            def sample = sample_from_library(f.getBaseName())
            tuple(sample,f)
        }.groupTuple()
        .map { k,files ->
            def meta = [:]
            meta.sample_id = k
            meta.single_end = false
            tuple(meta,files)
        },
        ch_bloomfilters
    )
    ch_versions = ch_versions.mix(CONTAMINATION.out.versions)
    multiqc_files = multiqc_files.mix(CONTAMINATION.out.qc)

    /*
    Perform basic read qc
    */
    FASTQC(
        ch_reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    multiqc_files = multiqc_files.mix(FASTQC.out.zip.map {m,z -> z})

    // Compute md5sum and store reads by sequencing project
    MD5SUM(
        ch_reads.map { m,reads -> 
            def meta = [:]
            def project = reads.getParent().getName()
            meta.sample_id = m.sample_id
            meta.project = project
            tuple(meta, reads)
        },
        true
    )

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

def sample_from_library(lib)  {
    def sample = lib.split("_L00")[0]
    return sample
}