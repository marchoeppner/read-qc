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
include { INTEROP_SUMMARY }                         from './../modules/interop/summary'
/*
Import sub workflows
*/
include { CONTAMINATION }                           from './../subworkflows/contamination'

workflow READQC {
    main:

    ch_versions = channel.from([])
    multiqc_files = channel.from([])
    
    ch_multiqc_config = params.multiqc_config   ? channel.fromPath(params.multiqc_config, checkIfExists: true).collect() : channel.value([])
    ch_multiqc_logo   = params.multiqc_logo     ? channel.fromPath(params.multiqc_logo, checkIfExists: true).collect() : channel.value([])

    ch_reads = channel.empty()

    illumina_folder = channel.fromPath(params.input, checkIfExists: true).map { f -> 
        [
            [sample_id: file(params.input).getBaseName()],
            f
        ]
    }
    samplesheet = params.samplesheet ? channel.fromPath(params.samplesheet, checkIfExists: true) : channel.fromPath("${params.input}/SampleSheet.csv", checkIfExists: true)

    // Find the demuxed reads if any
    channel.fromPath(params.input + "/**/*.fastq.gz").ifEmpty(false).map { f ->
        [
            [sample_id: file(params.input).getBaseName()],
            f
        ]
    }.set { ch_reads }
    
    // Check if there are reads
    illumina_folder.merge(ch_reads).branch { m,f,n,r ->
        has_reads: r
            return [m, file(f)]
        has_no_reads: !r
            return [m, file(f)]
    }.set { illumina_folder_by_status }

    // bloom filter - built-in plus installed
    bloomfilters = [ "$baseDir/assets/bloomfilters/phix.bf", "$baseDir/assets/bloomfilters/univec_core.bf"]
    params.bloomfilter.keySet().each { b ->
        bloomfilters << params.bloomfilter[b].bf
    }
    ch_bloomfilters = bloomfilters.join(" ")

    /*
    Read the sub folders and any fastQ files therein
    */

    /*
    Run Illumina interop
    */
    INTEROP_SUMMARY(
        illumina_folder
    )
    ch_versions = ch_versions.mix(INTEROP_SUMMARY.out.versions)
    multiqc_files = multiqc_files.mix(INTEROP_SUMMARY.out.csv.map { m,t -> t})

    // Demux reads from scratch to obtain relevant metrics
    BCL2FASTQ(
        illumina_folder_by_status.has_no_reads.combine(samplesheet)
    )
    ch_versions = ch_versions.mix(BCL2FASTQ.out.versions)
    multiqc_files = multiqc_files.mix(BCL2FASTQ.out.stats.map {m,s -> s } )

    BCL2FASTQ.out.fastq.flatMap { m, fastqs ->
        fastqs.collect { fastq ->
            def meta = [:]
            meta.id = m.sample_id
            meta.sample_id = fastq.getSimpleName()
            tuple(meta,fastq)
        }
    }.set { ch_reads_demuxed }
    
    ch_reads = ch_reads.mix(ch_reads_demuxed)

    ch_reads.branch { m, r ->
        valid: r
        invalid: !r
    }.set { reads_by_status }

    // forward the illumina folder after demuxing; we use one of the outputs of bcl2fastq to trigger this
    illumina_folder.merge(reads_by_status.valid).map { m,f,n,r ->
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
    reads_by_status.valid.map { m,f ->
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
        ch_bloomfilters
    )
    ch_versions = ch_versions.mix(CONTAMINATION.out.versions)
    multiqc_files = multiqc_files.mix(CONTAMINATION.out.qc)

    /*
    Perform basic read qc
    */
    FASTQC(
        ch_grouped_reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    multiqc_files = multiqc_files.mix(FASTQC.out.zip.map {m,z -> z})

    // Compute md5sum and store reads by sequencing project
    MD5SUM(
        reads_by_status.valid.map { m,reads -> 
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