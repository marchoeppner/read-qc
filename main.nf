#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
===============================
Read QC Pipeline
===============================

This Pipeline performs basic QC of sequencing runs and project data therein

### Homepage / git
git@github.com:marchoeppner/read-qc.git

*/

// Pipeline version
params.version = workflow.manifest.version

include { READQC }              from './workflows/readqc'
include { BUILD_REFERENCES }    from './workflows/build_references'
include { PIPELINE_COMPLETION } from './subworkflows/pipeline_completion'


workflow {

    WorkflowMain.initialise(workflow, params, log)
    WorkflowPipeline.initialise(params, log)

    multiqc_report = channel.from([])

    if (params.build_references) {
        BUILD_REFERENCES()
    } else {
        READQC()
        multiqc_report = multiqc_report.mix(READQC.out.qc).toList()
    }

    PIPELINE_COMPLETION()
}