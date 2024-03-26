//
// This file holds functions to validate user-supplied arguments
//

class WorkflowPipeline {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        if (!params.run_name) {
            log.info 'Must provide a run_name (--run_name)'
            System.exit(1)
        }
        if (!params.input) {
            log.info 'Pipeline requires a folder with sequencing data as input!'
            System.exit(1)
        }
    }

}

