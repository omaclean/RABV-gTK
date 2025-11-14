

params.tax_id="11292"
params.db_name="rabv-jul0425"
params.master_acc="NC_001542"
params.is_segmented="N"
params.skip_fill=true
params.publish_dir="results"
params.email="your_email@example.com"
params.ref_list="${projectDir}/generic/rabv/ref_list.txt"

scripts_dir="${projectDir}/scripts"

params.test=0



// 1. List your script's explicitly defined parameters (keep this in sync!)
def scriptDefinedParams = [
    'tax_id', 'db_name', 'master_acc', 'is_segmented', 'skip_fill', 'test',
    "scripts_dir", "publish_dir", "email", "ref_list"
    // Add all parameter names defined above
]

// 2. List common/allowed built-in Nextflow parameters (customize as needed)
//    You might need to add others depending on features you use (e.g., cloud options)
def knownNextflowParams = [
    'help', 'version', 'profile', 'config', 'workDir', 'resume',
    'entry', 'validateParams', 'cached', 'dumpHashes', 'listInputs',
    // Add cloud-specific ones if used: 'awsqueue', 'clusterOptions', etc.
    // Add tower related ones if used: 'tower', 'computeEnvId', etc.
    // Check Nextflow documentation for a more exhaustive list if required
]

// 3. Combine allowed parameters
def allowedParams = (scriptDefinedParams + knownNextflowParams).unique()

// 4. Get all parameters provided via command line and config files
//    The 'params' map holds the merged view of all parameters
def providedParams = params.collect { it.key }

// 5. Find parameters that were provided but are not in the allowed list
def unexpectedParams = providedParams.findAll { ! (it in allowedParams) }

// 6. Error out if any unexpected parameters were found
if (unexpectedParams) {
    // Construct a helpful error message
    def errorMsg = """
    ERROR: Unknown command line parameter(s) provided: ${unexpectedParams.join(', ')}

    Please check your command. Allowed parameters are:
    Script parameters: ${scriptDefinedParams.join(', ')}
    Nextflow options : ${knownNextflowParams.join(', ')} (subset shown)
    """.stripIndent() // Use stripIndent for cleaner multi-line string formatting

    error(errorMsg) // Stop the pipeline
}

process FETCH_GENBANK{
    input:
        val (TAX_ID)
    output:
        path 'GenBank-XML', type: 'dir'
    shell:
    '''
    extra=""
    if( [ !{params.test} -eq 1 ] )
    then
        extra="--test_run"
    fi
    python !{scripts_dir}/GenBankFetcher.py --taxid !{TAX_ID} -b 50 \
            --update tmp/GenBank-matrix/gB_matrix_raw.tsv ${extra} -e !{params.email}
    #what's update doing with a tmp dir?

    '''
}

process DOWNLOAD_GFF{
    input:
        val master_acc
    output:
         GFF_FILE = file("*.gff"), emit: gff_file
    shell:
    '''
    python "!{scripts_dir}/DownloadGFF.py" --accession_ids "!{master_acc}"

    '''
}

//python GenBankParser.py
//python "${scripts_dir}/GenBankParser.py" -r "generic/rabv/ref_list.txt"

process GENBANK_PARSER{
    publish_dir "${params.publish_dir}"
    input:
         ref_list_path
         gen_bank_XML
    output:

    shell:
    '''
        python !{scripts_dir}/GenBankParser.py -r ref_list_path -d GenBank-XML
        python !{scripts_dir}/ValidateMatrix.py
    '''
}
process skip_fill{
    input:
         
    output:
         
    when:
        params.skip_fill.toBoolean()
    shell:
    '''
    echo "Skipping fill step as per user request."
    '''
}


workflow {
    FETCH_GENBANK(params.tax_id)
    DOWNLOAD_GFF(params.master_acc)
    GENBANK_PARSER(ref_list_path: params.ref_list, gen_bank_XML: FETCH_GENBANK.out)
    if(skip_fill){
        data=SKIP_FILL(FETCH_GENBANK.out)
    }else{
        data=FETCH_GENBANK()
    }
}
