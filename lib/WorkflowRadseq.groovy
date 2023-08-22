//
// This file holds several functions specific to the workflow/radseq.nf in the nf-core/radseq pipeline
//
import nextflow.Nextflow
import groovy.json.JsonSlurper

class WorkflowRadseq {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        if (!params.method) {
            log.error "type of workflow to execute not specified with e.g. '--method denovo' or via a detectable config file."
            System.exit(1)
        }
        if (params.method == 'reference') {
            if (!params.genome || params.genome == null) {
            log.error "need to specify a genome file with e.g. '--genome fasta' or via a detectable config file."
            System.exit(1)
            }
        }
        if (params.method == 'denovo'){
            if (!params.sequence_type) {
                log.error "need to specify the sequencing method with e.g. '--sequence_type' or via a detectable config file"
                System.exit(1)
            }
            if (!params.minreaddepth_withinindividual || params.minreaddepth_withinindividual == null) {
                log.warn("using default range of values for minReadDepth_withinIndividual")
            }
            if (params.method == 'denovo' && !params.minreaddepth_betweenindividual || params.minreaddepth_betweenindividual == null) {
                log.warn("using default range of values for minReadDepth_BetweenIndividual")
            }   
        }  
    }

    public static ArrayList groupMetaID(channel) {
        def meta = channel[0]
        def paths = channel[1]
        def metaf = [:] // initialize groovy map
        metaf.id =  meta.id.split(/[^\p{L}]/)[0] // set id splits at the first number appearance and retains items to the left
        metaf.single_end = meta.single_end
        metaf.ref_id = meta.ref_id

        [metaf, paths]
    } 

    public static ArrayList addRefIdToChannels(params, channel) {
        def metaf = [:]
        def meta = channel[0]
        def file = channel[1]
        def array = []
        // swtiched for two scenarios large and small channel
        if (channel.size() > 3) {
            // pair bam with correct index
            def meta2 = channel[2]
            def index = channel[3]
            
            if (params.method == 'denovo') {
                metaf.id = meta.id 
                metaf.umi_barcodes = meta.umi_barcodes
                metaf.single_end = meta.single_end
                metaf.ref_id = meta2.ref_id.tokenize( '_' )[0] + '_' + meta2.ref_id.tokenize( '_' )[1]

            } else {
                metaf.id = meta.id
                metaf.umi_barcodes = meta.umi_barcodes
                metaf.single_end = meta.single_end
                metaf.ref_id = meta2.ref_id
            }
            array = [ metaf, file, index ] 
        } else {
            // impacts fasta, fai channels
            if (params.method == 'denovo') {
                // pair d fasta
                metaf.id = meta.id
                metaf.single_end = meta.single_end
                if (file.toString().endsWith('.fasta')) {
                    metaf.ref_id = file.baseName.toString().tokenize('_')[1] + '_' + file.baseName.toString().tokenize('_')[2]
                } else {
                    metaf.ref_id = meta.ref_id
                }
                array = [ metaf, file ]
            } else {
                // pair reference fasta
                // have to include id, single_end, ref_id in meta for fai and fasta file
                def meta2 = channel[2]
                metaf.id = meta2.id.split(/[^\p{L}]/)[0]
                metaf.single_end = meta2.single_end
                if (file.toString().endsWith('.fasta')) {
                    metaf.ref_id = file.baseName
                } else {
                    metaf.ref_id = meta.ref_id
                }
                array = [metaf, file]
            }
        }
    }

    //
    // Function to check alignment rate
    //
    public static ArrayList getBamPercentMapped(params, flagstat_log) {
        def percent_aligned = 0
        // 46 + 0 mapped (100.00% : N/A)
        def pattern = /primary\s+mapped\s+\((\d+\.\d+)%\s*:/
        flagstat_log[1].eachLine { line ->
            def matcher = line =~ pattern
            if (matcher) {
                percent_aligned = matcher[0][1].toFloat()
                }
            }
            def pass = false
            if (percent_aligned >= params.min_percent_aligned.toFloat()) {
                pass = true
            } 
            return [ flagstat_log[0], percent_aligned, pass ]
        }

    //
    // Function to generate error if no reads after filtering
    //
    public static ArrayList getFastpReadsAfterFiltering(params, log, json_file) {
        // read json file
        def json = new JsonSlurper().parseText(json_file[1].text)
        
        // extract particular element
        def after_filtering = json.summary.after_filtering.total_reads.toInteger()
        
        def pass = false
        if (after_filtering >= params.min_reads_after_fastp.toInteger()) {
            pass = true
        } else {
            //log.warn("low read count after fastp filtering ${json_file[0]}.id ${after_filtering}")
            log.warn("low read count after fastp filtering")
        }
        return [json_file[0], after_filtering ,pass]
    }

    public static ArrayList splitBedFile(params, log, bed) {
        def meta = bed[0]
        def inputLines = bed[1]

        def outputFilePath  = inputLines.getParent() + '/split_bed_files'
        outputFilePath.mkdir()
        def outputFilePaths = []

        inputLines.eachLine { line ->
            def fields = line.split('\t')
            def chrom = fields[0]
            def startPos = fields[1] as int
            def endPos = fields[2] as int

            def outputFile = new File("${outputFilePath}/${chrom}_${startPos}_${endPos}.bed")
            outputFile.write(line)

            outputFilePaths << outputFile
        }
        return [meta, outputFilePaths]
    } 

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "=============================================================================\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "==================================================================================="
            System.exit(1)
        }
    }
}
